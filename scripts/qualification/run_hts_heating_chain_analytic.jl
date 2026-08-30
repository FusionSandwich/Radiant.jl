#!/usr/bin/env julia

using Radiant
using TOML
using Dates
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT = joinpath(
    ROOT,"qualification-results","julia-$(VERSION)","hts-heating-chain-analytic.toml",
)
const EV_TO_J = 1.602176634e-19

function _build_geometry()
    geometry = Geometry()
    geometry.type = "cartesian"
    geometry.dimension = 1
    geometry.axis = ["x"]
    geometry.number_of_regions["x"] = 2
    geometry.voxels_per_region["x"] = [1,1]
    geometry.region_boundaries["x"] = [0.0,0.5,1.0]
    geometry.number_of_voxels["x"] = 2
    geometry.voxels_width["x"] = [0.5,0.5]
    geometry.voxels_position["x"] = [0.25,0.75]
    geometry.voxels_boundaries["x"] = [0.0,0.5,1.0]
    geometry.boundary_conditions["X-"] = "void"
    geometry.boundary_conditions["X+"] = "void"
    geometry.volume_per_voxel = reshape([0.5,0.5],2,1,1)
    geometry.material_per_voxel = reshape(Int64[1,2],2,1,1)
    geometry.is_build = true
    return geometry
end

function _build_cross_sections(photon::Particle)
    first_material = Material("analytic-layer-1")
    second_material = Material("analytic-layer-2")
    Radiant.set_density(first_material,1.0)
    Radiant.set_density(second_material,1.0)

    sigma_absorption = [0.10,0.20]
    cross_sections = Cross_Sections()
    Radiant.set_source(cross_sections,"custom")
    Radiant.set_materials(cross_sections,[first_material,second_material])
    Radiant.set_particles(cross_sections,photon)
    Radiant.set_absorption(cross_sections,sigma_absorption)
    Radiant.set_scattering(cross_sections,[0.0,0.0])
    Radiant.set_energy_boundaries(cross_sections,[[1.0,0.001]])
    Radiant.build(cross_sections)

    midpoint_energy_MeV = Radiant.get_energies(cross_sections,photon)[1]
    process_fractions = [0.75,0.35]
    for material_index in 1:2
        coefficient = sigma_absorption[material_index]*midpoint_energy_MeV
        mcs = cross_sections.multigroup_cross_sections[1,material_index]
        Radiant.set_energy_deposition(mcs,[coefficient,0.0])
        Radiant.add_response_channel!(
            mcs,
            "energy-deposition|photoelectric-absorption",
            [process_fractions[material_index]*coefficient,0.0],
        )
        Radiant.add_response_channel!(
            mcs,
            "energy-deposition|secondary-electron-thermalization",
            [(1.0-process_fractions[material_index])*coefficient,0.0],
        )
        reconstructed = Radiant.sum_response_channels(mcs,"energy-deposition")
        isapprox(reconstructed,[coefficient,0.0];rtol=1.0e-13,atol=1.0e-15) || error(
            "Process-resolved coefficients do not reconstruct the total deposition response.",
        )
    end
    return cross_sections,sigma_absorption,midpoint_energy_MeV
end

function _build_solver(photon::Particle)
    solver = SN()
    Radiant.set_particle(solver,photon)
    Radiant.set_solver_type(solver,"BTE")
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,1)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DD",1)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
    return solver,solvers
end

function _build_boundary_source(
    photon::Particle,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN,
)
    _,weights,directions,_ = Radiant._solver_quadrature(solver,geometry)
    normal = [-1.0 0.0 0.0]
    angular_flux = zeros(Float64,1,1,length(weights))
    incoming_direction = 0
    for direction in eachindex(weights)
        if dot(view(directions,direction,:),view(normal,1,:)) < 0.0
            incoming_direction = direction
            angular_flux[1,1,direction] = 2.0
        end
    end
    incoming_direction > 0 || error("No positive-x incoming ordinate was found.")
    normalization = Source_Normalization(
        basis=:per_history,
        source_rate_per_s=1.0,
        symmetry_factor=1.0,
        source_hash="analytic-two-layer-heating-v1",
        provenance=Dict(
            "classification" => "software-verification",
            "cross_section_area_cm2" => "1.0",
        ),
    )
    source = Boundary_Angular_Current_Source(
        photon,[1],[0.0 0.0 0.0],[1.0],normal,
        [0.0 1.0 0.0],[0.0 0.0 -1.0],
        [1.0e3,1.0e6],directions,weights,angular_flux,normalization,
    )
    projected,projection_receipt = project_boundary_source(
        source,cross_sections,geometry,solver,
    )
    projection_receipt.max_relative_error ≤ 1.0e-12 || error(
        "Boundary-current projection did not close.",
    )
    return source,weights,directions,incoming_direction,angular_flux,projection_receipt
end

function _diamond_difference_reference(
    sigma_absorption::Vector{Float64},
    widths::Vector{Float64},
    mu::Float64,
    weight::Float64,
    incident_angular_flux::Float64,
)
    cell_scalar_flux = zeros(Float64,length(widths))
    outgoing_angular_flux = incident_angular_flux
    for cell in eachindex(widths)
        incoming = outgoing_angular_flux
        streaming = 2.0*mu/widths[cell]
        cell_angular_flux = streaming/(sigma_absorption[cell]+streaming)*incoming
        outgoing_angular_flux = 2.0*cell_angular_flux-incoming
        outgoing_angular_flux ≥ -1.0e-14 || error(
            "Analytic fixture entered the negative diamond-difference regime.",
        )
        cell_scalar_flux[cell] = weight*cell_angular_flux
    end
    return cell_scalar_flux,outgoing_angular_flux
end

function _partition_layer_energy(label::String,available_energy_eV::Float64)
    partition = Energy_Partition(
        label=label,
        available_energy_eV=available_energy_eV,
        ionization_eV=0.15*available_energy_eV,
        electronic_excitation_eV=0.15*available_energy_eV,
        prompt_lattice_heat_eV=0.60*available_energy_eV,
        defect_stored_eV=0.05*available_energy_eV,
        escaping_particle_eV=0.03*available_energy_eV,
        cutoff_handoff_eV=0.02*available_energy_eV,
        metadata=Dict(
            "classification" => "analytic-partition-fixture",
            "instantaneous_heat_is_not_total_deposition" => "true",
        ),
    )
    is_energy_partition_closed(partition;rtol=1.0e-13,atol_eV=1.0e-9) || error(
        "Layer energy partition does not close.",
    )
    return partition
end

function _electrothermal_response(prompt_heat_J::Float64,source_hash::String)
    heat_capacity_J_K = 1.0e-12
    node = Cryogenic_Thermal_Node(
        "active-ybco","YBCO-analytic",
        constant_cryogenic_property(
            heat_capacity_J_K;units="J/K",property_hash="analytic-C",
        ),
        constant_cryogenic_property(
            0.0;units="W/K",property_hash="analytic-zero-G",allow_zero=true,
        );
        initial_temperature_K=20.0,
        metadata=Dict("classification" => "software-verification"),
    )
    resistance = HTS_Transition_Resistance_Model(
        normal_resistance_ohm=1.0,
        residual_resistance_ohm=0.0,
        critical_temperature_K=100.0,
        transition_width_K=1.0,
        model_hash="analytic-transition",
        qualification_status=:verification,
    )
    circuit = Electrothermal_Circuit(mode=:current_bias,bias_current_A=1.0e-12)
    model = Cryogenic_Electrothermal_Model(
        [node],Cryogenic_Thermal_Link[],"active-ybco",resistance,circuit;
        bath_temperature_K=20.0,
        magnetic_field_T=0.0,
        geometry_hash="analytic-two-layer-geometry",
        material_state_hash="analytic-ybco-state",
        model_hash="analytic-heating-chain",
        qualification_status=:verification,
    )
    impulse = Electrothermal_Energy_Impulse(
        0.0,"active-ybco",prompt_heat_J;
        source_hash=source_hash,
        channel=:prompt_lattice_heat,
    )
    result = simulate_cryogenic_electrothermal(
        model,[0.0,1.0e-6,2.0e-6];
        impulses=[impulse],energy_balance_rtol=1.0e-10,
        energy_balance_atol_J=1.0e-24,
    )
    expected_peak_K = 20.0+prompt_heat_J/heat_capacity_J_K
    observed_peak_K = peak_active_temperature(result,"active-ybco")
    isapprox(observed_peak_K,expected_peak_K;rtol=1.0e-10,atol=1.0e-12) || error(
        "Electrothermal temperature rise does not match Q/C for the isolated node.",
    )
    maximum(abs.(result.step_energy_residual_J)) ≤ 1.0e-22 || error(
        "Electrothermal energy balance did not close.",
    )
    return result,expected_peak_K,observed_peak_K
end

function main()
    photon = Photon()
    geometry = _build_geometry()
    cross_sections,sigma_absorption,midpoint_energy_MeV =
        _build_cross_sections(photon)
    solver,solvers = _build_solver(photon)
    source,weights,directions,incoming_direction,angular_flux,projection_receipt =
        _build_boundary_source(photon,cross_sections,geometry,solver)

    fixed_sources = Fixed_Sources(cross_sections,geometry,solvers)
    Radiant.add_source(fixed_sources,source)
    unit = Computation_Unit()
    Radiant.set_cross_sections(unit,cross_sections)
    Radiant.set_geometry(unit,geometry)
    Radiant.set_solvers(unit,solvers)
    Radiant.set_sources(unit,fixed_sources)
    Radiant.run(unit)

    calculated_flux = vec(Radiant.get_flux(unit,photon)[1,:])
    calculated_deposition = vec(Radiant.get_energy_deposition(unit,photon))
    mu = directions[incoming_direction,1]
    expected_flux,_ = _diamond_difference_reference(
        sigma_absorption,geometry.voxels_width["x"],mu,
        weights[incoming_direction],angular_flux[1,1,incoming_direction],
    )
    expected_deposition = sigma_absorption.*midpoint_energy_MeV.*expected_flux
    isapprox(calculated_flux,expected_flux;rtol=1.0e-10,atol=1.0e-12) || error(
        "Two-layer scalar flux does not match the diamond-difference reference.",
    )
    isapprox(calculated_deposition,expected_deposition;rtol=1.0e-10,atol=1.0e-12) || error(
        "Two-layer energy deposition does not match the analytic reference.",
    )

    process_score = score_process_responses(
        cross_sections,geometry,solvers,fixed_sources,unit.flux,photon;
        quantity="energy-deposition",mass_normalized=true,
    )
    assert_process_score_closure(
        process_score,reshape(calculated_deposition,2,1,1);
        rtol=1.0e-12,atol=1.0e-14,
    )

    densities = Radiant.get_densities(cross_sections)
    volumes = vec(geometry.volume_per_voxel)
    layer_energy_eV = calculated_deposition .* densities .* volumes .* 1.0e6
    partitions = [
        _partition_layer_energy("layer-$(index)",layer_energy_eV[index])
        for index in eachindex(layer_energy_eV)
    ]
    prompt_heat_J = sum(partition.prompt_lattice_heat_eV for partition in partitions)*EV_TO_J
    total_deposited_J = sum(layer_energy_eV)*EV_TO_J
    prompt_heat_J < total_deposited_J || error(
        "Analytic heating ledger incorrectly treated all deposited energy as prompt heat.",
    )
    electrothermal,expected_peak_K,observed_peak_K =
        _electrothermal_response(prompt_heat_J,source.normalization.source_hash)

    receipt = Dict{String,Any}(
        "schema" => "radiant.hts.heating_chain_analytic/v1",
        "classification" => "software-verification",
        "status" => "HTS_HEATING_CHAIN_ANALYTIC_PASS",
        "completed_at" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "source_hash" => source.normalization.source_hash,
        "boundary_projection_max_relative_error" => projection_receipt.max_relative_error,
        "sigma_absorption_cm-1" => sigma_absorption,
        "cell_width_cm" => geometry.voxels_width["x"],
        "cell_volume_cm3_for_unit_area" => volumes,
        "expected_scalar_flux" => expected_flux,
        "calculated_scalar_flux" => calculated_flux,
        "expected_deposition_MeV_per_g" => expected_deposition,
        "calculated_deposition_MeV_per_g" => calculated_deposition,
        "process_keys" => get_process_keys(process_score),
        "process_sum_closure" => true,
        "layer_available_energy_eV" => layer_energy_eV,
        "prompt_heat_fraction" => 0.60,
        "prompt_heat_J" => prompt_heat_J,
        "total_deposited_J" => total_deposited_J,
        "expected_peak_temperature_K" => expected_peak_K,
        "observed_peak_temperature_K" => observed_peak_K,
        "maximum_electrothermal_energy_residual_J" =>
            maximum(abs.(electrothermal.step_energy_residual_J)),
        "physical_material_data" => false,
        "physical_openmc_reference" => false,
        "physical_opensn_reference" => false,
        "physical_gate_promoted" => false,
    )
    mkpath(dirname(OUTPUT))
    open(OUTPUT,"w") do io
        TOML.print(io,receipt)
    end
    println("HTS_HEATING_CHAIN_ANALYTIC_PASS")
    println("Heating-chain receipt: $(OUTPUT)")
end

try
    main()
catch exception
    println(stderr,"HTS_HEATING_CHAIN_ANALYTIC_PASS not established.")
    showerror(stderr,exception,catch_backtrace())
    println(stderr)
    exit(1)
end
