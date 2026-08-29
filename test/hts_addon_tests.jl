using Radiant
using Test
using LinearAlgebra

function _process_scoring_fixture()
    particle = Photon()
    material = Material("fixture-material")
    Radiant.set_density(material,2.0)

    multigroup = Radiant.Multigroup_Cross_Sections(1)
    Radiant.set_total(multigroup,[0.0])
    Radiant.set_absorption(multigroup,[0.0])
    Radiant.set_boundary_stopping_powers(multigroup,[0.0,0.0])
    Radiant.set_stopping_powers(multigroup,[0.0])
    Radiant.set_momentum_transfer(multigroup,[0.0])
    Radiant.set_energy_deposition(multigroup,[5.0,0.0])
    Radiant.set_charge_deposition(multigroup,[0.0,0.0])
    Radiant.set_scattering(multigroup,zeros(Float64,1,1,1))
    add_response_channel!(
        multigroup,
        "energy-deposition|Photoelectric/fixture/in=photon/out=electron",
        [2.0,0.0],
    )
    add_response_channel!(
        multigroup,
        "energy-deposition|Compton/fixture/in=photon/out=photon",
        [3.0,0.0],
    )

    cross_sections = Cross_Sections()
    cross_sections.materials = [material]
    cross_sections.number_of_materials = 1
    cross_sections.particles = [particle]
    cross_sections.number_of_particles = 1
    cross_sections.number_of_groups = [1]
    cross_sections.energy_boundaries = [[1.0,0.001]]
    cross_sections.multigroup_cross_sections = reshape([multigroup],1,1)
    cross_sections.is_build = true

    geometry = Geometry()
    geometry.type = "cartesian"
    geometry.dimension = 1
    geometry.axis = ["x"]
    geometry.number_of_voxels["x"] = 1
    geometry.voxels_width["x"] = [1.0]
    geometry.voxels_position["x"] = [0.5]
    geometry.voxels_boundaries["x"] = [0.0,1.0]
    geometry.volume_per_voxel = [1.0]
    geometry.material_per_voxel = ones(Int64,1,1,1)
    geometry.is_build = true

    solver = SN()
    Radiant.set_particle(solver,particle)
    Radiant.set_solver_type(solver,"BTE")
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,1)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DD",1)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)

    sources = Fixed_Sources(cross_sections,geometry,solvers)
    sources.normalization_factor = 1.0

    flux_per_particle = Radiant.Flux_Per_Particle(particle)
    flux_array = zeros(Float64,1,1,1,1,1,1)
    flux_array[1,1,1,1,1,1] = 4.0
    Radiant.add_flux(flux_per_particle,flux_array)
    flux = Radiant.Flux()
    Radiant.add_flux(flux,flux_per_particle)
    return particle,cross_sections,geometry,solvers,sources,flux,multigroup
end

@testset "Temporary HTS add-on extraction boundary" begin
    manifest = default_hts_addon_extraction_manifest()
    @test validate_addon_extraction_manifest(manifest)
    @test manifest.intended_package == "RadiantHTS.jl"
    @test manifest.metadata["branch_policy"] == "branch-only-no-pr"
    @test get_addon_component(manifest,:microdosimetry).status == :ready_to_extract
    @test all(
        startswith(path,"src/hts_addon/")
        for component in manifest.components for path in component.files
    )
end

@testset "Native process-resolved scoring" begin
    particle,cross_sections,geometry,solvers,sources,flux,multigroup =
        _process_scoring_fixture()
    @test sum_response_channels(multigroup,"energy-deposition") == [5.0,0.0]
    @test length(get_response_channel_keys(multigroup;quantity="energy-deposition")) == 2
    score = score_process_responses(
        cross_sections,geometry,solvers,sources,flux,particle;
        quantity="energy-deposition",mass_normalized=true,
    )
    @test sort(get_process_keys(score)) == [
        "Compton/fixture/in=photon/out=photon",
        "Photoelectric/fixture/in=photon/out=electron",
    ]
    @test score.total[1,1,1] ≈ 10.0 atol=1.0e-12
    @test get_process_score(
        score,"Photoelectric/fixture/in=photon/out=electron",
    )[1,1,1] ≈ 4.0 atol=1.0e-12
    @test assert_process_score_closure(score,[10.0])
    aggregated = aggregate_process_scores(score)
    @test sort(get_process_keys(aggregated)) == ["Compton","Photoelectric"]
    @test aggregated.total == score.total
end

@testset "Spatial magnetic field maps" begin
    values = zeros(Float64,3,2,1,1)
    values[3,1,1,1] = 10.0
    values[3,2,1,1] = 12.0
    field = Cartesian_Magnetic_Field_Map(
        [0.0,1.0],[0.0],[0.0],values;
        interpolation=:trilinear,field_hash="analytic-field-map",
    )
    @test field_at(field,[0.5,0.0,0.0]) ≈ [0.0,0.0,11.0]
    local_frame = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test field_in_local_frame(field,[0.25,0.0,0.0],local_frame) ≈ [0.0,0.0,10.5]
    em = electromagnetic_field_at(field,[0.5,0.0,0.0])
    @test Radiant.get_magnetic_field(em) ≈ [0.0,0.0,11.0]
    @test field_map_receipt(field)["maximum_magnitude_T"] == 12.0
end

@testset "Piecewise-flat tape atlas" begin
    values = zeros(Float64,3,2,2,1)
    values[3,1,:,:] .= 10.0
    values[3,2,:,:] .= 12.0
    field = Cartesian_Magnetic_Field_Map(
        [0.0,2.0],[0.0,0.1],[0.0],values;
        field_hash="atlas-field-map",
    )
    centerline = [
        0.0 0.0 0.0
        1.0 0.0 0.0
        2.0 0.1 0.0
    ]
    atlas = build_piecewise_flat_tape_atlas(
        centerline;
        width_cm=0.4,thickness_cm=0.01632,initial_normal=[0.0,0.0,1.0],
        magnetic_field=field,source_patch_ids=[[1,2],[3,4]],
        geometry_hash="analytic-atlas",
    )
    @test length(atlas.patches) == 2
    @test atlas.total_length_cm > 2.0
    for patch in atlas.patches
        frame = patch_frame(patch)
        @test frame'*frame ≈ Matrix{Float64}(I,3,3) atol=1.0e-10
        local_point = [0.1,0.02,0.001]
        @test global_to_patch_coordinates(
            patch,patch_to_global_coordinates(patch,local_point),
        ) ≈ local_point atol=1.0e-12
        @test patch.magnetic_field_T[3] > 0.0
    end
    report = atlas_refinement_report(
        atlas;
        maximum_turn_angle_rad=1.0,
        maximum_sagitta_cm=1.0,
        maximum_field_variation_T=2.0,
    )
    @test report["screening_pass"]
    @test !report["production_pass"]
end

@testset "Gd prompt capture cascade source" begin
    cascade = synthetic_gd157_capture_fixture()
    rates = Capture_Rate_Field(
        "Gd-157",[1,2],[1.0,2.0],[2.0,4.0];
        variance_per_s2=[0.04,0.16],
        material_tags=["GdBCO","GdBCO"],
        transport_artifact_hash="self-shielded-openmc-fixture",
        provenance=Dict("capture_rates_are_self_shielded" => "true"),
    )
    edges = [0.0,1.0e4,1.0e5,1.0e6,1.0e7]
    bundle = build_gd_capture_source_bundle(
        cascade,rates;
        photon_energy_edges_eV=edges,electron_energy_edges_eV=edges,
    )
    @test Set(keys(bundle.sources)) == Set([
        :photon,:xray,:conversion_electron,:auger_electron,
    ])
    @test bundle.recoil_source !== nothing
    @test abs(bundle.energy_residual_eV_per_capture) ≤ 1.0e-8
    @test sum(get_volume_source_rate(bundle.sources[:photon])) ≈ 12.0
    @test sum(get_volume_source_rate(bundle.sources[:conversion_electron])) ≈ 3.0
    @test bundle.provenance["correlations_preserved_in_scalar_sources"] == "false"
end

@testset "Sub-keV and non-equilibrium thermalization" begin
    kernel = synthetic_subkev_kernel_fixture()
    result = thermalize_subkev_event(kernel,100.0)
    @test is_subkev_partition_closed(result)
    @test is_energy_partition_closed(result.base_partition)
    @test !subkev_model_is_production_ready(kernel)
    model = NonEquilibrium_Thermalization_Model(
        NonEquilibrium_Decay_Channel(
            :athermal_phonon,1.0e-3;
            heat_fraction=0.8,escape_fraction=0.2,
        ),
        NonEquilibrium_Decay_Channel(
            :quasiparticle,2.0e-3;
            heat_fraction=1.0,
        );
        model_hash="synthetic-nonequilibrium",qualification_status=:verification,
    )
    initial = non_equilibrium_release(result,model,0.0)
    late = non_equilibrium_release(result,model,1.0)
    @test initial.heat_eV == 0.0
    @test initial.remaining_non_equilibrium_eV ≈ get_non_equilibrium_energy(result)
    @test late.remaining_non_equilibrium_eV ≤ 1.0e-10
    @test late.heat_eV+late.escaped_eV+late.optical_eV ≈
          get_non_equilibrium_energy(result) atol=1.0e-8
end

@testset "Weighted statistical microdosimetry" begin
    kernel = synthetic_microdosimetry_kernel_fixture()
    events = sample_microdosimetry_events(
        kernel,100;time_window_s=(0.0,2.0),seed=17,
    )
    @test length(events) == 100
    @test sum(getfield.(events,:statistical_weight_events)) ≈ 200.0
    @test microdosimetry_effective_sample_size(events) ≈ 100.0
    @test all(event.deposited_energy_eV ≤ event.incident_energy_eV for event in events)
    @test all(is_energy_partition_closed(event.partition) for event in events)
    @test all(length(event.correlated_secondaries) == 1 for event in events)
    @test 0.0 ≤ detector_trigger_probability(events,500.0) ≤ 1.0
    @test weighted_deposited_energy_mean(events) > 0.0
    @test weighted_deposited_energy_variance(events) ≥ 0.0
    @test expected_specific_energy_Gy(events,1.0e-9) > 0.0
end

@testset "Cryogenic electrothermal coupling" begin
    heat_capacity = constant_cryogenic_property(
        1.0;units="J/K",property_hash="fixture-C",
    )
    bath_conductance = constant_cryogenic_property(
        0.1;units="W/K",property_hash="fixture-G",
    )
    node = Cryogenic_Thermal_Node(
        "active","YBCO",heat_capacity,bath_conductance;
        initial_temperature_K=10.0,
    )
    resistance = HTS_Transition_Resistance_Model(
        normal_resistance_ohm=1.0,
        residual_resistance_ohm=0.0,
        critical_temperature_K=20.0,
        transition_width_K=1.0,
        model_hash="fixture-R",
        qualification_status=:verification,
    )
    circuit = Electrothermal_Circuit(mode=:current_bias,bias_current_A=0.01)
    model = Cryogenic_Electrothermal_Model(
        [node],Cryogenic_Thermal_Link[],"active",resistance,circuit;
        bath_temperature_K=10.0,geometry_hash="fixture-geometry",
        material_state_hash="fixture-material",model_hash="fixture-electrothermal",
        qualification_status=:verification,
    )
    impulse = Electrothermal_Energy_Impulse(
        0.0,"active",1.0e-2;
        source_hash="fixture-process-partition",
    )
    result = simulate_cryogenic_electrothermal(
        model,[0.0,1.0e-3,1.0e-2,1.0e-1];impulses=[impulse],
    )
    @test peak_active_temperature(result,"active") > 10.0
    @test peak_voltage(result) ≥ 0.0
    @test maximum(abs.(result.step_energy_residual_J)) ≤ 1.0e-8
    @test result.source_hashes == ["fixture-process-partition"]
end
