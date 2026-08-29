#!/usr/bin/env julia

using Radiant
using TOML
using Dates
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT = joinpath(
    ROOT,"qualification-results","julia-$(VERSION)","radiant-em-analytic.toml",
)

function main()
    photon = Photon()
    material = Material("analytic-absorber")
    Radiant.set_density(material,1.0)

    cross_sections = Cross_Sections()
    Radiant.set_source(cross_sections,"custom")
    Radiant.set_materials(cross_sections,material)
    Radiant.set_particles(cross_sections,photon)
    sigma_absorption = 0.2
    sigma_scattering = 0.0
    Radiant.set_absorption(cross_sections,[sigma_absorption])
    Radiant.set_scattering(cross_sections,[sigma_scattering])
    Radiant.set_energy_boundaries(cross_sections,[[1.0,0.001]])
    Radiant.build(cross_sections)

    midpoint_energy_MeV = Radiant.get_energies(cross_sections,photon)[1]
    energy_deposition_coefficient = sigma_absorption*midpoint_energy_MeV
    mcs = cross_sections.multigroup_cross_sections[1,1]
    Radiant.set_energy_deposition(mcs,[energy_deposition_coefficient,0.0])

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
    Radiant.set_particle(solver,photon)
    Radiant.set_solver_type(solver,"BTE")
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,1)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DD",1)

    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
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
        source_hash="analytic-em-dd-v1",
        provenance=Dict("classification" => "software-verification"),
    )
    source = Boundary_Angular_Current_Source(
        photon,[1],[0.0 0.0 0.0],[1.0],normal,
        [0.0 1.0 0.0],[0.0 0.0 -1.0],
        [1.0e3,1.0e6],directions,weights,angular_flux,normalization,
    )
    fixed_sources = Fixed_Sources(cross_sections,geometry,solvers)
    Radiant.add_source(fixed_sources,source)

    unit = Computation_Unit()
    Radiant.set_cross_sections(unit,cross_sections)
    Radiant.set_geometry(unit,geometry)
    Radiant.set_solvers(unit,solvers)
    Radiant.set_sources(unit,fixed_sources)
    Radiant.run(unit)

    calculated_flux = Radiant.get_flux(unit,photon)[1,1]
    calculated_deposition = Radiant.get_energy_deposition(unit,photon)[1]

    mu = directions[incoming_direction,1]
    delta_x = geometry.voxels_width["x"][1]
    incoming_flux = angular_flux[1,1,incoming_direction]
    expected_cell_angular_flux =
        (2.0*mu/delta_x)/(sigma_absorption+2.0*mu/delta_x)*incoming_flux
    expected_scalar_flux = weights[incoming_direction]*expected_cell_angular_flux
    expected_deposition = energy_deposition_coefficient*expected_scalar_flux

    isapprox(calculated_flux,expected_scalar_flux;rtol=1.0e-10,atol=1.0e-12) || error(
        "Diamond-difference analytic scalar-flux comparison failed.",
    )
    isapprox(calculated_deposition,expected_deposition;rtol=1.0e-10,atol=1.0e-12) || error(
        "Analytic energy-deposition comparison failed.",
    )
    calculated_deposition ≥ 0.0 || error("Energy deposition must be nonnegative.")

    receipt = Dict{String,Any}(
        "schema" => "radiant.em_analytic_qualification/v1",
        "classification" => "software-verification",
        "status" => "RADIANT_EM_ANALYTIC_PASS",
        "completed_at" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "source_hash" => normalization.source_hash,
        "spatial_scheme" => "DD1",
        "quadrature" => "Gauss-Legendre-2",
        "mu" => mu,
        "sigma_absorption_cm-1" => sigma_absorption,
        "expected_scalar_flux" => expected_scalar_flux,
        "calculated_scalar_flux" => calculated_flux,
        "expected_deposition" => expected_deposition,
        "calculated_deposition" => calculated_deposition,
        "physical_compound_data" => false,
        "geant4_reference" => false,
    )
    mkpath(dirname(OUTPUT))
    open(OUTPUT,"w") do io
        TOML.print(io,receipt)
    end
    println("RADIANT_EM_ANALYTIC_PASS")
    println("EM qualification receipt: $(OUTPUT)")
end

try
    main()
catch exception
    println(stderr,"RADIANT_EM_ANALYTIC_PASS not established.")
    showerror(stderr,exception,catch_backtrace())
    println(stderr)
    exit(1)
end
