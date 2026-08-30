#!/usr/bin/env julia

using Radiant
using TOML
using SHA
using LinearAlgebra
using Dates

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const FIXTURE = joinpath(ROOT,"test","fixtures","opensn_radiant_boundary_replay_v1.toml")
const OUTPUT = joinpath(
    ROOT,"qualification-results","julia-$(VERSION)","opensn-radiant-replay.toml",
)

function file_sha256(path::AbstractString)
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function matrix_row(values)
    result = zeros(Float64,1,length(values))
    result[1,:] .= Float64.(values)
    return result
end

function main()
    data = TOML.parsefile(FIXTURE)
    data["classification"] == "software-verification" || error(
        "Only explicitly classified software-verification fixtures are accepted here.",
    )

    photon = Photon()
    material = Material("analytic-vacuum")
    Radiant.set_density(material,1.0)
    cross_sections = Cross_Sections()
    Radiant.set_source(cross_sections,"custom")
    Radiant.set_materials(cross_sections,material)
    Radiant.set_particles(cross_sections,photon)
    Radiant.set_absorption(cross_sections,[0.0])
    Radiant.set_scattering(cross_sections,[0.0])
    source_edges_eV = Float64.(data["energy_edges_eV"])
    Radiant.set_energy_boundaries(cross_sections,[reverse(source_edges_eV)./1.0e6])
    Radiant.build(cross_sections)

    geometry = Geometry()
    geometry.type = "cartesian"
    geometry.dimension = 1
    geometry.axis = ["x"]
    geometry.number_of_regions["x"] = 1
    geometry.voxels_per_region["x"] = [1]
    geometry.region_boundaries["x"] = [0.0,1.0]
    geometry.number_of_voxels["x"] = 1
    geometry.voxels_width["x"] = [1.0]
    geometry.voxels_position["x"] = [0.5]
    geometry.voxels_boundaries["x"] = [0.0,1.0]
    geometry.boundary_conditions["X-"] = "void"
    geometry.boundary_conditions["X+"] = "void"
    geometry.volume_per_voxel = [1.0]
    geometry.material_per_voxel = ones(Int64,1,1,1)
    geometry.is_build = true

    boundary_conditions = Radiant.get_boundary_conditions(geometry)
    boundary_conditions == [0,0] || error(
        "Replay fixture boundary ordering is inconsistent with the Radiant 1-D sweep.",
    )

    solver = SN()
    Radiant.set_particle(solver,photon)
    Radiant.set_solver_type(solver,"BTE")
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,1)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DD",1)

    solvers = Solvers()
    Radiant.add_solver(solvers,solver)

    direction_records = data["direction"]
    directions = zeros(Float64,length(direction_records),3)
    angular_flux = zeros(Float64,1,1,length(direction_records))
    for index in eachindex(direction_records)
        directions[index,:] .= Float64.(direction_records[index]["cosines"])
        angular_flux[1,1,index] = Float64(direction_records[index]["angular_flux"][1])
    end
    patch = data["patch"]
    normalization = Source_Normalization(
        basis=Symbol(data["normalization_basis"]),
        source_rate_per_s=Float64(data["source_rate_per_s"]),
        symmetry_factor=Float64(data["symmetry_factor"]),
        source_hash=String(data["source_hash"]),
        provenance=Dict(
            "producer" => String(data["producer"]),
            "fixture_sha256" => file_sha256(FIXTURE),
        ),
    )
    source = Boundary_Angular_Current_Source(
        photon,
        [Int64(patch["id"])],
        matrix_row(patch["centroid_cm"]),
        [Float64(patch["area_cm2"])],
        matrix_row(patch["normal"]),
        matrix_row(patch["tangent_1"]),
        matrix_row(patch["tangent_2"]),
        source_edges_eV,
        directions,
        Float64.(data["quadrature_weights"]),
        angular_flux,
        normalization,
    )

    expected = Float64(data["expected_total_incoming_current"])
    observed = get_total_incoming_current(source)
    isapprox(observed,expected;rtol=1.0e-12,atol=1.0e-14) || error(
        "Fixture incoming-current invariant failed.",
    )

    projected,projection_receipt = project_boundary_source(
        source,cross_sections,geometry,solver,
    )
    projection_receipt.max_relative_error ≤ 1.0e-10 || error(
        "OpenSn-to-Radiant source projection failed current closure.",
    )

    fixed_sources = Fixed_Sources(cross_sections,geometry,solvers)
    Radiant.add_source(fixed_sources,source)
    Radiant.build(fixed_sources)
    length(get_projection_receipts(fixed_sources)) == 1 || error(
        "Exactly one boundary projection receipt was expected.",
    )
    sum(projected[1,:,1]) > 0.0 || error("Projected source was not installed on X-.")

    receipt = Dict{String,Any}(
        "schema" => "opensn_radiant_replay_receipt/v1",
        "classification" => "software-verification",
        "status" => "OPENSN_RADIANT_INTERFACE_SOFTWARE_PASS",
        "completed_at" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "fixture_sha256" => file_sha256(FIXTURE),
        "source_hash" => normalization.source_hash,
        "target_total_incoming_current" => expected,
        "observed_total_incoming_current" => observed,
        "projection_max_relative_error" => projection_receipt.max_relative_error,
        "energy_group_map" => projection_receipt.energy_group_map,
        "direction_map" => projection_receipt.direction_map,
        "surface_index" => projection_receipt.surface_index,
        "boundary_conditions_encoded" => boundary_conditions,
        "physical_openmc_source" => false,
        "physical_opensn_solution" => false,
    )
    mkpath(dirname(OUTPUT))
    open(OUTPUT,"w") do io
        TOML.print(io,receipt)
    end
    println("OPENSN_RADIANT_INTERFACE_SOFTWARE_PASS")
    println("Replay receipt: $(OUTPUT)")
end

try
    main()
catch exception
    println(stderr,"OPENSN_RADIANT_INTERFACE_SOFTWARE_PASS not established.")
    showerror(stderr,exception,catch_backtrace())
    println(stderr)
    exit(1)
end
