#!/usr/bin/env julia

using Radiant
using SHA
using TOML

length(ARGS) == 2 || error(
    "Usage: validate_cross_language_source_fixture.jl <fixture-directory> <receipt.toml>",
)
fixture_directory = abspath(ARGS[1])
receipt_path = abspath(ARGS[2])

function file_sha256(path::AbstractString)
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

boundary_path = joinpath(fixture_directory,"python_boundary_source.h5")
volume_path = joinpath(fixture_directory,"python_volume_source.h5")

boundary_hash = file_sha256(boundary_path)
volume_hash = file_sha256(volume_path)
boundary = read_boundary_angular_current_hdf5(
    boundary_path,Photon();expected_file_sha256=boundary_hash,
)
volume = read_anisotropic_volume_source_hdf5(
    volume_path,Photon();expected_file_sha256=volume_hash,
)

boundary_current = get_incoming_current(boundary)
size(boundary_current) == (1,2) || error("Unexpected boundary-current shape.")
isapprox(boundary_current[1,1],8.0;rtol=1.0e-12,atol=1.0e-12) || error(
    "Boundary group-1 current closure failed.",
)
isapprox(boundary_current[1,2],4.0;rtol=1.0e-12,atol=1.0e-12) || error(
    "Boundary group-2 current closure failed.",
)
isapprox(get_total_incoming_current(boundary),12.0;rtol=1.0e-12,atol=1.0e-12) || error(
    "Total boundary-current closure failed.",
)
isapprox(get_total_incoming_current(boundary;physical=true),96.0;rtol=1.0e-12,atol=1.0e-12) || error(
    "Physical boundary-current normalization failed.",
)

volume_rate = get_volume_source_rate(volume)
length(volume_rate) == 2 || error("Unexpected volume-source group count.")
isapprox(volume_rate[1],7.0;rtol=1.0e-12,atol=1.0e-12) || error(
    "Volume group-1 source-rate closure failed.",
)
isapprox(volume_rate[2],10.0;rtol=1.0e-12,atol=1.0e-12) || error(
    "Volume group-2 source-rate closure failed.",
)

receipt = Dict{String,Any}(
    "schema" => "radiant.cross_language_source_replay/v2",
    "status" => "PASS",
    "julia_version" => string(VERSION),
    "storage_order" => "external-row-major",
    "boundary" => Dict{String,Any}(
        "path" => boundary_path,
        "sha256" => boundary_hash,
        "current_per_group" => vec(boundary_current),
        "total_current" => get_total_incoming_current(boundary),
        "physical_total_current_per_s" => get_total_incoming_current(boundary;physical=true),
    ),
    "volume" => Dict{String,Any}(
        "path" => volume_path,
        "sha256" => volume_hash,
        "rate_per_group" => volume_rate,
        "time_class" => string(volume.normalization.time_class),
        "parent_reaction" => volume.parent_reaction,
    ),
    "gates" => Dict{String,Any}(
        "ANALYTIC_CROSSLANGUAGE_REPLAY_PASS" => true,
        "PHYSICAL_OPENMC_REPLAY_PASS" => false,
        "PHYSICAL_OPENSN_REPLAY_PASS" => false,
        "OPENSN_RADIANT_COUPLING_PASS" => false,
    ),
)
mkpath(dirname(receipt_path))
open(receipt_path,"w") do io
    TOML.print(io,receipt)
end
println("Cross-language source replay PASS: $(receipt_path)")
