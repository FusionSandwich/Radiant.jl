#!/usr/bin/env julia

using Pkg
using SHA
using Dates
using TOML

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const RESULT_ROOT = joinpath(ROOT,"qualification-results","julia-$(VERSION)")
const RECEIPT_PATH = joinpath(RESULT_ROOT,"qualification.toml")

function file_sha256(path::AbstractString)
    isfile(path) || return "missing"
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function git_value(arguments...)
    try
        return readchomp(`git -C $ROOT $(arguments...)`)
    catch
        return "unavailable"
    end
end

function julia_commit()
    try
        return string(Base.GIT_VERSION_INFO.commit)
    catch
        return "unavailable"
    end
end

function write_receipt(receipt)
    mkpath(dirname(RECEIPT_PATH))
    open(RECEIPT_PATH,"w") do io
        TOML.print(io,receipt)
    end
end

function run_step!(operation::Function,receipt,name::String)
    started = time()
    receipt["steps"][name] = Dict(
        "status" => "RUNNING","started_at" => string(Dates.now()),
    )
    write_receipt(receipt)
    try
        operation()
        receipt["steps"][name] = Dict(
            "status" => "PASS","duration_s" => time()-started,
            "completed_at" => string(Dates.now()),
        )
        write_receipt(receipt)
        return true
    catch exception
        receipt["steps"][name] = Dict(
            "status" => "FAIL","duration_s" => time()-started,
            "completed_at" => string(Dates.now()),
            "error" => sprint(showerror,exception,catch_backtrace()),
        )
        receipt["overall_status"] = "FAIL"
        write_receipt(receipt)
        rethrow()
    end
end

function run_subqualification(script_name::String)
    script = joinpath(ROOT,"scripts","qualification",script_name)
    run(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $script`)
end

function main()
    cd(ROOT)
    mkpath(RESULT_ROOT)
    receipt = Dict{String,Any}(
        "schema" => "radiant.julia_qualification/v2",
        "overall_status" => "RUNNING",
        "started_at" => string(Dates.now()),
        "repository_root" => ROOT,
        "git_commit" => git_value("rev-parse","HEAD"),
        "git_branch" => git_value("rev-parse","--abbrev-ref","HEAD"),
        "julia_version" => string(VERSION),
        "julia_commit" => julia_commit(),
        "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH),
        "threads" => Threads.nthreads(),
        "project_sha256" => file_sha256(joinpath(ROOT,"Project.toml")),
        "manifest_sha256_before" => file_sha256(joinpath(ROOT,"Manifest.toml")),
        "branch_policy" => "branch-only-no-pr",
        "future_package" => "RadiantHTS.jl",
        "steps" => Dict{String,Any}(),
        "gates" => Dict{String,Any}(
            "RADIANT_BUILD_PASS" => false,
            "HTS_ADDON_SOFTWARE_PASS" => false,
            "SPATIAL_FIELD_SOFTWARE_PASS" => false,
            "FACETED_GEOMETRY_SOFTWARE_PASS" => false,
            "MATERIAL_RESPONSE_REGISTRY_PASS" => false,
            "ATOMISTIC_RESPONSE_PIPELINE_PASS" => false,
            "CLOSED_COUPLING_SOFTWARE_PASS" => false,
            "RADIANT_EM_PASS" => false,
            "PHYSICAL_OPENMC_REPLAY_PASS" => false,
            "PHYSICAL_OPENSN_REPLAY_PASS" => false,
            "OPENSN_RADIANT_COUPLING_PASS" => false,
            "PIECEWISE_FLAT_ATLAS_PASS" => false,
        ),
    )
    write_receipt(receipt)

    run_step!(receipt,"instantiate") do
        Pkg.instantiate()
    end
    receipt["manifest_sha256_after_instantiate"] = file_sha256(joinpath(ROOT,"Manifest.toml"))
    write_receipt(receipt)

    run_step!(receipt,"precompile") do
        Pkg.precompile()
    end
    run_step!(receipt,"package_import") do
        @eval using Radiant
        nothing
    end
    receipt["gates"]["RADIANT_BUILD_PASS"] = true
    write_receipt(receipt)

    run_step!(receipt,"package_tests") do
        Pkg.test(;coverage=false)
    end
    run_step!(receipt,"extended_hts_tests") do
        include(joinpath(ROOT,"test","hts_extended_tests.jl"))
    end
    run_step!(receipt,"hts_addon_tests") do
        include(joinpath(ROOT,"test","hts_addon_tests.jl"))
    end
    run_step!(receipt,"hts_production_foundations_tests") do
        include(joinpath(ROOT,"test","hts_production_foundations_tests.jl"))
    end
    for gate in (
        "HTS_ADDON_SOFTWARE_PASS","SPATIAL_FIELD_SOFTWARE_PASS",
        "FACETED_GEOMETRY_SOFTWARE_PASS","MATERIAL_RESPONSE_REGISTRY_PASS",
        "ATOMISTIC_RESPONSE_PIPELINE_PASS","CLOSED_COUPLING_SOFTWARE_PASS",
    )
        receipt["gates"][gate] = true
    end
    write_receipt(receipt)

    run_step!(receipt,"hdf5_interchange_tests") do
        include(joinpath(ROOT,"test","hdf5_interchange_tests.jl"))
    end
    run_step!(receipt,"radiant_em_analytic") do
        run_subqualification("run_radiant_em_analytic.jl")
    end
    run_step!(receipt,"opensn_radiant_fixture_replay") do
        run_subqualification("run_opensn_radiant_fixture_replay.jl")
    end

    receipt["overall_status"] = "PASS"
    receipt["completed_at"] = string(Dates.now())
    receipt["manifest_sha256_final"] = file_sha256(joinpath(ROOT,"Manifest.toml"))
    write_receipt(receipt)
    println("RADIANT_BUILD_PASS: Julia $(VERSION)")
    println("HTS_ADDON_SOFTWARE_PASS: Julia $(VERSION)")
    println("Qualification receipt: $(RECEIPT_PATH)")
end

try
    main()
catch exception
    println(stderr,"RADIANT_BUILD_PASS not established for Julia $(VERSION).")
    showerror(stderr,exception,catch_backtrace())
    println(stderr)
    exit(1)
end
