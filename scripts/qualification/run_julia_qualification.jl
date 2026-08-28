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

function run_step!(receipt,name::String,operation::Function)
    started = time()
    try
        operation()
        receipt["steps"][name] = Dict(
            "status" => "PASS",
            "duration_s" => time()-started,
        )
        write_receipt(receipt)
        return true
    catch exception
        message = sprint(showerror,exception,catch_backtrace())
        receipt["steps"][name] = Dict(
            "status" => "FAIL",
            "duration_s" => time()-started,
            "error" => message,
        )
        receipt["overall_status"] = "FAIL"
        write_receipt(receipt)
        rethrow()
    end
end

function main()
    cd(ROOT)
    mkpath(RESULT_ROOT)
    receipt = Dict{String,Any}(
        "schema" => "radiant.julia_qualification/v1",
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
        "steps" => Dict{String,Any}(),
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

    run_step!(receipt,"package_tests") do
        Pkg.test(;coverage=false)
    end

    run_step!(receipt,"extended_hts_tests") do
        include(joinpath(ROOT,"test","hts_extended_tests.jl"))
    end

    run_step!(receipt,"opensn_radiant_fixture_replay") do
        replay_script = joinpath(
            ROOT,"scripts","qualification","run_opensn_radiant_fixture_replay.jl",
        )
        run(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $replay_script`)
    end

    receipt["overall_status"] = "PASS"
    receipt["completed_at"] = string(Dates.now())
    receipt["manifest_sha256_final"] = file_sha256(joinpath(ROOT,"Manifest.toml"))
    write_receipt(receipt)
    println("RADIANT_BUILD_PASS: Julia $(VERSION)")
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
