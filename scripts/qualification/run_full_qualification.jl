#!/usr/bin/env julia

using Dates
using Pkg
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT_DIR = get(ENV,"RADIANT_QUALIFICATION_OUTPUT",joinpath(ROOT,"artifacts","qualification"))
mkpath(OUTPUT_DIR)

function file_sha256(path::AbstractString)
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function record_step!(receipt::Dict{String,Any},name::String,operation::Function)
    started = time()
    step = Dict{String,Any}(
        "started_utc" => string(Dates.now(Dates.UTC)),
        "status" => "RUNNING",
    )
    receipt["steps"][name] = step
    try
        operation()
        step["status"] = "PASS"
    catch error_value
        step["status"] = "FAIL"
        step["error_type"] = string(typeof(error_value))
        step["error"] = sprint(showerror,error_value,catch_backtrace())
        rethrow()
    finally
        step["elapsed_s"] = time() - started
        step["finished_utc"] = string(Dates.now(Dates.UTC))
    end
    return nothing
end

receipt = Dict{String,Any}(
    "schema" => "radiant.qualification_receipt/v1",
    "started_utc" => string(Dates.now(Dates.UTC)),
    "julia_version" => string(VERSION),
    "julia_commit" => string(Base.GIT_VERSION_INFO.commit),
    "platform" => string(Sys.MACHINE),
    "threads" => Threads.nthreads(),
    "project_sha256" => file_sha256(joinpath(ROOT,"Project.toml")),
    "steps" => Dict{String,Any}(),
    "gates" => Dict{String,Any}(
        "RADIANT_BUILD_PASS" => false,
        "RADIANT_EM_PASS" => false,
        "ANALYTIC_HDF5_REPLAY_PASS" => false,
        "PHYSICAL_OPENMC_REPLAY_PASS" => false,
        "PHYSICAL_OPENSN_REPLAY_PASS" => false,
        "OPENSN_RADIANT_COUPLING_PASS" => false,
    ),
)

receipt_path = joinpath(OUTPUT_DIR,"qualification-julia-$(VERSION.major).$(VERSION.minor).toml")
exit_code = 0

try
    cd(ROOT) do
        record_step!(receipt,"instantiate") do
            Pkg.activate(ROOT)
            Pkg.instantiate()
        end
        record_step!(receipt,"precompile") do
            Pkg.precompile()
        end
        record_step!(receipt,"package_import") do
            @eval using Radiant
        end
        receipt["gates"]["RADIANT_BUILD_PASS"] = true

        record_step!(receipt,"package_tests") do
            Pkg.test(;coverage=false)
        end
        receipt["gates"]["RADIANT_EM_PASS"] = true

        extended_test = joinpath(ROOT,"test","hts_extended_tests.jl")
        if isfile(extended_test)
            record_step!(receipt,"hts_extended_tests") do
                Base.include(Main,extended_test)
            end
        end

        hdf5_test = joinpath(ROOT,"test","source_hdf5_tests.jl")
        record_step!(receipt,"analytic_hdf5_replay") do
            isfile(hdf5_test) || error("Missing HDF5 source replay test: $(hdf5_test).")
            Base.include(Main,hdf5_test)
        end
        receipt["gates"]["ANALYTIC_HDF5_REPLAY_PASS"] = true

        # Physical artifacts are deliberately opt-in and hash-bound. A path alone is not enough:
        # the corresponding expected SHA-256 environment variable is required.
        physical_inputs = [
            ("openmc", "RADIANT_OPENMC_BOUNDARY_HDF5", "RADIANT_OPENMC_BOUNDARY_SHA256"),
            ("opensn", "RADIANT_OPENSN_VOLUME_HDF5", "RADIANT_OPENSN_VOLUME_SHA256"),
        ]
        for (producer,path_variable,hash_variable) in physical_inputs
            if haskey(ENV,path_variable)
                haskey(ENV,hash_variable) || error(
                    "$(path_variable) was supplied without required $(hash_variable).",
                )
                artifact_path = abspath(ENV[path_variable])
                expected_hash = ENV[hash_variable]
                record_step!(receipt,"physical_$(producer)_source_validation") do
                    source = Radiant.read_radiant_source_hdf5(
                        artifact_path;expected_sha256=expected_hash,
                    )
                    source_particle = Radiant.get_particle(source)
                    receipt["physical_$(producer)_source"] = Dict{String,Any}(
                        "path" => artifact_path,
                        "sha256" => Radiant.source_artifact_sha256(artifact_path),
                        "particle" => Radiant.get_tag(source_particle),
                        "source_type" => string(typeof(source)),
                    )
                end
                receipt["gates"][producer == "openmc" ?
                    "PHYSICAL_OPENMC_REPLAY_PASS" : "PHYSICAL_OPENSN_REPLAY_PASS"] = true
            end
        end

        # Full OpenSn-Radiant coupling requires both a physical OpenSn-produced source and a
        # transport/interface closure receipt. Source-file validation alone never establishes it.
        coupling_receipt = get(ENV,"RADIANT_OPENSN_COUPLING_RECEIPT","")
        if !isempty(coupling_receipt)
            isfile(coupling_receipt) || error("OpenSn coupling receipt does not exist.")
            expected = get(ENV,"RADIANT_OPENSN_COUPLING_RECEIPT_SHA256","")
            isempty(expected) && error("OpenSn coupling receipt SHA-256 is required.")
            lowercase(file_sha256(coupling_receipt)) == lowercase(expected) || error(
                "OpenSn coupling receipt hash mismatch.",
            )
            parsed = TOML.parsefile(coupling_receipt)
            get(parsed,"interface_current_closure_pass",false) || error(
                "OpenSn coupling receipt does not establish interface-current closure.",
            )
            get(parsed,"energy_balance_pass",false) || error(
                "OpenSn coupling receipt does not establish energy balance.",
            )
            receipt["gates"]["OPENSN_RADIANT_COUPLING_PASS"] = true
        end
    end
catch
    exit_code = 1
finally
    receipt["finished_utc"] = string(Dates.now(Dates.UTC))
    receipt["overall_status"] = exit_code == 0 ? "PASS" : "FAIL"
    open(receipt_path,"w") do io
        TOML.print(io,receipt;sorted=true)
    end
    println("Qualification receipt: $(receipt_path)")
    println("Overall status: $(receipt[\"overall_status\"])")
end

exit(exit_code)
