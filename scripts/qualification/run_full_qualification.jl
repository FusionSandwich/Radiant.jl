#!/usr/bin/env julia

using Dates
using Pkg
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT_DIR = get(
    ENV,"RADIANT_QUALIFICATION_OUTPUT",joinpath(ROOT,"artifacts","qualification"),
)
mkpath(OUTPUT_DIR)

function file_sha256(path::AbstractString)
    isfile(path) || error("Qualification artifact does not exist: $(path).")
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

# Julia do-block syntax passes the anonymous function first.
function record_step!(
    operation::Function,
    receipt::Dict{String,Any},
    name::String,
)
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
        step["elapsed_s"] = time()-started
        step["finished_utc"] = string(Dates.now(Dates.UTC))
    end
    return nothing
end

function verify_hash_bound_receipt(
    path_variable::String,
    hash_variable::String,
)
    path = get(ENV,path_variable,"")
    isempty(path) && return nothing
    expected = get(ENV,hash_variable,"")
    isempty(expected) && error("$(path_variable) requires $(hash_variable).")
    artifact_path = abspath(path)
    lowercase(file_sha256(artifact_path)) == lowercase(expected) || error(
        "Hash mismatch for $(path_variable).",
    )
    return artifact_path,TOML.parsefile(artifact_path)
end

receipt = Dict{String,Any}(
    "schema" => "radiant.qualification_receipt/v2",
    "started_utc" => string(Dates.now(Dates.UTC)),
    "julia_version" => string(VERSION),
    "julia_commit" => try string(Base.GIT_VERSION_INFO.commit) catch; "unavailable" end,
    "platform" => string(Sys.MACHINE),
    "threads" => Threads.nthreads(),
    "project_sha256" => file_sha256(joinpath(ROOT,"Project.toml")),
    "branch_policy" => "branch-only-no-pr",
    "future_package" => "RadiantHTS.jl",
    "steps" => Dict{String,Any}(),
    "gates" => Dict{String,Any}(
        "RADIANT_BUILD_PASS" => false,
        "HTS_ADDON_SOFTWARE_PASS" => false,
        "ANALYTIC_HDF5_REPLAY_PASS" => false,
        "OPENMC_SOURCE_ARTIFACT_VALID_PASS" => false,
        "OPENSN_SOURCE_ARTIFACT_VALID_PASS" => false,
        "RADIANT_EM_PASS" => false,
        "PHYSICAL_OPENMC_REPLAY_PASS" => false,
        "PHYSICAL_OPENSN_REPLAY_PASS" => false,
        "OPENSN_RADIANT_COUPLING_PASS" => false,
        "PIECEWISE_FLAT_ATLAS_PASS" => false,
    ),
)

receipt_path = joinpath(
    OUTPUT_DIR,"qualification-julia-$(VERSION.major).$(VERSION.minor).toml",
)
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

        for (step_name,test_name) in (
            ("hts_extended_tests","hts_extended_tests.jl"),
            ("hts_addon_tests","hts_addon_tests.jl"),
            ("hdf5_interchange_tests","hdf5_interchange_tests.jl"),
            ("source_hdf5_tests","source_hdf5_tests.jl"),
        )
            test_path = joinpath(ROOT,"test",test_name)
            if isfile(test_path)
                record_step!(receipt,step_name) do
                    Base.include(Main,test_path)
                end
            end
        end
        receipt["gates"]["HTS_ADDON_SOFTWARE_PASS"] = true
        receipt["gates"]["ANALYTIC_HDF5_REPLAY_PASS"] = true

        if haskey(ENV,"RADIANT_OPENMC_BOUNDARY_HDF5")
            haskey(ENV,"RADIANT_OPENMC_BOUNDARY_SHA256") || error(
                "RADIANT_OPENMC_BOUNDARY_HDF5 requires RADIANT_OPENMC_BOUNDARY_SHA256.",
            )
            artifact_path = abspath(ENV["RADIANT_OPENMC_BOUNDARY_HDF5"])
            expected_hash = ENV["RADIANT_OPENMC_BOUNDARY_SHA256"]
            record_step!(receipt,"openmc_boundary_source_validation") do
                source = Radiant.read_boundary_angular_current_hdf5(
                    artifact_path,Radiant.Photon();
                    expected_file_sha256=expected_hash,
                )
                receipt["openmc_boundary_source"] = Dict{String,Any}(
                    "path" => artifact_path,
                    "sha256" => file_sha256(artifact_path),
                    "particle" => Radiant.get_tag(Radiant.get_particle(source)),
                    "source_type" => string(typeof(source)),
                    "total_incoming_current" => Radiant.get_total_incoming_current(source),
                )
            end
            receipt["gates"]["OPENMC_SOURCE_ARTIFACT_VALID_PASS"] = true
        end

        if haskey(ENV,"RADIANT_OPENSN_VOLUME_HDF5")
            haskey(ENV,"RADIANT_OPENSN_VOLUME_SHA256") || error(
                "RADIANT_OPENSN_VOLUME_HDF5 requires RADIANT_OPENSN_VOLUME_SHA256.",
            )
            artifact_path = abspath(ENV["RADIANT_OPENSN_VOLUME_HDF5"])
            expected_hash = ENV["RADIANT_OPENSN_VOLUME_SHA256"]
            record_step!(receipt,"opensn_volume_source_validation") do
                source = Radiant.read_anisotropic_volume_source_hdf5(
                    artifact_path,Radiant.Photon();
                    expected_file_sha256=expected_hash,
                )
                receipt["opensn_volume_source"] = Dict{String,Any}(
                    "path" => artifact_path,
                    "sha256" => file_sha256(artifact_path),
                    "particle" => Radiant.get_tag(Radiant.get_particle(source)),
                    "source_type" => string(typeof(source)),
                    "integrated_rate" => Radiant.get_volume_source_rate(source),
                )
            end
            receipt["gates"]["OPENSN_SOURCE_ARTIFACT_VALID_PASS"] = true
        end

        em_evidence = verify_hash_bound_receipt(
            "RADIANT_EM_QUALIFICATION_RECEIPT",
            "RADIANT_EM_QUALIFICATION_SHA256",
        )
        if !isnothing(em_evidence)
            path,parsed = em_evidence
            get(parsed,"schema","") == "radiant.em_qualification/v1" || error(
                "Unsupported Radiant EM qualification receipt schema.",
            )
            get(parsed,"all_required_benchmarks_pass",false) || error(
                "Radiant EM receipt does not pass all required benchmarks.",
            )
            get(parsed,"physical_reference_data",false) || error(
                "Radiant EM receipt lacks physical reference data.",
            )
            get(parsed,"clipping_used",true) && error(
                "Radiant EM qualification cannot use clipping.",
            )
            receipt["radiant_em_receipt"] = Dict(
                "path" => path,"sha256" => file_sha256(path),
            )
            receipt["gates"]["RADIANT_EM_PASS"] = true
        end

        for (prefix,gate) in (
            ("RADIANT_OPENMC_REPLAY","PHYSICAL_OPENMC_REPLAY_PASS"),
            ("RADIANT_OPENSN_REPLAY","PHYSICAL_OPENSN_REPLAY_PASS"),
        )
            evidence = verify_hash_bound_receipt(
                string(prefix,"_RECEIPT"),string(prefix,"_SHA256"),
            )
            isnothing(evidence) && continue
            path,parsed = evidence
            get(parsed,"schema","") == "radiant.physical_replay/v1" || error(
                "Unsupported physical replay receipt schema for $(prefix).",
            )
            get(parsed,"response_convergence_pass",false) || error(
                "Physical replay response convergence failed for $(prefix).",
            )
            get(parsed,"particle_balance_pass",false) || error(
                "Physical replay particle balance failed for $(prefix).",
            )
            get(parsed,"energy_balance_pass",false) || error(
                "Physical replay energy balance failed for $(prefix).",
            )
            receipt[string(lowercase(prefix),"_receipt")] = Dict(
                "path" => path,"sha256" => file_sha256(path),
            )
            receipt["gates"][gate] = true
        end

        coupling_evidence = verify_hash_bound_receipt(
            "RADIANT_OPENSN_COUPLING_RECEIPT",
            "RADIANT_OPENSN_COUPLING_SHA256",
        )
        if !isnothing(coupling_evidence)
            path,parsed = coupling_evidence
            get(parsed,"schema","") == "opensn-radiant.coupling/v1" || error(
                "Unsupported OpenSn-Radiant coupling receipt schema.",
            )
            get(parsed,"forward_current_closure_pass",false) || error(
                "OpenSn-to-Radiant current closure failed.",
            )
            get(parsed,"return_current_closure_pass",false) || error(
                "Radiant-to-OpenSn return-current closure failed.",
            )
            get(parsed,"energy_current_closure_pass",false) || error(
                "OpenSn-Radiant energy-current closure failed.",
            )
            get(parsed,"response_convergence_pass",false) || error(
                "OpenSn-Radiant response convergence failed.",
            )
            get(parsed,"clipping_used",true) && error(
                "OpenSn-Radiant coupling cannot use clipping.",
            )
            receipt["opensn_radiant_coupling_receipt"] = Dict(
                "path" => path,"sha256" => file_sha256(path),
            )
            receipt["gates"]["OPENSN_RADIANT_COUPLING_PASS"] = true
        end

        atlas_evidence = verify_hash_bound_receipt(
            "RADIANT_CURVED_ATLAS_RECEIPT",
            "RADIANT_CURVED_ATLAS_SHA256",
        )
        if !isnothing(atlas_evidence)
            path,parsed = atlas_evidence
            get(parsed,"schema","") == "radiant.hts.curved_atlas_qualification/v1" || error(
                "Unsupported curved-atlas qualification receipt schema.",
            )
            get(parsed,"flat_curved_response_convergence_pass",false) || error(
                "Piecewise-flat versus curved response convergence failed.",
            )
            get(parsed,"source_mapping_closure_pass",false) || error(
                "Atlas source mapping closure failed.",
            )
            receipt["curved_atlas_receipt"] = Dict(
                "path" => path,"sha256" => file_sha256(path),
            )
            receipt["gates"]["PIECEWISE_FLAT_ATLAS_PASS"] = true
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
