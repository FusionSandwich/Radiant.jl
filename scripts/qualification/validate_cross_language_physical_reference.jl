#!/usr/bin/env julia

using Radiant
using SHA
using TOML

function file_sha256(path::AbstractString)
    isfile(path) || error("Physical reference fixture is missing: $(path).")
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function main(bundle_path::AbstractString,receipt_path::AbstractString)
    bundle_path = abspath(bundle_path)
    receipt_path = abspath(receipt_path)
    digest = file_sha256(bundle_path)
    bundle = read_physical_reference_bundle_hdf5(
        bundle_path;expected_file_sha256=digest,
    )

    expected_photoelectric = [1.0 2.0; 3.0 4.0]
    expected_photoelectric_uncertainty = [0.01 0.02; 0.03 0.04]
    expected_compton = [0.5 0.25; 0.125 0.0625]
    expected_compton_uncertainty = fill(0.005,2,2)

    photoelectric = bundle.responses["photoelectric-layer-field"]
    compton = bundle.responses["compton-layer-field"]
    photoelectric.values == expected_photoelectric || error(
        "Python-to-Julia photoelectric tensor ordering changed.",
    )
    photoelectric.standard_uncertainty == expected_photoelectric_uncertainty || error(
        "Python-to-Julia photoelectric uncertainty ordering changed.",
    )
    compton.values ≈ expected_compton || error(
        "Python-to-Julia Compton tensor ordering changed.",
    )
    compton.standard_uncertainty ≈ expected_compton_uncertainty || error(
        "Python-to-Julia Compton uncertainty ordering changed.",
    )
    all(response.classification == :synthetic for response in Base.values(bundle.responses)) ||
        error("Cross-language fixture must remain classified as synthetic.")
    all(response.source_artifact_hash == repeat("1",64)
        for response in Base.values(bundle.responses)) || error(
        "Source lineage was not preserved by the Python producer.",
    )

    raw_reference_by_process = references_by_process_key(bundle)
    reference_by_process = Dict{String,Physical_Reference_Response}()
    for (key,response) in raw_reference_by_process
        reference_by_process[key] = Physical_Reference_Response(
            response_id=response.response_id,
            process_key=response.process_key,
            particle_tag=response.particle_tag,
            units=response.units,
            values=reshape(response.values,2,2,1),
            standard_uncertainty=reshape(response.standard_uncertainty,2,2,1),
            classification=response.classification,
            producer=response.producer,
            source_artifact_hash=response.source_artifact_hash,
            result_artifact_hash=response.result_artifact_hash,
            geometry_hash=response.geometry_hash,
            material_state_hash=response.material_state_hash,
            normalization_hash=response.normalization_hash,
            metadata=response.metadata,
        )
    end

    score_channels = Dict(
        "photoelectric" => reshape(expected_photoelectric,2,2,1),
        "compton" => reshape(expected_compton,2,2,1),
    )
    score_total = score_channels["photoelectric"]+score_channels["compton"]
    score = Process_Resolved_Score(
        "photon","energy-deposition",score_channels,score_total;
        units="MeV/g per transport source basis",
        provenance=Dict(
            "classification" => "software-verification",
            "source_artifact_hash" => repeat("1",64),
        ),
    )
    qualification = qualify_process_resolved_score(
        score,reference_by_process;
        required_process_keys=["photoelectric","compton"],
        relative_tolerance=0.0,
        absolute_tolerance=1.0e-14,
        uncertainty_multiplier=0.0,
        physical_reference_only=false,
    )
    qualification.passed || error("Synthetic cross-language response replay did not close.")

    physical_only_rejected = false
    try
        qualify_process_resolved_score(
            score,reference_by_process;
            required_process_keys=["photoelectric","compton"],
            relative_tolerance=0.0,
            absolute_tolerance=1.0e-14,
            uncertainty_multiplier=0.0,
            physical_reference_only=true,
        )
    catch
        physical_only_rejected = true
    end
    physical_only_rejected || error(
        "Synthetic reference bundle was incorrectly accepted as physical qualification.",
    )

    bundle_receipt = physical_reference_bundle_receipt(bundle)
    bundle_receipt["physical_reference_only"] == false || error(
        "Synthetic fixture receipt incorrectly reports physical references.",
    )

    receipt = Dict{String,Any}(
        "schema" => "radiant.hts.cross_language_physical_reference/v1",
        "status" => "ANALYTIC_CROSSLANGUAGE_REFERENCE_PASS",
        "classification" => "software-verification",
        "julia_version" => string(VERSION),
        "bundle_path" => bundle_path,
        "bundle_sha256" => digest,
        "response_ids" => sort(collect(keys(bundle.responses))),
        "tensor_shape" => [2,2],
        "tensor_order_preserved" => true,
        "uncertainty_order_preserved" => true,
        "lineage_preserved" => true,
        "synthetic_process_qualification_pass" => qualification.passed,
        "physical_only_rejection_pass" => physical_only_rejected,
        "physical_reference_pass" => false,
        "physical_openmc_replay_pass" => false,
        "physical_opensn_replay_pass" => false,
    )
    mkpath(dirname(receipt_path))
    open(receipt_path,"w") do io
        TOML.print(io,receipt)
    end
    println("ANALYTIC_CROSSLANGUAGE_REFERENCE_PASS")
    println("Reference replay receipt: $(receipt_path)")
end

length(ARGS) == 2 || error(
    "Usage: validate_cross_language_physical_reference.jl <bundle.h5> <receipt.toml>",
)
main(ARGS[1],ARGS[2])
