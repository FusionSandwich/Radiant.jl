using Radiant
using Test

function _reference_response(
    response_id::String,
    process_key::String,
    value::Float64;
    classification::Symbol=:opensn,
    result_hash::String="result-hash",
)
    return Physical_Reference_Response(
        response_id=response_id,
        process_key=process_key,
        particle_tag="photon",
        units="MeV/g per transport source basis",
        values=reshape([value],1,1,1),
        standard_uncertainty=reshape([0.01],1,1,1),
        classification=classification,
        producer=classification == :opensn ? "OpenSn" : "OpenMC",
        source_artifact_hash="matched-source-hash",
        result_artifact_hash=result_hash,
        geometry_hash="matched-geometry-hash",
        material_state_hash="matched-material-state",
        normalization_hash="matched-normalization-hash",
        metadata=Dict(
            "scoring_semantics" => "layer-process-energy-deposition",
            "classification" => string(classification),
        ),
    )
end

@testset "Hash-bound physical reference bundle" begin
    directory = mktempdir()
    path = joinpath(directory,"physical-reference.h5")
    references = [
        _reference_response(
            "photoelectric-response","photoelectric",1.0;
            result_hash="opensn-photoelectric-result",
        ),
        _reference_response(
            "compton-response","compton",2.0;
            result_hash="opensn-compton-result",
        ),
    ]
    digest = write_physical_reference_bundle_hdf5(
        path,references;
        metadata=Dict(
            "case_id" => "matched-two-process-fixture",
            "physical_values_are_synthetic_for_test" => "true",
        ),
    )
    @test length(digest) == 64
    replayed = read_physical_reference_bundle_hdf5(
        path;expected_file_sha256=digest,
    )
    @test replayed.artifact_sha256 == digest
    @test Set(keys(replayed.responses)) ==
          Set(["photoelectric-response","compton-response"])
    @test replayed.responses["photoelectric-response"].values ==
          references[1].values
    @test replayed.responses["photoelectric-response"].standard_uncertainty ==
          references[1].standard_uncertainty
    @test replayed.responses["photoelectric-response"].classification == :opensn
    @test references_by_process_key(replayed)["compton"].values[1] == 2.0
    @test references_by_response_id(replayed)["photoelectric-response"].process_key ==
          "photoelectric"

    receipt = physical_reference_bundle_receipt(replayed)
    @test receipt["artifact_sha256"] == digest
    @test receipt["response_count"] == 2
    @test receipt["physical_reference_only"]
    @test receipt["classifications"] == ["opensn"]

    score = Process_Resolved_Score(
        "photon",
        "energy-deposition",
        Dict(
            "photoelectric" => reshape([1.005],1,1,1),
            "compton" => reshape([1.99],1,1,1),
        ),
        reshape([2.995],1,1,1);
        units="MeV/g per transport source basis",
        provenance=Dict("source_hash" => "matched-source-hash"),
    )
    qualification = qualify_process_score_from_reference_bundle(
        score,replayed;
        candidate_standard_uncertainty=Dict(
            "photoelectric" => reshape([0.005],1,1,1),
            "compton" => reshape([0.005],1,1,1),
        ),
        required_process_keys=["photoelectric","compton"],
        relative_tolerance=0.01,
        absolute_tolerance=0.0,
        uncertainty_multiplier=2.0,
        physical_reference_only=true,
    )
    @test qualification.passed
    @test isempty(qualification.missing_process_keys)

    @test_throws ErrorException read_physical_reference_bundle_hdf5(
        path;expected_file_sha256=repeat("0",64),
    )

    duplicate_process = Physical_Reference_Bundle([
        _reference_response("a","same",1.0;result_hash="a"),
        _reference_response("b","same",1.0;result_hash="b"),
    ])
    @test_throws ErrorException references_by_process_key(duplicate_process)

    synthetic_bundle = Physical_Reference_Bundle([
        _reference_response(
            "synthetic","synthetic-process",1.0;
            classification=:synthetic,
            result_hash="synthetic-result",
        ),
    ])
    @test !physical_reference_bundle_receipt(synthetic_bundle)["physical_reference_only"]
end
