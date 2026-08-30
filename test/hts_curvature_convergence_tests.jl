using Radiant
using Test

function _curvature_level(
    id::String,
    h::Float64,
    sagitta::Float64,
    heating::Float64;
    source_hash::String="analytic-curved-source",
    hotspot::String="facet-17",
)
    return Curvature_Refinement_Level(
        id,h,sagitta,Dict("rebco_peak_heating" => heating);
        standard_uncertainties=Dict("rebco_peak_heating" => 1.0e-4),
        hotspot_ids=Dict("rebco_peak_heating" => hotspot),
        geometry_hash="geometry-$(id)",
        source_hash=source_hash,
        metadata=Dict("classification" => "software-verification"),
    )
end

@testset "Three-level curvature convergence and physical reference gate" begin
    levels = [
        _curvature_level("coarse",0.4,0.004,1.04),
        _curvature_level("medium",0.2,0.001,1.01),
        _curvature_level("fine",0.1,0.00025,1.0025),
    ]
    software = evaluate_curvature_convergence(
        levels;
        response_relative_tolerance=0.01,
        response_absolute_tolerance=0.0,
        uncertainty_multiplier=2.0,
    )
    @test software.software_pass
    @test !software.physical_reference_present
    @test !software.production_pass
    result = software.responses["rebco_peak_heating"]
    @test result.monotonic
    @test result.hotspot_stable
    @test result.observed_order ≈ 2.0 rtol=1.0e-12
    @test result.richardson_limit ≈ 1.0 rtol=1.0e-12
    @test result.fine_grid_convergence_index > 0.0
    receipt = curvature_convergence_receipt(software)
    @test receipt["software_pass"]
    @test !receipt["production_pass"]

    reference = Physical_Reference_Response(
        response_id="rebco_peak_heating",
        process_key="layer-heating",
        particle_tag="photon",
        units="W/cm3",
        values=reshape([1.0],1),
        standard_uncertainty=reshape([1.0e-3],1),
        classification=:opensn,
        producer="OpenSn-curved-reference",
        source_artifact_hash="analytic-curved-source",
        result_artifact_hash="opensn-curved-result-hash",
        geometry_hash="curved-reference-geometry",
        material_state_hash="analytic-ybco-state",
        normalization_hash="analytic-normalization",
        metadata=Dict("classification" => "test-physical-reference-contract"),
    )
    physical = evaluate_curvature_convergence(
        levels;
        response_relative_tolerance=0.01,
        response_absolute_tolerance=0.0,
        uncertainty_multiplier=2.0,
        physical_references=Dict("rebco_peak_heating" => reference),
    )
    @test physical.software_pass
    @test physical.physical_reference_present
    @test physical.responses["rebco_peak_heating"].physical_pass
    @test physical.production_pass
    @test curvature_convergence_receipt(physical)["production_pass"]

    hotspot_change = [
        _curvature_level("coarse",0.4,0.004,1.04;hotspot="facet-3"),
        _curvature_level("medium",0.2,0.001,1.01;hotspot="facet-9"),
        _curvature_level("fine",0.1,0.00025,1.0025;hotspot="facet-17"),
    ]
    unstable = evaluate_curvature_convergence(
        hotspot_change;
        response_relative_tolerance=0.01,
    )
    @test !unstable.responses["rebco_peak_heating"].hotspot_stable
    @test !unstable.software_pass
    @test !unstable.production_pass

    mismatched_source = copy(levels)
    mismatched_source[3] = _curvature_level(
        "fine",0.1,0.00025,1.0025;source_hash="different-source",
    )
    @test_throws ErrorException evaluate_curvature_convergence(
        mismatched_source;
        response_relative_tolerance=0.01,
    )

    synthetic_reference = Physical_Reference_Response(
        response_id="rebco_peak_heating",
        process_key="layer-heating",
        particle_tag="photon",
        units="W/cm3",
        values=reshape([1.0],1),
        standard_uncertainty=reshape([1.0e-3],1),
        classification=:synthetic,
        producer="synthetic",
        source_artifact_hash="analytic-curved-source",
        result_artifact_hash="synthetic-result",
        geometry_hash="synthetic-geometry",
        material_state_hash="analytic-ybco-state",
        normalization_hash="analytic-normalization",
    )
    @test_throws ErrorException evaluate_curvature_convergence(
        levels;
        response_relative_tolerance=0.01,
        physical_references=Dict("rebco_peak_heating" => synthetic_reference),
    )
end
