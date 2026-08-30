using Radiant
using Test

function _uncertainty_component(
    id::String,
    category::Symbol,
    values;
    correlation_group::String=id,
    classification::Symbol=:analytic,
)
    return Response_Uncertainty_Component(
        component_id=id,
        category=category,
        standard_uncertainty=reshape(Float64.(values),2,1,1),
        correlation_group=correlation_group,
        artifact_hash="$(id)-artifact-hash",
        classification=classification,
        metadata=Dict("classification" => "software-verification"),
    )
end

@testset "Protected-response uncertainty budget" begin
    budget = Protected_Response_Uncertainty_Budget(
        "REBCO-prompt-heating",reshape([10.0,20.0],2,1,1);
        units="W/cm3",
        result_artifact_hash="heating-result-hash",
    )
    add_uncertainty_component!(budget,_uncertainty_component(
        "source-statistics",:source_bank_sampling,[1.0,2.0];
        correlation_group="source",
    ))
    add_uncertainty_component!(budget,_uncertainty_component(
        "source-projection",:source_projection,[0.5,1.0];
        correlation_group="source",
    ))
    add_uncertainty_component!(budget,_uncertainty_component(
        "material-property",:material_property,[2.0,1.0];
        correlation_group="material",
    ))

    combined = vec(combined_standard_uncertainty(budget))
    @test combined[1] ≈ 2.5 atol=1.0e-14
    @test combined[2] ≈ sqrt(10.0) atol=1.0e-14
    relative = vec(relative_combined_uncertainty(budget))
    @test relative[1] ≈ 0.25
    @test relative[2] ≈ sqrt(10.0)/20.0

    source_only = vec(combined_standard_uncertainty(
        budget;categories=(:source_bank_sampling,:source_projection),
    ))
    @test source_only == [1.5,3.0]
    receipt = uncertainty_budget_receipt(budget)
    @test receipt["component_count"] == 3
    @test receipt["correlation_groups"] == ["material","source"]
    @test !receipt["all_components_physical"]

    @test_throws ErrorException add_uncertainty_component!(
        budget,_uncertainty_component("source-statistics",:other,[0.1,0.1]),
    )
    bad_shape = Response_Uncertainty_Component(
        component_id="bad-shape",category=:other,
        standard_uncertainty=zeros(3,1,1),
        artifact_hash="bad-shape-hash",
    )
    @test_throws ErrorException add_uncertainty_component!(budget,bad_shape)
end

function _level(id::String,step::Float64,values;uncertainty=fill(0.0,2))
    return Protected_Response_Convergence_Level(
        level_id=id,
        characteristic_step=step,
        values=reshape(Float64.(values),2,1,1),
        standard_uncertainty=reshape(Float64.(uncertainty),2,1,1),
        artifact_hash="$(id)-result-hash",
    )
end

function _passing_study(axis::Symbol;physical::Bool=false)
    levels = [
        _level("coarse-$(axis)",1.0,[1.20,2.40]),
        _level("medium-$(axis)",0.5,[1.05,2.10]),
        _level("fine-$(axis)",0.25,[1.01,2.02]),
    ]
    return evaluate_response_convergence(
        levels;
        response_id="REBCO-prompt-heating",
        axis=axis,
        relative_tolerance=0.05,
        absolute_tolerance=1.0e-12,
        uncertainty_multiplier=2.0,
        physical_reference=physical,
        source_hash="boundary-source-hash",
        geometry_hash="tape-geometry-hash",
    )
end

@testset "Protected-response convergence" begin
    spatial = _passing_study(:spatial)
    @test spatial.passed
    @test spatial.maximum_finest_relative_change < 0.05
    @test spatial.observed_order !== nothing
    @test spatial.observed_order > 1.0
    receipt = convergence_result_receipt(spatial)
    @test receipt["axis"] == "spatial"
    @test receipt["passed"]
    @test !receipt["physical_reference"]

    statistically_consistent = evaluate_response_convergence(
        [
            _level("coarse",1.0,[1.0,2.0];uncertainty=[0.5,0.5]),
            _level("medium",0.5,[1.2,2.2];uncertainty=[0.5,0.5]),
            _level("fine",0.25,[1.4,2.4];uncertainty=[0.5,0.5]),
        ];
        response_id="statistically-limited-response",
        axis=:source_bank,
        relative_tolerance=0.0,
        absolute_tolerance=0.0,
        uncertainty_multiplier=1.0,
        source_hash="source-bank-hash",
        geometry_hash="geometry-hash",
    )
    @test statistically_consistent.passed

    failing = evaluate_response_convergence(
        [
            _level("coarse-fail",1.0,[1.0,2.0]),
            _level("medium-fail",0.5,[2.0,4.0]),
            _level("fine-fail",0.25,[3.0,6.0]),
        ];
        response_id="diverging-response",
        axis=:angular,
        relative_tolerance=0.01,
        absolute_tolerance=0.0,
        source_hash="source-hash",
        geometry_hash="geometry-hash",
    )
    @test !failing.passed

    @test_throws ErrorException evaluate_response_convergence(
        [
            _level("bad-coarse",1.0,[1.0,2.0]),
            _level("bad-medium",2.0,[1.0,2.0]),
            _level("bad-fine",0.5,[1.0,2.0]),
        ];
        response_id="bad-order",axis=:spatial,
        relative_tolerance=0.1,absolute_tolerance=0.0,
        source_hash="source",geometry_hash="geometry",
    )
end

@testset "Multi-axis response qualification" begin
    spatial = _passing_study(:spatial)
    angular = _passing_study(:angular)
    energy = _passing_study(:energy)
    qualification = qualify_response_across_axes(
        [spatial,angular,energy];
        response_id="REBCO-prompt-heating",
        required_axes=[:spatial,:angular,:energy],
        require_physical_reference=false,
        metadata=Dict("classification" => "software-verification"),
    )
    @test qualification.software_passed
    @test !qualification.physical_passed
    @test isempty(qualification.missing_axes)
    receipt = multi_axis_qualification_receipt(qualification)
    @test receipt["software_passed"]
    @test !receipt["physical_passed"]

    incomplete = qualify_response_across_axes(
        [spatial,energy];
        response_id="REBCO-prompt-heating",
        required_axes=[:spatial,:angular,:energy],
    )
    @test !incomplete.software_passed
    @test incomplete.missing_axes == [:angular]

    @test_throws ErrorException qualify_response_across_axes(
        [spatial,angular,energy];
        response_id="REBCO-prompt-heating",
        required_axes=[:spatial,:angular,:energy],
        require_physical_reference=true,
    )

    physical_qualification = qualify_response_across_axes(
        [
            _passing_study(:spatial;physical=true),
            _passing_study(:angular;physical=true),
            _passing_study(:energy;physical=true),
        ];
        response_id="REBCO-prompt-heating",
        required_axes=[:spatial,:angular,:energy],
        require_physical_reference=true,
    )
    @test physical_qualification.physical_passed
end
