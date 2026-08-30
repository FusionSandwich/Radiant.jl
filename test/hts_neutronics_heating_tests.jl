using Radiant
using Test

function _heating_contribution(
    id::String,
    response::Symbol,
    values;
    population::String="prompt-photon-population",
    owner::String="Radiant",
    time_class::Symbol=:prompt,
    source_class::Symbol=:external_prompt,
    additive::Bool=true,
    comparison::Bool=false,
)
    return Layer_Heating_Contribution(
        contribution_id=id,
        population_id=population,
        owner=owner,
        domain_id="tape-microdomain",
        layer_id="REBCO",
        particle_tag=response == :joule_heat ? "electrical" : "photon-electron",
        source_class=source_class,
        time_class=time_class,
        response_channel=response,
        values=reshape(Float64.(values),2,1,1),
        standard_uncertainty=reshape(fill(0.01,2),2,1,1),
        units="W/cm3",
        source_hash="source-$(population)",
        geometry_hash="two-cell-rebco",
        normalization_hash="absolute-rate-v1",
        production_additive=additive,
        comparison_only=comparison,
        metadata=Dict("classification" => "software-verification"),
    )
end

@testset "No-double-counting layer heating ledger" begin
    ledger = Layer_Heating_Ledger(metadata=Dict(
        "classification" => "software-verification",
    ))
    add_heating_contribution!(ledger,_heating_contribution(
        "radiant-total",:total_deposition,[10.0,20.0],
    ))
    add_heating_contribution!(ledger,_heating_contribution(
        "radiant-prompt-heat",:prompt_lattice_heat,[6.0,12.0],
    ))
    add_heating_contribution!(ledger,_heating_contribution(
        "radiant-nonequilibrium",:non_equilibrium,[2.0,4.0],
    ))
    add_heating_contribution!(ledger,_heating_contribution(
        "radiant-defect-storage",:defect_storage,[1.0,2.0],
    ))
    add_heating_contribution!(ledger,_heating_contribution(
        "radiant-escape",:escaped,[1.0,2.0],
    ))

    # OpenSn may provide a comparison for the same photon population, but it cannot be added.
    add_heating_contribution!(ledger,_heating_contribution(
        "opensn-comparison",:total_deposition,[9.8,20.2];
        owner="OpenSn",additive=false,comparison=true,
    ))
    @test validate_heating_ledger(ledger)
    closure = validate_population_energy_closure(ledger,"prompt-photon-population")
    @test closure.closed
    @test closure.maximum_absolute_residual <= 1.0e-14

    prompt = production_heating_total(ledger)
    @test vec(prompt.values) == [6.0,12.0]
    @test prompt.contribution_ids == ["radiant-prompt-heat"]
    comparison = heating_comparison_set(
        ledger,"prompt-photon-population",:total_deposition,
    )
    @test length(comparison) == 2
    @test count(value -> value.comparison_only,comparison) == 1

    # A second additive scorer for the same population/response/layer is rejected.
    @test_throws ErrorException add_heating_contribution!(ledger,_heating_contribution(
        "opensn-illegal-addition",:total_deposition,[9.8,20.2];owner="OpenSn",
    ))

    add_heating_contribution!(ledger,_heating_contribution(
        "activation-delayed",:delayed_lattice_heat,[1.0,1.0];
        population="activation-decay-photons",owner="Radiant",
        time_class=:delayed,source_class=:activation_delayed,
    ))
    add_heating_contribution!(ledger,_heating_contribution(
        "current-sharing-joule",:joule_heat,[0.5,0.5];
        population="electrothermal-current",owner="Electrothermal",
        time_class=:dynamic,source_class=:joule,
    ))
    all_heat = production_heating_total(ledger)
    @test vec(all_heat.values) == [7.5,13.5]
    @test production_heating_total(ledger;time_class=:delayed).values[:] == [1.0,1.0]
    receipt = heating_ledger_receipt(ledger)
    @test receipt["comparison_fields_excluded_from_production_sum"]
    @test receipt["one_additive_owner_per_population_domain_layer_response_time"]
    @test receipt["comparison_only_count"] == 1
end

@testset "Groupwise Gd self-shielding and capture depth" begin
    edges = [0.0,1.0e3,1.0e6]
    gd155 = Gd_Groupwise_Capture_Component(
        "Gd-155",1.0e22,edges,[10.0,1.0];
        data_hash="synthetic-gd155-xs",status=:verification,
    )
    gd157 = Gd_Groupwise_Capture_Component(
        "Gd-157",2.0e22,edges,[20.0,2.0];
        data_hash="synthetic-gd157-xs",status=:verification,
    )
    layer = Gd_Self_Shielding_Layer(
        "gdbco-film-stack","GdBCO",1.0,[gd155,gd157];
        subdivisions=2,background_absorption_cm_inv=[0.1,0.0],
        metadata=Dict("classification" => "analytic-thick-screening-fixture"),
    )
    incident = [100.0,200.0]
    result = solve_gd_groupwise_self_shielding(
        [layer],incident;
        direction_cosine=1.0,geometry_hash="analytic-gdbco-slab",
        transport_artifact_hash="analytic-gdbco-screening-v1",
    )

    @test length(result.cells) == 2
    @test result.transmitted_current_per_s ≈ [
        100.0*exp(-0.6),200.0*exp(-0.05),
    ] rtol=1.0e-13
    removed_group_1 = 100.0*(-expm1(-0.6))
    removed_group_2 = 200.0*(-expm1(-0.05))
    @test result.capture_rate_per_s["Gd-155"] ≈ [
        removed_group_1*(0.1/0.6),removed_group_2*(0.01/0.05),
    ] rtol=1.0e-13
    @test result.capture_rate_per_s["Gd-157"] ≈ [
        removed_group_1*(0.4/0.6),removed_group_2*(0.04/0.05),
    ] rtol=1.0e-13
    @test result.background_absorption_rate_per_s[1] ≈
          removed_group_1*(0.1/0.6) rtol=1.0e-13
    @test maximum(abs.(result.particle_balance_residual_per_s)) <= 1.0e-12
    @test result.cells[2].capture_rate_per_s["Gd-157"][1] <
          result.cells[1].capture_rate_per_s["Gd-157"][1]

    factors_157 = gd_groupwise_self_shielding_factor(result,[layer],"Gd-157")
    @test all(0.0 .< factors_157 .< 1.0)

    oblique = solve_gd_groupwise_self_shielding(
        [layer],incident;
        direction_cosine=0.5,geometry_hash="analytic-gdbco-oblique",
        transport_artifact_hash="analytic-gdbco-oblique-v1",
    )
    @test all(oblique.transmitted_current_per_s .< result.transmitted_current_per_s)
    @test all(oblique.capture_rate_per_s["Gd-157"] .>
              result.capture_rate_per_s["Gd-157"])

    capture_field = capture_rate_field_from_gd_self_shielding(
        result,"Gd-157",[101,102],[0.5,0.5],
    )
    @test capture_field.voxel_ids == [101,102]
    @test sum(capture_field.capture_rate_per_s) ≈
          sum(result.capture_rate_per_s["Gd-157"]) rtol=1.0e-13
    @test capture_field.provenance["capture_rates_are_self_shielded"] == "true"

    receipt = gd_self_shielding_receipt(result)
    @test receipt["particle_balance_pass"]
    @test receipt["screening_only"]
    @test !receipt["physical_transport_qualification"]
    @test receipt["input_data_status"] == "verification"
end
