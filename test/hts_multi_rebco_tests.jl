using Radiant
using Test

function family_packet(tag,symbol,status;packet_id=lowercase(tag),product_id=tag,
                       construction_hash="matched-construction",design_hash="matched-design",
                       packet_status=:candidate)
    return REBCO_Material_Family_Packet(
        packet_id,tag;rare_earth_symbol=symbol,identity_status=status,
        identity_evidence_hash="identity-$(packet_id)",composition_hash="composition-$(packet_id)",
        product_id=product_id,product_construction_hash=construction_hash,
        controlled_design_hash=design_hash,status=packet_status,
    )
end

@testset "Versioned multi-REBCO family packets" begin
    packets = Dict(
        "YBCO" => family_packet("YBCO","Y",:confirmed),
        "GdBCO" => family_packet("GdBCO","Gd",:confirmed),
        "EuBCO" => family_packet("EuBCO","Eu",:confirmed),
        "SaBCO" => family_packet("SaBCO","",:unresolved),
        "SmBCO" => family_packet("SmBCO","Sm",:confirmed),
    )
    @test all(validate_rebco_material_packet(packet) for packet in values(packets))
    @test all(!rebco_packet_is_physically_qualified(packet) for packet in values(packets))
    @test packets["SaBCO"].identity_status == :unresolved
    @test isempty(packets["SaBCO"].rare_earth_symbol)
    @test packets["SmBCO"].identity_status == :confirmed
    @test_throws ErrorException family_packet("SaBCO","Sm",:confirmed)
    @test_throws ErrorException family_packet(
        "SaBCO","",:unresolved;packet_status=:qualified_input,
    )

    adapted = adapt_rebco_material_packet(Dict(
        "packet_id" => "eu-adapted","material_tag" => "EuBCO",
        "rare_earth_symbol" => "Eu","identity_status" => "confirmed",
        "identity_evidence_hash" => "eu-identity","composition_hash" => "eu-composition",
        "product_id" => "eu-product","product_construction_hash" => "eu-construction",
        "controlled_design_hash" => "eu-design","status" => "candidate",
    ))
    @test adapted.material_tag == "EuBCO"
    @test adapted.schema == REBCO_MATERIAL_FAMILY_PACKET_SCHEMA
end

@testset "Controlled substitution and realistic products stay comparison-only" begin
    y = family_packet("YBCO","Y",:confirmed;packet_id="y-control",product_id="control-y")
    gd = family_packet("GdBCO","Gd",:confirmed;packet_id="gd-control",product_id="control-gd")
    controlled = REBCO_Comparison_Contract(
        "y-vs-gd-controlled",y,gd;mode=:controlled_substitution,
        response_ids=["deposition","PKA"],paired_axes=["thickness","turn","incidence"],
    )
    @test controlled.comparison_only
    @test !controlled.physical_qualification

    eu = family_packet(
        "EuBCO","Eu",:confirmed;packet_id="eu-product",product_id="vendor-eu",
        construction_hash="vendor-eu-build",design_hash="vendor-eu-design",
    )
    sm = family_packet(
        "SmBCO","Sm",:confirmed;packet_id="sm-product",product_id="vendor-sm",
        construction_hash="vendor-sm-build",design_hash="vendor-sm-design",
    )
    realistic = REBCO_Comparison_Contract(
        "eu-vs-sm-products",eu,sm;mode=:realistic_product,
        response_ids=["activation"],paired_axes=["source","temperature"],
    )
    @test "product differences cannot be attributed solely to rare-earth identity" in
        realistic.limitations
    @test_throws ErrorException REBCO_Comparison_Contract(
        "invalid-controlled",y,eu;mode=:controlled_substitution,
        response_ids=["deposition"],paired_axes=["thickness"],
    )
end

@testset "Rare-earth generic path and isotope preservation" begin
    gd = family_packet("GdBCO","Gd",:confirmed;packet_id="gd-shield")
    gd155 = REBCO_Isotope_Capture_Component(
        gd,"Gd-155",1.0e20,[0.0,1.0,2.0],[10.0,20.0];
        data_hash="gd155-verification",status=:verification,
    )
    gd157 = REBCO_Isotope_Capture_Component(
        gd,"Gd-157",2.0e20,[0.0,1.0,2.0],[30.0,40.0];
        data_hash="gd157-verification",status=:verification,
    )
    result = solve_rebco_self_shielding_axis_sweep(
        gd,[gd155,gd157],[100.0,50.0];thickness_axis_cm=[1.0e-4,2.0e-4],
        turn_axis=[1,9],incidence_axis=[[0.0,0.0,-1.0],[0.6,0.0,-0.8]],
        layer_normal=[0.0,0.0,2.0],geometry_hash="geometry-fixture",
        transport_artifact_hash="analytic-fixture",
    )
    @test length(result.cases) == 8
    @test result.cases[1].path_length_cm ≈ 1.0e-4
    @test result.cases[2].path_length_cm ≈ 1.25e-4
    @test result.cases[1].turn_index == 1
    @test result.cases[5].turn_index == 1
    @test sort(collect(keys(result.cases[1].capture_rate_per_s))) == ["Gd-155","Gd-157"]
    @test maximum(maximum(abs.(case.particle_balance_residual_per_s)) for case in result.cases) <
        1.0e-12
    receipt = rebco_self_shielding_receipt(result)
    @test receipt["screening_only"]
    @test !receipt["physical_transport_qualification"]
    @test receipt["thickness_axis_cm"] == [1.0e-4,2.0e-4]
    @test receipt["turn_axis"] == [1,9]
    @test !rebco_self_shielding_is_physically_qualified(result)
    @test receipt["material_packet_status"] == "candidate"
    @test receipt["isotope_data_status"] == "verification"
    @test_throws ErrorException REBCO_Isotope_Capture_Component(
        gd,"Sm-149",1.0,[0.0,1.0],[1.0];data_hash="wrong-family",
    )
    @test_throws ErrorException solve_rebco_self_shielding_axis_sweep(
        gd,[gd155],[1.0,1.0];thickness_axis_cm=[1.0e-4],turn_axis=[1],
        incidence_axis=[[1.0,0.0,0.0]],layer_normal=[0.0,0.0,1.0],
        geometry_hash="grazing-fixture",
    )

    sm = family_packet("SmBCO","Sm",:confirmed;packet_id="sm-shield")
    sm149 = REBCO_Isotope_Capture_Component(
        sm,"Sm-149",1.0e20,[0.0,1.0],[5.0];
        data_hash="sm149-verification",status=:verification,
    )
    sm_result = solve_rebco_self_shielding_axis_sweep(
        sm,[sm149],[10.0];thickness_axis_cm=[3.0e-5],turn_axis=[4],
        incidence_axis=[[0.0,0.0,1.0]],layer_normal=[0.0,0.0,1.0],
        geometry_hash="sm-geometry",
    )
    @test collect(keys(sm_result.cases[1].capture_rate_per_s)) == ["Sm-149"]
    @test !sm_result.physical_transport_qualification
end

@testset "Family-specific consumer bindings" begin
    sm = family_packet("SmBCO","Sm",:confirmed;packet_id="sm-consumers")
    for consumer_class in (:source,:deposition,:pka,:activation)
        binding = bind_rebco_family_consumer(
            sm,consumer_class;binding_id="sm-$(consumer_class)",
            producer_artifact_hash="producer-$(consumer_class)",
            consumer_artifact_hash="consumer-$(consumer_class)",
            data_hashes=["sm-data-$(consumer_class)"],status=:candidate,
        )
        @test binding.material_tag == "SmBCO"
        @test binding.rare_earth_symbol == "Sm"
        @test binding.consumer_class == consumer_class
        @test !binding.physical_qualification
    end
    sa = family_packet("SaBCO","",:unresolved;packet_id="sa-consumers")
    blocked = bind_rebco_family_consumer(
        sa,:deposition;binding_id="sa-deposition",producer_artifact_hash="sa-producer",
        consumer_artifact_hash="sa-consumer",data_hashes=["sa-unresolved"],
        status=:blocked_input,
    )
    @test blocked.status == :blocked_input
    @test_throws ErrorException bind_rebco_family_consumer(
        sa,:pka;binding_id="sa-pka",producer_artifact_hash="sa-producer",
        consumer_artifact_hash="sa-consumer",data_hashes=["sa-unresolved"],
        status=:blocked_input,
    )
    @test_throws ErrorException bind_rebco_family_consumer(
        sa,:source;binding_id="sa-source",producer_artifact_hash="sa-producer",
        consumer_artifact_hash="sa-consumer",data_hashes=["sa-unresolved"],
        status=:candidate,
    )
end

@testset "Multi-REBCO extraction ownership" begin
    manifest = default_hts_addon_extraction_manifest()
    component = get_addon_component(manifest,:multi_rebco_family_contracts)
    @test component.files == ["src/hts_addon/Multi_REBCO_Family.jl"]
    @test validate_addon_extraction_manifest(manifest)
end
