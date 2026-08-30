using Radiant
using Test
using LinearAlgebra

@testset "Response-preserving energy-group condensation" begin
    mapping = Energy_Group_Condensation_Map(
        [0.0,1.0,2.0,3.0,4.0],
        [0.0,2.0,4.0],
    )
    @test fine_group_count(mapping) == 4
    @test coarse_group_count(mapping) == 2
    @test mapping.fine_to_coarse == [1,1,2,2]
    @test mapping.coarse_members == [[1,2],[3,4]]
    @test condense_group_integrals(mapping,[1.0,2.0,3.0,4.0]) == [3.0,7.0]

    flux = [1.0,2.0,3.0,4.0]
    response = [10.0,20.0,30.0,40.0]
    condensed,receipt = condense_response_coefficients(
        mapping,flux,response;
        response_id="REBCO-prompt-heating",
        physical_reference=false,
        metadata=Dict("classification" => "analytic-verification"),
    )
    @test condensed ≈ [50.0/3.0,250.0/7.0] rtol=1.0e-14
    @test dot(flux,response) == 300.0
    @test dot(condense_group_integrals(mapping,flux),condensed) ≈ 300.0 atol=1.0e-13
    @test receipt.fine_response_integral == 300.0
    @test receipt.coarse_response_integral ≈ 300.0 atol=1.0e-13
    @test receipt.relative_residual <= 1.0e-15
    @test !receipt.physical_reference
    receipt_map = response_condensation_receipt(receipt)
    @test receipt_map["metadata"]["spectrum_specific"] == "true"

    covariance = Matrix(Diagonal([1.0,4.0,9.0,16.0]))
    condensed_covariance = condense_response_covariance(mapping,flux,covariance)
    @test condensed_covariance[1,1] ≈ 17.0/9.0 rtol=1.0e-14
    @test condensed_covariance[2,2] ≈ 337.0/49.0 rtol=1.0e-14
    @test condensed_covariance[1,2] == 0.0
    @test condensed_covariance == condensed_covariance'

    # A protected response cannot be assigned in an energy interval absent from the reference
    # spectrum. An explicit alternate weighting spectrum is required instead.
    @test_throws ErrorException condense_response_coefficients(
        mapping,[0.0,0.0,3.0,4.0],response;
        response_id="invalid-zero-flux-group",
    )
    @test_throws ErrorException Energy_Group_Condensation_Map(
        [0.0,1.0,2.0,3.0],[0.0,1.5,3.0],
    )
end

@testset "Neutron-to-photon transfer condensation" begin
    incoming_map = Energy_Group_Condensation_Map(
        [0.0,1.0,2.0,3.0,4.0],
        [0.0,2.0,4.0],
    )
    outgoing_map = Energy_Group_Condensation_Map(
        [0.0,10.0,20.0,30.0,40.0],
        [0.0,20.0,40.0],
    )
    neutron_flux = [1.0,2.0,3.0,4.0]
    production = [
        1.0  2.0  3.0  4.0
        0.5  1.0  1.5  2.0
        2.0  1.0  0.5  0.25
        4.0  3.0  2.0  1.0
    ]
    condensed,receipt = condense_transfer_matrix(
        incoming_map,outgoing_map,neutron_flux,production;
        response_id="neutron-to-photon-production",
        physical_reference=false,
        metadata=Dict(
            "units" => "photons/source-neutron",
            "classification" => "analytic-verification",
        ),
    )
    @test size(condensed) == (2,2)
    fine_source = vec(transpose(neutron_flux)*production)
    expected_outgoing = [sum(fine_source[1:2]),sum(fine_source[3:4])]
    coarse_flux = condense_group_integrals(incoming_map,neutron_flux)
    observed_outgoing = vec(transpose(coarse_flux)*condensed)
    @test observed_outgoing ≈ expected_outgoing rtol=1.0e-14 atol=1.0e-14
    @test receipt.aggregated_fine_outgoing_source ≈ expected_outgoing
    @test receipt.coarse_outgoing_source ≈ expected_outgoing
    @test receipt.maximum_relative_residual <= 1.0e-15
    @test !receipt.physical_reference
    receipt_map = transfer_matrix_condensation_receipt(receipt)
    @test receipt_map["metadata"]["matrix_orientation"] ==
          "incoming-group-by-outgoing-group"
    @test receipt_map["metadata"]["spectrum_specific"] == "true"
end

@testset "Independent protected responses use independent condensations" begin
    mapping = Energy_Group_Condensation_Map(
        [0.0,1.0,2.0,3.0,4.0],
        [0.0,2.0,4.0],
    )
    flux = [1.0,2.0,3.0,4.0]
    heating = [1.0,4.0,2.0,8.0]
    gd_capture = [100.0,10.0,1.0,0.1]
    heating_coarse,heating_receipt = condense_response_coefficients(
        mapping,flux,heating;response_id="REBCO-heating",
    )
    capture_coarse,capture_receipt = condense_response_coefficients(
        mapping,flux,gd_capture;response_id="Gd-157-capture",
    )
    @test heating_coarse != capture_coarse
    @test heating_receipt.coarse_response_integral ≈ dot(flux,heating)
    @test capture_receipt.coarse_response_integral ≈ dot(flux,gd_capture)
    @test heating_receipt.reference_flux_hash == capture_receipt.reference_flux_hash
    @test heating_receipt.group_map_hash == capture_receipt.group_map_hash
end
