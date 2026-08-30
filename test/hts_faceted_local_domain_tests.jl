using Radiant
using Test
using LinearAlgebra

@testset "Faceted local HTS tape domain" begin
    patch = Faceted_Tape_Patch(
        1,
        101,
        (0.0,0.0,0.0),
        0.4,
        Matrix{Float64}(I,3,3),
        0.0,
        "tape-surface",
    )
    layer = Layer_Definition(
        "rebco","YBCO",1.0e-4;
        role=:superconductor,
        transport_subdivisions=2,
        activation_enabled=true,
        atomistic_enabled=true,
    )
    stack = Tape_Stack_Definition(
        "one-layer-facet-fixture",0.4,[layer];
        metadata=Dict("classification" => "software-verification"),
    )
    domain = build_faceted_local_tape_domain(
        patch,stack,Dict("YBCO" => 1);
        along_cells=2,
        width_cells=2,
        layer_refinement=2,
        source_mesh_hash="analytic-facet-mesh",
    )

    @test domain.surrogate_length_cm ≈ 1.0 atol=1.0e-14
    @test domain.area_ratio ≈ 1.0 atol=1.0e-14
    @test size(domain.geometry.material_per_voxel) == (2,2,4)
    @test sum(domain.geometry.volume_per_voxel) ≈
          patch.area_cm2*get_total_thickness(stack) rtol=1.0e-12
    @test Radiant.get_boundary_conditions(domain.geometry) == zeros(Int64,6)
    @test domain.layer_voxel_ranges["rebco"] == 1:4
    @test local_layer_for_voxel(domain,3).name == "rebco"
    @test faceted_local_domain_receipt(domain)["physical_curvature_qualification"] == false

    local_position = [0.2,-0.1,2.0e-5]
    global_position = local_to_global(domain,local_position)
    @test global_to_local(domain,global_position) ≈ local_position atol=1.0e-14

    photon = Photon()
    normalization = Source_Normalization(
        basis=:per_history,
        source_rate_per_s=1.0,
        symmetry_factor=1.0,
        source_hash="analytic-facet-source",
        provenance=Dict("classification" => "software-verification"),
    )
    angular_flux = reshape([2.5],1,1,1)
    source = Boundary_Angular_Current_Source(
        photon,
        [11],
        [0.0 0.0 0.0],
        [patch.area_cm2],
        [0.0 0.0 1.0],
        [1.0 0.0 0.0],
        [0.0 1.0 0.0],
        [1.0e3,1.0e6],
        [0.0 0.0 -1.0],
        [1.0],
        angular_flux,
        normalization,
    )
    typed_mapping = Facet_Boundary_Map(
        source.normalization.source_hash,
        domain.source_mesh_hash,
        [Facet_Boundary_Mapping_Entry(
            source_patch_id=11,
            facet_index=patch.facet_index,
            surface_id=patch.surface_id,
            centroid_distance_cm=0.0,
            normal_dot=1.0,
            facet_area_cm2=patch.area_cm2,
            source_area_cm2=patch.area_cm2,
        )],
    )
    local_sources,receipts = localize_boundary_source_to_facet(
        source,typed_mapping,domain,
    )
    @test Set(keys(local_sources)) == Set([:back])
    @test Set(keys(receipts)) == Set([:back])
    @test get_incoming_current(local_sources[:back]) ≈ get_incoming_current(source)
    @test receipts[:back].maximum_current_relative_error ≤ 1.0e-12
    @test receipts[:back].area_relative_error ≤ 1.0e-12
    @test local_sources[:back].normals ≈ [0.0 0.0 1.0]
    @test local_sources[:back].directions ≈ [0.0 0.0 -1.0]

    dictionary_mapping = Dict{String,Any}(
        "schema" => "radiant.faceted_source_mapping/v1",
        "source_hash" => source.normalization.source_hash,
        "geometry_hash" => domain.source_mesh_hash,
        "mappings" => [Dict{String,Any}(
            "source_patch_id" => 11,
            "facet_index" => 1,
            "surface_id" => 101,
            "centroid_distance_cm" => 0.0,
            "normal_dot" => 1.0,
            "facet_area_cm2" => patch.area_cm2,
            "source_area_cm2" => patch.area_cm2,
        )],
    )
    replayed_mapping = facet_boundary_map(dictionary_mapping)
    replayed_sources,replayed_receipts = localize_boundary_source_to_facet(
        source,replayed_mapping,domain,
    )
    @test get_incoming_current(replayed_sources[:back]) ≈ get_incoming_current(source)
    @test replayed_receipts[:back].local_source_hash == receipts[:back].local_source_hash

    @test_throws ErrorException Facet_Boundary_Map(
        source.normalization.source_hash,
        "wrong-geometry-hash",
        typed_mapping.mappings,
    ) |> mapping -> localize_boundary_source_to_facet(source,mapping,domain)
end
