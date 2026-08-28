using Radiant
using Test

function tape_material_library(stack::Tape_Stack_Definition)
    tags = unique(layer.material_tag for layer in stack.layers)
    materials = Material[]
    for tag in tags
        material = Material(tag)
        Radiant.set_density(material,1.0)
        push!(materials,material)
    end
    cross_sections = Cross_Sections()
    Radiant.set_materials(cross_sections,materials)
    return cross_sections
end

@testset "Executable planar HTS tape geometry" begin
    stack = verification_eight_layer_stack()
    cross_sections = tape_material_library(stack)
    mapping = build_planar_tape_geometry(
        stack,cross_sections;
        dimension=3,
        length_cm=0.1,
        thickness_refinement=2,
        width_subdivisions=2,
        length_subdivisions=3,
    )
    geometry = mapping.geometry

    @test geometry.is_build
    @test Radiant.get_dimension(geometry) == 3
    @test geometry.number_of_voxels["x"] == 16
    @test geometry.number_of_voxels["y"] == 2
    @test geometry.number_of_voxels["z"] == 3
    @test length(mapping.layer_voxel_ranges_x) == 8
    @test mapping.layer_voxel_ranges_x[3] == 5:6
    @test mapping.metadata["curvature"] == "not-included"

    represented_volume = sum(geometry.volume_per_voxel)
    expected_volume = get_total_thickness(stack)*stack.width_cm*0.1
    @test represented_volume ≈ expected_volume atol=1.0e-15

    @test length(get_layer_voxel_ids(mapping,3)) == 12
    @test length(get_activation_voxel_ids(mapping)) == 12
    @test length(get_atomistic_voxel_ids(mapping)) == 12

    values = zeros(Float64,16,2,3)
    for ix in 1:16
        values[ix,:,:] .= mapping.voxel_layer_index_x[ix]
    end
    @test aggregate_layer_response(
        mapping,values;reduction=:volume_average,
    ) ≈ collect(1.0:8.0)
    @test aggregate_layer_response(
        mapping,values;reduction=:maximum,
    ) ≈ collect(1.0:8.0)
    @test aggregate_layer_response(
        mapping,ones(Float64,16,2,3);reduction=:sum,
    ) == fill(12.0,8)

    one_dimensional = build_planar_tape_geometry(
        stack,cross_sections;dimension=1,thickness_refinement=1,
    )
    @test one_dimensional.geometry.number_of_voxels["x"] == 8
    @test sum(one_dimensional.geometry.volume_per_voxel) ≈ get_total_thickness(stack)

    @test_throws ErrorException build_planar_tape_geometry(
        stack,cross_sections;dimension=3,
    )
end
