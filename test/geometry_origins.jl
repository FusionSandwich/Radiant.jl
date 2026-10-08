using Test
using Radiant

@testset "Native Cartesian shifted origins and unequal cells" begin
    mat = Water()
    cs = Cross_Sections()
    cs.set_materials([mat])
    for dimension in 1:3, shift in (0.0, -7.25, 3.125)
        geo = Geometry()
        geo.set_dimension(dimension)
        axes = ["x", "y", "z"][1:dimension]
        for (i,axis) in enumerate(axes)
            origin = shift + 0.25*(i-1)
            geo.set_number_of_regions(axis,2)
            geo.set_voxels_per_region(axis,[2,1])
            geo.set_region_boundaries(axis,[origin,origin+1.0,origin+3.5])
            geo.set_boundary_conditions(uppercase(axis)*"-","void")
            geo.set_boundary_conditions(uppercase(axis)*"+","void")
        end
        geo.set_material_per_region(fill(mat, ntuple(_ -> 2,dimension)))
        geo.build(cs)
        for (i,axis) in enumerate(axes)
            origin = shift + 0.25*(i-1)
            @test geo.get_voxels_boundaries(axis) ≈ origin .+ [0.0,0.5,1.0,3.5]
            @test geo.get_voxels_position(axis) ≈ origin .+ [0.25,0.75,2.25]
            @test geo.voxels_width[axis] ≈ [0.5,0.5,2.5]
        end
        @test sum(geo.volume_per_voxel) ≈ 3.5^dimension
        @test all(==(1),geo.material_per_voxel)
    end
end
