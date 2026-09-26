using Radiant, Test, LinearAlgebra

function origin_fixture(dimension, origins; multiple=false)
    mat = Material("origin-test")
    cs = Cross_Sections()
    Radiant.set_materials(cs,mat)
    geo = Geometry()
    Radiant.set_dimension(geo,dimension)
    for (i,axis) in enumerate(["x","y","z"][1:dimension])
        edges = multiple ? origins[i] .+ [0.0,0.25,1.0] : origins[i] .+ [0.0,1.0]
        Radiant.set_number_of_regions(geo,axis,length(edges)-1)
        Radiant.set_region_boundaries(geo,axis,edges)
        Radiant.set_voxels_per_region(geo,axis,multiple ? [1,3] : [1])
        for side in ("-","+")
            Radiant.set_boundary_conditions(geo,axis*side,"void")
        end
    end
    Radiant.set_material_per_region(geo,fill(mat,ntuple(_ -> multiple ? 2 : 1,dimension)))
    Radiant.build(geo,cs)
    return geo,cs
end

@testset "Cartesian origins in every active dimension" begin
    for dimension in 1:3, origins in ([0.0,0.0,0.0],[-0.4,-0.2,-0.3],[0.5,1.0,2.0]), multiple in (false,true)
        geo,_ = origin_fixture(dimension,origins;multiple=multiple)
        offsets = multiple ? [0.0,0.25,0.5,0.75,1.0] : [0.0,1.0]
        for (i,axis) in enumerate(["x","y","z"][1:dimension])
            edges = geo.voxels_boundaries[axis]
            @test edges ≈ origins[i] .+ offsets
            @test diff(edges) ≈ geo.voxels_width[axis]
            @test (edges[1:end-1]+edges[2:end])/2 ≈ geo.voxels_position[axis]
        end
        @test sum(geo.volume_per_voxel) ≈ 1.0
    end
end

@testset "Shifted proton projection on all six faces" begin
    geo,cs = origin_fixture(3,[-0.4,-0.2,-0.3])
    proton = Proton()
    Radiant.set_particles(cs,proton)
    Radiant.set_energy_boundaries(cs,[[1.0,0.001]])
    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,"CSD")
    Radiant.set_quadrature(solver,"gauss-legendre-chebychev",2,3)
    Radiant.set_legendre_order(solver,0)
    Radiant.set_angular_boltzmann(solver,"standard")
    omega,w = Radiant.quadrature(2,"gauss-legendre-chebychev",3,3)
    directions = hcat(omega...)
    centroids = [-0.4 0.3 0.2; 0.6 0.3 0.2; 0.1 -0.2 0.2; 0.1 0.8 0.2; 0.1 0.3 -0.3; 0.1 0.3 0.7]
    normals = [-1.0 0 0; 1.0 0 0; 0 -1.0 0; 0 1.0 0; 0 0 -1.0; 0 0 1.0]
    t1 = [0.0 1 0; 0 1 0; 0 0 1; 0 0 1; 1 0 0; 1 0 0]
    t2 = [0.0 0 -1; 0 0 1; -1 0 0; 1 0 0; 0 -1 0; 0 1 0]
    values = zeros(6,1,length(w))
    for p in 1:6, n in eachindex(w)
        values[p,1,n] = dot(directions[n,:],normals[p,:]) < 0 ? 1.0 : 0.0
    end
    source = Boundary_Angular_Current_Source(proton,collect(1:6),centroids,ones(6),normals,t1,t2,
        [1.0e3,1.0e6],directions,w,values,
        Source_Normalization(basis=:per_history,source_hash="six-shifted-faces"))
    _,receipt = project_boundary_source(source,cs,geo,solver)
    @test receipt.surface_index == collect(1:6)
    @test receipt.surface_cell == fill((1,1,1),6)
    @test all(receipt.target_current .> 0)
    @test receipt.target_current ≈ receipt.projected_current atol=1e-12
    @test receipt.max_relative_error <= 1e-10
end
