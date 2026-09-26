using Radiant, Test, LinearAlgebra

function fidelity_fixture(;basis="standard", order=0, scale=1.0, pencil=true,
                          edges=[1.0e3,1.0e5,1.0e6], permutation=nothing)
    proton = Proton()
    cs = Cross_Sections()
    cs.particles = [proton]
    cs.number_of_particles = 1
    cs.number_of_groups = [2]
    cs.energy_boundaries = [[1.0,0.1,0.001]]
    geo = Geometry()
    geo.type = "cartesian"
    geo.dimension = 3
    geo.axis = ["x","y","z"]
    for (axis,origin,width) in zip(geo.axis,[-0.4,1.2,-2.3],[2.0,3.0,4.0])
        geo.number_of_voxels[axis] = 1
        geo.voxels_boundaries[axis] = [origin,origin+width]
        geo.voxels_position[axis] = [origin+width/2]
        geo.voxels_width[axis] = [width]
    end
    geo.volume_per_voxel = fill(24.0,1,1,1)
    geo.material_per_voxel = ones(Int64,1,1,1)
    geo.is_build = true
    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,"CSD")
    Radiant.set_quadrature(solver,"gauss-legendre-chebychev",2,3)
    Radiant.set_legendre_order(solver,order)
    Radiant.set_angular_boltzmann(solver,basis)
    _,w,dirs,_ = Radiant._solver_quadrature(solver,geo)
    normals = [-1.0 0 0; 1.0 0 0; 0 -1.0 0; 0 1.0 0; 0 0 -1.0; 0 0 1.0]
    centroids = [-0.4 2.7 -0.3; 1.6 2.7 -0.3; 0.6 1.2 -0.3; 0.6 4.2 -0.3; 0.6 2.7 -2.3; 0.6 2.7 1.7]
    t1 = [0.0 1 0; 0 1 0; 0 0 1; 0 0 1; 1 0 0; 1 0 0]
    t2 = [0.0 0 -1; 0 0 1; -1 0 0; 1 0 0; 0 -1 0; 0 1 0]
    flux = zeros(6,2,length(w))
    selected = Int[]
    for p in 1:6
        incoming = findall(d -> dot(dirs[d,:],normals[p,:]) < 0,eachindex(w))
        push!(selected,incoming[1])
        for g in 1:2, d in incoming
            flux[p,g,d] = scale*g*(pencil ? (d == incoming[1] ? 1.0 : 0.0) :
                1.0+0.2*dirs[d,1]-0.1*dirs[d,2]+0.15*dirs[d,3])
        end
    end
    perm = isnothing(permutation) ? collect(eachindex(w)) : permutation
    source = Boundary_Angular_Current_Source(proton,collect(1:6),centroids,[12.,12.,8.,8.,6.,6.],
        normals,t1,t2,edges,dirs[perm,:],w[perm],flux[:,:,perm],
        Source_Normalization(basis=:per_history,source_hash="six-shifted-fidelity"))
    return source,cs,geo,solver,dirs,w,selected
end

@testset "Boundary angular source fidelity" begin
    source,cs,geo,solver,dirs,w,selected = fidelity_fixture()
    projected,receipt = project_boundary_source(source,cs,geo,solver)
    diagnostic = Radiant.boundary_source_fidelity(source,projected,cs,geo,solver)
    @test receipt.energy_group_map == [2,1]
    @test diagnostic.energy_group_map == [2,1]
    @test receipt.surface_index == collect(1:6)
    @test receipt.surface_cell == fill((1,1,1),6)
    @test receipt.max_relative_error ≤ 1e-10
    @test diagnostic.max_current_relative_error ≤ 1e-10
    @test all(diagnostic.distribution_l1_error .> 0.1)
    @test all(diagnostic.first_moment_error .> 0.1)
    @test all(diagnostic.second_moment_error .> 0.01)
    @test minimum(diagnostic.minimum_reconstructed_flux) ≥ 0
    @test_throws ErrorException Radiant.assert_boundary_source_fidelity(diagnostic)
    # Independent intended moments and extensive current, including unequal face areas.
    for p in 1:6, sg in 1:2
        g = 3-sg
        d = selected[p]
        expected = source.areas_cm2[p]*w[d]*abs(dot(dirs[d,:],source.normals[p,:]))*sg
        @test diagnostic.target_current[p,g] ≈ expected atol=1e-12
        @test diagnostic.target_first_moment[p,g,:] ≈ dirs[d,:] atol=1e-12
        @test diagnostic.target_second_moment[p,g,:,:] ≈ dirs[d,:]*dirs[d,:]' atol=1e-12
    end
    println("L0 direction erasure: max_current_error=",diagnostic.max_current_relative_error,
        "; min_distribution_L1=",minimum(diagnostic.distribution_l1_error),
        "; min_first_error=",minimum(diagnostic.first_moment_error),
        "; min_second_error=",minimum(diagnostic.second_moment_error))

    # Square half-range Galerkin inversion resolves the same selected ordinate on all faces.
    for pencil in (true,false)
        rich,cs,geo,solver,dirs,w,_ = fidelity_fixture(basis="galerkin-d",order=2,pencil=pencil,
            permutation=reverse(collect(eachindex(w))))
        installed,receipt = project_boundary_source(rich,cs,geo,solver)
        diag = Radiant.boundary_source_fidelity(rich,installed,cs,geo,solver)
        @test Radiant.assert_boundary_source_fidelity(diag;positivity_atol=1e-10)
        @test maximum(diag.distribution_l1_error) ≤ 1e-10
        @test maximum(diag.first_moment_error) ≤ 1e-10
        @test maximum(diag.second_moment_error) ≤ 1e-10
        println("Galerkin six-face ",pencil ? "pencil" : "smooth",": max_L1=",maximum(diag.distribution_l1_error))

        triple,_,_,_,_,_,_ = fidelity_fixture(basis="galerkin-d",order=2,pencil=pencil,scale=3,
            permutation=reverse(collect(eachindex(w))))
        installed3,receipt3 = project_boundary_source(triple,cs,geo,solver)
        diag3 = Radiant.boundary_source_fidelity(triple,installed3,cs,geo,solver)
        @test receipt3.target_current ≈ 3 .* receipt.target_current
        @test diag3.target_first_moment ≈ diag.target_first_moment
        @test diag3.target_second_moment ≈ diag.target_second_moment
        @test Radiant.assert_boundary_source_fidelity(diag3;positivity_atol=1e-10)
        # Additive installation of three unit-strength sources equals one strength-three source.
        summed = Radiant._empty_surface_source_array(2,size(installed,2),geo)
        for g in axes(installed,1), p in axes(installed,2), face in axes(installed,3)
            summed[g,p,face] = installed[g,p,face]+installed[g,p,face]+installed[g,p,face]
            @test summed[g,p,face] ≈ installed3[g,p,face] atol=1e-11
        end
        @test Radiant.assert_boundary_source_fidelity(
            Radiant.boundary_source_fidelity(triple,summed,cs,geo,solver);positivity_atol=1e-10)
        # Wrong magnitude, signed array and missing/wrong basis must not pass silently.
        @test_throws ErrorException Radiant.assert_boundary_source_fidelity(
            Radiant.boundary_source_fidelity(triple,installed,cs,geo,solver))
        negative = deepcopy(installed)
        for g in axes(negative,1), p in axes(negative,2), face in axes(negative,3)
            negative[g,p,face] .*= -1
        end
        @test_throws ErrorException Radiant.assert_boundary_source_fidelity(
            Radiant.boundary_source_fidelity(rich,negative,cs,geo,solver))
        @test_throws ErrorException Radiant.boundary_source_fidelity(rich,installed[:,1:1,:],cs,geo,solver)
        @test_throws ErrorException Radiant.assert_boundary_source_fidelity(diag;first_moment_atol=NaN)
    end

    # The diagnostic must reject installed nonfinite coefficients, even when the source itself is valid.
    rich,cs,geo,solver,_,_,_ = fidelity_fixture(basis="galerkin-d",order=2)
    installed,_ = project_boundary_source(rich,cs,geo,solver)
    nonfinite = deepcopy(installed)
    nonfinite[1,1,1][1,1] = NaN
    @test_throws ErrorException Radiant.boundary_source_fidelity(rich,nonfinite,cs,geo,solver)

    # An absent source group must not acquire a plausible-looking reconstructed current.
    subset = Boundary_Angular_Current_Source(
        rich.particle,rich.patch_ids,rich.centroids_cm,rich.areas_cm2,rich.normals,
        rich.tangent_1,rich.tangent_2,[1.0e3,1.0e5],rich.directions,
        rich.quadrature_weights,copy(rich.angular_flux[:,1:1,:]),rich.normalization;
        variance=isnothing(rich.variance) ? nothing : copy(rich.variance[:,1:1,:]),
        provenance=copy(rich.provenance),
    )
    subset_projection,_ = project_boundary_source(subset,cs,geo,solver)
    subset_diag = Radiant.boundary_source_fidelity(subset,subset_projection,cs,geo,solver)
    @test all(iszero,subset_diag.target_current[:,1])
    spurious = deepcopy(subset_projection)
    spurious[1,1,1][1,1] = 1.0
    spurious_diag = Radiant.boundary_source_fidelity(subset,spurious,cs,geo,solver)
    @test spurious_diag.reconstructed_current[1,1] > 0
    @test_throws ErrorException Radiant.assert_boundary_source_fidelity(spurious_diag;positivity_atol=1e-10)

    # Fidelity's patch-to-face mapping is explicitly Cartesian-only.
    unsupported_geo = deepcopy(geo)
    unsupported_geo.type = "cylindrical"
    @test_throws ErrorException Radiant.boundary_source_fidelity(
        rich,installed,cs,unsupported_geo,solver)

    wrong,cs,geo,solver,_,_,_ = fidelity_fixture(edges=[1.0e3,1.1e5,1.0e6])
    @test_throws ErrorException project_boundary_source(wrong,cs,geo,solver)
    @test_throws ErrorException Radiant.boundary_source_fidelity(wrong,projected,cs,geo,solver)
    empty,cs,geo,solver,_,_,_ = fidelity_fixture(scale=0)
    installed,_ = project_boundary_source(empty,cs,geo,solver)
    empty_diag = Radiant.boundary_source_fidelity(empty,installed,cs,geo,solver)
    @test Radiant.assert_boundary_source_fidelity(empty_diag)
    @test all(iszero,empty_diag.target_current)
    @test all(iszero,empty_diag.target_first_moment)
    @test all(iszero,empty_diag.target_second_moment)
end
