using Test, SHA, LinearAlgebra

# Kept separate from the normal runner until the root execution freeze.
function source_api_problem(;n=2, groups=[10.0,5.0,1.0], stopping=0.0,
                            dimension=3, angles=2, cap=100, origin=-0.7)
    p = Proton(); materials = Material[]; tables = Proton_Material_Data[]
    removals = Proton_Nonelastic_Data[]
    for (i,density) in enumerate((1.0,2.0))
        id = "source-api-manufactured-$i"; m = Material(id)
        Radiant.set_density(m,density); push!(materials,m)
        h = bytes2hex(sha256("$id/$stopping"))
        model = Tabulated_Ion_Transport_Model(species_id="proton",material_id=id,
            energy_MeV=[1.0,10.0],electronic_stopping_MeV_cm=fill(stopping,2),
            nuclear_stopping_MeV_cm=zeros(2),energy_straggling_variance_MeV2_cm=zeros(2),
            angular_variance_rad2_cm=zeros(2),data_hash=h,qualification_status=:synthetic)
        push!(tables,Proton_Material_Data(m,model;material_state="solid",
            density_g_cm3=density,source_sha256=h))
        push!(removals,Proton_Nonelastic_Data(material_id=id,energy_MeV=[1.0,10.0],
            removal_cm_inv=zeros(2),source_sha256=h))
    end
    binding = bind_proton_native(p,tables,removals,groups)
    cs = binding.cross_sections; Radiant.build(cs)
    geo = Geometry(); Radiant.set_dimension(geo,dimension)
    Radiant.set_number_of_regions(geo,"x",2)
    Radiant.set_region_boundaries(geo,"x",[origin,origin+0.5,origin+1.0])
    Radiant.set_voxels_per_region(geo,"x",[n÷2,n÷2])
    for axis in ("y","z")[1:dimension-1]
        Radiant.set_number_of_regions(geo,axis,1)
        a = axis == "y" ? 1.3 : -2.4
        Radiant.set_region_boundaries(geo,axis,[a,a+1.0])
        Radiant.set_voxels_per_region(geo,axis,[1])
    end
    for axis in ("x","y","z")[1:dimension], side in ("-","+")
        Radiant.set_boundary_conditions(geo,axis*side,"void")
    end
    Radiant.set_material_per_region(geo,dimension == 1 ? materials : reshape(materials,2,1,1))
    Radiant.build(geo,cs)
    sn = SN(); Radiant.set_particle(sn,p); Radiant.set_solver_type(sn,"CSD")
    if dimension == 1
        Radiant.set_quadrature(sn,"gauss-legendre",angles)
        Radiant.set_legendre_order(sn,angles-1)
    else
        Radiant.set_quadrature(sn,"gauss-legendre-chebychev",angles,3)
        Radiant.set_legendre_order(sn,angles)
    end
    Radiant.set_angular_boltzmann(sn,"galerkin-d")
    for axis in ("x","y","z")[1:dimension]
        Radiant.set_scheme(sn,axis,"DG",1)
    end
    Radiant.set_scheme(sn,"E","DG",1)
    Radiant.set_maximum_iteration(sn,cap); Radiant.set_convergence_criterion(sn,1e-10)
    solvers = Solvers(); Radiant.add_solver(solvers,sn)
    return (;p,binding,cs,geo,sn,solvers)
end

function source_api_convergence(data)
    all(r -> r.converged && r.iterations <= r.iteration_cap &&
        isfinite(r.solver_residual) && r.solver_residual <= r.tolerance &&
        isfinite(r.reconstruction_residual) && r.reconstruction_residual <= r.tolerance,
        data.convergence)
end

function source_api_solve(problem, specs)
    fs = Fixed_Sources(problem.cs,problem.geo,problem.solvers)
    for s in specs; Radiant.add_source(fs,s); end
    Radiant.build(fs)
    f = Radiant.transport(problem.cs,problem.geo,problem.solvers,fs;retain_boundary_flux=true)
    data = only(Radiant.get_boundary_flux(f,problem.p))
    @test source_api_convergence(data)
    @test all(isfinite,Radiant.get_outgoing_current(f,problem.p))
    @test minimum(minimum(v[:,:,1,:,:]) for v in data.faces) >= -1e-10
    return (;fs,f,data,current=Radiant.get_outgoing_current(f,problem.p))
end

@testset "Native explicit Fixed_Sources strength, moments and bins" begin
  for dimension in (1,3)
    problem = source_api_problem(;dimension=dimension)
    _,w,dirs,_ = Radiant._solver_quadrature(problem.sn,problem.geo)
    nd = length(w); norm = Source_Normalization(basis=:per_source_particle,
        source_hash="source-api-v2-manufactured")
    # Unequal widths, unequal bin values. Values are integrated over the bin,
    # not a per-MeV density; no extra width factor belongs in this integral.
    values = zeros(2,2,nd)
    for i in 1:2, g in 1:2, d in 1:nd
        values[i,g,d] = g*(1+0.2*dirs[d,1]-0.1*dirs[d,2]+0.05*dirs[d,3])
    end
    volumes = vec(problem.geo.volume_per_voxel)
    volume(scale) = Anisotropic_Volume_Source(problem.p,[1,2],volumes,[1e6,5e6,10e6],
        :ordinates,scale .* values,norm;directions=dirs,quadrature_weights=w)
    expected = [sum(volumes[i]*w[d]*values[i,g,d] for i in 1:2,d in 1:nd) for g in 1:2]
    @test get_volume_source_rate(volume(1)) ≈ expected atol=1e-11
    projected,receipt = project_volume_source(volume(1),problem.cs,problem.geo,problem.sn)
    @test receipt.energy_group_map == [2,1]
    @test receipt.target_rate ≈ reverse(expected) atol=1e-11
    Ω = [dirs[:,a] for a in 1:3]
    basis_order = Radiant.get_legendre_order(problem.sn)
    qdim = Radiant.get_quadrature_dimension(problem.sn,dimension)
    _,Mn,_,_,_ = Radiant.angular_polynomial_basis(Ω,w,basis_order,"galerkin-d",qdim)
    for i in 1:2,g in 1:2
        reconstructed = Mn*projected[3-g,:,1,i,1,1]
        target = values[i,g,:]
        @test reconstructed ≈ target atol=1e-10
        for a in 1:3
            @test sum(w .* dirs[:,a] .* reconstructed) ≈ sum(w .* dirs[:,a] .* target) atol=1e-10
            for b in 1:3
                @test sum(w .* dirs[:,a] .* dirs[:,b] .* reconstructed) ≈
                    sum(w .* dirs[:,a] .* dirs[:,b] .* target) atol=1e-10
            end
        end
    end
    # Native boundary install and native transport, on a shifted heterogeneous
    # material map. The zero coefficients deliberately isolate the source API.
    incoming = zeros(1,2,nd)
    for g in 1:2,d in 1:nd
        dirs[d,1] > 0 && (incoming[1,g,d] = g*(1+0.1*dirs[d,2]))
    end
    boundary(scale) = Boundary_Angular_Current_Source(problem.p,[1],
        reshape([-0.7,1.8,-1.9],1,3),[1.0],reshape([-1.,0.,0.],1,3),
        reshape([0.,1.,0.],1,3),reshape([0.,0.,-1.],1,3),[1e6,5e6,10e6],
        dirs,w,scale .* incoming,norm)
    expected_boundary = [sum(w[d]*max(dirs[d,1],0)*incoming[1,g,d] for d in 1:nd) for g in 1:2]
    @test vec(Radiant.get_incoming_current(boundary(1))) ≈ expected_boundary atol=1e-11
    for (kind,maker,truth) in (("volume",volume,expected),("boundary",boundary,expected_boundary))
        one = source_api_solve(problem,[maker(1)])
        three = source_api_solve(problem,[maker(3)])
        added = source_api_solve(problem,[maker(1),maker(3)])
        @test Radiant.get_normalization_factor(one.fs) == 1.0
        @test vec(sum(one.current;dims=1)) ≈ reverse(truth) atol=1e-8 rtol=1e-9
        @test three.current ≈ 3one.current atol=1e-8 rtol=1e-9
        @test added.current ≈ 4one.current atol=1e-8 rtol=1e-9
        @test length(get_projection_receipts(added.fs)) == 2
        old = deepcopy(Radiant.get_source(added.fs,problem.p).volume_sources)
        Radiant.build(added.fs)
        @test Radiant.get_source(added.fs,problem.p).volume_sources == old
        println("SOURCE_API kind=",kind," independent_group_source=",reverse(truth),
            " outgoing=",vec(sum(one.current;dims=1))," strength_cases=1,3,1+3")
        println("SOURCE_API_CONVERGENCE kind=",kind," rows=",one.data.convergence)
    end
    installed,_ = project_boundary_source(boundary(1),problem.cs,problem.geo,problem.sn)
    diag = Radiant.boundary_source_fidelity(boundary(1),installed,problem.cs,problem.geo,problem.sn)
    @test Radiant.assert_boundary_source_fidelity(diag;positivity_atol=1e-10)
    # Shape rejection is separate from current conservation. This volume pencil
    # demonstrates L0 erasure without claiming a continuous-pencil benchmark.
    pencil = zeros(2,2,nd); pencil[:,:,1] .= 1.0
    discrete = Anisotropic_Volume_Source(problem.p,[1,2],volumes,[1e6,5e6,10e6],
        :ordinates,pencil,norm;directions=dirs,quadrature_weights=w)
    low = deepcopy(problem.sn); Radiant.set_legendre_order(low,0); Radiant.set_angular_boltzmann(low,"standard")
    lowprojection,_ = project_volume_source(discrete,problem.cs,problem.geo,low)
    _,lowMn,_,_,_ = Radiant.angular_polynomial_basis(Ω,w,0,"standard",qdim)
    erased = lowMn*lowprojection[1,:,1,1,1,1]
    @test sum(w .* erased) ≈ w[1] atol=1e-10
    @test sum(w .* abs.(erased-pencil[1,2,:]))/w[1] > 0.5
    full,_ = project_volume_source(discrete,problem.cs,problem.geo,problem.sn)
    @test Mn*full[1,:,1,1,1,1] ≈ pencil[1,2,:] atol=1e-10
  end
end
