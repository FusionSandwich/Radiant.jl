using Test, SHA, LinearAlgebra

# Self-contained fixtures: only public source constructors/projectors install Qk.
# The synthetic stopping data qualify API/math behavior, not physical transport.
function public_energy_problem(;dimension=1, energy_order=4, full=true)
    p = Proton(); material = Material("public-energy-source-manufactured")
    Radiant.set_density(material,1.0)
    h = bytes2hex(sha256("public-energy-source-constant-S0.6"))
    model = Tabulated_Ion_Transport_Model(species_id="proton",material_id=material.tag,
        energy_MeV=[1.0,10.0],electronic_stopping_MeV_cm=fill(0.6,2),
        nuclear_stopping_MeV_cm=zeros(2),energy_straggling_variance_MeV2_cm=zeros(2),
        angular_variance_rad2_cm=zeros(2),data_hash=h,qualification_status=:synthetic)
    table = Proton_Material_Data(material,model;material_state="solid",density_g_cm3=1.0,source_sha256=h)
    removal = Proton_Nonelastic_Data(material_id=material.tag,energy_MeV=[1.0,10.0],
        removal_cm_inv=zeros(2),source_sha256=h)
    binding = bind_proton_native(p,[table],[removal],[10.0,5.0,1.0])
    cs = binding.cross_sections; Radiant.build(cs)
    geo = Geometry(); Radiant.set_dimension(geo,dimension)
    for axis in ("x","y","z")[1:dimension]
        origin = axis == "x" ? -0.7 : axis == "y" ? 1.3 : -2.4
        Radiant.set_number_of_regions(geo,axis,1)
        Radiant.set_region_boundaries(geo,axis,[origin,origin+1.0])
        Radiant.set_voxels_per_region(geo,axis,[axis == "x" ? 2 : 1])
        for side in ("-","+"); Radiant.set_boundary_conditions(geo,axis*side,"void"); end
    end
    Radiant.set_material_per_region(geo,dimension == 1 ? [material] : reshape([material],1,1,1))
    Radiant.build(geo,cs)
    sn = SN(); Radiant.set_particle(sn,p); Radiant.set_solver_type(sn,"CSD")
    if dimension == 1
        Radiant.set_quadrature(sn,"gauss-legendre",2); Radiant.set_legendre_order(sn,1)
    else
        Radiant.set_quadrature(sn,"gauss-legendre-chebychev",2,3); Radiant.set_legendre_order(sn,2)
    end
    Radiant.set_angular_boltzmann(sn,"galerkin-d")
    for axis in ("x","y","z")[1:dimension]; Radiant.set_scheme(sn,axis,"DG",2); end
    Radiant.set_scheme(sn,"E","DG",energy_order); Radiant.set_is_full_coupling(sn,full)
    Radiant.set_maximum_iteration(sn,100); Radiant.set_convergence_criterion(sn,1e-10)
    solvers = Solvers(); Radiant.add_solver(solvers,sn)
    return (;p,binding,cs,geo,sn,solvers)
end

function public_energy_base(problem, moments;representation=:isotropic, directions=zeros(0,3),
                            weights=Float64[], provenance=Dict{String,String}())
    norm = Source_Normalization(basis=:per_source_particle,source_rate_per_s=5.0,
        symmetry_factor=2.0,source_hash="public-energy-manufactured-v1")
    return Anisotropic_Volume_Source(problem.p,[1,2],vec(problem.geo.volume_per_voxel),
        [1e6,5e6,10e6],representation,moments[:,:,:,1],norm;
        directions=directions,quadrature_weights=weights,provenance=provenance,source_tolerance=0.0)
end

public_energy_basis(u) = [1.0,sqrt(3.0)*u,sqrt(5.0)*(3u^2-1)/2,sqrt(7.0)*(5u^3-3u)/2]

@testset "Public exact energy source integrals and reconstructed positivity" begin
    edges = [1e6,5e6,10e6]
    coefficients = zeros(2,2,1,4)
    for v in 1:2,g in 1:2
        coefficients[v,g,1,:] .= (v+g)*1e-6 .* [1.0,-0.2,0.3,0.1]
    end
    moments = source_energy_moments(edges,coefficients;energy_order=4)
    # Independent four-point integration, exact for cubic q times cubic basis.
    outer = sqrt((3+2sqrt(6/5))/7); inner = sqrt((3-2sqrt(6/5))/7)
    nodes = [-outer,-inner,inner,outer]
    weights = [(18-sqrt(30.0))/36,(18+sqrt(30.0))/36,(18+sqrt(30.0))/36,(18-sqrt(30.0))/36]
    for v in 1:2,g in 1:2,k in 1:4
        a = coefficients[v,g,1,:]; width = edges[g+1]-edges[g]
        integral = width/2*sum(weights[d]*sum(a[j]*nodes[d]^(j-1) for j in 1:4)*
            public_energy_basis(nodes[d])[k] for d in 1:4)
        @test moments[v,g,1,k] ≈ integral rtol=2e-14 atol=1e-14
    end
    for v in 1:2,g in 1:2,u in (-1.0,-0.71,0.0,0.23,1.0)
        q = sum(coefficients[v,g,1,j]*u^(j-1) for j in 1:4)
        @test sum(moments[v,g,1,:].*public_energy_basis(u))/(edges[g+1]-edges[g]) ≈ q rtol=2e-14
        @test q > 0
    end
    @test all(moments[:,:,:,2] .< 0) # signed coefficients are legitimate.
    padded = source_energy_moments(edges,ones(2,2,1,1);energy_order=4)
    @test all(iszero,padded[:,:,:,2:4]) # explicit degree-zero polynomial, declared DG4.
    @test_throws ErrorException source_energy_moments(edges,coefficients;energy_order=3)
    @test_throws ErrorException source_energy_moments(reverse(edges),coefficients)
    @test_throws ErrorException source_energy_moments([1e6,5e6,5e6],coefficients)
    @test_throws ErrorException source_energy_moments(edges,fill(NaN,2,2,1,4))
    @test_throws ErrorException source_energy_moments(edges,zeros(2,1,1,4))
    problem = public_energy_problem()
    base = public_energy_base(problem,moments)
    wrapper = Energy_Moment_Volume_Source(base,moments)
    @test get_volume_source_rate(wrapper) == get_volume_source_rate(base)
    @test get_volume_source_rate(wrapper;physical=true) == 10get_volume_source_rate(base)
    @test_throws ErrorException Energy_Moment_Volume_Source(base,moments[:,:,:,1])
    @test_throws ErrorException Energy_Moment_Volume_Source(base,zeros(2,2,2,4))
    @test_throws ErrorException Energy_Moment_Volume_Source(base,zeros(2,2,1,0))
    @test_throws ErrorException Energy_Moment_Volume_Source(base,cat(moments,moments[:,:,:,1:1];dims=4))
    @test_throws ErrorException Energy_Moment_Volume_Source(base,fill(Inf,size(moments)))
    changed = copy(moments); changed[1,1,1,1] += 1.0
    @test_throws ErrorException Energy_Moment_Volume_Source(base,changed)
    @test_throws ErrorException Energy_Moment_Volume_Source(base,moments;positivity_tolerance=1e-8)
    # Positive endpoints conceal a negative interior; checking means/endpoints alone is insufficient.
    negative = zeros(2,2,1,4); negative[:,:,:,1] .= 1.0; negative[:,:,:,3] .= 1.0
    @test sum(negative[1,1,1,:].*public_energy_basis(-1.0)) > 0
    @test sum(negative[1,1,1,:].*public_energy_basis(1.0)) > 0
    @test_throws ErrorException Energy_Moment_Volume_Source(public_energy_base(problem,negative),negative)
    # Tangential zero q=u^2: retain raw minimum and a reported numerical sign bound.
    tangent_power = zeros(2,2,1,3); tangent_power[:,:,:,3] .= 1e-6
    tangent = source_energy_moments(edges,tangent_power;energy_order=4)
    tangential = Energy_Moment_Volume_Source(public_energy_base(problem,tangent),tangent)
    _,receipt = project_volume_source(tangential,problem.cs,problem.geo,problem.sn)
    @test all(receipt.reconstructed_energy_minima .>= -receipt.reconstruction_roundoff_bounds)
    @test maximum(receipt.reconstruction_roundoff_bounds) < 1e-11
    genuinely_negative = copy(tangent_power); genuinely_negative[:,:,:,1] .= -1e-14
    bad = source_energy_moments(edges,genuinely_negative;energy_order=4)
    @test_throws ErrorException Energy_Moment_Volume_Source(public_energy_base(problem,bad),bad)
    # Finite huge cubic must not overflow the derivative discriminant and lose its minimum.
    huge = zeros(2,2,1,4)
    for v in 1:2,g in 1:2; huge[v,g,1,:] .= 1e200 .* [0.0292,-0.388,0.94,0.1]; end
    huge_moments = source_energy_moments([0.0,1.0,2.0],huge)
    diagnostic = Radiant._source_energy_polynomial_diagnostic(huge_moments[1,1,1,:])
    @test isfinite(diagnostic.minimum) && diagnostic.minimum < -1e197
    @test diagnostic.minimum < -diagnostic.roundoff_bound
    huge_base = Anisotropic_Volume_Source(problem.p,[1,2],vec(problem.geo.volume_per_voxel),
        [0.0,1.0,2.0],:isotropic,huge_moments[:,:,:,1],base.normalization;source_tolerance=0.0)
    @test_throws ErrorException Energy_Moment_Volume_Source(huge_base,huge_moments)
    # The remote derivative root of a positive subnormal cubic is irrelevant;
    # its finite in-interval stationary point must still be examined.
    nearly_quadratic = zeros(2,2,1,4)
    nearly_quadratic[:,:,:,1] .= 1.0
    nearly_quadratic[:,:,:,3] .= 1.0
    nearly_quadratic[:,:,:,4] .= 1e-320
    tiny_cubic = source_energy_moments(edges,nearly_quadratic)
    tiny_diagnostic = Radiant._source_energy_polynomial_diagnostic(tiny_cubic[1,1,1,:])
    @test isfinite(tiny_diagnostic.minimum) && tiny_diagnostic.minimum > 3.9e6
    @test Energy_Moment_Volume_Source(public_energy_base(problem,tiny_cubic),tiny_cubic) isa Energy_Moment_Volume_Source
end

@testset "Public Cartesian native energy mode placement, rates and angular fidelity" begin
    for order in 1:4
        problem = public_energy_problem(;energy_order=order)
        coefficients = zeros(2,2,1,order)
        for v in 1:2,g in 1:2
            coefficients[v,g,1,:] .= 1e-6 .* [1.0,-0.2,0.3,0.1][1:order]
        end
        moments = source_energy_moments([1e6,5e6,10e6],coefficients;energy_order=order)
        wrapper = Energy_Moment_Volume_Source(public_energy_base(problem,moments),moments)
        projected,receipt = project_volume_source(wrapper,problem.cs,problem.geo,problem.sn)
        @test receipt.native_energy_indices == collect(1:order)
        for v in 1:2,g in 1:2,k in 1:order
            @test projected[g,1,k,v,1,1] == moments[v,3-g,1,k]
        end
        if order == 1
            legacy,_ = project_volume_source(wrapper.source,problem.cs,problem.geo,problem.sn)
            @test projected ≈ legacy rtol=1e-12 atol=1e-12
        end
    end
    for dimension in (1,2,3), full in (false,true)
        problem = public_energy_problem(;dimension=dimension,full=full)
        coefficients = zeros(2,2,1,4)
        for v in 1:2,g in 1:2; coefficients[v,g,1,:] .= (v+g)*1e-6 .* [1.0,-0.2,0.3,0.1]; end
        moments = source_energy_moments([1e6,5e6,10e6],coefficients)
        wrapper = Energy_Moment_Volume_Source(public_energy_base(problem,moments),moments)
        projected,receipt = project_volume_source(wrapper,problem.cs,problem.geo,problem.sn)
        @test receipt.energy_group_map == [2,1]
        @test receipt.native_energy_indices == [1,2,3,4]
        @test receipt.source_hash == wrapper.source.normalization.source_hash
        expected = [sum(wrapper.source.voxel_volumes_cm3[v]*moments[v,3-g,1,k] for v in 1:2) for g in 1:2,k in 1:4]
        @test receipt.target_energy_moment_rates ≈ expected rtol=1e-14 atol=1e-14
        @test receipt.projected_energy_moment_rates ≈ expected rtol=1e-12 atol=1e-12
        @test all(projected[:,:,5:end,:,:,:] .== 0)
        @test all(projected[:,2:end,:,:,:,:] .== 0)
        for v in 1:2,g in 1:2,k in 1:4
            @test projected[g,1,k,v,1,1] == moments[v,3-g,1,k]
        end
        fs = Fixed_Sources(problem.cs,problem.geo,problem.solvers)
        Radiant.add_source(fs,wrapper); Radiant.add_source(fs,Energy_Moment_Volume_Source(
            public_energy_base(problem,3moments),3moments)); Radiant.build(fs)
        @test Radiant.get_normalization_factor(fs) == 1.0
        @test Radiant.get_source(fs,problem.p).volume_sources ≈ 4projected rtol=1e-14
        @test length(get_projection_receipts(fs)) == 2
        saved = copy(Radiant.get_source(fs,problem.p).volume_sources); Radiant.build(fs)
        @test Radiant.get_source(fs,problem.p).volume_sources == saved
        @test_throws ErrorException project_volume_source(wrapper,problem.cs,problem.geo,problem.sn;rate_rtol=NaN)
    end
    problem = public_energy_problem()
    zerosource = zeros(2,2,1,4)
    zerowrapper = Energy_Moment_Volume_Source(public_energy_base(problem,zerosource),zerosource)
    zeroprojected,zeroreceipt = project_volume_source(zerowrapper,problem.cs,problem.geo,problem.sn;rate_atol=0.0)
    @test all(iszero,zeroprojected)
    @test zeroreceipt.max_relative_error == 0.0
    Ω,w,dirs,qdim = Radiant._solver_quadrature(problem.sn,problem.geo)
    _,Mn,Dn,_,_ = Radiant.angular_polynomial_basis(Ω,w,1,"galerkin-d",qdim)
    coefficients = zeros(2,2,length(w),4)
    for v in 1:2,g in 1:2,d in eachindex(w)
        coefficients[v,g,d,:] .= d*1e-6 .* [1.0,-0.2,0.3,0.1]
    end
    moments = source_energy_moments([1e6,5e6,10e6],coefficients)
    base = public_energy_base(problem,moments;representation=:ordinates,directions=dirs,weights=w)
    source = Energy_Moment_Volume_Source(base,moments)
    projected,receipt = project_volume_source(source,problem.cs,problem.geo,problem.sn)
    for v in 1:2,g in 1:2,k in 1:4
        @test Mn*projected[3-g,:,k,v,1,1] ≈ moments[v,g,:,k] rtol=1e-12 atol=1e-12
        @test receipt.target_energy_moment_rates[3-g,k] ≈ sum(w .* moments[1,g,:,k]) rtol=1e-12
    end
    erased = deepcopy(problem.sn); Radiant.set_legendre_order(erased,0); Radiant.set_angular_boltzmann(erased,"standard")
    @test_throws ErrorException project_volume_source(source,problem.cs,problem.geo,erased)
    incomplete = Energy_Moment_Volume_Source(public_energy_base(problem,moments[:,:,:,1:3];
        representation=:ordinates,directions=dirs,weights=w),moments[:,:,:,1:3])
    @test_throws ErrorException project_volume_source(incomplete,problem.cs,problem.geo,problem.sn)
    bte = deepcopy(problem.sn); Radiant.set_solver_type(bte,"BTE")
    @test_throws ErrorException project_volume_source(source,problem.cs,problem.geo,bte)
    mutated = deepcopy(source); mutated.energy_moments[1,1,1,1] += 1
    @test_throws ErrorException project_volume_source(mutated,problem.cs,problem.geo,problem.sn)
    provenance = Dict("angular_basis"=>"radiant-volume-moments/v1",
        "zeroth_moment_index"=>"1","zeroth_moment_is_angle_integrated"=>"true")
    angular = zeros(2,2,size(Dn,1),4)
    for v in 1:2,g in 1:2,k in 1:4; angular[v,g,:,k] .= Dn*moments[v,g,:,k]; end
    native = Energy_Moment_Volume_Source(public_energy_base(problem,angular;
        representation=:moments,provenance=provenance),angular)
    native_projected,_ = project_volume_source(native,problem.cs,problem.geo,problem.sn)
    @test native_projected ≈ projected rtol=1e-12 atol=1e-12
    forbidden = copy(angular)
    for v in 1:2,g in 1:2
        forbidden[v,g,:,1] .= Dn*[-0.1,1.1]; forbidden[v,g,:,2:4] .= 0
    end
    badnative = Energy_Moment_Volume_Source(public_energy_base(problem,forbidden;
        representation=:moments,provenance=provenance),forbidden)
    @test_throws ErrorException project_volume_source(badnative,problem.cs,problem.geo,problem.sn)
    misaligned = deepcopy(source); misaligned.source.energy_edges_eV[2] += 1.0
    @test_throws ErrorException project_volume_source(misaligned,problem.cs,problem.geo,problem.sn)
    reversed = deepcopy(source); reverse!(reversed.source.energy_edges_eV)
    @test_throws ErrorException project_volume_source(reversed,problem.cs,problem.geo,problem.sn)
    malformed_edges = deepcopy(source); malformed_edges.source.energy_edges_eV[2] = NaN
    @test_throws ErrorException project_volume_source(malformed_edges,problem.cs,problem.geo,problem.sn)
    incoming = zeros(1,2,length(w))
    for d in eachindex(w); dirs[d,1] > 0 && (incoming[:,:,d] .= 1); end
    args = (problem.p,[1],reshape([-0.7,0.,0.],1,3),[1.0],reshape([-1.,0.,0.],1,3),
        reshape([0.,1.,0.],1,3),reshape([0.,0.,-1.],1,3),[1e6,5e6,10e6],dirs,w,incoming,base.normalization)
    @test Boundary_Angular_Current_Source(args...) isa Boundary_Angular_Current_Source
    @test_throws ErrorException Boundary_Angular_Current_Source(args...;energy_moments=zeros(1,2,length(w),4))
end

@testset "Public energy source bounded native transport and strength additivity" begin
    problem = public_energy_problem()
    coefficients = zeros(2,2,1,4)
    for v in 1:2,g in 1:2; coefficients[v,g,1,:] .= 1e-6 .* [1.0,-0.2,0.3,0.1]; end
    moments = source_energy_moments([1e6,5e6,10e6],coefficients)
    maker(scale) = Energy_Moment_Volume_Source(public_energy_base(problem,scale*moments),scale*moments)
    currents = Any[]
    for strengths in ((1,),(3,),(1,3))
        fs = Fixed_Sources(problem.cs,problem.geo,problem.solvers)
        for strength in strengths; Radiant.add_source(fs,maker(strength)); end
        Radiant.build(fs)
        flux = Radiant.transport(problem.cs,problem.geo,problem.solvers,fs;retain_boundary_flux=true)
        data = only(Radiant.get_boundary_flux(flux,problem.p))
        @test !isempty(data.convergence)
        @test all(r -> r.converged && r.iterations <= r.iteration_cap &&
            isfinite(r.solver_residual) && r.solver_residual <= r.tolerance &&
            isfinite(r.reconstruction_residual) && r.reconstruction_residual <= r.tolerance,data.convergence)
        current = Radiant.get_outgoing_current(flux,problem.p)
        @test all(isfinite,current)
        push!(currents,current)
        println("PUBLIC_ENERGY_SOURCE strengths=",strengths," outgoing=",current,
            " convergence=",data.convergence," physical_validation=false")
    end
    @test currents[2] ≈ 3currents[1] rtol=1e-9 atol=1e-10
    @test currents[3] ≈ 4currents[1] rtol=1e-9 atol=1e-10
end
