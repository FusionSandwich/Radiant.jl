using Test

# Deliberately tiny independent constant-source transport controls. Four diagonal
# directions have equal |cosine|, so reflection/periodic half-range averages are exact.
function boundary_balance_fixture(dimension;acceleration="none",cap=200,bc=zeros(Int64,2dimension))
    μ = inv(sqrt(Float64(dimension)))
    signs = vec(collect(Iterators.product(ntuple(_ -> (-1.0,1.0),dimension)...)))
    Ω = [Float64[a <= dimension ? s[a]*μ : 0.0 for s in signs] for a in 1:3]
    Nd = length(signs); w = ones(Nd)
    Ns = Int64[a <= dimension ? 2 : 1 for a in 1:3]
    widths = [a <= dimension ? [0.3,0.7] : [1.0] for a in 1:3]
    orders = ones(Int64,4); moments = ones(Int64,5)
    ω,C,is_adaptive,W = Radiant.scheme_weights(orders,fill("DD",4),Int64(dimension),false)
    Mn = ones(Nd,1); Dn = reshape(w,1,Nd)
    Mn_surf = Array{Float64}[]; Dn_surf = Array{Float64}[]; mapping = Vector{Int64}[]
    for axis in 1:dimension, side in 1:2
        # Surface half-range bases describe incoming directions; sweep uses the
        # same coefficient for its opposite outgoing face.
        ids = findall(side == 1 ? Ω[axis] .> 0 : Ω[axis] .< 0)
        push!(Mn_surf,ones(length(ids),1))
        push!(Dn_surf,fill(1.0/length(ids),1,length(ids)))
        map = zeros(Int64,Nd); map[ids] .= 1:length(ids); push!(mapping,map)
    end
    sources = Array{Union{Array{Float64},Float64}}(undef,1,2dimension)
    for f in 1:2dimension
        tangent = [a for a in 1:dimension if a != (f+1)÷2]
        sources[1,f] = isempty(tangent) ? 0.0 : zeros((Ns[a] for a in tangent)...)
    end
    flux = zeros(1,1,Ns...); Q = ones(size(flux))
    data = SN_Boundary_Flux(dimension,Ω,w,widths,[2.0,1.0],bc,1,moments)
    Radiant.sn_one_speed(flux,Q,[0.0],zeros(1,1),ones(Int64,Ns...),Int64(dimension),Int64(Nd),Int64(1),Ns,widths,Ω,Mn,Dn,Int64(1),Int64[0],Mn_surf,Dn_surf,Int64(1),mapping,orders,moments,false,C,ω,Int64(cap),1e-10,sources,is_adaptive,false,Int64(1),0.0,zeros(0),Float64[],Float64[],zeros(0),Float64[],zeros(0),acceleration,Int64(0),false,zeros(1,1),W,bc,Int64(1);boundary_flux=data)
    return data,Float64(Nd) # uniform Q=1, unit domain volume, weights all 1
end

@testset "Actual SN outgoing face current" begin
    for dimension in 1:3
        data,source = boundary_balance_fixture(dimension)
        @test sum(get_outgoing_current(data)) ≈ source rtol=1e-9
        @test data.convergence[1].converged
        @test data.convergence[1].reconstruction_residual < 1e-9
        @test length(data.faces) == 2dimension
    end
    # With an open x+ path a reflective x- face cannot trap source indefinitely.
    reflected,source = boundary_balance_fixture(1;bc=Int64[1,0])
    @test get_outgoing_current(reflected)[2,1] ≈ source rtol=1e-9
    @test get_outgoing_current(reflected)[1,1] > 0 # crossing, not escaped leakage
    # Transverse reflection and periodicity have an open longitudinal escape path.
    for transverse in (1,2)
        data,source = boundary_balance_fixture(2;bc=Int64[0,0,transverse,transverse])
        current = get_outgoing_current(data)
        @test sum(current[1:2,:]) ≈ source rtol=1e-8
        @test current[3,1] ≈ current[4,1] rtol=1e-8
        @test data.convergence[1].converged
    end
    # Krylov/Anderson scratch sweeps must not be mistaken for physical final faces.
    for acceleration in ("livolant","anderson","gmres","bicgstab")
        data,source = boundary_balance_fixture(1;acceleration=acceleration,bc=Int64[1,0])
        @test get_outgoing_current(data)[2,1] ≈ source rtol=1e-8
        @test data.convergence[1].converged
        @test data.convergence[1].reconstruction_residual < 1e-8
    end
    capped,_ = boundary_balance_fixture(1;cap=1,bc=Int64[1,0])
    @test !capped.convergence[1].converged
    @test capped.convergence[1].terminal == :solver_not_converged
    @test capped.convergence[1].reconstruction_residual > capped.convergence[1].tolerance
end

@testset "Boundary normal, quadrature and geometry discriminators" begin
    Ω = [[-0.5,0.25],[0.0,0.0],[0.0,0.0]]
    data = SN_Boundary_Flux(3,Ω,[2.0,3.0],[[7.0],[2.0,3.0],[5.0,11.0]],[4.0,2.0,1.0],zeros(Int64,6),2,ones(Int64,5))
    data.faces[1][1,1,1,:,:] .= [1.0 2.0;3.0 4.0]
    data.faces[1][1,2,1,:,:] .= 999.0 # incoming must not contribute on x-
    data.faces[2][2,2,1,:,:] .= 2.0
    current = get_outgoing_current(data)
    @test current[1,1] == 2.0*0.5*(1*2*5+2*2*11+3*3*5+4*3*11)
    @test current[2,2] == 3.0*0.25*2.0*(2+3)*(5+11)
    @test current[1,2] == 0.0
    # Counterexamples remain visible: replacing normals, weights or areas changes score.
    wrong_normal = deepcopy(data); wrong_normal.directions[1] .*= -1
    wrong_weight = deepcopy(data); wrong_weight.weights .*= 2
    wrong_area = deepcopy(data); wrong_area.widths[2] .*= 2
    @test get_outgoing_current(wrong_normal) != current
    @test get_outgoing_current(wrong_weight) != current
    @test get_outgoing_current(wrong_area) != current
    data.faces[1][1,1,1,1,1] = -1.0
    @test data.faces[1][1,1,1,1,1] < 0 # no clipping of adverse values
    @test_throws ErrorException SN_Boundary_Flux(4,Ω,[2.0,3.0],data.widths,[2.0,1.0],zeros(Int64,8),1,ones(Int64,5))
    @test_throws ErrorException SN_Boundary_Flux(1,Ω,[0.0,3.0],data.widths,[2.0,1.0],zeros(Int64,2),1,ones(Int64,5))
    particle = Radiant.Flux_Per_Particle(Photon())
    Radiant.add_flux(particle,zeros(2,1,1,1,1,1))
    @test_throws ErrorException get_outgoing_current(particle)
    push!(particle.boundary_flux,data)
    @test get_outgoing_current(particle) == get_outgoing_current(data)
    second = Radiant.Flux_Per_Particle(particle.particle)
    Radiant.add_flux(second,zeros(2,1,1,1,1,1)); push!(second.boundary_flux,deepcopy(data))
    Radiant.add_flux(particle,second)
    @test get_outgoing_current(particle) == 2*get_outgoing_current(data)
    # Native 1D Geometry returns a three-slot vector with inactive axes undefined.
    # Both the first generation and a subsequent compatibility check must be safe.
    native_widths = Vector{Vector{Float64}}(undef,3)
    native_widths[1] = [0.3,0.7]
    native = SN_Boundary_Flux(1,Ω,[2.0,3.0],native_widths,[2.0,1.0],zeros(Int64,2),1,ones(Int64,5))
    @test native.widths == [[0.3,0.7],Float64[],Float64[]]
    native.faces[1][1,1,1,1,1] = 4.0
    native_particle = Radiant.Flux_Per_Particle(Photon())
    Radiant.add_flux(native_particle,zeros(1,1,1,2,1,1)); push!(native_particle.boundary_flux,native)
    @test get_outgoing_current(native_particle) == reshape([4.0,0.0],2,1)
    native_second = Radiant.Flux_Per_Particle(native_particle.particle)
    Radiant.add_flux(native_second,zeros(1,1,1,2,1,1)); push!(native_second.boundary_flux,deepcopy(native))
    Radiant.add_flux(native_particle,native_second)
    @test get_outgoing_current(native_particle) == reshape([8.0,0.0],2,1)
    native_second.boundary_flux[1].directions[1][1] = -0.75
    @test_throws ErrorException get_outgoing_current(native_particle)
end
