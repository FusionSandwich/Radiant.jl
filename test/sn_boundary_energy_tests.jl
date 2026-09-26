using Test

# Independent four-point Gauss-Legendre rule on [-1,1]. The reconstruction and
# tuple indexing here are deliberately local to these tests.
const _energy_gl_x = [-0.8611363115940526, -0.3399810435848563,
                       0.3399810435848563,  0.8611363115940526]
const _energy_gl_w = [0.3478548451374538, 0.6521451548625461,
                      0.6521451548625461, 0.3478548451374538]
_energy_legendre(n, x) = n == 0 ? 1.0 : n == 1 ? sqrt(3.0)*x :
    n == 2 ? sqrt(5.0)*(3x^2-1)/2 : error("test rule only needs degrees 0–2")

function _energy_data(dim, orders, fc, bc=zeros(Int64,2dim))
    Ω = [[-0.8, 0.35], [0.0, 0.0], [0.0, 0.0]]
    # For each dimension, both signs of the active normal cosine occur; the
    # deliberately unequal weights exercise ordinate weighting.
    Ω[1] = [-0.8, 0.35]
    dim >= 2 && (Ω[2] = [0.2, -0.55])
    dim >= 3 && (Ω[3] = [-0.4, 0.7])
    widths = [[2.0, 3.0], [1.5, 4.0], [2.5, 5.0]]
    nm = Int64[]
    for a in 1:dim
        tangents = [b for b in 1:dim if b != a]
        push!(nm, fc ? prod((orders[b] for b in tangents); init=1)*orders[4] :
            1 + (orders[4]-1) + sum((orders[b]-1 for b in tangents); init=0))
    end
    SN_Boundary_Flux(dim, Ω, [1.25, 2.5], widths, [5.0, 1.0], bc, 1,
        nm; orders=Int64.(orders),
        fully_coupled=fc, is_csd=true)
end

function _face_index(e, a, b, oe, oa, fc)
    fc && return e + oe*(a-1) + oe*oa*(b-1)
    e > 1 && a > 1 && return 0
    e > 1 && b > 1 && return 0
    a > 1 && b > 1 && return 0
    e > 1 && return e
    a > 1 && return oe+a-1
    b > 1 && return oe+(oa-1)+b-1
    return 1
end

function _independent_energy_integral(coeff, orders, fc, eupper, elower)
    # Integrate the represented polynomial in all active coordinates, without
    # using any Radiant indexing or contraction routine.
    oe, oa, ob = orders
    width = eupper-elower
    total = 0.0
    for ie in eachindex(_energy_gl_x), ia in eachindex(_energy_gl_x), ib in eachindex(_energy_gl_x)
        ue, ua, ub = _energy_gl_x[ie], _energy_gl_x[ia], _energy_gl_x[ib]
        psi = 0.0
        for e in 1:oe, a in 1:oa, b in 1:ob
            idx = _face_index(e,a,b,oe,oa,fc)
            idx == 0 && continue
            psi += coeff[idx]*_energy_legendre(e-1,ue)*
                   _energy_legendre(a-1,ua)*_energy_legendre(b-1,ub)
        end
        E = (eupper+elower)/2 + width*ue/2
        total += _energy_gl_w[ie]*_energy_gl_w[ia]*_energy_gl_w[ib]*E*psi/8
    end
    return total
end

@testset "Represented energy polynomial counterexamples" begin
    for slope in (sqrt(3.0)/2, -sqrt(3.0)/2)
        d = _energy_data(1, [1,1,1,3], true)
        d.directions[1] .= 1.0; d.weights .= 1.0
        d.faces[2][1,1,1,1,1] = 2.0
        d.faces[2][1,1,2,1,1] = slope
        d.faces[2][1,1,3,1,1] = 91.0 # energy degree 2 integrates to zero
        @test get_outgoing_energy_current(d)[2,1] ≈ (slope > 0 ? 7.0 : 5.0)
    end
    # OE=1 means face moment 2 is a tangent mode, never an energy slope.
    d = _energy_data(2, [1,2,1,1], true)
    d.directions[1] .= 1.0; d.weights .= 1.0
    d.faces[2][1,1,1,1,1] = 2.0
    d.faces[2][1,1,2,1,1] = 37.0
    @test get_outgoing_energy_current(d)[2,1] ≈ 6.0*1.5
end

@testset "Independent coupled and reduced face polynomial quadrature" begin
    # Unequal orders and distinguishable tuple coefficients expose wrong energy
    # axes, slope signs, flattening, and accidental inclusion of mixed modes.
    spatial_orders = [2,3,2]
    for dim in 1:3, fc in (true,false), axis in 1:dim
        orders = [spatial_orders; 2]
        orders[dim+1:3] .= 1
        d = _energy_data(dim, orders, fc)
        f = 2axis
        nf = size(d.faces[f],2)
        tangent_axes = [a for a in 1:dim if a != axis]
        ncells = size(d.faces[f],4)*size(d.faces[f],5)
        for n in 1:nf, c in 1:ncells
            coeff = zeros(size(d.faces[f],3))
            # Constant and pure energy slope have independent signatures.
            coeff[1] = 0.8 + 0.3n + 0.07c
            ei = _face_index(2,1,1,orders[4],
                isempty(tangent_axes) ? 1 : orders[tangent_axes[1]],fc)
            orders[4] >= 2 && (coeff[ei] = (-1)^n*(0.17+0.02c))
            # Populate every representable higher/tangent mode so quadrature
            # confirms its zero integral, including large and negative values.
            for i in 2:length(coeff)
                i == ei && continue
                coeff[i] = (-1)^i*(0.11i + 0.03n + 0.01c)
            end
            j = mod1(c,size(d.faces[f],4)); k = div(c-1,size(d.faces[f],4))+1
            d.faces[f][1,n,:,j,k] .= coeff
            expected = _independent_energy_integral(coeff,
                (orders[4], isempty(tangent_axes) ? 1 : orders[tangent_axes[1]],
                 length(tangent_axes)<2 ? 1 : orders[tangent_axes[2]]),fc,5.0,1.0)
            area = isempty(tangent_axes) ? 1.0 : d.widths[tangent_axes[1]][j]
            length(tangent_axes)==2 && (area *= d.widths[tangent_axes[2]][k])
            outward = max(0.0, (isodd(f) ? -1 : 1)*d.directions[axis][n])
            score = get_outgoing_energy_current(d)[f,1]
            # Compare one ordinate/cell contribution by subtracting a clean copy
            # with that entry removed; this leaves the independent hand result.
            clean = deepcopy(d); clean.faces[f][1,n,:,j,k] .= 0
            isolated = score-get_outgoing_energy_current(clean)[f,1]
            @test isolated ≈ d.weights[n]*outward*area*expected rtol=2e-11 atol=2e-11
        end
    end
end

@testset "Crossing versus void escaped energy" begin
    labels = Int64[0,1,2,0,1,2]
    d = _energy_data(3,[2,3,2,2],true,labels)
    for f in 1:6, n in 1:2
        d.faces[f][1,n,1,1,1] = 1+n
        d.faces[f][1,n,2,1,1] = 0.1n
    end
    crossing = get_outgoing_energy_current(d)
    escaped = get_escaped_energy_current(d)
    @test size(escaped) == size(crossing)
    @test escaped[1,:] == crossing[1,:]
    @test all(iszero, escaped[2,:]) && all(iszero,escaped[3,:])
    @test escaped[4,:] == crossing[4,:]
    @test all(iszero, escaped[5,:]) && all(iszero,escaped[6,:])
    d.boundary_conditions[2] = 0
    @test get_outgoing_energy_current(d) == crossing
    @test get_escaped_energy_current(d)[2,:] == crossing[2,:]
    d.boundary_conditions[2] = 9
    @test_throws ErrorException get_escaped_energy_current(d)
end

@testset "Energy scorer metadata and fail-closed inputs" begin
    Ω = [[1.0], [0.0], [0.0]]
    legacy = SN_Boundary_Flux(1,Ω,[1.0],[[1.0],[1.0],[1.0]],[5.0,1.0],zeros(Int64,2),1,Int64[2])
    @test_throws ErrorException get_outgoing_energy_current(legacy)
    @test_throws ErrorException get_escaped_energy_current(legacy)
    d = _energy_data(1,[1,1,1,2],true)
    d.faces[2][1,2,1,1,1] = NaN
    @test_throws ErrorException get_outgoing_energy_current(d)
    d.faces[2][1,2,1,1,1] = Inf
    @test_throws ErrorException get_escaped_energy_current(d)
    d.faces[2][1,2,1,1,1] = -2.0
    @test get_outgoing_energy_current(d)[2,1] < 0 # diagnostics retain signed values
    @test_throws ErrorException SN_Boundary_Flux(1,Ω,[1.0],[[1.0],[1.0],[1.0]],
        [5.0,1.0],zeros(Int64,2),1,Int64[1]; orders=Int64[1,1,1,2],
        fully_coupled=true,is_csd=true)
    @test_throws ErrorException SN_Boundary_Flux(1,Ω,[1.0],[[1.0],[1.0],[1.0]],
        [5.0,1.0],zeros(Int64,2),1,Int64[2]; orders=Int64[1,1,1,0],
        fully_coupled=true,is_csd=true)
    d.faces[2][1,2,1,1,1] = floatmax(Float64)
    @test_throws ErrorException get_outgoing_energy_current(d)
    malformed = _energy_data(1,[1,1,1,2],true)
    malformed.faces[2] = zeros(1,1,1,1,1)
    @test_throws ErrorException get_outgoing_energy_current(malformed)
end
