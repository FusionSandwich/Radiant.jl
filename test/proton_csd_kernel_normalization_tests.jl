using Test

# Independent fixed six-node Gauss-Legendre integration, exact through degree11.
# No native polynomial, source-projection, triple-product or scheme helper is
# used to create the expected solution or quadrature tensors.
const _CSD_K_NODES = [-0.9324695142031521,-0.6612093864662645,-0.2386191860831969,
    0.2386191860831969,0.6612093864662645,0.9324695142031521]
const _CSD_K_WEIGHTS = [0.1713244923791704,0.3607615730481386,0.4679139345726910,
    0.4679139345726910,0.3607615730481386,0.1713244923791704]
function _csd_k_polynomials(u,oe)
    p = ones(oe)
    oe >= 2 && (p[2]=u)
    for k in 2:oe-1
        p[k+1]=((2*k-1)*u*p[k]-(k-1)*p[k-1])/k
    end
    return p
end

function _csd_k_call(dimension,oe,coupled,upper,lower,slope;old_factor=false)
    width = upper-lower
    midpoint = (upper+lower)/2
    space_widths = [2.0,3.5,5.0]
    direction = [1.0,-1.0,1.0]./sqrt(3.0)
    h = sum(abs(direction[a])/space_widths[a] for a in 1:dimension)
    phi(E) = upper-E
    stopping(E) = 2.0+slope*(E-lower)
    source(E) = h*phi(E)-slope*phi(E)+stopping(E)
    c = sqrt.(2 .* collect(1:oe) .- 1)
    expected = zeros(oe)
    q = zeros(oe)
    s = zeros(oe)
    w = zeros(oe,oe,oe)
    for (u,weight) in zip(_CSD_K_NODES,_CSD_K_WEIGHTS)
        p = _csd_k_polynomials(u,oe)
        E = midpoint+width*u/2
        expected .+= width*weight/2 .* phi(E) .* c .* p
        q .+= width*weight/2 .* source(E) .* c .* p
        # Native stopping coefficients are normalized energy averages and are
        # divided by width before the cell solve.
        s .+= weight/(2*width) .* stopping(E) .* c .* p
        for j in 1:oe,k in 1:oe,l in 1:oe
            w[j,k,l] += weight*p[j]*p[k]*p[l]/2
        end
    end
    # Stopping is analytically affine. Its higher coefficients vanish exactly;
    # suppressing quadrature roundoff here does not modify a flux observation.
    s[3:end] .= 0.0
    if old_factor
        # With the corrected native C[k]^2 coefficient, W/C[k] exactly
        # recreates the original wrong derivative normalization for each row.
        w ./= reshape(c,1,oe,1)
    end
    lower_trace = width*phi(lower)
    se = stopping(upper)/width
    sl = stopping(lower)/width
    if dimension == 1
        we = zeros(oe+1,1,1); we[2:end,1,1] .= 1.0
        wx = zeros(2,oe,oe)
        for k in 1:oe; wx[2,k,k] = 1.0; end
        result = Radiant.flux_1D_BFP(direction[1],0.0,space_widths[1],q,zeros(oe),
            se,sl,s,width,zeros(1),oe,1,c,we,wx,false,w,coupled)
    elseif dimension == 2
        we = zeros(oe+1,1,1); we[2:end,1,1] .= 1.0
        wx = zeros(2,1,oe); wx[2,1,:] .= 1.0
        wy = copy(wx)
        result = Radiant.flux_2D_BFP(direction[1],direction[2],0.0,se,sl,s,width,
            space_widths[1],space_widths[2],q,zeros(oe),zeros(oe),zeros(1),
            oe,1,1,c,we,wx,wy,false,w,coupled)
    else
        we = zeros(oe+1,1,1,1); we[2:end,1,1,1] .= 1.0
        wx = zeros(2,1,1,oe); wx[2,1,1,:] .= 1.0
        wy,wz = copy(wx),copy(wx)
        result = Radiant.flux_3D_BFP(direction[1],direction[2],direction[3],0.0,se,sl,s,width,
            space_widths[1],space_widths[2],space_widths[3],q,zeros(oe),zeros(oe),zeros(oe),zeros(1),
            oe,1,1,1,c,we,wx,wy,wz,false,w,coupled)
    end
    return (;result,expected,lower_trace)
end

@testset "CSD derivative normalization across native Cartesian kernels" begin
    # Spatial orders are1; coupling flags are exercised but no mixed spatial
    # degree qualification is claimed. Widths and h differ across dimensions.
    for dimension in 1:3, oe in (2,3,4), coupled in (false,true),
        slope in (0.0,0.1), (upper,lower) in ((10.0,1.0),(7.5,2.0))
        f = _csd_k_call(dimension,oe,coupled,upper,lower,slope)
        println("CSD_KERNEL_CASE dimension=",dimension," oe=",oe," coupled=",coupled,
            " slope=",slope," bounds_MeV=",[upper,lower],
            " coefficients=",f.result[1]," expected_coefficients=",f.expected,
            " lower_trace=",f.result[end]," expected_lower_trace=",f.lower_trace)
        @test f.result[1] ≈ f.expected atol=1e-9 rtol=1e-11
        for axis in 1:dimension
            @test f.result[axis+1] ≈ f.expected atol=1e-9 rtol=1e-11
        end
        @test only(f.result[end]) ≈ f.lower_trace atol=1e-9 rtol=1e-11
        wrong = _csd_k_call(dimension,oe,coupled,upper,lower,slope;old_factor=true)
        if oe >= 3
            @test maximum(abs,wrong.result[1]-f.expected) > 1e-4
        else
            @test wrong.result[1] ≈ f.expected atol=1e-9 rtol=1e-11
        end
    end
end
