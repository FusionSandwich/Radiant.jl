using Test
using LinearAlgebra

# GN has angular coefficient blocks, not ordinate values. In an orthonormal
# angular basis, energy differentiation acts as identity on each angular block.
# Prescribe phi_p(E)=v_p(U-E). The independently integrated source is
# sum_a N_a*v*(U-E)/dx_a + v*(S-S'*(U-E)). Non-diagonal streaming matrices ensure
# this is not a scalar SN fixture reused without the angular-block coupling.
# Spatial polynomial degree is zero; higher spatial degrees remain unqualified.
function _gn_csd_case(dim,oe,coupled,upper,lower,slope;old_factor=false)
    width = upper-lower
    midpoint = (upper+lower)/2
    widths = [2.0,3.5,5.0]
    v = [1.0,0.35]
    matrices = [[0.6 0.08;0.1 0.5], [0.5 -0.04;0.03 0.4], [0.7 0.02;-0.05 0.6]]
    hv = sum(matrices[a]*v/widths[a] for a in 1:dim)
    c = sqrt.(2 .* collect(1:oe) .- 1)
    phi(E) = upper-E
    stopping(E) = 2.0+slope*(E-lower)
    expected, qn = zeros(2,oe), zeros(2,oe)
    s = zeros(oe)
    w = zeros(oe,oe,oe)
    for (u,weight) in zip(_CSD_K_NODES,_CSD_K_WEIGHTS)
        p = _csd_k_polynomials(u,oe)
        E = midpoint+width*u/2
        basis = c .* p
        expected .+= width*weight/2 .* (v .* phi(E)) * transpose(basis)
        source = hv .* phi(E) .+ v .* (stopping(E)-slope*phi(E))
        qn .+= width*weight/2 .* source * transpose(basis)
        s .+= weight/(2*width) .* stopping(E) .* basis
        for j in 1:oe,k in 1:oe,l in 1:oe
            w[j,k,l] += weight*p[j]*p[k]*p[l]/2
        end
    end
    s[3:end] .= 0.0
    old_factor && (w ./= reshape(c,1,oe,1))
    se,sl = stopping(upper)/width,stopping(lower)/width
    a = zeros(2*oe,2*oe)
    q,solution = zeros(2*oe),zeros(2*oe)
    fn,fe = zeros(2,oe),zeros(2,1)
    faces = [zeros(2,oe) for _ in 1:dim]
    wx,wy,wz = [0.0,1.0],[0.0,1.0],[0.0,1.0]
    we = vcat(0.0,ones(oe))
    optimized = zeros(2*oe,2*oe)
    if dim == 1
        Radiant.gn_1D_BFP!(fn,faces[1],fe,1,0.0,se,sl,s,widths[1],qn,a,q,solution,1,oe,2,c,wx,we,matrices[1],w,coupled)
        Radiant.gn_1D_BFP_matrix!(optimized,1,0.0,sl,s,widths[1],1,oe,2,c,wx,we,matrices[1],w,coupled)
    elseif dim == 2
        Radiant.gn_2D_BFP!(fn,faces[1],faces[2],fe,1,1,0.0,se,sl,s,widths[1],widths[2],qn,a,q,solution,1,1,oe,2,c,wx,wy,we,matrices[1],matrices[2],w,coupled)
        Radiant.gn_2D_BFP_matrix!(optimized,1,1,0.0,sl,s,widths[1],widths[2],1,1,oe,2,c,wx,wy,we,matrices[1],matrices[2],w,coupled)
    else
        Radiant.gn_3D_BFP!(fn,faces[1],faces[2],faces[3],fe,1,1,1,0.0,se,sl,s,widths[1],widths[2],widths[3],qn,a,q,solution,1,1,1,oe,2,c,wx,wy,wz,we,matrices[1],matrices[2],matrices[3],w,coupled)
        Radiant.gn_3D_BFP_matrix!(optimized,1,1,1,0.0,sl,s,widths[1],widths[2],widths[3],1,1,1,oe,2,c,wx,wy,wz,we,matrices[1],matrices[2],matrices[3],w,coupled)
    end
    optimized_solution = transpose(reshape(optimized \ vec(transpose(qn)),oe,2))
    return (;fn,faces,fe,expected,lower_trace=v .* (width*phi(lower)),optimized_solution)
end

@testset "GN independent angular-block CSD normalization" begin
    for dim in 1:3, oe in (2,3,4), coupled in (false,true), slope in (0.0,0.1),
        (upper,lower) in ((10.0,1.0),(7.5,2.0))
        f = _gn_csd_case(dim,oe,coupled,upper,lower,slope)
        @test f.fn ≈ f.expected atol=1e-9 rtol=1e-11
        @test f.optimized_solution ≈ f.expected atol=1e-9 rtol=1e-11
        for face in f.faces
            @test face ≈ f.expected atol=1e-9 rtol=1e-11
        end
        @test vec(f.fe) ≈ f.lower_trace atol=1e-9 rtol=1e-11
        bad = _gn_csd_case(dim,oe,coupled,upper,lower,slope;old_factor=true)
        if oe >= 3
            @test maximum(abs,bad.fn-f.expected) > 1e-4
        else
            @test bad.fn ≈ f.expected atol=1e-9 rtol=1e-11
        end
    end
end
