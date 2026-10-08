using Test
using Radiant

# Independent monomial antiderivatives in 256-bit arithmetic. Float64 endpoints
# are promoted exactly so thin-interval comparisons test the supplied interval.
function hard_moment_monomial_oracle(cache,n,lo,hi)
    setprecision(BigFloat,256) do
        a,b = BigFloat(lo),BigFloat(hi)
        coeff_b = BigFloat(cache.β²)/BigFloat(cache.W_ridge)
        coeff_c = (1-BigFloat(cache.β²))/(2*BigFloat(cache.ratio_mass)^2)
        primitive(power) = power == -1 ? log(b/a) : (b^(power+1)-a^(power+1))/(power+1)
        return Float64(BigFloat(cache.A)*(primitive(n-2)-coeff_b*primitive(n-1)+coeff_c*primitive(n)))
    end
end

# Independent GL integration in log W. This resolves a broad dynamic range
# without concentrating almost all linear-W nodes away from the lower limit.
function hard_moment_log_quadrature(cache,n,lo,hi; nodes=96)
    x,w = Radiant.quadrature(nodes,"gauss-legendre")
    # Integrate t=log(W/lo); log1p retains a thin interval's actual width.
    length_log = log1p((hi-lo)/lo)
    half = length_log/2
    sumvalue = 0.0
    for j in eachindex(x)
        W = lo*exp(half*(x[j]+1))
        # Direct kernel polynomial, independently weighted by dW=W dt.
        polynomial = 1-cache.β²*W/cache.W_ridge+(1-cache.β²)*W^2/(2*cache.ratio_mass^2)
        sumvalue += w[j]*cache.A*W^(n-1)*polynomial
    end
    return half*sumvalue
end

@testset "Hard electronic collision raw moments, n=0/1/2" begin
    me = 0.51099895069
    for particle in (Alpha(),Proton()), energy_MeV in (1.0,5.0,20.0,500.0)
        Ei = energy_MeV/me
        cache = Radiant.HeavyInelasticCache(2.0,Ei,particle)
        m = Radiant.get_mass(particle)
        # Independent invariant transfer expression in MeV; no cache gamma/beta.
        expected_max_MeV = 2*me*energy_MeV*(energy_MeV+2*m)/
            (m*m+2*me*(energy_MeV+m)+me*me)
        @test isapprox(cache.W_ridge*me,expected_max_MeV; rtol=3e-14)
        @test 0 < cache.W_ridge <= Ei
        limit = min(Ei,cache.W_ridge)
        for n in 0:2
            lo,hi = limit*1e-5,limit
            actual = Radiant.heavy_raw_moment(cache,n,lo,hi)
            oracle = hard_moment_monomial_oracle(cache,n,lo,hi)
            numeric = hard_moment_log_quadrature(cache,n,lo,hi)
            @test actual > 0.0
            @test isapprox(actual,oracle; rtol=3e-13)
            @test isapprox(actual,numeric; rtol=3e-12)
            # Cutoff sweep: no negative partition or cancellation is hidden.
            for fraction in (1e-4,0.01,0.1,0.5,0.99,1-1e-8)
                split = limit*fraction
                lower = Radiant.heavy_raw_moment(cache,n,lo,split)
                upper = Radiant.heavy_raw_moment(cache,n,split,hi)
                @test lower >= 0 && upper >= 0
                @test isapprox(lower+upper,actual; rtol=3e-13)
            end
            @test Radiant.heavy_raw_moment(cache,n,limit,2*limit) == 0.0
            @test Radiant.heavy_raw_moment(cache,n,lo,2*limit) == actual
            @test Radiant.heavy_raw_moment(cache,n,hi,hi) == 0.0
            for width in (1e-4,1e-5,1e-6,1e-7,1e-8,1e-9)
                # Near the old branch at 1e-6, its n0 delta/hi² shortcut
                # has O(width) relative error and fails the 1e-9 requirement.
                thin_hi = limit*0.7
                thin_lo = thin_hi*(1-width)
                thin = Radiant.heavy_raw_moment(cache,n,thin_lo,thin_hi)
                exact = hard_moment_monomial_oracle(cache,n,thin_lo,thin_hi)
                @test isapprox(thin,exact; rtol=1e-12)
                @test isapprox(thin,hard_moment_log_quadrature(cache,n,thin_lo,thin_hi); rtol=3e-12)
            end
        end
        # A second raw moment is finite down to zero, unlike number/first moment.
        c = (1-cache.β²)/(2*cache.ratio_mass^2)
        exact_second = cache.A*(limit-cache.β²*limit^2/(2*cache.W_ridge)+c*limit^3/3)
        @test isapprox(Radiant.heavy_raw_moment(cache,2,0.0,limit),exact_second; rtol=3e-14)
        # Differential endpoint and signed regular partition.
        regular,singular = Radiant.heavy_sigma_components(cache,limit*0.5)
        @test regular < 0 < singular
        @test regular+singular > 0
        one = Radiant.HeavyInelasticCache(1.0,Ei,particle)
        @test Radiant.inelastic_collision_heavy_particle(1.0,Ei,limit,particle) > 0
        @test Radiant.inelastic_collision_heavy_particle(1.0,Ei,nextfloat(limit),particle) == 0.0
        for n in 0:2
            threshold = 0.1*limit
            expected = Radiant.heavy_raw_moment(one,n,threshold,limit)
            @test isapprox(Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(
                Ei,n,threshold,particle),expected; rtol=3e-14)
            @test isapprox(Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(
                Ei,n,0.01*limit,particle,threshold),expected; rtol=3e-14)
            @test Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(
                Ei,n,2*limit,particle) == 0.0
        end
        # Feed uses the same primitive and clips outgoing groups at Wmax.
        for width in (1e-4,1e-6,1e-9)
            upper_transfer = 0.7*limit
            lower_transfer = upper_transfer*(1-width)
            Efhi,Eflo = Ei-lower_transfer,Ei-upper_transfer
            realized_lo,realized_hi = Ei-Efhi,Ei-Eflo
            @test isapprox(Radiant.feed_analytical_heavy_particle(cache,Efhi,Eflo),
                hard_moment_monomial_oracle(cache,0,realized_lo,realized_hi); rtol=1e-12)
            @test isapprox(Radiant.feed_first_moment_heavy_particle(cache,Efhi,Eflo),
                hard_moment_monomial_oracle(cache,1,realized_lo,realized_hi); rtol=1e-12)
            expected_singular = cache.A*((realized_hi-realized_lo)/realized_hi)/realized_lo
            @test isapprox(Radiant.integrate_A_over_W2_per_subshell(cache,Efhi,Eflo),expected_singular; rtol=3e-14)
        end
        Efhi,Eflo = Ei-0.2*limit,Ei-2*limit
        @test isapprox(Radiant.feed_analytical_heavy_particle(cache,Efhi,Eflo),
            Radiant.heavy_raw_moment(cache,0,Ei-Efhi,limit); rtol=3e-14)
        @test Radiant.integrate_leading_1overW_per_subshell(cache,Efhi,Eflo,1) < 0
    end
end

@testset "Hard collision moments reject invalid and divergent domains" begin
    particle = Alpha()
    cache = Radiant.HeavyInelasticCache(1.0,20.0/0.51099895069,particle)
    hi = cache.W_ridge
    for n in (-1,3,8)
        @test_throws ArgumentError Radiant.heavy_raw_moment(cache,n,hi,hi)
        @test_throws ArgumentError Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(
            cache.Ei,n,2*hi,particle)
        @test_throws ArgumentError Radiant.integrate_inelastic_collision_heavy_particle(1,cache.Ei,n,particle)
    end
    for n in 0:2, (lo,upper) in ((-1.0,hi),(hi,0.0),(NaN,hi),(0.1*hi,Inf))
        @test_throws ArgumentError Radiant.heavy_raw_moment(cache,n,lo,upper)
    end
    for n in 0:1
        @test_throws DomainError Radiant.heavy_raw_moment(cache,n,0.0,hi)
        @test_throws DomainError Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(cache.Ei,n,0.0,particle)
    end
    for E in (0.0,-1.0,NaN,Inf)
        @test_throws ArgumentError Radiant.HeavyInelasticCache(1.0,E,particle)
    end
    for Zi in (-1.0,NaN,Inf)
        @test_throws ArgumentError Radiant.HeavyInelasticCache(Zi,cache.Ei,particle)
    end
    for unsupported in (Electron(),Muon(),Particle("charged fixture",1000.0,1.0))
        @test_throws ArgumentError Radiant.HeavyInelasticCache(1.0,cache.Ei,unsupported)
    end
    for invalid in (-1.0,NaN,Inf)
        @test_throws ArgumentError Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(cache.Ei,0,invalid,particle)
        @test_throws ArgumentError Radiant.integrate_inelastic_collision_heavy_particle_per_subshell(cache.Ei,0,0.1*hi,particle,invalid)
        @test_throws ArgumentError Radiant.heavy_sigma_components(cache,invalid)
    end
    @test_throws ArgumentError Radiant.heavy_sigma_components(cache,0.0)
    for f in (Radiant.feed_analytical_heavy_particle,Radiant.feed_first_moment_heavy_particle,
              Radiant.integrate_A_over_W2_per_subshell)
        @test_throws DomainError f(cache,cache.Ei,cache.Ei-0.1*hi)
        @test_throws ArgumentError f(cache,cache.Ei+1,cache.Ei-0.1*hi)
        @test_throws ArgumentError f(cache,cache.Ei-0.2*hi,cache.Ei-0.1*hi)
        @test_throws ArgumentError f(cache,NaN,0.0)
        @test f(cache,cache.Ei-0.1*hi,cache.Ei-0.1*hi) == 0.0
    end
end
