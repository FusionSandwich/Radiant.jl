using Test
using Radiant

@testset "Alpha electronic stopping foundations" begin
    me = 0.51099895069
    re = 2.81794092e-13
    fine = 1 / 137.035999177
    alpha = Alpha()
    @test Radiant.get_type(alpha) === Alpha
    @test Radiant.get_charge(alpha) == 2
    @test Radiant.get_mass(alpha) == 3727.3794118

    # An independent SI/MeV kinematic calculation detects using proton mass or
    # a nonrelativistic maximum transfer. The cache uses electron-rest-energy units.
    for p in (alpha, Proton()), T in (5.0, 200.0, 10000.0)
        M = Radiant.get_mass(p)
        gamma = 1 + T / M
        b2 = 1 - inv(gamma^2)
        tmax = 2 * me * b2 * gamma^2 / (1 + 2 * gamma * me / M + (me/M)^2)
        cache = Radiant.HeavyInelasticCache(1.0, T/me, p)
        @test cache.W_ridge * me ≈ tmax rtol=1e-11
        @test cache.W_ridge * me <= T
    end

    # Equal beta*gamma must give the same density correction for every projectile
    # mass. The old momentum/me expression fails this alpha/electron comparison.
    Ne, I = 3.0e23, 78e-6/me
    for bg in (0.1, 20.0, 2000.0)
        g = sqrt(1 + bg^2)
        d = Float64[]
        for p in (Electron(), Proton(), alpha)
            massratio = Radiant.get_mass(p) / me
            E = (g - 1) * massratio
            push!(d, Radiant.fermi_density_effect_from_electrons([1], [1.0], Ne,
                E, "solid", "sternheimer", massratio, I))
        end
        @test d[1] ≈ d[2] rtol=1e-11 atol=1e-12
        @test d[1] ≈ d[3] rtol=1e-11 atol=1e-12
        if bg == 2000.0
            omega = sqrt(4*pi*Ne*(3.86189603e-11)^2*re)
            @test d[1] ≈ 2*log(bg) - 2*log(I/omega) - 1 rtol=1e-12
        elseif bg == 0.1
            @test d[3] == 0.0
        end
    end

    # Independent unscreened Bethe bracket, with all correction switches off.
    # Changing isotope masses changes electron density, even at identical Z/w/rho.
    Z, w, rho = [1, 8], [0.2, 0.8], 1.2
    for A in ([1.0, 16.0], [2.0, 18.0]), p in (alpha, Proton(), Electron())
        E = 200.0 / me
        M = Radiant.get_mass(p) / me
        gamma = 1 + E/M
        b2 = 1 - inv(gamma^2)
        q = Radiant.get_charge(p)
        if p.type in (Alpha, Proton)
            q *= -expm1(-sqrt(b2)/(fine*q))
        end
        ne = rho * Radiant.Nₐ * sum(w .* Z ./ A)
        oracle = 4*pi*re^2/b2*q^2*ne*(log(2*b2*gamma^2/I)-b2)
        actual = Radiant.bethe(Z,w,rho,E,p,"none","solid",I;
            atomic_weights=A,is_bethe_correction=false,is_density_correction=false,
            is_shell_correction=false,is_lindhard_sorensen=false,is_barkas=false)
        @test actual ≈ oracle rtol=1e-11
        atoms = w ./ A
        atz = atoms ./ sum(atoms)
        N = rho * Radiant.Nₐ * sum(atoms)
        @test Radiant.bethe(Z,w,atz,rho,N,E,p,"none","solid",I;
            atomic_weights=A,is_bethe_correction=false,is_density_correction=false,
            is_shell_correction=false,is_lindhard_sorensen=false,is_barkas=false) ≈ oracle rtol=1e-11
    end

    # Hydrogen has one oscillator, W=I. Its closed-form Fano root gives an
    # independent oracle without reusing the implementation's root finder.
    gamma = 100.0
    omega2 = 4*pi*Ne*(3.86189603e-11)^2*re
    L2 = gamma^2*omega2 - I^2
    @test L2 > 0
    expected_fano = log1p(L2/I^2) - L2/(gamma^2*omega2)
    @test Radiant.fermi_density_effect_from_electrons([1],[1.0],Ne,
        gamma-1,"solid","fano",1.0,I) ≈ expected_fano rtol=1e-5

    # Duplicating an elemental row must preserve oscillator normalization.
    b1 = Radiant.barkas_correction([1],[1.0],Ne,2.0,alpha,I)
    b2 = Radiant.barkas_correction([1,1],[0.25,0.75],Ne,2.0,alpha,I)
    @test b1 > 0
    @test b2 ≈ b1 rtol=1e-13
    @test Radiant.barkas_correction([1],[1.0],Ne,2.0,alpha,2*I) != b1

    # Independent mixed-element oscillator sum. The old nested loop multiplied
    # the entire compound oscillator sum by sum(w/Z), introducing cross terms.
    za, aa = [1,8], [0.8,0.2]
    gamma = 2.0
    beta2 = 1 - inv(gamma^2)
    omega = sqrt(omega2)
    zbar = sum(aa .* za)
    cb = max(1.0,zbar/10)
    oscillator_sum = 0.0
    for j in eachindex(za)
        _,population,_,_,_,_ = Radiant.electron_subshells(za[j])
        frequencies = Radiant.resonance_energy(za[j],omega,I)
        for k in eachindex(population)
            xi = 0.5616*cb*frequencies[k]/(gamma*beta2)
            oscillator_sum += aa[j]/zbar*population[k]*frequencies[k]*
                (Radiant.I₁_barkas(xi)+Radiant.I₂_barkas(xi)/gamma^2)
        end
    end
    expected_barkas = 2*fine/(gamma^2*beta2^1.5)*oscillator_sum
    @test Radiant.barkas_correction(za,aa,Ne,gamma,alpha,I) ≈ expected_barkas rtol=1e-12

    # Analytic joins at Emax and Ecut, with numeric slopes from distinct sides.
    bare(E) = Radiant.bethe([1],[1.0],1.0,E,alpha,"none","solid",I;
        is_bethe_correction=false,is_density_correction=false,is_shell_correction=false,
        is_lindhard_sorensen=false,is_barkas=false)
    Ec, Em = 5/me, 0.8/me
    @test bare(0.0) == 0.0
    for join in (Ec,)
        step = join*1e-5
        left = (bare(join)-bare(join-step))/step
        right = (bare(join+step)-bare(join))/step
        @test bare(join-step) ≈ bare(join+step) rtol=2e-4
        @test left ≈ right rtol=3e-3 atol=bare(join)/join*1e-4
    end
    # The upper branch at Em has log(E/Em)^(3/2), hence a continuous zero
    # derivative with only O(sqrt(h)) convergence of a one-sided secant.
    # Checking two finite secants for near equality at one h is not a C1 test.
    errors = Float64[]
    for epsilon in (1e-3,1e-5,1e-7)
        step = Em*epsilon
        left = (bare(Em)-bare(Em-step))/step
        right = (bare(Em+step)-bare(Em))/step
        push!(errors,max(abs(left),abs(right))/(bare(Em)/Em))
        @test left >= 0 && right <= 0
    end
    @test errors[2] < 0.12*errors[1]
    @test errors[3] < 0.12*errors[2]
    @test errors[3] < 1e-3
    @test bare(Em) > bare(0.9Em)
    @test bare(Em) > bare(1.1Em)

    interaction = Inelastic_Collision()
    @test Radiant.get_ecut(interaction,Alpha) == 5.0
    @test Radiant.get_emax(interaction,Alpha) == 0.8
    @test Radiant.get_low_energy_model(interaction) == "analytic"
    @test !Radiant.get_is_bragg_additivity(interaction)
    Radiant.set_emax(interaction,Alpha,0.1)
    Radiant.set_ecut(interaction,Alpha,0.6)
    @test Radiant.get_emax(interaction,Alpha) == 0.1
    @test Radiant.get_ecut(interaction,Alpha) == 0.6
    # A supported positive custom model, and a lower cutoff rejected for its
    # increasing Bethe slope (DomainError), not a lost maximum (ArgumentError).
    custom = Radiant.bethe([1],[1.0],1.0,0.3/me,alpha,"none","solid",I;
        user_ecut=4.0,user_emax=0.6,is_bethe_correction=false,
        is_density_correction=false,is_shell_correction=false,
        is_lindhard_sorensen=false,is_barkas=false)
    @test isfinite(custom) && custom > 0
    @test_throws DomainError Radiant.bethe([1],[1.0],1.0,0.3/me,alpha,"none","solid",I;
        user_ecut=0.6,user_emax=0.1,is_bethe_correction=false,
        is_density_correction=false,is_shell_correction=false,
        is_lindhard_sorensen=false,is_barkas=false)
    interaction = Inelastic_Collision()
    for bad in (NaN,Inf,-Inf,0.0,-1.0)
        @test_throws ArgumentError Radiant.set_ecut(interaction,Alpha,bad)
        @test_throws ArgumentError Radiant.set_emax(interaction,Alpha,bad)
    end
    @test_throws ArgumentError Radiant.set_emax(interaction,Alpha,5.0)
    @test_throws ArgumentError Radiant.set_ecut(interaction,Alpha,0.5)
    for bad in (NaN,Inf,-Inf,-1.0)
        @test_throws ArgumentError bare(bad)
    end
    @test_throws ArgumentError Radiant.bethe([1],[1.0],1.0,20.0,alpha;
        user_ecut=Inf)
    @test_throws ArgumentError Radiant.bethe([1],[1.0],1.0,20.0,alpha;
        user_emax=5.0)
    # Deliberately inconsistent model partition must be rejected, never clipped.
    @test_throws DomainError Radiant.checked_soft_stopping(1.0,2.0)
    @test_throws DomainError Radiant.checked_soft_stopping(Inf,0.0)
    @test Radiant.checked_soft_stopping(3.0,2.0) == 1.0
    @test_throws ArgumentError Radiant.bethe([1],[1.0],[1.0],1.0,1.0,20.0,alpha)
    @test_throws ArgumentError Radiant.bethe([1],[1.0],1.0,20.0,alpha,"fano","solid",Inf)
    @test_throws ArgumentError Radiant.bethe([1],[1.0],1.0,20.0,Electron();low_energy_model="icru49")
    @test_throws ArgumentError Radiant.bethe([1],[1.0],1.0,0.1/me,alpha;low_energy_model="icru49")
end
