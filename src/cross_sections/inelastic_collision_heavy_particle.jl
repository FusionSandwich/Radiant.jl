"""
    HeavyInelasticCache

Stores reusable constants for heavy-particle inelastic collisions with one subshell
electron population at one incident energy.

# Field(s)
- `Zi::Float64` : number of electrons in the subshell.
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `mₑ::Float64` : electron rest energy [MeV].
- `rₑ::Float64` : classical electron radius [cm].
- `Zc::Float64` : incoming particle charge.
- `ratio_mass::Float64` : incoming particle mass divided by the electron rest mass.
- `γ::Float64` : Lorentz factor of the incoming particle.
- `β²::Float64` : squared reduced speed of the incoming particle.
- `R::Float64` : recoil factor used in the maximum energy transfer.
- `W_ridge::Float64` : maximum reduced energy transfer.
- `A::Float64` : leading cross-section coefficient.
- `mu_lin::Float64` : coefficient used in the angular small-energy-transfer expansion.
"""
struct HeavyInelasticCache
    Zi::Float64
    Ei::Float64
    mₑ::Float64
    rₑ::Float64
    Zc::Float64
    ratio_mass::Float64
    γ::Float64
    β²::Float64
    R::Float64
    W_ridge::Float64
    A::Float64  # 2*π*rₑ^2/β² * Zc^2 * Zi
    mu_lin::Float64  # For angular derivatives: 1 / (Ei * (Ei + 2*ratio_mass))
end

"""
    HeavyInelasticCache(Zi::Real, Ei::Float64, particle::Particle)

Builds a cache of heavy-particle inelastic-collision constants.

# Input Argument(s)
- `Zi::Real` : number of electrons in the subshell.
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `particle::Particle` : incoming heavy charged particle.

# Output Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for cross-section and feed integrations.
"""
function HeavyInelasticCache(Zi::Real, Ei::Float64, particle::Particle)
    (is_proton(particle) || is_alpha(particle)) ||
        throw(ArgumentError("The native heavy collision kernel supports only Proton and Alpha."))
    mₑ = 0.51099895069
    rₑ = 2.81794092e-13
    Zc = get_charge(particle)
    ratio_mass = get_mass(particle)/mₑ
    all(isfinite, (Zi, Ei, Zc, ratio_mass)) && Zi >= 0 && Ei > 0 && ratio_mass > 0 && Zc != 0 ||
        throw(ArgumentError("Heavy collision cache requires finite Zi >= 0, Ei > 0, mass > 0 and nonzero charge."))
    tau = Ei / ratio_mass
    γ = 1.0 + tau
    β² = tau * (tau + 2.0) / γ^2
    R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
    W_ridge = 2*tau*(tau + 2.0)*R
    A = 2*π*rₑ^2/β² * Zc^2 * Zi
    mu_lin = 1.0 / (Ei * (Ei + 2.0 * ratio_mass))
    return HeavyInelasticCache(Float64(Zi), Ei, mₑ, rₑ, Zc, ratio_mass, γ, β², R, W_ridge, A, mu_lin)
end

"""
    inelastic_collision_heavy_particle(Zi::Real, Ei::Float64, W::Float64,
    particle::Particle)

Gives the inelastic differential cross-sections with atomic electrons for heavy particles.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `W::Float64` : energy transferred to the atomic electron [in m_e*c^2].
- `particle::Particle` : incoming heavy charged particle.

# Output Argument(s)
- `σs::Float64` : inelastic differential cross-section with atomic electrons for heavy particles.

# Reference(s)
- Salvat and Heredia (2024), Electromagnetic interaction models for Monte-Carlo simulation
  of protons and alpha particles.

"""
function inelastic_collision_heavy_particle(Zi::Real,Ei::Float64,W::Float64,particle::Particle)
    cache = HeavyInelasticCache(Zi, Ei, particle)
    regular, singular = heavy_sigma_components(cache, W)
    return regular + singular
end

"""
    heavy_sigma_components(cache::HeavyInelasticCache, W::Float64)

Splits the heavy-particle inelastic differential cross-section into its singular
`A / W^2` part and the remaining regular part.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `W::Float64` : energy transferred to the atomic electron [in m_e*c^2].

# Output Argument(s)
- `sigma_reg::Float64` : regular contribution, equal to `sigma_full(W) - A / W^2`.
- `A_over_W2::Float64` : leading singular contribution.
"""
function heavy_sigma_components(cache::HeavyInelasticCache, W::Float64)
    isfinite(W) && W > 0.0 || throw(ArgumentError("Heavy collision transfer must be finite and positive."))
    if W > min(cache.Ei, cache.W_ridge)
        return 0.0, 0.0
    end
    A_over_W2 = cache.A / W^2
    # This partition is signed; do not clip the negative regular term.
    sigma_reg = cache.A * (-cache.β²/(cache.W_ridge*W) + (1-cache.β²)/(2*cache.ratio_mass^2))
    return sigma_reg, A_over_W2
end

"""Validate and clip an ordered transfer interval to the hard-kernel support.

Invalid/negative/reversed bounds are errors. An ordered interval above Wmax is
empty. Zero lower transfer is divergent for n=0/1 on a nonempty interval and
must not be replaced by an arbitrary small positive cutoff.
"""
function heavy_transfer_interval(cache::HeavyInelasticCache, lo::Float64, hi::Float64, n::Int)
    n in (0,1,2) || throw(ArgumentError("Heavy collision raw moments support only n=0,1,2."))
    all(isfinite, (lo,hi)) && 0.0 <= lo <= hi ||
        throw(ArgumentError("Heavy collision transfer bounds must be finite, nonnegative and ordered."))
    if lo == 0.0 && hi > 0.0 && n < 2
        throw(DomainError(lo, "The n=0/1 hard collision moment diverges at zero lower transfer."))
    end
    return lo, min(hi, cache.Ei, cache.W_ridge)
end

function heavy_interval_log(lo::Float64, hi::Float64)
    ratio = (hi-lo)/lo
    return isfinite(ratio) ? log1p(ratio) : log(hi)-log(lo)
end

"""Exact finite hard-kernel raw transfer moment ∫ W^n (dσ/dW) dW.

W is measured in electron rest energies; units are cm²*(m_e*c²)^n.
n=2 is a raw hard collision moment, not a central moment or soft straggling
coefficient. No diffusion factor 1/2 is included.
"""
function heavy_raw_moment(cache::HeavyInelasticCache, n::Int, lo::Float64, hi::Float64)
    lo, hi = heavy_transfer_interval(cache,lo,hi,n)
    hi <= lo && return 0.0
    delta = hi-lo
    b = cache.β²/cache.W_ridge
    c = (1-cache.β²)/(2*cache.ratio_mass^2)
    if n == 0
        # delta/(lo*hi) exactly, without the thin-interval delta/hi² shortcut.
        value = (delta/hi)/lo - b*heavy_interval_log(lo,hi) + c*delta
    elseif n == 1
        value = heavy_interval_log(lo,hi) - b*delta + c*delta*(hi+lo)/2
    else
        value = delta * (1 - b*(hi+lo)/2 + c*(hi^2+hi*lo+lo^2)/3)
    end
    return cache.A*value
end

function heavy_feed_interval(cache::HeavyInelasticCache, Ef_minus::Float64, Ef_plus::Float64, n::Int)
    all(isfinite, (Ef_minus,Ef_plus)) && 0.0 <= Ef_plus <= Ef_minus <= cache.Ei ||
        throw(ArgumentError("Heavy collision outgoing bounds must be finite, ordered and in [0,Ei]."))
    return heavy_transfer_interval(cache,cache.Ei-Ef_minus,cache.Ei-Ef_plus,n)
end

"""
    integrate_A_over_W2_per_subshell(cache::HeavyInelasticCache,
    Ef_minus::Float64, Ef_plus::Float64)

Integrates the leading `A / W^2` part analytically over an outgoing-energy group.
The integration bounds are converted with `W = Ei - Ef` and clipped to
`0 < W <= W_ridge`.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef_minus::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef_plus::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].

# Output Argument(s)
- `σ::Float64` : integral of `A / W^2` over the clipped group.
"""
function integrate_A_over_W2_per_subshell(cache::HeavyInelasticCache, Ef_minus::Float64, Ef_plus::Float64)
    lo,hi = heavy_feed_interval(cache,Ef_minus,Ef_plus,0)
    hi <= lo && return 0.0
    return cache.A*((hi-lo)/hi)/lo
end

"""
    heavy_leading_1overW_coeff(cache::HeavyInelasticCache, l::Int64)

Computes the coefficient of the leading `1 / W` angular correction for Legendre order `l`.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `l::Int64` : Legendre order.

# Output Argument(s)
- `coef::Float64` : coefficient such that the leading correction is `coef / W`.
"""
function heavy_leading_1overW_coeff(cache::HeavyInelasticCache, l::Int64)
    if l == 0
        return 0.0
    end
    # P_l'(1) = l*(l+1)/2
    Plp1 = l * (l + 1) / 2.0
    # coef for 1/W term: -A * P_l'(1) * mu_lin
    coef = -cache.A * Plp1 * cache.mu_lin
    return coef
end

"""
    integrate_leading_1overW_per_subshell(cache::HeavyInelasticCache,
    Ef_minus::Float64, Ef_plus::Float64, l::Int64)

Integrates the leading `1 / W` angular correction analytically over an outgoing-energy
group for one subshell and one Legendre order.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef_minus::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef_plus::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].
- `l::Int64` : Legendre order.

# Output Argument(s)
- `σl::Float64` : integral of the leading angular correction over the clipped group.
"""
function integrate_leading_1overW_per_subshell(cache::HeavyInelasticCache, Ef_minus::Float64, Ef_plus::Float64, l::Int64)
    l >= 0 || throw(ArgumentError("Legendre order must be nonnegative."))
    # n=2 permits a zero bound for the identically zero l=0 correction.
    W_lo,W_hi = heavy_feed_interval(cache,Ef_minus,Ef_plus,l == 0 ? 2 : 1)
    if l == 0
        return 0.0
    end
    if W_hi <= W_lo
        return 0.0
    end

    coef = heavy_leading_1overW_coeff(cache, l)
    return coef * heavy_interval_log(W_lo,W_hi)
end

"""
    angular_inelastic_collision_heavy_particle(Ei::Float64, Ef::Float64, L::Int64,
    particle::Particle)

Computes Legendre polynomial values for the inelastic angular deflection of a heavy
particle after an energy loss.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `Ef::Float64` : outgoing particle energy [in m_e*c^2].
- `L::Int64` : maximum Legendre order.
- `particle::Particle` : incoming heavy charged particle.

# Output Argument(s)
- `Pl::Vector{Float64}` : Legendre polynomial values from order 0 to `L`.
"""
function angular_inelastic_collision_heavy_particle(Ei::Float64, Ef::Float64, L::Int64, particle::Particle)
    mₑ = 0.51099895069
    M  = get_mass(particle) / mₑ  # ratio_mass

    W = Ei - Ef
    if W <= 0.0 || Ef <= 0.0 || Ei <= 0.0
        return legendre_polynomials_up_to_L(L, 1.0)  # forward
    end

    pi2 = Ei*(Ei + 2.0*M)
    pf2 = Ef*(Ef + 2.0*M)
    pe2 = W*(W + 2.0)

    denom = 2.0*sqrt(pi2*pf2)
    if denom <= 0.0
        μ = 1.0
    else
        μ = (pi2 + pf2 - pe2)/denom
        μ = clamp(μ, -1.0, 1.0)
    end

    return legendre_polynomials_up_to_L(L, μ)
end

"""
    integrate_inelastic_collision_heavy_particle(Z::Int64, Ei::Float64, n::Int64,
    particle::Particle, W⁻min::Float64=0.0, is_subshells::Bool=true)

Computes the subshell-summed analytical energy-loss moment of the heavy-particle
inelastic collision cross-section for one element.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `n::Int64` : raw energy-loss moment order, 0, 1 or 2.
- `particle::Particle` : incoming heavy charged particle.
- `W⁻min::Float64=0.0` : minimum energy loss W of the incoming particle [in m_e*c^2].
- `is_subshells::Bool=true` : boolean indicating if bounded subshells are resolved.

# Output Argument(s)
- `σn::Float64` : analytical moment summed over all subshells.
"""
function integrate_inelastic_collision_heavy_particle(Z::Int64,Ei::Float64,n::Int64,particle::Particle,W⁻min::Float64=0.0,is_subshells::Bool=true)
    n in (0,1,2) || throw(ArgumentError("Heavy collision raw moments support only n=0,1,2."))
    isfinite(W⁻min) && W⁻min >= 0.0 || throw(ArgumentError("Minimum transfer must be finite and nonnegative."))
    HeavyInelasticCache(1.0,Ei,particle) # validate even if no subshell contributes
    Nshells,Zi,Ui,_,_,_ = electron_subshells(Z,~is_subshells)
    σn = 0
    for δi in range(1,Nshells)
        σn += Zi[δi] * integrate_inelastic_collision_heavy_particle_per_subshell(Ei,n,Ui[δi],particle,W⁻min)
    end
    return σn
end

"""
    integrate_inelastic_collision_heavy_particle_per_subshell(Ei::Float64, n::Int64,
    Ui::Float64, particle::Particle, W⁻min::Float64=0.0)

Computes the analytical energy-loss moment of the heavy-particle inelastic collision
cross-section for one subshell.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `n::Int64` : raw energy-loss moment order, 0, 1 or 2.
- `Ui::Float64` : subshell binding energy [in m_e*c^2].
- `particle::Particle` : incoming heavy charged particle.
- `W⁻min::Float64=0.0` : minimum energy loss W of the incoming particle [in m_e*c^2].

# Output Argument(s)
- `σni::Float64` : analytical moment for the subshell.
"""
function integrate_inelastic_collision_heavy_particle_per_subshell(Ei::Float64,n::Int64,Ui::Float64,particle::Particle,W⁻min::Float64=0.0)
    n in (0,1,2) || throw(ArgumentError("Heavy collision raw moments support only n=0,1,2."))
    all(isfinite, (Ui,W⁻min)) && Ui >= 0.0 && W⁻min >= 0.0 ||
        throw(ArgumentError("Binding energy and minimum transfer must be finite and nonnegative."))
    cache = HeavyInelasticCache(1.0,Ei,particle)
    lower = max(W⁻min,Ui)
    # A threshold above the maximum has empty support; preserve API semantics.
    upper = max(lower,min(Ei,cache.W_ridge))
    return heavy_raw_moment(cache,n,lower,upper)
end

"""
    feed_analytical_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64,
    Ef⁺::Float64)

Computes the scalar analytical feed integral for heavy-particle inelastic collisions over
an outgoing-energy group.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef⁻::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef⁺::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].

# Output Argument(s)
- `σs::Float64` : scalar feed integral over the group.

"""
function feed_analytical_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64, Ef⁺::Float64)
    lo,hi = heavy_feed_interval(cache,Ef⁻,Ef⁺,0)
    hi <= lo && return 0.0
    return heavy_raw_moment(cache,0,lo,hi)
end

"""
    feed_first_moment_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64,
    Ef⁺::Float64)

Computes the first energy-loss moment of the heavy-particle inelastic feed integral over
an outgoing-energy group.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef⁻::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef⁺::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].

# Output Argument(s)
- `M₁::Float64` : first energy-loss moment over the group.

"""
function feed_first_moment_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64, Ef⁺::Float64)
    lo,hi = heavy_feed_interval(cache,Ef⁻,Ef⁺,1)
    hi <= lo && return 0.0
    return heavy_raw_moment(cache,1,lo,hi)
end
