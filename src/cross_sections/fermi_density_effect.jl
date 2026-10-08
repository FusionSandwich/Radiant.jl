# Adapted selectively from zmeskma/Radiant.jl PR #127, commit
# 7054f9ffc5ce107d29fde6a16249aebad1f5566f (CBienvenue/Radiant.jl).
# Input guards and compatibility wrappers were added locally. Empirical accuracy
# claims in the author's source are unverified here; this is not ASTAR qualification.

"""
    fermi_density_effect_from_electrons(Z::Vector{Int64},atz::Vector{Float64},𝒩ₑ_eff::Float64,
    Ei::Float64,state_of_matter::String,type::String,ratio_mass::Float64,I::Float64)

Compute the Fermi density effect correction.

# Input Argument(s)
- `Z::Vector{Int64}`: atomic number of the element(s) composing the material.
- `atz::Vector{Float64}`: atomic fraction of the elements composing the material.
- `𝒩ₑ_eff::Float64`: effective electron density of the material [in cm⁻³], already built
  by the caller. Passing it in rather than rebuilding it from the composition keeps one
  definition of the electron density: the Bethe bracket and the plasma energy cannot
  drift apart.
- `Ei::Float64`: energy of the incoming particle [in mₑc²].
- `type::String`: type of calculation.
- `ratio_mass::Float64`: mass of the incoming particle in unit of the electron mass.
- `I::Float64`: mean excitation energy of the material [in mₑc²], already resolved by the
  caller. It is passed in rather than looked up here so that every term of the Bethe
  formula uses one and the same value.

# Output Argument(s)
- `δF::Float64`: Fermi density effect correction.

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.
- Lorence et al. (1989) : Physics Guide to CEPXS: A Multigroup Coupled Electron-Photon
  Cross-Section Generating Code.

"""
function fermi_density_effect_from_electrons(Z::Vector{Int64},atz::Vector{Float64},𝒩ₑ_eff::Float64,Ei::Float64,state_of_matter::String,type::String,ratio_mass::Float64,I::Float64)

(isfinite(Ei) && Ei >= 0 && isfinite(ratio_mass) && ratio_mass > 0) ||
    throw(ArgumentError("Energy must be finite and nonnegative and mass ratio finite and positive."))
(isfinite(𝒩ₑ_eff) && 𝒩ₑ_eff > 0 && isfinite(I) && I > 0) ||
    throw(ArgumentError("Electron density and excitation energy must be finite and positive."))
validate_atom_fractions(Z,atz)
if type == "fano"

    # Initialization
    γ = (Ei + ratio_mass)/ratio_mass
    β² = (γ^2-1)/γ^2
    Ωp = plasma_energy(𝒩ₑ_eff)

    # Mean electron number per atom.
    Zbar = sum(atz .* Z)

    # Extract the subshell of all atoms of the compound
    Zi = zeros(0); Ui = zeros(0); Wi = zeros(0); Nz = length(Z)
    for i in range(1,Nz)
        _,Ze,Ue,_,_,_ = electron_subshells(Z[i])
        We = resonance_energy(Z[i],Ωp,I)
        append!(Zi,(atz[i]/Zbar).*Ze)
        append!(Ui,Ue)
        append!(Wi,We)
    end

    if (1-β²) < Ωp^2 * sum(Zi./Wi.^2)

        # Find the root L
        f(L) = Ωp^2 * sum(Zi./(Wi.^2 .+ L^2)) - (1-β²)
        dfdL(L) = -2*L*Ωp^2 * sum(Zi./((Wi.^2 .+ L.^2).^2))
        L = newton_bisection(f,dfdL,10*Ωp/sqrt(1-β²),0,1e-5,500)

        # Compute the Fermi density effect
        δF = sum(Zi.* log.(1 .+ L^2 ./Wi.^2)) - L^2/Ωp^2 * (1-β²)

    else
        # No root for L
        δF = 0
    end

elseif type == "sternheimer"

    mₑc² = 0.510999
    Ωp = plasma_energy(𝒩ₑ_eff)

    C = -2*log(I/Ωp) - 1
    # X = log10(βγ) = log10(p/Mc), momentum in units of the incident particle's own
    # mass. sqrt(2*Ei*ratio_mass + Ei^2) is p/(mₑc), so it must be divided by
    # ratio_mass = M/mₑ. Do not "simplify" this away: ratio_mass == 1 for electrons
    # and positrons, which is why the two forms look identical there.
    X = log10(sqrt(2*Ei*ratio_mass+Ei^2)/ratio_mass)

    if state_of_matter ∈ ["solid","liquid"]
        if I ≥ 0.0001/mₑc² && -C ≥ 5.215
            X0 = -0.326*C-1.5
        elseif I < 0.0001/mₑc² && -C ≥ 3.681
            X0 = -0.326*C-1
        else
            X0 = 0.2
        end
        if I ≥ 0.0001/mₑc²
            X1 = 3
        else
            X1 = 2
        end
    elseif state_of_matter ∈ ["gaz"]
        if 13.804 ≤ -C
            X0 = -0.326*C-2.5
        elseif 11.5 ≤ -C < 13.804
            X0 = 2
        elseif 11 ≤ -C < 11.5
            X0 = 1.9
        elseif 10.5 ≤ -C < 11
            X0 = 1.8
        elseif 10 ≤ -C < 10.5
            X0 = 1.7
        else
            X0 = 1.6
        end
        if 12.25 ≤ -C
            X1 = 5
        else
            X1 = 4
        end
    else
        error("The state of matter is either solid, liquid or gaz.")
    end
    B = (-C - (2*log(10.0))*X0)/(X1-X0)^3
    if X ≤ X0
        D = 0.0
    elseif X0 < X ≤ X1 
        D = (2*log(10.0))*X + C + B*(X1-X)^3
    else
        D = (2*log(10.0))*X + C
    end
    δF = D

elseif type == "none"
    δF = 0
else    
    error("Unknown `$type` of Fermi density effect calculations")
end

return δF
end

"""
    fermi_density_effect(Z,ωz,ρ,Ei,state,type,ratio_mass=1.0,I_eff=NaN;atomic_weights=nothing)

Compatibility entry point with mass fractions and mass density. The canonical
electron-density entry point is `fermi_density_effect_from_electrons`.
"""
function fermi_density_effect(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,
    state_of_matter::String,type::String,ratio_mass::Float64=1.0,I_eff::Float64=NaN;
    atomic_weights::Union{Nothing,Vector{Float64}}=nothing)
    atz,_,Ne = stopping_composition(Z,ωz,ρ; atomic_weights=atomic_weights)
    I = stopping_excitation_energy(Z,ωz,atz,I_eff)
    return fermi_density_effect_from_electrons(Z,atz,Ne,Ei,state_of_matter,type,ratio_mass,I)
end
