"""
    plasma_energy(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64;
    atomic_weights::Union{Nothing,Vector{Float64}}=nothing)

Calculate the plasma energy of a free-electron gas with the electron density of the medium.

# Input Argument(s)
- `Z::Vector{Int64}`: atomic number of the element(s) composing the material.
- `ωz::Vector{Float64}`: weight fraction of the element(s) composing the material.
- `ρ::Float64`: density of the material [in g/cm³].
- `atomic_weights::Union{Nothing,Vector{Float64}}`: element atomic weights [u].

# Output Argument(s)
- `Ωp::Float64`: plasma energy [in mₑc²].

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function plasma_energy(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64; atomic_weights::Union{Nothing,Vector{Float64}}=nothing)
    element_atomic_weights = isnothing(atomic_weights) ? standard_atomic_weight.(Z) : atomic_weights
    𝒩ₑ_eff = sum(ωz .* Z .* nuclei_density.(element_atomic_weights, ρ)) # (cm⁻³)
    return plasma_energy(𝒩ₑ_eff)
end

"""
    plasma_energy(𝒩ₑ_eff::Float64)

Calculate the plasma energy of a free-electron gas with the electron density of the medium.

# Input Argument(s)
- `𝒩ₑ_eff::Float64`: effective electron density in the medium.

# Output Argument(s)
- `Ωp::Float64`: plasma energy [in mₑc²].

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function plasma_energy(𝒩ₑ_eff::Float64)
    ħc = 3.86189603E-11                          # (mₑc²⋅cm)
    rₑ = 2.81794092E-13                          # (cm)
    Ωp = sqrt(4*π*𝒩ₑ_eff*ħc^2*rₑ)
    return Ωp
end

# Composition foundation shared by stopping and plasma-energy callers.
function validate_atom_fractions(Z::Vector{Int64},atz::Vector{Float64})
    (!isempty(Z) && length(Z) == length(atz) && all(z -> 1 <= z <= 100,Z)) ||
        throw(ArgumentError("Atomic numbers and fractions must have matching nonempty lengths, with Z in 1:100."))
    (all(x -> isfinite(x) && x >= 0,atz) && isapprox(sum(atz),1.0; atol=1e-8,rtol=1e-8)) ||
        throw(ArgumentError("Atomic fractions must be finite, nonnegative and sum to one."))
    return nothing
end

function stopping_composition(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64;
    atomic_weights::Union{Nothing,Vector{Float64}}=nothing)
    (!isempty(Z) && length(Z) == length(ωz) && all(z -> 1 <= z <= 100,Z)) ||
        throw(ArgumentError("Composition must have matching nonempty lengths and Z in 1:100."))
    (all(w -> isfinite(w) && w >= 0,ωz) && sum(ωz) > 0) ||
        throw(ArgumentError("Mass fractions must be finite, nonnegative and have positive sum."))
    (isfinite(ρ) && ρ > 0) || throw(ArgumentError("Mass density must be finite and positive."))
    A = isnothing(atomic_weights) ? standard_atomic_weight.(Z) : atomic_weights
    (length(A) == length(Z) && all(a -> isfinite(a) && a > 0,A)) ||
        throw(ArgumentError("Atomic weights must match the composition and be finite and positive."))
    atoms = ωz ./ A
    atz = atoms ./ sum(atoms)
    N_density = ρ * Nₐ * sum(atoms)
    Ne = N_density * sum(atz .* Z)
    return atz,N_density,Ne
end

function validate_stopping_inputs(Z,ωz,atz,ρ,N_density,Ei,I_eff; atomic_weights=nothing)
    a,N,_ = stopping_composition(Z,ωz,ρ; atomic_weights=atomic_weights)
    validate_atom_fractions(Z,atz)
    (isfinite(N_density) && N_density > 0) ||
        throw(ArgumentError("Nuclei density must be finite and positive."))
    (all(isapprox.(atz,a; atol=1e-8,rtol=1e-8)) && isapprox(N_density,N;rtol=1e-8)) ||
        throw(ArgumentError("Atomic fractions and nuclei density disagree with mass fractions, density or isotope masses."))
    (isfinite(Ei) && Ei >= 0) || throw(ArgumentError("Energy must be finite and nonnegative."))
    (isnan(I_eff) || (isfinite(I_eff) && I_eff > 0)) ||
        throw(ArgumentError("Excitation energy must be finite and positive, or NaN for the default."))
    return nothing
end

# Keep canonical compound I values, then use electron fractions for Bragg's log rule.
function stopping_excitation_energy(Z,ωz,atz,I_override::Float64=NaN)
    if !isnan(I_override)
        (isfinite(I_override) && I_override > 0) ||
            throw(ArgumentError("Excitation energy must be finite and positive."))
        return I_override
    end
    Itab = tabulated_mean_excitation_energy(Z,ωz)
    !isnan(Itab) && return Itab
    g = atz .* Z
    return exp(sum(g .* log.(mean_excitation_energy.(Z))) / sum(g))
end
