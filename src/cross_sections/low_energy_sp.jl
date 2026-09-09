"""
    low_energy_stopping(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},
    ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,
    density_correction::String,state_of_matter::String,I_eff::Float64,Ecut::Float64,
    Emax::Float64)

Gives the collisional stopping power below `Ecut`, where the corrected Bethe formula is no
longer valid.

Everything Radiant computes below `Ecut` lives in this file. `bethe`'s low-energy branch is
a single call into it and holds no low-energy physics of its own.

The analytical extension solves both coefficients of

    S = exp(A - B*(log t)^1.5)        t = Ei/Emax ≥ 1
    S = 1.5*exp(A)*sqrt(t)*(1 - t/3)  t < 1

from the value and the slope of the Bethe branch at `Ecut`, so the join is smooth. Both
branches have zero slope at `t = 1`, which puts the stopping-power maximum at `Emax`.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `N_density::Float64` : atom density of the material [in cm⁻³].
- `Ei::Float64` : incoming particle energy [in mₑc²].
- `particle::Particle` : incoming particle.
- `density_correction::String` : density-effect correction type.
- `state_of_matter::String` : material state.
- `I_eff::Float64` : user-defined mean excitation energy override [in mₑc²]; NaN ⟹ table.
- `Ecut::Float64` : energy below which this branch is used [in mₑc²].
- `Emax::Float64` : energy of the stopping-power maximum of the extension [in mₑc²].
- `atomic_weights::Union{Nothing,Vector{Float64}}` : element atomic weights [u].
- `is_bethe_correction::Bool` : boolean to enable or not the Bethe correction.
- `is_density_correction::Bool` : boolean to enable or not the density correction.
- `is_shell_correction::Bool` : boolean to enable or not the shell correction.
- `is_lindhard_sorensen::Bool` : boolean to enable or not the Lindhard-Sorensen correction.
- `is_barkas::Bool` : boolean to enable or not the Barkas correction.
- `ecut_MeV::Float64` : the cutoff `bethe` was given [in MeV], forwarded unchanged to the
  anchor calls so they resolve the same `Ecut` and take the corrected-Bethe branch.

# Output Argument(s)
- `S::Float64` : stopping power.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function low_energy_stopping(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,density_correction::String,state_of_matter::String,I_eff::Float64,Ecut::Float64,Emax::Float64; atomic_weights::Union{Nothing,Vector{Float64}}=nothing,is_bethe_correction::Bool=true,is_density_correction::Bool=true,is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true,ecut_MeV::Float64=NaN)

    #----
    # Estimate parameters A and B
    #----
    # The anchor energies are both ≥ Ecut, so these two calls take the corrected-Bethe
    # branch and never come back here: the recursion is one level deep by construction.
    # The correction switches are forwarded unchanged, or the extension would be anchored
    # on a different formula from the one it is joined to.
    h = 0.00005
    Sc = bethe(Z,ωz,atz,ρ,N_density,Ecut,particle,density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_bethe_correction=is_bethe_correction,
        is_density_correction=is_density_correction,
        is_shell_correction=is_shell_correction,
        is_lindhard_sorensen=is_lindhard_sorensen,
        is_barkas=is_barkas,
        ecut_MeV=ecut_MeV)
    Sc_plus_h = bethe(Z,ωz,atz,ρ,N_density,Ecut+h,particle,density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_bethe_correction=is_bethe_correction,
        is_density_correction=is_density_correction,
        is_shell_correction=is_shell_correction,
        is_lindhard_sorensen=is_lindhard_sorensen,
        is_barkas=is_barkas,
        ecut_MeV=ecut_MeV)
    dSdE = (Sc_plus_h - Sc)/h
    A = (3*log(Sc)*Sc - 2*dSdE*Ecut*log(Ecut/Emax))/(3*Sc)
    B = -2*Ecut*dSdE/(3*Sc*sqrt(log(Ecut/Emax)))

    #----
    # Compute total cross-sections using analytical formula
    #----
    t = Ei/Emax
    if t ≥ 1
        S = exp(A-B*(log(t))^1.5)
    else
        S = 1.5*exp(A)*sqrt(t)*(1-t/3)
    end
    return S
end
