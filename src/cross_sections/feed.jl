"""
    feed(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,
    is_elastic::Bool,is_subshells::Bool)

Calculate the feed function 𝓕 (normalized probability of scattering from Ei into each
group gf) for each Legendre moment up to order L. Also calculate the energy weighted
feed function 𝓕ₑ for energy-deposition cross section.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `atz::Vector{Float64}` : atomic percent of the element(s) composing the material.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction informations.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` :  number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of the particles imply in the interaction.
- `type::String` : type of interaction (scattering or production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic::Bool` : boolean indicating if the scattering is elastic.
- `is_subshells::Bool` : boolean indicating if the cross-sections are subshells dependant.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function (per unit nuclei density).
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function (per unit nuclei density).

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool)

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)
ΔQ = get_mass_energy_variation(interaction,type,true)

# Outgoing particle energy spectrum
is_dirac, Np, q_type = out_distribution_dispatch(interaction,type)
if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

#----
# Feed function over all groups and under the cutoff energy
#----

# Heavy inelastic S (same scattered particle)
is_heavy_inelastic_S = (is_proton(incoming_particle) || is_alpha(incoming_particle)) && (incoming_particle == scattered_particle)

# Loop over the compound elements
Nz = length(Z)
for i in range(1,Nz)

    # Loop over subshells and outgoing groups
    Nshells,Zi,Ui,Ti,ri,_ = electron_subshells(Z[i],~is_subshells)
    for gf in range(1,Ng), δi in range(1,Nshells)

        # Final energy group
        Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
        Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui[δi],Ec,incoming_particle)
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺

        # Integration over the energy group
        𝓕i = zeros(L+1)
        𝓕iₑ = 0

        # For heavy particles compute cache once and reuse for all analytic integrals
        analytic_A = 0.0
        M₁ = 0.0
        cache = nothing
        if is_heavy_inelastic_S
            cache = HeavyInelasticCache(Zi[δi], Ei, incoming_particle)
            analytic_A = integrate_A_over_W2_per_subshell(cache, Ef⁻, Ef⁺) * atz[i]
            M₁ = feed_first_moment_heavy_particle(cache, Ef⁻, Ef⁺) * atz[i]
        end

        # Use quadrature integration for all particles
        for n in range(1,Np)
            # Outgoing particle energy group
            if (is_elastic) Ef = Ei else Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2 end

            # Compute Legendre angular flux moments
            Σsᵢ = ΔEf .* w[n]/2 .* dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui[δi],Zi[δi],Ti[δi],ri[δi],Ec,incoming_particle) * atz[i]
            if is_dirac Σsᵢ /= ΔEf end
            𝓕i .+= Σsᵢ
            𝓕iₑ += Σsᵢ[1] * (Ef + ΔQ)
        end

        # Add analytic singular contribution
        if is_heavy_inelastic_S
            𝓕i .+= analytic_A .* ones(L+1)
            for l in range(0,L)
                lead_log = integrate_leading_1overW_per_subshell(cache, Ef⁻, Ef⁺, l) * atz[i]
                𝓕i[l+1] += lead_log
            end
            σ_analytic = feed_analytical_heavy_particle(cache, Ef⁻, Ef⁺) * atz[i]
            𝓕iₑ = Ei * σ_analytic - M₁

            if ~isapprox(𝓕i[1], σ_analytic; rtol=1e-3, atol=1e-12)
                rel = abs(𝓕i[1] - σ_analytic) / max(abs(σ_analytic), 1e-300)
                print("FEED_WARN: Z=$(Z[i]), δi=$(δi), gf=$(gf), σ_analytic=$(σ_analytic), numeric_l0=$(𝓕i[1]), rel_diff=$(rel)\n")
                𝓕i[1] = σ_analytic
            end
        end
        𝓕[gf,:] .+= 𝓕i
        𝓕ₑ[gf] += 𝓕iₑ
    end
end
return 𝓕, 𝓕ₑ
end

"""
    feed_elastic_scattering(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,
    is_subshells::Bool,A::Vector{Vector{Int64}},
    atpercentA::Vector{Vector{Float64}})

Calculate the elastic-scattering feed function 𝓕 from incident energy `Ei` into each
outgoing energy group and Legendre moment up to order `L`. Also calculate the
energy-weighted feed function 𝓕ₑ for energy-deposition cross sections, including isotope
fractions when isotope-resolved data are provided.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `atz::Vector{Float64}` : atomic percent of the element(s) composing the material.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction information.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` : number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of particles involved in the interaction.
- `type::String` : type of interaction (`"S"` for scattering or `"P"` for production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic::Bool` : boolean indicating if the outgoing particle energy is equal to `Ei`.
- `is_subshells::Bool` : boolean indicating if subshell-dependent cross sections are used.
- `A::Vector{Vector{Int64}}` : isotope mass numbers per element.
- `atpercentA::Vector{Vector{Float64}}` : isotope atomic fractions per element.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy-weighted feed function.

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed_elastic_scattering(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool,A::Vector{Vector{Int64}},atpercentA::Vector{Vector{Float64}})

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)

#----
# Feed function over all groups and under the cutoff energy
#----

# Loop over the compound elements
Nz = length(Z)
for i in range(1,Nz)
    # Loop over isotopes
    for (Ai, atai) in zip(A[i], atpercentA[i])
        if type == "P" && !(Z[i] == 1 && Ai == 1)
            continue
        end
        δi = 0
        Ui = 0.0
        Zi = Z[i]
        Ti = 0.0
        ri = 0.0
        for gf in range(1,Ng)

            # Final energy group
            Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
            M_target = get_mass(Z[i], Ai)
            Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui,Ec,incoming_particle,M_target)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺

            # Integration over the energy group
            𝓕i = zeros(L+1)
            𝓕iₑ = 0
            Ef = is_elastic ? Ei : (Ef⁻ + Ef⁺) / 2
            Σsᵢ = dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui,Zi,Ti,ri,Ec,incoming_particle,Ai) * atz[i] * atai
            𝓕i .+= Σsᵢ
            𝓕iₑ += Σsᵢ[1] * Ef
            𝓕[gf,:] .+= 𝓕i
            𝓕ₑ[gf] += 𝓕iₑ
        end
    end
end
return 𝓕, 𝓕ₑ
end

"""
    feed_nuclear_reaction(Z::Vector{Int64}, atz::Vector{Float64},
    ωz::Vector{Float64}, ρ::Float64, N_density::Float64, state_of_matter::String,
    L::Int64, Ei::Float64, Eout::Vector{Float64}, Ng::Int64,
    interaction::Nuclear_Reaction, gi::Int64, Ngi::Int64,
    particles::Vector{Particle}, type::String,
    incoming_particle::Particle, Ein::Vector{Float64},
    I_eff::Float64,
    A::Vector{Vector{Int64}},
    atpercentA::Vector{Vector{Float64}})

Calculate the nuclear-reaction feed function 𝓕 for type `"P"` (equivalent-proton
production following Salvat & Quesada 2020 §4.2). The charged particles a reaction emits
are replaced by the protons of equal CSDA range, weighted by w = E_b/E_eq so that the
energy they carry is preserved.

For each isotope, `equivalent_proton_spectra` interpolates the outgoing spectrum of every
production channel at the incident energy and carries its energy points onto the
equivalent-proton scale; each outgoing group is then integrated over those spectra by
`integrate_equivalent_proton_group`. A group narrower than the tabulation receives its
share of the interval covering it and one wider than it receives the whole, so both the
number of equivalent protons and the energy they carry are conserved whatever the group
structure.

The angular distribution is isotropic, so only the Legendre moment L=0 is filled.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of each element.
- `ωz::Vector{Float64}` : weight fraction of each element.
- `ρ::Float64` : material density [g/cm³].
- `state_of_matter::String` : material state (`"solid"`, `"liquid"`, or `"gas"`).
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy [mₑc²].
- `Eout::Vector{Float64}` : outgoing energy group boundaries [mₑc²], decreasing.
- `Ng::Int64` : number of outgoing energy groups.
- `interaction::Nuclear_Reaction` : nuclear-reaction interaction structure.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` : number of incoming particle groups.
- `particles::Vector{Particle}` : particles in the simulation.
- `type::String` : interaction type string (`"P"`).
- `incoming_particle::Particle` : incoming particle.
- `Ein::Vector{Float64}` : incoming energy group boundaries [mₑc²].
- `I_eff::Float64` : effective mean excitation energy override [mₑc²]; NaN ⟹ tables.
- `A::Vector{Vector{Int64}}` : isotope mass numbers per element.
- `atpercentA::Vector{Vector{Float64}}` : isotope atomic fractions per element.

# Output Argument(s)
- `𝓕::Array{Float64,2}` : feed function [cm²], size (Ng+1, L+1).
- `𝓕ₑ::Vector{Float64}` : energy-weighted feed function [cm² × mₑc²], size (Ng+1).

# Reference(s)
- Salvat & Quesada (2020), NIMB 475, 49–62, §4.2 equivalent-proton method.
- MacFarlane et al. (2021), The NJOY Nuclear Data Processing System, Version 2012.
"""
function feed_nuclear_reaction(Z::Vector{Int64}, atz::Vector{Float64},
        ωz::Vector{Float64}, ρ::Float64, N_density::Float64, state_of_matter::String,
        L::Int64, Ei::Float64, Eout::Vector{Float64}, Ng::Int64,
        interaction::Nuclear_Reaction, _gi::Int64, _Ngi::Int64,
        _particles::Vector{Particle}, _type::String,
        incoming_particle::Particle, _Ein::Vector{Float64},
        I_eff::Float64,
        A::Vector{Vector{Int64}},
        atpercentA::Vector{Vector{Float64}})

#----
# Initialization
#----
𝓕  = zeros(Ng+1, L+1)
𝓕ₑ = zeros(Ng+1)

ptype = get_type(incoming_particle)
haskey(interaction.production_db, ptype) || error("Requested nuclear production feed database is missing for $ptype.")
E_in = Ei * M_E_C2_MEV * 1e6   # mₑc² -> eV, the unit of the ENDF distributions

#----
# Feed function over all groups
#----

# Loop over the compound elements
Nz = length(Z)
for i in range(1,Nz)

    # Loop over isotopes
    for (Ai, atai) in zip(A[i], atpercentA[i])
        haskey(interaction.production_db[ptype], (Z[i], Ai)) || error("Requested nuclear production feed isotope ($(Z[i]),$Ai) is missing.")
        channels = interaction.production_db[ptype][(Z[i], Ai)]
        isempty(channels) && error("Requested nuclear production feed has no evaluated channels for ($(Z[i]),$Ai).")

        # Energy leaving with the neutrons and the photons, which nothing transports. It
        # is added to the energy-weighted feed without any matching entry in 𝓕, so that
        # Σe = Σtₑ - ΣΣsₑ loses it while no particle is put in a group for it. Σsₑ enters
        # the balance only through its sum over the outgoing groups, so the group it is
        # filed under is immaterial; it goes in the first one.
        if interaction.get_is_neutral_escape()
            𝓕ₑ[1] += neutral_energy_release(channels, E_in) * atz[i] * atai
        end

        # Outgoing spectra of this isotope, carried onto the equivalent-proton scale
        spectra = equivalent_proton_spectra(interaction, channels, E_in, Z[i], Ai,
                                            Z, ωz, atz, ρ, N_density, state_of_matter,
                                            I_eff)
        isempty(spectra) && continue

        # Loop over outgoing groups
        for gf in range(1,Ng)

            # Final energy group
            Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]

            # Integration over the energy group
            𝓕i, 𝓕iₑ = integrate_equivalent_proton_group(spectra, Ef⁺, Ef⁻, L)
            𝓕i[1] == 0.0 && continue

            𝓕[gf,:] .+= 𝓕i  .* (atz[i] * atai)
            𝓕ₑ[gf]  += 𝓕iₑ  * (atz[i] * atai)
        end
    end
end
return 𝓕, 𝓕ₑ
end
