# Optional stopping controls adapted from zmeskma/Radiant.jl PR #127,
# commit 7054f9ffc5ce107d29fde6a16249aebad1f5566f. Guards added locally.
"""
Inelastic_Collision

Structure used to define parameters for production of multigroup inelastic collisional cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

"""
mutable struct Inelastic_Collision <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_subshells_dependant::Bool
    is_shell_correction::Bool
    is_focusing_møller::Bool
    is_hydrogenic_distribution_term::Bool
    is_distant_collision::Bool
    density_correction::String
    scattering_model::String
    user_ecut::Dict{Type,Float64}
    user_emax::Dict{Type,Float64}
    user_emax_fit::Set{Type}
    low_energy_model::String
    low_energy_source::String
    is_bragg_additivity::Bool

    # Constructor(s)
    function Inelastic_Collision()
        this = new()
        this.name = "Inelastic_Collision"
        this.interaction_types = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"],(Proton,Proton) => ["S"],(Proton,Electron) => ["P"],(Alpha,Alpha) => ["S"],(Alpha,Electron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.set_density_correction("fano")
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_focusing_møller = false
        this.is_hydrogenic_distribution_term = false
        this.set_is_subshells_dependant(true)
        this.set_is_shell_correction(true)
        this.set_scattering_model("BFP")
        this.user_ecut = Dict{Type,Float64}()
        this.user_emax = Dict{Type,Float64}()
        this.user_emax_fit = Set{Type}()
        this.set_low_energy_model("analytic")
        this.set_low_energy_source("srim")
        this.set_is_bragg_additivity(false)

        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Inelastic_Collision,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for inelastic collision processes.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # No knock-on electrons.
```
"""
function set_interaction_types(this::Inelastic_Collision,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_density_correction(this::Inelastic_Collision,density_correction::String)

Set the Fermi density correction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `density_correction::String` : type of density effect:
    - `fano` : Fano density effect.
    - `sternheimer` : Sternheimer semi-empirical density effect.
    - `none` : no density effect.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_density_correction("sternheimer")
```
"""
function set_density_correction(this::Inelastic_Collision,density_correction::String)
    if lowercase(density_correction) ∉ ["fano","sternheimer","none"] error("Unknown '$density_correction' density correction.") end
    this.density_correction = lowercase(density_correction)
end

"""
    set_is_shell_correction(this::Inelastic_Collision,is_shell_correction::Bool)

Activate or desactivate shell correction for stopping powers.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_shell_correction::Bool` : activate (true) or desactivate (false) the shell correction.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_is_shell_correction(false)
```
"""
function set_is_shell_correction(this::Inelastic_Collision,is_shell_correction::Bool)
    this.is_shell_correction = is_shell_correction
end

"""
    set_is_subshells_dependant(this::Inelastic_Collision,is_subshells_dependant::Bool)

Compute the inelastic cross-sections assuming bounded or unbounded electrons.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_subshells_dependant::Bool` : bounded (true) or unbounded (false) electrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_is_subshells_dependant(false)
```
"""
function set_is_subshells_dependant(this::Inelastic_Collision,is_subshells_dependant::Bool)
    this.is_subshells_dependant = is_subshells_dependant
end

"""
    scattering_model(this::Inelastic_Collision,scattering_model::String)

To define the solver for inelastic scattering.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic structure.
- `solver::String` : solver for inelastic scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Inelastic_Collision,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    set_ecut(this::Inelastic_Collision,particle::Type,energy::Real)

To set the energy above which the corrected Bethe formula is used for a given particle,
and below which the stopping power is taken from the low-energy branch.

Lowering it extends the Bethe formula into a region where its corrections are less well
determined; the shell-correction tables stop at a tabulated γ and raise an error naming
this setter when they run out.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type, one of `Electron`, `Positron`, `Muon`, `Antimuon`,
  `Proton`, `Antiproton` or `Alpha`.
- `energy::Real` : cutoff energy [in MeV].

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.set_ecut(Proton,1.0)
```
"""
function set_ecut(this::Inelastic_Collision,particle::Type,energy::Real)
    if !haskey(ECUT_DEFAULTS,particle) error("The cutoff energy is undefined for the '$particle' particle. It is defined for $(collect(keys(ECUT_DEFAULTS))).") end
    (isfinite(energy) && energy > 0) ||
        throw(ArgumentError("The cutoff energy must be finite and positive."))
    emax = get_emax(this,particle)
    (emax == "f" || energy > emax) ||
        throw(ArgumentError("The cutoff energy must exceed the analytic stopping maximum."))
    this.user_ecut[particle] = float(energy)
end

"""
    get_ecut(this::Inelastic_Collision,particle::Type)

To get the energy above which the corrected Bethe formula is used for a given particle:
the value `set_ecut` stored, or `ECUT_DEFAULTS[particle]` if none was set.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type.

# Output Argument(s)
- `energy::Float64` : cutoff energy [in MeV].

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_ecut(Proton)
0.75
```
"""
function get_ecut(this::Inelastic_Collision,particle::Type)
    if !haskey(ECUT_DEFAULTS,particle) error("The cutoff energy is undefined for the '$particle' particle. It is defined for $(collect(keys(ECUT_DEFAULTS))).") end
    return get(this.user_ecut,particle,ECUT_DEFAULTS[particle])
end

"""
    set_emax(this::Inelastic_Collision,particle::Type,energy::Real)
    set_emax(this::Inelastic_Collision,particle::Type,mode::String)

Set a finite positive analytic maximum energy [MeV] below Ecut, or "f" for the author-supplied empirical maximum fit for Proton and Alpha. The fit is unqualified here.
"""
function set_emax(this::Inelastic_Collision,particle::Type,energy::Real)
    if !haskey(EMAX_DEFAULTS,particle) error("The stopping-power maximum is undefined for the '$particle' particle. It is defined for $(collect(keys(EMAX_DEFAULTS))).") end
    (isfinite(energy) && 0 < energy < get_ecut(this,particle)) ||
        throw(ArgumentError("The analytic stopping maximum must be finite and satisfy 0 < Emax < Ecut."))
    delete!(this.user_emax_fit,particle)
    this.user_emax[particle] = float(energy)
end

function set_emax(this::Inelastic_Collision,particle::Type,mode::String)
    if lowercase(mode) != "f" error("Unknown '$mode' mode for the stopping-power maximum. Give an energy [in MeV], or \"f\" for the Z-dependent fit.") end
    if particle ∉ [Proton,Alpha] error("The \"f\" mode is only supported for Proton and Alpha; there is no Z-dependent fit for the '$particle' particle.") end
    delete!(this.user_emax,particle)
    push!(this.user_emax_fit,particle)
end

"""
    get_emax(this::Inelastic_Collision,particle::Type)

To get the energy of the stopping-power maximum of the analytic low-energy extension.

Gives the value `set_emax` stored, or `EMAX_DEFAULTS[particle]` if none was set. Returns
`"f"` unchanged when the Z-dependent fit is selected: the resolved value depends on the
material and cannot be known here.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type.

# Output Argument(s)
- `energy::Union{Float64,String}` : energy of the stopping-power maximum [in MeV], or `"f"`.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_emax(Proton)
0.05
```
"""
function get_emax(this::Inelastic_Collision,particle::Type)
    if !haskey(EMAX_DEFAULTS,particle) error("The stopping-power maximum is undefined for the '$particle' particle. It is defined for $(collect(keys(EMAX_DEFAULTS))).") end
    if particle in this.user_emax_fit return "f" end
    return get(this.user_emax,particle,EMAX_DEFAULTS[particle])
end

"""
    set_low_energy_model(this::Inelastic_Collision,model::String)

Select "analytic" (default) or optional "icru49" stopping. The empirical branch needs separately qualified coefficient data, uses a smooth join above Ecut, and is unqualified below 10 keV. Its soft/catastrophic split uses a separate analytic kernel.
"""
function set_low_energy_model(this::Inelastic_Collision,model::String)
    if lowercase(model) ∉ ["analytic","icru49"] error("Unknown '$model' low-energy model.") end
    this.low_energy_model = lowercase(model)
end

"""
    get_low_energy_model(this::Inelastic_Collision)

To get the model used for the stopping power below `Ecut`.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `model::String` : low-energy model.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_low_energy_model()
"analytic"
```
"""
function get_low_energy_model(this::Inelastic_Collision)
    return this.low_energy_model
end

"""
    set_low_energy_source(this::Inelastic_Collision,source::String)

Select "srim" or "iaea" coefficients for optional empirical stopping. IAEA entries override SRIM per element and phase; missing entries retain SRIM. Dataset coverage and accuracy are unverified here.
"""
function set_low_energy_source(this::Inelastic_Collision,source::String)
    if lowercase(source) ∉ ["srim","iaea"] error("Unknown '$source' low-energy coefficient source.") end
    this.low_energy_source = lowercase(source)
end

"""
    get_low_energy_source(this::Inelastic_Collision)

To get the coefficient source of the fitted low-energy branch. Inert unless the `icru49`
low-energy model is selected.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `source::String` : coefficient source.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_low_energy_source()
"srim"
```

# Reference(s)
- IAEA, Electronic Stopping Power of Matter for Ions (database),
  https://www-nds.iaea.org/stopping/. Sole source of the "iaea" coefficients.
"""
function get_low_energy_source(this::Inelastic_Collision)
    return this.low_energy_source
end

"""
    set_is_bragg_additivity(this::Inelastic_Collision,is_bragg_additivity::Bool)

Enable optional mass-fraction additivity of isolated elemental stopping models evaluated at the same mass density. This overrides compound I and uses each elemental electron density; it is an alternative approximation.
"""
function set_is_bragg_additivity(this::Inelastic_Collision,is_bragg_additivity::Bool)
    this.is_bragg_additivity = is_bragg_additivity
end

"""
    get_is_bragg_additivity(this::Inelastic_Collision)

To know whether compounds are built from their elements by Bragg additivity.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_bragg_additivity::Bool` : whether Bragg additivity is enabled.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_is_bragg_additivity()
false
```
"""
function get_is_bragg_additivity(this::Inelastic_Collision)
    return this.is_bragg_additivity
end

"""
    in_distribution(this::Inelastic_Collision)

Describe the energy discretization method for the incoming particle in the inelastic collision
interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Inelastic_Collision)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Inelastic_Collision)

Describe the energy discretization method for the outgoing particle in the inelastic collision
interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Inelastic_Collision)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Inelastic_Collision,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64,Ui::Float64,particle::Particle)

Gives the integration energy bounds for the outgoing particle for inelastic collision
interaction. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.
- `Ec::Float64` : energy cutoff between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Inelastic_Collision,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64,Ui::Float64,particle::Particle)
    if is_electron(particle)
        # Scattered electron
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,(Ei-Ui)/2)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,(Ei-Ui)/2)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    elseif is_positron(particle)
        # Scattered positron/proton/alpha
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    elseif is_proton(particle) || is_alpha(particle)
        mₑ = 0.51099895069
        ratio_mass = get_mass(particle)/mₑ
        γ = (Ei + ratio_mass)/ratio_mass
        β² = (γ^2-1)/γ^2
        R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
        W_ridge = 2*β²*γ^2*R
        # Scattered heavy particle
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0,Ei-W_ridge)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui,W_ridge-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    else
        error("Unknown particle")
    end
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Inelastic_Collision,L::Int64,Ei::Float64,Ef::Float64,type::String,
    particle::Particle,Ui::Float64,Zi::Real,Ti::Float64)

Gives the Legendre moments of the scattering cross-sections for inelastic collision
interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `Z::Int64` : atomic number.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `type::String` : type of interaction.
- `particle::Particle` : incoming particle.
- `δi::Int64` : subshell index.
- `Ui::Float64` : binding energy of the subshell.
- `Zi::Int64` : number of electrons in the subshell.
- `Ti::Float64` : mean kinetic energy of electrons in the subshell.
- `ri::Float64` : average radius per subshell.

# Output Argument(s)
- `σl::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Inelastic_Collision,Z::Int64,L::Int64,Ei::Float64,Ef::Float64,type::String,particle::Particle,δi::Int64,Ui::Float64,Zi::Real,Ti::Float64,ri::Float64)

    #----
    # Initialization
    #----
    σs = 0.0
    σl = zeros(L+1)
    if type == "S"
        E⁺ = Ef
        E⁻ = Ei-Ui-E⁺
        W = E⁻+Ui
    elseif type == "P"
        E⁻ = Ef
        E⁺ = Ei-Ui-E⁻
        W = E⁻+Ui
    else
        error("Unknown type")
    end
    if W < 0 return σl end

    #----
    # Close collisions
    #----
    is_heavy = false
    sigma_reg = 0.0
    A_over_W2 = 0.0
    cache = nothing
    if is_electron(particle)
        σs += moller(Zi,Ei,W,Ui,Ti,this.is_focusing_møller,this.is_hydrogenic_distribution_term)
    elseif is_positron(particle)
        σs += bhabha(Zi,Ei,W)
    elseif is_proton(particle) || is_alpha(particle)
        is_heavy = true
        # Create cache once for all Legendre orders
        cache = HeavyInelasticCache(Zi, Ei, particle)
        sigma_reg, A_over_W2 = heavy_sigma_components(cache, W)
    else
        error("Unknown particle")
    end

    #----
    # Compute the angular distribution
    #----
    if is_heavy
        Wl = angular_inelastic_collision_heavy_particle(Ei,Ef,L,particle)
    else
        Wl = angular_moller(Ei,Ef,L)
    end

    #----
    # Compute the Legendre moments of the cross-section
    #----
    if is_heavy
        for l in range(0,L)
            # Subtract analytic leading 1/W piece to regularize the numeric integrand
            if W > 0
                coef = heavy_leading_1overW_coeff(cache, l)
                leading_term = coef / W
            else
                leading_term = 0.0
            end
            σl[l+1] += Wl[l+1] * sigma_reg + (Wl[l+1] - 1.0) * A_over_W2 - leading_term
        end
    else
        for l in range(0,L) σl[l+1] += Wl[l+1] * σs end
    end
    return σl
end

"""
    tcs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the total cross-section for inelastic collision interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

    #----
    # Initialization
    #----
    σt = 0

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σt += integrate_moller(Z,Ei,0,Ei-Ec,this.is_focusing_møller,this.is_hydrogenic_distribution_term,this.is_subshells_dependant)
    elseif is_positron(particle)
        σt += integrate_bhabha(Z,Ei,0,Ei-Ec,this.is_subshells_dependant)
    elseif is_proton(particle) || is_alpha(particle)
        σt += integrate_inelastic_collision_heavy_particle(Z,Ei,0,particle,Ei-Ec,this.is_subshells_dependant)
    else
        error("Unknown particle")
    end

    return σt
end

"""
    acs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the absorption cross-section for inelastic collision interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `Ecutoff::Float64` : cutoff energy.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64,Ecutoff::Float64)

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σa = integrate_moller(Z,Ei,0,Ei-min(Ec,Ecutoff),this.is_focusing_møller,this.is_hydrogenic_distribution_term,this.is_subshells_dependant)
    elseif is_positron(particle)
        σa = integrate_bhabha(Z,Ei,0,Ei-min(Ec,Ecutoff),this.is_subshells_dependant)
    else
        error("Unknown particle")
    end
    return σa
    
end

"""
    sp(this::Inelastic_Collision,Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},
    ρ::Float64,N_density::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,
    particle::Particle,I_eff::Float64,atomic_weights::Union{Nothing,Vector{Float64}})

Gives the stopping power for inelastic collision interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `N_density::Float64` : nuclei density of the material [in cm⁻³].
- `state_of_matter::String` : state of matter.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `I_eff::Float64` : effective mean excitation energy override [in mₑc²]; NaN to use tables.
- `atomic_weights::Union{Nothing,Vector{Float64}}` : element atomic weights [u].

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Inelastic_Collision,Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,N_density::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle,I_eff::Float64=NaN,atomic_weights::Union{Nothing,Vector{Float64}}=nothing)

    # Compute the total stopping power
    Stot = bethe(Z,ωz,atz,ρ,N_density,Ei,particle,this.density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_shell_correction=this.is_shell_correction,
        user_ecut=get(this.user_ecut,get_type(particle),NaN),
        user_emax=get(this.user_emax,get_type(particle),NaN),
        is_emax_fit=(get_type(particle) in this.user_emax_fit),
        low_energy_model=this.low_energy_model,
        low_energy_source=this.low_energy_source,
        is_bragg_additivity=this.is_bragg_additivity)

    # Compute the catastrophic Møller- or Bhabha- derived stopping power
    Sc = 0
    Nz = length(Z)
    for i in range(1,Nz)

        # Close collision
        if is_electron(particle)
            Sc += atz[i] * N_density * integrate_moller(Z[i],Ei,1,Ei-Ec,this.is_focusing_møller,this.is_hydrogenic_distribution_term,this.is_subshells_dependant)
        elseif is_positron(particle)
            Sc += atz[i] * N_density * integrate_bhabha(Z[i],Ei,1,Ei-Ec,this.is_subshells_dependant)
        elseif is_proton(particle) || is_alpha(particle)
            Sc += atz[i] * N_density * integrate_inelastic_collision_heavy_particle(Z[i],Ei,1,particle,Ei-Ec,this.is_subshells_dependant)
        end
    end

    # Compute the soft stopping power
    Sr = checked_soft_stopping(Stot,Sc)

    return Sr
end

"""
    mt(this::Inelastic_Collision)

Gives the momentum transfer for inelastic collision interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Inelastic_Collision)
    T = 0.0
    return T
end



"""
    checked_soft_stopping(total,catastrophic)

Compute the residual electronic stopping. A negative or nonfinite partition raises
`DomainError`, including small negative residuals; no clipping conceals a mismatch.
"""
function checked_soft_stopping(total::Real,catastrophic::Real)
    (isfinite(total) && isfinite(catastrophic) && total >= 0 && catastrophic >= 0) ||
        throw(DomainError((total,catastrophic),"Electronic stopping components must be finite and nonnegative."))
    soft = total - catastrophic
    soft >= 0 || throw(DomainError(soft,
        "Catastrophic electronic stopping exceeds total stopping; this model/cutoff partition is inconsistent."))
    return soft
end
