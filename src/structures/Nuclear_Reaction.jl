"""
    IsotopeNuclearReactionDB

Stores the ENDF MF=3 total cross-section curve for one target isotope and one reaction
(MT number).

# Field(s)
- `E::Vector{Float64}` : incident energy grid [eV].
- `S::Vector{Float64}` : tabulated cross-section [barn].
- `NBT::Vector{Int}` : interpolation-region breakpoints.
- `INT::Vector{Int}` : interpolation laws.
"""
struct IsotopeNuclearReactionDB
    E::Vector{Float64}
    S::Vector{Float64}
    NBT::Vector{Int}
    INT::Vector{Int}
end

"""
    KalbachMannTable

Stores the Kalbach-Mann (ENDF LAW=1 LANG=2) energy distribution for one charged-particle
production channel at one target isotope.

# Field(s)
- `E_in::Vector{Float64}` : incident-energy grid [eV].
- `NBT_Ei::Vector{Int}` : interpolation-region breakpoints for the E_in direction.
- `INT_Ei::Vector{Int}` : interpolation laws for the E_in direction.
- `E_out::Vector{Vector{Float64}}` : outgoing-energy grid per incident energy [eV].
- `f::Vector{Vector{Float64}}` : probability-density spectrum per incident energy [eV⁻¹],
  normalized so that ∫f dE_out = 1.
- `r::Vector{Vector{Float64}}` : Kalbach precompound fraction per incident energy
  (dimensionless, 0 ≤ r ≤ 1).
- `LEP::Int` : ENDF interpolation law of the outgoing energy, 1 for a histogram spectrum
  and 2 for a linear-linear one. It sets the quadrature weights with which the spectrum is
  integrated over each outgoing energy bin.
"""
struct KalbachMannTable
    E_in    :: Vector{Float64}
    NBT_Ei  :: Vector{Int}
    INT_Ei  :: Vector{Int}
    E_out   :: Vector{Vector{Float64}}
    f       :: Vector{Vector{Float64}}
    r       :: Vector{Vector{Float64}}
    LEP     :: Int
end

"""
    IsotopeProductionChannelDB

Stores the cross-section and energy-distribution data for one charged-particle production
channel (one MT reaction, one product ZAP) at one target isotope.

# Field(s)
- `mt::Int` : ENDF reaction MT number.
- `ZAP::Int` : ENDF product identifier (1001=p, 1002=d, 1003=t, 2003=He3, 2004=α).
- `AWP::Float64` : product atomic-mass ratio (product mass / neutron mass).
- `E_xs::Vector{Float64}` : MF=3 incident-energy grid [eV].
- `S_xs::Vector{Float64}` : MF=3 total cross-section [barn].
- `NBT_xs::Vector{Int}` : interpolation breakpoints for MF=3.
- `INT_xs::Vector{Int}` : interpolation laws for MF=3.
- `yield_E::Vector{Float64}` : MF=6 yield TAB1 incident-energy grid [eV].
- `yield_y::Vector{Float64}` : MF=6 yield (number of this product per reaction).
- `kalbach::Union{KalbachMannTable,Nothing}` : Kalbach-Mann energy distribution (LAW=1
  LANG=2), or `nothing` when the product distribution is not given as Kalbach-Mann.
"""
struct IsotopeProductionChannelDB
    mt      :: Int
    ZAP     :: Int
    AWP     :: Float64
    E_xs    :: Vector{Float64}
    S_xs    :: Vector{Float64}
    NBT_xs  :: Vector{Int}
    INT_xs  :: Vector{Int}
    yield_E :: Vector{Float64}
    yield_y :: Vector{Float64}
    kalbach :: Union{KalbachMannTable, Nothing}
end

"""
    NuclearReactionENDFDB

Stores the nuclear-reaction ENDF database for one incoming particle type.

# Field(s)
- `mt::Int` : ENDF reaction number (MT) this database was built for.
- `isotopes::Dict{Tuple{Int,Int},IsotopeNuclearReactionDB}` : isotope data indexed by
  `(Z, A)`.
"""
struct NuclearReactionENDFDB
    mt::Int
    isotopes::Dict{Tuple{Int,Int},IsotopeNuclearReactionDB}
end

"""
    Nuclear_Reaction

Structure for multigroup nuclear-reaction cross-sections read from ENDF MF=3 (total XS)
and optionally MF=6 (product energy distributions).

Two interaction modes are supported via `interaction_types`:
- `"A"` (absorption): reads MF=3/MT=`mt` → contributes to Σt and Σa. Represents
  complete removal of the incoming proton by all non-elastic processes.
- `"P"` (equivalent-proton production): reads MF=3 per charged-particle channel and
  MF=6 Kalbach-Mann distributions → contributes to Σsl and Σsₑ. Secondary charged
  particles (p, d, t, He3, α) from the reactions selected by `production_mts` are converted
  to equivalent protons (Salvat & Quesada 2020, §4.2): each product at kinetic energy
  E_b is replaced by a proton at energy E_eq such that their CSDA ranges match, with
  an energy-conservation weight w = E_b/E_eq. This correctly reduces the local energy
  deposition Σe = Σtₑ − Σsₑ by the energy carried away by secondary particles.

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{Type,Type},Vector{String}}` : default
  `Dict((Proton,Proton) => ["A"])` for absorption only; use `["A","P"]` to also generate
  equivalent-proton production cross-sections.
- `library::String = "TENDL2023"` : nuclear-reaction ENDF library name.
- `data_path::String = "../../data"` : root directory for the ENDF data.
- `mt::Int64 = 3` : ENDF MT number for MF=3 total XS (default: 3 = all non-elastic).
- `production_mts::Vector{Int} = Int[]` : MT numbers searched for charged-particle
  products. Empty, the default, lets each evaluation decide: the non-elastic reactions it
  tabulates in both MF=3 and MF=6, reduced by the ENDF sum rules so that no reaction is
  counted twice. A non-empty list overrides that choice.
- `angular_distribution::String = "isotropic"` : how the equivalent protons are emitted,
  either `"isotropic"`, filling the Legendre moment of order zero alone, or
  `"true-angular"`, following the Kalbach-Mann distribution of the product and carrying it
  from the centre of mass to the laboratory.
- `is_neutral_escape::Bool = false` : whether the neutrons and photons emitted by the
  reactions are taken to leave the control volume. They are transported by nothing, so
  their energy is deposited where the reaction happened by default; set to `true`, that
  energy is removed from the local deposition instead.

# Reference(s)
- Salvat & Quesada (2020), NIMB 475, 49–62.
"""
mutable struct Nuclear_Reaction <: Interaction

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
    scattering_model::String
    library::String
    data_path::String
    mt::Int64
    endf_db::Dict{Type,NuclearReactionENDFDB}
    production_mts::Vector{Int}
    is_neutral_escape::Bool
    angular_distribution::String
    production_db::Dict{Type,Dict{Tuple{Int,Int},Vector{IsotopeProductionChannelDB}}}
    production_eq_cache::Dict{Tuple,Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}

    # Constructor(s)
    function Nuclear_Reaction()
        this = new()
        this.name = "Nuclear_Reaction"
        this.interaction_types = Dict((Proton,Proton) => ["A"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.scattering_model = "BTE"
        this.set_library("TENDL2023")
        this.set_data_path("../../data")
        this.set_mt(3)
        this.set_endf_db(Dict{Type,NuclearReactionENDFDB}())
        this.set_production_mts(Int[])
        this.set_is_neutral_escape(false)
        this.set_angular_distribution("isotropic")
        this.production_db = Dict{Type,Dict{Tuple{Int,Int},Vector{IsotopeProductionChannelDB}}}()
        this.production_eq_cache = Dict{Tuple,Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}()
        return this
    end
end

# Method(s)
# Constructor-safe invalidation: setters are also used before all fields exist.
function invalidate_nuclear_caches!(this::Nuclear_Reaction)
    for field in (:endf_db, :production_db, :production_eq_cache)
        isdefined(this,field) && empty!(getfield(this,field))
    end
    return nothing
end

"""
    set_interaction_types(this::Nuclear_Reaction,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

Define the interaction types for nuclear-reaction processes.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : interaction process
  dictionary of the form `(incident particle, outgoing particle) => [type, ...]`:
    - `(Proton,Proton) => ["A"]` : absorption-only — contributes Σt/Σa from MF=3/MT=`mt`.
    - `(Proton,Proton) => ["A","P"]` : absorption + equivalent-proton production — also
      reads per-channel MF=6 Kalbach-Mann data and contributes to Σsl/Σsₑ via the
      Salvat & Quesada (2020) CSDA range-matching method.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> nr = Nuclear_Reaction()
julia> nr.set_interaction_types( Dict((Proton,Proton) => ["A"]) )         # absorption only
julia> nr.set_interaction_types( Dict((Proton,Proton) => ["A","P"]) )     # + production
```
"""
function set_interaction_types(this::Nuclear_Reaction,interaction_types)
    all(v -> length(unique(v)) == length(v),values(interaction_types)) ||
        error("Duplicate nuclear interaction process strings are unsupported.")
    all(kv -> all(t -> t in ("A","P"),kv.second),interaction_types) || error("Unsupported nuclear interaction type.")
    all(kv -> !("A" in kv.second) || kv.first[1] == kv.first[2],interaction_types) ||
        error("Nuclear absorption requires same-species input/output pairs.")
    if any(v -> "P" in v,values(interaction_types))
        all(kv -> !("P" in kv.second) || kv.first == (Proton,Proton),interaction_types) ||
            error("Equivalent-proton production requires Proton-to-Proton.")
        "A" in get(interaction_types,(Proton,Proton),String[]) ||
            error("Proton production requires the Proton-to-Proton absorption/removal channel.")
    end
    this.interaction_types = Dict(k => copy(v) for (k,v) in interaction_types)
    this.incoming_particle = unique([k[1] for k in keys(interaction_types)])
    this.interaction_particles = unique([k[2] for k in keys(interaction_types)])
    invalidate_nuclear_caches!(this)
end

"""
    set_library(this::Nuclear_Reaction,library::String)

Define the nuclear-reaction database library name (e.g., TENDL2023).

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `library::String` : library name.

# Output Argument(s)
N/A
"""
function set_library(this::Nuclear_Reaction,library::String)
    this.library = library
    invalidate_nuclear_caches!(this)
end

"""
    get_library(this::Nuclear_Reaction)

Get the nuclear-reaction database library name.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `library::String` : library name.
"""
function get_library(this::Nuclear_Reaction)
    return this.library
end

"""
    set_data_path(this::Nuclear_Reaction,data_path::String)

Define the root directory containing the nuclear-reaction ENDF data files. Relative
paths are resolved from the location of `Nuclear_Reaction.jl`.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `data_path::String` : data root path. The default `../../data` resolves to the
  package data directory.

# Output Argument(s)
N/A
"""
function set_data_path(this::Nuclear_Reaction,data_path::String)
    path = isabspath(data_path) ? normpath(data_path) : normpath(joinpath(@__DIR__, data_path))
    this.data_path = path
    invalidate_nuclear_caches!(this)
end

"""
    get_data_path(this::Nuclear_Reaction)

Get the root directory containing the nuclear-reaction ENDF data files.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `data_path::String` : data root path.
"""
function get_data_path(this::Nuclear_Reaction)
    return this.data_path
end

"""
    set_mt(this::Nuclear_Reaction,mt::Int64)

Define the ENDF reaction number (MT) to read from MF=3.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `mt::Int64` : ENDF MT number (e.g. 3 for non-elastic total).

# Output Argument(s)
N/A
"""
function set_mt(this::Nuclear_Reaction,mt::Int64)
    this.mt = mt
    invalidate_nuclear_caches!(this)
end

"""
    get_mt(this::Nuclear_Reaction)

Get the ENDF reaction number (MT) read from MF=3.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `mt::Int64` : ENDF MT number.
"""
function get_mt(this::Nuclear_Reaction)
    return this.mt
end

"""
    set_endf_db(this::Nuclear_Reaction,endf_db::Dict{Type,NuclearReactionENDFDB})

Set the nuclear-reaction ENDF database.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `endf_db::Dict{Type,NuclearReactionENDFDB}` : per-particle ENDF database.

# Output Argument(s)
N/A
"""
function set_endf_db(this::Nuclear_Reaction,endf_db::Dict{Type,NuclearReactionENDFDB})
    this.endf_db = endf_db
end

"""
    get_endf_db(this::Nuclear_Reaction)

Get the nuclear-reaction ENDF database.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `endf_db::Dict{Type,NuclearReactionENDFDB}` : per-particle ENDF database.
"""
function get_endf_db(this::Nuclear_Reaction)
    return this.endf_db
end

"""
    set_production_mts(this::Nuclear_Reaction, mts::Vector{Int})

Set the list of ENDF MT reaction numbers for which charged-particle production data are
read from MF=6. Reactions absent from the ENDF file for a given isotope are silently
skipped.

An empty list, the default, leaves the choice to each evaluation: the non-elastic reactions
it tabulates in both MF=3 and MF=6, reduced by the ENDF sum rules so that a reaction
already contained in a coarser one is not counted twice. That is the recommended setting,
since a fixed list cannot know whether an evaluation gives its production explicitly, as
MT=103 or MT=107, or lumped into MT=5.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `mts::Vector{Int}` : list of MT numbers to include, or `Int[]` to decide per evaluation.

# Output Argument(s)
N/A
"""
function set_production_mts(this::Nuclear_Reaction, mts::Vector{Int})
    validate_production_mts(mts)
    this.production_mts = copy(mts)
    invalidate_nuclear_caches!(this)
end

"""
    get_production_mts(this::Nuclear_Reaction)

Get the list of ENDF MT reaction numbers used for charged-particle production.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `mts::Vector{Int}` : list of MT numbers, empty when they are chosen per evaluation.
"""
function get_production_mts(this::Nuclear_Reaction)
    return this.production_mts
end

"""
    initialize(this::Nuclear_Reaction, particles, isotopes::Vector{Tuple{Int,Int}};
    energy_boundaries::Union{Nothing,Vector{Vector{Float64}}}=nothing)

Initializes the nuclear-reaction ENDF databases required by the incoming particles, and
the charged-particle production databases when an equivalent-proton production type `"P"`
is requested. Already initialized particle databases are reused. When the energy group
structure is provided, the tabulated energy range of each isotope is compared with the highest energy
of the corresponding particle and a warning is issued for the evaluations that stop below
it.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `particles` : particle objects used in the simulation.
- `isotopes::Vector{Tuple{Int,Int}}` : target isotope pairs `(Z, A)` required by the materials.
- `energy_boundaries::Union{Nothing,Vector{Vector{Float64}}}` : energy group boundaries of
  each particle [MeV], in the same order as `particles`.

# Output Argument(s)
N/A
"""
function initialize(this::Nuclear_Reaction, particles, isotopes::Vector{Tuple{Int,Int}};
    energy_boundaries::Union{Nothing,Vector{Vector{Float64}}}=nothing)
    endf_db = this.get_endf_db()
    all(v -> length(unique(v)) == length(v),values(this.interaction_types)) ||
        error("Duplicate nuclear interaction process strings are unsupported.")
    all(kv -> !("A" in kv.second) || kv.first[1] == kv.first[2],this.interaction_types) ||
        error("Nuclear absorption requires same-species input/output pairs.")
    validate_production_mts(this.production_mts)
    needs_production = any(v -> "P" ∈ v, values(this.interaction_types))
    isnothing(energy_boundaries) && error("Nuclear admission requires explicit energy boundaries.")
    if needs_production
        "A" in get(this.interaction_types,(Proton,Proton),String[]) ||
            error("Proton production requires its absorption/removal channel.")
        this.angular_distribution == "true-angular" ||
            error("Equivalent-proton production requires explicit true-angular mode and CM data; isotropic production is unsupported.")
        all(kv -> !("P" in kv.second) || (kv.first[1] == Proton && kv.first[2] == Proton), this.interaction_types) ||
            error("Equivalent-proton production supports only Proton-to-Proton; explicit alpha/products are unsupported.")
    end
    for (n, particle) in enumerate(particles)
        ptype = get_type(particle)
        needs_particle_production = any(kv -> kv.first[1] == ptype && "P" in kv.second,
            this.interaction_types)
        if ptype ∈ this.get_in_particles()
            if !haskey(endf_db, ptype)
                endf_db[ptype] = nuclear_reaction_endf(
                    this.get_library(), isotopes, particle, this.get_mt();
                    data_root=this.get_data_path())
            end
            if needs_particle_production && !haskey(this.production_db, ptype)
                this.production_db[ptype] = nuclear_production_endf(
                    this.get_library(), isotopes, particle, this.get_production_mts();
                    data_root=this.get_data_path())
            end
            if !isnothing(energy_boundaries)
                for ((Zi,Ai),iso) in endf_db[ptype].isotopes
                    author_nuclear_table(iso.E,iso.S,"removal isotope ($Zi,$Ai)")
                    author_nuclear_energy(minimum(energy_boundaries[n])*1e6,iso.E,
                        "removal lower energy ($Zi,$Ai)";threshold=iso.S[1] == 0.0)
                end
                if needs_particle_production
                    for ((Zi,Ai),channels) in this.production_db[ptype], channel in channels
                        Emax = maximum(energy_boundaries[n])*1e6
                        author_nuclear_energy(Emax,channel.E_xs,"production cross-section maximum")
                        author_nuclear_energy(Emax,channel.yield_E,"production yield maximum")
                        author_nuclear_energy(Emax,channel.kalbach.E_in,"production incident-spectrum maximum")
                    end
                end
                warn_energy_range(endf_db[ptype], maximum(energy_boundaries[n]),
                    this.get_library(), particle)
            end
        end
    end
    this.set_endf_db(endf_db)
    return nothing
end

"""
    in_distribution(this::Nuclear_Reaction)

Describe the energy discretization method for the incoming particle in the nuclear
reaction interaction.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.
"""
function in_distribution(this::Nuclear_Reaction)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Nuclear_Reaction)

Describe the energy discretization method for the outgoing equivalent proton in the
nuclear-reaction production interaction.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.
"""
function out_distribution(this::Nuclear_Reaction)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    tcs(this::Nuclear_Reaction,Ei::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

Gives the total cross-section for the nuclear-reaction interaction.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `Ei::Float64` : incoming particle energy.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `A::Union{Nothing,Vector{Int64}}` : isotope mass numbers for this element.
- `atpercentA::Union{Nothing,Vector{Float64}}` : isotope atomic fractions.

# Output Argument(s)
- `σt::Float64` : total cross-section.
"""
function tcs(this::Nuclear_Reaction,Ei::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

    if isnothing(A) || isnothing(atpercentA)
        pairs = isotopic_composition(Z)
        A = [p[1] for p in pairs]
        atpercentA = [p[2] for p in pairs]
    end

    ptype = get_type(particle)
    db = this.get_endf_db()[ptype]
    E_eV = Ei * M_E_C2_MEV * 1.0e6 # mₑc² -> MeV -> eV (ENDF native unit)

    σt = 0.0
    for (Ai, atai) in zip(A, atpercentA)
        iso = db.isotopes[(Z,Ai)]
        author_nuclear_energy(E_eV,iso.E,"removal isotope ($Z,$Ai)";threshold=iso.S[1] == 0.0) || continue
        σt += atai * interp_TAB1(E_eV, iso.E, iso.S, iso.NBT, iso.INT)
    end

    return σt * BARN_TO_CM2
end

"""
    set_is_neutral_escape(this::Nuclear_Reaction, is_neutral_escape::Bool)

Set whether the neutrons and photons emitted by the nuclear reactions leave the control
volume.

No interaction type transports them. By default their energy is deposited where the
reaction happened, which overestimates the local dose by whatever they carry out of it: at
100 MeV the neutrons alone take between 5% and 19% of the energy the reaction removes,
more on the heavy elements than on the light ones. Set to `true`, that energy is removed
from the energy-deposition cross-section, which amounts to assuming they escape without
interacting.

The option acts through the `"P"` production feed, so it requires that type to be active;
with `["A"]` alone no feed is computed and nothing is removed.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `is_neutral_escape::Bool` : true to remove the energy of the neutrons and photons from
  the local deposition, false to deposit it.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> nr = Nuclear_Reaction()
julia> nr.set_interaction_types( Dict((Proton,Proton) => ["A","P"]) )
julia> nr.set_is_neutral_escape(true)   # neutrons and photons leave the volume
```
"""
function set_is_neutral_escape(this::Nuclear_Reaction, is_neutral_escape::Bool)
    this.is_neutral_escape = is_neutral_escape
end

"""
    get_is_neutral_escape(this::Nuclear_Reaction)

Get whether the energy of the neutrons and photons is removed from the local deposition.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `is_neutral_escape::Bool` : true when their energy leaves the control volume.
"""
function get_is_neutral_escape(this::Nuclear_Reaction)
    return this.is_neutral_escape
end

"""
    set_angular_distribution(this::Nuclear_Reaction, angular_distribution::String)

Set how the equivalent protons are emitted.

- `"isotropic"` : the emission is spread evenly over the sphere, so only the Legendre
  moment of order zero is filled. The cheapest choice, and a first approximation only.
- `"true-angular"` : the Kalbach-Mann distribution of the product is followed. Its slope
  comes from the Kalbach systematics, the evaluations tabulating the precompound fraction
  alone, and the emission is carried from the centre of mass, where the distribution is
  given, to the laboratory, where the transport happens. The energies move with the angles:
  a product emitted backward comes out slower and one emitted forward faster, which the
  isotropic treatment ignores altogether.

The transformation matters most on the light elements, the centre of mass moving faster the
lighter the target: a 10 MeV alpha emitted sideways in the centre of mass leaves at 26
degrees from the beam on carbon and at 2 degrees on lead.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `angular_distribution::String` : `"isotropic"` or `"true-angular"`.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> nr = Nuclear_Reaction()
julia> nr.set_interaction_types( Dict((Proton,Proton) => ["A","P"]) )
julia> nr.set_angular_distribution("true-angular")
```

# Reference(s)
- Kalbach (1988), Phys. Rev. C 37, 2350.
"""
function set_angular_distribution(this::Nuclear_Reaction, angular_distribution::String)
    if angular_distribution ∉ ("isotropic", "true-angular")
        error("Unknown angular distribution: $(angular_distribution). It should be " *
              "either isotropic or true-angular.")
    end
    this.angular_distribution = angular_distribution
end

"""
    get_angular_distribution(this::Nuclear_Reaction)

Get how the equivalent protons are emitted.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `angular_distribution::String` : `"isotropic"` or `"true-angular"`.
"""
function get_angular_distribution(this::Nuclear_Reaction)
    return this.angular_distribution
end
