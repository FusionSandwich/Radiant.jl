"""
    Flux_Per_Particle

Structure used to contain flux information per particle.

"""
mutable struct Flux_Per_Particle

    # Variable(s)
    particle                        ::Particle
    flux                            ::Vector{Array{Float64,6}}
    flux_cutoff                     ::Vector{Array{Float64,5}}
    spectral_radius                 ::Vector{Vector{Float64}}
    boundary_flux                   ::Vector{SN_Boundary_Flux}

    # Constructor(s)
    function Flux_Per_Particle(particle::Particle)

        this = new()
        this.particle = particle
        this.flux = Vector{Array{Float64,6}}()
        this.flux_cutoff = Vector{Array{Float64,5}}()
        this.spectral_radius = Vector{Vector{Float64}}()
        this.boundary_flux = SN_Boundary_Flux[]

        return this
    end
end

# Method(s)
"""
    add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})

Add flux solution to the list of flux solutions.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `flux::Array{Float64,6}` : flux solution.

# Output Argument(s)
N/A

"""
function add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})
    push!(this.flux,flux)
end

"""
    add_flux(this::Flux_Per_Particle,flux_per_particle::Flux_Per_Particle)

Add flux solutions from another Flux_Per_Particle structure.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `flux_per_particle::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
N/A

"""
function add_flux(this::Flux_Per_Particle,flux_per_particle::Flux_Per_Particle)
    if get_tag(this.particle) != get_tag(flux_per_particle.particle) error("Flux particle don't fit.") end
    append!(this.flux,flux_per_particle.flux)
    append!(this.flux_cutoff,flux_per_particle.flux_cutoff)
    append!(this.spectral_radius,flux_per_particle.spectral_radius)
    append!(this.boundary_flux,flux_per_particle.boundary_flux)
end

"""
    add_flux_cutoff(this::Flux_Per_Particle,flux_cutoff::Array{Float64,5})

Add flux at cutoff solutions to the list of flux at cutoff solutions.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `flux_cutoff::Array{Float64,5}` : flux at cutoff solution.

# Output Argument(s)
N/A

"""
function add_flux_cutoff(this::Flux_Per_Particle,flux_cutoff::Array{Float64,5})
    push!(this.flux_cutoff,flux_cutoff)
end

"""
    add_spectral_radius(this::Flux_Per_Particle,spectral_radius::Vector{Float64})

Add the per-energy-group in-group spectral-radius estimate of one generation to the list.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `spectral_radius::Vector{Float64}` : estimated in-group spectral radius per energy group.

# Output Argument(s)
N/A

"""
function add_spectral_radius(this::Flux_Per_Particle,spectral_radius::Vector{Float64})
    push!(this.spectral_radius,spectral_radius)
end

"""
    get_spectral_radius(this::Flux_Per_Particle)

Get the per-energy-group in-group spectral-radius estimate of the last generation solved for the
particle (see the convergence-acceleration section of the documentation for its definition).

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
- `spectral_radius::Vector{Float64}` : estimated in-group spectral radius per energy group.

"""
function get_spectral_radius(this::Flux_Per_Particle)
    if length(this.spectral_radius) == 0 error("No spectral-radius data available.") end
    return this.spectral_radius[end]
end

"""
    get_flux(this::Flux_Per_Particle)

Get the total flux solution for the particle.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
- `flux::Array{Float64}` : flux solution.

"""
function get_flux(this::Flux_Per_Particle)
    return sum(this.flux)
end

"""
    get_flux_cutoff(this::Flux_Per_Particle)

Get the total flux solution at cutoff for the particle.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
- `flux_cutoff::Array{Float64}` : flux solution at cutoff.

"""
function get_flux_cutoff(this::Flux_Per_Particle)
    return sum(this.flux_cutoff)
end

"""Return retained SN boundary generations. Missing generations explicitly reject scoring."""
function get_boundary_flux(this::Flux_Per_Particle)
    length(this.boundary_flux) == length(this.flux) && !isempty(this.boundary_flux) ||
        error("Boundary flux was not retained for every generation; solve with retain_boundary_flux=true.")
    return this.boundary_flux
end

"""Sum outward face/group crossing currents over all retained particle generations."""
function _compatible_boundary_generations(this::Flux_Per_Particle)
    generations = get_boundary_flux(this)
    first = generations[1]
    for data in generations
        data.dimension == first.dimension && data.directions == first.directions &&
            data.weights == first.weights && data.widths == first.widths &&
            data.energy_boundaries == first.energy_boundaries &&
            data.boundary_conditions == first.boundary_conditions &&
            data.orders == first.orders && data.fully_coupled == first.fully_coupled &&
            data.is_csd == first.is_csd && data.basis == first.basis &&
            size.(data.faces) == size.(first.faces) ||
            error("Cannot accumulate boundary generations with different geometry, quadrature or basis metadata.")
    end
    return generations
end

function get_outgoing_current(this::Flux_Per_Particle)
    return sum(get_outgoing_current(data) for data in _compatible_boundary_generations(this))
end

"""Sum represented crossing kinetic energy over compatible captured generations."""
function get_outgoing_energy_current(this::Flux_Per_Particle)
    score = sum(get_outgoing_energy_current(data) for data in _compatible_boundary_generations(this))
    all(isfinite,score) || error("Accumulated boundary energy score overflowed.")
    return score
end

"""Sum void-only represented escaped kinetic energy over compatible generations."""
function get_escaped_energy_current(this::Flux_Per_Particle)
    score = sum(get_escaped_energy_current(data) for data in _compatible_boundary_generations(this))
    all(isfinite,score) || error("Accumulated escaped energy score overflowed.")
    return score
end
