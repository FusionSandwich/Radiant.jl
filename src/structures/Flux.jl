"""
    Flux

Structure used to contain flux information for all particles.

"""
mutable struct Flux

    # Variable(s)
    number_of_particles             ::Int64
    particles                       ::Vector{Particle}
    flux_per_particle               ::Vector{Flux_Per_Particle}

    # Constructor(s)
    function Flux()

        this = new()
        this.number_of_particles = 0
        this.particles = Vector{Particle}()
        this.flux_per_particle = Vector{Flux_Per_Particle}()

        return this
    end
end

# Method(s)
"""
    add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})

Add flux solution for a particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `flux_per_particle::Flux_Per_Particle` : flux solution for a given particle.

# Output Argument(s)
N/A

"""
function add_flux(this::Flux,flux_per_particle::Flux_Per_Particle)
    particle = flux_per_particle.particle
    if get_tag(particle) ∈ get_tag.(this.particles) 
        index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
        this.flux_per_particle[index].add_flux(flux_per_particle)
    else
        this.number_of_particles += 1
        push!(this.particles,particle)
        push!(this.flux_per_particle,flux_per_particle)
    end
end

"""
    get_particles(this::Flux)

Get particle list.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.

# Output Argument(s)
- `particles::Vector{Particle}` : particle list.

"""
function get_particles(this::Flux)
    return this.particles
end

"""
    get_flux(this::Flux,particle::Particle)

Get the total flux solution for a given particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `particle::Particle` : particle

# Output Argument(s)
- `flux_per_particle::Array{Float64}` : flux solution.

"""
function get_flux(this::Flux,particle::Particle)
    if get_tag(particle) ∉ get_tag.(this.particles) error("No data for the specified particle.") end
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    return this.flux_per_particle[index].get_flux()
end

"""
    get_flux_cutoff(this::Flux,particle::Particle)

Get the total flux solution at cutoff for a given particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `particle::Particle` : particle

# Output Argument(s)
- `flux_per_particle::Array{Float64}` : flux solution at cutoff.

"""
function get_flux_cutoff(this::Flux,particle::Particle)
    if get_tag(particle) ∉ get_tag.(this.particles) error("No data for the specified particle.") end
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    return this.flux_per_particle[index].get_flux_cutoff()
end

"""
    get_spectral_radius(this::Flux,particle::Particle)

Get the per-energy-group in-group spectral-radius estimate for a given particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `particle::Particle` : particle.

# Output Argument(s)
- `spectral_radius::Vector{Float64}` : estimated in-group spectral radius per energy group.

"""
function get_spectral_radius(this::Flux,particle::Particle)
    if get_tag(particle) ∉ get_tag.(this.particles) error("No data for the specified particle.") end
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    return this.flux_per_particle[index].get_spectral_radius()
end

"""Get retained outgoing SN boundary generations for the specified particle."""
function get_boundary_flux(this::Flux,particle::Particle)
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    index === nothing && error("No data for the specified particle.")
    return get_boundary_flux(this.flux_per_particle[index])
end

"""Get generation-summed outward crossing current by face and energy group."""
function get_outgoing_current(this::Flux,particle::Particle)
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    index === nothing && error("No data for the specified particle.")
    return get_outgoing_current(this.flux_per_particle[index])
end

"""Get generation-summed represented crossing kinetic energy by face/group."""
function get_outgoing_energy_current(this::Flux,particle::Particle)
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    index === nothing && error("No data for the specified particle.")
    return get_outgoing_energy_current(this.flux_per_particle[index])
end

"""Get generation-summed void-only escaped kinetic energy by face/group."""
function get_escaped_energy_current(this::Flux,particle::Particle)
    index = findfirst(x -> get_tag(x) == get_tag(particle),this.particles)
    index === nothing && error("No data for the specified particle.")
    return get_escaped_energy_current(this.flux_per_particle[index])
end
