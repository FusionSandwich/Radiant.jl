const Fixed_Source_Input = Union{
    Surface_Source,
    Volume_Source,
    Boundary_Angular_Current_Source,
    Anisotropic_Volume_Source,
}

"""
    Fixed_Sources

Collection of fixed sources for a Radiant calculation.

Legacy `Surface_Source` and `Volume_Source` objects retain their historical normalization. New
HTS coupling sources carry an explicit `Source_Normalization`; legacy and explicit source
contracts cannot be mixed in one calculation because their normalization semantics differ.
"""
mutable struct Fixed_Sources
    number_of_particles        ::Int64
    particles                  ::Vector{Particle}
    normalization_factor       ::Float64
    sources_names              ::Vector{String}
    sources_list               ::Vector{Source}
    cross_sections             ::Cross_Sections
    geometry                   ::Geometry
    solvers                    ::Solvers
    source_collection          ::Vector{Fixed_Source_Input}
    explicit_normalization     ::Union{Nothing,Source_Normalization}
    projection_receipts        ::Vector{Any}
    is_build                   ::Bool

    function Fixed_Sources(cross_sections,geometry,solvers)
        this = new()
        this.number_of_particles = 0
        this.normalization_factor = 0.0
        this.particles = Vector{Particle}()
        this.sources_names = Vector{String}()
        this.sources_list = Vector{Source}()
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.solvers = solvers
        this.source_collection = Vector{Fixed_Source_Input}()
        this.explicit_normalization = nothing
        this.projection_receipts = Any[]
        this.is_build = false
        return this
    end
end

"""
    add_source(this::Fixed_Sources, fixed_source)

Add a legacy or explicit deterministic source specification. The source is copied so later
projection/build steps cannot mutate the caller's object. Adding a source invalidates a prior
build.
"""
function add_source(this::Fixed_Sources,fixed_source::Fixed_Source_Input)
    push!(this.source_collection,deepcopy(fixed_source))
    this.is_build = false
    return this
end

function _normalizations_compatible(
    first::Source_Normalization,
    second::Source_Normalization;
    rtol::Real = 1.0e-12,
    atol::Real = 0.0,
)
    return first.basis == second.basis &&
           first.time_class == second.time_class &&
           first.time_interval_s == second.time_interval_s &&
           isapprox(first.source_rate_per_s,second.source_rate_per_s;rtol=rtol,atol=atol) &&
           isapprox(first.symmetry_factor,second.symmetry_factor;rtol=rtol,atol=atol)
end

function _validate_explicit_normalizations(source_collection)
    explicit_sources = [source for source in source_collection if _is_explicit_source(source)]
    if isempty(explicit_sources)
        return nothing
    end
    reference = get_source_normalization(explicit_sources[1])
    for source in explicit_sources[2:end]
        candidate = get_source_normalization(source)
        if !_normalizations_compatible(reference,candidate)
            error("Explicit fixed sources have incompatible basis, source rate, symmetry, time interval, or prompt/delayed class.")
        end
    end
    return reference
end

"""
    build(this::Fixed_Sources)

Build and merge source arrays. The operation is idempotent: repeated calls after a successful
build return immediately, and rebuilding after a new source starts from empty derived state.

For explicit HTS coupling sources the legacy normalization divisor is set to one, preserving the
declared source basis (`per_history`, `per_source_particle`, or `per_second`). Physical source-rate
and symmetry scaling are applied later through `Source_Normalization` rather than hidden inside the
transport solve.
"""
function build(this::Fixed_Sources)
    if this.is_build
        return this
    end
    if isempty(this.source_collection)
        error("At least one fixed source is required.")
    end

    has_explicit = any(_is_explicit_source, this.source_collection)
    has_legacy = any(source -> !_is_explicit_source(source), this.source_collection)
    if has_explicit && has_legacy
        error("Legacy and explicit source-normalization contracts cannot be mixed in one Fixed_Sources object.")
    end

    this.number_of_particles = 0
    this.normalization_factor = has_explicit ? 1.0 : 0.0
    empty!(this.particles)
    empty!(this.sources_names)
    empty!(this.sources_list)
    empty!(this.projection_receipts)
    this.explicit_normalization = has_explicit ?
        _validate_explicit_normalizations(this.source_collection) : nothing

    for source_specification in this.source_collection
        fixed_source = deepcopy(source_specification)
        particle = get_particle(fixed_source)
        method = get_method(this.solvers,particle)
        source = Source(particle,this.cross_sections,this.geometry,method)
        receipt = add_source(source,fixed_source)
        if receipt isa Boundary_Projection_Receipt || receipt isa Volume_Projection_Receipt
            push!(this.projection_receipts,receipt)
        end

        index = findfirst(
            candidate -> get_tag(candidate) == get_tag(particle),
            this.particles,
        )
        if isnothing(index)
            this.number_of_particles += 1
            push!(this.particles,particle)
            push!(this.sources_list,source)
        else
            this.sources_list[index] += source
        end

        if !has_explicit
            this.normalization_factor += get_normalization_factor(fixed_source)
        end
    end

    if !has_explicit && (!isfinite(this.normalization_factor) || this.normalization_factor ≤ 0.0)
        error("Legacy fixed-source normalization must be finite and positive.")
    end
    this.is_build = true
    return this
end

"""
    get_source(this::Fixed_Sources, particle::Particle)

Return the merged source for a particle, or an initialized zero source when no fixed source was
provided for that particle.
"""
function get_source(this::Fixed_Sources,particle::Particle)
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        this.particles,
    )
    method = get_method(this.solvers,particle)
    if isnothing(index)
        return Source(particle,this.cross_sections,this.geometry,method)
    end
    return this.sources_list[index]
end

get_particles(this::Fixed_Sources) = this.particles
get_normalization_factor(this::Fixed_Sources) = this.normalization_factor
get_projection_receipts(this::Fixed_Sources) = this.projection_receipts

"""
    get_source_normalization(this::Fixed_Sources)

Return the compatible explicit normalization shared by the new deterministic coupling sources.
An error is raised for legacy source collections.
"""
function get_source_normalization(this::Fixed_Sources)
    if isnothing(this.explicit_normalization)
        error("This Fixed_Sources object uses the legacy normalization contract.")
    end
    return this.explicit_normalization
end
