const Fixed_Source_Input = Union{
    Surface_Source,
    Volume_Source,
    Boundary_Angular_Current_Source,
    Anisotropic_Volume_Source,
    Energy_Moment_Volume_Source
}

_is_explicit_source(::Energy_Moment_Volume_Source) = true

"""Receipt for the public energy-source projection; zeroth rate retains its existing units."""
struct Energy_Source_Projection_Receipt
    energy_group_map::Vector{Int64}
    direction_map::Vector{Int64}
    native_energy_indices::Vector{Int64}
    target_rate::Vector{Float64}
    projected_rate::Vector{Float64}
    target_energy_moment_rates::Matrix{Float64}
    projected_energy_moment_rates::Matrix{Float64}
    reconstructed_energy_minima::Array{Float64,3}
    reconstruction_roundoff_bounds::Array{Float64,3}
    max_relative_error::Float64
    source_hash::String
end

"""
    project_volume_source(source::Energy_Moment_Volume_Source, cross_sections, geometry, solver)

Project all declared energy modes through the selected angular basis, without angular rescaling.
Only Cartesian SN CSD-family solvers with DG energy order exactly equal to the explicit source
mode count are supported. Pure energy indices are `1:energy_order` for both full tensor coupling
and reduced coupling; all spatial source modes are zero. Unsupported mode counts are rejected,
never discarded or padded. Every represented ordinate energy polynomial is checked over the
complete interval, and every energy moment's scalar rate must close. The receipt reports
group-integrated moment rates in native descending group order, plus unmodified reconstructed
minima and scale-dependent floating-point sign bounds. A value below the negative bound fails;
a zero within the bound has a numerically unresolved sign, not an exact positivity proof.
"""
function project_volume_source(wrapped::Energy_Moment_Volume_Source,
    cross_sections::Cross_Sections, geometry::Geometry, solver::SN;
    rate_rtol::Real=1.0e-10, rate_atol::Real=1.0e-12,
    direction_atol::Real=1.0e-10, volume_rtol::Real=1.0e-10,
    volume_atol_cm3::Real=1.0e-14)
    source = wrapped.source
    edges = source.energy_edges_eV
    length(edges) == size(source.values,2)+1 && all(isfinite,edges) &&
        all(edges .>= 0.0) && all(diff(edges) .> 0.0) ||
        error("Energy-source boundaries must retain their finite, nonnegative, ascending group convention.")
    all(t -> isfinite(t) && t >= 0.0,(rate_rtol,rate_atol,direction_atol,volume_rtol,volume_atol_cm3)) ||
        error("Energy-source projection tolerances must be finite and nonnegative.")
    size(wrapped.energy_moments)[1:3] == size(source.values) &&
        1 <= size(wrapped.energy_moments,4) <= 4 && all(isfinite,wrapped.energy_moments) &&
        wrapped.energy_moments[:,:,:,1] == source.values ||
        error("Energy-source storage no longer satisfies its shape, finite or exact Q0 contract.")
    geometry.is_build || error("Geometry must be built before energy-source projection.")
    get_type(geometry) == "cartesian" || error("Energy sources require Cartesian geometry.")
    get_tag(source.particle) == get_tag(get_particle(solver)) || error("Energy source and solver particles differ.")
    _,is_csd = get_solver_type(solver)
    is_csd || error("Higher energy source representation requires a CSD-family solver.")
    schemes,orders,Nm = get_schemes(solver,geometry,get_is_full_coupling(solver))
    modes = size(wrapped.energy_moments,4)
    schemes[4] == "DG" && orders[4] == modes || error(
        "Source energy mode count must exactly match the selected DG energy order.")
    group_map = _energy_group_map(source.energy_edges_eV,cross_sections,source.particle)
    Ω,w,directions,Qdims = _solver_quadrature(solver,geometry)
    Np,Mn,Dn,_,_ = angular_polynomial_basis(Ω,w,get_legendre_order(solver),get_angular_boltzmann(solver),Qdims)
    direction_map = source.angular_representation == :ordinates ?
        _direction_map(source.directions,source.quadrature_weights,directions,w;direction_atol=direction_atol) : Int64[]
    if source.angular_representation == :moments
        get(source.provenance,"angular_basis","") == "radiant-volume-moments/v1" &&
        get(source.provenance,"zeroth_moment_is_angle_integrated","") == "true" &&
        get(source.provenance,"zeroth_moment_index","") == "1" || error(
            "Energy/angular moment sources require the native angle-integrated zeroth-index-one convention.")
        size(source.values,3) == Np || error("Angular moment count does not match the native basis.")
    end
    stored = _stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored)
    Nx = geometry.number_of_voxels["x"]
    Ny = get_dimension(geometry) >= 2 ? geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(geometry) >= 3 ? geometry.number_of_voxels["z"] : 1
    projected = zeros(Float64,Ng,Np,Nm[5],Nx,Ny,Nz)
    target_rates = zeros(Float64,Ng,modes)
    installed_rates = zeros(Float64,Ng,modes)
    minima = zeros(Float64,length(source.voxel_ids),size(source.values,2),length(w))
    roundoff_bounds = similar(minima)
    for voxel in eachindex(source.voxel_ids)
        ix,iy,iz = _linear_voxel_index(geometry,source.voxel_ids[voxel])
        volume = _geometry_voxel_volume(geometry,ix,iy,iz)
        isapprox(source.voxel_volumes_cm3[voxel],volume;rtol=volume_rtol,atol=volume_atol_cm3) ||
            error("Energy-source voxel volume does not match geometry.")
        for source_group in axes(source.values,2)
            group = group_map[source_group]
            angular = zeros(Float64,Np,modes)
            target = zeros(Float64,modes)
            for mode in 1:modes
                data = view(wrapped.energy_moments,voxel,source_group,:,mode)
                if source.angular_representation == :isotropic
                    angular[1,mode] = data[1]
                    target[mode] = data[1]
                elseif source.angular_representation == :ordinates
                    discrete = zeros(Float64,length(w))
                    discrete[direction_map] .= data
                    angular[:,mode] .= Dn * discrete
                    reconstructed = Mn * view(angular,:,mode)
                    all(isapprox(reconstructed[d],discrete[d];rtol=rate_rtol,atol=rate_atol) for d in eachindex(w)) ||
                        error("Energy source angular representation is not faithfully represented by the selected basis.")
                    target[mode] = sum(w .* discrete)
                else
                    angular[:,mode] .= data
                    target[mode] = data[1]
                end
            end
            reconstructed = Mn * angular
            all(isfinite,reconstructed) || error("Reconstructed energy source is nonfinite.")
            for d in axes(reconstructed,1)
                diagnostic = _source_energy_polynomial_diagnostic(view(reconstructed,d,:))
                minima[voxel,source_group,d] = diagnostic.minimum
                roundoff_bounds[voxel,source_group,d] = diagnostic.roundoff_bound
                diagnostic.minimum >= -diagnostic.roundoff_bound ||
                    error("Reconstructed energy/angular source is negative inside an energy interval.")
            end
            for mode in 1:modes
                installed = sum(w .* view(reconstructed,:,mode))
                isapprox(installed,target[mode];rtol=rate_rtol,atol=rate_atol) ||
                    error("Energy source angular projection failed moment-rate closure.")
                projected[group,:,mode,ix,iy,iz] .+= view(angular,:,mode)
                target_rates[group,mode] += volume*target[mode]
                installed_rates[group,mode] += volume*installed
            end
        end
    end
    all(isfinite,target_rates) && all(isfinite,installed_rates) || error("Energy-source integrated rates overflowed.")
    maximum_relative_error = 0.0
    for index in eachindex(target_rates)
        denominator = max(abs(target_rates[index]),Float64(rate_atol))
        difference = abs(installed_rates[index]-target_rates[index])
        relative = denominator == 0.0 ? (difference == 0.0 ? 0.0 : Inf) : difference/denominator
        maximum_relative_error = max(maximum_relative_error,relative)
    end
    isfinite(maximum_relative_error) || error("Energy-source closure diagnostic is nonfinite.")
    receipt = Energy_Source_Projection_Receipt(group_map,direction_map,Int64.(1:modes),
        target_rates[:,1],installed_rates[:,1],target_rates,installed_rates,minima,roundoff_bounds,
        maximum_relative_error,source.normalization.source_hash)
    return projected,receipt
end

function add_source(this::Source,source::Energy_Moment_Volume_Source)
    this.solver isa SN || error("Energy volume-source projection requires an SN solver.")
    projected,receipt = project_volume_source(source,this.cross_sections,this.geometry,this.solver)
    size(projected) == size(this.volume_sources) || error("Projected energy-source shape is incompatible with native source storage.")
    this.volume_sources .+= projected
    return receipt
end

"""
    Fixed_Sources

Collection of fixed sources for a Radiant calculation.

Legacy `Surface_Source` and `Volume_Source` objects retain their historical normalization. New
HTS coupling sources carry an explicit `Source_Normalization`. The two contracts cannot be mixed
in one calculation because their source-rate semantics differ.
"""
mutable struct Fixed_Sources
    number_of_particles    ::Int64
    particles              ::Vector{Particle}
    normalization_factor   ::Float64
    sources_names          ::Vector{String}
    sources_list           ::Vector{Source}
    cross_sections         ::Cross_Sections
    geometry               ::Geometry
    solvers                ::Solvers
    source_collection      ::Vector{Fixed_Source_Input}
    explicit_normalization ::Union{Nothing,Source_Normalization}
    projection_receipts    ::Vector{Any}
    is_build               ::Bool

    function Fixed_Sources(cross_sections,geometry,solvers)
        this = new()
        this.number_of_particles = 0
        this.particles = Particle[]
        this.normalization_factor = 0.0
        this.sources_names = String[]
        this.sources_list = Source[]
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.solvers = solvers
        this.source_collection = Fixed_Source_Input[]
        this.explicit_normalization = nothing
        this.projection_receipts = Any[]
        this.is_build = false
        return this
    end
end

"""
    add_source(this::Fixed_Sources, fixed_source)

Copy a legacy or explicit source specification into the collection. Adding a source invalidates a
prior build.
"""
function add_source(this::Fixed_Sources,fixed_source::Fixed_Source_Input)
    push!(this.source_collection,deepcopy(fixed_source))
    this.is_build = false
    return this
end

function _normalizations_compatible(
    first::Source_Normalization,
    second::Source_Normalization;
    rtol::Real=1.0e-12,
    atol::Real=0.0,
)
    return first.basis == second.basis &&
           first.time_class == second.time_class &&
           first.time_interval_s == second.time_interval_s &&
           isapprox(first.source_rate_per_s,second.source_rate_per_s;rtol=rtol,atol=atol) &&
           isapprox(first.symmetry_factor,second.symmetry_factor;rtol=rtol,atol=atol)
end

function _validate_explicit_normalizations(source_collection)
    explicit_sources = [source for source in source_collection if _is_explicit_source(source)]
    isempty(explicit_sources) && return nothing
    reference = get_source_normalization(explicit_sources[1])
    for source in explicit_sources[2:end]
        _normalizations_compatible(reference,get_source_normalization(source)) || error(
            "Explicit sources have incompatible basis, rate, symmetry, time interval, or time class.",
        )
    end
    return reference
end

"""
    build(this::Fixed_Sources)

Build and merge source arrays. The operation is idempotent. Rebuilding after a newly added source
starts from empty derived state.

For explicit HTS sources the legacy normalization divisor is exactly one, preserving the declared
per-history, per-source-particle, or per-second basis. Physical source-rate and symmetry scaling
remain in `Source_Normalization`.
"""
function build(this::Fixed_Sources)
    this.is_build && return this
    isempty(this.source_collection) && error("At least one fixed source is required.")

    has_explicit = any(_is_explicit_source,this.source_collection)
    has_legacy = any(source -> !_is_explicit_source(source),this.source_collection)
    has_explicit && has_legacy && error(
        "Legacy and explicit source-normalization contracts cannot be mixed.",
    )

    this.number_of_particles = 0
    this.normalization_factor = has_explicit ? 1.0 : 0.0
    empty!(this.particles)
    empty!(this.sources_names)
    empty!(this.sources_list)
    empty!(this.projection_receipts)
    if has_explicit
        this.explicit_normalization = _validate_explicit_normalizations(this.source_collection)
    else
        this.explicit_normalization = nothing
    end

    for source_specification in this.source_collection
        fixed_source = deepcopy(source_specification)
        particle = get_particle(fixed_source)
        method = get_method(this.solvers,particle)
        formatted_source = Source(particle,this.cross_sections,this.geometry,method)
        receipt = add_source(formatted_source,fixed_source)
        if receipt isa Boundary_Projection_Receipt || receipt isa Volume_Projection_Receipt || receipt isa Energy_Source_Projection_Receipt
            push!(this.projection_receipts,receipt)
        end

        index = findfirst(
            candidate -> get_tag(candidate) == get_tag(particle),
            this.particles,
        )
        if isnothing(index)
            this.number_of_particles += 1
            push!(this.particles,particle)
            push!(this.sources_list,formatted_source)
        else
            this.sources_list[index] += formatted_source
        end

        if !has_explicit
            this.normalization_factor += get_normalization_factor(fixed_source)
        end
    end

    if !has_explicit
        isfinite(this.normalization_factor) && this.normalization_factor > 0.0 || error(
            "Legacy fixed-source normalization must be finite and positive.",
        )
    end
    this.is_build = true
    return this
end

"""Return the merged source for a particle, or a compatible initialized zero source."""
function get_source(this::Fixed_Sources,particle::Particle)
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        this.particles,
    )
    method = get_method(this.solvers,particle)
    return isnothing(index) ?
        Source(particle,this.cross_sections,this.geometry,method) :
        this.sources_list[index]
end

get_particles(this::Fixed_Sources) = this.particles
get_normalization_factor(this::Fixed_Sources) = this.normalization_factor
get_projection_receipts(this::Fixed_Sources) = this.projection_receipts

"""Return the shared explicit normalization; legacy source collections raise an error."""
function get_source_normalization(this::Fixed_Sources)
    isnothing(this.explicit_normalization) && error(
        "This Fixed_Sources object uses the legacy normalization contract.",
    )
    return this.explicit_normalization
end
