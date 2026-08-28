"""
    Cross_Sections

Container for a coupled multigroup cross-section library.

The structure supports Radiant physics-model generation, FMAC-M, MATXS, and a bounded one-group
custom source used by analytic verification. Per-particle values accept either one value (broadcast
to all currently registered particles) or one value per particle. Setters replace prior state and
are therefore idempotent.
"""
mutable struct Cross_Sections
    source                    ::Union{Missing,String}
    file                      ::Union{Nothing,String}
    number_of_materials       ::Int64
    materials                 ::Vector{Material}
    number_of_particles       ::Int64
    particles                 ::Vector{Particle}
    energy                    ::Union{Missing,Vector{Float64}}
    cutoff                    ::Union{Missing,Vector{Float64}}
    number_of_groups          ::Union{Missing,Vector{Int64}}
    group_structure           ::Union{Missing,Dict{Any,Vector{Float64}}}
    interactions              ::Union{Missing,Vector{Interaction}}
    legendre_order            ::Union{Missing,Int64}
    energy_boundaries         ::Union{Missing,Vector{Vector{Float64}}}
    multigroup_cross_sections ::Union{Missing,Array{Multigroup_Cross_Sections}}
    is_build                  ::Bool
    custom_absorption         ::Vector{Float64}
    custom_scattering         ::Vector{Float64}

    function Cross_Sections()
        this = new()
        this.source = missing
        this.file = nothing
        this.number_of_materials = 0
        this.materials = Material[]
        this.number_of_particles = 0
        this.particles = Particle[]
        this.energy = missing
        this.cutoff = missing
        this.number_of_groups = missing
        this.group_structure = missing
        this.interactions = missing
        this.legendre_order = missing
        this.energy_boundaries = missing
        this.multigroup_cross_sections = missing
        this.is_build = false
        this.custom_absorption = Float64[]
        this.custom_scattering = Float64[]
        return this
    end
end

function _invalidate!(this::Cross_Sections)
    this.is_build = false
    return this
end

function _expand_float_values(values::AbstractVector{<:Real},count::Int64,label::String)
    result = Float64.(values)
    isempty(result) && error("At least one $(label) value is required.")
    any(x -> !isfinite(x),result) && error("$(label) values must be finite.")
    if count > 0
        if length(result) == 1 && count > 1
            result = fill(result[1],count)
        elseif length(result) != count
            error("$(label) requires one value or one value per particle.")
        end
    end
    return result
end

function _expand_integer_values(values::AbstractVector{<:Integer},count::Int64,label::String)
    result = Int64.(values)
    isempty(result) && error("At least one $(label) value is required.")
    if count > 0
        if length(result) == 1 && count > 1
            result = fill(result[1],count)
        elseif length(result) != count
            error("$(label) requires one value or one value per particle.")
        end
    end
    return result
end

function _expand_existing_particle_values!(this::Cross_Sections)
    count = this.number_of_particles
    if !ismissing(this.energy)
        this.energy = _expand_float_values(this.energy,count,"energy")
    end
    if !ismissing(this.cutoff)
        this.cutoff = _expand_float_values(this.cutoff,count,"cutoff")
    end
    if !ismissing(this.number_of_groups)
        this.number_of_groups = _expand_integer_values(this.number_of_groups,count,"group count")
    end
    if !ismissing(this.energy_boundaries)
        if length(this.energy_boundaries) == 1 && count > 1
            this.energy_boundaries = [copy(this.energy_boundaries[1]) for _ in 1:count]
        elseif count > 0 && length(this.energy_boundaries) != count
            error("Energy boundaries require one vector or one vector per particle.")
        end
    end
    return this
end

function _normalize_group_boundaries(boundaries::AbstractVector{<:Real})
    result = Float64.(boundaries)
    length(result) ≥ 2 || error("At least two energy boundaries are required.")
    any(x -> !isfinite(x) || x < 0.0,result) && error(
        "Energy boundaries must be finite and nonnegative.",
    )
    descending = all(result[index] > result[index+1] for index in 1:length(result)-1)
    ascending = all(result[index] < result[index+1] for index in 1:length(result)-1)
    if ascending
        reverse!(result)
    elseif !descending
        error("Energy boundaries must be strictly monotonic.")
    end
    return result
end

function _check_unique_tags(items,kind::String)
    tags = String[]
    for (index,item) in enumerate(items)
        tag = get_tag(item)
        isempty(tag) && error("The $(kind) at position $(index) has an empty tag.")
        tag in tags && error("Duplicate $(kind) tag \"$(tag)\".")
        push!(tags,tag)
    end
    return true
end

function _particle_index(this::Cross_Sections,particle::Particle)
    index = findfirst(candidate -> get_tag(candidate) == get_tag(particle),this.particles)
    isnothing(index) && error("Cross sections do not contain particle $(get_tag(particle)).")
    return index
end

function _require_multigroup_library(this::Cross_Sections)
    ismissing(this.multigroup_cross_sections) && error("Multigroup cross sections are missing.")
    return this.multigroup_cross_sections
end

"""Validate the state required by the selected cross-section source."""
function is_ready_to_build(this::Cross_Sections)
    ismissing(this.source) && error("The cross-section source is not specified.")
    this.number_of_materials > 0 || error("At least one material is required.")
    this.number_of_particles > 0 || error("At least one particle is required.")
    source = lowercase(this.source)

    if source == "fmac-m" || source == "matxs"
        isnothing(this.file) && error("A cross-section file is required for $(source).")
        isempty(this.file) && error("The cross-section file path cannot be empty.")
    elseif source == "physics-models"
        ismissing(this.group_structure) && error("An energy-group structure is required.")
        ismissing(this.interactions) && error("At least one interaction is required.")
        ismissing(this.legendre_order) && error("A Legendre order is required.")
    elseif source == "custom"
        this.number_of_particles == 1 || error(
            "The analytic custom library currently supports exactly one particle.",
        )
        length(this.custom_absorption) == this.number_of_materials || error(
            "Custom absorption requires one value per material.",
        )
        length(this.custom_scattering) == this.number_of_materials || error(
            "Custom scattering requires one value per material.",
        )
        any(x -> !isfinite(x) || x < 0.0,this.custom_absorption) && error(
            "Custom absorption cross sections must be finite and nonnegative.",
        )
        any(x -> !isfinite(x) || x < 0.0,this.custom_scattering) && error(
            "Custom scattering cross sections must be finite and nonnegative.",
        )
    else
        error("Unknown cross-section source: $(this.source).")
    end
    return true
end

"""Build or read the selected multigroup library."""
function build(this::Cross_Sections)
    this.is_build && return this
    _expand_existing_particle_values!(this)
    is_ready_to_build(this)
    source = lowercase(this.source)
    if source == "fmac-m"
        read_fmac_m(this)
    elseif source == "matxs"
        read_matxs(this)
    elseif source == "physics-models"
        try
            generate_cross_sections(this)
        finally
            empty!(cache_radiant[])
        end
    elseif source == "custom"
        custom_cross_sections(this)
    end
    this.is_build = true
    return this
end

function write(this::Cross_Sections,file::String,format::String="fmac-m",binary::Bool=false)
    this.is_build || build(this)
    normalized_format = lowercase(format)
    if normalized_format == "fmac-m"
        write_fmac_m(this,file)
    elseif normalized_format == "matxs"
        write_matxs(this,file,binary)
    else
        error("Unknown cross-section output format: $(format).")
    end
end

function set_source(this::Cross_Sections,source::String)
    normalized = lowercase(source)
    normalized in ("fmac-m","matxs","physics-models","custom") || error(
        "Unknown cross-section source: $(source).",
    )
    this.source = normalized
    return _invalidate!(this)
end

function set_file(this::Cross_Sections,file::String)
    isempty(file) && error("The cross-section file path cannot be empty.")
    this.file = file
    return _invalidate!(this)
end

function set_materials(this::Cross_Sections,material::Material)
    return set_materials(this,Material[material])
end

function set_materials(this::Cross_Sections,materials::AbstractVector{<:Material})
    isempty(materials) && error("At least one material is required.")
    material_list = Material[materials...]
    _check_unique_tags(material_list,"material")
    this.materials = material_list
    this.number_of_materials = length(material_list)
    return _invalidate!(this)
end

function set_particles(this::Cross_Sections,particle::Particle)
    return set_particles(this,Particle[particle])
end

function set_particles(this::Cross_Sections,particles::AbstractVector{<:Particle})
    isempty(particles) && error("At least one particle is required.")
    particle_list = Particle[particles...]
    _check_unique_tags(particle_list,"particle")
    this.particles = particle_list
    this.number_of_particles = length(particle_list)
    _expand_existing_particle_values!(this)
    return _invalidate!(this)
end

function set_energy(this::Cross_Sections,energy::Real)
    return set_energy(this,[energy])
end

function set_energy(this::Cross_Sections,energy::AbstractVector{<:Real})
    result = _expand_float_values(energy,this.number_of_particles,"energy")
    any(x -> x ≤ 0.0,result) && error("Energy values must be positive.")
    this.energy = result
    return _invalidate!(this)
end

function set_cutoff(this::Cross_Sections,cutoff::Real)
    return set_cutoff(this,[cutoff])
end

function set_cutoff(this::Cross_Sections,cutoff::AbstractVector{<:Real})
    result = _expand_float_values(cutoff,this.number_of_particles,"cutoff")
    any(x -> x < 0.0,result) && error("Cutoff energies must be nonnegative.")
    this.cutoff = result
    return _invalidate!(this)
end

function set_number_of_groups(this::Cross_Sections,number_of_groups::Integer)
    return set_number_of_groups(this,[number_of_groups])
end

function set_number_of_groups(this::Cross_Sections,number_of_groups::AbstractVector{<:Integer})
    result = _expand_integer_values(number_of_groups,this.number_of_particles,"group count")
    any(x -> x ≤ 0,result) && error("Each particle must have at least one energy group.")
    this.number_of_groups = result
    return _invalidate!(this)
end

function set_group_structure(this::Cross_Sections,boundaries::AbstractVector{<:Real})
    ismissing(this.group_structure) && (this.group_structure = Dict{Any,Vector{Float64}}())
    this.group_structure["default"] = _normalize_group_boundaries(boundaries)
    return _invalidate!(this)
end

function set_group_structure(
    this::Cross_Sections,
    particle::Particle,
    boundaries::AbstractVector{<:Real},
)
    ismissing(this.group_structure) && (this.group_structure = Dict{Any,Vector{Float64}}())
    this.group_structure[particle] = _normalize_group_boundaries(boundaries)
    return _invalidate!(this)
end

function set_group_structure(this::Cross_Sections,type::String,Ng::Integer,E::Real,Ec::Real)
    return _set_generated_group_structure(this,"default",type,Ng,E,Ec)
end

function set_group_structure(
    this::Cross_Sections,
    particle::Particle,
    type::String,
    Ng::Integer,
    E::Real,
    Ec::Real,
)
    return _set_generated_group_structure(this,particle,type,Ng,E,Ec)
end

function _set_generated_group_structure(
    this::Cross_Sections,
    key,
    type::String,
    Ng::Integer,
    E::Real,
    Ec::Real,
)
    Ng ≥ 1 || error("The number of groups must be positive.")
    isfinite(E) && isfinite(Ec) || error("Energy limits must be finite.")
    E > Ec ≥ 0.0 || error("The upper energy must exceed the nonnegative cutoff.")
    ismissing(this.group_structure) && (this.group_structure = Dict{Any,Vector{Float64}}())
    normalized = lowercase(type)
    if normalized == "linear"
        this.group_structure[key] = linear_energy_group_structure(Int64(Ng),E,Ec)
    elseif normalized == "log"
        Ec > 0.0 || error("A logarithmic group structure requires a positive cutoff.")
        this.group_structure[key] = log_energy_group_structure(Int64(Ng),E,Ec)
    else
        error("Unknown energy-group structure: $(type).")
    end
    return _invalidate!(this)
end

function set_energy_boundaries(
    this::Cross_Sections,
    energy_boundaries::AbstractVector{<:AbstractVector{<:Real}},
)
    result = [_normalize_group_boundaries(boundaries) for boundaries in energy_boundaries]
    isempty(result) && error("At least one particle energy grid is required.")
    if this.number_of_particles > 0
        if length(result) == 1 && this.number_of_particles > 1
            result = [copy(result[1]) for _ in 1:this.number_of_particles]
        elseif length(result) != this.number_of_particles
            error("Energy boundaries require one grid or one grid per particle.")
        end
    end
    this.energy_boundaries = result
    this.number_of_groups = Int64[length(boundaries)-1 for boundaries in result]
    return _invalidate!(this)
end

function set_interactions(this::Cross_Sections,interaction::Interaction)
    return set_interactions(this,Interaction[interaction])
end

function set_interactions(this::Cross_Sections,interactions::AbstractVector{<:Interaction})
    isempty(interactions) && error("At least one interaction is required.")
    this.interactions = Interaction[interactions...]
    return _invalidate!(this)
end

function set_legendre_order(this::Cross_Sections,legendre_order::Integer)
    legendre_order ≥ 0 || error("The Legendre order must be nonnegative.")
    this.legendre_order = Int64(legendre_order)
    return _invalidate!(this)
end

function set_multigroup_cross_sections(
    this::Cross_Sections,
    multigroup_cross_sections::Array{Multigroup_Cross_Sections},
)
    if this.number_of_particles > 0 && this.number_of_materials > 0
        size(multigroup_cross_sections) == (this.number_of_particles,this.number_of_materials) || error(
            "The multigroup-library dimensions must be particle × material.",
        )
    end
    this.multigroup_cross_sections = multigroup_cross_sections
    return this
end

function set_absorption(this::Cross_Sections,Σa::AbstractVector{<:Real})
    values = Float64.(Σa)
    any(x -> !isfinite(x) || x < 0.0,values) && error(
        "Absorption cross sections must be finite and nonnegative.",
    )
    this.custom_absorption = values
    return _invalidate!(this)
end

function set_scattering(this::Cross_Sections,Σs::AbstractVector{<:Real})
    values = Float64.(Σs)
    any(x -> !isfinite(x) || x < 0.0,values) && error(
        "Scattering cross sections must be finite and nonnegative.",
    )
    this.custom_scattering = values
    return _invalidate!(this)
end

get_file(this::Cross_Sections) = this.file
get_number_of_materials(this::Cross_Sections) = this.number_of_materials
get_materials(this::Cross_Sections) = this.materials
get_number_of_particles(this::Cross_Sections) = this.number_of_particles
get_particles(this::Cross_Sections) = this.particles
get_number_of_groups(this::Cross_Sections) = this.number_of_groups
get_energy(this::Cross_Sections) = this.energy
get_cutoff(this::Cross_Sections) = this.cutoff
get_legendre_order(this::Cross_Sections) = this.legendre_order
get_group_structure(this::Cross_Sections) = this.group_structure
get_interactions(this::Cross_Sections) = this.interactions

function get_particle(this::Cross_Sections,tag::String)
    index = findfirst(particle -> get_tag(particle) == tag,this.particles)
    isnothing(index) && error("No particle with tag \"$(tag)\" exists in this library.")
    return this.particles[index]
end

function get_particle(this::Cross_Sections,type::Type)
    index = findfirst(particle -> get_type(particle) == type,this.particles)
    isnothing(index) && error("No particle of type $(type) exists in this library.")
    return this.particles[index]
end

function get_material(this::Cross_Sections,tag::String)
    index = findfirst(material -> get_tag(material) == tag,this.materials)
    isnothing(index) && error("No material with tag \"$(tag)\" exists in this library.")
    return this.materials[index]
end

function get_number_of_groups(this::Cross_Sections,particle::Particle)
    ismissing(this.number_of_groups) && error("Energy-group counts are missing.")
    return this.number_of_groups[_particle_index(this,particle)]
end

function get_energy_boundaries(this::Cross_Sections,particle::Particle)
    ismissing(this.energy_boundaries) && error("Energy boundaries are missing.")
    return this.energy_boundaries[_particle_index(this,particle)]
end

function get_energies(this::Cross_Sections,particle::Particle)
    boundaries = get_energy_boundaries(this,particle)
    return (boundaries[1:end-1] .+ boundaries[2:end]) ./ 2.0
end

function get_energy_width(this::Cross_Sections,particle::Particle)
    return abs.(diff(get_energy_boundaries(this,particle)))
end

function get_absorption(this::Cross_Sections,particle::Particle)
    library = _require_multigroup_library(this)
    particle_index = _particle_index(this,particle)
    Ng = get_number_of_groups(this,particle)
    result = zeros(Float64,Ng,this.number_of_materials)
    for material in 1:this.number_of_materials
        result[:,material] .= get_absorption(library[particle_index,material])
    end
    return result
end

function get_total(this::Cross_Sections,particle::Particle)
    library = _require_multigroup_library(this)
    particle_index = _particle_index(this,particle)
    Ng = get_number_of_groups(this,particle)
    result = zeros(Float64,Ng,this.number_of_materials)
    for material in 1:this.number_of_materials
        result[:,material] .= get_total(library[particle_index,material])
    end
    return result
end

function get_scattering(
    this::Cross_Sections,
    particle_in::Particle,
    particle_out::Particle,
    legendre_order::Integer,
)
    legendre_order ≥ 0 || error("The requested Legendre order must be nonnegative.")
    library = _require_multigroup_library(this)
    incoming_index = _particle_index(this,particle_in)
    outgoing_index = _particle_index(this,particle_out)
    Ngi = get_number_of_groups(this,particle_in)
    Ngf = get_number_of_groups(this,particle_out)
    result = zeros(Float64,this.number_of_materials,Ngi,Ngf,Int64(legendre_order)+1)
    for material in 1:this.number_of_materials
        scattering = get_scattering(library[incoming_index,material],outgoing_index)
        order_count = min(size(scattering,3),Int64(legendre_order)+1)
        result[material,:,:,1:order_count] .= scattering[:,:,1:order_count]
    end
    return result
end

function _material_response(
    this::Cross_Sections,
    particle::Particle,
    getter::Function,
    boundary_value::Bool=false,
)
    library = _require_multigroup_library(this)
    particle_index = _particle_index(this,particle)
    Ng = get_number_of_groups(this,particle) + (boundary_value ? 1 : 0)
    result = zeros(Float64,Ng,this.number_of_materials)
    for material in 1:this.number_of_materials
        result[:,material] .= getter(library[particle_index,material])
    end
    return result
end

get_boundary_stopping_powers(this::Cross_Sections,particle::Particle) =
    _material_response(this,particle,get_boundary_stopping_powers,true)
get_stopping_powers(this::Cross_Sections,particle::Particle) =
    _material_response(this,particle,get_stopping_powers,false)
get_momentum_transfer(this::Cross_Sections,particle::Particle) =
    _material_response(this,particle,get_momentum_transfer,false)
get_energy_deposition(this::Cross_Sections,particle::Particle) =
    _material_response(this,particle,get_energy_deposition,true)
get_charge_deposition(this::Cross_Sections,particle::Particle) =
    _material_response(this,particle,get_charge_deposition,true)

function get_densities(this::Cross_Sections)
    result = zeros(Float64,this.number_of_materials)
    for material in 1:this.number_of_materials
        result[material] = get_density(this.materials[material])
    end
    return result
end
