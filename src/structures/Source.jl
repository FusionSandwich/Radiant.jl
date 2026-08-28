"""Formatted fixed source arrays for one particle and one solver."""
mutable struct Source
    name                 ::Union{Missing,String}
    particle             ::Particle
    volume_sources       ::Array{Float64}
    surface_sources      ::Array{Union{Array{Float64},Float64}}
    normalization_factor ::Float64
    cross_sections       ::Cross_Sections
    geometry             ::Geometry
    solver               ::Solver

    function Source(
        particle::Particle,
        cross_sections::Cross_Sections,
        geometry::Geometry,
        solver::Solver,
    )
        this = new()
        this.name = missing
        this.particle = particle
        this.normalization_factor = 0.0
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.solver = solver
        initalize_sources(this,cross_sections,geometry,solver)
        return this
    end
end

"""Initialize zero volume and boundary arrays in the selected solver representation."""
function initalize_sources(
    this::Source,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::Solver,
)
    particle = this.particle
    get_tag(particle) in get_tag.(cross_sections.particles) || error(
        "No cross sections are available for particle $(get_tag(particle)).",
    )
    get_tag(particle) == get_tag(get_particle(solver)) || error(
        "Source particle and solver particle differ.",
    )
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        cross_sections.particles,
    )
    Ng = cross_sections.number_of_groups[index]
    Nx = geometry.number_of_voxels["x"]
    Ndims = geometry.dimension
    Ny = Ndims ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = Ndims ≥ 3 ? geometry.number_of_voxels["z"] : 1
    _,_,Nm = get_schemes(solver,geometry,get_is_full_coupling(solver))

    if solver isa SN
        Qdims = get_quadrature_dimension(solver,Ndims)
        Ω,w = quadrature(
            get_quadrature_order(solver),
            get_quadrature_type(solver),
            Ndims,
            Qdims,
        )
        if Ω isa Vector{Float64}
            Ω = [Ω,zeros(Float64,length(Ω)),zeros(Float64,length(Ω))]
        end
        P,_,_,_,_ = angular_polynomial_basis(
            Ω,w,get_legendre_order(solver),get_angular_boltzmann(solver),Qdims,
        )
    elseif solver isa GN
        polynomial_basis = get_polynomial_basis(solver,Ndims)
        Lp = get_legendre_order(solver)
        P = polynomial_basis == "legendre" ? Lp+1 : (Lp+1)^2
    elseif solver isa CP
        P = get_legendre_order(solver)+1
    else
        error("No source representation is available for solver $(typeof(solver)).")
    end

    this.volume_sources = zeros(Float64,Ng,P,Nm[5],Nx,Ny,Nz)
    this.surface_sources = Array{Union{Array{Float64},Float64}}(undef,Ng,1,2*Ndims)
    for igroup in 1:Ng
        if Ndims == 1
            this.surface_sources[igroup,1,1] = 0.0
            this.surface_sources[igroup,1,2] = 0.0
        elseif Ndims == 2
            this.surface_sources[igroup,1,1] = zeros(Float64,Ny)
            this.surface_sources[igroup,1,2] = zeros(Float64,Ny)
            this.surface_sources[igroup,1,3] = zeros(Float64,Nx)
            this.surface_sources[igroup,1,4] = zeros(Float64,Nx)
        elseif Ndims == 3
            this.surface_sources[igroup,1,1] = zeros(Float64,Ny,Nz)
            this.surface_sources[igroup,1,2] = zeros(Float64,Ny,Nz)
            this.surface_sources[igroup,1,3] = zeros(Float64,Nx,Nz)
            this.surface_sources[igroup,1,4] = zeros(Float64,Nx,Nz)
            this.surface_sources[igroup,1,5] = zeros(Float64,Nx,Ny)
            this.surface_sources[igroup,1,6] = zeros(Float64,Nx,Ny)
        else
            error("Only one-, two-, and three-dimensional sources are supported.")
        end
    end
    return this
end

"""Replace the source particle only with a solver-compatible `Particle` object."""
function set_particle(this::Source,particle::Particle)
    get_tag(particle) == get_tag(get_particle(this.solver)) || error(
        "Replacement source particle does not match the configured solver.",
    )
    get_tag(particle) in get_tag.(get_particles(this.cross_sections)) || error(
        "Replacement source particle is absent from the cross-section library.",
    )
    this.particle = particle
    return this
end

function set_particle(::Source,particle::String)
    error(
        "Source.set_particle requires a Particle object, not the string \"$(particle)\".",
    )
end

"""Add a legacy isotropic rectangular volume source."""
function add_source(this::Source,source::Volume_Source)
    get_tag(source.particle) == get_tag(this.particle) || error(
        "Volume-source particle does not match the formatted source.",
    )
    Qv,norm = volume_source(source.particle,source,this.cross_sections,this.geometry)
    size(Qv) == size(this.volume_sources) || error(
        "Legacy volume-source tensor is incompatible with the selected solver.",
    )
    this.volume_sources .+= Qv
    source.normalization_factor += norm
    return nothing
end

"""Add a legacy monodirectional rectangular surface source."""
function add_source(this::Source,source::Surface_Source)
    particle = source.particle
    get_tag(particle) in get_tag.(this.cross_sections.particles) || error(
        "No cross sections are available for particle $(get_tag(particle)).",
    )
    get_tag(particle) == get_tag(this.particle) || error(
        "Surface-source particle does not match the formatted source.",
    )
    Q_new,norm = surface_source(
        particle,source,this.cross_sections,this.geometry,this.solver,
    )
    _merge_legacy_surface_source!(this,Q_new)
    source.normalization_factor += norm
    return nothing
end

function _zero_surface_block(geometry::Geometry,Ng::Int64,Np::Int64)
    Ndims = geometry.dimension
    Nx = geometry.number_of_voxels["x"]
    Ny = Ndims ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = Ndims ≥ 3 ? geometry.number_of_voxels["z"] : 1
    block = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2*Ndims)
    for igroup in 1:Ng, coefficient in 1:Np
        if Ndims == 1
            block[igroup,coefficient,1] = 0.0
            block[igroup,coefficient,2] = 0.0
        elseif Ndims == 2
            block[igroup,coefficient,1] = zeros(Float64,Ny)
            block[igroup,coefficient,2] = zeros(Float64,Ny)
            block[igroup,coefficient,3] = zeros(Float64,Nx)
            block[igroup,coefficient,4] = zeros(Float64,Nx)
        else
            block[igroup,coefficient,1] = zeros(Float64,Ny,Nz)
            block[igroup,coefficient,2] = zeros(Float64,Ny,Nz)
            block[igroup,coefficient,3] = zeros(Float64,Nx,Nz)
            block[igroup,coefficient,4] = zeros(Float64,Nx,Nz)
            block[igroup,coefficient,5] = zeros(Float64,Nx,Ny)
            block[igroup,coefficient,6] = zeros(Float64,Nx,Ny)
        end
    end
    return block
end

function _merge_legacy_surface_source!(this::Source,additional_source)
    existing = this.surface_sources
    size(existing,1) == size(additional_source,1) || error(
        "Surface-source energy-group dimensions differ.",
    )
    size(existing,3) == size(additional_source,3) || error(
        "Surface-source boundary dimensions differ.",
    )
    Ng = size(existing,1)
    Np = max(size(existing,2),size(additional_source,2))
    merged = _zero_surface_block(this.geometry,Ng,Np)
    for igroup in 1:Ng, coefficient in axes(existing,2), boundary in axes(existing,3)
        merged[igroup,coefficient,boundary] += existing[igroup,coefficient,boundary]
    end
    for igroup in 1:Ng, coefficient in axes(additional_source,2), boundary in axes(additional_source,3)
        merged[igroup,coefficient,boundary] += additional_source[igroup,coefficient,boundary]
    end
    this.surface_sources = merged
    return this
end

function add_volume_source(this::Source,source::Array{Float64})
    size(source) == size(this.volume_sources) || error(
        "Volume-source tensor shape does not match the formatted source.",
    )
    this.volume_sources .+= source
    return this
end

get_surface_sources(this::Source) = this.surface_sources
get_volume_sources(this::Source) = this.volume_sources
get_normalization_factor(this::Source) = this.normalization_factor
get_particle(this::Source) = this.particle

function Base.:+(source1::Source,source2::Source)
    get_tag(get_particle(source1)) == get_tag(get_particle(source2)) || error(
        "Sources for different particles cannot be combined.",
    )
    source1.cross_sections === source2.cross_sections || error(
        "Combined sources must reference the same cross-section object.",
    )
    source1.geometry === source2.geometry || error(
        "Combined sources must reference the same geometry object.",
    )
    source1.volume_sources .+= source2.volume_sources
    _merge_legacy_surface_source!(source1,source2.surface_sources)
    source1.normalization_factor += source2.normalization_factor
    return source1
end
