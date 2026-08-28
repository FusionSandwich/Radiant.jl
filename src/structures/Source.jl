"""
    Source

Structure used to describe volume and boundary sources for one transported particle.
"""
mutable struct Source
    name                 ::Union{Missing,String}
    particle             ::Union{Missing,Particle}
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
        initialize_sources(this,cross_sections,geometry,solver)
        return this
    end
end

"""
    initialize_sources(this, cross_sections, geometry, solver)

Initialize zero volume- and surface-source arrays in the angular and spatial basis selected by the
solver. `initalize_sources` remains as a compatibility alias for the historical misspelling.
"""
function initialize_sources(
    this::Source,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::Solver,
)
    particle = this.particle
    if get_tag(particle) ∉ get_tag.(cross_sections.particles)
        error("No cross sections are available for particle $(get_tag(particle)).")
    end
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
            Ω,
            w,
            get_legendre_order(solver),
            get_angular_boltzmann(solver),
            Qdims,
        )
    elseif solver isa GN
        polynomial_basis = get_polynomial_basis(solver,Ndims)
        Lp = get_legendre_order(solver)
        P = polynomial_basis == "legendre" ? Lp+1 : (Lp+1)^2
    elseif solver isa CP
        P = get_legendre_order(solver)+1
    else
        error("No source initialization is available for the selected solver.")
    end

    this.volume_sources = zeros(Float64,Ng,P,Nm[5],Nx,Ny,Nz)
    this.surface_sources = _initialize_legacy_surface_array(Ng,1,Ndims,Nx,Ny,Nz)
    return this
end

initalize_sources(args...) = initialize_sources(args...)

function _initialize_legacy_surface_array(
    Ng::Int64,
    Np::Int64,
    Ndims::Int64,
    Nx::Int64,
    Ny::Int64,
    Nz::Int64,
)
    sources = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2*Ndims)
    for igroup in 1:Ng, coefficient in 1:Np
        if Ndims == 1
            sources[igroup,coefficient,1] = 0.0
            sources[igroup,coefficient,2] = 0.0
        elseif Ndims == 2
            sources[igroup,coefficient,1] = zeros(Float64,Ny)
            sources[igroup,coefficient,2] = zeros(Float64,Ny)
            sources[igroup,coefficient,3] = zeros(Float64,Nx)
            sources[igroup,coefficient,4] = zeros(Float64,Nx)
        elseif Ndims == 3
            sources[igroup,coefficient,1] = zeros(Float64,Ny,Nz)
            sources[igroup,coefficient,2] = zeros(Float64,Ny,Nz)
            sources[igroup,coefficient,3] = zeros(Float64,Nx,Nz)
            sources[igroup,coefficient,4] = zeros(Float64,Nx,Nz)
            sources[igroup,coefficient,5] = zeros(Float64,Nx,Ny)
            sources[igroup,coefficient,6] = zeros(Float64,Nx,Ny)
        else
            error("Cartesian source arrays are available only in one, two, or three dimensions.")
        end
    end
    return sources
end

"""
    set_particle(this::Source, particle::String)

Legacy string setter retained for compatibility.
"""
function set_particle(this::Source,particle::String)
    if lowercase(particle) ∉ ["photons","electrons","positrons"]
        error("Unknown particle type.")
    end
    this.particle = particle
end

"""
    add_source(this::Source, source::Volume_Source)

Add a legacy isotropic, single-group volume source.
"""
function add_source(this::Source,source::Volume_Source)
    particle = source.particle
    Qv,norm = volume_source(particle,source,this.cross_sections,this.geometry)
    this.volume_sources[:,1,1,:,:,:] .+= Qv[:,1,1,:,:,:]
    source.normalization_factor += norm
    return nothing
end

"""
    add_source(this::Source, source::Surface_Source)

Add a legacy monodirectional surface source.
"""
function add_source(this::Source,surface_source_specification::Surface_Source)
    particle = surface_source_specification.particle
    if get_tag(particle) ∉ get_tag.(this.cross_sections.particles)
        error("No cross sections are available for $(get_type(particle)).")
    end
    Ndims = this.geometry.dimension
    Nx = this.geometry.number_of_voxels["x"]
    Ny = Ndims ≥ 2 ? this.geometry.number_of_voxels["y"] : 1
    Nz = Ndims ≥ 3 ? this.geometry.number_of_voxels["z"] : 1
    if get_tag(particle) != get_tag(this.solver.particle)
        error("The source particle does not match the selected solver.")
    end

    Q_new,norm = surface_source(
        particle,
        surface_source_specification,
        this.cross_sections,
        this.geometry,
        this.solver,
    )
    _merge_legacy_surface_source!(this,Q_new,Nx,Ny,Nz)
    surface_source_specification.normalization_factor += norm
    return nothing
end

function _merge_legacy_surface_source!(
    this::Source,
    Q_new,
    Nx::Int64,
    Ny::Int64,
    Nz::Int64,
)
    Q_old = this.surface_sources
    dims_new = size(Q_new)
    dims_old = size(Q_old)
    if dims_new[1] != dims_old[1] || dims_new[3] != dims_old[3]
        error("Surface sources have incompatible energy-group or boundary dimensions.")
    end
    Ng = dims_old[1]
    Ndims = div(dims_old[3],2)
    Np = max(dims_old[2],dims_new[2])
    Q_merged = _initialize_legacy_surface_array(Ng,Np,Ndims,Nx,Ny,Nz)

    for igroup in 1:Ng, coefficient in 1:dims_old[2], boundary in 1:dims_old[3]
        Q_merged[igroup,coefficient,boundary] += Q_old[igroup,coefficient,boundary]
    end
    for igroup in 1:Ng, coefficient in 1:dims_new[2], boundary in 1:dims_new[3]
        Q_merged[igroup,coefficient,boundary] += Q_new[igroup,coefficient,boundary]
    end
    this.surface_sources = Q_merged
    return this
end

"""
    add_volume_source(this::Source, source::Array{Float64})

Replace the volume-source tensor with a preformatted source.
"""
function add_volume_source(this::Source,source::Array{Float64})
    this.volume_sources = source
    return this
end

get_surface_sources(this::Source) = this.surface_sources
get_volume_sources(this::Source) = this.volume_sources
get_normalization_factor(this::Source) = this.normalization_factor
get_particle(this::Source) = this.particle

"""
    Base.:+(source1::Source, source2::Source)

Merge two formatted sources for the same particle.
"""
function Base.:+(source1::Source,source2::Source)
    if get_tag(get_particle(source1)) != get_tag(get_particle(source2))
        error("Sources for different particles cannot be added.")
    end
    if size(source1.volume_sources) != size(source2.volume_sources)
        error("Volume-source arrays are incompatible.")
    end
    source1.volume_sources .+= source2.volume_sources

    Ndims = source1.geometry.dimension
    Nx = source1.geometry.number_of_voxels["x"]
    Ny = Ndims ≥ 2 ? source1.geometry.number_of_voxels["y"] : 1
    Nz = Ndims ≥ 3 ? source1.geometry.number_of_voxels["z"] : 1
    _merge_legacy_surface_source!(source1,source2.surface_sources,Nx,Ny,Nz)
    return source1
end
