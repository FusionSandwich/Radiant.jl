"""Receipt for a conservative boundary-source projection into Radiant's SN basis."""
struct Boundary_Projection_Receipt
    energy_group_map::Vector{Int64}
    direction_map::Vector{Int64}
    surface_index::Vector{Int64}
    surface_cell::Vector{NTuple{3,Int64}}
    target_current::Matrix{Float64}
    projected_current::Matrix{Float64}
    max_relative_error::Float64
    source_hash::String
end

"""Receipt for a conservative volume-source projection into Radiant's SN basis."""
struct Volume_Projection_Receipt
    energy_group_map::Vector{Int64}
    direction_map::Vector{Int64}
    target_rate::Vector{Float64}
    projected_rate::Vector{Float64}
    max_relative_error::Float64
    source_hash::String
end

get_particle(source::Boundary_Angular_Current_Source) = source.particle
get_source_normalization(source::Boundary_Angular_Current_Source) = source.normalization

_is_explicit_source(::Boundary_Angular_Current_Source) = true
_is_explicit_source(::Anisotropic_Volume_Source) = true
_is_explicit_source(::Any) = false

function _stored_particle(cross_sections::Cross_Sections,particle::Particle)
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        get_particles(cross_sections),
    )
    isnothing(index) && error("No cross sections are available for particle $(get_tag(particle)).")
    return get_particles(cross_sections)[index]
end

function _energy_group_map(
    source_edges_eV::AbstractVector{<:Real},
    cross_sections::Cross_Sections,
    particle::Particle;
    rtol::Real=1.0e-10,
    atol_eV::Real=1.0e-8,
)
    ismissing(cross_sections.energy_boundaries) && error(
        "Cross-section energy boundaries are required for source projection.",
    )
    stored_particle = _stored_particle(cross_sections,particle)
    target_edges_eV = 1.0e6 .* Float64.(
        get_energy_boundaries(cross_sections,stored_particle),
    )
    source_edges = Float64.(source_edges_eV)
    mapping = zeros(Int64,length(source_edges)-1)
    used = Set{Int64}()

    for source_group in eachindex(mapping)
        source_low = min(source_edges[source_group],source_edges[source_group+1])
        source_high = max(source_edges[source_group],source_edges[source_group+1])
        matches = Int64[]
        for target_group in 1:length(target_edges_eV)-1
            target_low = min(target_edges_eV[target_group],target_edges_eV[target_group+1])
            target_high = max(target_edges_eV[target_group],target_edges_eV[target_group+1])
            if isapprox(source_low,target_low;rtol=rtol,atol=atol_eV) &&
               isapprox(source_high,target_high;rtol=rtol,atol=atol_eV)
                push!(matches,target_group)
            end
        end
        length(matches) == 1 || error(
            "Source energy group $(source_group) does not match exactly one Radiant group.",
        )
        matches[1] in used && error("Multiple source groups map to one Radiant group.")
        mapping[source_group] = matches[1]
        push!(used,matches[1])
    end
    return mapping
end

"""
    _solver_quadrature(solver, geometry)

Return the native Radiant quadrature together with a three-component unit-vector embedding used by
the coupling schemas. A one-dimensional SN quadrature stores only the transport direction cosine
`μ`; its unused transverse component is chosen deterministically as `sqrt(1-μ^2)`. This embedding
preserves the physically relevant normal cosine while satisfying the source schema's unit-vector
invariant. No azimuthal information is implied for a one-dimensional solve.
"""
function _solver_quadrature(solver::SN,geometry::Geometry)
    Ndims = get_dimension(geometry)
    Qdims = get_quadrature_dimension(solver,Ndims)
    Ω,w = quadrature(
        get_quadrature_order(solver),
        get_quadrature_type(solver),
        Ndims,
        Qdims,
    )
    if Ω isa Vector{Float64}
        μ = Float64.(Ω)
        η = sqrt.(max.(0.0,1.0 .- μ.^2))
        Ω = [μ,η,zeros(Float64,length(μ))]
    end
    return Ω,Float64.(w),hcat(Ω[1],Ω[2],Ω[3]),Qdims
end

function _direction_map(
    source_directions::AbstractMatrix{<:Real},
    source_weights::AbstractVector{<:Real},
    target_directions::AbstractMatrix{<:Real},
    target_weights::AbstractVector{<:Real};
    direction_atol::Real=1.0e-10,
    weight_rtol::Real=1.0e-10,
    weight_atol::Real=1.0e-14,
)
    Nsource = size(source_directions,1)
    size(source_directions,2) == 3 || error("Source directions must be three-vectors.")
    length(source_weights) == Nsource || error("Source directions and weights are inconsistent.")
    mapping = zeros(Int64,Nsource)
    used = Set{Int64}()

    for source_index in 1:Nsource
        distances = [
            norm(
                view(source_directions,source_index,:) .-
                view(target_directions,target_index,:),
            )
            for target_index in axes(target_directions,1)
        ]
        target_index = argmin(distances)
        distances[target_index] ≤ direction_atol || error(
            "Source direction $(source_index) is not an ordinate of the Radiant quadrature.",
        )
        target_index in used && error("Multiple source directions map to one Radiant ordinate.")
        isapprox(
            Float64(source_weights[source_index]),
            Float64(target_weights[target_index]);
            rtol=weight_rtol,
            atol=weight_atol,
        ) || error("Source and Radiant quadrature weights do not match.")
        mapping[source_index] = target_index
        push!(used,target_index)
    end
    return mapping
end

function _empty_surface_source_array(Ng::Int64,Np::Int64,geometry::Geometry)
    Ndims = get_dimension(geometry)
    Nx = geometry.number_of_voxels["x"]
    Ny = Ndims ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = Ndims ≥ 3 ? geometry.number_of_voxels["z"] : 1
    output = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2*Ndims)
    for igroup in 1:Ng, coefficient in 1:Np
        if Ndims == 1
            output[igroup,coefficient,1] = 0.0
            output[igroup,coefficient,2] = 0.0
        elseif Ndims == 2
            output[igroup,coefficient,1] = zeros(Float64,Ny)
            output[igroup,coefficient,2] = zeros(Float64,Ny)
            output[igroup,coefficient,3] = zeros(Float64,Nx)
            output[igroup,coefficient,4] = zeros(Float64,Nx)
        elseif Ndims == 3
            output[igroup,coefficient,1] = zeros(Float64,Ny,Nz)
            output[igroup,coefficient,2] = zeros(Float64,Ny,Nz)
            output[igroup,coefficient,3] = zeros(Float64,Nx,Nz)
            output[igroup,coefficient,4] = zeros(Float64,Nx,Nz)
            output[igroup,coefficient,5] = zeros(Float64,Nx,Ny)
            output[igroup,coefficient,6] = zeros(Float64,Nx,Ny)
        else
            error("Only one-, two-, and three-dimensional Cartesian sources are supported.")
        end
    end
    return output
end

function _merge_surface_source!(source::Source,additional_source)
    existing = source.surface_sources
    size(existing,1) == size(additional_source,1) || error("Surface-source group dimensions differ.")
    size(existing,3) == size(additional_source,3) || error("Surface-source boundary dimensions differ.")
    Ng = size(existing,1)
    Np = max(size(existing,2),size(additional_source,2))
    merged = _empty_surface_source_array(Ng,Np,source.geometry)
    for igroup in 1:Ng, coefficient in axes(existing,2), boundary in axes(existing,3)
        merged[igroup,coefficient,boundary] += existing[igroup,coefficient,boundary]
    end
    for igroup in 1:Ng, coefficient in axes(additional_source,2), boundary in axes(additional_source,3)
        merged[igroup,coefficient,boundary] += additional_source[igroup,coefficient,boundary]
    end
    source.surface_sources = merged
    return source
end

function _nearest_voxel_index(
    geometry::Geometry,
    axis::String,
    coordinate::Real,
    tolerance_cm::Real,
)
    distances = abs.(geometry.voxels_position[axis] .- Float64(coordinate))
    index = argmin(distances)
    distances[index] ≤ tolerance_cm || error(
        "Patch centroid does not coincide with a boundary-cell centroid along $(axis).",
    )
    return index
end

function _cartesian_patch_location(
    geometry::Geometry,
    centroid::AbstractVector{<:Real},
    normal::AbstractVector{<:Real},
    patch_measure::Real;
    normal_atol::Real=1.0e-8,
    coordinate_atol_cm::Real=1.0e-10,
    measure_rtol::Real=1.0e-8,
    measure_atol::Real=1.0e-14,
)
    get_type(geometry) == "cartesian" || error("Only Cartesian projection is implemented.")
    canonical_normals = (
        (-1.0,0.0,0.0),(1.0,0.0,0.0),
        (0.0,-1.0,0.0),(0.0,1.0,0.0),
        (0.0,0.0,-1.0),(0.0,0.0,1.0),
    )
    distances = [norm(Float64.(normal) .- collect(candidate)) for candidate in canonical_normals]
    boundary = argmin(distances)
    Ndims = get_dimension(geometry)
    (boundary ≤ 2*Ndims && distances[boundary] ≤ normal_atol) || error(
        "Patch normal is not an outward Cartesian boundary normal.",
    )

    axis_index = boundary ≤ 2 ? 1 : (boundary ≤ 4 ? 2 : 3)
    axis = ("x","y","z")[axis_index]
    axis_boundaries = geometry.voxels_boundaries[axis]
    expected_coordinate = isodd(boundary) ? axis_boundaries[1] : axis_boundaries[end]
    abs(Float64(centroid[axis_index])-expected_coordinate) ≤ coordinate_atol_cm || error(
        "Patch centroid is not on the boundary selected by its normal.",
    )

    cell = (1,1,1)
    expected_measure = 1.0
    if Ndims == 2
        if boundary ≤ 2
            iy = _nearest_voxel_index(geometry,"y",centroid[2],coordinate_atol_cm)
            cell = (iy,1,1)
            expected_measure = geometry.voxels_width["y"][iy]
        else
            ix = _nearest_voxel_index(geometry,"x",centroid[1],coordinate_atol_cm)
            cell = (ix,1,1)
            expected_measure = geometry.voxels_width["x"][ix]
        end
    elseif Ndims == 3
        if boundary ≤ 2
            iy = _nearest_voxel_index(geometry,"y",centroid[2],coordinate_atol_cm)
            iz = _nearest_voxel_index(geometry,"z",centroid[3],coordinate_atol_cm)
            cell = (iy,iz,1)
            expected_measure = geometry.voxels_width["y"][iy] * geometry.voxels_width["z"][iz]
        elseif boundary ≤ 4
            ix = _nearest_voxel_index(geometry,"x",centroid[1],coordinate_atol_cm)
            iz = _nearest_voxel_index(geometry,"z",centroid[3],coordinate_atol_cm)
            cell = (ix,iz,1)
            expected_measure = geometry.voxels_width["x"][ix] * geometry.voxels_width["z"][iz]
        else
            ix = _nearest_voxel_index(geometry,"x",centroid[1],coordinate_atol_cm)
            iy = _nearest_voxel_index(geometry,"y",centroid[2],coordinate_atol_cm)
            cell = (ix,iy,1)
            expected_measure = geometry.voxels_width["x"][ix] * geometry.voxels_width["y"][iy]
        end
    end

    isapprox(Float64(patch_measure),expected_measure;rtol=measure_rtol,atol=measure_atol) || error(
        "Patch measure does not match the corresponding boundary-cell face.",
    )
    return boundary,cell,expected_measure
end

function _add_surface_moment!(array,igroup,coefficient,boundary,cell,value,Ndims)
    if Ndims == 1
        array[igroup,coefficient,boundary] += value
    elseif Ndims == 2
        array[igroup,coefficient,boundary][cell[1]] += value
    else
        array[igroup,coefficient,boundary][cell[1],cell[2]] += value
    end
end

function _maximum_relative_error(reference,calculated,atol::Real)
    size(reference) == size(calculated) || error("Closure arrays have different shapes.")
    result = 0.0
    for index in eachindex(reference)
        denominator = max(abs(Float64(reference[index])),Float64(atol))
        result = max(result,abs(Float64(calculated[index])-Float64(reference[index]))/denominator)
    end
    return result
end

"""
    project_boundary_source(source, cross_sections, geometry, solver)

Project a patch/group/direction source into Radiant's half-range SN boundary basis. Source energy
intervals and ordinates must map exactly. Each patch must equal one Cartesian boundary-cell face.
The installed moments remain angular-flux densities; patch measure is used only to verify the
extensive incoming current. Negative reconstructed angular flux is rejected, never clipped.
"""
function project_boundary_source(
    source::Boundary_Angular_Current_Source,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN;
    current_rtol::Real=1.0e-10,
    current_atol::Real=1.0e-12,
    positivity_atol::Real=1.0e-12,
    direction_atol::Real=1.0e-10,
)
    geometry.is_build || error("Geometry must be built before boundary-source projection.")
    get_tag(source.particle) == get_tag(get_particle(solver)) || error("Source and solver particles differ.")

    group_map = _energy_group_map(source.energy_edges_eV,cross_sections,source.particle)
    stored_particle = _stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    Ω,w,target_directions,Qdims = _solver_quadrature(solver,geometry)
    direction_map = _direction_map(
        source.directions,source.quadrature_weights,target_directions,w;
        direction_atol=direction_atol,
    )
    Np,Mn,Dn,incoming_indices,_,_,_ = surface_angular_polynomial_basis(
        Ω,w,get_legendre_order(solver),get_angular_boltzmann(solver),Qdims,
        get_dimension(geometry),get_type(geometry),
    )

    projected_source = _empty_surface_source_array(Ng,Np,geometry)
    Npatch = length(source.patch_ids)
    target_current = zeros(Float64,Npatch,Ng)
    projected_current = zeros(Float64,Npatch,Ng)
    surface_index = zeros(Int64,Npatch)
    surface_cell = Vector{NTuple{3,Int64}}(undef,Npatch)
    occupied = Set{Tuple{Int64,NTuple{3,Int64}}}()
    Ndims = get_dimension(geometry)

    for patch in 1:Npatch
        boundary,cell,patch_measure = _cartesian_patch_location(
            geometry,view(source.centroids_cm,patch,:),view(source.normals,patch,:),
            source.areas_cm2[patch],
        )
        key = (boundary,cell)
        key in occupied && error("Multiple source patches map to one boundary-cell face.")
        push!(occupied,key)
        surface_index[patch] = boundary
        surface_cell[patch] = cell
        incoming = incoming_indices[boundary]

        for source_group in axes(source.angular_flux,2)
            target_group = group_map[source_group]
            discrete_flux = zeros(Float64,size(target_directions,1))
            for source_direction in axes(source.angular_flux,3)
                discrete_flux[direction_map[source_direction]] =
                    source.angular_flux[patch,source_group,source_direction]
            end

            moments = Dn[boundary] * discrete_flux[incoming]
            reconstruction = Mn[boundary] * moments
            target_density = 0.0
            reconstructed_density = 0.0
            for (half_index,direction) in enumerate(incoming)
                cosine = abs(dot(
                    view(target_directions,direction,:),view(source.normals,patch,:),
                ))
                target_density += w[direction] * cosine * discrete_flux[direction]
                reconstructed_density += w[direction] * cosine * reconstruction[half_index]
            end

            if target_density > current_atol
                reconstructed_density > 0.0 || error("Projected current is nonpositive for a nonzero source.")
                moments .*= target_density/reconstructed_density
                reconstruction = Mn[boundary] * moments
            elseif abs(reconstructed_density) > current_atol
                error("Projection generated current from a zero-current group.")
            end

            scale = max(1.0,maximum(abs.(reconstruction)))
            minimum(reconstruction) ≥ -Float64(positivity_atol)*scale || error(
                "Boundary projection produced negative angular-flux lobes.",
            )
            reconstructed_density = 0.0
            for (half_index,direction) in enumerate(incoming)
                cosine = abs(dot(
                    view(target_directions,direction,:),view(source.normals,patch,:),
                ))
                reconstructed_density += w[direction] * cosine * reconstruction[half_index]
            end
            isapprox(reconstructed_density,target_density;rtol=current_rtol,atol=current_atol) || error(
                "Boundary-current density failed closure after projection.",
            )

            target_current[patch,target_group] += patch_measure * target_density
            projected_current[patch,target_group] += patch_measure * reconstructed_density
            for coefficient in 1:Np
                _add_surface_moment!(
                    projected_source,target_group,coefficient,boundary,cell,
                    moments[coefficient],Ndims,
                )
            end
        end
    end

    receipt = Boundary_Projection_Receipt(
        group_map,direction_map,surface_index,surface_cell,target_current,projected_current,
        _maximum_relative_error(target_current,projected_current,current_atol),
        source.normalization.source_hash,
    )
    return projected_source,receipt
end

function _linear_voxel_index(geometry::Geometry,voxel_id::Int64)
    Nx = geometry.number_of_voxels["x"]
    Ny = get_dimension(geometry) ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(geometry) ≥ 3 ? geometry.number_of_voxels["z"] : 1
    1 ≤ voxel_id ≤ Nx*Ny*Nz || error("Volume-source voxel identifier is outside the geometry.")
    return Tuple(CartesianIndices((Nx,Ny,Nz))[voxel_id])
end

function _geometry_voxel_volume(geometry::Geometry,ix::Int64,iy::Int64,iz::Int64)
    get_dimension(geometry) == 1 && return geometry.volume_per_voxel[ix]
    get_dimension(geometry) == 2 && return geometry.volume_per_voxel[ix,iy]
    return geometry.volume_per_voxel[ix,iy,iz]
end

"""
    project_volume_source(source, cross_sections, geometry, solver)

Project isotropic, ordinate, or native Radiant moments into the SN volume-source tensor. Voxel
identifiers use one-based Julia linear indexing into `(Nx,Ny,Nz)`. Scalar source rate and
nonnegative reconstructed angular source are protected.
"""
function project_volume_source(
    source::Anisotropic_Volume_Source,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN;
    rate_rtol::Real=1.0e-10,
    rate_atol::Real=1.0e-12,
    positivity_atol::Real=1.0e-12,
    direction_atol::Real=1.0e-10,
    volume_rtol::Real=1.0e-10,
    volume_atol_cm3::Real=1.0e-14,
)
    geometry.is_build || error("Geometry must be built before projecting a volume source.")
    get_tag(source.particle) == get_tag(get_particle(solver)) || error(
        "Volume-source particle does not match the selected solver.",
    )

    group_map = _energy_group_map(source.energy_edges_eV,cross_sections,source.particle)
    stored_particle = _stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    Ω,w,target_directions,Qdims = _solver_quadrature(solver,geometry)
    Np,Mn,Dn,_,_ = angular_polynomial_basis(
        Ω,w,get_legendre_order(solver),get_angular_boltzmann(solver),Qdims,
    )
    direction_map = source.angular_representation == :ordinates ?
        _direction_map(
            source.directions,source.quadrature_weights,target_directions,w;
            direction_atol=direction_atol,
        ) : Int64[]

    if source.angular_representation == :moments
        get(source.provenance,"angular_basis","") == "radiant-volume-moments/v1" || error(
            "Moment source must declare angular_basis=radiant-volume-moments/v1.",
        )
        size(source.values,3) == Np || error(
            "Moment-source coefficient count does not match the selected Radiant basis.",
        )
    end

    _,_,Nm = get_schemes(solver,geometry,get_is_full_coupling(solver))
    Nx = geometry.number_of_voxels["x"]
    Ny = get_dimension(geometry) ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(geometry) ≥ 3 ? geometry.number_of_voxels["z"] : 1
    projected_source = zeros(Float64,Ng,Np,Nm[5],Nx,Ny,Nz)
    target_rate = zeros(Float64,Ng)
    projected_rate = zeros(Float64,Ng)

    for source_voxel in eachindex(source.voxel_ids)
        ix,iy,iz = _linear_voxel_index(geometry,source.voxel_ids[source_voxel])
        expected_volume = _geometry_voxel_volume(geometry,ix,iy,iz)
        isapprox(
            source.voxel_volumes_cm3[source_voxel],expected_volume;
            rtol=volume_rtol,atol=volume_atol_cm3,
        ) || error("Volume-source voxel volume does not match the Radiant geometry.")

        for source_group in axes(source.values,2)
            target_group = group_map[source_group]
            moments = zeros(Float64,Np)
            target_density = 0.0

            if source.angular_representation == :isotropic
                target_density = source.values[source_voxel,source_group,1]
                moments[1] = target_density
            elseif source.angular_representation == :ordinates
                discrete_source = zeros(Float64,size(target_directions,1))
                for source_direction in axes(source.values,3)
                    discrete_source[direction_map[source_direction]] =
                        source.values[source_voxel,source_group,source_direction]
                end
                moments .= Dn * discrete_source
                target_density = sum(w .* discrete_source)
            else
                moments .= view(source.values,source_voxel,source_group,:)
                get(source.provenance,"zeroth_moment_is_angle_integrated","false") == "true" || error(
                    "Native moment projection requires zeroth_moment_is_angle_integrated=\"true\".",
                )
                index = try
                    parse(Int,get(source.provenance,"zeroth_moment_index",""))
                catch
                    error("Native moment projection requires a valid zeroth_moment_index.")
                end
                1 ≤ index ≤ Np || error("zeroth_moment_index lies outside the Radiant basis.")
                target_density = moments[index]
            end

            reconstruction = Mn * moments
            projected_density = sum(w .* reconstruction)
            if target_density > rate_atol
                projected_density > 0.0 || error(
                    "Volume projection has nonpositive reconstructed scalar source for a nonzero source.",
                )
                moments .*= target_density/projected_density
                reconstruction = Mn * moments
            elseif abs(projected_density) > rate_atol
                error("Volume projection generated scalar source from a zero source group.")
            end

            scale = max(1.0,maximum(abs.(reconstruction)))
            minimum(reconstruction) ≥ -Float64(positivity_atol)*scale || error(
                "Volume projection produced negative angular-source lobes; increase angular order or change basis.",
            )
            projected_density = sum(w .* reconstruction)
            isapprox(projected_density,target_density;rtol=rate_rtol,atol=rate_atol) || error(
                "Volume-source scalar-rate closure failed after angular projection.",
            )

            for coefficient in 1:Np
                projected_source[target_group,coefficient,1,ix,iy,iz] += moments[coefficient]
            end
            target_rate[target_group] += target_density * expected_volume
            projected_rate[target_group] += projected_density * expected_volume
        end
    end

    receipt = Volume_Projection_Receipt(
        group_map,direction_map,target_rate,projected_rate,
        _maximum_relative_error(target_rate,projected_rate,rate_atol),
        source.normalization.source_hash,
    )
    return projected_source,receipt
end

"""
    add_source(this::Source, source::Boundary_Angular_Current_Source)

Project and add a conservative tabulated boundary source to an existing Radiant particle source.
Returns a `Boundary_Projection_Receipt`.
"""
function add_source(this::Source,source::Boundary_Angular_Current_Source)
    projected,receipt = project_boundary_source(
        source,this.cross_sections,this.geometry,this.solver,
    )
    _merge_surface_source!(this,projected)
    return receipt
end

"""
    add_source(this::Source, source::Anisotropic_Volume_Source)

Project and add a voxel-resolved angular volume source to an existing Radiant particle source.
Returns a `Volume_Projection_Receipt`.
"""
function add_source(this::Source,source::Anisotropic_Volume_Source)
    projected,receipt = project_volume_source(
        source,this.cross_sections,this.geometry,this.solver,
    )
    size(projected) == size(this.volume_sources) || error(
        "Projected volume-source array is incompatible with the initialized Radiant source.",
    )
    this.volume_sources .+= projected
    return receipt
end