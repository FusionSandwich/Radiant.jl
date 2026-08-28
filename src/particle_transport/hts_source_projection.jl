"""
    Boundary_Projection_Receipt

Audit receipt for projection of a patch/group/direction boundary source into Radiant's Cartesian
half-range angular basis.
"""
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

"""
    Volume_Projection_Receipt

Audit receipt for projection of a voxel/group/angular volume source into Radiant's volume-source
basis.
"""
struct Volume_Projection_Receipt
    energy_group_map::Vector{Int64}
    direction_map::Vector{Int64}
    target_rate::Vector{Float64}
    projected_rate::Vector{Float64}
    max_relative_error::Float64
    source_hash::String
end

get_particle(this::Boundary_Angular_Current_Source) = this.particle
get_source_normalization(this::Boundary_Angular_Current_Source) = this.normalization

_is_explicit_source(::Boundary_Angular_Current_Source) = true
_is_explicit_source(::Anisotropic_Volume_Source) = true
_is_explicit_source(::Any) = false

function _stored_particle(cross_sections::Cross_Sections,particle::Particle)
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        get_particles(cross_sections),
    )
    if isnothing(index)
        error("Cross-sections do not contain particle tag $(get_tag(particle)).")
    end
    return get_particles(cross_sections)[index]
end

function _energy_group_map(
    source_edges_eV::AbstractVector{<:Real},
    cross_sections::Cross_Sections,
    particle::Particle;
    rtol::Real = 1.0e-10,
    atol_eV::Real = 1.0e-8,
)
    stored_particle = _stored_particle(cross_sections,particle)
    if ismissing(cross_sections.energy_boundaries)
        error("Cross-section energy boundaries are required for deterministic source projection.")
    end
    radiant_edges_eV = 1.0e6 .* Float64.(get_energy_boundaries(cross_sections,stored_particle))
    source_edges = Float64.(source_edges_eV)
    mapping = zeros(Int64,length(source_edges)-1)
    used = Set{Int64}()

    for source_group in eachindex(mapping)
        source_low = min(source_edges[source_group],source_edges[source_group+1])
        source_high = max(source_edges[source_group],source_edges[source_group+1])
        matches = Int64[]
        for radiant_group in 1:length(radiant_edges_eV)-1
            radiant_low = min(radiant_edges_eV[radiant_group],radiant_edges_eV[radiant_group+1])
            radiant_high = max(radiant_edges_eV[radiant_group],radiant_edges_eV[radiant_group+1])
            if isapprox(source_low,radiant_low;rtol=rtol,atol=atol_eV) &&
               isapprox(source_high,radiant_high;rtol=rtol,atol=atol_eV)
                push!(matches,radiant_group)
            end
        end
        if length(matches) != 1
            error("Source energy group $(source_group) does not match exactly one Radiant group.")
        end
        if matches[1] in used
            error("Multiple source groups map to the same Radiant energy group.")
        end
        mapping[source_group] = matches[1]
        push!(used,matches[1])
    end
    return mapping
end

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
        Ω = [Ω,zeros(Float64,length(Ω)),zeros(Float64,length(Ω))]
    end
    directions = hcat(Ω[1],Ω[2],Ω[3])
    return Ω,Float64.(w),directions,Qdims
end

function _direction_map(
    source_directions::AbstractMatrix{<:Real},
    source_weights::AbstractVector{<:Real},
    solver_directions::AbstractMatrix{<:Real},
    solver_weights::AbstractVector{<:Real};
    direction_atol::Real = 1.0e-10,
    weight_rtol::Real = 1.0e-10,
    weight_atol::Real = 1.0e-14,
)
    Nsource = size(source_directions,1)
    if size(source_directions,2) != 3 || length(source_weights) != Nsource
        error("Source directions and weights are inconsistent.")
    end
    mapping = zeros(Int64,Nsource)
    used = Set{Int64}()

    for source_index in 1:Nsource
        best_index = 0
        best_distance = Inf
        for solver_index in 1:size(solver_directions,1)
            distance = norm(
                view(source_directions,source_index,:) .-
                view(solver_directions,solver_index,:),
            )
            if distance < best_distance
                best_distance = distance
                best_index = solver_index
            end
        end
        if best_index == 0 || best_distance > direction_atol
            error("Source direction $(source_index) is not an ordinate of the selected Radiant quadrature.")
        end
        if best_index in used
            error("Multiple source directions map to the same Radiant ordinate.")
        end
        if !isapprox(
            Float64(source_weights[source_index]),
            Float64(solver_weights[best_index]);
            rtol=weight_rtol,
            atol=weight_atol,
        )
            error("Source and Radiant quadrature weights differ for direction $(source_index).")
        end
        mapping[source_index] = best_index
        push!(used,best_index)
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
        else
            output[igroup,coefficient,1] = zeros(Float64,Ny,Nz)
            output[igroup,coefficient,2] = zeros(Float64,Ny,Nz)
            output[igroup,coefficient,3] = zeros(Float64,Nx,Nz)
            output[igroup,coefficient,4] = zeros(Float64,Nx,Nz)
            output[igroup,coefficient,5] = zeros(Float64,Nx,Ny)
            output[igroup,coefficient,6] = zeros(Float64,Nx,Ny)
        end
    end
    return output
end

function _merge_surface_source!(this::Source,new_source)
    old_source = this.surface_sources
    if size(old_source,1) != size(new_source,1) || size(old_source,3) != size(new_source,3)
        error("Surface source arrays have incompatible energy-group or boundary dimensions.")
    end
    Ng = size(old_source,1)
    Np = max(size(old_source,2),size(new_source,2))
    merged = _empty_surface_source_array(Ng,Np,this.geometry)

    for igroup in 1:Ng, coefficient in 1:size(old_source,2), boundary in 1:size(old_source,3)
        merged[igroup,coefficient,boundary] += old_source[igroup,coefficient,boundary]
    end
    for igroup in 1:Ng, coefficient in 1:size(new_source,2), boundary in 1:size(new_source,3)
        merged[igroup,coefficient,boundary] += new_source[igroup,coefficient,boundary]
    end
    this.surface_sources = merged
    return this
end

function _nearest_voxel_index(
    geometry::Geometry,
    axis::String,
    coordinate::Real,
    tolerance_cm::Real,
)
    positions = geometry.voxels_position[axis]
    distances = abs.(positions .- Float64(coordinate))
    index = argmin(distances)
    if distances[index] > tolerance_cm
        error("Patch centroid does not coincide with a Radiant boundary-voxel centroid along $(axis).")
    end
    return index
end

function _cartesian_patch_location(
    geometry::Geometry,
    centroid::AbstractVector{<:Real},
    normal::AbstractVector{<:Real},
    area_cm2::Real;
    normal_atol::Real = 1.0e-8,
    coordinate_atol_cm::Real = 1.0e-10,
    area_rtol::Real = 1.0e-8,
    area_atol_cm2::Real = 1.0e-14,
)
    if get_type(geometry) != "cartesian"
        error("Boundary phase-space projection currently supports Cartesian Radiant geometry only.")
    end
    canonical_normals = (
        (-1.0,0.0,0.0),
        ( 1.0,0.0,0.0),
        (0.0,-1.0,0.0),
        (0.0, 1.0,0.0),
        (0.0,0.0,-1.0),
        (0.0,0.0, 1.0),
    )
    distances = [norm(Float64.(normal) .- collect(candidate)) for candidate in canonical_normals]
    boundary = argmin(distances)
    Ndims = get_dimension(geometry)
    if boundary > 2*Ndims || distances[boundary] > normal_atol
        error("Patch normal is not an outward normal of the selected Cartesian geometry.")
    end

    axis_index = boundary ≤ 2 ? 1 : (boundary ≤ 4 ? 2 : 3)
    axis = ("x","y","z")[axis_index]
    external_boundaries = geometry.voxels_boundaries[axis]
    expected_coordinate = isodd(boundary) ? external_boundaries[1] : external_boundaries[end]
    if abs(Float64(centroid[axis_index])-expected_coordinate) > coordinate_atol_cm
        error("Patch centroid is not on the external boundary selected by its normal.")
    end

    cell = (1,1,1)
    expected_area = 1.0
    if Ndims == 2
        if boundary ≤ 2
            iy = _nearest_voxel_index(geometry,"y",centroid[2],coordinate_atol_cm)
            cell = (iy,1,1)
            expected_area = geometry.voxels_width["y"][iy]
        else
            ix = _nearest_voxel_index(geometry,"x",centroid[1],coordinate_atol_cm)
            cell = (ix,1,1)
            expected_area = geometry.voxels_width["x"][ix]
        end
    elseif Ndims == 3
        if boundary ≤ 2
            iy = _nearest_voxel_index(geometry,"y",centroid[2],coordinate_atol_cm)
            iz = _nearest_voxel_index(geometry,"z",centroid[3],coordinate_atol_cm)
            cell = (iy,iz,1)
            expected_area = geometry.voxels_width["y"][iy] * geometry.voxels_width["z"][iz]
        elseif boundary ≤ 4
            ix = _nearest_voxel_index(geometry,"x",centroid[1],coordinate_atol_cm)
            iz = _nearest_voxel_index(geometry,"z",centroid[3],coordinate_atol_cm)
            cell = (ix,iz,1)
            expected_area = geometry.voxels_width["x"][ix] * geometry.voxels_width["z"][iz]
        else
            ix = _nearest_voxel_index(geometry,"x",centroid[1],coordinate_atol_cm)
            iy = _nearest_voxel_index(geometry,"y",centroid[2],coordinate_atol_cm)
            cell = (ix,iy,1)
            expected_area = geometry.voxels_width["x"][ix] * geometry.voxels_width["y"][iy]
        end
    end

    if !isapprox(Float64(area_cm2),expected_area;rtol=area_rtol,atol=area_atol_cm2)
        error("Patch area does not match the corresponding Radiant boundary-voxel face.")
    end
    return boundary,cell,expected_area
end

function _add_surface_moment!(source_array,igroup,coefficient,boundary,cell,value,Ndims)
    if Ndims == 1
        source_array[igroup,coefficient,boundary] += value
    elseif Ndims == 2
        source_array[igroup,coefficient,boundary][cell[1]] += value
    else
        source_array[igroup,coefficient,boundary][cell[1],cell[2]] += value
    end
end

function _maximum_relative_error(reference::AbstractArray,calculated::AbstractArray,atol::Real)
    if size(reference) != size(calculated)
        error("Cannot compare arrays with different shapes.")
    end
    maximum_error = 0.0
    for index in eachindex(reference)
        denominator = max(abs(Float64(reference[index])),Float64(atol))
        maximum_error = max(
            maximum_error,
            abs(Float64(calculated[index])-Float64(reference[index]))/denominator,
        )
    end
    return maximum_error
end

"""
    project_boundary_source(source, cross_sections, geometry, solver)

Project a `Boundary_Angular_Current_Source` onto the selected Radiant SN quadrature and half-range
surface basis. The initial production-safe implementation requires:

- exact energy-group interval matches;
- source directions that are ordinates of the Radiant quadrature with matching weights;
- one source patch per Cartesian boundary voxel face;
- no duplicate patch-to-face mapping;
- nonnegative reconstructed incoming angular flux;
- patch/group current closure after projection.

The half-range moments are rescaled per patch and group only to remove truncation error in the
protected current. Negative lobes are never clipped.
"""
function project_boundary_source(
    source::Boundary_Angular_Current_Source,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN;
    current_rtol::Real = 1.0e-10,
    current_atol::Real = 1.0e-12,
    positivity_atol::Real = 1.0e-12,
    direction_atol::Real = 1.0e-10,
)
    if !geometry.is_build
        error("Geometry must be built before projecting a boundary source.")
    end
    if get_tag(source.particle) != get_tag(get_particle(solver))
        error("Boundary-source particle does not match the selected solver.")
    end

    group_map = _energy_group_map(source.energy_edges_eV,cross_sections,source.particle)
    stored_particle = _stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    Ω,w,solver_directions,Qdims = _solver_quadrature(solver,geometry)
    direction_map = _direction_map(
        source.directions,
        source.quadrature_weights,
        solver_directions,
        w;
        direction_atol=direction_atol,
    )

    Np,Mn,Dn,nplus_to_n,n_to_nplus,_,_ = surface_angular_polynomial_basis(
        Ω,
        w,
        get_legendre_order(solver),
        get_angular_boltzmann(solver),
        Qdims,
        get_dimension(geometry),
        get_type(geometry),
    )
    projected_source = _empty_surface_source_array(Ng,Np,geometry)
    target_current = zeros(Float64,length(source.patch_ids),Ng)
    projected_current = zeros(Float64,length(source.patch_ids),Ng)
    surface_index = zeros(Int64,length(source.patch_ids))
    surface_cell = Vector{NTuple{3,Int64}}(undef,length(source.patch_ids))
    occupied_cells = Set{Tuple{Int64,NTuple{3,Int64}}}()
    Ndims = get_dimension(geometry)

    for patch in eachindex(source.patch_ids)
        boundary,cell,_ = _cartesian_patch_location(
            geometry,
            view(source.centroids_cm,patch,:),
            view(source.normals,patch,:),
            source.areas_cm2[patch],
        )
        key = (boundary,cell)
        if key in occupied_cells
            error("Multiple source patches map to the same Radiant boundary-voxel face.")
        end
        push!(occupied_cells,key)
        surface_index[patch] = boundary
        surface_cell[patch] = cell
        incoming_solver_indices = nplus_to_n[boundary]
        local_direction_index = n_to_nplus[boundary]

        for source_group in axes(source.angular_flux,2)
            radiant_group = group_map[source_group]
            discrete_flux = zeros(Float64,size(solver_directions,1))
            for source_direction in axes(source.angular_flux,3)
                solver_direction = direction_map[source_direction]
                discrete_flux[solver_direction] = source.angular_flux[patch,source_group,source_direction]
            end

            moments = Dn[boundary] * discrete_flux[incoming_solver_indices]
            moments .*= source.areas_cm2[patch]
            reconstruction = Mn[boundary] * moments
            target = 0.0
            reconstructed = 0.0
            for (half_index,solver_direction) in enumerate(incoming_solver_indices)
                mu = abs(dot(
                    view(solver_directions,solver_direction,:),
                    view(source.normals,patch,:),
                ))
                target += w[solver_direction] * mu *
                          discrete_flux[solver_direction] * source.areas_cm2[patch]
                reconstructed += w[solver_direction] * mu * reconstruction[half_index]
            end

            if target > current_atol
                if reconstructed ≤ 0.0
                    error("Boundary projection has nonpositive reconstructed current for a nonzero source.")
                end
                moments .*= target/reconstructed
                reconstruction = Mn[boundary] * moments
            elseif abs(reconstructed) > current_atol
                error("Boundary projection generated current from a zero-current source group.")
            end

            scale = max(1.0,maximum(abs.(reconstruction)))
            if minimum(reconstruction) < -Float64(positivity_atol)*scale
                error("Boundary projection produced negative incoming angular-flux lobes; increase angular order or change basis.")
            end

            reconstructed = 0.0
            for (half_index,solver_direction) in enumerate(incoming_solver_indices)
                mu = abs(dot(
                    view(solver_directions,solver_direction,:),
                    view(source.normals,patch,:),
                ))
                reconstructed += w[solver_direction] * mu * reconstruction[half_index]
            end
            if !isapprox(reconstructed,target;rtol=current_rtol,atol=current_atol)
                error("Boundary-current closure failed after half-range projection.")
            end

            target_current[patch,radiant_group] += target
            projected_current[patch,radiant_group] += reconstructed
            for coefficient in 1:Np
                _add_surface_moment!(
                    projected_source,
                    radiant_group,
                    coefficient,
                    boundary,
                    cell,
                    moments[coefficient],
                    Ndims,
                )
            end
        end
    end

    receipt = Boundary_Projection_Receipt(
        group_map,
        direction_map,
        surface_index,
        surface_cell,
        target_current,
        projected_current,
        _maximum_relative_error(target_current,projected_current,current_atol),
        source.normalization.source_hash,
    )
    return projected_source,receipt
end

function _linear_voxel_index(geometry::Geometry,voxel_id::Int64)
    Nx = geometry.number_of_voxels["x"]
    Ny = get_dimension(geometry) ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(geometry) ≥ 3 ? geometry.number_of_voxels["z"] : 1
    if voxel_id < 1 || voxel_id > Nx*Ny*Nz
        error("Volume-source voxel identifier is outside the Radiant geometry.")
    end
    cartesian = CartesianIndices((Nx,Ny,Nz))[voxel_id]
    return Tuple(cartesian)
end

"""
    project_volume_source(source, cross_sections, geometry, solver)

Project isotropic, ordinate, or native Radiant moment sources into the SN volume-source array.
Voxel identifiers are one-based Julia linear indices into `(Nx,Ny,Nz)`. Source voxel volumes must
match the built Radiant geometry. Scalar source rate and nonnegative reconstructed angular source
are protected; negative moment coefficients themselves are allowed.
"""
function project_volume_source(
    source::Anisotropic_Volume_Source,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN;
    rate_rtol::Real = 1.0e-10,
    rate_atol::Real = 1.0e-12,
    positivity_atol::Real = 1.0e-12,
    direction_atol::Real = 1.0e-10,
    volume_rtol::Real = 1.0e-10,
    volume_atol_cm3::Real = 1.0e-14,
)
    if !geometry.is_build
        error("Geometry must be built before projecting a volume source.")
    end
    if get_tag(source.particle) != get_tag(get_particle(solver))
        error("Volume-source particle does not match the selected solver.")
    end

    group_map = _energy_group_map(source.energy_edges_eV,cross_sections,source.particle)
    stored_particle = _stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    Ω,w,solver_directions,Qdims = _solver_quadrature(solver,geometry)
    Np,Mn,Dn,_,_ = angular_polynomial_basis(
        Ω,
        w,
        get_legendre_order(solver),
        get_angular_boltzmann(solver),
        Qdims,
    )
    direction_map = source.angular_representation == :ordinates ?
        _direction_map(
            source.directions,
            source.quadrature_weights,
            solver_directions,
            w;
            direction_atol=direction_atol,
        ) : Int64[]

    if source.angular_representation == :moments
        if get(source.provenance,"angular_basis","") != "radiant-volume-moments/v1"
            error("Moment source must declare angular_basis=radiant-volume-moments/v1.")
        end
        if size(source.values,3) != Np
            error("Moment source coefficient count does not match the selected Radiant basis.")
        end
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
        expected_volume = geometry.volume_per_voxel[ix,iy,iz]
        if get_dimension(geometry) == 1
            expected_volume = geometry.volume_per_voxel[ix]
        elseif get_dimension(geometry) == 2
            expected_volume = geometry.volume_per_voxel[ix,iy]
        end
        if !isapprox(
            source.voxel_volumes_cm3[source_voxel],
            expected_volume;
            rtol=volume_rtol,
            atol=volume_atol_cm3,
        )
            error("Volume-source voxel volume does not match the Radiant geometry.")
        end

        for source_group in axes(source.values,2)
            radiant_group = group_map[source_group]
            moments = zeros(Float64,Np)
            target_density = 0.0

            if source.angular_representation == :isotropic
                target_density = source.values[source_voxel,source_group,1]
                moments[1] = target_density
            elseif source.angular_representation == :ordinates
                discrete_source = zeros(Float64,size(solver_directions,1))
                for source_direction in axes(source.values,3)
                    discrete_source[direction_map[source_direction]] =
                        source.values[source_voxel,source_group,source_direction]
                end
                moments .= Dn * discrete_source
                target_density = sum(w .* discrete_source)
            else
                moments .= view(source.values,source_voxel,source_group,:)
                if get(source.provenance,"zeroth_moment_is_angle_integrated","false") != "true"
                    error("Native moment projection requires zeroth_moment_is_angle_integrated=\"true\".")
                end
                index = try
                    parse(Int,get(source.provenance,"zeroth_moment_index",""))
                catch
                    error("Native moment projection requires a valid zeroth_moment_index.")
                end
                if index < 1 || index > Np
                    error("zeroth_moment_index lies outside the Radiant basis.")
                end
                target_density = moments[index]
            end

            reconstruction = Mn * moments
            projected_density = sum(w .* reconstruction)
            if target_density > rate_atol
                if projected_density ≤ 0.0
                    error("Volume projection has nonpositive reconstructed scalar source for a nonzero source.")
                end
                moments .*= target_density/projected_density
                reconstruction = Mn * moments
            elseif abs(projected_density) > rate_atol
                error("Volume projection generated scalar source from a zero source group.")
            end

            scale = max(1.0,maximum(abs.(reconstruction)))
            if minimum(reconstruction) < -Float64(positivity_atol)*scale
                error("Volume projection produced negative angular-source lobes; increase angular order or change basis.")
            end
            projected_density = sum(w .* reconstruction)
            if !isapprox(projected_density,target_density;rtol=rate_rtol,atol=rate_atol)
                error("Volume-source scalar-rate closure failed after angular projection.")
            end

            for coefficient in 1:Np
                projected_source[radiant_group,coefficient,1,ix,iy,iz] += moments[coefficient]
            end
            target_rate[radiant_group] += target_density * expected_volume
            projected_rate[radiant_group] += projected_density * expected_volume
        end
    end

    receipt = Volume_Projection_Receipt(
        group_map,
        direction_map,
        target_rate,
        projected_rate,
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
        source,
        this.cross_sections,
        this.geometry,
        this.solver,
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
        source,
        this.cross_sections,
        this.geometry,
        this.solver,
    )
    if size(projected) != size(this.volume_sources)
        error("Projected volume-source array is incompatible with the initialized Radiant source.")
    end
    this.volume_sources .+= projected
    return receipt
end
