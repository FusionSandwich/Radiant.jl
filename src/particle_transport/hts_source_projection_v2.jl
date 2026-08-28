"""Receipt for one boundary phase-space projection into Radiant source moments."""
struct Boundary_Projection_Receipt
    energy_group_map::Vector{Int64}
    direction_map::Vector{Int64}
    surface_index::Vector{Int64}
    boundary_cell_index::Vector{Int64}
    target_current::Matrix{Float64}
    projected_current::Matrix{Float64}
    max_relative_error::Float64
    source_hash::String
end

"""Receipt for one anisotropic volume-source projection into Radiant moments."""
struct Volume_Projection_Receipt
    energy_group_map::Vector{Int64}
    direction_map::Vector{Int64}
    target_rate::Vector{Float64}
    projected_rate::Vector{Float64}
    max_relative_error::Float64
    source_hash::String
end

_is_explicit_source(::Boundary_Angular_Current_Source) = true
_is_explicit_source(::Anisotropic_Volume_Source) = true
_is_explicit_source(::Any) = false

function _projection_maximum_relative_error(reference,calculated,atol::Real)
    size(reference) == size(calculated) || error("Projection closure arrays have different shapes.")
    error_value = 0.0
    for index in eachindex(reference)
        denominator = max(abs(reference[index]),Float64(atol))
        error_value = max(error_value,abs(calculated[index]-reference[index])/denominator)
    end
    return error_value
end

function _projection_stored_particle(
    cross_sections::Cross_Sections,
    particle::Particle,
)
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        get_particles(cross_sections),
    )
    isnothing(index) && error(
        "Particle $(get_tag(particle)) is absent from the cross-section library.",
    )
    return get_particles(cross_sections)[index]
end

function _projection_energy_group_map(
    source_edges_eV::AbstractVector{<:Real},
    cross_sections::Cross_Sections,
    particle::Particle;
    rtol::Real=1.0e-10,
    atol_eV::Real=1.0e-8,
)
    stored_particle = _projection_stored_particle(cross_sections,particle)
    radiant_edges_eV = 1.0e6 .* get_energy_boundaries(cross_sections,stored_particle)
    source_edges = Float64.(source_edges_eV)
    source_groups = length(source_edges)-1
    radiant_groups = length(radiant_edges_eV)-1
    source_groups > 0 || error("At least one source energy group is required.")

    mapping = zeros(Int64,source_groups)
    for source_group in 1:source_groups
        source_low = min(source_edges[source_group],source_edges[source_group+1])
        source_high = max(source_edges[source_group],source_edges[source_group+1])
        matches = Int64[]
        for radiant_group in 1:radiant_groups
            radiant_low = min(radiant_edges_eV[radiant_group],radiant_edges_eV[radiant_group+1])
            radiant_high = max(radiant_edges_eV[radiant_group],radiant_edges_eV[radiant_group+1])
            if isapprox(source_low,radiant_low;rtol=rtol,atol=atol_eV) &&
               isapprox(source_high,radiant_high;rtol=rtol,atol=atol_eV)
                push!(matches,radiant_group)
            end
        end
        length(matches) == 1 || error(
            "Source energy group $(source_group) does not map uniquely to a Radiant group.",
        )
        mapping[source_group] = matches[1]
    end
    length(unique(mapping)) == length(mapping) || error(
        "Multiple source groups map to one Radiant group.",
    )
    return mapping
end

function _projection_solver_quadrature(solver::SN,geometry::Geometry)
    dimension = get_dimension(geometry)
    quadrature_dimension = get_quadrature_dimension(solver,dimension)
    directions_raw,weights_raw = quadrature(
        get_quadrature_order(solver),
        get_quadrature_type(solver),
        dimension,
        quadrature_dimension,
    )
    directions = if isa(directions_raw,Vector{Float64})
        [
            Float64.(directions_raw),
            zeros(Float64,length(directions_raw)),
            zeros(Float64,length(directions_raw)),
        ]
    else
        [Float64.(component) for component in directions_raw]
    end
    matrix = hcat(directions[1],directions[2],directions[3])
    return (directions,Float64.(weights_raw),matrix,quadrature_dimension)
end

function _projection_direction_map(
    source_directions::AbstractMatrix{<:Real},
    source_weights::AbstractVector{<:Real},
    solver_directions::AbstractMatrix{<:Real},
    solver_weights::AbstractVector{<:Real};
    direction_atol::Real=1.0e-10,
    weight_rtol::Real=1.0e-10,
    weight_atol::Real=1.0e-14,
)
    size(source_directions,2) == 3 || error("Source directions must have three components.")
    size(solver_directions,2) == 3 || error("Solver directions must have three components.")
    size(source_directions,1) == length(source_weights) || error(
        "Source direction and weight counts differ.",
    )
    size(solver_directions,1) == length(solver_weights) || error(
        "Solver direction and weight counts differ.",
    )

    mapping = zeros(Int64,size(source_directions,1))
    used = Set{Int64}()
    for source_direction in axes(source_directions,1)
        candidates = Int64[]
        for solver_direction in axes(solver_directions,1)
            distance = norm(
                view(source_directions,source_direction,:) .-
                view(solver_directions,solver_direction,:),
            )
            if distance ≤ direction_atol && isapprox(
                source_weights[source_direction],solver_weights[solver_direction];
                rtol=weight_rtol,atol=weight_atol,
            )
                push!(candidates,solver_direction)
            end
        end
        length(candidates) == 1 || error(
            "Source ordinate $(source_direction) does not map uniquely to the Radiant quadrature.",
        )
        candidates[1] in used && error("Two source ordinates map to one Radiant ordinate.")
        mapping[source_direction] = candidates[1]
        push!(used,candidates[1])
    end
    return mapping
end

function _projection_surface_name(surface_index::Int64)
    names = ("X-","X+","Y-","Y+","Z-","Z+")
    1 ≤ surface_index ≤ length(names) || error("Unknown Cartesian surface index.")
    return names[surface_index]
end

function _projection_axis_data(geometry::Geometry,axis::String)
    return (
        geometry.voxels_position[axis],
        geometry.voxels_width[axis],
        geometry.voxels_boundaries[axis],
        geometry.number_of_voxels[axis],
    )
end

function _projection_find_cell(
    positions::Vector{Float64},
    value::Float64;
    atol::Real=1.0e-10,
)
    candidates = findall(position -> isapprox(position,value;rtol=0.0,atol=atol),positions)
    length(candidates) == 1 || error(
        "Patch centroid does not map uniquely to a boundary-cell centroid.",
    )
    return Int64(candidates[1])
end

function _projection_patch_mapping(
    source::Boundary_Angular_Current_Source,
    geometry::Geometry;
    coordinate_atol_cm::Real=1.0e-10,
    area_rtol::Real=1.0e-10,
    area_atol_cm2::Real=1.0e-14,
    normal_atol::Real=1.0e-10,
)
    get_type(geometry) == "cartesian" || error(
        "Boundary source projection currently requires Cartesian geometry.",
    )
    dimension = get_dimension(geometry)
    surface_indices = zeros(Int64,length(source.patch_ids))
    cell_indices = zeros(Int64,length(source.patch_ids))
    occupied = Set{Tuple{Int64,Int64}}()

    x,x_width,x_boundaries,Nx = _projection_axis_data(geometry,"x")
    y = dimension ≥ 2 ? geometry.voxels_position["y"] : [0.0]
    y_width = dimension ≥ 2 ? geometry.voxels_width["y"] : [1.0]
    Ny = dimension ≥ 2 ? geometry.number_of_voxels["y"] : 1
    z = dimension ≥ 3 ? geometry.voxels_position["z"] : [0.0]
    z_width = dimension ≥ 3 ? geometry.voxels_width["z"] : [1.0]
    Nz = dimension ≥ 3 ? geometry.number_of_voxels["z"] : 1

    for patch in eachindex(source.patch_ids)
        normal = view(source.normals,patch,:)
        centroid = view(source.centroids_cm,patch,:)
        surface = 0
        cell = 0
        expected_area = 0.0

        if norm(normal .- [-1.0,0.0,0.0]) ≤ normal_atol
            isapprox(centroid[1],x_boundaries[1];rtol=0.0,atol=coordinate_atol_cm) || error(
                "X- patch is not on the X- boundary.",
            )
            iy = dimension ≥ 2 ? _projection_find_cell(y,centroid[2];atol=coordinate_atol_cm) : 1
            iz = dimension ≥ 3 ? _projection_find_cell(z,centroid[3];atol=coordinate_atol_cm) : 1
            surface = 1
            cell = LinearIndices((Ny,Nz))[iy,iz]
            expected_area = y_width[iy]*z_width[iz]
        elseif norm(normal .- [1.0,0.0,0.0]) ≤ normal_atol
            isapprox(centroid[1],x_boundaries[end];rtol=0.0,atol=coordinate_atol_cm) || error(
                "X+ patch is not on the X+ boundary.",
            )
            iy = dimension ≥ 2 ? _projection_find_cell(y,centroid[2];atol=coordinate_atol_cm) : 1
            iz = dimension ≥ 3 ? _projection_find_cell(z,centroid[3];atol=coordinate_atol_cm) : 1
            surface = 2
            cell = LinearIndices((Ny,Nz))[iy,iz]
            expected_area = y_width[iy]*z_width[iz]
        elseif dimension ≥ 2 && norm(normal .- [0.0,-1.0,0.0]) ≤ normal_atol
            y_boundaries = geometry.voxels_boundaries["y"]
            isapprox(centroid[2],y_boundaries[1];rtol=0.0,atol=coordinate_atol_cm) || error(
                "Y- patch is not on the Y- boundary.",
            )
            ix = _projection_find_cell(x,centroid[1];atol=coordinate_atol_cm)
            iz = dimension ≥ 3 ? _projection_find_cell(z,centroid[3];atol=coordinate_atol_cm) : 1
            surface = 3
            cell = LinearIndices((Nx,Nz))[ix,iz]
            expected_area = x_width[ix]*z_width[iz]
        elseif dimension ≥ 2 && norm(normal .- [0.0,1.0,0.0]) ≤ normal_atol
            y_boundaries = geometry.voxels_boundaries["y"]
            isapprox(centroid[2],y_boundaries[end];rtol=0.0,atol=coordinate_atol_cm) || error(
                "Y+ patch is not on the Y+ boundary.",
            )
            ix = _projection_find_cell(x,centroid[1];atol=coordinate_atol_cm)
            iz = dimension ≥ 3 ? _projection_find_cell(z,centroid[3];atol=coordinate_atol_cm) : 1
            surface = 4
            cell = LinearIndices((Nx,Nz))[ix,iz]
            expected_area = x_width[ix]*z_width[iz]
        elseif dimension ≥ 3 && norm(normal .- [0.0,0.0,-1.0]) ≤ normal_atol
            z_boundaries = geometry.voxels_boundaries["z"]
            isapprox(centroid[3],z_boundaries[1];rtol=0.0,atol=coordinate_atol_cm) || error(
                "Z- patch is not on the Z- boundary.",
            )
            ix = _projection_find_cell(x,centroid[1];atol=coordinate_atol_cm)
            iy = _projection_find_cell(y,centroid[2];atol=coordinate_atol_cm)
            surface = 5
            cell = LinearIndices((Nx,Ny))[ix,iy]
            expected_area = x_width[ix]*y_width[iy]
        elseif dimension ≥ 3 && norm(normal .- [0.0,0.0,1.0]) ≤ normal_atol
            z_boundaries = geometry.voxels_boundaries["z"]
            isapprox(centroid[3],z_boundaries[end];rtol=0.0,atol=coordinate_atol_cm) || error(
                "Z+ patch is not on the Z+ boundary.",
            )
            ix = _projection_find_cell(x,centroid[1];atol=coordinate_atol_cm)
            iy = _projection_find_cell(y,centroid[2];atol=coordinate_atol_cm)
            surface = 6
            cell = LinearIndices((Nx,Ny))[ix,iy]
            expected_area = x_width[ix]*y_width[iy]
        else
            error("Patch normal is not an outward Cartesian boundary normal.")
        end

        isapprox(source.areas_cm2[patch],expected_area;rtol=area_rtol,atol=area_atol_cm2) || error(
            "Patch area does not equal the matched Cartesian boundary-cell measure.",
        )
        key = (surface,cell)
        key in occupied && error("Multiple source patches map to the same boundary-cell face.")
        push!(occupied,key)
        surface_indices[patch] = surface
        cell_indices[patch] = cell
    end
    return surface_indices,cell_indices
end

function _projection_empty_surface_array(
    geometry::Geometry,
    energy_groups::Int64,
    coefficients::Int64,
)
    dimension = get_dimension(geometry)
    Nx = geometry.number_of_voxels["x"]
    Ny = dimension ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = dimension ≥ 3 ? geometry.number_of_voxels["z"] : 1
    result = Array{Union{Array{Float64},Float64}}(
        undef,energy_groups,coefficients,2*dimension,
    )
    for group in 1:energy_groups, coefficient in 1:coefficients
        if dimension == 1
            result[group,coefficient,1] = 0.0
            result[group,coefficient,2] = 0.0
        elseif dimension == 2
            result[group,coefficient,1] = zeros(Float64,Ny)
            result[group,coefficient,2] = zeros(Float64,Ny)
            result[group,coefficient,3] = zeros(Float64,Nx)
            result[group,coefficient,4] = zeros(Float64,Nx)
        else
            result[group,coefficient,1] = zeros(Float64,Ny,Nz)
            result[group,coefficient,2] = zeros(Float64,Ny,Nz)
            result[group,coefficient,3] = zeros(Float64,Nx,Nz)
            result[group,coefficient,4] = zeros(Float64,Nx,Nz)
            result[group,coefficient,5] = zeros(Float64,Nx,Ny)
            result[group,coefficient,6] = zeros(Float64,Nx,Ny)
        end
    end
    return result
end

function _projection_add_boundary_value!(
    destination,
    group::Int64,
    coefficient::Int64,
    surface::Int64,
    cell::Int64,
    value::Float64,
)
    current = destination[group,coefficient,surface]
    if current isa Float64
        cell == 1 || error("One-dimensional boundary cell index must be one.")
        destination[group,coefficient,surface] = current+value
    else
        current[cell] += value
    end
    return destination
end

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
    get_tag(source.particle) == get_tag(get_particle(solver)) || error(
        "Boundary-source particle does not match the selected solver.",
    )

    group_map = _projection_energy_group_map(
        source.energy_edges_eV,cross_sections,source.particle,
    )
    stored_particle = _projection_stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    directions,weights,direction_matrix,Qdims = _projection_solver_quadrature(solver,geometry)
    direction_map = _projection_direction_map(
        source.directions,source.quadrature_weights,direction_matrix,weights;
        direction_atol=direction_atol,
    )
    surface_indices,cell_indices = _projection_patch_mapping(source,geometry)

    basis_by_surface = Dict{Int64,Any}()
    coefficient_count = 0
    for surface in unique(surface_indices)
        basis = surface_angular_polynomial_basis(
            directions,weights,get_legendre_order(solver),
            get_angular_boltzmann(solver),Qdims,_projection_surface_name(surface),
        )
        basis_by_surface[surface] = basis
        Np = basis[1]
        coefficient_count == 0 && (coefficient_count = Np)
        Np == coefficient_count || error(
            "Cartesian boundaries produced incompatible source-basis sizes.",
        )
    end
    projected = _projection_empty_surface_array(geometry,Ng,coefficient_count)
    target_current = get_incoming_current(source)
    projected_current = zeros(Float64,size(target_current))

    for patch in eachindex(source.patch_ids)
        surface = surface_indices[patch]
        Np,Mn,Dn,incoming_directions,full_to_half,_,_ = basis_by_surface[surface]
        normal = view(source.normals,patch,:)
        for source_group in axes(source.angular_flux,2)
            radiant_group = group_map[source_group]
            discrete = zeros(Float64,length(incoming_directions))
            for source_direction in axes(source.angular_flux,3)
                solver_direction = direction_map[source_direction]
                half_direction = full_to_half[solver_direction]
                value = source.angular_flux[patch,source_group,source_direction]
                if value > 0.0
                    half_direction > 0 || error(
                        "A nonzero boundary source maps to an outgoing solver ordinate.",
                    )
                    discrete[half_direction] += value
                end
            end
            moments = Dn*discrete
            reconstruction = Mn*moments
            scale = max(1.0,maximum(abs.(reconstruction)))
            minimum(reconstruction) ≥ -Float64(positivity_atol)*scale || error(
                "Boundary projection produced negative angular-source lobes.",
            )

            target_density = target_current[patch,source_group]/source.areas_cm2[patch]
            calculated_density = 0.0
            for half_direction in eachindex(incoming_directions)
                full_direction = incoming_directions[half_direction]
                calculated_density += weights[full_direction]*
                    abs(dot(view(direction_matrix,full_direction,:),normal))*
                    reconstruction[half_direction]
            end
            if target_density > current_atol
                calculated_density > 0.0 || error(
                    "Boundary projection has nonpositive current for a nonzero source.",
                )
                moments .*= target_density/calculated_density
                reconstruction = Mn*moments
                calculated_density = 0.0
                for half_direction in eachindex(incoming_directions)
                    full_direction = incoming_directions[half_direction]
                    calculated_density += weights[full_direction]*
                        abs(dot(view(direction_matrix,full_direction,:),normal))*
                        reconstruction[half_direction]
                end
            elseif abs(calculated_density) > current_atol
                error("Boundary projection generated current from a zero source group.")
            end
            isapprox(calculated_density,target_density;rtol=current_rtol,atol=current_atol) || error(
                "Boundary current-density closure failed after projection.",
            )

            for coefficient in 1:Np
                _projection_add_boundary_value!(
                    projected,radiant_group,coefficient,surface,cell_indices[patch],
                    moments[coefficient],
                )
            end
            projected_current[patch,source_group] =
                calculated_density*source.areas_cm2[patch]
        end
    end

    receipt = Boundary_Projection_Receipt(
        group_map,direction_map,surface_indices,cell_indices,target_current,projected_current,
        _projection_maximum_relative_error(target_current,projected_current,current_atol),
        source.normalization.source_hash,
    )
    return projected,receipt
end

function _projection_merge_surface_source!(this::Source,additional)
    size(this.surface_sources,1) == size(additional,1) || error(
        "Projected and initialized surface sources have different energy dimensions.",
    )
    size(this.surface_sources,3) == size(additional,3) || error(
        "Projected and initialized surface sources have different boundary dimensions.",
    )
    Ng = size(additional,1)
    coefficients = max(size(this.surface_sources,2),size(additional,2))
    merged = _projection_empty_surface_array(this.geometry,Ng,coefficients)
    for group in 1:Ng, coefficient in axes(this.surface_sources,2), surface in axes(this.surface_sources,3)
        merged[group,coefficient,surface] += this.surface_sources[group,coefficient,surface]
    end
    for group in 1:Ng, coefficient in axes(additional,2), surface in axes(additional,3)
        merged[group,coefficient,surface] += additional[group,coefficient,surface]
    end
    this.surface_sources = merged
    return this
end

function _projection_linear_voxel_index(geometry::Geometry,voxel_id::Int64)
    dimension = get_dimension(geometry)
    Nx = geometry.number_of_voxels["x"]
    Ny = dimension ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = dimension ≥ 3 ? geometry.number_of_voxels["z"] : 1
    total = Nx*Ny*Nz
    1 ≤ voxel_id ≤ total || error("Volume-source voxel ID is outside the Radiant geometry.")
    cartesian = CartesianIndices((Nx,Ny,Nz))[voxel_id]
    return Tuple(cartesian)
end

function _projection_voxel_volume(geometry::Geometry,ix::Int64,iy::Int64,iz::Int64)
    dimension = get_dimension(geometry)
    dimension == 1 && return geometry.volume_per_voxel[ix]
    dimension == 2 && return geometry.volume_per_voxel[ix,iy]
    return geometry.volume_per_voxel[ix,iy,iz]
end

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
    geometry.is_build || error("Geometry must be built before volume-source projection.")
    get_tag(source.particle) == get_tag(get_particle(solver)) || error(
        "Volume-source particle does not match the selected solver.",
    )

    group_map = _projection_energy_group_map(
        source.energy_edges_eV,cross_sections,source.particle,
    )
    stored_particle = _projection_stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    directions,weights,direction_matrix,Qdims = _projection_solver_quadrature(solver,geometry)
    Np,Mn,Dn,_,_ = angular_polynomial_basis(
        directions,weights,get_legendre_order(solver),
        get_angular_boltzmann(solver),Qdims,
    )
    direction_map = source.angular_representation == :ordinates ?
        _projection_direction_map(
            source.directions,source.quadrature_weights,direction_matrix,weights;
            direction_atol=direction_atol,
        ) : Int64[]

    if source.angular_representation == :moments
        get(source.provenance,"angular_basis","") == "radiant-volume-moments/v1" || error(
            "Moment projection requires angular_basis=radiant-volume-moments/v1.",
        )
        size(source.values,3) == Np || error(
            "Moment source coefficient count differs from the selected Radiant basis.",
        )
    end

    _,_,Nm = get_schemes(solver,geometry,get_is_full_coupling(solver))
    dimension = get_dimension(geometry)
    Nx = geometry.number_of_voxels["x"]
    Ny = dimension ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = dimension ≥ 3 ? geometry.number_of_voxels["z"] : 1
    projected = zeros(Float64,Ng,Np,Nm[5],Nx,Ny,Nz)
    target_rate = zeros(Float64,Ng)
    projected_rate = zeros(Float64,Ng)

    for source_voxel in eachindex(source.voxel_ids)
        ix,iy,iz = _projection_linear_voxel_index(geometry,source.voxel_ids[source_voxel])
        expected_volume = _projection_voxel_volume(geometry,ix,iy,iz)
        isapprox(
            source.voxel_volumes_cm3[source_voxel],expected_volume;
            rtol=volume_rtol,atol=volume_atol_cm3,
        ) || error("Volume-source voxel volume does not match the Radiant geometry.")

        for source_group in axes(source.values,2)
            radiant_group = group_map[source_group]
            moments = zeros(Float64,Np)
            target_density = 0.0
            if source.angular_representation == :isotropic
                target_density = source.values[source_voxel,source_group,1]
                moments[1] = target_density
            elseif source.angular_representation == :ordinates
                discrete = zeros(Float64,size(direction_matrix,1))
                for source_direction in axes(source.values,3)
                    discrete[direction_map[source_direction]] +=
                        source.values[source_voxel,source_group,source_direction]
                end
                moments .= Dn*discrete
                target_density = sum(weights .* discrete)
            else
                moments .= view(source.values,source_voxel,source_group,:)
                get(source.provenance,"zeroth_moment_is_angle_integrated","false") == "true" || error(
                    "Moment projection requires zeroth_moment_is_angle_integrated=true.",
                )
                zeroth_index = try
                    parse(Int,get(source.provenance,"zeroth_moment_index",""))
                catch
                    error("Moment projection requires a valid zeroth_moment_index.")
                end
                1 ≤ zeroth_index ≤ Np || error("Zeroth moment index is outside the basis.")
                target_density = moments[zeroth_index]
            end

            reconstruction = Mn*moments
            calculated_density = sum(weights .* reconstruction)
            if target_density > rate_atol
                calculated_density > 0.0 || error(
                    "Volume projection has nonpositive scalar source for a nonzero source.",
                )
                moments .*= target_density/calculated_density
                reconstruction = Mn*moments
                calculated_density = sum(weights .* reconstruction)
            elseif abs(calculated_density) > rate_atol
                error("Volume projection generated scalar source from a zero source group.")
            end
            scale = max(1.0,maximum(abs.(reconstruction)))
            minimum(reconstruction) ≥ -Float64(positivity_atol)*scale || error(
                "Volume projection produced negative angular-source lobes.",
            )
            isapprox(calculated_density,target_density;rtol=rate_rtol,atol=rate_atol) || error(
                "Volume-source scalar-rate closure failed after projection.",
            )

            for coefficient in 1:Np
                projected[radiant_group,coefficient,1,ix,iy,iz] += moments[coefficient]
            end
            target_rate[radiant_group] += target_density*expected_volume
            projected_rate[radiant_group] += calculated_density*expected_volume
        end
    end

    receipt = Volume_Projection_Receipt(
        group_map,direction_map,target_rate,projected_rate,
        _projection_maximum_relative_error(target_rate,projected_rate,rate_atol),
        source.normalization.source_hash,
    )
    return projected,receipt
end

function add_source(this::Source,source::Boundary_Angular_Current_Source)
    projected,receipt = project_boundary_source(
        source,this.cross_sections,this.geometry,this.solver,
    )
    _projection_merge_surface_source!(this,projected)
    return receipt
end

function add_source(this::Source,source::Anisotropic_Volume_Source)
    projected,receipt = project_volume_source(
        source,this.cross_sections,this.geometry,this.solver,
    )
    size(projected) == size(this.volume_sources) || error(
        "Projected volume source is incompatible with the initialized Radiant source.",
    )
    this.volume_sources .+= projected
    return receipt
end
