"""
    Boundary_Projection_Receipt

Audit receipt for projection of a patch/group/direction boundary source into Radiant's Cartesian
half-range angular basis. Currents are extensive over each patch; the values installed in the
transport sweep remain angular-flux densities.
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
    isnothing(index) && error(
        "Cross sections do not contain particle tag $(get_tag(particle)).",
    )
    return get_particles(cross_sections)[index]
end

function _energy_group_map(
    source_edges_eV::AbstractVector{<:Real},
    cross_sections::Cross_Sections,
    particle::Particle;
    rtol::Real=1.0e-10,
    atol_eV::Real=1.0e-8,
)
    stored_particle = _stored_particle(cross_sections,particle)
    ismissing(cross_sections.energy_boundaries) && error(
        "Cross-section energy boundaries are required for source projection.",
    )
    radiant_edges_eV = 1.0e6 .* Float64.(
        get_energy_boundaries(cross_sections,stored_particle),
    )
    source_edges = Float64.(source_edges_eV)
    mapping = zeros(Int64,length(source_edges)-1)
    used = Set{Int64}()

    for source_group in eachindex(mapping)
        source_interval = extrema(source_edges[source_group:source_group+1])
        matches = Int64[]
        for radiant_group in 1:length(radiant_edges_eV)-1
            radiant_interval = extrema(
                radiant_edges_eV[radiant_group:radiant_group+1],
            )
            if isapprox(source_interval[1],radiant_interval[1];rtol=rtol,atol=atol_eV) &&
               isapprox(source_interval[2],radiant_interval[2];rtol=rtol,atol=atol_eV)
                push!(matches,radiant_group)
            end
        end
        length(matches) == 1 || error(
            "Source energy group $(source_group) does not match exactly one Radiant group.",
        )
        matches[1] in used && error(
            "Multiple source groups map to the same Radiant group.",
        )
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
    return Ω,Float64.(w),hcat(Ω[1],Ω[2],Ω[3]),Qdims
end

function _direction_map(
    source_directions::AbstractMatrix{<:Real},
    source_weights::AbstractVector{<:Real},
    solver_directions::AbstractMatrix{<:Real},
    solver_weights::AbstractVector{<:Real};
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
                view(solver_directions,solver_index,:),
            )
            for solver_index in axes(solver_directions,1)
        ]
        solver_index = argmin(distances)
        distances[solver_index] ≤ direction_atol || error(
            "Source direction $(source_index) is not an ordinate of the Radiant quadrature.",
        )
        solver_index in used && error(
            "Multiple source directions map to one Radiant ordinate.",
        )
        isapprox(
            Float64(source_weights[source_index]),
            Float64(solver_weights[solver_index]);
            rtol=weight_rtol,
            atol=weight_atol,
        ) || error(
            "Source and Radiant quadrature weights differ for direction $(source_index).",
        )
        mapping[source_index] = solver_index
        push!(used,solver_index)
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
            error("Cartesian sources support one, two, or three dimensions.")
        end
    end
    return output
end

function _merge_surface_source!(this::Source,new_source)
    old_source = this.surface_sources
    size(old_source,1) == size(new_source,1) || error(
        "Surface-source group dimensions are incompatible.",
    )
    size(old_source,3) == size(new_source,3) || error(
        "Surface-source boundary dimensions are incompatible.",
    )
    Ng = size(old_source,1)
    Np = max(size(old_source,2),size(new_source,2))
    merged = _empty_surface_source_array(Ng,Np,this.geometry)
    for igroup in 1:Ng, coefficient in axes(old_source,2), boundary in axes(old_source,3)
        merged[igroup,coefficient,boundary] += old_source[igroup,coefficient,boundary]
    end
    for igroup in 1:Ng, coefficient in axes(new_source,2), boundary in axes(new_source,3)
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
    distances = abs.(geometry.voxels_position[axis] .- Float64(coordinate))
    index = argmin(distances)
    distances[index] ≤ tolerance_cm || error(
        "Patch centroid does not coincide with a boundary-voxel centroid along $(axis).",
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
    get_type(geometry) == "cartesian" || error(
        "Boundary projection currently supports Cartesian Radiant geometry only.",
    )
    canonical_normals = (
        (-1.0,0.0,0.0),(1.0,0.0,0.0),
        (0.0,-1.0,0.0),(0.0,1.0,0.0),
        (0.0,0.0,-1.0),(0.0,0.0,1.0),
    )
    distances = [norm(Float64.(normal) .- collect(candidate)) for candidate in canonical_normals]
    boundary = argmin(distances)
    Ndims = get_dimension(geometry)
    boundary ≤ 2*Ndims && distances[boundary] ≤ normal_atol || error(
        "Patch normal is not an outward normal of the Cartesian geometry.",
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
        "Patch measure does not match the corresponding boundary-voxel face.",
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

Project a `Boundary_Angular_Current_Source` into Radiant's SN half-range boundary basis.
The source directions and weights must match a subset of the selected Radiant quadrature, each
patch must equal one Cartesian boundary-voxel face, and energy intervals must match exactly.

The installed moments describe angular-flux density. Patch area is used only in the closure
receipt:

```math
J^-_{k,g}=A_k\sum_{m\in-}w_m|\Omega_m\cdot n_k|\psi^-_{k,g,m}.
```

A current-preserving scalar rescale may remove basis truncation error. Negative reconstructed
angular flux is rejected and never clipped.
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
    get_tag(source.particle) == get_tag(get_particle(solver)) || error(
        "Boundary-source particle does not match the solver.",
    )

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
    Np,Mn,Dn,incoming_indices,_,_,_ = surface_angular_polynomial_basis(
        Ω,w,get_legendre_order(solver),get_angular_boltzmann(solver),Qdims,
        get_dimension(geometry),get_type(geometry),
    )

    projected_source = _empty_surface_source_array(Ng,Np,geometry)
    target_current = zeros(Float64,length(source.patch_ids),Ng)
    projected_current = similar(target_current)
    surface_index = zeros(Int64,length(source.patch_ids))
    surface_cell = Vector{NTuple{3,Int64}}(undef,length(source.patch_ids))
    occupied_cells = Set{Tuple{Int64,NTuple{3,Int64}}}()
    Ndims = get_dimension(geometry)

    for patch in eachindex(source.patch_ids)
        boundary,cell,patch_measure = _cartesian_patch_location(
            geometry,view(source.centroids_cm,patch,:),view(source.normals,patch,:),
            source.areas_cm2[patch],
        )
        key = (boundary,cell)
        key in occupied_cells && error(
            "Multiple source patches map to one boundary-voxel face.",
        )
        push!(occupied_cells,key)
        surface_index[patch] = boundary
        surface_cell[patch] = cell
        incoming = incoming_indices[boundary]

        for source_group in axes(source.angular_flux,2)
            radiant_group = group_map[source_group]
            discrete_flux = zeros(Float64,size(solver_directions,1))
            for source_direction in axes(source.angular_flux,3)
                discrete_flux[direction_map[source_direction]] =
                    source.angular_flux[patch,source_group,source_direction]
            end

            moments = Dn[boundary] * discrete_flux[incoming]
            reconstruction = Mn[boundary] * moments
            target_density = 0.0
            projected_density = 0.0
            for (half_index,solver_direction) in enumerate(incoming)
                cosine = abs(dot(
                    view(solver_directions,solver_direction,:),
                    view(source.normals,patch,:),
                ))
                target_density += w[solver_direction] * cosine * discrete_flux[solver_direction]
                projected_density += w[solver_direction] * cosine * reconstruction[half_index]
            end

            if target_density > current_atol
                projected_density > 0.0 || error(
                    "Boundary projection has nonpositive current for a nonzero source.",
                )
                moments .*= target_density/projected_density
                reconstruction = Mn[boundary] * moments
            elseif abs(projected_density) > current_atol
                error("Boundary projection generated current from a zero-current source group.")
            end

            scale = max(1.0,maximum(abs.(reconstruction)))
            minimum(reconstruction) ≥ -Float64(positivity_atol)*scale || error(
                "Boundary projection produced negative incoming angular-flux lobes.",
            )
            projected_density = 0.0
            for (half_index,solver_direction) in enumerate(incoming)
                cosine = abs(dot(
                    view(solver_directions,solver_direction,:),
                    view(source.normals,patch,:),
                ))
                projected_density += w[solver_direction] * cosine * reconstruction[half_index]
            end
            isapprox(projected_density,target_density;rtol=current_rtol,atol=current_atol) || error(
                "Boundary-current density failed closure after projection.",
            )

            target_current[patch,radiant_group] += patch_measure * target_density
            projected_current[patch,radiant_group] += patch_measure * projected_density
            for coefficient in 1:Np
                _add_surface_moment!(
                    projected_source,radiant_group,coefficient,boundary,cell,
                    moments[coefficient],Ndims,
                )
            end
        end
    end

    return projected_source,Boundary_Projection_Receipt(
        group_map,direction_map,surface_index,surface_cell,target_current,projected_current,
        _maximum_relative_error(target_current,projected_current,current_atol),
        source.normalization.source_hash,
    )
end

function _linear_voxel_index(geometry::Geometry,voxel_id::Int64)
    Nx = geometry.number_of_voxels["x"]
    Ny = get_dimension(geometry) ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(geometry) ≥ 3 ? geometry.number_of_voxels["z"] : 1
    1 ≤ voxel_id ≤ Nx*Ny*Nz || error(
        "Volume-source voxel identifier is outside the Radiant geometry.",
    )
    return Tuple(CartesianIndices((Nx,Ny,Nz))[voxel_id])
end

function _geometry_voxel_volume(geometry::Geometry,ix::Int64,iy::Int64,iz::Int64)
    if get_dimension(geometry) == 1
        return geometry.volume_per_voxel[ix]
    elseif get_dimension(geometry) == 2
        return geometry.volume_per_voxel[ix,iy]
    end
    return geometry.volume_per_voxel[ix,iy,iz]
end

"""
    project_volume_source(source, cross_sections, geometry, solver)

Project isotropic, ordinate, or native Radiant moment sources into the SN volume-source tensor.
Voxel identifiers use Julia's one-based linear indexing into `(Nx,Ny,Nz)`. The scalar source rate
and nonnegative reconstructed angular source are protected. Signed higher moments are allowed.
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
    geometry.is_build || error("Geometry must be built before volume-source projection.")
    get_tag(source.particle) == get_tag(get_particle(solver)) || error(
        "Volume-source particle does not match the solver.",
    )

    group_map = _energy_group_map(source.energy_edges_eV,cross_sections,source.particle)
    stored_particle = _stored_particle(cross_sections,source.particle)
    Ng = get_number_of_groups(cross_sections,stored_particle)
    Ω,w,solver_directions,Qdims = _solver_quadrature(solver,geometry)
    Np,Mn,Dn,_,_ = angular_polynomial_basis(
        Ω,w,get_legendre_order(solver),get_angular_boltzmann(solver),Qdims,
    )
    direction_map = source.angular_representation == :ordinates ? _direction_map(
        source.directions,source.quadrature_weights,solver_directions,w;
        direction_atol=direction_atol,
    ) : Int64[]

    if source.angular_representation == :moments
        get(source.provenance,"angular_basis","") == "radiant-volume-moments/v1" || error(
            "Moment source must declare angular_basis=radiant-volume-moments/v1.",
        )
        size(source.values,3) == Np || error(
            "Moment-source coefficient count does not match the Radiant basis.",
        )
        get(source.provenance,"zeroth_moment_index","") == "1" || error(
            "Radiant moment coefficient 1 must be declared as the zeroth moment.",
        )
        get(source.provenance,"zeroth_moment_is_angle_integrated","false") == "true" || error(
            "Radiant moment coefficient 1 must be declared angle-integrated.",
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
        voxel_volume = _geometry_voxel_volume(geometry,ix,iy,iz)
        isapprox(
            source.voxel_volumes_cm3[source_voxel],voxel_volume;
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
                discrete_source = zeros(Float64,size(solver_directions,1))
                for source_direction in axes(source.values,3)
                    discrete_source[direction_map[source_direction]] =
                        source.values[source_voxel,source_group,source_direction]
                end
                moments .= Dn * discrete_source
                target_density = sum(w .* discrete_source)
            else
                moments .= view(source.values,source_voxel,source_group,:)
                target_density = moments[1]
            end

            target_density ≥ -rate_atol || error(
                "Angle-integrated volume-source density must be nonnegative.",
            )
            abs(target_density) ≤ rate_atol && (target_density = 0.0)
            reconstruction = Mn * moments
            projected_density = sum(w .* reconstruction)
            if target_density > rate_atol
                projected_density > 0.0 || error(
                    "Volume projection has nonpositive scalar source for a nonzero source.",
                )
                moments .*= target_density/projected_density
                reconstruction = Mn * moments
            elseif abs(projected_density) > rate_atol
                error("Volume projection generated scalar source from a zero source group.")
            end

            scale = max(1.0,maximum(abs.(reconstruction)))
            minimum(reconstruction) ≥ -Float64(positivity_atol)*scale || error(
                "Volume projection produced negative angular-source lobes.",
            )
            projected_density = sum(w .* reconstruction)
            isapprox(projected_density,target_density;rtol=rate_rtol,atol=rate_atol) || error(
                "Volume-source scalar rate failed closure after projection.",
            )

            for coefficient in 1:Np
                projected_source[radiant_group,coefficient,1,ix,iy,iz] += moments[coefficient]
            end
            target_rate[radiant_group] += target_density * voxel_volume
            projected_rate[radiant_group] += projected_density * voxel_volume
        end
    end

    return projected_source,Volume_Projection_Receipt(
        group_map,direction_map,target_rate,projected_rate,
        _maximum_relative_error(target_rate,projected_rate,rate_atol),
        source.normalization.source_hash,
    )
end

"""
    add_source(this::Source, source::Boundary_Angular_Current_Source)

Project and add a conservative tabulated boundary source. Returns a closure receipt.
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

Project and add a voxel-resolved angular volume source. Returns a closure receipt.
"""
function add_source(this::Source,source::Anisotropic_Volume_Source)
    projected,receipt = project_volume_source(
        source,this.cross_sections,this.geometry,this.solver,
    )
    size(projected) == size(this.volume_sources) || error(
        "Projected volume-source tensor is incompatible with the initialized source.",
    )
    this.volume_sources .+= projected
    return receipt
end
