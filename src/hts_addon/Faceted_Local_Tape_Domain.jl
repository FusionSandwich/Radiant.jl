"""One source-patch-to-CAD-facet association with geometric closure metadata."""
struct Facet_Boundary_Mapping_Entry
    source_patch_id::Int64
    facet_index::Int64
    surface_id::Int64
    centroid_distance_cm::Float64
    normal_dot::Float64
    facet_area_cm2::Float64
    source_area_cm2::Float64

    function Facet_Boundary_Mapping_Entry(;
        source_patch_id::Integer,
        facet_index::Integer,
        surface_id::Integer,
        centroid_distance_cm::Real,
        normal_dot::Real,
        facet_area_cm2::Real,
        source_area_cm2::Real,
    )
        source_patch_id >= 1 || error("Source patch identifiers must be positive.")
        facet_index >= 1 || error("Facet indices must be positive.")
        surface_id >= 1 || error("Surface identifiers must be positive.")
        distance = Float64(centroid_distance_cm)
        alignment = Float64(normal_dot)
        facet_area = Float64(facet_area_cm2)
        source_area = Float64(source_area_cm2)
        isfinite(distance) && distance >= 0.0 || error(
            "Source-to-facet centroid distance must be finite and nonnegative.",
        )
        isfinite(alignment) && -1.0 <= alignment <= 1.0 || error(
            "Source-to-facet normal dot product must lie in [-1,1].",
        )
        all(value -> isfinite(value) && value > 0.0,(facet_area,source_area)) || error(
            "Facet and source areas must be finite and positive.",
        )
        return new(
            Int64(source_patch_id),Int64(facet_index),Int64(surface_id),distance,
            alignment,facet_area,source_area,
        )
    end
end

"""Typed form of the dictionary returned by `map_boundary_patches_to_facets`."""
struct Facet_Boundary_Map
    schema::String
    source_hash::String
    geometry_hash::String
    mappings::Vector{Facet_Boundary_Mapping_Entry}

    function Facet_Boundary_Map(
        source_hash::AbstractString,
        geometry_hash::AbstractString,
        mappings::AbstractVector{Facet_Boundary_Mapping_Entry};
        schema::AbstractString="radiant.faceted_source_mapping/v1",
    )
        isempty(source_hash) && error("Facet boundary map source hash cannot be empty.")
        isempty(geometry_hash) && error("Facet boundary map geometry hash cannot be empty.")
        entries = Facet_Boundary_Mapping_Entry[mappings...]
        isempty(entries) && error("Facet boundary map cannot be empty.")
        source_ids = getfield.(entries,:source_patch_id)
        length(unique(source_ids)) == length(source_ids) || error(
            "Each source patch may occur only once in a facet boundary map.",
        )
        return new(String(schema),String(source_hash),String(geometry_hash),entries)
    end
end

function facet_boundary_map(mapping::AbstractDict)
    string(get(mapping,"schema","")) == "radiant.faceted_source_mapping/v1" || error(
        "Unsupported faceted source mapping schema.",
    )
    entries = Facet_Boundary_Mapping_Entry[]
    for raw in get(mapping,"mappings",Any[])
        push!(entries,Facet_Boundary_Mapping_Entry(
            source_patch_id=Int(raw["source_patch_id"]),
            facet_index=Int(raw["facet_index"]),
            surface_id=Int(raw["surface_id"]),
            centroid_distance_cm=Float64(raw["centroid_distance_cm"]),
            normal_dot=Float64(raw["normal_dot"]),
            facet_area_cm2=Float64(raw["facet_area_cm2"]),
            source_area_cm2=Float64(raw["source_area_cm2"]),
        ))
    end
    return Facet_Boundary_Map(
        string(mapping["source_hash"]),string(mapping["geometry_hash"]),entries;
        schema=string(mapping["schema"]),
    )
end

"""
Area-equivalent local Cartesian domain constructed from one triangular CAD/DAGMC tape facet.

The facet is flattened into a rectangular face with exactly the same area. The transverse width is
the physical tape width and the along-tape length is `facet_area/tape_width`. Physical layer
thicknesses and requested subdivisions are preserved exactly through the local `z` grid.
"""
struct Faceted_Local_Tape_Domain
    facet_patch::Faceted_Tape_Patch
    tape_stack::Tape_Stack_Definition
    geometry::Geometry
    material_indices::Dict{String,Int64}
    layer_voxel_ranges::Dict{String,UnitRange{Int64}}
    surrogate_length_cm::Float64
    area_ratio::Float64
    geometry_hash::String
    source_mesh_hash::String
    metadata::Dict{String,String}
end

struct Faceted_Local_Source_Receipt
    facet_index::Int64
    face::Symbol
    source_patch_ids::Vector{Int64}
    source_area_cm2::Float64
    facet_area_cm2::Float64
    area_relative_error::Float64
    source_current::Matrix{Float64}
    local_current::Matrix{Float64}
    maximum_current_relative_error::Float64
    source_hash::String
    local_source_hash::String
end

_facet_frame(patch::Faceted_Tape_Patch) = patch.local_to_global
_facet_normal(patch::Faceted_Tape_Patch) = vec(patch.local_to_global[:,3])

function _uniform_boundaries(lower::Real,upper::Real,count::Integer)
    count >= 1 || error("A local geometry axis requires at least one cell.")
    lower_value = Float64(lower)
    upper_value = Float64(upper)
    isfinite(lower_value) && isfinite(upper_value) && upper_value > lower_value || error(
        "Local geometry boundaries must be finite and increasing.",
    )
    return collect(range(lower_value,stop=upper_value,length=Int(count)+1))
end

function _cell_widths_and_centres(boundaries::Vector{Float64})
    widths = diff(boundaries)
    all(widths .> 0.0) || error("Local geometry contains a non-positive cell width.")
    centres = 0.5 .* (boundaries[1:end-1] .+ boundaries[2:end])
    return widths,centres
end

function _faceted_domain_hash(
    patch::Faceted_Tape_Patch,
    stack::Tape_Stack_Definition,
    material_indices::Dict{String,Int64},
    x_boundaries::Vector{Float64},
    y_boundaries::Vector{Float64},
    z_boundaries::Vector{Float64},
    source_mesh_hash::String,
)
    io = IOBuffer()
    print(io,"radiant.hts.faceted_local_tape_domain/v2\n")
    print(io,patch.facet_index,"\n",patch.surface_id,"\n",source_mesh_hash,"\n")
    print(io,stack.name,"\n",stack.width_cm,"\n")
    for layer in stack.layers
        print(io,layer.name,"|",layer.material_tag,"|",layer.thickness_cm,"|",
              layer.transport_subdivisions,"\n")
    end
    for key in sort(collect(keys(material_indices)))
        print(io,key,"=",material_indices[key],"\n")
    end
    print(io,join(x_boundaries,','),"\n",join(y_boundaries,','),"\n",
          join(z_boundaries,','),"\n")
    return bytes2hex(SHA.sha256(take!(io)))
end

"""
    build_faceted_local_tape_domain(patch, stack, material_indices; ...)

Create a fully built three-dimensional Radiant `Geometry` in the facet-local frame. Local `x` is
the facet tangent, `y` is the tape-width direction, and `z` is the facet normal. The equivalent
rectangular face has the exact triangular area, so integrated boundary current and layer volume are
conserved under flattening.
"""
function build_faceted_local_tape_domain(
    patch::Faceted_Tape_Patch,
    stack::Tape_Stack_Definition,
    material_indices::AbstractDict;
    along_cells::Integer=1,
    width_cells::Integer=1,
    layer_refinement::Integer=1,
    boundary_condition::AbstractString="void",
    source_mesh_hash::AbstractString="unbound",
    metadata::AbstractDict=Dict{String,String}(),
)
    along_cells >= 1 || error("along_cells must be positive.")
    width_cells >= 1 || error("width_cells must be positive.")
    layer_refinement >= 1 || error("layer_refinement must be positive.")
    boundary_condition in ("void","reflective","periodic") || error(
        "Unsupported local boundary condition $(boundary_condition).",
    )
    isempty(source_mesh_hash) && error("Source mesh hash cannot be empty.")

    frame = _facet_frame(patch)
    size(frame) == (3,3) || error("Facet local frame must be 3x3.")
    isapprox(frame'*frame,Matrix{Float64}(I,3,3);rtol=1.0e-9,atol=1.0e-9) || error(
        "Facet local frame must be orthonormal.",
    )
    det(frame) > 0.0 || error("Facet local frame must be right-handed.")

    indices = Dict{String,Int64}()
    for (key,value) in material_indices
        index = Int64(value)
        index >= 1 || error("Material indices must be positive and one-based.")
        indices[string(key)] = index
    end
    for layer in stack.layers
        haskey(indices,layer.material_tag) || error(
            "No material index was provided for layer material $(layer.material_tag).",
        )
    end

    surrogate_length = patch.area_cm2/stack.width_cm
    isfinite(surrogate_length) && surrogate_length > 0.0 || error(
        "Facet area and tape width do not define a positive surrogate length.",
    )
    thickness = get_total_thickness(stack)
    x_boundaries = _uniform_boundaries(
        -0.5*surrogate_length,0.5*surrogate_length,along_cells,
    )
    y_boundaries = _uniform_boundaries(
        -0.5*stack.width_cm,0.5*stack.width_cm,width_cells,
    )

    z_boundaries = Float64[-0.5*thickness]
    layer_ranges = Dict{String,UnitRange{Int64}}()
    z_material_indices = Int64[]
    z_cell_index = 0
    cursor = -0.5*thickness
    for layer in stack.layers
        subdivisions = layer.transport_subdivisions*Int(layer_refinement)
        local_boundaries = _uniform_boundaries(
            cursor,cursor+layer.thickness_cm,subdivisions,
        )
        append!(z_boundaries,local_boundaries[2:end])
        first_cell = z_cell_index+1
        append!(z_material_indices,fill(indices[layer.material_tag],subdivisions))
        z_cell_index += subdivisions
        layer_ranges[layer.name] = first_cell:z_cell_index
        cursor += layer.thickness_cm
    end
    isapprox(z_boundaries[end],0.5*thickness;rtol=1.0e-12,atol=1.0e-14) || error(
        "Layer boundaries do not reconstruct the tape thickness.",
    )

    x_widths,x_centres = _cell_widths_and_centres(x_boundaries)
    y_widths,y_centres = _cell_widths_and_centres(y_boundaries)
    z_widths,z_centres = _cell_widths_and_centres(z_boundaries)
    nx,ny,nz = length(x_widths),length(y_widths),length(z_widths)
    material_map = zeros(Int64,nx,ny,nz)
    volumes = zeros(Float64,nx,ny,nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        material_map[ix,iy,iz] = z_material_indices[iz]
        volumes[ix,iy,iz] = x_widths[ix]*y_widths[iy]*z_widths[iz]
    end
    expected_volume = patch.area_cm2*thickness
    isapprox(sum(volumes),expected_volume;rtol=1.0e-12,atol=1.0e-16) || error(
        "Flattened local-domain volume does not conserve facet area times tape thickness.",
    )

    geometry = Geometry()
    geometry.name = "faceted-tape-facet-$(patch.facet_index)"
    geometry.type = "cartesian"
    geometry.dimension = 3
    geometry.axis = ["x","y","z"]
    for (axis,boundaries,widths,centres) in (
        ("x",x_boundaries,x_widths,x_centres),
        ("y",y_boundaries,y_widths,y_centres),
        ("z",z_boundaries,z_widths,z_centres),
    )
        geometry.number_of_voxels[axis] = length(widths)
        geometry.voxels_boundaries[axis] = boundaries
        geometry.voxels_width[axis] = widths
        geometry.voxels_position[axis] = centres
        geometry.number_of_regions[axis] = 1
        geometry.voxels_per_region[axis] = [length(widths)]
        geometry.region_boundaries[axis] = [boundaries[1],boundaries[end]]
    end
    for boundary in ("X-","X+","Y-","Y+","Z-","Z+")
        geometry.boundary_conditions[boundary] = String(boundary_condition)
    end
    geometry.material_per_voxel = material_map
    geometry.volume_per_voxel = volumes
    geometry.is_build = true

    geometry_hash = _faceted_domain_hash(
        patch,stack,indices,x_boundaries,y_boundaries,z_boundaries,
        String(source_mesh_hash),
    )
    metadata_string = Dict{String,String}(
        "schema" => "radiant.hts.faceted_local_tape_domain/v2",
        "coordinate_x" => "facet-tangent",
        "coordinate_y" => "facet-width",
        "coordinate_z" => "through-thickness",
        "flattening" => "area-equivalent-rectangular-facet",
        "source_mesh_hash" => String(source_mesh_hash),
        "classification" => "software-qualified-physical-curvature-pending",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    area_ratio = surrogate_length*stack.width_cm/patch.area_cm2
    return Faceted_Local_Tape_Domain(
        patch,stack,geometry,indices,layer_ranges,surrogate_length,area_ratio,
        geometry_hash,String(source_mesh_hash),metadata_string,
    )
end

local_domain_frame(domain::Faceted_Local_Tape_Domain) = _facet_frame(domain.facet_patch)

function global_to_local(
    domain::Faceted_Local_Tape_Domain,
    position_cm::AbstractVector{<:Real},
)
    length(position_cm) == 3 || error("Global position must be a three-vector.")
    return local_domain_frame(domain)' *
           (Float64.(position_cm)-collect(domain.facet_patch.centroid_cm))
end

function local_to_global(
    domain::Faceted_Local_Tape_Domain,
    position_cm::AbstractVector{<:Real},
)
    length(position_cm) == 3 || error("Local position must be a three-vector.")
    return collect(domain.facet_patch.centroid_cm)+
           local_domain_frame(domain)*Float64.(position_cm)
end

function local_layer_for_voxel(domain::Faceted_Local_Tape_Domain,iz::Integer)
    1 <= iz <= domain.geometry.number_of_voxels["z"] || error(
        "Through-thickness voxel index is out of range.",
    )
    for layer in domain.tape_stack.layers
        iz in domain.layer_voxel_ranges[layer.name] && return layer
    end
    error("Through-thickness voxel was not assigned to a tape layer.")
end

function _relative_array_error(first::AbstractArray,second::AbstractArray;atol::Real=1.0e-30)
    size(first) == size(second) || error("Array comparison shapes differ.")
    denominator = max.(abs.(Float64.(first)),abs.(Float64.(second)),Float64(atol))
    return maximum(abs.(Float64.(first).-Float64.(second))./denominator)
end

function localize_boundary_source_to_facet(
    source::Boundary_Angular_Current_Source,
    mapping::AbstractDict,
    domain::Faceted_Local_Tape_Domain;
    kwargs...,
)
    return localize_boundary_source_to_facet(source,facet_boundary_map(mapping),domain;kwargs...)
end

"""
    localize_boundary_source_to_facet(source, mapping, domain; ...)

Aggregate every global source patch mapped to one CAD facet into area-equivalent local `Z-` and/or
`Z+` boundary patches. Directional integrated current is conserved before inversion to local
angular flux. Opposite-facing source populations remain separate so they cannot be mixed into one
boundary normal.
"""
function localize_boundary_source_to_facet(
    source::Boundary_Angular_Current_Source,
    mapping::Facet_Boundary_Map,
    domain::Faceted_Local_Tape_Domain;
    area_rtol::Real=1.0e-6,
    normal_alignment_minimum::Real=0.999,
)
    mapping.source_hash == source.normalization.source_hash || error(
        "Facet mapping source hash does not match the boundary source.",
    )
    mapping.geometry_hash == domain.source_mesh_hash || error(
        "Facet mapping geometry hash does not match the local domain source mesh.",
    )
    facet_index = domain.facet_patch.facet_index
    entries = [entry for entry in mapping.mappings if entry.facet_index == facet_index]
    isempty(entries) && error("No boundary patches map to the requested facet.")
    source_index = Dict(id => index for (index,id) in enumerate(source.patch_ids))
    selected = Int64[]
    for entry in entries
        haskey(source_index,entry.source_patch_id) || error(
            "Facet mapping references an unknown source patch.",
        )
        push!(selected,source_index[entry.source_patch_id])
    end
    length(unique(selected)) == length(selected) || error(
        "The same source patch was mapped to the facet more than once.",
    )

    frame = local_domain_frame(domain)
    local_directions = source.directions*frame
    facet_normal = _facet_normal(domain.facet_patch)
    face_indices = Dict(:front => Int64[],:back => Int64[])
    for index in selected
        alignment = dot(view(source.normals,index,:),facet_normal)
        abs(alignment) >= normal_alignment_minimum || error(
            "Source patch normal is not aligned with the mapped CAD facet.",
        )
        push!(face_indices[alignment > 0.0 ? :back : :front],index)
    end

    outputs = Dict{Symbol,Boundary_Angular_Current_Source}()
    receipts = Dict{Symbol,Faceted_Local_Source_Receipt}()
    for face in (:front,:back)
        indices_face = face_indices[face]
        isempty(indices_face) && continue
        sign_z = face == :back ? 1.0 : -1.0
        local_normal = [0.0 0.0 sign_z]
        local_tangent_1 = [1.0 0.0 0.0]
        local_tangent_2 = face == :back ? [0.0 1.0 0.0] : [0.0 -1.0 0.0]
        source_area = sum(source.areas_cm2[indices_face])
        area_error = abs(source_area-domain.facet_patch.area_cm2)/
                     max(domain.facet_patch.area_cm2,eps(Float64))
        area_error <= area_rtol || error(
            "Mapped source-patch area does not close the CAD facet area.",
        )

        ngroup = length(source.energy_edges_eV)-1
        ndirection = length(source.quadrature_weights)
        integrated_current_by_direction = zeros(Float64,ngroup,ndirection)
        integrated_variance_by_direction = isnothing(source.variance) ? nothing :
            zeros(Float64,ngroup,ndirection)
        source_current = zeros(Float64,1,ngroup)
        source_current_all = get_incoming_current(source)
        for index in indices_face
            area = source.areas_cm2[index]
            normal = view(source.normals,index,:)
            for direction in 1:ndirection
                mu = dot(view(source.directions,direction,:),normal)
                if mu < 0.0
                    factor = area*source.quadrature_weights[direction]*abs(mu)
                    integrated_current_by_direction[:,direction] .+=
                        factor .* view(source.angular_flux,index,:,direction)
                    if !isnothing(integrated_variance_by_direction)
                        integrated_variance_by_direction[:,direction] .+=
                            factor^2 .* view(source.variance,index,:,direction)
                    end
                end
            end
            source_current .+= reshape(source_current_all[index,:],1,ngroup)
        end
        directional_current_density = reshape(
            integrated_current_by_direction./domain.facet_patch.area_cm2,
            1,ngroup,ndirection,
        )
        directional_variance_density = isnothing(integrated_variance_by_direction) ? nothing :
            reshape(
                integrated_variance_by_direction./domain.facet_patch.area_cm2^2,
                1,ngroup,ndirection,
            )
        local_centroid = [0.0 0.0 sign_z*0.5*get_total_thickness(domain.tape_stack)]
        local_hash = bytes2hex(SHA.sha256(codeunits(
            string(source.normalization.source_hash,"|facet=",facet_index,"|face=",face),
        )))
        normalization = Source_Normalization(
            basis=source.normalization.basis,
            source_rate_per_s=source.normalization.source_rate_per_s,
            symmetry_factor=source.normalization.symmetry_factor,
            time_interval_s=source.normalization.time_interval_s,
            time_class=source.normalization.time_class,
            source_hash=local_hash,
            provenance=merge(copy(source.normalization.provenance),Dict(
                "parent_source_hash" => source.normalization.source_hash,
                "facet_geometry_hash" => domain.geometry_hash,
            )),
        )
        local_source = boundary_source_from_directional_current(
            source.particle,[facet_index],local_centroid,[domain.facet_patch.area_cm2],
            local_normal,local_tangent_1,local_tangent_2,source.energy_edges_eV,
            local_directions,source.quadrature_weights,directional_current_density,
            normalization;
            variance=directional_variance_density,
            provenance=merge(copy(source.provenance),Dict(
                "coordinate_frame" => "facet-local",
                "facet_index" => string(facet_index),
                "local_face" => string(face),
                "local_geometry_hash" => domain.geometry_hash,
            )),
        )
        local_current = get_incoming_current(local_source)
        maximum_error = _relative_array_error(source_current,local_current)
        maximum_error <= 1.0e-10 || error(
            "Global-to-local faceted source mapping failed current closure.",
        )
        outputs[face] = local_source
        receipts[face] = Faceted_Local_Source_Receipt(
            facet_index,face,source.patch_ids[indices_face],source_area,
            domain.facet_patch.area_cm2,area_error,source_current,local_current,
            maximum_error,source.normalization.source_hash,local_hash,
        )
    end
    return outputs,receipts
end

function faceted_local_domain_receipt(domain::Faceted_Local_Tape_Domain)
    layer_receipt = Dict{String,Any}()
    for layer in domain.tape_stack.layers
        layer_receipt[layer.name] = Dict(
            "material_tag" => layer.material_tag,
            "material_index" => domain.material_indices[layer.material_tag],
            "thickness_cm" => layer.thickness_cm,
            "voxel_range" => [first(domain.layer_voxel_ranges[layer.name]),
                               last(domain.layer_voxel_ranges[layer.name])],
        )
    end
    return Dict{String,Any}(
        "schema" => "radiant.hts.faceted_local_tape_domain_receipt/v2",
        "facet_index" => domain.facet_patch.facet_index,
        "surface_id" => domain.facet_patch.surface_id,
        "geometry_hash" => domain.geometry_hash,
        "source_mesh_hash" => domain.source_mesh_hash,
        "facet_area_cm2" => domain.facet_patch.area_cm2,
        "surrogate_length_cm" => domain.surrogate_length_cm,
        "tape_width_cm" => domain.tape_stack.width_cm,
        "tape_thickness_cm" => get_total_thickness(domain.tape_stack),
        "area_ratio" => domain.area_ratio,
        "voxel_shape" => collect(size(domain.geometry.material_per_voxel)),
        "layers" => layer_receipt,
        "physical_curvature_qualification" => false,
    )
end
