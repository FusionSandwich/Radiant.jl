"""
    Tape_Geometry_Map

Mapping between a planar `Tape_Stack_Definition` and a built Cartesian Radiant geometry. The x-axis
is always the front-to-back layer-normal direction; y is tape width and z is local arc length when
those dimensions are enabled.
"""
struct Tape_Geometry_Map
    stack::Tape_Stack_Definition
    geometry::Geometry
    layer_voxel_ranges_x::Vector{UnitRange{Int64}}
    voxel_layer_index_x::Vector{Int64}
    activation_voxel_ids::Vector{Int64}
    atomistic_voxel_ids::Vector{Int64}
    metadata::Dict{String,String}
end

function _material_for_tag(cross_sections::Cross_Sections,tag::String)
    index = findfirst(material -> get_tag(material) == tag,get_materials(cross_sections))
    isnothing(index) && error("No Radiant material matches tape layer tag $(tag).")
    return get_materials(cross_sections)[index]
end

function _layer_material_array(
    stack::Tape_Stack_Definition,
    cross_sections::Cross_Sections,
    dimension::Int64,
)
    layer_materials = Material[
        _material_for_tag(cross_sections,layer.material_tag) for layer in stack.layers
    ]
    if dimension == 1
        return layer_materials
    elseif dimension == 2
        result = Array{Material}(undef,length(layer_materials),1)
        for layer in eachindex(layer_materials)
            result[layer,1] = layer_materials[layer]
        end
        return result
    elseif dimension == 3
        result = Array{Material}(undef,length(layer_materials),1,1)
        for layer in eachindex(layer_materials)
            result[layer,1,1] = layer_materials[layer]
        end
        return result
    end
    error("Tape geometry dimension must be one, two, or three.")
end

function _layer_voxel_ranges(stack::Tape_Stack_Definition,refinement::Int64)
    ranges = Vector{UnitRange{Int64}}(undef,length(stack.layers))
    voxel_layers = Int64[]
    first_voxel = 1
    for (layer_index,layer) in enumerate(stack.layers)
        count = layer.transport_subdivisions*refinement
        last_voxel = first_voxel+count-1
        ranges[layer_index] = first_voxel:last_voxel
        append!(voxel_layers,fill(layer_index,count))
        first_voxel = last_voxel+1
    end
    return ranges,voxel_layers
end

function _linear_voxel_ids_for_layers(
    voxel_layers_x::Vector{Int64},
    Ny::Int64,
    Nz::Int64,
    selected_layers::Vector{Int64},
)
    selected = Set(selected_layers)
    ids = Int64[]
    Nx = length(voxel_layers_x)
    for cartesian in CartesianIndices((Nx,Ny,Nz))
        voxel_layers_x[cartesian[1]] in selected || continue
        push!(ids,LinearIndices((Nx,Ny,Nz))[cartesian])
    end
    return ids
end

"""
    build_planar_tape_geometry(stack, cross_sections; ...)

Build a one-, two-, or three-dimensional Cartesian tape segment. Physical layer boundaries are
exact decimal inputs from the stack. Every layer receives
`layer.transport_subdivisions*thickness_refinement` x-cells. Two-dimensional cases add the finite
tape width, while three-dimensional cases additionally require a positive local segment length.

This is a straight local geometry. Curvature is handled by a separate piecewise-flat atlas and must
be convergence-tested before a curved-device response is accepted.
"""
function build_planar_tape_geometry(
    stack::Tape_Stack_Definition,
    cross_sections::Cross_Sections;
    dimension::Integer=1,
    length_cm::Union{Nothing,Real}=nothing,
    thickness_refinement::Integer=1,
    width_subdivisions::Integer=1,
    length_subdivisions::Integer=1,
    boundary_condition::AbstractString="void",
    origin_cm::Real=0.0,
)
    Ndims = Int64(dimension)
    Ndims in (1,2,3) || error("Tape geometry dimension must be one, two, or three.")
    refinement = Int64(thickness_refinement)
    refinement ≥ 1 || error("Thickness refinement must be at least one.")
    width_subdivisions ≥ 1 || error("Width subdivisions must be at least one.")
    length_subdivisions ≥ 1 || error("Length subdivisions must be at least one.")
    boundary = lowercase(String(boundary_condition))
    boundary in ("void","reflective","periodic") || error(
        "Unknown tape boundary condition $(boundary_condition).",
    )
    origin = Float64(origin_cm)
    isfinite(origin) || error("Tape origin must be finite.")

    if Ndims == 3
        isnothing(length_cm) && error("Three-dimensional tape geometry requires length_cm.")
        length_value = Float64(length_cm)
        isfinite(length_value) && length_value > 0.0 || error(
            "Tape segment length must be finite and positive.",
        )
    else
        length_value = isnothing(length_cm) ? 0.0 : Float64(length_cm)
    end

    geometry = Geometry()
    set_type(geometry,"cartesian")
    set_dimension(geometry,Ndims)
    for axis in ("x","y","z")[1:Ndims]
        set_boundary_conditions(geometry,string(axis,"-"),boundary)
        set_boundary_conditions(geometry,string(axis,"+"),boundary)
    end

    number_of_layers = length(stack.layers)
    set_number_of_regions(geometry,"x",number_of_layers)
    x_voxels = Int64[
        layer.transport_subdivisions*refinement for layer in stack.layers
    ]
    set_voxels_per_region(geometry,"x",x_voxels)
    set_region_boundaries(geometry,"x",get_layer_boundaries(stack;origin_cm=origin))

    if Ndims ≥ 2
        set_number_of_regions(geometry,"y",1)
        set_voxels_per_region(geometry,"y",[Int64(width_subdivisions)])
        set_region_boundaries(geometry,"y",[-stack.width_cm/2.0,stack.width_cm/2.0])
    end
    if Ndims == 3
        set_number_of_regions(geometry,"z",1)
        set_voxels_per_region(geometry,"z",[Int64(length_subdivisions)])
        set_region_boundaries(geometry,"z",[0.0,length_value])
    end

    set_material_per_region(
        geometry,_layer_material_array(stack,cross_sections,Ndims),
    )
    build(geometry,cross_sections)
    geometry.is_build = true

    layer_ranges,voxel_layers = _layer_voxel_ranges(stack,refinement)
    Ny = Ndims ≥ 2 ? geometry.number_of_voxels["y"] : 1
    Nz = Ndims == 3 ? geometry.number_of_voxels["z"] : 1
    activation_layers = findall(layer -> layer.activation_enabled,stack.layers)
    atomistic_layers = findall(layer -> layer.atomistic_enabled,stack.layers)
    activation_ids = _linear_voxel_ids_for_layers(
        voxel_layers,Ny,Nz,Int64.(activation_layers),
    )
    atomistic_ids = _linear_voxel_ids_for_layers(
        voxel_layers,Ny,Nz,Int64.(atomistic_layers),
    )

    metadata = Dict(
        "geometry_type" => "planar-cartesian",
        "thickness_axis" => "x",
        "width_axis" => Ndims ≥ 2 ? "y" : "not-resolved",
        "length_axis" => Ndims == 3 ? "z" : "not-resolved",
        "curvature" => "not-included",
        "thickness_refinement" => string(refinement),
    )
    return Tape_Geometry_Map(
        stack,geometry,layer_ranges,voxel_layers,activation_ids,atomistic_ids,metadata,
    )
end

get_activation_voxel_ids(this::Tape_Geometry_Map) = copy(this.activation_voxel_ids)
get_atomistic_voxel_ids(this::Tape_Geometry_Map) = copy(this.atomistic_voxel_ids)

function get_layer_voxel_ids(this::Tape_Geometry_Map,layer_index::Integer)
    1 ≤ layer_index ≤ length(this.stack.layers) || error("Layer index is outside the tape stack.")
    Nx = this.geometry.number_of_voxels["x"]
    Ny = get_dimension(this.geometry) ≥ 2 ? this.geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(this.geometry) == 3 ? this.geometry.number_of_voxels["z"] : 1
    ids = Int64[]
    for cartesian in CartesianIndices((Nx,Ny,Nz))
        this.voxel_layer_index_x[cartesian[1]] == layer_index || continue
        push!(ids,LinearIndices((Nx,Ny,Nz))[cartesian])
    end
    return ids
end

function _canonical_voxel_array(this::Tape_Geometry_Map,values)
    Nx = this.geometry.number_of_voxels["x"]
    Ny = get_dimension(this.geometry) ≥ 2 ? this.geometry.number_of_voxels["y"] : 1
    Nz = get_dimension(this.geometry) == 3 ? this.geometry.number_of_voxels["z"] : 1
    if size(values) == (Nx,Ny,Nz)
        return Float64.(values)
    elseif get_dimension(this.geometry) == 1 && size(values) == (Nx,)
        return reshape(Float64.(values),Nx,1,1)
    elseif get_dimension(this.geometry) == 2 && size(values) == (Nx,Ny)
        return reshape(Float64.(values),Nx,Ny,1)
    end
    error("Response array shape $(size(values)) does not match tape geometry $(Nx,Ny,Nz).")
end

function _voxel_volume(this::Tape_Geometry_Map,ix::Int64,iy::Int64,iz::Int64)
    dimension = get_dimension(this.geometry)
    dimension == 1 && return this.geometry.volume_per_voxel[ix]
    dimension == 2 && return this.geometry.volume_per_voxel[ix,iy]
    return this.geometry.volume_per_voxel[ix,iy,iz]
end

"""
    aggregate_layer_response(map, values; reduction=:volume_average)

Reduce a voxel response to one value per physical tape layer. `:sum` preserves extensive values,
`:maximum` protects local hotspots, and `:volume_average` computes a volume-weighted average of an
intensive field.
"""
function aggregate_layer_response(
    this::Tape_Geometry_Map,
    values;
    reduction::Symbol=:volume_average,
)
    reduction in (:sum,:maximum,:volume_average) || error(
        "Unknown layer response reduction $(reduction).",
    )
    array = _canonical_voxel_array(this,values)
    any(x -> !isfinite(x),array) && error("Layer response contains non-finite values.")
    result = zeros(Float64,length(this.stack.layers))
    denominator = zeros(Float64,length(this.stack.layers))
    if reduction == :maximum
        result .= -Inf
    end

    for cartesian in CartesianIndices(array)
        ix,iy,iz = Tuple(cartesian)
        layer = this.voxel_layer_index_x[ix]
        value = array[cartesian]
        if reduction == :sum
            result[layer] += value
        elseif reduction == :maximum
            result[layer] = max(result[layer],value)
        else
            volume = _voxel_volume(this,ix,iy,iz)
            result[layer] += value*volume
            denominator[layer] += volume
        end
    end
    if reduction == :volume_average
        all(denominator .> 0.0) || error("A tape layer has zero represented volume.")
        result ./= denominator
    end
    return result
end
