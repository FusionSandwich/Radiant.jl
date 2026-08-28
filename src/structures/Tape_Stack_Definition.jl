"""
    Tape_Stack_Definition

Ordered front-to-back definition of an HTS tape or layered local structure.
"""
struct Tape_Stack_Definition
    name::String
    width_cm::Float64
    layers::Vector{Layer_Definition}
    metadata::Dict{String,String}

    function Tape_Stack_Definition(
        name::AbstractString,
        width_cm::Real,
        layers::AbstractVector{Layer_Definition};
        metadata::AbstractDict = Dict{String,String}(),
    )
        width = Float64(width_cm)
        layer_vector = collect(layers)
        if isempty(name)
            error("Tape-stack name cannot be empty.")
        end
        if !isfinite(width) || width ≤ 0.0
            error("Tape-stack width must be finite and greater than zero.")
        end
        if isempty(layer_vector)
            error("Tape stack must contain at least one layer.")
        end
        names = [layer.name for layer in layer_vector]
        if length(unique(names)) != length(names)
            error("Layer names must be unique within a tape stack.")
        end

        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end

        return new(String(name),width,layer_vector,metadata_string)
    end
end

"""
    get_total_thickness(this::Tape_Stack_Definition)

Return total stack thickness in centimetres.
"""
function get_total_thickness(this::Tape_Stack_Definition)
    return sum(layer.thickness_cm for layer in this.layers)
end

"""
    get_layer_boundaries(this::Tape_Stack_Definition; origin_cm=0.0)

Return front-to-back layer boundaries in centimetres.
"""
function get_layer_boundaries(this::Tape_Stack_Definition; origin_cm::Real=0.0)
    boundaries = zeros(Float64,length(this.layers)+1)
    boundaries[1] = Float64(origin_cm)
    for index in eachindex(this.layers)
        boundaries[index+1] = boundaries[index] + this.layers[index].thickness_cm
    end
    return boundaries
end

"""
    get_layer_index(this::Tape_Stack_Definition, depth_cm; origin_cm=0.0)

Return the one-based layer index containing `depth_cm`. The rear-most external boundary belongs to
the last layer. An error is raised for points outside the stack.
"""
function get_layer_index(this::Tape_Stack_Definition,depth_cm::Real; origin_cm::Real=0.0)
    depth = Float64(depth_cm)
    boundaries = get_layer_boundaries(this;origin_cm=origin_cm)
    if depth < boundaries[1] || depth > boundaries[end]
        error("Depth lies outside the tape stack.")
    end
    if depth == boundaries[end]
        return length(this.layers)
    end
    index = searchsortedlast(boundaries,depth)
    return min(index,length(this.layers))
end

"""
    verification_eight_layer_stack(; rebco_material_tag="YBCO")

Return the canonical 4 mm wide, 163.2 micrometre eight-layer verification fixture. This object is
explicitly marked verification-only and must not be substituted for a measured production tape.
"""
function verification_eight_layer_stack(; rebco_material_tag::AbstractString="YBCO")
    micrometre_to_cm = 1.0e-4
    layers = Layer_Definition[
        Layer_Definition("copper-front","Cu",20.0*micrometre_to_cm;role=:stabilizer,transport_subdivisions=2),
        Layer_Definition("silver-cap","Ag",2.0*micrometre_to_cm;role=:cap,transport_subdivisions=2),
        Layer_Definition("rebco",rebco_material_tag,1.0*micrometre_to_cm;role=:superconductor,transport_subdivisions=4,activation_enabled=true,atomistic_enabled=true),
        Layer_Definition("buffer","buffer",0.2*micrometre_to_cm;role=:buffer,transport_subdivisions=2),
        Layer_Definition("hastelloy","Hastelloy",50.0*micrometre_to_cm;role=:substrate,transport_subdivisions=2),
        Layer_Definition("copper-rear","Cu",20.0*micrometre_to_cm;role=:stabilizer,transport_subdivisions=2),
        Layer_Definition("solder","solder",20.0*micrometre_to_cm;role=:joint,transport_subdivisions=2),
        Layer_Definition("insulation","insulation",50.0*micrometre_to_cm;role=:insulation,transport_subdivisions=2),
    ]
    return Tape_Stack_Definition(
        "verification-eight-layer-163p2um",
        0.4,
        layers;
        metadata=Dict(
            "classification" => "verification-only",
            "coordinate_order" => "front-to-back",
            "nominal_total_thickness_um" => "163.2",
        ),
    )
end
