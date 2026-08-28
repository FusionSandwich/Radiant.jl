"""
    Layer_Definition

Physical layer in an HTS tape, insulation stack, sensor, or other local heterogeneous structure.

`transport_subdivisions` is a minimum requested subdivision count; convergence studies may refine
it further. Activation and atomistic flags identify layers for which downstream reaction-rate or
event-kernel products are required.
"""
struct Layer_Definition
    name::String
    material_tag::String
    thickness_cm::Float64
    role::Symbol
    transport_subdivisions::Int64
    activation_enabled::Bool
    atomistic_enabled::Bool
    metadata::Dict{String,String}

    function Layer_Definition(
        name::AbstractString,
        material_tag::AbstractString,
        thickness_cm::Real;
        role::Symbol = :structural,
        transport_subdivisions::Integer = 1,
        activation_enabled::Bool = false,
        atomistic_enabled::Bool = false,
        metadata::AbstractDict = Dict{String,String}(),
    )
        thickness = Float64(thickness_cm)
        subdivisions = Int64(transport_subdivisions)
        if isempty(name) || isempty(material_tag)
            error("Layer name and material tag cannot be empty.")
        end
        if !isfinite(thickness) || thickness ≤ 0.0
            error("Layer thickness must be finite and greater than zero.")
        end
        if subdivisions < 1
            error("Layer transport subdivisions must be at least one.")
        end

        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end

        return new(
            String(name),
            String(material_tag),
            thickness,
            role,
            subdivisions,
            activation_enabled,
            atomistic_enabled,
            metadata_string,
        )
    end
end
