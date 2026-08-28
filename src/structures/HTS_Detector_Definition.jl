"""Neutron or photon converter layer placed next to an HTS detector active region."""
struct Detector_Converter_Layer
    name::String
    material_tag::String
    thickness_cm::Float64
    reaction_channel::String
    reaction_products::Vector{String}
    transport_owner::String

    function Detector_Converter_Layer(
        name::AbstractString,
        material_tag::AbstractString,
        thickness_cm::Real,
        reaction_channel::AbstractString;
        reaction_products::AbstractVector{<:AbstractString}=String[],
        transport_owner::AbstractString="OpenMC/Geant4",
    )
        isempty(name) && error("Detector-converter name cannot be empty.")
        isempty(material_tag) && error("Detector-converter material tag cannot be empty.")
        thickness = Float64(thickness_cm)
        isfinite(thickness) && thickness > 0.0 || error(
            "Detector-converter thickness must be finite and positive.",
        )
        isempty(reaction_channel) && error("Detector reaction channel cannot be empty.")
        isempty(transport_owner) && error("Detector converter requires a transport owner.")
        products = String[string(product) for product in reaction_products]
        return new(
            String(name),String(material_tag),thickness,String(reaction_channel),products,
            String(transport_owner),
        )
    end
end

"""
    Detector_Thermal_Model

Optional lumped electrothermal parameters for a bolometric or transition-edge detector response.
The deterministic transport solve supplies deposited-energy impulses or rates; this model carries
the separate thermal/electrical state needed to predict pulse amplitude, decay time, or resistance.
"""
struct Detector_Thermal_Model
    active_heat_capacity_J_K::Union{Nothing,Float64}
    substrate_heat_capacity_J_K::Union{Nothing,Float64}
    active_substrate_conductance_W_K::Union{Nothing,Float64}
    substrate_bath_conductance_W_K::Union{Nothing,Float64}
    transition_temperature_K::Union{Nothing,Float64}
    transition_width_K::Union{Nothing,Float64}
    normal_resistance_ohm::Union{Nothing,Float64}
    inductance_H::Union{Nothing,Float64}

    function Detector_Thermal_Model(;
        active_heat_capacity_J_K::Union{Nothing,Real}=nothing,
        substrate_heat_capacity_J_K::Union{Nothing,Real}=nothing,
        active_substrate_conductance_W_K::Union{Nothing,Real}=nothing,
        substrate_bath_conductance_W_K::Union{Nothing,Real}=nothing,
        transition_temperature_K::Union{Nothing,Real}=nothing,
        transition_width_K::Union{Nothing,Real}=nothing,
        normal_resistance_ohm::Union{Nothing,Real}=nothing,
        inductance_H::Union{Nothing,Real}=nothing,
    )
        values = (
            active_heat_capacity_J_K,
            substrate_heat_capacity_J_K,
            active_substrate_conductance_W_K,
            substrate_bath_conductance_W_K,
            transition_temperature_K,
            transition_width_K,
            normal_resistance_ohm,
            inductance_H,
        )
        converted = map(value -> isnothing(value) ? nothing : Float64(value),values)
        for value in converted
            if !isnothing(value)
                isfinite(value) && value > 0.0 || error(
                    "Specified detector electrothermal parameters must be finite and positive.",
                )
            end
        end
        return new(converted...)
    end
end

"""
    HTS_Detector_Definition

Geometry, converter, operating state, and response contract for using an HTS film or tape as a
radiation detector. Transport and electrothermal response remain separate: Radiant supplies
layer- and process-resolved deposition, while a transient response model consumes the deposition
impulse together with the thermal and electrical parameters.
"""
struct HTS_Detector_Definition
    name::String
    mode::Symbol
    active_material_tag::String
    operating_temperature_K::Float64
    bias_current_A::Float64
    active_width_cm::Float64
    active_length_cm::Float64
    active_thickness_cm::Float64
    substrate_material_tag::String
    converter_layers::Vector{Detector_Converter_Layer}
    thermal_model::Detector_Thermal_Model
    responses::Vector{Symbol}
    metadata::Dict{String,String}

    function HTS_Detector_Definition(
        name::AbstractString,
        mode::Symbol,
        active_material_tag::AbstractString;
        operating_temperature_K::Real,
        bias_current_A::Real=0.0,
        active_width_cm::Real,
        active_length_cm::Real,
        active_thickness_cm::Real,
        substrate_material_tag::AbstractString,
        converter_layers::AbstractVector{Detector_Converter_Layer}=Detector_Converter_Layer[],
        thermal_model::Detector_Thermal_Model=Detector_Thermal_Model(),
        responses::AbstractVector{Symbol}=Symbol[
            :converter_reaction_rate,
            :reaction_product_escape_probability,
            :active_layer_deposition_distribution,
            :prompt_lattice_heat,
            :ionization_excitation,
            :pulse_height,
            :pulse_decay_time,
            :resistance_change,
        ],
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(name) && error("HTS detector name cannot be empty.")
        mode in (:transition_edge,:bolometric,:resistive,:critical_current_shift) || error(
            "Unknown HTS detector mode: $(mode).",
        )
        isempty(active_material_tag) && error("Active HTS material tag cannot be empty.")
        isempty(substrate_material_tag) && error("Detector substrate tag cannot be empty.")
        temperature = Float64(operating_temperature_K)
        bias = Float64(bias_current_A)
        dimensions = Float64[active_width_cm,active_length_cm,active_thickness_cm]
        isfinite(temperature) && temperature > 0.0 || error(
            "Detector operating temperature must be finite and positive.",
        )
        isfinite(bias) && bias ≥ 0.0 || error("Detector bias current must be finite and nonnegative.")
        any(value -> !isfinite(value) || value ≤ 0.0,dimensions) && error(
            "Detector active dimensions must be finite and positive.",
        )
        response_vector = Symbol[responses...]
        isempty(response_vector) && error("At least one detector response is required.")
        length(unique(response_vector)) == length(response_vector) || error(
            "Detector response identifiers must be unique.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(name),mode,String(active_material_tag),temperature,bias,dimensions[1],
            dimensions[2],dimensions[3],String(substrate_material_tag),
            Detector_Converter_Layer[converter_layers...],thermal_model,response_vector,
            metadata_string,
        )
    end
end

function is_ready_for_transient_response(this::HTS_Detector_Definition)
    model = this.thermal_model
    required = (
        model.active_heat_capacity_J_K,
        model.substrate_heat_capacity_J_K,
        model.active_substrate_conductance_W_K,
        model.substrate_bath_conductance_W_K,
        model.transition_temperature_K,
        model.transition_width_K,
        model.normal_resistance_ohm,
    )
    return all(value -> !isnothing(value),required)
end

get_active_volume(this::HTS_Detector_Definition) =
    this.active_width_cm*this.active_length_cm*this.active_thickness_cm

"""Create a geometry/schema fixture for a B-10-enriched B4C/YBCO transition-edge detector."""
function verification_b10_ybco_detector()
    converter = Detector_Converter_Layer(
        "B10 converter","B10-enriched-B4C",1.0e-4,"B-10(n,alpha)Li-7";
        reaction_products=["alpha","Li-7","gamma"],
        transport_owner="OpenMC/Geant4",
    )
    return HTS_Detector_Definition(
        "B10-YBCO verification detector",:transition_edge,"YBCO";
        operating_temperature_K=77.0,
        bias_current_A=0.0,
        active_width_cm=0.1,
        active_length_cm=0.1,
        active_thickness_cm=1.0e-5,
        substrate_material_tag="sapphire",
        converter_layers=[converter],
        metadata=Dict(
            "classification" => "verification-only",
            "thermal_parameters" => "unbound",
        ),
    )
end
