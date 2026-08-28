"""One concentric material region in a fibre-optic diagnostic package."""
struct Fiber_Optic_Layer
    role::Symbol
    material_tag::String
    outer_radius_cm::Float64
    transport_subdivisions::Int64

    function Fiber_Optic_Layer(
        role::Symbol,
        material_tag::AbstractString,
        outer_radius_cm::Real;
        transport_subdivisions::Integer=1,
    )
        role ∈ (:core,:cladding,:coating,:adhesive,:capillary,:jacket,:other) || error(
            "Unsupported fibre layer role $(role).",
        )
        isempty(material_tag) && error("Fibre layer material tag cannot be empty.")
        radius = Float64(outer_radius_cm)
        isfinite(radius) && radius > 0.0 || error("Fibre layer radius must be finite and positive.")
        subdivisions = Int64(transport_subdivisions)
        subdivisions ≥ 1 || error("Fibre transport subdivisions must be at least one.")
        return new(role,String(material_tag),radius,subdivisions)
    end
end

"""
    Fiber_Optic_Definition

Explicit concentric geometry and response contract for fibre diagnostics embedded near an HTS
tape. Radiant supplies ionizing-particle transport responses. Optical attenuation, colour-centre
kinetics, radioluminescence, Bragg-wavelength drift, and annealing remain owned by an external
optical-material model and must not be inferred directly from deposited energy.
"""
struct Fiber_Optic_Definition
    tag::String
    layers::Vector{Fiber_Optic_Layer}
    diagnostic_modes::Vector{Symbol}
    protected_responses::Vector{Symbol}
    optical_response_owner::String
    metadata::Dict{String,String}

    function Fiber_Optic_Definition(
        tag::AbstractString,
        layers::AbstractVector{Fiber_Optic_Layer};
        diagnostic_modes::AbstractVector{Symbol}=Symbol[:fbg],
        protected_responses::AbstractVector{Symbol}=Symbol[
            :ionizing_energy_deposition,
            :electron_spectrum,
            :charge_deposition,
            :nuclear_recoil_source,
            :activation_source,
            :temperature_source,
        ],
        optical_response_owner::AbstractString="external-optical-material-model",
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(tag) && error("Fibre definition tag cannot be empty.")
        isempty(layers) && error("A fibre definition requires at least one material layer.")
        roles = getfield.(layers,:role)
        count(==( :core),roles) == 1 || error("A fibre definition requires exactly one core layer.")
        radii = getfield.(layers,:outer_radius_cm)
        all(radii[index] < radii[index+1] for index in 1:length(radii)-1) || error(
            "Fibre layer radii must increase strictly from the core outward.",
        )
        allowed_modes = (:fbg,:distributed_rayleigh,:distributed_raman,:distributed_brillouin,
                         :radioluminescence,:dosimetry,:other)
        all(mode in allowed_modes for mode in diagnostic_modes) || error("Unsupported fibre diagnostic mode.")
        isempty(optical_response_owner) && error("The optical-response owner must be explicit.")
        metadata_string = Dict{String,String}(string(key)=>string(value) for (key,value) in metadata)
        return new(
            String(tag),collect(layers),collect(diagnostic_modes),collect(protected_responses),
            String(optical_response_owner),metadata_string,
        )
    end
end

get_fiber_outer_radius(this::Fiber_Optic_Definition) = this.layers[end].outer_radius_cm

function get_fiber_layer_index(this::Fiber_Optic_Definition,radius_cm::Real)
    radius = Float64(radius_cm)
    isfinite(radius) && radius ≥ 0.0 || error("Fibre radius query must be finite and nonnegative.")
    radius ≤ get_fiber_outer_radius(this) || error("Fibre radius query is outside the diagnostic package.")
    index = findfirst(layer -> radius ≤ layer.outer_radius_cm,this.layers)
    isnothing(index) && error("Unable to identify a fibre layer.")
    return index
end
