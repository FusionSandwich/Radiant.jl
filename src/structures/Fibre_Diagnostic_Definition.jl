"""One concentric material layer of a fibre-optic diagnostic assembly."""
struct Fibre_Radial_Layer
    name::String
    material_tag::String
    inner_radius_cm::Float64
    outer_radius_cm::Float64
    transport_subdivisions::Int64
    optical_role::Symbol

    function Fibre_Radial_Layer(
        name::AbstractString,
        material_tag::AbstractString,
        inner_radius_cm::Real,
        outer_radius_cm::Real;
        transport_subdivisions::Integer=1,
        optical_role::Symbol=:passive,
    )
        isempty(name) && error("Fibre-layer name cannot be empty.")
        isempty(material_tag) && error("Fibre-layer material tag cannot be empty.")
        inner = Float64(inner_radius_cm)
        outer = Float64(outer_radius_cm)
        isfinite(inner) && isfinite(outer) || error("Fibre radii must be finite.")
        inner ≥ 0.0 || error("Fibre inner radius must be nonnegative.")
        outer > inner || error("Fibre outer radius must exceed its inner radius.")
        transport_subdivisions ≥ 1 || error("Each fibre layer requires at least one subdivision.")
        optical_role in (:core,:cladding,:coating,:adhesive,:capillary,:passive) || error(
            "Unknown fibre optical role: $(optical_role).",
        )
        return new(
            String(name),String(material_tag),inner,outer,Int64(transport_subdivisions),
            optical_role,
        )
    end
end

"""
    Fibre_Diagnostic_Definition

Concentric fibre, coating, adhesive, and capillary definition used to generate explicit local
transport regions and a protected-response list. Optical degradation or wavelength shift is not
silently inferred from ionizing dose; those effects remain named downstream response models.
"""
struct Fibre_Diagnostic_Definition
    name::String
    mode::Symbol
    layers::Vector{Fibre_Radial_Layer}
    operating_temperature_K::Float64
    interrogation_wavelength_nm::Union{Nothing,Float64}
    responses::Vector{Symbol}
    metadata::Dict{String,String}

    function Fibre_Diagnostic_Definition(
        name::AbstractString,
        mode::Symbol,
        layers::AbstractVector{Fibre_Radial_Layer};
        operating_temperature_K::Real=293.15,
        interrogation_wavelength_nm::Union{Nothing,Real}=nothing,
        responses::AbstractVector{Symbol}=Symbol[
            :ionizing_dose,
            :neutron_scalar_flux,
            :capture_rate,
            :displacement_source,
            :secondary_electron_spectrum,
            :hydrogen_helium_production,
            :optical_emission_source,
            :radiation_induced_attenuation_state,
        ],
        metadata::AbstractDict=Dict{String,String}(),
        interface_tolerance_cm::Real=1.0e-12,
    )
        isempty(name) && error("Fibre diagnostic name cannot be empty.")
        mode in (:dosimetry,:fbg,:luminescence,:distributed_temperature,:ria) || error(
            "Unknown fibre diagnostic mode: $(mode).",
        )
        isempty(layers) && error("At least one fibre radial layer is required.")
        temperature = Float64(operating_temperature_K)
        isfinite(temperature) && temperature > 0.0 || error(
            "Fibre operating temperature must be finite and positive.",
        )
        tolerance = Float64(interface_tolerance_cm)
        tolerance ≥ 0.0 || error("Fibre interface tolerance must be nonnegative.")

        layer_vector = Fibre_Radial_Layer[layers...]
        abs(layer_vector[1].inner_radius_cm) ≤ tolerance || error(
            "The innermost fibre layer must start at radius zero.",
        )
        for index in 2:length(layer_vector)
            isapprox(
                layer_vector[index-1].outer_radius_cm,
                layer_vector[index].inner_radius_cm;
                rtol=0.0,atol=tolerance,
            ) || error("Fibre radial layers must be contiguous and non-overlapping.")
        end
        names = getfield.(layer_vector,:name)
        length(unique(names)) == length(names) || error("Fibre layer names must be unique.")

        wavelength = isnothing(interrogation_wavelength_nm) ? nothing :
            Float64(interrogation_wavelength_nm)
        if !isnothing(wavelength)
            isfinite(wavelength) && wavelength > 0.0 || error(
                "Interrogation wavelength must be finite and positive.",
            )
        end
        isempty(responses) && error("At least one fibre response must be protected.")
        response_vector = Symbol[responses...]
        length(unique(response_vector)) == length(response_vector) || error(
            "Fibre response identifiers must be unique.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(name),mode,layer_vector,temperature,wavelength,response_vector,
            metadata_string,
        )
    end
end

get_outer_radius(this::Fibre_Diagnostic_Definition) = this.layers[end].outer_radius_cm
get_cross_sectional_area(this::Fibre_Diagnostic_Definition) = π*get_outer_radius(this)^2

function get_fibre_layer_index(
    this::Fibre_Diagnostic_Definition,
    radius_cm::Real;
    tolerance_cm::Real=1.0e-12,
)
    radius = Float64(radius_cm)
    radius ≥ -tolerance_cm || error("Fibre radius lies outside the assembly.")
    radius ≤ get_outer_radius(this)+tolerance_cm || error("Fibre radius lies outside the assembly.")
    for (index,layer) in enumerate(this.layers)
        if radius ≥ layer.inner_radius_cm-tolerance_cm &&
           radius ≤ layer.outer_radius_cm+tolerance_cm
            return index
        end
    end
    error("Fibre radius could not be assigned to a radial layer.")
end

function requires_optical_response_model(this::Fibre_Diagnostic_Definition)
    return any(
        response -> response in (
            :optical_emission_source,
            :radiation_induced_attenuation_state,
            :fbg_wavelength_shift,
        ),
        this.responses,
    )
end
