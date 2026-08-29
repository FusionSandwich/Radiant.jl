abstract type Abstract_Material_Property_Model end

"""NIST log-polynomial fit: log10(y) = Σ aᵢ(log10(T))^(i-1)."""
struct NIST_Log10_Polynomial <: Abstract_Material_Property_Model
    coefficients::Vector{Float64}
end

"""NIST OFHC-copper conductivity rational fit in powers of sqrt(T)."""
struct NIST_Copper_Conductivity <: Abstract_Material_Property_Model
    coefficients::NTuple{9,Float64}
    rrr::Int64
end

"""Reference exists but its numerical table has not been digitized or generated."""
struct Source_Only_Material_Property <: Abstract_Material_Property_Model
    required_generation_path::String
end

const MATERIAL_PROPERTY_STATUSES = (
    :literature_fit,
    :literature_table_pending,
    :atomistic_candidate,
    :qualified,
)

struct Material_Property_Record
    record_id::String
    material_tag::String
    property::Symbol
    direction::Symbol
    state_description::String
    units::String
    temperature_min_K::Float64
    temperature_max_K::Float64
    nominal_relative_uncertainty::Float64
    model::Abstract_Material_Property_Model
    source_title::String
    source_identifier::String
    source_url::String
    status::Symbol
    metadata::Dict{String,String}

    function Material_Property_Record(;
        record_id::AbstractString,
        material_tag::AbstractString,
        property::Symbol,
        direction::Symbol=:isotropic,
        state_description::AbstractString,
        units::AbstractString,
        temperature_min_K::Real,
        temperature_max_K::Real,
        nominal_relative_uncertainty::Real,
        model::Abstract_Material_Property_Model,
        source_title::AbstractString,
        source_identifier::AbstractString,
        source_url::AbstractString,
        status::Symbol,
        metadata::AbstractDict=Dict{String,String}(),
    )
        for (label,value) in (
            ("record ID",record_id),("material tag",material_tag),
            ("state description",state_description),("units",units),
            ("source title",source_title),("source identifier",source_identifier),
            ("source URL",source_url),
        )
            isempty(value) && error("Material-property $(label) cannot be empty.")
        end
        property in (:thermal_conductivity,:specific_heat,:volumetric_heat_capacity,
                     :interface_conductance,:relaxation_time) || error(
            "Unsupported material property $(property).",
        )
        lower = Float64(temperature_min_K)
        upper = Float64(temperature_max_K)
        uncertainty = Float64(nominal_relative_uncertainty)
        isfinite(lower) && isfinite(upper) && lower > 0.0 && upper > lower || error(
            "Material-property temperature range must be finite, positive, and increasing.",
        )
        isfinite(uncertainty) && uncertainty >= 0.0 || error(
            "Material-property uncertainty must be finite and nonnegative.",
        )
        status in MATERIAL_PROPERTY_STATUSES || error(
            "Unknown material-property status $(status).",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(record_id),String(material_tag),property,direction,String(state_description),
            String(units),lower,upper,uncertainty,model,String(source_title),
            String(source_identifier),String(source_url),status,metadata_string,
        )
    end
end

struct Material_Response_Registry
    schema::String
    records::Vector{Material_Property_Record}

    function Material_Response_Registry(records::AbstractVector{Material_Property_Record})
        record_vector = Material_Property_Record[records...]
        isempty(record_vector) && error("Material response registry cannot be empty.")
        identifiers = getfield.(record_vector,:record_id)
        length(unique(identifiers)) == length(identifiers) || error(
            "Material-property record identifiers must be unique.",
        )
        return new("radiant.hts.material_response_registry/v1",record_vector)
    end
end

function material_property_value(
    model::NIST_Log10_Polynomial,
    temperature_K::Real,
)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error(
        "Material-property temperature must be finite and positive.",
    )
    x = log10(temperature)
    exponent = 0.0
    for (index,coefficient) in enumerate(model.coefficients)
        exponent += coefficient*x^(index-1)
    end
    return 10.0^exponent
end

function material_property_value(
    model::NIST_Copper_Conductivity,
    temperature_K::Real,
)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error(
        "Material-property temperature must be finite and positive.",
    )
    a,b,c,d,e,f,g,h,i = model.coefficients
    root = sqrt(temperature)
    numerator = a+c*root+e*temperature+g*temperature*root+i*temperature^2
    denominator = 1.0+b*root+d*temperature+f*temperature*root+h*temperature^2
    denominator != 0.0 || error("NIST copper-conductivity fit has a zero denominator.")
    return 10.0^(numerator/denominator)
end

function material_property_value(
    model::Source_Only_Material_Property,
    temperature_K::Real,
)
    error(
        "This material property is source-only. Required generation path: " *
        model.required_generation_path,
    )
end

function material_property_value(
    record::Material_Property_Record,
    temperature_K::Real;
    allow_extrapolation::Bool=false,
)
    temperature = Float64(temperature_K)
    if !allow_extrapolation &&
       (temperature < record.temperature_min_K || temperature > record.temperature_max_K)
        error(
            "Temperature $(temperature) K lies outside the supported range for " *
            record.record_id * ".",
        )
    end
    return material_property_value(record.model,temperature)
end

function get_material_property(
    registry::Material_Response_Registry;
    material_tag::AbstractString,
    property::Symbol,
    direction::Symbol=:isotropic,
    state_contains::Union{Nothing,AbstractString}=nothing,
)
    matches = Material_Property_Record[]
    for record in registry.records
        record.material_tag == material_tag || continue
        record.property == property || continue
        record.direction == direction || continue
        if !isnothing(state_contains)
            occursin(lowercase(String(state_contains)),lowercase(record.state_description)) ||
                continue
        end
        push!(matches,record)
    end
    isempty(matches) && error("No material-property record matches the requested state.")
    length(matches) == 1 || error(
        "Multiple material-property records match; specify state_contains more precisely.",
    )
    return matches[1]
end

function tabulate_material_property(
    record::Material_Property_Record,
    temperatures_K::AbstractVector{<:Real},
)
    temperatures = Float64.(temperatures_K)
    values = [material_property_value(record,temperature) for temperature in temperatures]
    uncertainty = abs.(values).*record.nominal_relative_uncertainty
    return Dict{String,Any}(
        "schema" => "radiant.hts.material_property_table/v1",
        "record_id" => record.record_id,
        "material_tag" => record.material_tag,
        "property" => string(record.property),
        "direction" => string(record.direction),
        "state_description" => record.state_description,
        "units" => record.units,
        "temperature_K" => temperatures,
        "value" => values,
        "standard_uncertainty" => uncertainty,
        "status" => string(record.status),
        "source_title" => record.source_title,
        "source_identifier" => record.source_identifier,
        "source_url" => record.source_url,
    )
end

function cryogenic_property_from_record(
    record::Material_Property_Record,
    temperatures_K::AbstractVector{<:Real},
)
    record.status != :literature_table_pending || error(
        "Cannot build a numerical cryogenic property from an undigitized source-only record.",
    )
    table = tabulate_material_property(record,temperatures_K)
    return Tabulated_Cryogenic_Property(
        table["temperature_K"],table["value"];
        units=record.units,
        property_hash=record.source_identifier * ":" * record.record_id,
        qualification_status=record.status == :qualified ? :qualified : :candidate,
        allow_zero=false,
    )
end

_nist_log(coefficients) = NIST_Log10_Polynomial(Float64.(coefficients))
_nist_copper(rrr,coefficients) = NIST_Copper_Conductivity(Tuple(Float64.(coefficients)),Int64(rrr))

function _source_only_record(;
    record_id,
    material_tag,
    property,
    direction=:isotropic,
    state_description,
    units,
    temperature_min_K,
    temperature_max_K,
    source_title,
    source_identifier,
    source_url,
    required_generation_path,
    metadata=Dict{String,String}(),
)
    return Material_Property_Record(
        record_id=record_id,material_tag=material_tag,property=property,
        direction=direction,state_description=state_description,units=units,
        temperature_min_K=temperature_min_K,temperature_max_K=temperature_max_K,
        nominal_relative_uncertainty=0.0,
        model=Source_Only_Material_Property(required_generation_path),
        source_title=source_title,source_identifier=source_identifier,source_url=source_url,
        status=:literature_table_pending,metadata=metadata,
    )
end

"""
Return source-backed candidate material models. Numerical fits are included only where an official
NIST equation and coefficients are directly available. YBCO, GdBCO, Hastelloy, Ag, and complete
tape records remain source-only until tables are digitized or generated and hash-bound.
"""
function default_material_response_registry()
    records = Material_Property_Record[]
    copper_coefficients = Dict(
        50 => [1.8743,-0.41538,-0.6018,0.13294,0.26426,-0.0219,-0.051276,0.0014871,0.003723],
        100 => [2.2154,-0.47461,-0.88068,0.13871,0.29505,-0.02043,-0.04831,0.001281,0.003207],
        150 => [2.3797,-0.4918,-0.98615,0.13942,0.30475,-0.019713,-0.046897,0.0011969,0.0029988],
        300 => [1.357,0.3981,2.669,-0.1346,-0.6683,0.01342,0.05773,0.0002147,0.0],
        500 => [2.8075,-0.54074,-1.2777,0.15362,0.36444,-0.02105,-0.051727,0.0012226,0.0030964],
    )
    copper_errors = Dict(50 => 0.02,100 => 0.01,150 => 0.02,300 => 0.01,500 => 0.02)
    for rrr in sort(collect(keys(copper_coefficients)))
        push!(records,Material_Property_Record(
            record_id="NIST-OFHC-Cu-k-RRR$(rrr)",material_tag="Cu-OFHC",
            property=:thermal_conductivity,direction=:isotropic,
            state_description="OFHC copper, RRR=$(rrr)",units="W/(m*K)",
            temperature_min_K=4.0,temperature_max_K=300.0,
            nominal_relative_uncertainty=copper_errors[rrr],
            model=_nist_copper(rrr,copper_coefficients[rrr]),
            source_title="NIST cryogenic material properties: OFHC Copper",
            source_identifier="NIST-MONO-177:OFHC-k-RRR$(rrr)",
            source_url="https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm",
            status=:literature_fit,
            metadata=Dict("fit_form" => "NIST copper rational log10 fit"),
        ))
    end
    push!(records,Material_Property_Record(
        record_id="NIST-OFHC-Cu-cp",material_tag="Cu-OFHC",property=:specific_heat,
        state_description="OFHC copper",units="J/(kg*K)",temperature_min_K=4.0,
        temperature_max_K=300.0,nominal_relative_uncertainty=0.10,
        model=_nist_log([-1.91844,-0.15973,8.61013,-18.996,21.9661,-12.7328,3.54322,-0.3797,0.0]),
        source_title="NIST cryogenic material properties: OFHC Copper",
        source_identifier="NIST-MONO-177:OFHC-cp",
        source_url="https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm",
        status=:literature_fit,
        metadata=Dict("uncertainty_note" => "10% below 15 K; 5% at and above 15 K"),
    ))
    for (record_id,material_tag,property,direction,coefficients,error_fraction) in (
        ("NIST-Kapton-k","Kapton",:thermal_conductivity,:isotropic,[5.73101,-39.5199,79.9313,-83.8572,50.9157,-17.9835,3.42413,-0.27133,0.0],0.02),
        ("NIST-Kapton-cp","Kapton",:specific_heat,:isotropic,[-1.3684,0.65892,2.8719,0.42651,-3.0088,1.9558,-0.51998,0.051574,0.0],0.03),
        ("NIST-G10CR-k-normal","G10-CR",:thermal_conductivity,:normal,[-4.1236,13.788,-26.068,26.272,-14.663,4.4954,-0.6905,0.0397,0.0],0.05),
        ("NIST-G10CR-k-warp","G10-CR",:thermal_conductivity,:warp,[-2.64827,8.80228,-24.8998,41.1625,-39.8754,23.1778,-7.95635,1.48806,-0.11701],0.05),
        ("NIST-G10CR-cp","G10-CR",:specific_heat,:isotropic,[-2.4083,7.6006,-8.2982,7.3301,-4.2386,1.4294,-0.24396,0.015236,0.0],0.02),
    )
        lower = record_id == "NIST-G10CR-k-normal" ? 10.0 :
                record_id == "NIST-G10CR-k-warp" ? 12.0 : 4.0
        units = property == :specific_heat ? "J/(kg*K)" : "W/(m*K)"
        source_url = material_tag == "Kapton" ?
            "https://trc.nist.gov/cryogenics/materials/Polyimide%20Kapton/PolyimideKapton_rev.htm" :
            "https://trc.nist.gov/cryogenics/materials/G-10%20CR%20Fiberglass%20Epoxy/G10CRFiberglassEpoxy_rev.htm"
        push!(records,Material_Property_Record(
            record_id=record_id,material_tag=material_tag,property=property,
            direction=direction,state_description=material_tag,units=units,
            temperature_min_K=lower,temperature_max_K=300.0,
            nominal_relative_uncertainty=error_fraction,model=_nist_log(coefficients),
            source_title="NIST cryogenic material properties: $(material_tag)",
            source_identifier=record_id,source_url=source_url,status=:literature_fit,
        ))
    end

    append!(records,[
        _source_only_record(
            record_id="YBCO-coated-conductor-cp",material_tag="YBCO-coated-conductor",
            property=:specific_heat,state_description="Ag/YBCO/Hastelloy and Cu-reinforced coated conductors",
            units="J/(kg*K)",temperature_min_K=4.0,temperature_max_K=300.0,
            source_title="Specific Heat and Thermal Diffusivity of YBCO Coated Conductors",
            source_identifier="doi:10.1016/j.phpro.2012.06.316",
            source_url="https://doi.org/10.1016/j.phpro.2012.06.316",
            required_generation_path="digitize the open-access C(T) curves or reconstruct the layer mass-weighted sum using source-backed constituent tables",
        ),
        _source_only_record(
            record_id="Hastelloy-C276-cryogenic-cp",material_tag="Hastelloy-C276",
            property=:specific_heat,state_description="Hastelloy C-276, 2-300 K measured specimen",
            units="J/(kg*K)",temperature_min_K=2.0,temperature_max_K=300.0,
            source_title="Physical properties of Hastelloy C-276 at cryogenic temperatures",
            source_identifier="doi:10.1063/1.2899058",
            source_url="https://doi.org/10.1063/1.2899058",
            required_generation_path="digitize the published Cp curve and retain specimen/composition metadata",
        ),
        _source_only_record(
            record_id="Hastelloy-C276-cryogenic-k",material_tag="Hastelloy-C276",
            property=:thermal_conductivity,state_description="Hastelloy C-276, 2-200 K measured specimen",
            units="W/(m*K)",temperature_min_K=2.0,temperature_max_K=200.0,
            source_title="Physical properties of Hastelloy C-276 at cryogenic temperatures",
            source_identifier="doi:10.1063/1.2899058",
            source_url="https://doi.org/10.1063/1.2899058",
            required_generation_path="digitize the published thermal-conductivity curve and retain specimen/composition metadata",
        ),
        _source_only_record(
            record_id="REBCO-tape-longitudinal-k-2024",material_tag="REBCO-coated-conductor",
            property=:thermal_conductivity,direction=:longitudinal,
            state_description="Cu, Ag, or Ag-Au stabilizer; 4.2-200 K; RRR-dependent",
            units="W/(m*K)",temperature_min_K=4.2,temperature_max_K=200.0,
            source_title="Thermal conductivity of REBCO tapes with different stabilizers from 4.2 to 200 K",
            source_identifier="arXiv:2410.09915",
            source_url="https://arxiv.org/abs/2410.09915",
            required_generation_path="digitize each tape/stabilizer curve and bind Cu/Ag thickness plus measured RRR",
        ),
        _source_only_record(
            record_id="YBCO-anisotropic-k",material_tag="YBCO",
            property=:thermal_conductivity,direction=:crystal_tensor,
            state_description="YBa2Cu3O7-y single-crystal anisotropic conductivity",
            units="W/(m*K)",temperature_min_K=7.0,temperature_max_K=300.0,
            source_title="Anisotropy of the thermal conductivity of YBa2Cu3O7-y",
            source_identifier="doi:10.1103/PhysRevB.40.9389",
            source_url="https://doi.org/10.1103/PhysRevB.40.9389",
            required_generation_path="digitize a-, b-, and c-axis curves or replace with converged phono3py/Green-Kubo tensors",
        ),
        _source_only_record(
            record_id="GdBCO-thermal-response",material_tag="GdBCO",
            property=:thermal_conductivity,direction=:crystal_tensor,
            state_description="oxygen- and texture-resolved GdBCO film",
            units="W/(m*K)",temperature_min_K=4.0,temperature_max_K=300.0,
            source_title="No transferable production table selected",
            source_identifier="atomistic-generation-required:GdBCO",
            source_url="https://phonopy.github.io/phono3py/",
            required_generation_path="generate oxygen-, texture-, isotope-, and defect-state-resolved tensors with DFT/MLIP plus phono3py or Green-Kubo MD",
        ),
        _source_only_record(
            record_id="YBCO-electron-phonon-relaxation",material_tag="YBCO",
            property=:relaxation_time,state_description="epitaxial microbridge; nonequilibrium optical excitation",
            units="s",temperature_min_K=20.0,temperature_max_K=80.0,
            source_title="Intrinsic picosecond response times of Y-Ba-Cu-O superconducting photodetectors",
            source_identifier="doi:10.1063/1.123388",
            source_url="https://doi.org/10.1063/1.123388",
            required_generation_path="use published values only as detector-model priors; generate state-specific electron-phonon kernels with EPW and calibrate to the actual film",
        ),
    ])
    return Material_Response_Registry(records)
end
