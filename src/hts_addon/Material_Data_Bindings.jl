"""Source-backed tabulated material property with bounded interpolation and uncertainty."""
struct Tabulated_Material_Property <: Abstract_Material_Property_Model
    temperatures_K::Vector{Float64}
    values::Vector{Float64}
    standard_uncertainties::Vector{Float64}
    interpolation::Symbol

    function Tabulated_Material_Property(
        temperatures_K::AbstractVector{<:Real},
        values::AbstractVector{<:Real};
        standard_uncertainties::AbstractVector{<:Real}=zeros(length(values)),
        interpolation::Symbol=:linear,
        allow_zero::Bool=false,
    )
        temperatures = Float64.(temperatures_K)
        data = Float64.(values)
        uncertainties = Float64.(standard_uncertainties)
        length(temperatures) >= 2 || error(
            "A tabulated material property requires at least two temperatures.",
        )
        length(data) == length(temperatures) || error(
            "Tabulated material-property values do not match the temperature axis.",
        )
        length(uncertainties) == length(temperatures) || error(
            "Tabulated material-property uncertainties do not match the temperature axis.",
        )
        all(isfinite,temperatures) && all(temperatures .> 0.0) &&
            all(diff(temperatures) .> 0.0) || error(
            "Tabulated material-property temperatures must be positive and increasing.",
        )
        predicate = allow_zero ? (value -> isfinite(value) && value >= 0.0) :
                                 (value -> isfinite(value) && value > 0.0)
        all(predicate,data) || error("Tabulated material-property values are invalid.")
        all(value -> isfinite(value) && value >= 0.0,uncertainties) || error(
            "Tabulated material-property uncertainties must be nonnegative.",
        )
        interpolation in (:linear,:loglog) || error(
            "Material-property interpolation must be :linear or :loglog.",
        )
        interpolation == :loglog && any(data .<= 0.0) && error(
            "Log-log interpolation requires strictly positive values.",
        )
        return new(temperatures,data,uncertainties,interpolation)
    end
end

function _tabulated_bracket(temperatures::Vector{Float64},temperature::Float64)
    temperature < temperatures[1] && error(
        "Temperature lies below the tabulated material-property domain.",
    )
    temperature > temperatures[end] && error(
        "Temperature lies above the tabulated material-property domain.",
    )
    upper = searchsortedfirst(temperatures,temperature)
    upper == 1 && return 1,1,0.0
    temperature == temperatures[upper] && return upper,upper,0.0
    lower = upper-1
    fraction = (temperature-temperatures[lower])/(temperatures[upper]-temperatures[lower])
    return lower,upper,fraction
end

function material_property_value(model::Tabulated_Material_Property,temperature_K::Real)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error(
        "Material-property query temperature must be finite and positive.",
    )
    lower,upper,fraction = _tabulated_bracket(model.temperatures_K,temperature)
    lower == upper && return model.values[lower]
    if model.interpolation == :linear
        return (1.0-fraction)*model.values[lower]+fraction*model.values[upper]
    end
    return exp(
        (1.0-fraction)*log(model.values[lower])+fraction*log(model.values[upper]),
    )
end

function material_property_standard_uncertainty(
    model::Tabulated_Material_Property,
    temperature_K::Real,
)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error(
        "Material-property uncertainty query temperature must be finite and positive.",
    )
    lower,upper,fraction = _tabulated_bracket(model.temperatures_K,temperature)
    lower == upper && return model.standard_uncertainties[lower]
    return (1.0-fraction)*model.standard_uncertainties[lower]+
           fraction*model.standard_uncertainties[upper]
end

"""Insert or replace one record without changing the registry's existing vector contract."""
function register_material_property!(
    registry::Material_Response_Registry,
    record::Material_Property_Record;
    replace::Bool=false,
)
    matches = findall(existing -> existing.record_id == record.record_id,registry.records)
    length(matches) <= 1 || error(
        "Material response registry already contains duplicate record IDs.",
    )
    if isempty(matches)
        push!(registry.records,record)
    elseif replace
        registry.records[first(matches)] = record
    else
        error("Material-property record already exists: $(record.record_id).")
    end
    identifiers = getfield.(registry.records,:record_id)
    length(unique(identifiers)) == length(identifiers) || error(
        "Material-property record identifiers must remain unique.",
    )
    return registry
end

function _atomistic_temperature_axis(table::Atomistic_Response_Table)
    "temperature_K" in table.axis_order || error(
        "Atomistic response table has no temperature_K axis.",
    )
    return findfirst(==("temperature_K"),table.axis_order),table.axes["temperature_K"]
end

function _atomistic_component_axis(table::Atomistic_Response_Table)
    candidates = ("component_index","tensor_component","component")
    matches = [name for name in candidates if name in table.axis_order]
    length(matches) == 1 || error(
        "A tensor response table must declare exactly one component axis.",
    )
    name = first(matches)
    return findfirst(==(name),table.axis_order),name,table.axes[name]
end

function _atomistic_temperature_series(
    table::Atomistic_Response_Table;
    component::Union{Nothing,Integer}=nothing,
)
    temperature_index,temperatures = _atomistic_temperature_axis(table)
    if length(table.axis_order) == 1
        temperature_index == 1 || error("Scalar response-table axis ordering is invalid.")
        return temperatures,vec(table.values),
            isnothing(table.standard_uncertainty) ? nothing :
                vec(table.standard_uncertainty)
    end
    length(table.axis_order) == 2 || error(
        "Cryogenic material binding supports one temperature axis and at most one component axis.",
    )
    component_index,_,component_axis = _atomistic_component_axis(table)
    isnothing(component) && error(
        "An anisotropic atomistic property requires an explicit component index.",
    )
    selected = Int(component)
    1 <= selected <= length(component_axis) || error(
        "Requested atomistic tensor component is out of range.",
    )
    values = if temperature_index == 1 && component_index == 2
        vec(table.values[:,selected])
    elseif component_index == 1 && temperature_index == 2
        vec(table.values[selected,:])
    else
        error("Response-table temperature/component axes are inconsistent.")
    end
    uncertainty = if isnothing(table.standard_uncertainty)
        nothing
    elseif temperature_index == 1
        vec(table.standard_uncertainty[:,selected])
    else
        vec(table.standard_uncertainty[selected,:])
    end
    return temperatures,values,uncertainty
end

"""
    cryogenic_property_from_atomistic_table(table; component=nothing, ...)

Convert a current-format, hash-bound heat-capacity or thermal-conductivity response table into the
property object consumed by the electrothermal solver. Tensor averaging is never implicit.
"""
function cryogenic_property_from_atomistic_table(
    table::Atomistic_Response_Table;
    component::Union{Nothing,Integer}=nothing,
    property_hash::AbstractString=table.table_hash,
    qualification_status::Symbol=table.status,
)
    table.quantity in (:heat_capacity,:thermal_conductivity_tensor) || error(
        "Cryogenic binding requires a heat-capacity or thermal-conductivity response table.",
    )
    temperatures,values,_ = _atomistic_temperature_series(table;component=component)
    all(value -> isfinite(value) && value >= 0.0,values) || error(
        "Atomistic thermal response contains invalid values.",
    )
    isempty(property_hash) && error("Cryogenic property hash cannot be empty.")
    return Tabulated_Cryogenic_Property(
        temperatures,values;
        units=table.units,
        property_hash=property_hash,
        qualification_status=qualification_status,
        allow_zero=any(values .== 0.0),
    )
end

function _atomistic_record_status(status::Symbol)
    status == :qualified && return :qualified
    status in (:verification,:candidate) && return :atomistic_candidate
    error("Unknown atomistic response-table status: $(status).")
end

function _relative_uncertainty_bound(
    values::Vector{Float64},
    uncertainty::Union{Nothing,Vector{Float64}},
)
    isnothing(uncertainty) && return 0.0
    bounds = Float64[]
    for index in eachindex(values)
        scale = abs(values[index])
        if scale > 0.0
            push!(bounds,uncertainty[index]/scale)
        elseif uncertainty[index] > 0.0
            return Inf
        end
    end
    return isempty(bounds) ? 0.0 : maximum(bounds)
end

"""Convert one current-format atomistic response table into a source-backed registry record."""
function material_record_from_atomistic_table(
    table::Atomistic_Response_Table;
    record_id::AbstractString=table.table_id,
    property::Symbol,
    direction::Symbol=:isotropic,
    component::Union{Nothing,Integer}=nothing,
    state_description::AbstractString=table.material_state_hash,
    source_title::AbstractString=string(
        table.calculation.code," ",table.calculation.method," response calculation",
    ),
    source_identifier::AbstractString=string(
        table.calculation.calculation_id,":",table.table_hash,
    ),
    source_url::AbstractString=get(
        table.metadata,"source_url",string("urn:sha256:",table.table_hash),
    ),
    notes::AbstractString="Generated from a hash-bound DFT/MD response table.",
)
    temperatures,values,uncertainty = _atomistic_temperature_series(
        table;component=component,
    )
    model = Tabulated_Material_Property(
        temperatures,values;
        standard_uncertainties=isnothing(uncertainty) ? zeros(length(values)) : uncertainty,
        interpolation=:linear,
        allow_zero=any(values .== 0.0),
    )
    relative_uncertainty = _relative_uncertainty_bound(values,uncertainty)
    isfinite(relative_uncertainty) || error(
        "A zero response with nonzero uncertainty cannot be represented by one relative uncertainty.",
    )
    metadata = merge(copy(table.metadata),Dict(
        "response_table_id" => table.table_id,
        "response_table_hash" => table.table_hash,
        "material_state_hash" => table.material_state_hash,
        "calculation_id" => table.calculation.calculation_id,
        "calculation_input_hash" => table.calculation.input_hash,
        "calculation_structure_hash" => table.calculation.structure_hash,
        "component" => isnothing(component) ? "scalar" : string(component),
        "notes" => String(notes),
    ))
    return Material_Property_Record(
        record_id=record_id,
        material_tag=table.material_tag,
        property=property,
        direction=direction,
        state_description=state_description,
        units=table.units,
        temperature_min_K=first(temperatures),
        temperature_max_K=last(temperatures),
        nominal_relative_uncertainty=relative_uncertainty,
        model=model,
        source_title=source_title,
        source_identifier=source_identifier,
        source_url=source_url,
        status=_atomistic_record_status(table.status),
        metadata=metadata,
    )
end

"""
Compatibility wrapper around the current `subkev_kernel_from_response_table` implementation.
The result remains tied to the response-table hash, material state, and requested temperature.
"""
function subkev_kernel_from_atomistic_table(
    table::Atomistic_Response_Table,
    temperature_K::Real;
    model_id::AbstractString=table.table_id,
    qualification_status::Symbol=table.status,
)
    binding = subkev_kernel_from_response_table(
        table,temperature_K;qualification_status=qualification_status,
    )
    kernel = binding.kernel
    model_id == kernel.model_id && return kernel
    return SubkeV_Thermalization_Kernel(
        kernel.material_tag,kernel.temperature_K,kernel.energy_grid_eV,kernel.fractions;
        model_id=model_id,
        model_hash=kernel.model_hash,
        material_state_hash=kernel.material_state_hash,
        qualification_status=kernel.qualification_status,
        metadata=merge(copy(kernel.metadata),Dict(
            "compatibility_wrapper" => "subkev_kernel_from_atomistic_table",
            "source_response_table_hash" => table.table_hash,
        )),
    )
end

"""Return source-backed acquisition and generation routes for missing HTS material responses."""
function hts_material_response_source_map()
    return Dict{String,Any}(
        "Cu" => Dict(
            "measured_properties" => ["specific_heat","thermal_conductivity"],
            "source" => "NIST Cryogenic Materials database",
            "implementation" => "default_material_response_registry",
        ),
        "Kapton" => Dict(
            "measured_properties" => ["specific_heat","thermal_conductivity"],
            "source" => "NIST Cryogenic Materials database",
            "implementation" => "default_material_response_registry",
        ),
        "G10" => Dict(
            "measured_properties" => ["specific_heat","normal/warp thermal conductivity"],
            "source" => "NIST Cryogenic Materials database",
            "implementation" => "default_material_response_registry",
        ),
        "Hastelloy-C276" => Dict(
            "measured_properties" => ["specific_heat","thermal_conductivity"],
            "source" => "NIST report / OSTI 5410732, 4.2--300 K",
            "next_step" => "digitize tables with uncertainty and bind exact alloy heat treatment",
        ),
        "YBCO" => Dict(
            "measured_properties" => [
                "coated-conductor specific heat to 2 K",
                "coated-conductor thermal diffusivity/specific heat",
                "single-crystal thermal-conductivity anisotropy",
                "ultrafast carrier/phonon relaxation",
            ],
            "generated_properties" => [
                "phonon heat capacity from Phonopy",
                "lattice thermal conductivity from phono3py or Green--Kubo MD",
                "electron-phonon spectral function and rates from EPW",
                "electronic stopping from real-time TDDFT",
                "defect storage and displacement thresholds from cascade MD",
            ],
        ),
        "GdBCO" => Dict(
            "measured_properties" => [
                "compound/tape heat capacity and thermal conductivity where available",
                "isotope-resolved capture rates and prompt cascade data",
            ],
            "generated_properties" => [
                "phonon and thermal tables using GdBCO-specific DFT/MLIP",
                "sub-keV partition from GdBCO dielectric/track-structure calculations",
            ],
            "prohibited_shortcut" =>
                "Do not reuse YBCO values without an explicit surrogate uncertainty.",
        ),
    )
end

export Tabulated_Material_Property,material_property_standard_uncertainty
export register_material_property!,cryogenic_property_from_atomistic_table
export material_record_from_atomistic_table,subkev_kernel_from_atomistic_table
export hts_material_response_source_map
