"""Source-backed tabulated material property with bounded interpolation."""
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
        length(data) == length(temperatures) == length(uncertainties) || error(
            "Tabulated material-property arrays have inconsistent lengths.",
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

function material_property_value(model::Tabulated_Material_Property,temperature_K::Real)
    temperature = Float64(temperature_K)
    temperature < model.temperatures_K[1] && error(
        "Temperature lies below the tabulated material-property domain.",
    )
    temperature > model.temperatures_K[end] && error(
        "Temperature lies above the tabulated material-property domain.",
    )
    temperature == model.temperatures_K[end] && return model.values[end]
    upper = searchsortedfirst(model.temperatures_K,temperature)
    upper == 1 && return model.values[1]
    model.temperatures_K[upper] == temperature && return model.values[upper]
    lower = upper-1
    if model.interpolation == :linear
        fraction = (temperature-model.temperatures_K[lower])/
                   (model.temperatures_K[upper]-model.temperatures_K[lower])
        return (1.0-fraction)*model.values[lower]+fraction*model.values[upper]
    end
    fraction = (log(temperature)-log(model.temperatures_K[lower]))/
               (log(model.temperatures_K[upper])-log(model.temperatures_K[lower]))
    return exp((1.0-fraction)*log(model.values[lower])+
               fraction*log(model.values[upper]))
end

function material_property_standard_uncertainty(
    model::Tabulated_Material_Property,
    temperature_K::Real,
)
    temperature = Float64(temperature_K)
    temperature < model.temperatures_K[1] && error(
        "Temperature lies below the tabulated uncertainty domain.",
    )
    temperature > model.temperatures_K[end] && error(
        "Temperature lies above the tabulated uncertainty domain.",
    )
    temperature == model.temperatures_K[end] && return model.standard_uncertainties[end]
    upper = searchsortedfirst(model.temperatures_K,temperature)
    upper == 1 && return model.standard_uncertainties[1]
    model.temperatures_K[upper] == temperature && return model.standard_uncertainties[upper]
    lower = upper-1
    fraction = (temperature-model.temperatures_K[lower])/
               (model.temperatures_K[upper]-model.temperatures_K[lower])
    return (1.0-fraction)*model.standard_uncertainties[lower]+
           fraction*model.standard_uncertainties[upper]
end

function register_material_property!(
    registry::Material_Response_Registry,
    record::Material_Property_Record;
    replace::Bool=false,
)
    key = (record.material_tag,record.property_id,record.orientation)
    haskey(registry.records,key) && !replace && error(
        "Material property already exists: $(key).",
    )
    registry.records[key] = record
    return registry
end

function _atomistic_temperature_axis(table::Atomistic_Response_Table)
    index = findfirst(==("temperature_K"),table.axis_names)
    isnothing(index) && error("Atomistic table has no temperature_K axis.")
    return index,table.axes[index]
end

"""
    cryogenic_property_from_atomistic_table(table; component=nothing, ...)

Convert a qualified or explicitly candidate atomistic response table into the interpolation object
used by the electrothermal solver. Scalar tables require only a temperature axis. Tensor tables
must declare a component axis and caller-selected component; implicit tensor averaging is rejected.
"""
function cryogenic_property_from_atomistic_table(
    table::Atomistic_Response_Table;
    component::Union{Nothing,Integer}=nothing,
    property_hash::AbstractString=table.table_hash,
    qualification_status::Symbol=table.qualification_status,
)
    temperature_index,temperatures = _atomistic_temperature_axis(table)
    values = if length(table.axis_names) == 1
        temperature_index == 1 || error("Scalar atomistic table axis ordering is invalid.")
        vec(table.values)
    elseif length(table.axis_names) == 2
        component_index = findfirst(name -> name in ("tensor_component","component"),
                                    table.axis_names)
        isnothing(component_index) && error(
            "Two-dimensional thermal table must declare a tensor_component axis.",
        )
        isnothing(component) && error(
            "An anisotropic atomistic property requires an explicit component index.",
        )
        component_value = Int(component)
        1 <= component_value <= length(table.axes[component_index]) || error(
            "Requested atomistic tensor component is out of range.",
        )
        if temperature_index == 1
            vec(table.values[:,component_value])
        else
            vec(table.values[component_value,:])
        end
    else
        error("Electrothermal conversion supports one- or two-axis thermal tables only.")
    end
    all(value -> isfinite(value) && value >= 0.0,values) || error(
        "Atomistic thermal table contains invalid property values.",
    )
    units = table.units
    occursin("J",units) || occursin("W",units) || error(
        "Atomistic thermal table units are not recognized as heat capacity or conductance.",
    )
    return Tabulated_Cryogenic_Property(
        temperatures,values;
        units=units,property_hash=property_hash,
        qualification_status=qualification_status,
        allow_zero=any(values .== 0.0),
    )
end

function material_record_from_atomistic_table(
    table::Atomistic_Response_Table;
    material_tag::AbstractString=table.material_tag,
    property_id::Symbol,
    orientation::Symbol=:isotropic,
    component::Union{Nothing,Integer}=nothing,
    source_id::AbstractString=table.manifest_hash,
    citation::AbstractString=get(table.metadata,"citation","atomistic-generated"),
    notes::AbstractString="Generated from a hash-bound DFT/MD response table.",
)
    temperature_index,temperatures = _atomistic_temperature_axis(table)
    values = if length(table.axis_names) == 1
        vec(table.values)
    else
        component_index = findfirst(name -> name in ("tensor_component","component"),
                                    table.axis_names)
        isnothing(component_index) && error(
            "Atomistic material record requires a declared component axis.",
        )
        isnothing(component) && error("Explicit tensor component is required.")
        temperature_index == 1 ? vec(table.values[:,Int(component)]) :
                                 vec(table.values[Int(component),:])
    end
    uncertainties = zeros(Float64,length(values))
    if !isnothing(table.standard_uncertainty)
        uncertainties = if length(table.axis_names) == 1
            vec(table.standard_uncertainty)
        else
            temperature_index == 1 ? vec(table.standard_uncertainty[:,Int(component)]) :
                                     vec(table.standard_uncertainty[Int(component),:])
        end
    end
    model = Tabulated_Material_Property(
        temperatures,values;
        standard_uncertainties=uncertainties,
        interpolation=:linear,
        allow_zero=any(values .== 0.0),
    )
    return Material_Property_Record(
        String(material_tag),property_id,model;
        orientation=orientation,units=table.units,source_id=source_id,
        citation=citation,data_hash=table.table_hash,
        qualification_status=table.qualification_status,notes=notes,
    )
end

function _subkev_channel_names(table::Atomistic_Response_Table)
    raw = get(table.metadata,"channel_names","")
    isempty(raw) && error("Sub-keV partition table metadata requires channel_names.")
    return Symbol.(strip.(split(raw,',')))
end

"""
    subkev_kernel_from_atomistic_table(table, temperature_K; ...)

Build a `SubkeV_Thermalization_Kernel` from an atomistic table with axes
`energy_eV`, `temperature_K`, and `channel_index`. The table values are channel fractions and must
close to one at every sampled energy/temperature point. No unlisted energy channel is sent to heat.
"""
function subkev_kernel_from_atomistic_table(
    table::Atomistic_Response_Table,
    temperature_K::Real;
    model_id::AbstractString=table.model_id,
    qualification_status::Symbol=table.qualification_status,
)
    Set(table.axis_names) == Set(["energy_eV","temperature_K","channel_index"]) || error(
        "Sub-keV table must have energy_eV, temperature_K, and channel_index axes.",
    )
    energy_index = findfirst(==("energy_eV"),table.axis_names)
    temperature_index = findfirst(==("temperature_K"),table.axis_names)
    channel_index = findfirst(==("channel_index"),table.axis_names)
    energies = table.axes[energy_index]
    temperatures = table.axes[temperature_index]
    channels = _subkev_channel_names(table)
    length(channels) == length(table.axes[channel_index]) || error(
        "Sub-keV channel_names count does not match channel_index axis.",
    )
    Set(channels) == Set(SUBKEV_PARTITION_CHANNELS) || error(
        "Sub-keV atomistic table channels do not match the required partition ledger.",
    )
    temperature = Float64(temperature_K)
    temperature < temperatures[1] && error("Temperature lies below the sub-keV table.")
    temperature > temperatures[end] && error("Temperature lies above the sub-keV table.")
    upper = searchsortedfirst(temperatures,temperature)
    lower = upper
    fraction = 0.0
    if upper > 1 && (upper > length(temperatures) || temperatures[upper] != temperature)
        lower = upper-1
        fraction = (temperature-temperatures[lower])/(temperatures[upper]-temperatures[lower])
    end

    fractions = Dict{Symbol,Vector{Float64}}()
    for (channel_position,channel) in enumerate(channels)
        values = zeros(Float64,length(energies))
        for energy_position in eachindex(energies)
            index_low = ntuple(dimension -> begin
                dimension == energy_index ? energy_position :
                dimension == temperature_index ? lower : channel_position
            end,3)
            index_high = ntuple(dimension -> begin
                dimension == energy_index ? energy_position :
                dimension == temperature_index ? upper : channel_position
            end,3)
            values[energy_position] = (1.0-fraction)*table.values[index_low...]+
                                      fraction*table.values[index_high...]
        end
        fractions[channel] = values
    end
    for energy_position in eachindex(energies)
        total = sum(fractions[channel][energy_position]
                    for channel in SUBKEV_PARTITION_CHANNELS)
        isapprox(total,1.0;rtol=1.0e-8,atol=1.0e-10) || error(
            "Atomistic sub-keV fractions fail closure at energy index $(energy_position).",
        )
    end
    return SubkeV_Thermalization_Kernel(
        table.material_tag,temperature,energies,fractions;
        model_id=model_id,model_hash=table.table_hash,
        material_state_hash=table.material_state_hash,
        qualification_status=qualification_status,
        metadata=merge(copy(table.metadata),Dict(
            "manifest_hash" => table.manifest_hash,
            "source_table_hash" => table.table_hash,
            "interpolated_temperature_K" => string(temperature),
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
            "prohibited_shortcut" => "Do not reuse YBCO values without an explicit surrogate uncertainty.",
        ),
    )
end

export Tabulated_Material_Property,material_property_standard_uncertainty
export register_material_property!,cryogenic_property_from_atomistic_table
export material_record_from_atomistic_table,subkev_kernel_from_atomistic_table
export hts_material_response_source_map
