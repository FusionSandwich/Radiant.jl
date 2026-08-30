const SUBKEV_RESPONSE_CHANNEL_NAMES = String[
    "ionization",
    "electronic_excitation",
    "prompt_phonon",
    "athermal_phonon",
    "quasiparticle",
    "defect_storage",
    "optical_emission",
    "electron_escape",
]

function _response_file_sha256(path::AbstractString)
    isfile(path) || error("Response-table artifact does not exist: $(path).")
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _sample_mean_standard_error(values::AbstractVector{<:Real})
    data = Float64.(values)
    length(data) >= 2 || error("At least two independent replicas are required.")
    all(isfinite,data) || error("Replica values must be finite.")
    mean_value = sum(data)/length(data)
    variance = sum((value-mean_value)^2 for value in data)/(length(data)-1)
    return mean_value,sqrt(max(0.0,variance)/length(data))
end

"""
    read_green_kubo_replicas_csv(path, calculation; ...)

Read independent Green–Kubo replica results with columns
`temperature_K,replica_id,kxx_W_mK,kyy_W_mK,kzz_W_mK`. Means and standard errors are calculated
independently at every temperature. The result is a current-format `Atomistic_Response_Table`, not
a legacy response-table object.
"""
function read_green_kubo_replicas_csv(
    path::AbstractString,
    calculation::Atomistic_Calculation_Manifest;
    table_id::AbstractString="green-kubo-thermal-conductivity",
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
    table_hash::Union{Nothing,AbstractString}=nothing,
    status::Symbol=:candidate,
    metadata::AbstractDict=Dict{String,String}(),
)
    calculation.method == :green_kubo_md || error(
        "Green–Kubo replica ingestion requires calculation.method=:green_kubo_md.",
    )
    actual_hash = _response_file_sha256(path)
    if !isnothing(expected_file_sha256)
        lowercase(actual_hash) == lowercase(String(expected_file_sha256)) || error(
            "Green–Kubo source SHA-256 mismatch.",
        )
    end
    header,rows = _simple_csv_rows(path)
    required = ["temperature_K","replica_id","kxx_W_mK","kyy_W_mK","kzz_W_mK"]
    _csv_required_columns(header,required)

    temperatures = sort(unique(_csv_float(row,"temperature_K") for row in rows))
    all(value -> value > 0.0,temperatures) || error(
        "Green–Kubo temperatures must be positive.",
    )
    values = zeros(Float64,length(temperatures),3)
    uncertainty = zeros(Float64,size(values))
    replica_counts = Int64[]
    seen = Set{Tuple{Float64,String}}()
    component_columns = ("kxx_W_mK","kyy_W_mK","kzz_W_mK")
    for (temperature_index,temperature) in enumerate(temperatures)
        temperature_rows = [
            row for row in rows if _csv_float(row,"temperature_K") == temperature
        ]
        replica_ids = String[]
        for row in temperature_rows
            replica_id = strip(row["replica_id"])
            isempty(replica_id) && error("Green–Kubo replica identifiers cannot be empty.")
            key = (temperature,replica_id)
            key in seen && error(
                "Duplicate Green–Kubo temperature/replica pair: $(temperature), $(replica_id).",
            )
            push!(seen,key)
            push!(replica_ids,replica_id)
        end
        length(unique(replica_ids)) >= 2 || error(
            "Every Green–Kubo temperature requires at least two independent replicas.",
        )
        push!(replica_counts,length(replica_ids))
        for (component_index,column) in enumerate(component_columns)
            component_values = [_csv_float(row,column) for row in temperature_rows]
            all(value -> value >= 0.0,component_values) || error(
                "Thermal-conductivity components must be nonnegative.",
            )
            values[temperature_index,component_index],
                uncertainty[temperature_index,component_index] =
                _sample_mean_standard_error(component_values)
        end
    end

    metadata_string = Dict{String,String}(
        "source_file_sha256" => actual_hash,
        "source_format" => "green-kubo-replicas-csv/v1",
        "component_convention" => "xx,yy,zz",
        "replica_count_per_temperature" => join(replica_counts,','),
        "replicas_declared_independent" => "true",
        "finite_size_and_model_discrepancy_included" => "false",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    bound_table_hash = isnothing(table_hash) ? actual_hash : String(table_hash)
    return Atomistic_Response_Table(
        table_id=table_id,
        material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:thermal_conductivity_tensor,
        axis_order=["temperature_K","component_index"],
        axes=Dict(
            "temperature_K" => temperatures,
            "component_index" => [1.0,2.0,3.0],
        ),
        values=values,
        standard_uncertainty=uncertainty,
        units="W/(m*K)",
        calculation=calculation,
        table_hash=bound_table_hash,
        status=status,
        metadata=metadata_string,
    )
end

"""
    read_subkev_partition_csv(path, calculation; ...)

Read a complete tensor product over `energy_eV` and `temperature_K`. The eight required fraction
columns are the names in `SUBKEV_RESPONSE_CHANNEL_NAMES`. Either all corresponding `_sigma`
columns are present or none are. Every row must close to unity without renormalization.
"""
function read_subkev_partition_csv(
    path::AbstractString,
    calculation::Atomistic_Calculation_Manifest;
    table_id::AbstractString="subkev-energy-partition",
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
    table_hash::Union{Nothing,AbstractString}=nothing,
    status::Symbol=:candidate,
    closure_tolerance::Real=1.0e-8,
    metadata::AbstractDict=Dict{String,String}(),
)
    actual_hash = _response_file_sha256(path)
    if !isnothing(expected_file_sha256)
        lowercase(actual_hash) == lowercase(String(expected_file_sha256)) || error(
            "Sub-keV response source SHA-256 mismatch.",
        )
    end
    header,rows = _simple_csv_rows(path)
    required = vcat(["energy_eV","temperature_K"],SUBKEV_RESPONSE_CHANNEL_NAMES)
    _csv_required_columns(header,required)
    sigma_columns = [string(channel,"_sigma") for channel in SUBKEV_RESPONSE_CHANNEL_NAMES]
    sigma_present = [column in header for column in sigma_columns]
    (all(sigma_present) || !any(sigma_present)) || error(
        "Sub-keV uncertainty columns must be either complete or absent.",
    )
    has_uncertainty = all(sigma_present)

    energies = sort(unique(_csv_float(row,"energy_eV") for row in rows))
    temperatures = sort(unique(_csv_float(row,"temperature_K") for row in rows))
    all(value -> value > 0.0,energies) || error("Sub-keV energies must be positive.")
    all(value -> value > 0.0,temperatures) || error(
        "Sub-keV temperatures must be positive.",
    )
    expected_row_count = length(energies)*length(temperatures)
    length(rows) == expected_row_count || error(
        "Sub-keV response table must contain the complete energy-temperature tensor product.",
    )
    energy_index = Dict(value => index for (index,value) in enumerate(energies))
    temperature_index = Dict(value => index for (index,value) in enumerate(temperatures))
    values = zeros(Float64,length(energies),length(temperatures),
                   length(SUBKEV_RESPONSE_CHANNEL_NAMES))
    uncertainty = has_uncertainty ? zeros(Float64,size(values)) : nothing
    seen = Set{Tuple{Float64,Float64}}()
    tolerance = Float64(closure_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 || error(
        "Sub-keV closure tolerance must be finite and nonnegative.",
    )

    for row in rows
        energy = _csv_float(row,"energy_eV")
        temperature = _csv_float(row,"temperature_K")
        key = (energy,temperature)
        key in seen && error("Duplicate sub-keV energy-temperature row: $(key).")
        push!(seen,key)
        ie = energy_index[energy]
        it = temperature_index[temperature]
        row_sum = 0.0
        for (channel_index,channel) in enumerate(SUBKEV_RESPONSE_CHANNEL_NAMES)
            fraction = _csv_float(row,channel)
            0.0 <= fraction <= 1.0 || error(
                "Sub-keV partition fractions must lie in [0,1].",
            )
            values[ie,it,channel_index] = fraction
            row_sum += fraction
            if has_uncertainty
                sigma = _csv_float(row,string(channel,"_sigma"))
                sigma >= 0.0 || error("Sub-keV standard uncertainties cannot be negative.")
                uncertainty[ie,it,channel_index] = sigma
            end
        end
        isapprox(row_sum,1.0;rtol=tolerance,atol=tolerance) || error(
            "Sub-keV energy partition does not close at E=$(energy) eV, T=$(temperature) K: $(row_sum).",
        )
    end

    metadata_string = Dict{String,String}(
        "source_file_sha256" => actual_hash,
        "source_format" => "subkev-partition-csv/v1",
        "channel_names" => join(SUBKEV_RESPONSE_CHANNEL_NAMES,','),
        "row_closure_required" => "true",
        "renormalization_used" => "false",
        "uncertainty_columns_present" => string(has_uncertainty),
        "physical_qualification" => "false",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    bound_table_hash = isnothing(table_hash) ? actual_hash : String(table_hash)
    return Atomistic_Response_Table(
        table_id=table_id,
        material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:energy_partition,
        axis_order=["energy_eV","temperature_K","channel_index"],
        axes=Dict(
            "energy_eV" => energies,
            "temperature_K" => temperatures,
            "channel_index" => collect(1.0:length(SUBKEV_RESPONSE_CHANNEL_NAMES)),
        ),
        values=values,
        standard_uncertainty=uncertainty,
        units="fraction",
        calculation=calculation,
        table_hash=bound_table_hash,
        status=status,
        metadata=metadata_string,
    )
end

function _temperature_bracket(temperatures::Vector{Float64},temperature::Float64)
    temperature < temperatures[1] && error(
        "Requested temperature lies below the response-table range.",
    )
    temperature > temperatures[end] && error(
        "Requested temperature lies above the response-table range.",
    )
    length(temperatures) == 1 && return 1,1,0.0
    upper = searchsortedfirst(temperatures,temperature)
    upper == 1 && return 1,1,0.0
    temperatures[upper] == temperature && return upper,upper,0.0
    lower = upper-1
    fraction = (temperature-temperatures[lower])/(temperatures[upper]-temperatures[lower])
    return lower,upper,fraction
end

struct SubkeV_Response_Binding
    kernel::SubkeV_Thermalization_Kernel
    channel_standard_uncertainty::Union{Nothing,Dict{Symbol,Vector{Float64}}}
    response_table_hash::String
    interpolation::String
end

"""Bind a current-format energy-partition response table to a fixed-temperature transport kernel."""
function subkev_kernel_from_response_table(
    table::Atomistic_Response_Table,
    temperature_K::Real;
    qualification_status::Symbol=table.status,
)
    table.quantity == :energy_partition || error(
        "Sub-keV kernel binding requires an :energy_partition response table.",
    )
    table.axis_order == ["energy_eV","temperature_K","channel_index"] || error(
        "Sub-keV response-table axis order is unsupported.",
    )
    table.units == "fraction" || error("Sub-keV response-table units must be fraction.")
    channels = split(get(table.metadata,"channel_names",""),',')
    channels == SUBKEV_RESPONSE_CHANNEL_NAMES || error(
        "Sub-keV response-table channel ordering does not match the kernel contract.",
    )
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error(
        "Sub-keV binding temperature must be finite and positive.",
    )
    temperatures = table.axes["temperature_K"]
    lower,upper,fraction = _temperature_bracket(temperatures,temperature)
    energy_grid = table.axes["energy_eV"]
    fractions = Dict{Symbol,Vector{Float64}}()
    uncertainty = isnothing(table.standard_uncertainty) ? nothing :
        Dict{Symbol,Vector{Float64}}()
    for (channel_index,channel) in enumerate(channels)
        values_lower = vec(table.values[:,lower,channel_index])
        values_upper = vec(table.values[:,upper,channel_index])
        values = (1.0-fraction).*values_lower+fraction.*values_upper
        fractions[Symbol(channel)] = values
        if !isnothing(uncertainty)
            sigma_lower = vec(table.standard_uncertainty[:,lower,channel_index])
            sigma_upper = vec(table.standard_uncertainty[:,upper,channel_index])
            uncertainty[Symbol(channel)] = sqrt.(
                ((1.0-fraction).*sigma_lower).^2+(fraction.*sigma_upper).^2,
            )
        end
    end
    kernel = SubkeV_Thermalization_Kernel(
        table.material_tag,temperature,energy_grid,fractions;
        model_id=table.table_id,
        model_hash=table.table_hash,
        material_state_hash=table.material_state_hash,
        qualification_status=qualification_status,
        metadata=merge(copy(table.metadata),Dict(
            "source_response_table_hash" => table.table_hash,
            "temperature_interpolation" => lower == upper ? "exact" : "linear",
            "temperature_bracket_K" => string(temperatures[lower],",",temperatures[upper]),
        )),
    )
    return SubkeV_Response_Binding(
        kernel,uncertainty,table.table_hash,lower == upper ? "exact" : "linear",
    )
end

"""
    read_phonopy_heat_capacity_yaml(path, calculation; ...)

Read `temperature` and `heat_capacity` from Phonopy's `thermal_properties.yaml`. Heat capacity is
left in J/(mol*K) unless both density and molar mass are supplied, in which case it is converted to
J/(cm^3*K). No uncertainty is invented.
"""
function read_phonopy_heat_capacity_yaml(
    path::AbstractString,
    calculation::Atomistic_Calculation_Manifest;
    table_id::AbstractString="phonopy-harmonic-heat-capacity",
    density_g_cm3::Union{Nothing,Real}=nothing,
    molar_mass_g_mol::Union{Nothing,Real}=nothing,
    standard_uncertainty::Union{Nothing,AbstractVector{<:Real}}=nothing,
    table_hash::Union{Nothing,AbstractString}=nothing,
    status::Symbol=:candidate,
    metadata::AbstractDict=Dict{String,String}(),
)
    isfile(path) || error("Phonopy thermal-properties file does not exist: $(path).")
    entries = Tuple{Float64,Float64}[]
    current_temperature = nothing
    for raw_line in eachline(path)
        line = strip(raw_line)
        if startswith(line,"- temperature:")
            current_temperature = parse(Float64,strip(split(line,':';limit=2)[2]))
        elseif startswith(line,"temperature:")
            current_temperature = parse(Float64,strip(split(line,':';limit=2)[2]))
        elseif startswith(line,"heat_capacity:")
            isnothing(current_temperature) && error(
                "Phonopy heat_capacity appeared before its temperature.",
            )
            value = parse(Float64,split(strip(split(line,':';limit=2)[2]))[1])
            push!(entries,(Float64(current_temperature),value))
            current_temperature = nothing
        end
    end
    isempty(entries) && error("No temperature/heat-capacity pairs were found in Phonopy YAML.")
    sort!(entries;by=first)
    temperatures = first.(entries)
    heat_capacity = last.(entries)
    all(diff(temperatures) .> 0.0) || error(
        "Phonopy temperatures must be strictly increasing and unique.",
    )
    all(value -> isfinite(value) && value >= 0.0,heat_capacity) || error(
        "Phonopy heat capacities must be finite and nonnegative.",
    )
    convert_to_volumetric = !isnothing(density_g_cm3) || !isnothing(molar_mass_g_mol)
    convert_to_volumetric && (isnothing(density_g_cm3) || isnothing(molar_mass_g_mol)) && error(
        "Density and molar mass must be supplied together for volumetric conversion.",
    )
    units = "J/(mol*K)"
    if convert_to_volumetric
        density = Float64(density_g_cm3)
        molar_mass = Float64(molar_mass_g_mol)
        all(value -> isfinite(value) && value > 0.0,(density,molar_mass)) || error(
            "Density and molar mass must be finite and positive.",
        )
        heat_capacity .*= density/molar_mass
        units = "J/(cm^3*K)"
    end
    uncertainty = isnothing(standard_uncertainty) ? nothing :
        reshape(Float64.(standard_uncertainty),length(temperatures))
    if !isnothing(uncertainty)
        length(uncertainty) == length(temperatures) || error(
            "Phonopy heat-capacity uncertainty length is inconsistent.",
        )
        all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
            "Heat-capacity uncertainties must be finite and nonnegative.",
        )
        if convert_to_volumetric
            uncertainty .*= Float64(density_g_cm3)/Float64(molar_mass_g_mol)
        end
    end
    source_hash = _response_file_sha256(path)
    metadata_string = Dict{String,String}(
        "source_file_sha256" => source_hash,
        "source_format" => "phonopy-thermal_properties-yaml",
        "harmonic_heat_capacity" => "true",
        "density_g_cm3" => isnothing(density_g_cm3) ? "not-applied" : string(density_g_cm3),
        "molar_mass_g_mol" => isnothing(molar_mass_g_mol) ? "not-applied" : string(molar_mass_g_mol),
        "anharmonic_and_electronic_heat_capacity_included" => "false",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Atomistic_Response_Table(
        table_id=table_id,
        material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:heat_capacity,
        axis_order=["temperature_K"],
        axes=Dict("temperature_K" => temperatures),
        values=heat_capacity,
        standard_uncertainty=uncertainty,
        units=units,
        calculation=calculation,
        table_hash=isnothing(table_hash) ? source_hash : String(table_hash),
        status=status,
        metadata=metadata_string,
    )
end

function cryogenic_heat_capacity_from_response_table(
    table::Atomistic_Response_Table,
    volume_cm3::Real,
)
    table.quantity == :heat_capacity || error(
        "Cryogenic heat-capacity binding requires a :heat_capacity table.",
    )
    table.units == "J/(cm^3*K)" || error(
        "Cryogenic node binding requires volumetric J/(cm^3*K) heat capacity.",
    )
    volume = Float64(volume_cm3)
    isfinite(volume) && volume > 0.0 || error("Thermal-node volume must be positive.")
    values = vec(table.values).*volume
    return Tabulated_Cryogenic_Property(
        table.axes["temperature_K"],values;
        units="J/K",
        property_hash=table.table_hash,
        qualification_status=table.status,
    )
end

function cryogenic_conductance_from_conductivity_table(
    table::Atomistic_Response_Table;
    component_index::Integer,
    area_cm2::Real,
    length_cm::Real,
)
    table.quantity == :thermal_conductivity_tensor || error(
        "Conductance binding requires a thermal-conductivity tensor table.",
    )
    table.units == "W/(m*K)" || error("Thermal conductivity must use W/(m*K).")
    1 <= component_index <= size(table.values,2) || error(
        "Thermal-conductivity component index is out of range.",
    )
    area = Float64(area_cm2)
    length_value = Float64(length_cm)
    all(value -> isfinite(value) && value > 0.0,(area,length_value)) || error(
        "Conductance geometry must be finite and positive.",
    )
    conversion = area*1.0e-2/length_value
    conductance = vec(table.values[:,component_index]).*conversion
    return Tabulated_Cryogenic_Property(
        table.axes["temperature_K"],conductance;
        units="W/K",
        property_hash=table.table_hash,
        qualification_status=table.status,
    )
end

export SUBKEV_RESPONSE_CHANNEL_NAMES
export read_green_kubo_replicas_csv,read_subkev_partition_csv
export SubkeV_Response_Binding,subkev_kernel_from_response_table
export read_phonopy_heat_capacity_yaml,cryogenic_heat_capacity_from_response_table
export cryogenic_conductance_from_conductivity_table
