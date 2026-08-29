function _simple_csv_rows(path::AbstractString)
    isfile(path) || error("Response-table CSV does not exist: $(path).")
    lines = [strip(line) for line in readlines(path) if !isempty(strip(line)) &&
             !startswith(strip(line),'#')]
    isempty(lines) && error("Response-table CSV is empty.")
    header = strip.(split(lines[1],','))
    isempty(header) && error("Response-table CSV header is empty.")
    length(unique(header)) == length(header) || error("Response-table CSV headers must be unique.")
    rows = Vector{Dict{String,String}}()
    for (line_number,line) in enumerate(lines[2:end],2)
        values = strip.(split(line,',';keepempty=true))
        length(values) == length(header) || error(
            "CSV row $(line_number) has $(length(values)) fields; expected $(length(header)).",
        )
        push!(rows,Dict(header[index] => values[index] for index in eachindex(header)))
    end
    isempty(rows) && error("Response-table CSV contains no data rows.")
    return header,rows
end

function _csv_float(row::Dict{String,String},key::String)
    haskey(row,key) || error("Response-table CSV is missing column $(key).")
    value = try
        parse(Float64,row[key])
    catch
        error("Response-table CSV value for $(key) is not numeric: $(row[key]).")
    end
    isfinite(value) || error("Response-table CSV value for $(key) must be finite.")
    return value
end

"""
    read_phonopy_thermal_properties_yaml(path; material_tag, ...)

Read the stable scalar fields written by Phonopy's `thermal_properties.yaml`. The adapter extracts
`temperature`, `heat_capacity`, `free_energy`, and `entropy` without requiring a YAML dependency.
Phonopy reports heat capacity and entropy in J/K/mol and free energy in kJ/mol. The resulting tables
remain molar until a separately source-backed molar mass and density conversion is applied.
"""
function read_phonopy_thermal_properties_yaml(
    path::AbstractString;
    material_tag::AbstractString,
    model_id::AbstractString="phonopy-harmonic",
    material_state_hash::AbstractString,
    manifest_hash::AbstractString,
    table_hash_prefix::AbstractString="phonopy",
    qualification_status::Symbol=:candidate,
    metadata::AbstractDict=Dict{String,String}(),
)
    isfile(path) || error("Phonopy thermal-properties file does not exist: $(path).")
    entries = Dict{String,Float64}[]
    current = Dict{String,Float64}()
    active = false
    for raw_line in readlines(path)
        line = strip(raw_line)
        startswith(line,"thermal_properties:") && (active = true; continue)
        active || continue
        startswith(line,"-") && begin
            !isempty(current) && push!(entries,current)
            current = Dict{String,Float64}()
            line = strip(line[2:end])
        end
        isempty(line) && continue
        occursin(':',line) || continue
        key,value_string = strip.(split(line,':';limit=2))
        key in ("temperature","heat_capacity","free_energy","entropy") || continue
        value = try
            parse(Float64,split(value_string)[1])
        catch
            error("Unable to parse Phonopy field $(key) from line: $(raw_line).")
        end
        current[key] = value
    end
    !isempty(current) && push!(entries,current)
    isempty(entries) && error("No Phonopy thermal-property entries were found.")
    required = ("temperature","heat_capacity","free_energy","entropy")
    all(all(haskey(entry,key) for key in required) for entry in entries) || error(
        "Every Phonopy thermal-property entry must contain temperature, heat capacity, free energy, and entropy.",
    )
    sort!(entries;by=entry -> entry["temperature"])
    temperatures = [entry["temperature"] for entry in entries]
    all(diff(temperatures) .> 0.0) || error("Phonopy temperatures must be strictly increasing.")
    metadata_string = Dict{String,String}(
        "source_format" => "phonopy-thermal_properties.yaml",
        "source_path" => abspath(path),
        "temperature_units" => "K",
        "heat_capacity_units" => "J/K/mol",
        "free_energy_units" => "kJ/mol",
        "entropy_units" => "J/K/mol",
        "molar_to_volumetric_conversion_applied" => "false",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    output = Dict{Symbol,Atomistic_Response_Table}()
    for (property_id,key,units) in (
        (:heat_capacity,"heat_capacity","J/K/mol"),
        (:free_energy,"free_energy","kJ/mol"),
        (:entropy,"entropy","J/K/mol"),
    )
        values = [entry[key] for entry in entries]
        output[property_id] = Atomistic_Response_Table(
            "$(table_hash_prefix)-$(property_id)",material_tag,property_id,
            ["temperature_K"],[temperatures],reshape(values,length(values));
            units=units,model_id=model_id,material_state_hash=material_state_hash,
            manifest_hash=manifest_hash,table_hash="$(table_hash_prefix)-$(property_id)",
            qualification_status=qualification_status,metadata=metadata_string,
        )
    end
    return output
end

function molar_heat_capacity_to_volumetric(
    table::Atomistic_Response_Table;
    molar_mass_kg_per_mol::Real,
    density_kg_per_m3::Real,
    conversion_hash::AbstractString,
)
    table.property_id == :heat_capacity || error("Input table is not heat capacity.")
    table.units == "J/K/mol" || error("Molar heat-capacity conversion expects J/K/mol.")
    molar_mass = Float64(molar_mass_kg_per_mol)
    density = Float64(density_kg_per_m3)
    isfinite(molar_mass) && molar_mass > 0.0 || error("Molar mass must be positive.")
    isfinite(density) && density > 0.0 || error("Density must be positive.")
    isempty(conversion_hash) && error("Molar conversion requires a provenance hash.")
    values = table.values .* (density/molar_mass)
    uncertainty = isnothing(table.standard_uncertainty) ? nothing :
        table.standard_uncertainty .* (density/molar_mass)
    metadata = merge(copy(table.metadata),Dict(
        "molar_to_volumetric_conversion_applied" => "true",
        "molar_mass_kg_per_mol" => string(molar_mass),
        "density_kg_per_m3" => string(density),
        "conversion_hash" => String(conversion_hash),
    ))
    return Atomistic_Response_Table(
        table.table_id*"-volumetric",table.material_tag,:volumetric_heat_capacity,
        table.axis_names,table.axes,values;
        standard_uncertainty=uncertainty,units="J/m3/K",model_id=table.model_id,
        material_state_hash=table.material_state_hash,manifest_hash=table.manifest_hash,
        table_hash=bytes2hex(SHA.sha256(codeunits(table.table_hash*"|"*conversion_hash))),
        qualification_status=table.qualification_status,metadata=metadata,
    )
end

"""Read replicated Green--Kubo conductivity results and return ensemble means and uncertainties."""
function read_green_kubo_conductivity_csv(
    path::AbstractString;
    material_tag::AbstractString,
    model_id::AbstractString,
    material_state_hash::AbstractString,
    manifest_hash::AbstractString,
    table_hash::AbstractString,
    qualification_status::Symbol=:candidate,
)
    header,rows = _simple_csv_rows(path)
    required = ["temperature_K","replica","kxx_W_m_K","kyy_W_m_K","kzz_W_m_K"]
    all(key in header for key in required) || error(
        "Green--Kubo CSV requires columns $(join(required,',' )).",
    )
    temperatures = sort(unique(_csv_float(row,"temperature_K") for row in rows))
    components = ("kxx_W_m_K","kyy_W_m_K","kzz_W_m_K")
    values = zeros(Float64,length(temperatures),3)
    uncertainties = zeros(Float64,length(temperatures),3)
    replica_counts = zeros(Int64,length(temperatures))
    for (temperature_index,temperature) in enumerate(temperatures)
        selected = [row for row in rows if _csv_float(row,"temperature_K") == temperature]
        replica_ids = [row["replica"] for row in selected]
        length(unique(replica_ids)) == length(replica_ids) || error(
            "Green--Kubo replica identifiers must be unique per temperature.",
        )
        replica_counts[temperature_index] = length(selected)
        for component_index in 1:3
            samples = [_csv_float(row,components[component_index]) for row in selected]
            all(samples .>= 0.0) || error("Thermal conductivity cannot be negative.")
            values[temperature_index,component_index] = sum(samples)/length(samples)
            if length(samples) > 1
                variance = sum((sample-values[temperature_index,component_index])^2
                               for sample in samples)/(length(samples)-1)
                uncertainties[temperature_index,component_index] =
                    sqrt(variance/length(samples))
            end
        end
    end
    return Atomistic_Response_Table(
        "green-kubo-kappa",material_tag,:thermal_conductivity,
        ["temperature_K","tensor_component"],
        [Float64.(temperatures),[1.0,2.0,3.0]],values;
        standard_uncertainty=uncertainties,units="W/m/K",model_id=model_id,
        material_state_hash=material_state_hash,manifest_hash=manifest_hash,
        table_hash=table_hash,qualification_status=qualification_status,
        metadata=Dict(
            "source_format" => "green-kubo-replicas-csv/v1",
            "tensor_component_names" => "xx,yy,zz",
            "replica_counts" => join(replica_counts,','),
            "uncertainty" => "standard-error-across-independent-replicas",
        ),
    )
end

"""Read a complete energy-temperature partition grid for a sub-keV response model."""
function read_subkev_partition_csv(
    path::AbstractString;
    material_tag::AbstractString,
    model_id::AbstractString,
    material_state_hash::AbstractString,
    manifest_hash::AbstractString,
    table_hash::AbstractString,
    qualification_status::Symbol=:candidate,
)
    header,rows = _simple_csv_rows(path)
    channel_names = String[string(channel) for channel in SUBKEV_PARTITION_CHANNELS]
    required = vcat(["energy_eV","temperature_K"],channel_names)
    all(key in header for key in required) || error(
        "Sub-keV partition CSV is missing required energy, temperature, or channel columns.",
    )
    energies = sort(unique(_csv_float(row,"energy_eV") for row in rows))
    temperatures = sort(unique(_csv_float(row,"temperature_K") for row in rows))
    expected_rows = length(energies)*length(temperatures)
    length(rows) == expected_rows || error(
        "Sub-keV partition CSV must contain one row for every energy-temperature pair.",
    )
    values = zeros(Float64,length(energies),length(temperatures),length(channel_names))
    occupied = falses(length(energies),length(temperatures))
    for row in rows
        energy_index = findfirst(==(_csv_float(row,"energy_eV")),energies)
        temperature_index = findfirst(==(_csv_float(row,"temperature_K")),temperatures)
        occupied[energy_index,temperature_index] && error(
            "Duplicate sub-keV energy-temperature row.",
        )
        occupied[energy_index,temperature_index] = true
        fractions = [_csv_float(row,channel) for channel in channel_names]
        all(fraction -> 0.0 <= fraction <= 1.0,fractions) || error(
            "Sub-keV partition fractions must lie in [0,1].",
        )
        isapprox(sum(fractions),1.0;rtol=1.0e-8,atol=1.0e-10) || error(
            "Sub-keV partition fractions do not close at energy $(energies[energy_index]) eV and temperature $(temperatures[temperature_index]) K.",
        )
        values[energy_index,temperature_index,:] .= fractions
    end
    all(occupied) || error("Sub-keV partition grid has missing points.")
    return Atomistic_Response_Table(
        "subkev-partition",material_tag,:subkev_partition,
        ["energy_eV","temperature_K","channel_index"],
        [Float64.(energies),Float64.(temperatures),collect(1.0:length(channel_names))],
        values;
        units="fraction",model_id=model_id,material_state_hash=material_state_hash,
        manifest_hash=manifest_hash,table_hash=table_hash,
        qualification_status=qualification_status,
        metadata=Dict(
            "source_format" => "subkev-partition-csv/v1",
            "channel_names" => join(channel_names,','),
            "energy_partition_closed" => "true",
        ),
    )
end

"""Read cascade-MD response replicas and return one table per requested response column."""
function read_cascade_md_response_csv(
    path::AbstractString;
    material_tag::AbstractString,
    response_columns::AbstractVector{<:AbstractString},
    model_id::AbstractString,
    material_state_hash::AbstractString,
    manifest_hash::AbstractString,
    table_hash_prefix::AbstractString,
    qualification_status::Symbol=:candidate,
)
    header,rows = _simple_csv_rows(path)
    required = ["recoil_energy_eV","orientation_id","replica"]
    all(key in header for key in required) || error(
        "Cascade-MD CSV requires recoil_energy_eV, orientation_id, and replica.",
    )
    all(String(column) in header for column in response_columns) || error(
        "Cascade-MD CSV is missing at least one requested response column.",
    )
    energies = sort(unique(_csv_float(row,"recoil_energy_eV") for row in rows))
    orientations = sort(unique(row["orientation_id"] for row in rows))
    output = Dict{String,Atomistic_Response_Table}()
    for response_column_value in response_columns
        response_column = String(response_column_value)
        values = zeros(Float64,length(energies),length(orientations))
        uncertainties = zeros(Float64,size(values))
        for (energy_index,energy) in enumerate(energies),
            (orientation_index,orientation) in enumerate(orientations)
            selected = [row for row in rows
                        if _csv_float(row,"recoil_energy_eV") == energy &&
                           row["orientation_id"] == orientation]
            isempty(selected) && error(
                "Cascade-MD grid lacks energy $(energy), orientation $(orientation).",
            )
            replicas = [row["replica"] for row in selected]
            length(unique(replicas)) == length(replicas) || error(
                "Cascade-MD replica identifiers must be unique per grid point.",
            )
            samples = [_csv_float(row,response_column) for row in selected]
            mean_value = sum(samples)/length(samples)
            values[energy_index,orientation_index] = mean_value
            if length(samples) > 1
                variance = sum((sample-mean_value)^2 for sample in samples)/
                           (length(samples)-1)
                uncertainties[energy_index,orientation_index] =
                    sqrt(variance/length(samples))
            end
        end
        output[response_column] = Atomistic_Response_Table(
            "cascade-md-$(response_column)",material_tag,Symbol(response_column),
            ["recoil_energy_eV","orientation_index"],
            [Float64.(energies),collect(1.0:length(orientations))],values;
            standard_uncertainty=uncertainties,units="declared-in-metadata",
            model_id=model_id,material_state_hash=material_state_hash,
            manifest_hash=manifest_hash,
            table_hash="$(table_hash_prefix)-$(response_column)",
            qualification_status=qualification_status,
            metadata=Dict(
                "source_format" => "cascade-md-replicas-csv/v1",
                "orientation_names" => join(orientations,','),
                "uncertainty" => "standard-error-across-independent-replicas",
            ),
        )
    end
    return output
end

function response_table_generation_commands()
    return Dict{String,Vector{String}}(
        "harmonic_heat_capacity" => [
            "phonopy --fc FORCE_CONSTANTS --mesh 30 30 15 --tprop --tmin 2 --tmax 300 --tstep 2",
            "Convert thermal_properties.yaml with read_phonopy_thermal_properties_yaml.",
        ],
        "anharmonic_conductivity" => [
            "phono3py --fc2 --fc3 --mesh 30 30 15 --br",
            "Retain kappa-*.hdf5, displacement sets, force outputs, and convergence receipts.",
        ],
        "green_kubo_conductivity" => [
            "Run independent equilibrated LAMMPS replicas per temperature and orientation.",
            "Export temperature_K,replica,kxx_W_m_K,kyy_W_m_K,kzz_W_m_K.",
        ],
        "electron_phonon_relaxation" => [
            "Run EPW on converged Wannier/electron/phonon meshes and export alpha2F or scattering-rate tables.",
        ],
        "electronic_stopping" => [
            "Run real-time TDDFT trajectories across projectile species, energies, directions, and impact parameters.",
            "Average stopping after equilibration and retain trajectory-level uncertainty.",
        ],
        "subkev_partition" => [
            "Combine dielectric-loss/track-structure event sampling with TDDFT carrier excitation, phonon, escape, and defect channels.",
            "Export a closed energy_eV x temperature_K x channel fraction table.",
        ],
        "cascade_damage" => [
            "Run MLIP/DFT-qualified cascade MD across recoil energy, species, crystal direction, temperature, and independent velocity seeds.",
            "Export surviving defects, antisites, clusters, stored energy, strain, and interface crossing.",
        ],
    )
end

export read_phonopy_thermal_properties_yaml,molar_heat_capacity_to_volumetric
export read_green_kubo_conductivity_csv,read_subkev_partition_csv
export read_cascade_md_response_csv,response_table_generation_commands
