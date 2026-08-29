const ATOMISTIC_RESPONSE_TABLE_STATUSES = (
    :verification,
    :candidate,
    :qualified,
)

const ATOMISTIC_RESPONSE_QUANTITIES = (
    :thermal_conductivity_tensor,
    :heat_capacity,
    :phonon_lifetime,
    :electron_phonon_relaxation,
    :electronic_stopping_power,
    :energy_partition,
    :displacement_threshold,
    :defect_stored_energy,
)

"""Provenance and convergence requirements for one DFT, DFPT, TDDFT, or MD calculation."""
struct Atomistic_Calculation_Manifest
    calculation_id::String
    material_tag::String
    material_state_hash::String
    method::Symbol
    code::String
    code_version::String
    input_hash::String
    structure_hash::String
    potential_or_functional::String
    convergence_parameters::Dict{String,String}
    references::Vector{String}
    status::Symbol
    metadata::Dict{String,String}

    function Atomistic_Calculation_Manifest(;
        calculation_id::AbstractString,
        material_tag::AbstractString,
        material_state_hash::AbstractString,
        method::Symbol,
        code::AbstractString,
        code_version::AbstractString,
        input_hash::AbstractString,
        structure_hash::AbstractString,
        potential_or_functional::AbstractString,
        convergence_parameters::AbstractDict,
        references::AbstractVector{<:AbstractString}=String[],
        status::Symbol=:candidate,
        metadata::AbstractDict=Dict{String,String}(),
    )
        for (label,value) in (
            ("calculation ID",calculation_id),("material tag",material_tag),
            ("material state hash",material_state_hash),("code",code),
            ("code version",code_version),("input hash",input_hash),
            ("structure hash",structure_hash),("potential or functional",potential_or_functional),
        )
            isempty(value) && error("Atomistic $(label) cannot be empty.")
        end
        method in (:phono3py,:green_kubo_md,:epw,:rt_tddft,:cascade_md,:dft_static) ||
            error("Unsupported atomistic response method $(method).")
        status in ATOMISTIC_RESPONSE_TABLE_STATUSES || error(
            "Unknown atomistic calculation status $(status).",
        )
        convergence = Dict{String,String}()
        for (key,value) in convergence_parameters
            convergence[string(key)] = string(value)
        end
        isempty(convergence) && error(
            "Atomistic calculation manifest requires explicit convergence parameters.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(calculation_id),String(material_tag),String(material_state_hash),method,
            String(code),String(code_version),String(input_hash),String(structure_hash),
            String(potential_or_functional),convergence,
            String[string(value) for value in references],status,metadata_string,
        )
    end
end

"""
    Atomistic_Response_Table

N-dimensional response table with named numeric axes. Values are stored in the declared
`axis_order`; the length of every axis must exactly match the corresponding array dimension.
Uncertainty is mandatory for `:qualified` tables.
"""
struct Atomistic_Response_Table
    table_id::String
    material_tag::String
    material_state_hash::String
    quantity::Symbol
    axis_order::Vector{String}
    axes::Dict{String,Vector{Float64}}
    values::Array{Float64}
    standard_uncertainty::Union{Nothing,Array{Float64}}
    units::String
    calculation::Atomistic_Calculation_Manifest
    table_hash::String
    status::Symbol
    metadata::Dict{String,String}

    function Atomistic_Response_Table(;
        table_id::AbstractString,
        material_tag::AbstractString,
        material_state_hash::AbstractString,
        quantity::Symbol,
        axis_order::AbstractVector{<:AbstractString},
        axes::AbstractDict,
        values::AbstractArray{<:Real},
        standard_uncertainty::Union{Nothing,AbstractArray{<:Real}}=nothing,
        units::AbstractString,
        calculation::Atomistic_Calculation_Manifest,
        table_hash::AbstractString="unbound",
        status::Symbol=:candidate,
        metadata::AbstractDict=Dict{String,String}(),
    )
        for (label,value) in (
            ("table ID",table_id),("material tag",material_tag),
            ("material state hash",material_state_hash),("units",units),
            ("table hash",table_hash),
        )
            isempty(value) && error("Atomistic response $(label) cannot be empty.")
        end
        quantity in ATOMISTIC_RESPONSE_QUANTITIES || error(
            "Unsupported atomistic response quantity $(quantity).",
        )
        status in ATOMISTIC_RESPONSE_TABLE_STATUSES || error(
            "Unknown atomistic response status $(status).",
        )
        order = String[string(value) for value in axis_order]
        isempty(order) && error("Atomistic response table requires at least one axis.")
        length(unique(order)) == length(order) || error(
            "Atomistic response axis names must be unique.",
        )
        axis_dictionary = Dict{String,Vector{Float64}}()
        for name in order
            haskey(axes,name) || error("Atomistic response axis $(name) is missing.")
            axis = Float64.(axes[name])
            isempty(axis) && error("Atomistic response axis $(name) cannot be empty.")
            all(isfinite,axis) || error("Atomistic response axes must be finite.")
            length(axis) == 1 || all(diff(axis) .> 0.0) || error(
                "Atomistic response axis $(name) must be strictly increasing.",
            )
            axis_dictionary[name] = axis
        end
        data = Float64.(values)
        ndims(data) == length(order) || error(
            "Atomistic response rank does not match axis_order.",
        )
        expected_shape = Tuple(length(axis_dictionary[name]) for name in order)
        size(data) == expected_shape || error(
            "Atomistic response shape " * string(size(data)) *
            " does not match axes " * string(expected_shape) * ".",
        )
        all(isfinite,data) || error("Atomistic response values must be finite.")
        uncertainty = isnothing(standard_uncertainty) ? nothing :
            Float64.(standard_uncertainty)
        if !isnothing(uncertainty)
            size(uncertainty) == size(data) || error(
                "Atomistic response uncertainty must match the value array.",
            )
            all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
                "Atomistic response uncertainty must be finite and nonnegative.",
            )
        end
        status == :qualified && isnothing(uncertainty) && error(
            "Qualified atomistic response tables require uncertainty.",
        )
        calculation.material_tag == material_tag || error(
            "Calculation and response-table material tags do not match.",
        )
        calculation.material_state_hash == material_state_hash || error(
            "Calculation and response-table material-state hashes do not match.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(table_id),String(material_tag),String(material_state_hash),quantity,order,
            axis_dictionary,data,uncertainty,String(units),calculation,String(table_hash),status,
            metadata_string,
        )
    end
end

function _atomistic_row_major_flat(array::AbstractArray)
    ndims(array) == 1 && return vec(Array(array))
    permutation = Tuple(reverse(collect(1:ndims(array))))
    return vec(permutedims(Array(array),permutation))
end

function _atomistic_restore_row_major(flat::AbstractVector,shape::AbstractVector{<:Integer})
    dimensions = Int64.(shape)
    isempty(dimensions) && error("Atomistic response shape cannot be empty.")
    prod(dimensions) == length(flat) || error(
        "Atomistic response flat data length and shape disagree.",
    )
    length(dimensions) == 1 && return reshape(Float64.(flat),dimensions[1])
    reverse_shape = Tuple(reverse(dimensions))
    permutation = Tuple(reverse(collect(1:length(dimensions))))
    return permutedims(reshape(Float64.(flat),reverse_shape),permutation)
end

function _write_atomistic_dictionary(parent,name::AbstractString,dictionary::AbstractDict)
    group = HDF5.create_group(parent,name)
    pairs = sort(collect(dictionary);by=entry -> string(first(entry)))
    group["count"] = Int64(length(pairs))
    for (index,(key,value)) in enumerate(pairs)
        entry = HDF5.create_group(group,@sprintf("entry_%06d",index))
        entry["key"] = string(key)
        entry["value"] = string(value)
    end
    return nothing
end

function _read_atomistic_dictionary(parent,name::AbstractString)
    haskey(parent,name) || return Dict{String,String}()
    group = parent[name]
    count = Int(read(group["count"]))
    output = Dict{String,String}()
    for index in 1:count
        entry = group[@sprintf("entry_%06d",index)]
        output[string(read(entry["key"]))] = string(read(entry["value"]))
    end
    return output
end

function write_atomistic_response_hdf5(
    path::AbstractString,
    table::Atomistic_Response_Table;
    overwrite::Bool=false,
)
    isfile(path) && !overwrite && error("Refusing to overwrite atomistic response $(path).")
    mkpath(dirname(abspath(path)))
    HDF5.h5open(path,"w") do file
        meta = HDF5.create_group(file,"meta")
        axes_group = HDF5.create_group(file,"axes")
        data_group = HDF5.create_group(file,"data")
        calculation = HDF5.create_group(file,"calculation")
        meta["schema"] = "radiant.hts.atomistic_response_table/v1"
        meta["table_id"] = table.table_id
        meta["material_tag"] = table.material_tag
        meta["material_state_hash"] = table.material_state_hash
        meta["quantity"] = string(table.quantity)
        meta["axis_order"] = table.axis_order
        meta["units"] = table.units
        meta["table_hash"] = table.table_hash
        meta["status"] = string(table.status)
        for name in table.axis_order
            axes_group[name] = table.axes[name]
        end
        data_group["values_flat"] = _atomistic_row_major_flat(table.values)
        data_group["shape"] = Int64[size(table.values)...]
        data_group["uncertainty_present"] = Int8(!isnothing(table.standard_uncertainty))
        if !isnothing(table.standard_uncertainty)
            data_group["standard_uncertainty_flat"] =
                _atomistic_row_major_flat(table.standard_uncertainty)
        end
        calculation["calculation_id"] = table.calculation.calculation_id
        calculation["method"] = string(table.calculation.method)
        calculation["code"] = table.calculation.code
        calculation["code_version"] = table.calculation.code_version
        calculation["input_hash"] = table.calculation.input_hash
        calculation["structure_hash"] = table.calculation.structure_hash
        calculation["potential_or_functional"] = table.calculation.potential_or_functional
        calculation["references"] = table.calculation.references
        calculation["status"] = string(table.calculation.status)
        _write_atomistic_dictionary(
            calculation,"convergence_parameters",table.calculation.convergence_parameters,
        )
        _write_atomistic_dictionary(calculation,"metadata",table.calculation.metadata)
        _write_atomistic_dictionary(meta,"metadata",table.metadata)
    end
    return _source_file_sha256(path)
end

function read_atomistic_response_hdf5(
    path::AbstractString;
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
)
    artifact_hash = _verify_file_hash(path,expected_file_sha256)
    return HDF5.h5open(path,"r") do file
        meta = file["meta"]
        axes_group = file["axes"]
        data_group = file["data"]
        calculation_group = file["calculation"]
        string(read(meta["schema"])) == "radiant.hts.atomistic_response_table/v1" || error(
            "Unsupported atomistic response table schema.",
        )
        axis_order = String[string(value) for value in vec(read(meta["axis_order"]))]
        axes = Dict(name => Float64.(vec(read(axes_group[name]))) for name in axis_order)
        shape = Int64.(vec(read(data_group["shape"])))
        values = _atomistic_restore_row_major(
            vec(read(data_group["values_flat"])),shape,
        )
        uncertainty_present = Bool(Int(read(data_group["uncertainty_present"])))
        uncertainty = uncertainty_present ? _atomistic_restore_row_major(
            vec(read(data_group["standard_uncertainty_flat"])),shape,
        ) : nothing
        material_tag = string(read(meta["material_tag"]))
        material_state_hash = string(read(meta["material_state_hash"]))
        calculation = Atomistic_Calculation_Manifest(
            calculation_id=string(read(calculation_group["calculation_id"])),
            material_tag=material_tag,
            material_state_hash=material_state_hash,
            method=Symbol(string(read(calculation_group["method"]))),
            code=string(read(calculation_group["code"])),
            code_version=string(read(calculation_group["code_version"])),
            input_hash=string(read(calculation_group["input_hash"])),
            structure_hash=string(read(calculation_group["structure_hash"])),
            potential_or_functional=string(read(calculation_group["potential_or_functional"])),
            convergence_parameters=_read_atomistic_dictionary(
                calculation_group,"convergence_parameters",
            ),
            references=String[string(value) for value in vec(read(calculation_group["references"]))],
            status=Symbol(string(read(calculation_group["status"]))),
            metadata=_read_atomistic_dictionary(calculation_group,"metadata"),
        )
        metadata = _read_atomistic_dictionary(meta,"metadata")
        metadata["artifact_sha256"] = artifact_hash
        return Atomistic_Response_Table(
            table_id=string(read(meta["table_id"])),material_tag=material_tag,
            material_state_hash=material_state_hash,
            quantity=Symbol(string(read(meta["quantity"]))),axis_order=axis_order,axes=axes,
            values=values,standard_uncertainty=uncertainty,units=string(read(meta["units"])),
            calculation=calculation,table_hash=string(read(meta["table_hash"])),
            status=Symbol(string(read(meta["status"]))),metadata=metadata,
        )
    end
end

function atomistic_table_is_production_ready(table::Atomistic_Response_Table)
    return table.status == :qualified && table.calculation.status == :qualified &&
           table.table_hash != "unbound" && table.material_state_hash != "unbound" &&
           !isnothing(table.standard_uncertainty)
end

"""
Read a phono3py `kappa-*.hdf5` tensor. Six Voigt components are retained as a component axis.
The file remains a candidate until force constants, q-mesh, isotope scattering, boundary model,
and convergence evidence are bound in `calculation`.
"""
function read_phono3py_kappa_hdf5(
    path::AbstractString,
    calculation::Atomistic_Calculation_Manifest;
    dataset::AbstractString="kappa",
    temperature_dataset::AbstractString="temperature",
    table_id::AbstractString="phono3py-kappa",
    units::AbstractString="W/(m*K)",
    table_hash::AbstractString="unbound",
    status::Symbol=:candidate,
)
    isfile(path) || error("phono3py HDF5 file does not exist: $(path).")
    temperatures,kappa = HDF5.h5open(path,"r") do file
        haskey(file,temperature_dataset) || error("phono3py temperature dataset is missing.")
        haskey(file,dataset) || error("phono3py conductivity dataset $(dataset) is missing.")
        Float64.(vec(read(file[temperature_dataset]))),Float64.(read(file[dataset]))
    end
    ndims(kappa) == 2 || error("phono3py conductivity dataset must have rank two.")
    if size(kappa,1) == length(temperatures)
        canonical = kappa
    elseif size(kappa,2) == length(temperatures)
        canonical = permutedims(kappa,(2,1))
    else
        error("phono3py conductivity dimensions do not match temperature count.")
    end
    component_count = size(canonical,2)
    component_count in (3,6,9) || error(
        "Expected 3, 6, or 9 thermal-conductivity components from phono3py.",
    )
    metadata = Dict(
        "source_file_sha256" => _source_file_sha256(path),
        "component_convention" => component_count == 6 ? "xx,yy,zz,yz,xz,xy" :
                                  "producer-declared",
        "physical_qualification" => "false",
    )
    return Atomistic_Response_Table(
        table_id=table_id,material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:thermal_conductivity_tensor,
        axis_order=["temperature_K","component_index"],
        axes=Dict(
            "temperature_K" => temperatures,
            "component_index" => collect(1.0:component_count),
        ),
        values=canonical,units=units,calculation=calculation,table_hash=table_hash,
        status=status,metadata=metadata,
    )
end

function read_epw_spectral_function(
    path::AbstractString,
    calculation::Atomistic_Calculation_Manifest;
    frequency_column::Integer=1,
    alpha2f_column::Integer=2,
    frequency_units::AbstractString="meV",
    table_id::AbstractString="epw-alpha2f",
    table_hash::AbstractString="unbound",
    status::Symbol=:candidate,
)
    isfile(path) || error("EPW spectral-function file does not exist: $(path).")
    frequency = Float64[]
    alpha2f = Float64[]
    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line,"#") && continue
        fields = split(replace(line,',' => ' '))
        maximum((frequency_column,alpha2f_column)) <= length(fields) || continue
        values = try
            parse.(Float64,fields)
        catch
            continue
        end
        push!(frequency,values[frequency_column])
        push!(alpha2f,values[alpha2f_column])
    end
    length(frequency) >= 2 || error("EPW spectral-function parser found too few rows.")
    permutation = sortperm(frequency)
    frequency = frequency[permutation]
    alpha2f = alpha2f[permutation]
    all(diff(frequency) .> 0.0) || error("EPW frequency axis must be strictly increasing.")
    return Atomistic_Response_Table(
        table_id=table_id,material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:electron_phonon_relaxation,axis_order=["frequency"],
        axes=Dict("frequency" => frequency),values=alpha2f,units="alpha2F($(frequency_units))",
        calculation=calculation,table_hash=table_hash,status=status,
        metadata=Dict(
            "source_file_sha256" => _source_file_sha256(path),
            "frequency_units" => String(frequency_units),
            "interpretation" => "spectral-input-not-direct-lifetime",
        ),
    )
end

function read_electronic_stopping_table(
    path::AbstractString,
    calculation::Atomistic_Calculation_Manifest;
    table_id::AbstractString="electronic-stopping",
    table_hash::AbstractString="unbound",
    status::Symbol=:candidate,
)
    isfile(path) || error("Electronic stopping table does not exist: $(path).")
    energy = Float64[]
    stopping = Float64[]
    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line,"#") && continue
        fields = split(line,',')
        length(fields) >= 2 || continue
        parsed = try
            (parse(Float64,strip(fields[1])),parse(Float64,strip(fields[2])))
        catch
            continue
        end
        push!(energy,parsed[1])
        push!(stopping,parsed[2])
    end
    length(energy) >= 2 || error("Stopping table parser found too few numerical rows.")
    permutation = sortperm(energy)
    energy = energy[permutation]
    stopping = stopping[permutation]
    all(diff(energy) .> 0.0) || error("Stopping-energy axis must be strictly increasing.")
    all(value -> isfinite(value) && value >= 0.0,stopping) || error(
        "Electronic stopping values must be finite and nonnegative.",
    )
    return Atomistic_Response_Table(
        table_id=table_id,material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:electronic_stopping_power,axis_order=["kinetic_energy_eV"],
        axes=Dict("kinetic_energy_eV" => energy),values=stopping,
        units="eV/Angstrom",calculation=calculation,table_hash=table_hash,status=status,
        metadata=Dict(
            "source_file_sha256" => _source_file_sha256(path),
            "required_columns" => "kinetic_energy_eV,stopping_eV_per_Angstrom",
        ),
    )
end

function default_ybco_gdbco_atomistic_plan()
    return Dict{String,Any}(
        "schema" => "radiant.hts.atomistic_response_plan/v1",
        "materials" => ["YBCO6+x","GdBCO6+x"],
        "state_axes" => [
            "oxygen content and ordered chain configuration",
            "temperature",
            "crystal orientation and texture",
            "isotope composition",
            "point-defect and cascade-damage state",
            "strain and magnetic field where relevant",
        ],
        "calculations" => [
            Dict(
                "quantity" => "thermal-conductivity tensor and phonon lifetime",
                "primary_method" => "DFT/DFPT third-order force constants plus phono3py",
                "cross_check" => "equilibrium Green-Kubo MD with a validated MLIP",
                "required_convergence" => [
                    "supercell","q mesh","third-order cutoff","displacement amplitude",
                    "isotope scattering","boundary mean free path","MLIP ensemble",
                ],
                "references" => [
                    "https://phonopy.github.io/phono3py/",
                    "https://docs.lammps.org/compute_heat_flux.html",
                ],
            ),
            Dict(
                "quantity" => "electron-phonon spectral function and relaxation",
                "primary_method" => "DFPT/Wannier interpolation with EPW",
                "cross_check" => "time-resolved optical-response literature and experiment",
                "required_convergence" => [
                    "electronic k mesh","phonon q mesh","Wannier window","smearing",
                    "Coulomb parameter","oxygen order","carrier concentration",
                ],
                "references" => ["https://docs.epw-code.org/"],
            ),
            Dict(
                "quantity" => "electronic stopping and sub-keV energy partition",
                "primary_method" => "real-time TDDFT trajectories across crystallographic directions",
                "cross_check" => "track-structure calculation and measured stopping where available",
                "required_convergence" => [
                    "time step","cell size","trajectory impact parameter","projectile charge state",
                    "velocity grid","exchange-correlation functional","pseudopotential",
                ],
                "references" => ["https://octopus-code.org/"],
            ),
            Dict(
                "quantity" => "defect production, stored energy, and threshold displacement",
                "primary_method" => "directional DFT/MLIP recoil calculations and cascade MD",
                "cross_check" => "independent-cascade and overlap benchmarks",
                "required_convergence" => [
                    "cell size","recoil direction","thermostat exclusion","time horizon",
                    "potential ensemble","electronic stopping coupling","defect classifier",
                ],
                "references" => [
                    "https://doi.org/10.1088/1361-6668/ac47dc",
                ],
            ),
        ],
        "qualification_rule" =>
            "No generated table is production-qualified until independent method, finite-size, " *
            "sampling, and material-state uncertainties are quantified and hash-bound.",
    )
end
