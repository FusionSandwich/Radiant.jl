const PHYSICAL_REFERENCE_BUNDLE_SCHEMA = "radiant.hts.physical_reference_bundle/v1"
const PHYSICAL_REFERENCE_CLASSIFICATIONS = (
    :physical,
    :continuous_energy_monte_carlo,
    :deterministic_reference,
    :analytic,
    :software_fixture,
)

"""One array-valued protected response from an independent calculation or experiment."""
struct Physical_Response_Field
    response_id::String
    values::Array{Float64}
    standard_uncertainty::Union{Nothing,Array{Float64}}
    units::String
    classification::Symbol
    scoring_semantics::String
    metadata::Dict{String,String}

    function Physical_Response_Field(
        response_id::AbstractString,
        values::AbstractArray{<:Real};
        standard_uncertainty::Union{Nothing,AbstractArray{<:Real}}=nothing,
        units::AbstractString,
        classification::Symbol,
        scoring_semantics::AbstractString,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(response_id) && error("Physical response identifier cannot be empty.")
        value_array = Array{Float64}(values)
        isempty(value_array) && error("Physical response values cannot be empty.")
        all(isfinite,value_array) || error("Physical response values must be finite.")
        uncertainty = isnothing(standard_uncertainty) ? nothing :
            Array{Float64}(standard_uncertainty)
        if !isnothing(uncertainty)
            size(uncertainty) == size(value_array) || error(
                "Physical response uncertainty must match the response shape.",
            )
            all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
                "Physical response standard uncertainties must be nonnegative.",
            )
        end
        isempty(units) && error("Physical response units cannot be empty.")
        classification in PHYSICAL_REFERENCE_CLASSIFICATIONS || error(
            "Unknown physical response classification $(classification).",
        )
        isempty(scoring_semantics) && error("Physical response scoring semantics cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(response_id),value_array,uncertainty,String(units),classification,
            String(scoring_semantics),metadata_string,
        )
    end
end

"""Hash-lineage container for matched OpenMC, OpenSn, Geant4, or measured responses."""
struct Physical_Reference_Bundle
    schema::String
    case_id::String
    producer_code::String
    producer_version::String
    source_hash::String
    geometry_hash::String
    material_state_hash::String
    model_hash::String
    responses::Dict{String,Physical_Response_Field}
    metadata::Dict{String,String}

    function Physical_Reference_Bundle(
        case_id::AbstractString,
        producer_code::AbstractString,
        producer_version::AbstractString,
        responses::AbstractVector{Physical_Response_Field};
        source_hash::AbstractString,
        geometry_hash::AbstractString,
        material_state_hash::AbstractString,
        model_hash::AbstractString,
        metadata::AbstractDict=Dict{String,String}(),
    )
        for (label,value) in (
            ("case",case_id),("producer code",producer_code),
            ("producer version",producer_version),("source hash",source_hash),
            ("geometry hash",geometry_hash),("material-state hash",material_state_hash),
            ("model hash",model_hash),
        )
            isempty(value) && error("Physical reference $(label) cannot be empty.")
        end
        dictionary = Dict{String,Physical_Response_Field}()
        for response in responses
            haskey(dictionary,response.response_id) && error(
                "Duplicate physical response identifier $(response.response_id).",
            )
            dictionary[response.response_id] = response
        end
        isempty(dictionary) && error("Physical reference bundle cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            PHYSICAL_REFERENCE_BUNDLE_SCHEMA,String(case_id),String(producer_code),
            String(producer_version),String(source_hash),String(geometry_hash),
            String(material_state_hash),String(model_hash),dictionary,metadata_string,
        )
    end
end

function _reference_c_flat(array::AbstractArray)
    ndims(array) == 1 && return collect(array)
    permutation = Tuple(reverse(collect(1:ndims(array))))
    return vec(permutedims(Array(array),permutation))
end

function _reference_restore_c(flat::AbstractVector,shape::AbstractVector{<:Integer})
    dimensions = Int64.(shape)
    isempty(dimensions) && error("Reference response shape cannot be empty.")
    prod(dimensions) == length(flat) || error(
        "Reference response data length does not match its declared shape.",
    )
    length(dimensions) == 1 && return reshape(Float64.(flat),dimensions[1])
    reverse_shape = Tuple(reverse(dimensions))
    permutation = Tuple(reverse(collect(1:length(dimensions))))
    return permutedims(reshape(Float64.(flat),reverse_shape),permutation)
end

function _write_string_dict_group(parent,name::AbstractString,values::Dict{String,String})
    group = HDF5.create_group(parent,name)
    keys_sorted = sort(collect(keys(values)))
    group["count"] = Int64(length(keys_sorted))
    for (index,key) in enumerate(keys_sorted)
        entry = HDF5.create_group(group,@sprintf("entry_%06d",index))
        entry["key"] = key
        entry["value"] = values[key]
    end
    return group
end

function _read_string_dict_group(parent,name::AbstractString)
    haskey(parent,name) || return Dict{String,String}()
    group = parent[name]
    count = Int(read(group["count"]))
    output = Dict{String,String}()
    for index in 1:count
        entry = group[@sprintf("entry_%06d",index)]
        key = string(read(entry["key"]))
        value = string(read(entry["value"]))
        haskey(output,key) && error("Duplicate reference metadata key $(key).")
        output[key] = value
    end
    return output
end

"""Write a cross-language, row-major flattened reference bundle and return its SHA-256."""
function write_physical_reference_bundle_hdf5(
    path::AbstractString,
    bundle::Physical_Reference_Bundle;
    overwrite::Bool=false,
)
    ispath(path) && !overwrite && error("Refusing to overwrite reference bundle $(path).")
    ispath(path) && rm(path;force=true)
    mkpath(dirname(abspath(path)))
    HDF5.h5open(path,"w") do file
        meta = HDF5.create_group(file,"meta")
        meta["schema"] = bundle.schema
        meta["storage_order"] = "canonical-c-row-major-flat/v1"
        meta["case_id"] = bundle.case_id
        meta["producer_code"] = bundle.producer_code
        meta["producer_version"] = bundle.producer_version
        meta["source_hash"] = bundle.source_hash
        meta["geometry_hash"] = bundle.geometry_hash
        meta["material_state_hash"] = bundle.material_state_hash
        meta["model_hash"] = bundle.model_hash
        _write_string_dict_group(meta,"metadata",bundle.metadata)

        responses = HDF5.create_group(file,"responses")
        response_ids = sort(collect(keys(bundle.responses)))
        responses["count"] = Int64(length(response_ids))
        for (index,response_id) in enumerate(response_ids)
            response = bundle.responses[response_id]
            group = HDF5.create_group(responses,@sprintf("response_%06d",index))
            group["response_id"] = response.response_id
            group["units"] = response.units
            group["classification"] = string(response.classification)
            group["scoring_semantics"] = response.scoring_semantics
            group["values_flat"] = _reference_c_flat(response.values)
            group["values_shape"] = Int64[size(response.values)...]
            present = !isnothing(response.standard_uncertainty)
            group["uncertainty_present"] = Int8(present)
            if present
                group["uncertainty_flat"] = _reference_c_flat(response.standard_uncertainty)
                group["uncertainty_shape"] = Int64[size(response.standard_uncertainty)...]
            end
            _write_string_dict_group(group,"metadata",response.metadata)
        end
    end
    return _source_file_sha256(path)
end

function _read_reference_string(parent,name::AbstractString)
    haskey(parent,name) || error("Reference bundle is missing $(name).")
    return string(read(parent[name]))
end

"""Read and hash-verify a physical reference bundle."""
function read_physical_reference_bundle_hdf5(
    path::AbstractString;
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
)
    artifact_hash = _verify_file_hash(path,expected_file_sha256)
    return HDF5.h5open(path,"r") do file
        meta = file["meta"]
        _read_reference_string(meta,"schema") == PHYSICAL_REFERENCE_BUNDLE_SCHEMA || error(
            "Unsupported physical reference bundle schema.",
        )
        _read_reference_string(meta,"storage_order") ==
            "canonical-c-row-major-flat/v1" || error(
            "Unsupported physical reference storage order.",
        )
        response_group = file["responses"]
        count = Int(read(response_group["count"]))
        responses = Physical_Response_Field[]
        for index in 1:count
            group = response_group[@sprintf("response_%06d",index)]
            values = _reference_restore_c(vec(read(group["values_flat"])),
                                          vec(Int64.(read(group["values_shape"]))))
            uncertainty = Bool(read(group["uncertainty_present"])) ?
                _reference_restore_c(vec(read(group["uncertainty_flat"])),
                                     vec(Int64.(read(group["uncertainty_shape"])))) : nothing
            push!(responses,Physical_Response_Field(
                _read_reference_string(group,"response_id"),values;
                standard_uncertainty=uncertainty,
                units=_read_reference_string(group,"units"),
                classification=Symbol(_read_reference_string(group,"classification")),
                scoring_semantics=_read_reference_string(group,"scoring_semantics"),
                metadata=_read_string_dict_group(group,"metadata"),
            ))
        end
        metadata = _read_string_dict_group(meta,"metadata")
        metadata["artifact_sha256"] = artifact_hash
        metadata["artifact_path"] = abspath(path)
        Physical_Reference_Bundle(
            _read_reference_string(meta,"case_id"),
            _read_reference_string(meta,"producer_code"),
            _read_reference_string(meta,"producer_version"),responses;
            source_hash=_read_reference_string(meta,"source_hash"),
            geometry_hash=_read_reference_string(meta,"geometry_hash"),
            material_state_hash=_read_reference_string(meta,"material_state_hash"),
            model_hash=_read_reference_string(meta,"model_hash"),metadata=metadata,
        )
    end
end

function get_physical_response_field(
    bundle::Physical_Reference_Bundle,
    response_id::AbstractString,
)
    key = String(response_id)
    haskey(bundle.responses,key) || error("Physical response $(key) is unavailable.")
    return bundle.responses[key]
end

function reduce_physical_response(
    field::Physical_Response_Field;
    reduction::Symbol=:sum,
    index::Union{Nothing,CartesianIndex}=nothing,
)
    if reduction == :sum
        value = sum(field.values)
        uncertainty = isnothing(field.standard_uncertainty) ? 0.0 :
            sqrt(sum(abs2,field.standard_uncertainty))
    elseif reduction == :mean
        value = sum(field.values)/length(field.values)
        uncertainty = isnothing(field.standard_uncertainty) ? 0.0 :
            sqrt(sum(abs2,field.standard_uncertainty))/length(field.values)
    elseif reduction == :maximum
        location = argmax(field.values)
        value = field.values[location]
        uncertainty = isnothing(field.standard_uncertainty) ? 0.0 :
            field.standard_uncertainty[location]
    elseif reduction == :index
        isnothing(index) && error("Index reduction requires a CartesianIndex.")
        value = field.values[index]
        uncertainty = isnothing(field.standard_uncertainty) ? 0.0 :
            field.standard_uncertainty[index]
    else
        error("Unknown physical response reduction $(reduction).")
    end
    return Physical_Reference_Response(
        field.response_id,value;
        standard_uncertainty=uncertainty,units=field.units,
        classification=field.classification,
        source_artifact_hash=get(field.metadata,"source_artifact_hash","unbound"),
        geometry_hash=get(field.metadata,"geometry_hash","unbound"),
        material_state_hash=get(field.metadata,"material_state_hash","unbound"),
        model_hash=get(field.metadata,"model_hash","unbound"),
        scoring_semantics=field.scoring_semantics,
        metadata=copy(field.metadata),
    )
end

function reference_bundle_is_physical(bundle::Physical_Reference_Bundle)
    return all(field.classification in (:physical,:continuous_energy_monte_carlo,
                                        :deterministic_reference)
               for field in values(bundle.responses)) &&
           all(value != "unbound" for value in (
               bundle.source_hash,bundle.geometry_hash,bundle.material_state_hash,
               bundle.model_hash,
           ))
end

"""Compare each process-resolved score channel to its independently named reference field."""
function qualify_process_score_from_bundle(
    score::Process_Resolved_Score,
    bundle::Physical_Reference_Bundle;
    response_map::AbstractDict,
    relative_tolerances::AbstractDict,
    absolute_tolerances::AbstractDict=Dict{String,Float64}(),
    reductions::AbstractDict=Dict{String,Symbol}(),
    candidate_standard_uncertainties::AbstractDict=Dict{String,Float64}(),
)
    comparisons = Dict{String,Matched_Response_Comparison}()
    for process_key in get_process_keys(score)
        haskey(response_map,process_key) || error(
            "No physical-reference mapping exists for process $(process_key).",
        )
        response_id = string(response_map[process_key])
        reference_field = get_physical_response_field(bundle,response_id)
        reduction = Symbol(get(reductions,process_key,:sum))
        reference = reduce_physical_response(reference_field;reduction=reduction)
        candidate_field = get_process_score(score,process_key)
        candidate = reduction == :sum ? sum(candidate_field) :
                    reduction == :mean ? sum(candidate_field)/length(candidate_field) :
                    reduction == :maximum ? maximum(candidate_field) : error(
                        "Process score bundle qualification supports sum, mean, or maximum.",
                    )
        haskey(relative_tolerances,process_key) || error(
            "No relative tolerance exists for process $(process_key).",
        )
        comparisons[process_key] = compare_matched_response(
            candidate,reference;
            candidate_standard_uncertainty=Float64(get(
                candidate_standard_uncertainties,process_key,0.0,
            )),
            relative_tolerance=Float64(relative_tolerances[process_key]),
            absolute_tolerance=Float64(get(absolute_tolerances,process_key,0.0)),
        )
    end
    all_pass = all(comparison.pass for comparison in values(comparisons))
    physical_pass = all_pass && reference_bundle_is_physical(bundle)
    return Dict{String,Any}(
        "schema" => "radiant.hts.process_reference_bundle_qualification/v1",
        "all_comparisons_pass" => all_pass,
        "reference_bundle_physical" => reference_bundle_is_physical(bundle),
        "production_pass" => physical_pass,
        "source_hash" => bundle.source_hash,
        "geometry_hash" => bundle.geometry_hash,
        "material_state_hash" => bundle.material_state_hash,
        "model_hash" => bundle.model_hash,
        "comparisons" => comparisons,
    )
end

export PHYSICAL_REFERENCE_BUNDLE_SCHEMA,Physical_Response_Field
export Physical_Reference_Bundle,write_physical_reference_bundle_hdf5
export read_physical_reference_bundle_hdf5,get_physical_response_field
export reduce_physical_response,reference_bundle_is_physical
export qualify_process_score_from_bundle
