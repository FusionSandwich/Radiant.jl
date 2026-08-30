const PHYSICAL_REFERENCE_BUNDLE_SCHEMA = "radiant.hts.physical_reference_bundle/v2"
const PHYSICAL_REFERENCE_STORAGE_ORDER = "canonical-c-row-major-flat/v1"

"""Hash-bound collection of current-format physical reference responses."""
struct Physical_Reference_Bundle
    schema::String
    responses::Dict{String,Physical_Reference_Response}
    artifact_sha256::Union{Nothing,String}
    metadata::Dict{String,String}

    function Physical_Reference_Bundle(
        responses_input::AbstractVector{Physical_Reference_Response};
        artifact_sha256::Union{Nothing,AbstractString}=nothing,
        metadata::AbstractDict=Dict{String,String}(),
    )
        responses = Dict{String,Physical_Reference_Response}()
        for response in responses_input
            haskey(responses,response.response_id) && error(
                "Physical reference response IDs must be unique: $(response.response_id).",
            )
            responses[response.response_id] = response
        end
        isempty(responses) && error("Physical reference bundle cannot be empty.")
        digest = isnothing(artifact_sha256) ? nothing : lowercase(String(artifact_sha256))
        if !isnothing(digest)
            occursin(r"^[0-9a-f]{64}$",digest) || error(
                "Physical reference bundle SHA-256 must contain 64 hexadecimal characters.",
            )
        end
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(PHYSICAL_REFERENCE_BUNDLE_SCHEMA,responses,digest,metadata_string)
    end
end

function _physical_reference_source_sha256(path::AbstractString)
    isfile(path) || error("Physical reference artifact does not exist: $(path).")
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _write_physical_response(group,response::Physical_Reference_Response)
    group["response_id"] = response.response_id
    group["process_key"] = response.process_key
    group["particle_tag"] = response.particle_tag
    group["units"] = response.units
    group["classification"] = string(response.classification)
    group["producer"] = response.producer
    group["source_artifact_hash"] = response.source_artifact_hash
    group["result_artifact_hash"] = response.result_artifact_hash
    group["geometry_hash"] = response.geometry_hash
    group["material_state_hash"] = response.material_state_hash
    group["normalization_hash"] = response.normalization_hash
    group["shape"] = Int64[size(response.values)...]
    group["values_flat"] = _atomistic_row_major_flat(response.values)
    group["standard_uncertainty_flat"] =
        _atomistic_row_major_flat(response.standard_uncertainty)
    _write_atomistic_dictionary(group,"metadata",response.metadata)
    return nothing
end

"""
    write_physical_reference_bundle_hdf5(path, responses; ...)

Write OpenMC, OpenSn, experiment, or other independently produced responses with explicit units,
uncertainty, scoring/process identity, and source/geometry/material/normalization lineage. This
format stores reference results only; it never promotes a classification or fabricates uncertainty.
"""
function write_physical_reference_bundle_hdf5(
    path::AbstractString,
    responses::AbstractVector{Physical_Reference_Response};
    metadata::AbstractDict=Dict{String,String}(),
    overwrite::Bool=false,
)
    ispath(path) && !overwrite && error(
        "Refusing to overwrite physical reference artifact: $(path).",
    )
    ispath(path) && rm(path;force=true)
    mkpath(dirname(abspath(path)))
    bundle = Physical_Reference_Bundle(responses;metadata=metadata)
    HDF5.h5open(path,"w") do file
        meta = HDF5.create_group(file,"meta")
        response_group = HDF5.create_group(file,"responses")
        meta["schema"] = PHYSICAL_REFERENCE_BUNDLE_SCHEMA
        meta["storage_order"] = PHYSICAL_REFERENCE_STORAGE_ORDER
        meta["response_count"] = Int64(length(bundle.responses))
        _write_atomistic_dictionary(meta,"metadata",bundle.metadata)
        ordered = sort(
            collect(Base.values(bundle.responses));by=response -> response.response_id,
        )
        for (index,response) in enumerate(ordered)
            entry = HDF5.create_group(response_group,@sprintf("entry_%06d",index))
            _write_physical_response(entry,response)
        end
    end
    return _physical_reference_source_sha256(path)
end

function _read_physical_string(group,key::AbstractString)
    haskey(group,key) || error("Physical reference field is missing: $(key).")
    value = read(group[key])
    if value isa AbstractArray
        length(value) == 1 || error("Physical reference field $(key) must be scalar.")
        value = first(value)
    end
    return string(value)
end

function _read_physical_response(group)
    shape = vec(Int64.(read(group["shape"])))
    isempty(shape) && error("Physical response shape cannot be empty.")
    values = _atomistic_restore_row_major(vec(read(group["values_flat"])),shape)
    uncertainty = _atomistic_restore_row_major(
        vec(read(group["standard_uncertainty_flat"])),shape,
    )
    return Physical_Reference_Response(
        response_id=_read_physical_string(group,"response_id"),
        process_key=_read_physical_string(group,"process_key"),
        particle_tag=_read_physical_string(group,"particle_tag"),
        units=_read_physical_string(group,"units"),
        values=values,
        standard_uncertainty=uncertainty,
        classification=Symbol(_read_physical_string(group,"classification")),
        producer=_read_physical_string(group,"producer"),
        source_artifact_hash=_read_physical_string(group,"source_artifact_hash"),
        result_artifact_hash=_read_physical_string(group,"result_artifact_hash"),
        geometry_hash=_read_physical_string(group,"geometry_hash"),
        material_state_hash=_read_physical_string(group,"material_state_hash"),
        normalization_hash=_read_physical_string(group,"normalization_hash"),
        metadata=_read_atomistic_dictionary(group,"metadata"),
    )
end

function read_physical_reference_bundle_hdf5(
    path::AbstractString;
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
)
    digest = _physical_reference_source_sha256(path)
    if !isnothing(expected_file_sha256)
        lowercase(digest) == lowercase(String(expected_file_sha256)) || error(
            "Physical reference bundle SHA-256 mismatch.",
        )
    end
    responses = Physical_Reference_Response[]
    metadata = Dict{String,String}()
    HDF5.h5open(path,"r") do file
        haskey(file,"meta") && haskey(file,"responses") || error(
            "Physical reference bundle is missing meta or responses group.",
        )
        meta = file["meta"]
        _read_physical_string(meta,"schema") == PHYSICAL_REFERENCE_BUNDLE_SCHEMA || error(
            "Unsupported physical reference bundle schema.",
        )
        _read_physical_string(meta,"storage_order") == PHYSICAL_REFERENCE_STORAGE_ORDER || error(
            "Unsupported physical reference tensor storage order.",
        )
        count = Int(read(meta["response_count"]))
        count >= 1 || error("Physical reference bundle response count must be positive.")
        metadata = _read_atomistic_dictionary(meta,"metadata")
        for index in 1:count
            key = @sprintf("entry_%06d",index)
            haskey(file["responses"],key) || error(
                "Physical reference bundle is missing response entry $(index).",
            )
            push!(responses,_read_physical_response(file["responses"][key]))
        end
    end
    metadata["artifact_sha256"] = digest
    metadata["artifact_path"] = abspath(path)
    return Physical_Reference_Bundle(
        responses;artifact_sha256=digest,metadata=metadata,
    )
end

function references_by_process_key(bundle::Physical_Reference_Bundle)
    output = Dict{String,Physical_Reference_Response}()
    for response in Base.values(bundle.responses)
        haskey(output,response.process_key) && error(
            "More than one physical reference uses process key $(response.process_key).",
        )
        output[response.process_key] = response
    end
    return output
end

references_by_response_id(bundle::Physical_Reference_Bundle) = copy(bundle.responses)

function qualify_process_score_from_reference_bundle(
    score::Process_Resolved_Score,
    bundle::Physical_Reference_Bundle;
    candidate_standard_uncertainty::AbstractDict=Dict{String,Any}(),
    required_process_keys::AbstractVector{<:AbstractString}=get_process_keys(score),
    relative_tolerance::Real,
    absolute_tolerance::Real,
    uncertainty_multiplier::Real=2.0,
    physical_reference_only::Bool=true,
)
    return qualify_process_resolved_score(
        score,references_by_process_key(bundle);
        candidate_standard_uncertainty=candidate_standard_uncertainty,
        required_process_keys=required_process_keys,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        uncertainty_multiplier=uncertainty_multiplier,
        physical_reference_only=physical_reference_only,
    )
end

function physical_reference_bundle_receipt(bundle::Physical_Reference_Bundle)
    all_responses = collect(Base.values(bundle.responses))
    return Dict{String,Any}(
        "schema" => bundle.schema,
        "storage_order" => PHYSICAL_REFERENCE_STORAGE_ORDER,
        "artifact_sha256" => isnothing(bundle.artifact_sha256) ?
            "in-memory-not-written" : bundle.artifact_sha256,
        "response_count" => length(bundle.responses),
        "response_ids" => sort(collect(keys(bundle.responses))),
        "classifications" => sort(unique(
            string(response.classification) for response in all_responses
        )),
        "producers" => sort(unique(response.producer for response in all_responses)),
        "source_artifact_hashes" => sort(unique(
            response.source_artifact_hash for response in all_responses
        )),
        "geometry_hashes" => sort(unique(
            response.geometry_hash for response in all_responses
        )),
        "physical_reference_only" => all(reference_is_physical,all_responses),
        "metadata" => copy(bundle.metadata),
    )
end

export PHYSICAL_REFERENCE_BUNDLE_SCHEMA,PHYSICAL_REFERENCE_STORAGE_ORDER
export Physical_Reference_Bundle
export write_physical_reference_bundle_hdf5,read_physical_reference_bundle_hdf5
export references_by_process_key,references_by_response_id
export qualify_process_score_from_reference_bundle,physical_reference_bundle_receipt
