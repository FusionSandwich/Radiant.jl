const PHYSICAL_REFERENCE_CLASSIFICATIONS = (
    :analytic,
    :synthetic,
    :continuous_energy_openmc,
    :geant4,
    :opensn,
    :experiment,
)

struct Physical_Reference_Response
    response_id::String
    process_key::String
    particle_tag::String
    units::String
    values::Array{Float64}
    standard_uncertainty::Array{Float64}
    classification::Symbol
    producer::String
    source_artifact_hash::String
    result_artifact_hash::String
    geometry_hash::String
    material_state_hash::String
    normalization_hash::String
    metadata::Dict{String,String}

    function Physical_Reference_Response(;
        response_id::AbstractString,
        process_key::AbstractString,
        particle_tag::AbstractString,
        units::AbstractString,
        values::AbstractArray{<:Real},
        standard_uncertainty::AbstractArray{<:Real},
        classification::Symbol,
        producer::AbstractString,
        source_artifact_hash::AbstractString,
        result_artifact_hash::AbstractString,
        geometry_hash::AbstractString,
        material_state_hash::AbstractString,
        normalization_hash::AbstractString,
        metadata::AbstractDict=Dict{String,String}(),
    )
        for (label,value) in (
            ("response ID",response_id),("process key",process_key),
            ("particle tag",particle_tag),("units",units),("producer",producer),
            ("source artifact hash",source_artifact_hash),
            ("result artifact hash",result_artifact_hash),("geometry hash",geometry_hash),
            ("material-state hash",material_state_hash),
            ("normalization hash",normalization_hash),
        )
            isempty(value) && error("Physical reference $(label) cannot be empty.")
        end
        classification in PHYSICAL_REFERENCE_CLASSIFICATIONS || error(
            "Unknown physical-reference classification $(classification).",
        )
        data = Float64.(values)
        uncertainty = Float64.(standard_uncertainty)
        size(data) == size(uncertainty) || error(
            "Physical-reference value and uncertainty arrays must match.",
        )
        all(isfinite,data) || error("Physical-reference values must be finite.")
        all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
            "Physical-reference uncertainties must be finite and nonnegative.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(response_id),String(process_key),String(particle_tag),String(units),data,
            uncertainty,classification,String(producer),String(source_artifact_hash),
            String(result_artifact_hash),String(geometry_hash),String(material_state_hash),
            String(normalization_hash),metadata_string,
        )
    end
end

struct Matched_Response_Comparison
    response_id::String
    process_key::String
    candidate_values::Array{Float64}
    reference_values::Array{Float64}
    candidate_standard_uncertainty::Array{Float64}
    reference_standard_uncertainty::Array{Float64}
    absolute_difference::Array{Float64}
    allowed_difference::Array{Float64}
    normalized_residual::Array{Float64}
    pass_mask::BitArray
    maximum_absolute_difference::Float64
    maximum_relative_difference::Float64
    maximum_normalized_residual::Float64
    passed::Bool
    metadata::Dict{String,String}
end

function compare_matched_response(
    candidate_values::AbstractArray{<:Real},
    candidate_standard_uncertainty::AbstractArray{<:Real},
    reference::Physical_Reference_Response;
    response_id::AbstractString=reference.response_id,
    process_key::AbstractString=reference.process_key,
    relative_tolerance::Real,
    absolute_tolerance::Real,
    uncertainty_multiplier::Real=2.0,
    allow_shape_flattening::Bool=false,
)
    candidate = Float64.(candidate_values)
    candidate_uncertainty = Float64.(candidate_standard_uncertainty)
    if allow_shape_flattening && length(candidate) == length(reference.values)
        candidate = reshape(candidate,size(reference.values))
        candidate_uncertainty = reshape(candidate_uncertainty,size(reference.values))
    end
    size(candidate) == size(reference.values) || error(
        "Candidate and physical-reference response shapes do not match.",
    )
    size(candidate_uncertainty) == size(candidate) || error(
        "Candidate response uncertainty shape does not match candidate values.",
    )
    all(isfinite,candidate) || error("Candidate response values must be finite.")
    all(value -> isfinite(value) && value >= 0.0,candidate_uncertainty) || error(
        "Candidate response uncertainty must be finite and nonnegative.",
    )
    rtol = Float64(relative_tolerance)
    atol = Float64(absolute_tolerance)
    multiplier = Float64(uncertainty_multiplier)
    all(value -> isfinite(value) && value >= 0.0,(rtol,atol,multiplier)) || error(
        "Response-comparison tolerances must be finite and nonnegative.",
    )
    absolute_difference = abs.(candidate-reference.values)
    combined_uncertainty = sqrt.(candidate_uncertainty.^2+reference.standard_uncertainty.^2)
    allowed = atol .+ rtol.*abs.(reference.values) .+ multiplier.*combined_uncertainty
    pass_mask = absolute_difference .<= allowed
    scale = max.(combined_uncertainty,atol,eps(Float64))
    normalized = absolute_difference./scale
    relative = absolute_difference./max.(abs.(reference.values),atol,eps(Float64))
    metadata = Dict(
        "relative_tolerance" => string(rtol),
        "absolute_tolerance" => string(atol),
        "uncertainty_multiplier" => string(multiplier),
        "reference_classification" => string(reference.classification),
        "reference_producer" => reference.producer,
        "source_artifact_hash" => reference.source_artifact_hash,
        "result_artifact_hash" => reference.result_artifact_hash,
        "geometry_hash" => reference.geometry_hash,
        "material_state_hash" => reference.material_state_hash,
        "normalization_hash" => reference.normalization_hash,
    )
    return Matched_Response_Comparison(
        String(response_id),String(process_key),candidate,copy(reference.values),
        candidate_uncertainty,copy(reference.standard_uncertainty),absolute_difference,allowed,
        normalized,pass_mask,maximum(absolute_difference),maximum(relative),maximum(normalized),
        all(pass_mask),metadata,
    )
end

function reference_is_physical(reference::Physical_Reference_Response)
    return reference.classification in (
        :continuous_energy_openmc,:geant4,:opensn,:experiment,
    )
end

struct Process_Score_Qualification
    schema::String
    score_particle::String
    score_quantity::String
    comparisons::Dict{String,Matched_Response_Comparison}
    required_process_keys::Vector{String}
    missing_process_keys::Vector{String}
    physical_reference_only::Bool
    passed::Bool
    source_hashes::Vector{String}
    result_hashes::Vector{String}
end

function qualify_process_resolved_score(
    score::Process_Resolved_Score,
    references::AbstractDict;
    candidate_standard_uncertainty::AbstractDict=Dict{String,Any}(),
    required_process_keys::AbstractVector{<:AbstractString}=get_process_keys(score),
    relative_tolerance::Real,
    absolute_tolerance::Real,
    uncertainty_multiplier::Real=2.0,
    physical_reference_only::Bool=true,
)
    required = String[string(value) for value in required_process_keys]
    length(unique(required)) == length(required) || error(
        "Required process keys must be unique.",
    )
    comparisons = Dict{String,Matched_Response_Comparison}()
    missing = String[]
    source_hashes = String[]
    result_hashes = String[]
    for process_key in required
        if !haskey(score.channels,process_key) || !haskey(references,process_key)
            push!(missing,process_key)
            continue
        end
        reference = references[process_key]
        reference isa Physical_Reference_Response || error(
            "Reference for $(process_key) has the wrong type.",
        )
        reference.process_key == process_key || error(
            "Reference process key does not match its dictionary key.",
        )
        reference.particle_tag == score.particle_tag || error(
            "Reference and process score particle tags do not match.",
        )
        reference.units == score.units || error(
            "Reference and process score units do not match.",
        )
        if physical_reference_only && !reference_is_physical(reference)
            error("Process qualification requires a physical reference for $(process_key).")
        end
        candidate_uncertainty = if haskey(candidate_standard_uncertainty,process_key)
            Float64.(candidate_standard_uncertainty[process_key])
        else
            zeros(Float64,size(score.channels[process_key]))
        end
        comparison = compare_matched_response(
            score.channels[process_key],candidate_uncertainty,reference;
            relative_tolerance=relative_tolerance,absolute_tolerance=absolute_tolerance,
            uncertainty_multiplier=uncertainty_multiplier,
        )
        comparisons[process_key] = comparison
        push!(source_hashes,reference.source_artifact_hash)
        push!(result_hashes,reference.result_artifact_hash)
    end
    passed = isempty(missing) && length(comparisons) == length(required) &&
             all(comparison.passed for comparison in values(comparisons))
    return Process_Score_Qualification(
        "radiant.hts.process_score_qualification/v1",score.particle_tag,score.quantity,
        comparisons,required,missing,physical_reference_only,passed,unique(source_hashes),
        unique(result_hashes),
    )
end

struct Curved_Reference_Qualification
    response_ids::Vector{String}
    comparisons::Dict{String,Matched_Response_Comparison}
    atlas_geometry_hash::String
    curved_geometry_hash::String
    source_artifact_hash::String
    reference_classification::Symbol
    source_mapping_closure_pass::Bool
    frame_mapping_closure_pass::Bool
    passed::Bool
end

function qualify_atlas_against_curved_reference(
    atlas_responses::AbstractDict,
    atlas_uncertainties::AbstractDict,
    curved_references::AbstractDict;
    atlas_geometry_hash::AbstractString,
    curved_geometry_hash::AbstractString,
    source_artifact_hash::AbstractString,
    relative_tolerance::Real,
    absolute_tolerance::Real,
    uncertainty_multiplier::Real=2.0,
    source_mapping_closure_pass::Bool,
    frame_mapping_closure_pass::Bool,
)
    for value in (atlas_geometry_hash,curved_geometry_hash,source_artifact_hash)
        isempty(value) && error("Curved-reference qualification hashes cannot be empty.")
    end
    response_ids = sort(String[string(key) for key in keys(curved_references)])
    isempty(response_ids) && error("Curved-reference qualification cannot be empty.")
    comparisons = Dict{String,Matched_Response_Comparison}()
    classifications = Symbol[]
    for response_id in response_ids
        haskey(atlas_responses,response_id) || error(
            "Atlas response $(response_id) is missing.",
        )
        haskey(atlas_uncertainties,response_id) || error(
            "Atlas uncertainty $(response_id) is missing.",
        )
        reference = curved_references[response_id]
        reference isa Physical_Reference_Response || error(
            "Curved reference $(response_id) has the wrong type.",
        )
        reference_is_physical(reference) || error(
            "Curved atlas qualification requires OpenSn, Geant4, OpenMC, or experiment.",
        )
        reference.geometry_hash == curved_geometry_hash || error(
            "Curved reference geometry hash mismatch.",
        )
        reference.source_artifact_hash == source_artifact_hash || error(
            "Curved reference source hash mismatch.",
        )
        comparisons[response_id] = compare_matched_response(
            atlas_responses[response_id],atlas_uncertainties[response_id],reference;
            response_id=response_id,process_key=reference.process_key,
            relative_tolerance=relative_tolerance,absolute_tolerance=absolute_tolerance,
            uncertainty_multiplier=uncertainty_multiplier,
        )
        push!(classifications,reference.classification)
    end
    length(unique(classifications)) == 1 || error(
        "Curved-reference responses must use one reference classification.",
    )
    passed = source_mapping_closure_pass && frame_mapping_closure_pass &&
             all(comparison.passed for comparison in values(comparisons))
    return Curved_Reference_Qualification(
        response_ids,comparisons,String(atlas_geometry_hash),String(curved_geometry_hash),
        String(source_artifact_hash),first(classifications),source_mapping_closure_pass,
        frame_mapping_closure_pass,passed,
    )
end

function physical_qualification_receipt(qualification::Process_Score_Qualification)
    return Dict{String,Any}(
        "schema" => qualification.schema,
        "score_particle" => qualification.score_particle,
        "score_quantity" => qualification.score_quantity,
        "required_process_keys" => qualification.required_process_keys,
        "missing_process_keys" => qualification.missing_process_keys,
        "physical_reference_only" => qualification.physical_reference_only,
        "passed" => qualification.passed,
        "source_hashes" => qualification.source_hashes,
        "result_hashes" => qualification.result_hashes,
        "comparisons" => Dict(
            key => Dict(
                "passed" => comparison.passed,
                "maximum_absolute_difference" => comparison.maximum_absolute_difference,
                "maximum_relative_difference" => comparison.maximum_relative_difference,
                "maximum_normalized_residual" => comparison.maximum_normalized_residual,
            ) for (key,comparison) in qualification.comparisons
        ),
    )
end

function curved_qualification_receipt(qualification::Curved_Reference_Qualification)
    return Dict{String,Any}(
        "schema" => "radiant.hts.curved_atlas_qualification/v1",
        "atlas_geometry_hash" => qualification.atlas_geometry_hash,
        "curved_geometry_hash" => qualification.curved_geometry_hash,
        "source_artifact_hash" => qualification.source_artifact_hash,
        "reference_classification" => string(qualification.reference_classification),
        "source_mapping_closure_pass" => qualification.source_mapping_closure_pass,
        "frame_mapping_closure_pass" => qualification.frame_mapping_closure_pass,
        "flat_curved_response_convergence_pass" => qualification.passed,
        "passed" => qualification.passed,
        "comparisons" => Dict(
            key => Dict(
                "passed" => comparison.passed,
                "maximum_relative_difference" => comparison.maximum_relative_difference,
                "maximum_normalized_residual" => comparison.maximum_normalized_residual,
            ) for (key,comparison) in qualification.comparisons
        ),
    )
end
