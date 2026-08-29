"""One nested faceting level and its protected local-response results."""
struct Curvature_Refinement_Level
    level_id::String
    characteristic_facet_size_cm::Float64
    maximum_sagitta_cm::Float64
    responses::Dict{String,Float64}
    standard_uncertainties::Dict{String,Float64}
    hotspot_ids::Dict{String,String}
    geometry_hash::String
    source_hash::String
    metadata::Dict{String,String}

    function Curvature_Refinement_Level(
        level_id::AbstractString,
        characteristic_facet_size_cm::Real,
        maximum_sagitta_cm::Real,
        responses::AbstractDict;
        standard_uncertainties::AbstractDict=Dict{String,Float64}(),
        hotspot_ids::AbstractDict=Dict{String,String}(),
        geometry_hash::AbstractString,
        source_hash::AbstractString,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(level_id) && error("Curvature refinement level identifier cannot be empty.")
        h = Float64(characteristic_facet_size_cm)
        sagitta = Float64(maximum_sagitta_cm)
        isfinite(h) && h > 0.0 || error("Characteristic facet size must be positive.")
        isfinite(sagitta) && sagitta >= 0.0 || error("Maximum sagitta must be nonnegative.")
        isempty(geometry_hash) && error("Curvature level geometry hash cannot be empty.")
        isempty(source_hash) && error("Curvature level source hash cannot be empty.")
        response_values = Dict{String,Float64}()
        for (key,value) in responses
            numeric = Float64(value)
            isfinite(numeric) || error("Curvature response values must be finite.")
            response_values[string(key)] = numeric
        end
        isempty(response_values) && error("Curvature refinement level requires responses.")
        uncertainty_values = Dict{String,Float64}()
        for key in keys(response_values)
            numeric = Float64(get(standard_uncertainties,key,0.0))
            isfinite(numeric) && numeric >= 0.0 || error(
                "Curvature response uncertainties must be nonnegative.",
            )
            uncertainty_values[key] = numeric
        end
        hotspots = Dict{String,String}()
        for (key,value) in hotspot_ids
            haskey(response_values,string(key)) || error(
                "Hotspot identity was supplied for an unknown response.",
            )
            hotspots[string(key)] = string(value)
        end
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(level_id),h,sagitta,response_values,uncertainty_values,hotspots,
            String(geometry_hash),String(source_hash),metadata_string,
        )
    end
end

struct Curvature_Response_Estimate
    response_id::String
    coarse_to_medium_change::Float64
    medium_to_fine_change::Float64
    observed_order::Union{Nothing,Float64}
    richardson_limit::Union{Nothing,Float64}
    fine_grid_gci::Union{Nothing,Float64}
    combined_standard_uncertainty::Float64
    hotspot_stable::Bool
    curved_reference_comparison::Union{Nothing,Matched_Response_Comparison}
    pass::Bool
    notes::String
end

struct Curvature_Convergence_Result
    schema::String
    levels::Vector{Curvature_Refinement_Level}
    estimates::Dict{String,Curvature_Response_Estimate}
    common_source_hash::String
    all_responses_pass::Bool
    physical_curved_reference_present::Bool
    production_pass::Bool
    metadata::Dict{String,String}
end

function _symmetric_relative_change(first::Float64,second::Float64;atol::Float64=1.0e-30)
    return abs(second-first)/max(abs(first),abs(second),atol)
end

function _observed_order_and_limit(
    values::NTuple{3,Float64},
    sizes::NTuple{3,Float64},
)
    y1,y2,y3 = values
    h1,h2,h3 = sizes
    d12 = y1-y2
    d23 = y2-y3
    if d12 == 0.0 || d23 == 0.0 || sign(d12) != sign(d23)
        return nothing,nothing,nothing
    end
    r12 = h1/h2
    r23 = h2/h3
    r12 > 1.0 && r23 > 1.0 || return nothing,nothing,nothing
    if !isapprox(r12,r23;rtol=0.05,atol=0.0)
        return nothing,nothing,nothing
    end
    ratio = abs(d12/d23)
    ratio > 0.0 || return nothing,nothing,nothing
    p = log(ratio)/log(sqrt(r12*r23))
    isfinite(p) && p > 0.0 || return nothing,nothing,nothing
    r = r23
    denominator = r^p-1.0
    denominator > 0.0 || return nothing,nothing,nothing
    limit = y3+(y3-y2)/denominator
    gci = 1.25*abs(y3-y2)/denominator/max(abs(y3),eps(Float64))
    return p,limit,gci
end

"""
    qualify_curvature_convergence(levels; ...)

Estimate faceting-response uncertainty from at least three nested levels. Levels are sorted from
coarse to fine by characteristic facet size. The same source hash is required at every level.
A Richardson/GCI estimate is reported only for monotone, approximately uniform refinement.

A software convergence pass is not promoted to `production_pass` unless every protected response
also has a physical curved-reference comparison that passes its declared tolerances.
"""
function qualify_curvature_convergence(
    levels::AbstractVector{Curvature_Refinement_Level};
    relative_tolerances::AbstractDict,
    absolute_tolerances::AbstractDict=Dict{String,Float64}(),
    curved_references::AbstractDict=Dict{String,Physical_Reference_Response}(),
    reference_uncertainty_multiplier::Real=2.0,
    require_hotspot_stability::Bool=true,
    metadata::AbstractDict=Dict{String,String}(),
)
    level_vector = sort(Curvature_Refinement_Level[levels...];
                        by=level -> level.characteristic_facet_size_cm,rev=true)
    length(level_vector) >= 3 || error(
        "Curvature convergence requires at least three nested facet levels.",
    )
    source_hashes = unique(getfield.(level_vector,:source_hash))
    length(source_hashes) == 1 || error(
        "Curvature levels must use the same hash-bound source.",
    )
    for index in 2:length(level_vector)
        level_vector[index].characteristic_facet_size_cm <
            level_vector[index-1].characteristic_facet_size_cm || error(
            "Curvature facet sizes must strictly decrease after sorting.",
        )
    end
    response_sets = [Set(keys(level.responses)) for level in level_vector]
    all(response_sets[index] == response_sets[1] for index in 2:length(response_sets)) || error(
        "Every curvature level must contain the same response identifiers.",
    )
    response_ids = sort(collect(response_sets[1]))
    for response_id in response_ids
        haskey(relative_tolerances,response_id) || error(
            "No relative curvature tolerance was supplied for $(response_id).",
        )
    end

    estimates = Dict{String,Curvature_Response_Estimate}()
    fine = level_vector[end]
    medium = level_vector[end-1]
    coarse = level_vector[end-2]
    for response_id in response_ids
        coarse_change = _symmetric_relative_change(
            coarse.responses[response_id],medium.responses[response_id],
        )
        fine_change = _symmetric_relative_change(
            medium.responses[response_id],fine.responses[response_id],
        )
        p,limit,gci = _observed_order_and_limit(
            (coarse.responses[response_id],medium.responses[response_id],
             fine.responses[response_id]),
            (coarse.characteristic_facet_size_cm,medium.characteristic_facet_size_cm,
             fine.characteristic_facet_size_cm),
        )
        uncertainty = sqrt(
            medium.standard_uncertainties[response_id]^2+
            fine.standard_uncertainties[response_id]^2,
        )
        hotspot_values = String[]
        for level in level_vector
            haskey(level.hotspot_ids,response_id) && push!(hotspot_values,
                                                           level.hotspot_ids[response_id])
        end
        hotspot_stable = isempty(hotspot_values) || length(unique(hotspot_values)) == 1
        relative_tolerance = Float64(relative_tolerances[response_id])
        absolute_tolerance = Float64(get(absolute_tolerances,response_id,0.0))
        convergence_pass = fine_change <= relative_tolerance ||
            abs(fine.responses[response_id]-medium.responses[response_id]) <= absolute_tolerance
        if !isnothing(gci)
            convergence_pass &= gci <= relative_tolerance
        end
        require_hotspot_stability && (convergence_pass &= hotspot_stable)

        reference_comparison = nothing
        if haskey(curved_references,response_id)
            candidate_uncertainty = fine.standard_uncertainties[response_id]
            reference_comparison = compare_matched_response(
                fine.responses[response_id],curved_references[response_id];
                candidate_standard_uncertainty=candidate_uncertainty,
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                uncertainty_multiplier=reference_uncertainty_multiplier,
            )
        end
        pass_value = convergence_pass &&
            (isnothing(reference_comparison) || reference_comparison.pass)
        notes = if isnothing(p)
            "Refinement was non-monotone or not approximately uniform; no Richardson/GCI claim."
        else
            "Monotone approximately uniform refinement supported Richardson/GCI estimation."
        end
        estimates[response_id] = Curvature_Response_Estimate(
            response_id,coarse_change,fine_change,p,limit,gci,uncertainty,
            hotspot_stable,reference_comparison,pass_value,notes,
        )
    end

    all_pass = all(estimate.pass for estimate in values(estimates))
    physical_reference_present = all(
        haskey(curved_references,response_id) &&
        reference_is_physical(curved_references[response_id])
        for response_id in response_ids
    )
    production_pass = all_pass && physical_reference_present
    metadata_string = Dict{String,String}(
        "classification" => production_pass ? "physical-curvature-qualified" :
            "software-convergence-only",
        "minimum_levels" => "3",
        "gci_safety_factor" => "1.25",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Curvature_Convergence_Result(
        "radiant.hts.curvature_convergence/v1",level_vector,estimates,
        source_hashes[1],all_pass,physical_reference_present,production_pass,
        metadata_string,
    )
end

function curvature_convergence_receipt(result::Curvature_Convergence_Result)
    responses = Dict{String,Any}()
    for (key,estimate) in result.estimates
        responses[key] = Dict{String,Any}(
            "coarse_to_medium_relative_change" => estimate.coarse_to_medium_change,
            "medium_to_fine_relative_change" => estimate.medium_to_fine_change,
            "observed_order" => isnothing(estimate.observed_order) ? "unavailable" :
                                estimate.observed_order,
            "richardson_limit" => isnothing(estimate.richardson_limit) ? "unavailable" :
                                  estimate.richardson_limit,
            "fine_grid_gci" => isnothing(estimate.fine_grid_gci) ? "unavailable" :
                               estimate.fine_grid_gci,
            "hotspot_stable" => estimate.hotspot_stable,
            "curved_reference_present" => !isnothing(estimate.curved_reference_comparison),
            "pass" => estimate.pass,
            "notes" => estimate.notes,
        )
    end
    return Dict{String,Any}(
        "schema" => result.schema,
        "source_hash" => result.common_source_hash,
        "level_geometry_hashes" => getfield.(result.levels,:geometry_hash),
        "all_responses_pass" => result.all_responses_pass,
        "physical_curved_reference_present" => result.physical_curved_reference_present,
        "production_pass" => result.production_pass,
        "responses" => responses,
        "metadata" => copy(result.metadata),
    )
end

export Curvature_Refinement_Level,Curvature_Response_Estimate
export Curvature_Convergence_Result,qualify_curvature_convergence
export curvature_convergence_receipt
