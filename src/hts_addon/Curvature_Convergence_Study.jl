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
            key_string = string(key)
            haskey(response_values,key_string) || error(
                "Hotspot identity was supplied for an unknown response.",
            )
            isempty(string(value)) && error("Hotspot identifiers cannot be empty.")
            hotspots[key_string] = string(value)
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

struct Curvature_Response_Convergence
    response_id::String
    level_values::Vector{Float64}
    level_standard_uncertainties::Vector{Float64}
    coarse_to_medium_relative_change::Float64
    medium_to_fine_relative_change::Float64
    monotonic::Bool
    observed_order::Union{Nothing,Float64}
    richardson_limit::Union{Nothing,Float64}
    fine_grid_convergence_index::Union{Nothing,Float64}
    hotspot_stable::Bool
    software_pass::Bool
    physical_comparison::Union{Nothing,Matched_Response_Comparison}
    physical_pass::Bool
end

struct Curvature_Convergence_Study
    schema::String
    levels::Vector{Curvature_Refinement_Level}
    responses::Dict{String,Curvature_Response_Convergence}
    source_hash::String
    response_relative_tolerance::Float64
    response_absolute_tolerance::Float64
    uncertainty_multiplier::Float64
    software_pass::Bool
    physical_reference_present::Bool
    production_pass::Bool
    metadata::Dict{String,String}
end

function _curvature_relative_change(current::Float64,previous::Float64)
    return abs(current-previous)/max(abs(current),abs(previous),eps(Float64))
end

function _hotspot_stable(levels::Vector{Curvature_Refinement_Level},response_id::String)
    supplied = [
        level.hotspot_ids[response_id] for level in levels
        if haskey(level.hotspot_ids,response_id)
    ]
    isempty(supplied) && return true
    length(supplied) == length(levels) || return false
    return length(unique(supplied)) == 1
end

function _observed_order_and_limit(
    h_coarse::Float64,
    h_medium::Float64,
    h_fine::Float64,
    coarse::Float64,
    medium::Float64,
    fine::Float64;
    ratio_tolerance::Real=0.05,
)
    ratio_coarse = h_coarse/h_medium
    ratio_fine = h_medium/h_fine
    ratio_coarse > 1.0 && ratio_fine > 1.0 || return nothing,nothing,nothing
    isapprox(ratio_coarse,ratio_fine;rtol=ratio_tolerance,atol=0.0) ||
        return nothing,nothing,nothing
    difference_coarse = medium-coarse
    difference_fine = fine-medium
    difference_coarse == 0.0 || difference_fine == 0.0 ||
        signbit(difference_coarse) != signbit(difference_fine) &&
        return nothing,nothing,nothing
    ratio = sqrt(ratio_coarse*ratio_fine)
    order = log(abs(difference_coarse/difference_fine))/log(ratio)
    isfinite(order) && order > 0.0 || return nothing,nothing,nothing
    denominator = ratio^order-1.0
    denominator > 0.0 || return nothing,nothing,nothing
    limit = fine+difference_fine/denominator
    gci = 1.25*abs(difference_fine)/max(abs(fine),eps(Float64))/denominator
    return order,limit,gci
end

"""
    evaluate_curvature_convergence(levels; ...)

Evaluate at least three nested facet levels. Levels are sorted coarse-to-fine by characteristic
facet size. A software pass requires stable response keys, source identity, decreasing sagitta,
medium-to-fine response closure, and hotspot stability. Production pass additionally requires a
matched physical curved reference (OpenSn, continuous-energy OpenMC, or experiment are sufficient;
Geant4 is not mandatory).
"""
function evaluate_curvature_convergence(
    levels_input::AbstractVector{Curvature_Refinement_Level};
    response_relative_tolerance::Real,
    response_absolute_tolerance::Real=0.0,
    uncertainty_multiplier::Real=2.0,
    physical_references::Union{Nothing,AbstractDict}=nothing,
    require_monotonic::Bool=true,
    metadata::AbstractDict=Dict{String,String}(),
)
    levels = sort(Curvature_Refinement_Level[levels_input...];
                  by=level -> level.characteristic_facet_size_cm,rev=true)
    length(levels) >= 3 || error("Curvature convergence requires at least three facet levels.")
    level_ids = getfield.(levels,:level_id)
    length(unique(level_ids)) == length(level_ids) || error(
        "Curvature refinement level identifiers must be unique.",
    )
    h = getfield.(levels,:characteristic_facet_size_cm)
    all(diff(h) .< 0.0) || error("Facet sizes must decrease strictly from coarse to fine.")
    sagittas = getfield.(levels,:maximum_sagitta_cm)
    all(diff(sagittas) .<= 0.0) || error(
        "Maximum sagitta must not increase under facet refinement.",
    )
    source_hashes = unique(getfield.(levels,:source_hash))
    length(source_hashes) == 1 || error(
        "Every curvature level must use the same source artifact.",
    )
    response_keys = sort(collect(keys(levels[1].responses)))
    for level in levels[2:end]
        sort(collect(keys(level.responses))) == response_keys || error(
            "Curvature levels do not expose identical protected responses.",
        )
    end
    rtol = Float64(response_relative_tolerance)
    atol = Float64(response_absolute_tolerance)
    multiplier = Float64(uncertainty_multiplier)
    all(value -> isfinite(value) && value >= 0.0,(rtol,atol,multiplier)) || error(
        "Curvature convergence tolerances must be finite and nonnegative.",
    )

    response_results = Dict{String,Curvature_Response_Convergence}()
    physical_present = !isnothing(physical_references)
    for response_id in response_keys
        values = [level.responses[response_id] for level in levels]
        uncertainty = [level.standard_uncertainties[response_id] for level in levels]
        coarse_to_medium = _curvature_relative_change(values[2],values[1])
        medium_to_fine = _curvature_relative_change(values[end],values[end-1])
        difference_coarse = values[end-1]-values[end-2]
        difference_fine = values[end]-values[end-1]
        monotonic = difference_coarse == 0.0 || difference_fine == 0.0 ||
                    signbit(difference_coarse) == signbit(difference_fine)
        order,limit,gci = _observed_order_and_limit(
            h[end-2],h[end-1],h[end],values[end-2],values[end-1],values[end],
        )
        hotspot_stable = _hotspot_stable(levels,response_id)
        combined_uncertainty = sqrt(
            uncertainty[end]^2+uncertainty[end-1]^2,
        )
        allowed_difference = atol+rtol*abs(values[end])+multiplier*combined_uncertainty
        software_pass = abs(values[end]-values[end-1]) <= allowed_difference &&
                        hotspot_stable && (!require_monotonic || monotonic)

        physical_comparison = nothing
        physical_pass = false
        if physical_present
            haskey(physical_references,response_id) || error(
                "Physical curved reference is missing response $(response_id).",
            )
            reference = physical_references[response_id]
            reference isa Physical_Reference_Response || error(
                "Curved physical reference has the wrong type for $(response_id).",
            )
            reference_is_physical(reference) || error(
                "Curvature production qualification requires a physical reference.",
            )
            reference.source_artifact_hash == first(source_hashes) || error(
                "Curved reference source hash does not match refinement levels.",
            )
            length(reference.values) == 1 || error(
                "Scalar curvature convergence requires a scalar physical reference.",
            )
            physical_comparison = compare_matched_response(
                reshape([values[end]],1),reshape([uncertainty[end]],1),reference;
                response_id=response_id,
                process_key=reference.process_key,
                relative_tolerance=rtol,
                absolute_tolerance=atol,
                uncertainty_multiplier=multiplier,
                allow_shape_flattening=true,
            )
            physical_pass = physical_comparison.passed
        end
        response_results[response_id] = Curvature_Response_Convergence(
            response_id,values,uncertainty,coarse_to_medium,medium_to_fine,monotonic,
            order,limit,gci,hotspot_stable,software_pass,physical_comparison,physical_pass,
        )
    end
    software_pass = all(result.software_pass for result in values(response_results))
    production_pass = software_pass && physical_present &&
        all(result.physical_pass for result in values(response_results))
    metadata_string = Dict{String,String}(
        "classification" => "facet-refinement-study",
        "physical_reference_required_for_production" => "true",
        "geant4_required" => "false",
        "accepted_reference_classes" => "opensn,continuous_energy_openmc,experiment,geant4",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Curvature_Convergence_Study(
        "radiant.hts.curvature_convergence/v2",levels,response_results,
        first(source_hashes),rtol,atol,multiplier,software_pass,physical_present,
        production_pass,metadata_string,
    )
end

function curvature_convergence_receipt(study::Curvature_Convergence_Study)
    responses = Dict{String,Any}()
    for (response_id,result) in study.responses
        responses[response_id] = Dict{String,Any}(
            "level_values" => result.level_values,
            "level_standard_uncertainties" => result.level_standard_uncertainties,
            "coarse_to_medium_relative_change" =>
                result.coarse_to_medium_relative_change,
            "medium_to_fine_relative_change" => result.medium_to_fine_relative_change,
            "monotonic" => result.monotonic,
            "observed_order" => isnothing(result.observed_order) ?
                "not-estimated" : result.observed_order,
            "richardson_limit" => isnothing(result.richardson_limit) ?
                "not-estimated" : result.richardson_limit,
            "fine_grid_convergence_index" =>
                isnothing(result.fine_grid_convergence_index) ?
                "not-estimated" : result.fine_grid_convergence_index,
            "hotspot_stable" => result.hotspot_stable,
            "software_pass" => result.software_pass,
            "physical_pass" => result.physical_pass,
            "physical_reference_result_hash" => isnothing(result.physical_comparison) ?
                "not-supplied" : result.physical_comparison.metadata["result_artifact_hash"],
        )
    end
    return Dict{String,Any}(
        "schema" => study.schema,
        "source_hash" => study.source_hash,
        "level_ids" => getfield.(study.levels,:level_id),
        "characteristic_facet_size_cm" =>
            getfield.(study.levels,:characteristic_facet_size_cm),
        "maximum_sagitta_cm" => getfield.(study.levels,:maximum_sagitta_cm),
        "responses" => responses,
        "software_pass" => study.software_pass,
        "physical_reference_present" => study.physical_reference_present,
        "production_pass" => study.production_pass,
        "metadata" => copy(study.metadata),
    )
end

export Curvature_Refinement_Level,Curvature_Response_Convergence
export Curvature_Convergence_Study,evaluate_curvature_convergence
export curvature_convergence_receipt
