const HTS_UNCERTAINTY_CATEGORIES = (
    :global_monte_carlo,
    :source_bank_sampling,
    :source_projection,
    :energy_condensation,
    :angular_discretization,
    :spatial_discretization,
    :facet_curvature,
    :field_map,
    :nuclear_atomic_data,
    :material_property,
    :subkev_model,
    :atomistic_model,
    :coupling_iteration,
    :electrothermal_model,
    :time_discretization,
    :other,
)

const HTS_CONVERGENCE_AXES = (
    :source_bank,
    :spatial,
    :angular,
    :energy,
    :scattering_order,
    :facet,
    :field_map,
    :material_table,
    :subkev_table,
    :coupling,
    :time,
    :other,
)

"""
    Response_Uncertainty_Component

One standard-uncertainty field for one protected response. Components sharing a nonempty
`correlation_group` are combined linearly within that group; different groups are combined in
quadrature. Use a unique group for independent components.
"""
struct Response_Uncertainty_Component
    component_id::String
    category::Symbol
    standard_uncertainty::Array{Float64}
    correlation_group::String
    artifact_hash::String
    classification::Symbol
    metadata::Dict{String,String}

    function Response_Uncertainty_Component(;
        component_id::AbstractString,
        category::Symbol,
        standard_uncertainty::AbstractArray{<:Real},
        correlation_group::AbstractString=component_id,
        artifact_hash::AbstractString,
        classification::Symbol=:software,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(component_id) && error("Uncertainty component ID cannot be empty.")
        category in HTS_UNCERTAINTY_CATEGORIES || error(
            "Unknown HTS uncertainty category: $(category).",
        )
        uncertainty = Float64.(standard_uncertainty)
        isempty(uncertainty) && error("Uncertainty component field cannot be empty.")
        all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
            "Standard uncertainty must be finite and nonnegative.",
        )
        isempty(correlation_group) && error("Uncertainty correlation group cannot be empty.")
        isempty(artifact_hash) && error("Uncertainty artifact hash cannot be empty.")
        classification in (:software,:analytic,:physical,:experimental) || error(
            "Unknown uncertainty classification: $(classification).",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(component_id),category,uncertainty,String(correlation_group),
            String(artifact_hash),classification,metadata_string,
        )
    end
end

mutable struct Protected_Response_Uncertainty_Budget
    schema::String
    response_id::String
    nominal_values::Array{Float64}
    units::String
    components::Dict{String,Response_Uncertainty_Component}
    result_artifact_hash::String
    metadata::Dict{String,String}

    function Protected_Response_Uncertainty_Budget(
        response_id::AbstractString,
        nominal_values::AbstractArray{<:Real};
        units::AbstractString,
        result_artifact_hash::AbstractString,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(response_id) && error("Protected response ID cannot be empty.")
        nominal = Float64.(nominal_values)
        isempty(nominal) && error("Protected response field cannot be empty.")
        all(isfinite,nominal) || error("Protected response values must be finite.")
        isempty(units) && error("Protected response units cannot be empty.")
        isempty(result_artifact_hash) && error("Protected response result hash cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            "radiant.hts.protected_response_uncertainty/v1",String(response_id),nominal,
            String(units),Dict{String,Response_Uncertainty_Component}(),
            String(result_artifact_hash),metadata_string,
        )
    end
end

function add_uncertainty_component!(
    budget::Protected_Response_Uncertainty_Budget,
    component::Response_Uncertainty_Component,
)
    haskey(budget.components,component.component_id) && error(
        "Duplicate uncertainty component ID: $(component.component_id).",
    )
    size(component.standard_uncertainty) == size(budget.nominal_values) || error(
        "Uncertainty component shape does not match the protected response.",
    )
    budget.components[component.component_id] = component
    return budget
end

function combined_standard_uncertainty(
    budget::Protected_Response_Uncertainty_Budget;
    categories::Union{Nothing,Tuple,AbstractVector}=nothing,
)
    isempty(budget.components) && error("Protected response uncertainty budget is empty.")
    grouped = Dict{String,Array{Float64}}()
    for component in values(budget.components)
        !isnothing(categories) && !(component.category in categories) && continue
        if haskey(grouped,component.correlation_group)
            grouped[component.correlation_group] .+= component.standard_uncertainty
        else
            grouped[component.correlation_group] = copy(component.standard_uncertainty)
        end
    end
    isempty(grouped) && error("No uncertainty components match the requested categories.")
    variance = zeros(Float64,size(budget.nominal_values))
    for group_uncertainty in values(grouped)
        variance .+= group_uncertainty.^2
    end
    return sqrt.(variance)
end

function relative_combined_uncertainty(
    budget::Protected_Response_Uncertainty_Budget;
    absolute_floor::Real=eps(Float64),
)
    floor_value = Float64(absolute_floor)
    isfinite(floor_value) && floor_value > 0.0 || error(
        "Relative-uncertainty floor must be finite and positive.",
    )
    denominator = max.(abs.(budget.nominal_values),floor_value)
    return combined_standard_uncertainty(budget)./denominator
end

function uncertainty_budget_receipt(budget::Protected_Response_Uncertainty_Budget)
    combined = combined_standard_uncertainty(budget)
    return Dict{String,Any}(
        "schema" => budget.schema,
        "response_id" => budget.response_id,
        "units" => budget.units,
        "field_shape" => collect(size(budget.nominal_values)),
        "component_count" => length(budget.components),
        "categories" => sort(string.(unique(getfield.(collect(values(budget.components)),:category)))),
        "correlation_groups" => sort(unique(
            getfield.(collect(values(budget.components)),:correlation_group),
        )),
        "maximum_combined_standard_uncertainty" => maximum(combined),
        "maximum_relative_combined_uncertainty" => maximum(
            combined./max.(abs.(budget.nominal_values),eps(Float64)),
        ),
        "result_artifact_hash" => budget.result_artifact_hash,
        "all_components_physical" => all(
            component.classification in (:physical,:experimental)
            for component in values(budget.components)
        ),
        "metadata" => copy(budget.metadata),
    )
end

"""One member of a coarse-to-fine protected-response sequence."""
struct Protected_Response_Convergence_Level
    level_id::String
    characteristic_step::Float64
    values::Array{Float64}
    standard_uncertainty::Array{Float64}
    artifact_hash::String
    metadata::Dict{String,String}

    function Protected_Response_Convergence_Level(;
        level_id::AbstractString,
        characteristic_step::Real,
        values::AbstractArray{<:Real},
        standard_uncertainty::Union{Nothing,AbstractArray{<:Real}}=nothing,
        artifact_hash::AbstractString,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(level_id) && error("Convergence level ID cannot be empty.")
        step = Float64(characteristic_step)
        isfinite(step) && step > 0.0 || error(
            "Convergence characteristic step must be finite and positive.",
        )
        data = Float64.(values)
        isempty(data) && error("Convergence response field cannot be empty.")
        all(isfinite,data) || error("Convergence response values must be finite.")
        uncertainty = isnothing(standard_uncertainty) ? zeros(Float64,size(data)) :
            Float64.(standard_uncertainty)
        size(uncertainty) == size(data) || error(
            "Convergence uncertainty shape does not match the response field.",
        )
        all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
            "Convergence standard uncertainty must be finite and nonnegative.",
        )
        isempty(artifact_hash) && error("Convergence-level artifact hash cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(level_id),step,data,uncertainty,String(artifact_hash),metadata_string,
        )
    end
end

struct Protected_Response_Convergence_Result
    schema::String
    response_id::String
    axis::Symbol
    level_ids::Vector{String}
    characteristic_steps::Vector{Float64}
    coarse_to_fine_changes::Vector{Float64}
    finest_absolute_difference::Array{Float64}
    finest_allowed_difference::Array{Float64}
    finest_pass_mask::BitArray
    maximum_finest_relative_change::Float64
    observed_order::Union{Nothing,Float64}
    passed::Bool
    physical_reference::Bool
    source_hash::String
    geometry_hash::String
    metadata::Dict{String,String}
end

function _maximum_relative_array_change(
    first::AbstractArray{<:Real},
    second::AbstractArray{<:Real};
    floor_value::Real=eps(Float64),
)
    size(first) == size(second) || error("Response convergence array shapes differ.")
    denominator = max.(max.(abs.(Float64.(first)),abs.(Float64.(second))),Float64(floor_value))
    return maximum(abs.(Float64.(first).-Float64.(second))./denominator)
end

function _observed_convergence_order(
    levels::Vector{Protected_Response_Convergence_Level},
)
    length(levels) < 3 && return nothing
    coarse,medium,fine = levels[end-2],levels[end-1],levels[end]
    ratio_1 = coarse.characteristic_step/medium.characteristic_step
    ratio_2 = medium.characteristic_step/fine.characteristic_step
    isapprox(ratio_1,ratio_2;rtol=0.05,atol=0.0) || return nothing
    difference_1 = norm(coarse.values-medium.values)
    difference_2 = norm(medium.values-fine.values)
    difference_1 > 0.0 && difference_2 > 0.0 || return nothing
    order = log(difference_1/difference_2)/log(0.5*(ratio_1+ratio_2))
    return isfinite(order) ? order : nothing
end

"""
    evaluate_response_convergence(levels; ...)

Evaluate the last two levels with absolute, relative, and combined statistical-uncertainty
allowances. `characteristic_step` must strictly decrease from coarse to fine. The function may
record a physical reference classification but never infers it from numerical convergence alone.
"""
function evaluate_response_convergence(
    levels::AbstractVector{Protected_Response_Convergence_Level};
    response_id::AbstractString,
    axis::Symbol,
    relative_tolerance::Real,
    absolute_tolerance::Real,
    uncertainty_multiplier::Real=2.0,
    minimum_levels::Integer=3,
    physical_reference::Bool=false,
    source_hash::AbstractString,
    geometry_hash::AbstractString,
    metadata::AbstractDict=Dict{String,String}(),
)
    isempty(response_id) && error("Convergence response ID cannot be empty.")
    axis in HTS_CONVERGENCE_AXES || error("Unknown convergence axis: $(axis).")
    level_vector = Protected_Response_Convergence_Level[levels...]
    length(level_vector) >= minimum_levels || error(
        "Protected response convergence requires at least $(minimum_levels) levels.",
    )
    identifiers = getfield.(level_vector,:level_id)
    length(unique(identifiers)) == length(identifiers) || error(
        "Convergence level IDs must be unique.",
    )
    steps = getfield.(level_vector,:characteristic_step)
    all(diff(steps) .< 0.0) || error(
        "Convergence characteristic steps must strictly decrease from coarse to fine.",
    )
    reference_shape = size(level_vector[1].values)
    all(level -> size(level.values) == reference_shape,level_vector) || error(
        "All convergence response fields must have the same shape.",
    )
    for value in (source_hash,geometry_hash)
        isempty(value) && error("Convergence lineage hashes cannot be empty.")
    end
    rtol = Float64(relative_tolerance)
    atol = Float64(absolute_tolerance)
    multiplier = Float64(uncertainty_multiplier)
    all(value -> isfinite(value) && value >= 0.0,(rtol,atol,multiplier)) || error(
        "Convergence tolerances must be finite and nonnegative.",
    )
    changes = Float64[]
    for index in 2:length(level_vector)
        push!(changes,_maximum_relative_array_change(
            level_vector[index-1].values,level_vector[index].values,
        ))
    end
    previous = level_vector[end-1]
    finest = level_vector[end]
    difference = abs.(finest.values-previous.values)
    combined_uncertainty = sqrt.(
        finest.standard_uncertainty.^2+previous.standard_uncertainty.^2,
    )
    allowed = atol .+ rtol.*abs.(finest.values) .+ multiplier.*combined_uncertainty
    pass_mask = difference .<= allowed
    metadata_string = Dict{String,String}(
        "relative_tolerance" => string(rtol),
        "absolute_tolerance" => string(atol),
        "uncertainty_multiplier" => string(multiplier),
        "numerical_convergence_does_not_imply_physical_validation" => "true",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Protected_Response_Convergence_Result(
        "radiant.hts.protected_response_convergence/v1",String(response_id),axis,
        identifiers,steps,changes,difference,allowed,BitArray(pass_mask),last(changes),
        _observed_convergence_order(level_vector),all(pass_mask),physical_reference,
        String(source_hash),String(geometry_hash),metadata_string,
    )
end

function convergence_result_receipt(result::Protected_Response_Convergence_Result)
    return Dict{String,Any}(
        "schema" => result.schema,
        "response_id" => result.response_id,
        "axis" => string(result.axis),
        "level_ids" => result.level_ids,
        "characteristic_steps" => result.characteristic_steps,
        "coarse_to_fine_relative_changes" => result.coarse_to_fine_changes,
        "maximum_finest_relative_change" => result.maximum_finest_relative_change,
        "maximum_finest_absolute_difference" => maximum(result.finest_absolute_difference),
        "minimum_finest_allowed_difference" => minimum(result.finest_allowed_difference),
        "observed_order" => isnothing(result.observed_order) ? "unresolved" :
                            result.observed_order,
        "passed" => result.passed,
        "physical_reference" => result.physical_reference,
        "source_hash" => result.source_hash,
        "geometry_hash" => result.geometry_hash,
        "metadata" => copy(result.metadata),
    )
end

struct Multi_Axis_Response_Qualification
    schema::String
    response_id::String
    studies::Dict{Symbol,Protected_Response_Convergence_Result}
    required_axes::Vector{Symbol}
    missing_axes::Vector{Symbol}
    software_passed::Bool
    physical_passed::Bool
    metadata::Dict{String,String}
end

"""Require independent convergence evidence for every response-protection axis requested."""
function qualify_response_across_axes(
    studies::AbstractVector{Protected_Response_Convergence_Result};
    response_id::AbstractString,
    required_axes::AbstractVector{Symbol},
    require_physical_reference::Bool=false,
    metadata::AbstractDict=Dict{String,String}(),
)
    study_dictionary = Dict{Symbol,Protected_Response_Convergence_Result}()
    for study in studies
        study.response_id == response_id || error(
            "All convergence studies must protect the same response ID.",
        )
        haskey(study_dictionary,study.axis) && error(
            "Duplicate convergence study for axis $(study.axis).",
        )
        study_dictionary[study.axis] = study
    end
    required = Symbol[required_axes...]
    all(axis -> axis in HTS_CONVERGENCE_AXES,required) || error(
        "Required response-convergence axis is invalid.",
    )
    length(unique(required)) == length(required) || error(
        "Required response-convergence axes must be unique.",
    )
    missing = Symbol[axis for axis in required if !haskey(study_dictionary,axis)]
    software_passed = isempty(missing) && all(
        study_dictionary[axis].passed for axis in required
    )
    physical_passed = software_passed && all(
        study_dictionary[axis].physical_reference for axis in required
    )
    require_physical_reference && !physical_passed && error(
        "Multi-axis response qualification lacks physical-reference evidence.",
    )
    metadata_string = Dict{String,String}()
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Multi_Axis_Response_Qualification(
        "radiant.hts.multi_axis_response_qualification/v1",String(response_id),
        study_dictionary,required,missing,software_passed,physical_passed,metadata_string,
    )
end

function multi_axis_qualification_receipt(result::Multi_Axis_Response_Qualification)
    return Dict{String,Any}(
        "schema" => result.schema,
        "response_id" => result.response_id,
        "required_axes" => string.(result.required_axes),
        "missing_axes" => string.(result.missing_axes),
        "software_passed" => result.software_passed,
        "physical_passed" => result.physical_passed,
        "studies" => Dict(
            string(axis) => convergence_result_receipt(study)
            for (axis,study) in result.studies
        ),
        "metadata" => copy(result.metadata),
    )
end
