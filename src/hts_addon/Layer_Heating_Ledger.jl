const HTS_HEATING_SOURCE_CLASSES = (
    :external_prompt,
    :neutron_induced_prompt,
    :gd_capture_prompt,
    :activation_delayed,
    :thermalization_release,
    :joule,
    :other,
)

const HTS_HEATING_TIME_CLASSES = (:prompt,:delayed,:dynamic)

const HTS_HEATING_RESPONSE_CHANNELS = (
    :total_deposition,
    :prompt_lattice_heat,
    :delayed_lattice_heat,
    :non_equilibrium,
    :defect_storage,
    :chemical_storage,
    :escaped,
    :cutoff_handoff,
    :recoil_handoff,
    :unresolved,
    :joule_heat,
)

"""
    Layer_Heating_Contribution

One response field from one uniquely identified physical particle population. `population_id`
identifies the population independently of the code that scored it. This is the key used to stop
OpenSn and Radiant from both contributing additive heating for the same photons.

Comparison-only fields may coexist with a production field, but they are never included in
`production_heating_total`.
"""
struct Layer_Heating_Contribution
    contribution_id::String
    population_id::String
    owner::String
    domain_id::String
    layer_id::String
    particle_tag::String
    source_class::Symbol
    time_class::Symbol
    response_channel::Symbol
    values::Array{Float64}
    standard_uncertainty::Array{Float64}
    units::String
    source_hash::String
    geometry_hash::String
    normalization_hash::String
    production_additive::Bool
    comparison_only::Bool
    metadata::Dict{String,String}

    function Layer_Heating_Contribution(;
        contribution_id::AbstractString,
        population_id::AbstractString,
        owner::AbstractString,
        domain_id::AbstractString,
        layer_id::AbstractString,
        particle_tag::AbstractString,
        source_class::Symbol,
        time_class::Symbol,
        response_channel::Symbol,
        values::AbstractArray{<:Real},
        standard_uncertainty::Union{Nothing,AbstractArray{<:Real}}=nothing,
        units::AbstractString,
        source_hash::AbstractString,
        geometry_hash::AbstractString,
        normalization_hash::AbstractString,
        production_additive::Bool=true,
        comparison_only::Bool=false,
        metadata::AbstractDict=Dict{String,String}(),
    )
        for (label,value) in (
            ("contribution ID",contribution_id),("population ID",population_id),
            ("owner",owner),("domain ID",domain_id),("layer ID",layer_id),
            ("particle tag",particle_tag),("units",units),("source hash",source_hash),
            ("geometry hash",geometry_hash),("normalization hash",normalization_hash),
        )
            isempty(value) && error("Layer-heating $(label) cannot be empty.")
        end
        source_class in HTS_HEATING_SOURCE_CLASSES || error(
            "Unknown layer-heating source class: $(source_class).",
        )
        time_class in HTS_HEATING_TIME_CLASSES || error(
            "Unknown layer-heating time class: $(time_class).",
        )
        response_channel in HTS_HEATING_RESPONSE_CHANNELS || error(
            "Unknown layer-heating response channel: $(response_channel).",
        )
        comparison_only && production_additive && error(
            "A comparison-only heating response cannot be production-additive.",
        )
        data = Float64.(values)
        isempty(data) && error("Layer-heating response field cannot be empty.")
        all(isfinite,data) || error("Layer-heating response values must be finite.")
        all(value -> value >= 0.0,data) || error(
            "Layer-heating response values must be nonnegative.",
        )
        uncertainty = isnothing(standard_uncertainty) ? zeros(Float64,size(data)) :
            Float64.(standard_uncertainty)
        size(uncertainty) == size(data) || error(
            "Layer-heating uncertainty must match the response field shape.",
        )
        all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
            "Layer-heating uncertainty must be finite and nonnegative.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(contribution_id),String(population_id),String(owner),String(domain_id),
            String(layer_id),String(particle_tag),source_class,time_class,response_channel,data,
            uncertainty,String(units),String(source_hash),String(geometry_hash),
            String(normalization_hash),production_additive,comparison_only,metadata_string,
        )
    end
end

"""Mutable collection enforcing one additive owner per population/domain/layer/response/time key."""
mutable struct Layer_Heating_Ledger
    schema::String
    contributions::Dict{String,Layer_Heating_Contribution}
    field_shape::Union{Nothing,Tuple}
    units::Union{Nothing,String}
    metadata::Dict{String,String}

    function Layer_Heating_Ledger(;
        metadata::AbstractDict=Dict{String,String}(),
    )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            "radiant.hts.layer_heating_ledger/v1",
            Dict{String,Layer_Heating_Contribution}(),nothing,nothing,metadata_string,
        )
    end
end

function heating_ownership_key(contribution::Layer_Heating_Contribution)
    return (
        contribution.population_id,
        contribution.domain_id,
        contribution.layer_id,
        contribution.response_channel,
        contribution.time_class,
    )
end

function add_heating_contribution!(
    ledger::Layer_Heating_Ledger,
    contribution::Layer_Heating_Contribution,
)
    haskey(ledger.contributions,contribution.contribution_id) && error(
        "Duplicate layer-heating contribution ID: $(contribution.contribution_id).",
    )
    if isnothing(ledger.field_shape)
        ledger.field_shape = size(contribution.values)
        ledger.units = contribution.units
    else
        size(contribution.values) == ledger.field_shape || error(
            "All layer-heating fields in one ledger must have the same shape.",
        )
        contribution.units == ledger.units || error(
            "All layer-heating fields in one ledger must have identical units.",
        )
    end
    if contribution.production_additive
        key = heating_ownership_key(contribution)
        for existing in values(ledger.contributions)
            if existing.production_additive && heating_ownership_key(existing) == key
                error(
                    "Heating population $(contribution.population_id) already has additive " *
                    "owner $(existing.owner) for $(contribution.response_channel) in layer " *
                    "$(contribution.layer_id); attempted second owner $(contribution.owner).",
                )
            end
        end
    end
    ledger.contributions[contribution.contribution_id] = contribution
    return ledger
end

function get_heating_contribution(
    ledger::Layer_Heating_Ledger,
    contribution_id::AbstractString,
)
    key = String(contribution_id)
    haskey(ledger.contributions,key) || error(
        "Unknown layer-heating contribution: $(key).",
    )
    return ledger.contributions[key]
end

function _heating_filter_match(
    contribution::Layer_Heating_Contribution;
    response_channels::Union{Nothing,Tuple,AbstractVector}=nothing,
    time_class::Union{Nothing,Symbol}=nothing,
    layer_id::Union{Nothing,AbstractString}=nothing,
    population_id::Union{Nothing,AbstractString}=nothing,
)
    !isnothing(response_channels) &&
        !(contribution.response_channel in response_channels) && return false
    !isnothing(time_class) && contribution.time_class != time_class && return false
    !isnothing(layer_id) && contribution.layer_id != String(layer_id) && return false
    !isnothing(population_id) &&
        contribution.population_id != String(population_id) && return false
    return true
end

"""
    production_heating_total(ledger; ...)

Sum only production-additive, non-comparison contributions. By default this includes prompt and
delayed lattice heat plus Joule heat. Total deposition, defect storage, escape, and handoff fields
are excluded unless requested explicitly.
"""
function production_heating_total(
    ledger::Layer_Heating_Ledger;
    response_channels::Union{Tuple,AbstractVector}=(
        :prompt_lattice_heat,:delayed_lattice_heat,:joule_heat,
    ),
    time_class::Union{Nothing,Symbol}=nothing,
    layer_id::Union{Nothing,AbstractString}=nothing,
    population_id::Union{Nothing,AbstractString}=nothing,
)
    isnothing(ledger.field_shape) && error("Cannot sum an empty layer-heating ledger.")
    total = zeros(Float64,ledger.field_shape)
    variance = zeros(Float64,ledger.field_shape)
    selected = String[]
    for contribution in values(ledger.contributions)
        contribution.production_additive || continue
        contribution.comparison_only && error(
            "Ledger invariant broken: comparison-only contribution is additive.",
        )
        _heating_filter_match(
            contribution;response_channels=response_channels,time_class=time_class,
            layer_id=layer_id,population_id=population_id,
        ) || continue
        total .+= contribution.values
        variance .+= contribution.standard_uncertainty.^2
        push!(selected,contribution.contribution_id)
    end
    return (
        values=total,
        standard_uncertainty=sqrt.(variance),
        units=ledger.units,
        contribution_ids=sort(selected),
    )
end

function heating_comparison_set(
    ledger::Layer_Heating_Ledger,
    population_id::AbstractString,
    response_channel::Symbol;
    layer_id::Union{Nothing,AbstractString}=nothing,
)
    output = Layer_Heating_Contribution[]
    for contribution in values(ledger.contributions)
        contribution.population_id == String(population_id) || continue
        contribution.response_channel == response_channel || continue
        !isnothing(layer_id) && contribution.layer_id != String(layer_id) && continue
        push!(output,contribution)
    end
    sort!(output;by=value -> (value.comparison_only,value.owner,value.contribution_id))
    return output
end

"""Check that deposition is partitioned without silently converting storage or handoff to heat."""
function validate_population_energy_closure(
    ledger::Layer_Heating_Ledger,
    population_id::AbstractString;
    component_channels::Tuple=(
        :prompt_lattice_heat,:delayed_lattice_heat,:non_equilibrium,:defect_storage,
        :chemical_storage,:escaped,:cutoff_handoff,:recoil_handoff,:unresolved,
    ),
    rtol::Real=1.0e-10,
    atol::Real=1.0e-12,
)
    population = String(population_id)
    totals = [
        contribution for contribution in values(ledger.contributions)
        if contribution.population_id == population &&
           contribution.response_channel == :total_deposition &&
           contribution.production_additive
    ]
    length(totals) == 1 || error(
        "Population $(population) requires exactly one additive total-deposition field.",
    )
    total = first(totals)
    reconstructed = zeros(Float64,size(total.values))
    component_ids = String[]
    for contribution in values(ledger.contributions)
        contribution.population_id == population || continue
        contribution.production_additive || continue
        contribution.response_channel in component_channels || continue
        heating_ownership_key(contribution)[1:3] == heating_ownership_key(total)[1:3] ||
            error("Population energy components do not share domain and layer identity.")
        reconstructed .+= contribution.values
        push!(component_ids,contribution.contribution_id)
    end
    isempty(component_ids) && error(
        "Population $(population) has no additive energy-partition components.",
    )
    isapprox(reconstructed,total.values;rtol=rtol,atol=atol) || error(
        "Population $(population) energy partition does not reconstruct total deposition.",
    )
    return (
        closed=true,
        total_contribution_id=total.contribution_id,
        component_contribution_ids=sort(component_ids),
        maximum_absolute_residual=maximum(abs.(reconstructed-total.values)),
    )
end

function validate_heating_ledger(ledger::Layer_Heating_Ledger)
    ids = collect(keys(ledger.contributions))
    length(ids) == length(unique(ids)) || error("Layer-heating IDs are not unique.")
    production_keys = Set{Tuple}()
    for contribution in values(ledger.contributions)
        contribution.comparison_only && contribution.production_additive && error(
            "Comparison-only layer heating cannot be additive.",
        )
        if contribution.production_additive
            key = heating_ownership_key(contribution)
            key in production_keys && error("Duplicate additive heating ownership key.")
            push!(production_keys,key)
        end
    end
    return true
end

function heating_ledger_receipt(ledger::Layer_Heating_Ledger)
    validate_heating_ledger(ledger)
    contributions = collect(values(ledger.contributions))
    return Dict{String,Any}(
        "schema" => ledger.schema,
        "contribution_count" => length(contributions),
        "production_additive_count" => count(value -> value.production_additive,contributions),
        "comparison_only_count" => count(value -> value.comparison_only,contributions),
        "population_ids" => sort(unique(getfield.(contributions,:population_id))),
        "owners" => sort(unique(getfield.(contributions,:owner))),
        "layers" => sort(unique(getfield.(contributions,:layer_id))),
        "response_channels" => sort(string.(unique(getfield.(contributions,:response_channel)))),
        "field_shape" => isnothing(ledger.field_shape) ? Int64[] : collect(ledger.field_shape),
        "units" => isnothing(ledger.units) ? "unbound" : ledger.units,
        "comparison_fields_excluded_from_production_sum" => true,
        "one_additive_owner_per_population_domain_layer_response_time" => true,
        "metadata" => copy(ledger.metadata),
    )
end
