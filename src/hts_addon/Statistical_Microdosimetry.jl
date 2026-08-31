const MICRODOSIMETRY_PARTITION_CHANNELS = (
    :ionization,
    :electronic_excitation,
    :prompt_lattice_heat,
    :nuclear_recoil_handoff,
    :defect_storage,
    :optical_emission,
    :escaping_particle,
    :cutoff_handoff,
    :unresolved,
)

struct Event_Energy_Partition_Fractions
    fractions::Dict{Symbol,Float64}

    function Event_Energy_Partition_Fractions(
        fractions::AbstractDict;
        closure_tolerance::Real=1.0e-10,
    )
        output = Dict{Symbol,Float64}()
        for channel in MICRODOSIMETRY_PARTITION_CHANNELS
            haskey(fractions,channel) || error(
                "Microdosimetry partition is missing channel $(channel).",
            )
            value = Float64(fractions[channel])
            isfinite(value) && value ≥ 0.0 && value ≤ 1.0 || error(
                "Microdosimetry partition fractions must lie in [0,1].",
            )
            output[channel] = value
        end
        isapprox(
            sum(values(output)),1.0;
            rtol=closure_tolerance,atol=closure_tolerance,
        ) || error("Microdosimetry partition fractions must sum to one.")
        return new(output)
    end
end

function energy_partition_from_fractions(
    fractions::Event_Energy_Partition_Fractions,
    deposited_energy_eV::Real;
    label::AbstractString="microdosimetry-event",
    metadata::AbstractDict=Dict{String,String}(),
)
    energy = Float64(deposited_energy_eV)
    isfinite(energy) && energy ≥ 0.0 || error("Deposited event energy must be nonnegative.")
    f = fractions.fractions
    return Energy_Partition(
        label=label,
        available_energy_eV=energy,
        ionization_eV=energy*f[:ionization],
        electronic_excitation_eV=energy*f[:electronic_excitation],
        prompt_lattice_heat_eV=energy*f[:prompt_lattice_heat],
        nuclear_recoil_handoff_eV=energy*f[:nuclear_recoil_handoff],
        defect_stored_eV=energy*f[:defect_storage],
        optical_emission_eV=energy*f[:optical_emission],
        escaping_particle_eV=energy*f[:escaping_particle],
        cutoff_handoff_eV=energy*f[:cutoff_handoff],
        unresolved_eV=energy*f[:unresolved],
        metadata=metadata,
    )
end

const MICRODOSIMETRY_DIRECTION_MODELS = (
    :fixed,
    :isotropic,
    :tabulated_angular,
    :joint_energy_angle,
    :parent_correlated,
)

function _validate_direction_support(
    direction_model::Symbol,
    fixed_direction::AbstractVector{<:Real},
    directions::AbstractVector,
    probabilities::AbstractVector{<:Real},
    energies_eV::AbstractVector{<:Real};
    allow_parent_correlated::Bool,
    label::AbstractString,
)
    direction_model in MICRODOSIMETRY_DIRECTION_MODELS || error(
        "Unknown $(label) direction model: $(direction_model).",
    )
    direction_model == :parent_correlated && !allow_parent_correlated && error(
        "$(label) cannot use :parent_correlated without an external parent vector.",
    )
    fixed = _atlas_normalize(fixed_direction,"$(label) fixed direction")
    support = NTuple{3,Float64}[]
    for (index,value) in enumerate(directions)
        normalized = _atlas_normalize(value,"$(label) direction support $(index)")
        push!(support,(normalized[1],normalized[2],normalized[3]))
    end
    weights = Float64.(probabilities)
    energies = Float64.(energies_eV)
    tabulated = direction_model in (:tabulated_angular,:joint_energy_angle)
    if tabulated
        !isempty(support) && length(support) == length(weights) || error(
            "$(label) tabulated directions and probabilities must be nonempty and equal length.",
        )
        all(value -> isfinite(value) && value ≥ 0.0,weights) || error(
            "$(label) direction probabilities must be finite and nonnegative.",
        )
        isapprox(sum(weights),1.0;rtol=1.0e-12,atol=1.0e-12) || error(
            "$(label) direction probabilities must sum to one.",
        )
        if direction_model == :joint_energy_angle
            length(energies) == length(support) && all(value -> isfinite(value) && value ≥ 0.0,energies) || error(
                "$(label) joint energy-angle support requires one nonnegative energy per direction.",
            )
        elseif !isempty(energies)
            error("$(label) tabulated angular model cannot carry joint energy support.")
        end
    elseif !isempty(support) || !isempty(weights) || !isempty(energies)
        error("$(label) non-tabulated direction model cannot carry tabulated support.")
    end
    return (
        (fixed[1],fixed[2],fixed[3]),support,weights,energies,
    )
end

"""One secondary emitted in correlation with the parent event prototype."""
struct Correlated_Secondary
    particle_tag::String
    kinetic_energy_eV::Float64
    direction_model::Symbol
    correlation_id::String
    metadata::Dict{String,String}
    fixed_direction::NTuple{3,Float64}
    direction_support::Vector{NTuple{3,Float64}}
    direction_probabilities::Vector{Float64}
    joint_energy_support_eV::Vector{Float64}

    function Correlated_Secondary(
        particle_tag::AbstractString,
        kinetic_energy_eV::Real;
        direction_model::Symbol=:missing,
        correlation_id::AbstractString="uncorrelated",
        metadata::AbstractDict=Dict{String,String}(),
        fixed_direction::AbstractVector{<:Real}=[1.0,0.0,0.0],
        direction_support::AbstractVector=NTuple{3,Float64}[],
        direction_probabilities::AbstractVector{<:Real}=Float64[],
        joint_energy_support_eV::AbstractVector{<:Real}=Float64[],
    )
        isempty(particle_tag) && error("Secondary particle tag cannot be empty.")
        energy = Float64(kinetic_energy_eV)
        isfinite(energy) && energy ≥ 0.0 || error("Secondary energy must be nonnegative.")
        isempty(correlation_id) && error("Secondary correlation identifier cannot be empty.")
        direction_model == :parent_correlated && correlation_id == "uncorrelated" && error(
            "Parent-correlated secondaries require an explicit correlation identifier.",
        )
        fixed,support,weights,energies = _validate_direction_support(
            direction_model,fixed_direction,direction_support,direction_probabilities,
            joint_energy_support_eV;allow_parent_correlated=true,label="secondary",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(particle_tag),energy,direction_model,String(correlation_id),metadata_string,
            fixed,support,weights,energies,
        )
    end
end

"""Positive energy-loss straggling model used by one event prototype."""
struct Energy_Straggling_Model
    distribution::Symbol
    relative_sigma::Float64
    maximum_deposited_energy_eV::Float64
    model_hash::String
    qualification_status::Symbol

    function Energy_Straggling_Model(
        distribution::Symbol;
        relative_sigma::Real=0.0,
        maximum_deposited_energy_eV::Real=Inf,
        model_hash::AbstractString="unbound",
        qualification_status::Symbol=:verification,
    )
        distribution in (:deterministic,:truncated_gaussian,:lognormal) || error(
            "Straggling distribution must be deterministic, truncated_gaussian, or lognormal.",
        )
        sigma = Float64(relative_sigma)
        maximum_energy = Float64(maximum_deposited_energy_eV)
        isfinite(sigma) && sigma ≥ 0.0 || error("Relative straggling sigma must be nonnegative.")
        (isfinite(maximum_energy) || isinf(maximum_energy)) && maximum_energy > 0.0 || error(
            "Maximum deposited energy must be positive.",
        )
        distribution == :deterministic && sigma != 0.0 && error(
            "Deterministic straggling requires zero relative sigma.",
        )
        qualification_status in (:verification,:candidate,:qualified) || error(
            "Straggling qualification status is invalid.",
        )
        isempty(model_hash) && error("Straggling model hash cannot be empty.")
        return new(distribution,sigma,maximum_energy,String(model_hash),qualification_status)
    end
end

"""Representative event outcome with an expected physical rate."""
struct Microdosimetry_Event_Prototype
    prototype_id::String
    particle_tag::String
    process_id::String
    material_tag::String
    layer_id::String
    position_cm::NTuple{3,Float64}
    direction_model::Symbol
    fixed_direction::NTuple{3,Float64}
    direction_support::Vector{NTuple{3,Float64}}
    direction_probabilities::Vector{Float64}
    joint_energy_support_eV::Vector{Float64}
    incident_energy_eV::Float64
    mean_deposited_energy_eV::Float64
    event_rate_per_s::Float64
    straggling::Energy_Straggling_Model
    partition_fractions::Event_Energy_Partition_Fractions
    correlated_secondaries::Vector{Correlated_Secondary}
    correlation_id::String
    provenance::Dict{String,String}

    function Microdosimetry_Event_Prototype(;
        prototype_id::AbstractString,
        particle_tag::AbstractString,
        process_id::AbstractString,
        material_tag::AbstractString,
        layer_id::AbstractString,
        position_cm::AbstractVector{<:Real},
        direction_model::Symbol=:missing,
        fixed_direction::AbstractVector{<:Real}=[1.0,0.0,0.0],
        direction_support::AbstractVector=NTuple{3,Float64}[],
        direction_probabilities::AbstractVector{<:Real}=Float64[],
        joint_energy_support_eV::AbstractVector{<:Real}=Float64[],
        incident_energy_eV::Real,
        mean_deposited_energy_eV::Real,
        event_rate_per_s::Real,
        straggling::Energy_Straggling_Model=Energy_Straggling_Model(:deterministic),
        partition_fractions::Event_Energy_Partition_Fractions,
        correlated_secondaries::AbstractVector{Correlated_Secondary}=Correlated_Secondary[],
        correlation_id::AbstractString="uncorrelated",
        provenance::AbstractDict=Dict{String,String}(),
    )
        for (name,value) in (
            ("prototype",prototype_id),("particle",particle_tag),("process",process_id),
            ("material",material_tag),("layer",layer_id),("correlation",correlation_id),
        )
            isempty(value) && error("Microdosimetry $(name) identifier cannot be empty.")
        end
        position = Float64.(position_cm)
        length(position) == 3 && all(isfinite,position) || error(
            "Microdosimetry position must be a finite three-vector.",
        )
        direction,support,weights,energies = _validate_direction_support(
            direction_model,fixed_direction,direction_support,direction_probabilities,
            joint_energy_support_eV;allow_parent_correlated=false,label="prototype",
        )
        incident = Float64(incident_energy_eV)
        deposited = Float64(mean_deposited_energy_eV)
        rate = Float64(event_rate_per_s)
        isfinite(incident) && incident > 0.0 || error("Incident energy must be positive.")
        isfinite(deposited) && deposited ≥ 0.0 && deposited ≤ incident || error(
            "Mean deposited energy must lie between zero and incident energy.",
        )
        isfinite(rate) && rate ≥ 0.0 || error("Event rate must be nonnegative.")
        straggling.maximum_deposited_energy_eV < deposited && error(
            "Straggling maximum lies below the mean deposited energy.",
        )
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        correlation = String(correlation_id)
        for secondary in correlated_secondaries
            if secondary.direction_model == :parent_correlated
                correlation == "uncorrelated" && error(
                    "A prototype with parent-correlated secondaries requires an explicit correlation identifier.",
                )
                secondary.correlation_id == correlation || error(
                    "Parent-correlated secondary and parent correlation identifiers must match.",
                )
            end
        end
        return new(
            String(prototype_id),String(particle_tag),String(process_id),String(material_tag),
            String(layer_id),(position[1],position[2],position[3]),direction_model,
            direction,support,weights,energies,incident,deposited,rate,straggling,
            partition_fractions,Correlated_Secondary[correlated_secondaries...],
            correlation,provenance_string,
        )
    end
end

struct Microdosimetry_Kernel
    schema::String
    prototypes::Vector{Microdosimetry_Event_Prototype}
    source_artifact_hash::String
    geometry_hash::String
    material_state_hash::String
    normalization_basis::Symbol
    metadata::Dict{String,String}

    function Microdosimetry_Kernel(
        prototypes::AbstractVector{Microdosimetry_Event_Prototype};
        source_artifact_hash::AbstractString,
        geometry_hash::AbstractString,
        material_state_hash::AbstractString,
        normalization_basis::Symbol=:physical_rate,
        metadata::AbstractDict=Dict{String,String}(),
    )
        prototype_vector = Microdosimetry_Event_Prototype[prototypes...]
        isempty(prototype_vector) && error("Microdosimetry kernel cannot be empty.")
        identifiers = getfield.(prototype_vector,:prototype_id)
        length(unique(identifiers)) == length(identifiers) || error(
            "Microdosimetry prototype identifiers must be unique.",
        )
        for value in (source_artifact_hash,geometry_hash,material_state_hash)
            isempty(value) && error("Microdosimetry lineage hashes cannot be empty.")
        end
        normalization_basis in (:physical_rate,:per_history) || error(
            "Microdosimetry normalization must be :physical_rate or :per_history.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            "radiant.hts.microdosimetry_kernel/v1",prototype_vector,
            String(source_artifact_hash),String(geometry_hash),String(material_state_hash),
            normalization_basis,metadata_string,
        )
    end
end

"""Actual sampled secondary state, retaining the parent correlation and vector."""
struct Sampled_Correlated_Secondary
    particle_tag::String
    kinetic_energy_eV::Float64
    direction::NTuple{3,Float64}
    direction_model::Symbol
    correlation_id::String
    metadata::Dict{String,String}
end

struct Weighted_Microdosimetry_Event
    sample_id::Int64
    prototype_id::String
    particle_tag::String
    process_id::String
    material_tag::String
    layer_id::String
    time_s::Float64
    position_cm::NTuple{3,Float64}
    direction::NTuple{3,Float64}
    direction_model::Symbol
    incident_energy_eV::Float64
    deposited_energy_eV::Float64
    partition::Energy_Partition
    statistical_weight_events::Float64
    correlated_secondaries::Vector{Sampled_Correlated_Secondary}
    correlation_id::String
    provenance::Dict{String,String}
end

function _sample_isotropic_direction(rng::Random.AbstractRNG)
    cosine = 2.0*rand(rng)-1.0
    azimuth = 2.0*π*rand(rng)
    sine = sqrt(max(0.0,1.0-cosine^2))
    return (sine*cos(azimuth),sine*sin(azimuth),cosine)
end

function _sample_directional_state(
    rng::Random.AbstractRNG,
    direction_model::Symbol,
    fixed_direction::NTuple{3,Float64},
    direction_support::Vector{NTuple{3,Float64}},
    direction_probabilities::Vector{Float64},
    joint_energy_support_eV::Vector{Float64},
    default_energy_eV::Float64;
    parent_direction::Union{Nothing,NTuple{3,Float64}}=nothing,
)
    direction_model == :fixed && return (default_energy_eV,fixed_direction)
    direction_model == :isotropic && return (default_energy_eV,_sample_isotropic_direction(rng))
    if direction_model == :parent_correlated
        parent_direction === nothing && error(
            "Parent-correlated direction sampling requires the actual sampled parent vector.",
        )
        return (default_energy_eV,parent_direction)
    end
    target = rand(rng)
    cumulative = 0.0
    selected = length(direction_probabilities)
    for index in eachindex(direction_probabilities)
        cumulative += direction_probabilities[index]
        if target ≤ cumulative
            selected = index
            break
        end
    end
    energy = direction_model == :joint_energy_angle ?
        joint_energy_support_eV[selected] : default_energy_eV
    return (energy,direction_support[selected])
end

function _sample_correlated_secondaries(
    rng::Random.AbstractRNG,
    secondaries::AbstractVector{Correlated_Secondary},
    parent_direction::NTuple{3,Float64},
)
    output = Sampled_Correlated_Secondary[]
    for secondary in secondaries
        energy,direction = _sample_directional_state(
            rng,secondary.direction_model,secondary.fixed_direction,
            secondary.direction_support,secondary.direction_probabilities,
            secondary.joint_energy_support_eV,secondary.kinetic_energy_eV;
            parent_direction=parent_direction,
        )
        push!(output,Sampled_Correlated_Secondary(
            secondary.particle_tag,energy,direction,secondary.direction_model,
            secondary.correlation_id,copy(secondary.metadata),
        ))
    end
    return output
end

function _sample_deposited_energy(
    rng::Random.AbstractRNG,
    mean_energy_eV::Float64,
    model::Energy_Straggling_Model;
    maximum_attempts::Integer=10_000,
)
    model.distribution == :deterministic && return mean_energy_eV
    mean_energy_eV == 0.0 && return 0.0
    sigma = model.relative_sigma*mean_energy_eV
    sigma == 0.0 && return mean_energy_eV
    for _ in 1:maximum_attempts
        value = if model.distribution == :truncated_gaussian
            mean_energy_eV+sigma*randn(rng)
        else
            variance = sigma^2
            log_sigma2 = log1p(variance/mean_energy_eV^2)
            log_mu = log(mean_energy_eV)-0.5*log_sigma2
            exp(log_mu+sqrt(log_sigma2)*randn(rng))
        end
        if value ≥ 0.0 && value ≤ model.maximum_deposited_energy_eV
            return value
        end
    end
    error(
        "Straggling sampler could not draw an admissible energy without clipping; " *
        "revise the distribution or physical upper bound.",
    )
end

"""
    sample_microdosimetry_events(kernel, sample_count; time_window_s=(0,1), seed=0)

Systematic weighted sampling of deterministic event rates. Samples are representative weighted
events, not analog histories. Their weights sum to the expected number of physical events in the
requested time window. Correlated secondaries remain attached to their parent prototype.
"""
function sample_microdosimetry_events(
    kernel::Microdosimetry_Kernel,
    sample_count::Integer;
    time_window_s::Tuple{<:Real,<:Real}=(0.0,1.0),
    seed::Integer=0,
)
    sample_count ≥ 1 || error("Microdosimetry sample count must be positive.")
    kernel.normalization_basis == :physical_rate || error(
        "Time-window sampling requires physical event rates.",
    )
    t0 = Float64(time_window_s[1])
    t1 = Float64(time_window_s[2])
    isfinite(t0) && isfinite(t1) && t1 > t0 || error(
        "Microdosimetry time window must be finite and increasing.",
    )
    rates = getfield.(kernel.prototypes,:event_rate_per_s)
    total_rate = sum(rates)
    total_rate > 0.0 || error("Microdosimetry kernel has zero total event rate.")
    cumulative = cumsum(rates)
    rng = Random.MersenneTwister(seed)
    spacing = total_rate/sample_count
    offset = rand(rng)*spacing
    expected_events = total_rate*(t1-t0)
    event_weight = expected_events/sample_count
    output = Weighted_Microdosimetry_Event[]

    prototype_index = 1
    for sample_id in 1:sample_count
        target = offset+(sample_id-1)*spacing
        while prototype_index < length(cumulative) && target > cumulative[prototype_index]
            prototype_index += 1
        end
        prototype = kernel.prototypes[prototype_index]
        deposited = _sample_deposited_energy(
            rng,prototype.mean_deposited_energy_eV,prototype.straggling,
        )
        deposited ≤ prototype.incident_energy_eV || error(
            "Sampled deposited energy exceeds incident energy; no clipping is permitted.",
        )
        sampled_incident,direction = _sample_directional_state(
            rng,prototype.direction_model,prototype.fixed_direction,
            prototype.direction_support,prototype.direction_probabilities,
            prototype.joint_energy_support_eV,prototype.incident_energy_eV,
        )
        deposited ≤ sampled_incident || error(
            "Sampled deposited energy exceeds sampled incident energy; no clipping is permitted.",
        )
        sampled_secondaries = _sample_correlated_secondaries(
            rng,prototype.correlated_secondaries,direction,
        )
        event_time = t0+(t1-t0)*rand(rng)
        partition = energy_partition_from_fractions(
            prototype.partition_fractions,deposited;
            label=prototype.process_id,
            metadata=Dict(
                "prototype_id" => prototype.prototype_id,
                "straggling_model_hash" => prototype.straggling.model_hash,
                "source_artifact_hash" => kernel.source_artifact_hash,
                "geometry_hash" => kernel.geometry_hash,
                "material_state_hash" => kernel.material_state_hash,
            ),
        )
        push!(output,Weighted_Microdosimetry_Event(
            Int64(sample_id),prototype.prototype_id,prototype.particle_tag,
            prototype.process_id,prototype.material_tag,prototype.layer_id,event_time,
            prototype.position_cm,direction,prototype.direction_model,
            sampled_incident,deposited,partition,
            event_weight,sampled_secondaries,prototype.correlation_id,
            merge(
                copy(prototype.provenance),
                Dict(
                    "kernel_schema" => kernel.schema,
                    "sampling" => "systematic-weighted",
                    "seed" => string(seed),
                ),
            ),
        ))
    end
    return output
end

function microdosimetry_effective_sample_size(events::AbstractVector{Weighted_Microdosimetry_Event})
    isempty(events) && return 0.0
    weights = getfield.(events,:statistical_weight_events)
    denominator = sum(abs2,weights)
    denominator == 0.0 && return 0.0
    return sum(weights)^2/denominator
end

function weighted_deposited_energy_mean(events::AbstractVector{Weighted_Microdosimetry_Event})
    isempty(events) && error("Cannot calculate a microdosimetry mean from no events.")
    weights = getfield.(events,:statistical_weight_events)
    total_weight = sum(weights)
    total_weight > 0.0 || error("Microdosimetry event weights must sum to a positive value.")
    return sum(
        event.statistical_weight_events*event.deposited_energy_eV for event in events
    )/total_weight
end

function weighted_deposited_energy_variance(events::AbstractVector{Weighted_Microdosimetry_Event})
    mean_value = weighted_deposited_energy_mean(events)
    weights = getfield.(events,:statistical_weight_events)
    total_weight = sum(weights)
    return sum(
        event.statistical_weight_events*(event.deposited_energy_eV-mean_value)^2
        for event in events
    )/total_weight
end

function detector_trigger_probability(
    events::AbstractVector{Weighted_Microdosimetry_Event},
    threshold_eV::Real,
)
    threshold = Float64(threshold_eV)
    isfinite(threshold) && threshold ≥ 0.0 || error("Detector threshold must be nonnegative.")
    isempty(events) && return 0.0
    total_weight = sum(event.statistical_weight_events for event in events)
    total_weight > 0.0 || return 0.0
    triggered = sum(
        event.statistical_weight_events for event in events
        if event.deposited_energy_eV ≥ threshold
    )
    return triggered/total_weight
end

function expected_specific_energy_Gy(
    events::AbstractVector{Weighted_Microdosimetry_Event},
    sensitive_mass_kg::Real,
)
    mass = Float64(sensitive_mass_kg)
    isfinite(mass) && mass > 0.0 || error("Sensitive mass must be finite and positive.")
    energy_eV = sum(
        event.statistical_weight_events*event.deposited_energy_eV for event in events
    )
    return energy_eV*1.602176634e-19/mass
end

const WEIGHTED_MICRODOSIMETRY_EVENT_BANK_HDF5_SCHEMA =
    "radiant.hts.weighted_microdosimetry_event_bank/v2"

function _require_sha256_hex(value::AbstractString,label::AbstractString)
    occursin(r"^[0-9a-fA-F]{64}$",value) || error("$(label) must be a SHA-256 digest.")
    return lowercase(String(value))
end

"""
    write_weighted_microdosimetry_event_bank_hdf5(path, events; source_hash, kernel_hash)

Serialize actual sampled parent and secondary direction vectors.  The bank is a weighted
representative sample of a deterministic rate kernel, never an analog-history claim.
"""
function write_weighted_microdosimetry_event_bank_hdf5(
    path::AbstractString,
    events::AbstractVector{Weighted_Microdosimetry_Event};
    source_hash::AbstractString,
    kernel_hash::AbstractString,
)
    isempty(events) && error("Weighted microdosimetry HDF5 bank cannot be empty.")
    source_digest = _require_sha256_hex(source_hash,"Weighted event bank source hash")
    kernel_digest = _require_sha256_hex(kernel_hash,"Weighted event bank kernel hash")
    mkpath(dirname(abspath(path)))
    HDF5.h5open(path,"w") do handle
        HDF5.attributes(handle)["schema_id"] = WEIGHTED_MICRODOSIMETRY_EVENT_BANK_HDF5_SCHEMA
        HDF5.attributes(handle)["representative_not_analog"] = true
        HDF5.attributes(handle)["source_hash"] = source_digest
        HDF5.attributes(handle)["kernel_hash"] = kernel_digest
        count = length(events)
        event_group = HDF5.create_group(handle,"events")
        partition_group = HDF5.create_group(event_group,"partition")
        secondary_group = HDF5.create_group(handle,"secondaries")
        event_group["sample_id"] = getfield.(events,:sample_id)
        event_group["prototype_id"] = getfield.(events,:prototype_id)
        event_group["particle_tag"] = getfield.(events,:particle_tag)
        event_group["process_id"] = getfield.(events,:process_id)
        event_group["material_tag"] = getfield.(events,:material_tag)
        event_group["layer_id"] = getfield.(events,:layer_id)
        event_group["correlation_id"] = getfield.(events,:correlation_id)
        event_group["time_s"] = getfield.(events,:time_s)
        event_group["position_cm"] = reduce(hcat,collect.(getfield.(events,:position_cm)))
        event_group["direction"] = reduce(hcat,collect.(getfield.(events,:direction)))
        event_group["direction_model"] = String.(getfield.(events,:direction_model))
        event_group["incident_energy_eV"] = getfield.(events,:incident_energy_eV)
        event_group["deposited_energy_eV"] = getfield.(events,:deposited_energy_eV)
        event_group["statistical_weight_events"] = getfield.(events,:statistical_weight_events)
        partitions = getfield.(events,:partition)
        for field in fieldnames(Energy_Partition)
            field in (:label,:metadata) && continue
            partition_group[string(field)] = getfield.(partitions,field)
        end
        parent_index = Int64[]
        particle = String[]
        energy = Float64[]
        direction = NTuple{3,Float64}[]
        model = String[]
        correlation = String[]
        for (index,event) in enumerate(events)
            for secondary in event.correlated_secondaries
                push!(parent_index,Int64(index))
                push!(particle,secondary.particle_tag)
                push!(energy,secondary.kinetic_energy_eV)
                push!(direction,secondary.direction)
                push!(model,String(secondary.direction_model))
                push!(correlation,secondary.correlation_id)
            end
        end
        secondary_group["parent_event_index_1based"] = parent_index
        secondary_group["particle_tag"] = particle
        secondary_group["kinetic_energy_eV"] = energy
        secondary_group["direction"] = isempty(direction) ? zeros(Float64,3,0) : reduce(hcat,collect.(direction))
        secondary_group["direction_model"] = model
        secondary_group["correlation_id"] = correlation
        # HDF5's language-neutral on-disk dimension order is event,component;
        # HDF5.jl presents the reversed (component,event) view on read.
        HDF5.attributes(event_group)["direction_axis_order"] = "event,component"
        HDF5.attributes(secondary_group)["direction_axis_order"] = "secondary,component"
        HDF5.attributes(event_group)["event_count"] = count
    end
    return String(path)
end

function synthetic_microdosimetry_kernel_fixture()
    fractions = Event_Energy_Partition_Fractions(Dict(
        :ionization => 0.20,
        :electronic_excitation => 0.10,
        :prompt_lattice_heat => 0.40,
        :nuclear_recoil_handoff => 0.10,
        :defect_storage => 0.05,
        :optical_emission => 0.05,
        :escaping_particle => 0.05,
        :cutoff_handoff => 0.05,
        :unresolved => 0.0,
    ))
    prototype = Microdosimetry_Event_Prototype(
        prototype_id="synthetic-event",
        particle_tag="electron",
        process_id="synthetic-ionization",
        material_tag="YBCO",
        layer_id="REBCO",
        position_cm=[0.0,0.0,0.0],
        direction_model=:isotropic,
        incident_energy_eV=1000.0,
        mean_deposited_energy_eV=500.0,
        event_rate_per_s=100.0,
        straggling=Energy_Straggling_Model(
            :lognormal;
            relative_sigma=0.2,
            maximum_deposited_energy_eV=1000.0,
            model_hash="synthetic-straggling",
            qualification_status=:verification,
        ),
        partition_fractions=fractions,
        correlated_secondaries=[
            Correlated_Secondary(
                "photon",50.0;
                direction_model=:isotropic,
                correlation_id="synthetic-bundle",
            ),
        ],
        correlation_id="synthetic-bundle",
    )
    return Microdosimetry_Kernel(
        [prototype];
        source_artifact_hash="synthetic-source",
        geometry_hash="synthetic-geometry",
        material_state_hash="synthetic-material-state",
        metadata=Dict("classification" => "verification-only"),
    )
end
