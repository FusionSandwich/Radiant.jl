const DELAYED_EMISSION_SPECIES = (
    :photon,
    :electron,
    :positron,
    :alpha,
    :recoil,
    :neutrino,
)

const DELAYED_DECAY_DATA_STATUSES = (:verification,:candidate,:qualified)

"""One monoenergetic line or representative continuous-spectrum bin per parent decay."""
struct Delayed_Emission_Bin
    species::Symbol
    representative_energy_eV::Float64
    yield_per_decay::Float64
    correlation_id::String
    transition_id::String
    metadata::Dict{String,String}

    function Delayed_Emission_Bin(
        species::Symbol,
        representative_energy_eV::Real,
        yield_per_decay::Real;
        correlation_id::AbstractString="uncorrelated",
        transition_id::AbstractString="unresolved",
        metadata::AbstractDict=Dict{String,String}(),
    )
        species in DELAYED_EMISSION_SPECIES || error(
            "Unsupported delayed-emission species: $(species).",
        )
        energy = Float64(representative_energy_eV)
        yield_value = Float64(yield_per_decay)
        isfinite(energy) && energy >= 0.0 || error(
            "Delayed-emission representative energy must be finite and nonnegative.",
        )
        isfinite(yield_value) && yield_value >= 0.0 || error(
            "Delayed-emission yield must be finite and nonnegative.",
        )
        yield_value > 0.0 && energy == 0.0 && species != :recoil && error(
            "A non-recoil delayed emission with positive yield requires positive energy.",
        )
        isempty(correlation_id) && error("Delayed-emission correlation ID cannot be empty.")
        isempty(transition_id) && error("Delayed-emission transition ID cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            species,energy,yield_value,String(correlation_id),String(transition_id),
            metadata_string,
        )
    end
end

"""
    Delayed_Decay_Scheme

Parent-resolved mean decay-energy ledger. Neutrino energy, recoil energy, and non-Radiant charged
products remain explicit. `unresolved_energy_eV_per_decay` is allowed for candidate data but blocks
production readiness. No unresolved energy is implicitly assigned to heat.
"""
struct Delayed_Decay_Scheme
    parent_nuclide::String
    daughter_nuclide::String
    decay_mode::String
    q_value_eV::Float64
    emissions::Vector{Delayed_Emission_Bin}
    unresolved_energy_eV_per_decay::Float64
    evaluation_id::String
    evaluation_hash::String
    status::Symbol
    provenance::Dict{String,String}

    function Delayed_Decay_Scheme(
        parent_nuclide::AbstractString,
        daughter_nuclide::AbstractString,
        decay_mode::AbstractString,
        q_value_eV::Real,
        emissions::AbstractVector{Delayed_Emission_Bin};
        unresolved_energy_eV_per_decay::Real=0.0,
        evaluation_id::AbstractString,
        evaluation_hash::AbstractString,
        status::Symbol=:candidate,
        provenance::AbstractDict=Dict{String,String}(),
        closure_rtol::Real=5.0e-3,
        closure_atol_eV::Real=1.0,
    )
        for (label,value) in (
            ("parent nuclide",parent_nuclide),("daughter nuclide",daughter_nuclide),
            ("decay mode",decay_mode),("evaluation ID",evaluation_id),
            ("evaluation hash",evaluation_hash),
        )
            isempty(value) && error("Delayed-decay $(label) cannot be empty.")
        end
        q_value = Float64(q_value_eV)
        unresolved = Float64(unresolved_energy_eV_per_decay)
        isfinite(q_value) && q_value > 0.0 || error("Decay Q value must be positive.")
        isfinite(unresolved) && unresolved >= 0.0 || error(
            "Unresolved decay energy must be finite and nonnegative.",
        )
        emission_vector = Delayed_Emission_Bin[emissions...]
        isempty(emission_vector) && error("Delayed-decay scheme requires at least one emission.")
        status in DELAYED_DECAY_DATA_STATUSES || error(
            "Unknown delayed-decay data status: $(status).",
        )
        represented = sum(
            emission.representative_energy_eV*emission.yield_per_decay
            for emission in emission_vector
        )+unresolved
        isapprox(represented,q_value;rtol=closure_rtol,atol=closure_atol_eV) || error(
            "Delayed-decay energy ledger does not close its Q value: represented=$(represented), Q=$(q_value).",
        )
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            String(parent_nuclide),String(daughter_nuclide),String(decay_mode),q_value,
            emission_vector,unresolved,String(evaluation_id),String(evaluation_hash),status,
            provenance_string,
        )
    end
end

"""Activity and uncertainty for one parent nuclide on explicit activation voxels."""
struct Delayed_Activity_Field
    parent_nuclide::String
    voxel_ids::Vector{Int64}
    voxel_volumes_cm3::Vector{Float64}
    activity_Bq::Vector{Float64}
    activity_standard_uncertainty_Bq::Vector{Float64}
    material_tags::Vector{String}
    cooling_interval_s::Tuple{Float64,Float64}
    inventory_hash::String
    activation_artifact_hash::String
    provenance::Dict{String,String}

    function Delayed_Activity_Field(
        parent_nuclide::AbstractString,
        voxel_ids::AbstractVector{<:Integer},
        voxel_volumes_cm3::AbstractVector{<:Real},
        activity_Bq::AbstractVector{<:Real};
        activity_standard_uncertainty_Bq::AbstractVector{<:Real}=zeros(length(activity_Bq)),
        material_tags::AbstractVector{<:AbstractString}=fill("unbound",length(activity_Bq)),
        cooling_interval_s::Tuple{<:Real,<:Real},
        inventory_hash::AbstractString,
        activation_artifact_hash::AbstractString,
        provenance::AbstractDict=Dict{String,String}(),
    )
        isempty(parent_nuclide) && error("Delayed-activity parent nuclide cannot be empty.")
        ids = Int64.(voxel_ids)
        volumes = Float64.(voxel_volumes_cm3)
        activity = Float64.(activity_Bq)
        uncertainty = Float64.(activity_standard_uncertainty_Bq)
        tags = String[string(value) for value in material_tags]
        count = length(ids)
        count > 0 || error("Delayed-activity field cannot be empty.")
        length(unique(ids)) == count || error("Delayed-activity voxel IDs must be unique.")
        length(volumes) == count && length(activity) == count &&
            length(uncertainty) == count && length(tags) == count || error(
            "Delayed-activity arrays must have matching lengths.",
        )
        all(value -> isfinite(value) && value > 0.0,volumes) || error(
            "Delayed-activity voxel volumes must be finite and positive.",
        )
        all(value -> isfinite(value) && value >= 0.0,activity) || error(
            "Delayed activities must be finite and nonnegative.",
        )
        all(value -> isfinite(value) && value >= 0.0,uncertainty) || error(
            "Delayed-activity uncertainties must be finite and nonnegative.",
        )
        all(value -> !isempty(value),tags) || error(
            "Delayed-activity material tags cannot be empty.",
        )
        t0 = Float64(cooling_interval_s[1])
        t1 = Float64(cooling_interval_s[2])
        isfinite(t0) && isfinite(t1) && t0 >= 0.0 && t1 >= t0 || error(
            "Cooling interval must be finite, nonnegative, and ordered.",
        )
        isempty(inventory_hash) && error("Delayed inventory hash cannot be empty.")
        isempty(activation_artifact_hash) && error(
            "Activation artifact hash cannot be empty.",
        )
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            String(parent_nuclide),ids,volumes,activity,uncertainty,tags,(t0,t1),
            String(inventory_hash),String(activation_artifact_hash),provenance_string,
        )
    end
end

struct Delayed_NonEM_Handoff
    species::Symbol
    parent_nuclide::String
    daughter_nuclide::String
    representative_energy_eV::Float64
    voxel_ids::Vector{Int64}
    event_rate_per_s::Vector{Float64}
    event_rate_standard_uncertainty_per_s::Vector{Float64}
    transport_owner::String
    correlation_id::String
    transition_id::String
    evaluation_hash::String
    inventory_hash::String
end

struct Activation_Delayed_Source_Bundle
    schema::String
    parent_nuclide::String
    daughter_nuclide::String
    decay_mode::String
    em_sources::Dict{Symbol,Anisotropic_Volume_Source}
    non_em_handoffs::Vector{Delayed_NonEM_Handoff}
    total_activity_Bq::Float64
    q_value_energy_rate_eV_per_s::Float64
    em_source_energy_rate_eV_per_s::Float64
    non_em_energy_rate_eV_per_s::Float64
    neutrino_energy_rate_eV_per_s::Float64
    unresolved_energy_rate_eV_per_s::Float64
    energy_balance_residual_eV_per_s::Float64
    evaluation_hash::String
    inventory_hash::String
    status::Symbol
    provenance::Dict{String,String}
end

function _delayed_particle(species::Symbol)
    species == :photon && return Photon()
    species == :electron && return Electron()
    species == :positron && return Positron()
    error("Species $(species) is not transported by the delayed EM-source builder.")
end

function _delayed_energy_group(edges::Vector{Float64},energy_eV::Float64)
    energy_eV < edges[1] && error("Delayed emission lies below the source energy grid.")
    energy_eV > edges[end] && error("Delayed emission lies above the source energy grid.")
    energy_eV == edges[end] && return length(edges)-1
    index = searchsortedlast(edges,energy_eV)
    1 <= index <= length(edges)-1 || error(
        "Delayed emission could not be assigned to an energy group.",
    )
    return index
end

function _delayed_em_source(
    species::Symbol,
    emissions::Vector{Delayed_Emission_Bin},
    scheme::Delayed_Decay_Scheme,
    activity::Delayed_Activity_Field,
    energy_edges_eV::AbstractVector{<:Real},
)
    edges = Float64.(energy_edges_eV)
    length(edges) >= 2 && all(isfinite,edges) && all(edges .>= 0.0) &&
        all(diff(edges) .> 0.0) || error(
        "Delayed-source energy boundaries must be finite and increasing.",
    )
    values = zeros(Float64,length(activity.voxel_ids),length(edges)-1,1)
    variance = zeros(Float64,size(values))
    yield_by_group = zeros(Float64,length(edges)-1)
    for emission in emissions
        group = _delayed_energy_group(edges,emission.representative_energy_eV)
        yield_by_group[group] += emission.yield_per_decay
    end
    for voxel in eachindex(activity.voxel_ids), group in eachindex(yield_by_group)
        source_density = activity.activity_Bq[voxel]*yield_by_group[group]/
                         activity.voxel_volumes_cm3[voxel]
        sigma_density = activity.activity_standard_uncertainty_Bq[voxel]*
                        yield_by_group[group]/activity.voxel_volumes_cm3[voxel]
        values[voxel,group,1] = source_density
        variance[voxel,group,1] = sigma_density^2
    end
    source_hash = bytes2hex(SHA.sha256(codeunits(string(
        activity.inventory_hash,"|",scheme.evaluation_hash,"|",species,
    ))))
    normalization = Source_Normalization(
        basis=:per_second,
        source_rate_per_s=1.0,
        symmetry_factor=1.0,
        time_interval_s=activity.cooling_interval_s,
        time_class=:delayed,
        source_hash=source_hash,
        provenance=Dict(
            "inventory_hash" => activity.inventory_hash,
            "activation_artifact_hash" => activity.activation_artifact_hash,
            "decay_evaluation_hash" => scheme.evaluation_hash,
        ),
    )
    return Anisotropic_Volume_Source(
        _delayed_particle(species),activity.voxel_ids,activity.voxel_volumes_cm3,
        edges,:isotropic,values,normalization;
        variance=variance,
        parent_reaction=string(scheme.parent_nuclide," ",scheme.decay_mode),
        provenance=Dict(
            "schema" => "radiant.hts.activation_delayed_em_source/v1",
            "parent_nuclide" => scheme.parent_nuclide,
            "daughter_nuclide" => scheme.daughter_nuclide,
            "decay_mode" => scheme.decay_mode,
            "emission_species" => string(species),
            "activity_unit" => "Bq=decays/s",
            "values_unit" => "particles/(cm3*s) integrated over energy group",
            "cooling_interval_s" => string(
                activity.cooling_interval_s[1],",",activity.cooling_interval_s[2],
            ),
            "source_energy_is_not_deposited_heat" => "true",
            "inventory_hash" => activity.inventory_hash,
            "activation_artifact_hash" => activity.activation_artifact_hash,
            "evaluation_hash" => scheme.evaluation_hash,
        ),
    )
end

"""
    build_activation_delayed_source_bundle(scheme, activity; energy_grids)

Create parent-resolved delayed photon/electron/positron volume sources and explicit alpha, recoil,
and neutrino handoffs. `energy_grids` maps each transported EM species to ascending group edges in
eV. The activity field is already a physical decay rate, so generated sources use `:per_second`
normalization with unit source-rate scaling.
"""
function build_activation_delayed_source_bundle(
    scheme::Delayed_Decay_Scheme,
    activity::Delayed_Activity_Field;
    energy_grids::AbstractDict,
    energy_closure_rtol::Real=1.0e-10,
    energy_closure_atol_eV_per_s::Real=1.0e-6,
)
    scheme.parent_nuclide == activity.parent_nuclide || error(
        "Delayed-decay scheme and activity-field parent nuclides do not match.",
    )
    grouped = Dict(
        species => Delayed_Emission_Bin[] for species in DELAYED_EMISSION_SPECIES
    )
    for emission in scheme.emissions
        push!(grouped[emission.species],emission)
    end
    em_sources = Dict{Symbol,Anisotropic_Volume_Source}()
    for species in (:photon,:electron,:positron)
        isempty(grouped[species]) && continue
        haskey(energy_grids,species) || error(
            "Delayed source requires an energy grid for species $(species).",
        )
        em_sources[species] = _delayed_em_source(
            species,grouped[species],scheme,activity,energy_grids[species],
        )
    end

    handoffs = Delayed_NonEM_Handoff[]
    for species in (:alpha,:recoil,:neutrino)
        owner = species == :alpha ? "BCA/Geant4/ion-transport" :
                (species == :recoil ? "SPECTRA-PKA/BCA/MD" : "escape-energy-ledger")
        for emission in grouped[species]
            rates = activity.activity_Bq .* emission.yield_per_decay
            uncertainty = activity.activity_standard_uncertainty_Bq .* emission.yield_per_decay
            push!(handoffs,Delayed_NonEM_Handoff(
                species,scheme.parent_nuclide,scheme.daughter_nuclide,
                emission.representative_energy_eV,copy(activity.voxel_ids),rates,uncertainty,
                owner,emission.correlation_id,emission.transition_id,scheme.evaluation_hash,
                activity.inventory_hash,
            ))
        end
    end

    total_activity = sum(activity.activity_Bq)
    q_rate = total_activity*scheme.q_value_eV
    em_rate = total_activity*sum(Float64[
        emission.representative_energy_eV*emission.yield_per_decay
        for emission in scheme.emissions if emission.species in (:photon,:electron,:positron)
    ])
    non_em_rate = total_activity*sum(Float64[
        emission.representative_energy_eV*emission.yield_per_decay
        for emission in scheme.emissions if emission.species in (:alpha,:recoil)
    ])
    neutrino_rate = total_activity*sum(Float64[
        emission.representative_energy_eV*emission.yield_per_decay
        for emission in scheme.emissions if emission.species == :neutrino
    ])
    unresolved_rate = total_activity*scheme.unresolved_energy_eV_per_decay
    residual = q_rate-em_rate-non_em_rate-neutrino_rate-unresolved_rate
    isapprox(
        em_rate+non_em_rate+neutrino_rate+unresolved_rate,q_rate;
        rtol=energy_closure_rtol,atol=energy_closure_atol_eV_per_s,
    ) || error("Delayed-source energy-rate ledger does not close the decay Q value.")

    production_ready = scheme.status == :qualified &&
        scheme.unresolved_energy_eV_per_decay == 0.0 &&
        length(scheme.evaluation_hash) == 64 &&
        length(activity.inventory_hash) == 64 &&
        length(activity.activation_artifact_hash) == 64
    provenance = Dict{String,String}(
        "decay_evaluation_id" => scheme.evaluation_id,
        "decay_evaluation_hash" => scheme.evaluation_hash,
        "inventory_hash" => activity.inventory_hash,
        "activation_artifact_hash" => activity.activation_artifact_hash,
        "activity_unit" => "Bq=decays/s",
        "prompt_delayed_separation" => "delayed-only",
        "source_energy_is_not_deposited_heat" => "true",
        "neutrino_energy_is_local_heat" => "false",
        "production_ready" => string(production_ready),
    )
    return Activation_Delayed_Source_Bundle(
        "radiant.hts.activation_delayed_source_bundle/v1",scheme.parent_nuclide,
        scheme.daughter_nuclide,scheme.decay_mode,em_sources,handoffs,total_activity,q_rate,
        em_rate,non_em_rate,neutrino_rate,unresolved_rate,residual,scheme.evaluation_hash,
        activity.inventory_hash,scheme.status,provenance,
    )
end

function activation_delayed_bundle_is_production_ready(
    bundle::Activation_Delayed_Source_Bundle,
)
    return bundle.status == :qualified &&
           bundle.unresolved_energy_rate_eV_per_s == 0.0 &&
           length(bundle.evaluation_hash) == 64 && length(bundle.inventory_hash) == 64 &&
           get(bundle.provenance,"production_ready","false") == "true"
end

function activation_delayed_source_receipt(bundle::Activation_Delayed_Source_Bundle)
    return Dict{String,Any}(
        "schema" => bundle.schema,
        "parent_nuclide" => bundle.parent_nuclide,
        "daughter_nuclide" => bundle.daughter_nuclide,
        "decay_mode" => bundle.decay_mode,
        "transported_em_species" => sort(string.(collect(keys(bundle.em_sources)))),
        "non_em_handoff_species" => sort(string.(unique(
            getfield.(bundle.non_em_handoffs,:species),
        ))),
        "total_activity_Bq" => bundle.total_activity_Bq,
        "q_value_energy_rate_eV_per_s" => bundle.q_value_energy_rate_eV_per_s,
        "em_source_energy_rate_eV_per_s" => bundle.em_source_energy_rate_eV_per_s,
        "non_em_energy_rate_eV_per_s" => bundle.non_em_energy_rate_eV_per_s,
        "neutrino_energy_rate_eV_per_s" => bundle.neutrino_energy_rate_eV_per_s,
        "unresolved_energy_rate_eV_per_s" => bundle.unresolved_energy_rate_eV_per_s,
        "energy_balance_residual_eV_per_s" => bundle.energy_balance_residual_eV_per_s,
        "source_energy_is_not_deposited_heat" => true,
        "neutrino_energy_is_local_heat" => false,
        "evaluation_hash" => bundle.evaluation_hash,
        "inventory_hash" => bundle.inventory_hash,
        "status" => string(bundle.status),
        "production_ready" => activation_delayed_bundle_is_production_ready(bundle),
        "provenance" => copy(bundle.provenance),
    )
end

"""Synthetic beta-minus decay used only for software verification."""
function synthetic_activation_decay_fixture()
    emissions = Delayed_Emission_Bin[
        Delayed_Emission_Bin(
            :photon,1.0e6,1.0;
            correlation_id="synthetic-decay",transition_id="gamma-1",
        ),
        Delayed_Emission_Bin(
            :electron,5.0e5,1.0;
            correlation_id="synthetic-decay",transition_id="beta-bin-mean",
        ),
        Delayed_Emission_Bin(
            :recoil,1.0e5,1.0;
            correlation_id="synthetic-decay",transition_id="daughter-recoil",
        ),
        Delayed_Emission_Bin(
            :neutrino,1.4e6,1.0;
            correlation_id="synthetic-decay",transition_id="antineutrino-bin-mean",
        ),
    ]
    return Delayed_Decay_Scheme(
        "X-100","Y-100","beta-minus",3.0e6,emissions;
        evaluation_id="synthetic-decay-verification",
        evaluation_hash="synthetic-decay-data",
        status=:verification,
        provenance=Dict("classification" => "software-verification"),
    )
end
