const GD_CAPTURE_EMISSION_SPECIES = (
    :photon,
    :xray,
    :conversion_electron,
    :auger_electron,
)

"""One evaluated prompt emission line or grouped continuum bin per neutron capture."""
struct Capture_Emission_Line
    species::Symbol
    energy_eV::Float64
    yield_per_capture::Float64
    delay_s::Float64
    correlation_id::String
    transition_id::String
    metadata::Dict{String,String}

    function Capture_Emission_Line(
        species::Symbol,
        energy_eV::Real,
        yield_per_capture::Real;
        delay_s::Real=0.0,
        correlation_id::AbstractString="uncorrelated",
        transition_id::AbstractString="unresolved",
        metadata::AbstractDict=Dict{String,String}(),
    )
        species in GD_CAPTURE_EMISSION_SPECIES || error(
            "Unsupported Gd capture emission species: $(species).",
        )
        energy = Float64(energy_eV)
        yield_value = Float64(yield_per_capture)
        delay = Float64(delay_s)
        isfinite(energy) && energy > 0.0 || error("Capture emission energy must be positive.")
        isfinite(yield_value) && yield_value ≥ 0.0 || error(
            "Capture emission yield must be finite and nonnegative.",
        )
        isfinite(delay) && delay ≥ 0.0 || error("Capture emission delay must be nonnegative.")
        isempty(correlation_id) && error("Capture emission correlation identifier cannot be empty.")
        isempty(transition_id) && error("Capture transition identifier cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            species,energy,yield_value,delay,String(correlation_id),String(transition_id),
            metadata_string,
        )
    end
end

"""
    Gd_Prompt_Capture_Cascade

Isotope-resolved evaluated capture cascade. No physical line list is built into RadiantHTS: the
caller must supply an evaluated or measured cascade with provenance. This prevents natural Gd,
Gd-155, and Gd-157 from being silently represented by one generic capture gamma.
"""
struct Gd_Prompt_Capture_Cascade
    nuclide::String
    residual_nuclide::String
    reaction_id::String
    q_value_eV::Float64
    emissions::Vector{Capture_Emission_Line}
    mean_recoil_energy_eV::Union{Nothing,Float64}
    evaluation_id::String
    evaluation_hash::String
    qualification_status::Symbol
    metadata::Dict{String,String}

    function Gd_Prompt_Capture_Cascade(
        nuclide::AbstractString,
        residual_nuclide::AbstractString,
        q_value_eV::Real,
        emissions::AbstractVector{Capture_Emission_Line};
        reaction_id::AbstractString="n,gamma",
        mean_recoil_energy_eV::Union{Nothing,Real}=nothing,
        evaluation_id::AbstractString,
        evaluation_hash::AbstractString,
        qualification_status::Symbol=:candidate,
        metadata::AbstractDict=Dict{String,String}(),
    )
        startswith(String(nuclide),"Gd-") || error(
            "Gd capture cascades require an isotope-resolved Gd nuclide tag.",
        )
        isempty(residual_nuclide) && error("Residual nuclide cannot be empty.")
        isempty(reaction_id) && error("Capture reaction identifier cannot be empty.")
        q_value = Float64(q_value_eV)
        isfinite(q_value) && q_value > 0.0 || error("Capture Q value must be positive.")
        emission_vector = Capture_Emission_Line[emissions...]
        isempty(emission_vector) && error("Capture cascade requires at least one emission.")
        recoil = isnothing(mean_recoil_energy_eV) ? nothing : Float64(mean_recoil_energy_eV)
        if !isnothing(recoil)
            isfinite(recoil) && recoil ≥ 0.0 || error("Mean capture recoil energy must be nonnegative.")
        end
        isempty(evaluation_id) && error("Capture cascade evaluation identifier cannot be empty.")
        isempty(evaluation_hash) && error("Capture cascade evaluation hash cannot be empty.")
        qualification_status in (:verification,:candidate,:qualified) || error(
            "Capture cascade status must be :verification, :candidate, or :qualified.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(nuclide),String(residual_nuclide),String(reaction_id),q_value,
            emission_vector,recoil,String(evaluation_id),String(evaluation_hash),
            qualification_status,metadata_string,
        )
    end
end

"""Isotope-resolved capture rates already corrected for transport and self-shielding."""
struct Capture_Rate_Field
    nuclide::String
    voxel_ids::Vector{Int64}
    voxel_volumes_cm3::Vector{Float64}
    capture_rate_per_s::Vector{Float64}
    variance_per_s2::Union{Nothing,Vector{Float64}}
    material_tags::Vector{String}
    transport_artifact_hash::String
    provenance::Dict{String,String}

    function Capture_Rate_Field(
        nuclide::AbstractString,
        voxel_ids::AbstractVector{<:Integer},
        voxel_volumes_cm3::AbstractVector{<:Real},
        capture_rate_per_s::AbstractVector{<:Real};
        variance_per_s2::Union{Nothing,AbstractVector{<:Real}}=nothing,
        material_tags::AbstractVector{<:AbstractString}=String[],
        transport_artifact_hash::AbstractString,
        provenance::AbstractDict=Dict{String,String}(),
    )
        ids = Int64.(voxel_ids)
        volumes = Float64.(voxel_volumes_cm3)
        rates = Float64.(capture_rate_per_s)
        count = length(ids)
        count > 0 || error("Capture-rate field cannot be empty.")
        length(unique(ids)) == count || error("Capture-rate voxel identifiers must be unique.")
        length(volumes) == count && length(rates) == count || error(
            "Capture-rate arrays must have matching lengths.",
        )
        all(value -> isfinite(value) && value > 0.0,volumes) || error(
            "Capture-rate voxel volumes must be finite and positive.",
        )
        all(value -> isfinite(value) && value ≥ 0.0,rates) || error(
            "Capture rates must be finite and nonnegative.",
        )
        variance = isnothing(variance_per_s2) ? nothing : Float64.(variance_per_s2)
        if !isnothing(variance)
            length(variance) == count || error("Capture-rate variance length is inconsistent.")
            all(value -> isfinite(value) && value ≥ 0.0,variance) || error(
                "Capture-rate variances must be finite and nonnegative.",
            )
        end
        tags = isempty(material_tags) ? fill("unbound",count) :
            String[string(value) for value in material_tags]
        length(tags) == count || error("Capture-rate material tags must match voxel count.")
        isempty(transport_artifact_hash) && error("Capture-rate transport hash cannot be empty.")
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            String(nuclide),ids,volumes,rates,variance,tags,
            String(transport_artifact_hash),provenance_string,
        )
    end
end

struct Gd_Capture_Recoil_Source
    nuclide::String
    residual_nuclide::String
    voxel_ids::Vector{Int64}
    recoil_energy_eV::Vector{Float64}
    event_rate_per_s::Vector{Float64}
    transport_owner::String
    evaluation_hash::String
end

struct Gd_Capture_Source_Bundle
    nuclide::String
    sources::Dict{Symbol,Anisotropic_Volume_Source}
    recoil_source::Union{Nothing,Gd_Capture_Recoil_Source}
    prompt_emitted_energy_eV_per_capture::Float64
    recoil_energy_eV_per_capture::Union{Nothing,Float64}
    energy_residual_eV_per_capture::Union{Nothing,Float64}
    correlation_records::Vector{Dict{String,String}}
    qualification_status::Symbol
    provenance::Dict{String,String}
end

function _capture_energy_bin(edges::Vector{Float64},energy::Float64)
    energy < edges[1] && error("Capture emission lies below the source group structure.")
    energy > edges[end] && error("Capture emission lies above the source group structure.")
    energy == edges[end] && return length(edges)-1
    index = searchsortedlast(edges,energy)
    index == 0 && error("Capture emission could not be assigned to an energy group.")
    return min(index,length(edges)-1)
end

function _capture_source_for_species(
    species::Symbol,
    lines::Vector{Capture_Emission_Line},
    cascade::Gd_Prompt_Capture_Cascade,
    rates::Capture_Rate_Field,
    energy_edges_eV::AbstractVector{<:Real},
)
    edges = Float64.(energy_edges_eV)
    length(edges) ≥ 2 && all(diff(edges) .> 0.0) || error(
        "Capture-source energy boundaries must be strictly increasing.",
    )
    values = zeros(Float64,length(rates.voxel_ids),length(edges)-1,1)
    variance = isnothing(rates.variance_per_s2) ? nothing : zeros(Float64,size(values))
    for line in lines
        group = _capture_energy_bin(edges,line.energy_eV)
        for voxel in eachindex(rates.voxel_ids)
            source_density = rates.capture_rate_per_s[voxel]*line.yield_per_capture/
                             rates.voxel_volumes_cm3[voxel]
            values[voxel,group,1] += source_density
            if !isnothing(variance)
                variance[voxel,group,1] += rates.variance_per_s2[voxel]*
                                           line.yield_per_capture^2/
                                           rates.voxel_volumes_cm3[voxel]^2
            end
        end
    end
    particle = species in (:photon,:xray) ? Photon() : Electron()
    normalization = Source_Normalization(
        basis=:per_second,
        source_rate_per_s=1.0,
        symmetry_factor=1.0,
        time_interval_s=(0.0,0.0),
        time_class=:prompt,
        source_hash=rates.transport_artifact_hash,
        provenance=Dict(
            "producer" => "Gd_Prompt_Capture_Cascade",
            "cascade_evaluation_hash" => cascade.evaluation_hash,
        ),
    )
    return Anisotropic_Volume_Source(
        particle,rates.voxel_ids,rates.voxel_volumes_cm3,edges,:isotropic,values,
        normalization;
        variance=variance,
        parent_reaction=string(cascade.nuclide,"(",cascade.reaction_id,")"),
        provenance=Dict(
            "emission_species" => string(species),
            "nuclide" => cascade.nuclide,
            "residual_nuclide" => cascade.residual_nuclide,
            "evaluation_id" => cascade.evaluation_id,
            "evaluation_hash" => cascade.evaluation_hash,
            "transport_artifact_hash" => rates.transport_artifact_hash,
            "capture_rates_are_self_shielded" => get(
                rates.provenance,"capture_rates_are_self_shielded","unbound",
            ),
        ),
    )
end

"""
    build_gd_capture_source_bundle(cascade, rates; photon_energy_edges_eV,
                                   electron_energy_edges_eV)

Convert isotope-resolved capture rates to separate prompt photon, x-ray, conversion-electron,
Auger-electron, and recoil products. Capture rates—not an infinitely dilute cross section—are the
input, so film-scale and repeated-tape self-shielding remain owned by the upstream neutral solve.

Correlation identifiers are retained in the bundle for event-level microdosimetry. The scalar
volume sources are expectation values and do not preserve event-by-event cascade correlations.
"""
function build_gd_capture_source_bundle(
    cascade::Gd_Prompt_Capture_Cascade,
    rates::Capture_Rate_Field;
    photon_energy_edges_eV::AbstractVector{<:Real},
    electron_energy_edges_eV::AbstractVector{<:Real},
    energy_closure_rtol::Real=5.0e-3,
    energy_closure_atol_eV::Real=1.0e3,
)
    cascade.nuclide == rates.nuclide || error(
        "Capture cascade nuclide and capture-rate nuclide do not match.",
    )
    grouped = Dict{Symbol,Vector{Capture_Emission_Line}}(
        species => Capture_Emission_Line[] for species in GD_CAPTURE_EMISSION_SPECIES
    )
    for line in cascade.emissions
        push!(grouped[line.species],line)
    end
    sources = Dict{Symbol,Anisotropic_Volume_Source}()
    for species in GD_CAPTURE_EMISSION_SPECIES
        isempty(grouped[species]) && continue
        edges = species in (:photon,:xray) ? photon_energy_edges_eV :
                                             electron_energy_edges_eV
        sources[species] = _capture_source_for_species(
            species,grouped[species],cascade,rates,edges,
        )
    end

    emitted_energy = sum(
        line.energy_eV*line.yield_per_capture for line in cascade.emissions
    )
    recoil_source = nothing
    energy_residual = nothing
    if !isnothing(cascade.mean_recoil_energy_eV)
        recoil_energy = cascade.mean_recoil_energy_eV
        recoil_source = Gd_Capture_Recoil_Source(
            cascade.nuclide,cascade.residual_nuclide,copy(rates.voxel_ids),
            fill(recoil_energy,length(rates.voxel_ids)),copy(rates.capture_rate_per_s),
            "SPECTRA-PKA/Geant4/BCA",cascade.evaluation_hash,
        )
        energy_residual = cascade.q_value_eV-emitted_energy-recoil_energy
        isapprox(
            emitted_energy+recoil_energy,cascade.q_value_eV;
            rtol=energy_closure_rtol,atol=energy_closure_atol_eV,
        ) || error(
            "Evaluated Gd capture cascade does not close its declared Q value. " *
            "Do not route the residual to heat implicitly.",
        )
    end

    correlation_records = Dict{String,String}[]
    for line in cascade.emissions
        push!(correlation_records,Dict(
            "correlation_id" => line.correlation_id,
            "transition_id" => line.transition_id,
            "species" => string(line.species),
            "energy_eV" => string(line.energy_eV),
            "yield_per_capture" => string(line.yield_per_capture),
        ))
    end
    provenance = Dict(
        "schema" => "radiant.hts.gd_capture_source_bundle/v1",
        "evaluation_id" => cascade.evaluation_id,
        "evaluation_hash" => cascade.evaluation_hash,
        "transport_artifact_hash" => rates.transport_artifact_hash,
        "correlations_preserved_in_scalar_sources" => "false",
        "recoil_status" => isnothing(recoil_source) ? "unresolved" : "evaluated-mean",
    )
    return Gd_Capture_Source_Bundle(
        cascade.nuclide,sources,recoil_source,emitted_energy,
        cascade.mean_recoil_energy_eV,energy_residual,correlation_records,
        cascade.qualification_status,provenance,
    )
end

"""Synthetic cascade used only to verify source construction and energy closure."""
function synthetic_gd157_capture_fixture()
    emissions = Capture_Emission_Line[
        Capture_Emission_Line(
            :photon,6.0e6,1.0;
            correlation_id="fixture-cascade-1",transition_id="gamma-primary",
        ),
        Capture_Emission_Line(
            :photon,1.8e6,1.0;
            correlation_id="fixture-cascade-1",transition_id="gamma-secondary",
        ),
        Capture_Emission_Line(
            :conversion_electron,1.5e5,0.5;
            correlation_id="fixture-cascade-1",transition_id="conversion-electron",
        ),
        Capture_Emission_Line(
            :auger_electron,2.0e4,0.5;
            correlation_id="fixture-cascade-1",transition_id="auger-electron",
        ),
        Capture_Emission_Line(
            :xray,4.0e4,0.5;
            correlation_id="fixture-cascade-1",transition_id="characteristic-xray",
        ),
    ]
    emitted = sum(line.energy_eV*line.yield_per_capture for line in emissions)
    recoil = 100.0
    return Gd_Prompt_Capture_Cascade(
        "Gd-157","Gd-158",emitted+recoil,emissions;
        mean_recoil_energy_eV=recoil,
        evaluation_id="synthetic-verification-only",
        evaluation_hash="synthetic-gd157-cascade",
        qualification_status=:verification,
        metadata=Dict("classification" => "verification-only"),
    )
end
