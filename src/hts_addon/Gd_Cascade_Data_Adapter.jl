const GD_CASCADE_DATA_STATUSES = (
    :discovered,
    :downloaded,
    :candidate,
    :qualified,
)

struct Gd_Cascade_File_Manifest
    nuclide::String
    residual_nuclide::String
    q_value_eV::Float64
    mean_recoil_energy_eV::Union{Nothing,Float64}
    source_title::String
    source_identifier::String
    source_url::String
    expected_file_sha256::String
    energy_multiplier_to_eV::Float64
    skipped_numeric_columns::Int64
    status::Symbol
    metadata::Dict{String,String}

    function Gd_Cascade_File_Manifest(;
        nuclide::AbstractString,
        residual_nuclide::AbstractString,
        q_value_eV::Real,
        mean_recoil_energy_eV::Union{Nothing,Real}=nothing,
        source_title::AbstractString,
        source_identifier::AbstractString,
        source_url::AbstractString,
        expected_file_sha256::AbstractString="unbound",
        energy_multiplier_to_eV::Real=1.0e6,
        skipped_numeric_columns::Integer=0,
        status::Symbol=:discovered,
        metadata::AbstractDict=Dict{String,String}(),
    )
        startswith(String(nuclide),"Gd-") || error("Gd cascade file requires an isotope tag.")
        for (label,value) in (
            ("residual nuclide",residual_nuclide),("source title",source_title),
            ("source identifier",source_identifier),("source URL",source_url),
            ("expected SHA-256",expected_file_sha256),
        )
            isempty(value) && error("Gd cascade $(label) cannot be empty.")
        end
        q_value = Float64(q_value_eV)
        isfinite(q_value) && q_value > 0.0 || error("Gd capture Q value must be positive.")
        recoil = isnothing(mean_recoil_energy_eV) ? nothing : Float64(mean_recoil_energy_eV)
        if !isnothing(recoil)
            isfinite(recoil) && recoil >= 0.0 || error("Mean recoil energy must be nonnegative.")
        end
        multiplier = Float64(energy_multiplier_to_eV)
        isfinite(multiplier) && multiplier > 0.0 || error(
            "Cascade energy multiplier must be positive.",
        )
        skipped_numeric_columns >= 0 || error("Skipped numeric columns must be nonnegative.")
        status in GD_CASCADE_DATA_STATUSES || error("Unknown Gd cascade data status $(status).")
        status in (:candidate,:qualified) && expected_file_sha256 == "unbound" && error(
            "Candidate or qualified Gd cascade files require a bound SHA-256.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(nuclide),String(residual_nuclide),q_value,recoil,String(source_title),
            String(source_identifier),String(source_url),String(expected_file_sha256),multiplier,
            Int64(skipped_numeric_columns),status,metadata_string,
        )
    end
end

struct Gd_Cascade_Spectrum_Summary
    manifest::Gd_Cascade_File_Manifest
    source_file_sha256::String
    event_count::Int64
    photon_count::Int64
    energy_edges_eV::Vector{Float64}
    photon_yield_per_capture::Vector{Float64}
    photon_yield_standard_uncertainty::Vector{Float64}
    mean_multiplicity::Float64
    multiplicity_standard_deviation::Float64
    mean_gamma_energy_sum_eV::Float64
    gamma_energy_sum_standard_deviation_eV::Float64
    energy_residual_eV::Union{Nothing,Float64}
    line_parse_failures::Int64
    status::Symbol
    metadata::Dict{String,String}
end

function _parse_cascade_numeric_line(line::AbstractString)
    cleaned = strip(first(split(line,'#';limit=2)))
    isempty(cleaned) && return nothing
    tokens = split(replace(replace(cleaned,',' => ' '),';' => ' '))
    values = Float64[]
    for token in tokens
        isempty(token) && continue
        value = try
            parse(Float64,token)
        catch
            return nothing
        end
        push!(values,value)
    end
    return isempty(values) ? nothing : values
end

"""
    summarize_gd_gamma_cascade_file(path, manifest, energy_edges_eV; ...)

Stream a large event file without retaining all cascades. Every numerical input line is interpreted
as one capture event after `skipped_numeric_columns`; remaining values are photon energies. This is
suited to event repositories containing millions of cascades. The parser is deliberately strict:
malformed data lines are counted and fail the candidate/qualified result unless explicitly allowed.
"""
function summarize_gd_gamma_cascade_file(
    path::AbstractString,
    manifest::Gd_Cascade_File_Manifest,
    energy_edges_eV::AbstractVector{<:Real};
    maximum_events::Union{Nothing,Integer}=nothing,
    allow_parse_failures::Bool=false,
    closure_relative_tolerance::Real=5.0e-3,
    closure_absolute_tolerance_eV::Real=5.0e3,
)
    isfile(path) || error("Gd cascade event file does not exist: $(path).")
    actual_hash = _source_file_sha256(path)
    if manifest.expected_file_sha256 != "unbound"
        lowercase(actual_hash) == lowercase(manifest.expected_file_sha256) || error(
            "Gd cascade event-file SHA-256 mismatch.",
        )
    end
    edges = Float64.(energy_edges_eV)
    length(edges) >= 2 && all(isfinite,edges) && all(diff(edges) .> 0.0) || error(
        "Gd cascade energy boundaries must be finite and strictly increasing.",
    )
    edges[1] >= 0.0 || error("Gd cascade energy boundaries cannot be negative.")
    maximum_events_value = isnothing(maximum_events) ? typemax(Int64) : Int64(maximum_events)
    maximum_events_value >= 1 || error("maximum_events must be positive.")
    counts = zeros(Int64,length(edges)-1)
    counts_squared_per_event = zeros(Float64,length(edges)-1)
    event_count = 0
    photon_count = 0
    multiplicity_sum = 0.0
    multiplicity_sum2 = 0.0
    energy_sum = 0.0
    energy_sum2 = 0.0
    parse_failures = 0
    for raw_line in eachline(path)
        values = _parse_cascade_numeric_line(raw_line)
        if isnothing(values)
            stripped = strip(first(split(raw_line,'#';limit=2)))
            isempty(stripped) || (parse_failures += 1)
            continue
        end
        length(values) > manifest.skipped_numeric_columns || begin
            parse_failures += 1
            continue
        end
        energies = values[manifest.skipped_numeric_columns+1:end].*
                   manifest.energy_multiplier_to_eV
        all(value -> isfinite(value) && value > 0.0,energies) || begin
            parse_failures += 1
            continue
        end
        event_histogram = zeros(Int64,length(counts))
        for energy in energies
            energy <= manifest.q_value_eV*(1.0+1.0e-6) || error(
                "A cascade photon exceeds the capture Q value.",
            )
            energy >= edges[1] && energy <= edges[end] || error(
                "A cascade photon lies outside the requested energy structure.",
            )
            group = energy == edges[end] ? length(counts) : searchsortedlast(edges,energy)
            group = min(max(group,1),length(counts))
            event_histogram[group] += 1
        end
        event_count += 1
        multiplicity = length(energies)
        event_energy = sum(energies)
        photon_count += multiplicity
        multiplicity_sum += multiplicity
        multiplicity_sum2 += multiplicity^2
        energy_sum += event_energy
        energy_sum2 += event_energy^2
        counts .+= event_histogram
        counts_squared_per_event .+= Float64.(event_histogram).^2
        event_count >= maximum_events_value && break
    end
    event_count >= 2 || error("Gd cascade summary requires at least two parsed events.")
    parse_failures == 0 || allow_parse_failures || error(
        "Gd cascade file contains $(parse_failures) malformed non-comment lines.",
    )
    mean_multiplicity = multiplicity_sum/event_count
    multiplicity_variance = max(
        0.0,(multiplicity_sum2-event_count*mean_multiplicity^2)/(event_count-1),
    )
    mean_energy = energy_sum/event_count
    energy_variance = max(0.0,(energy_sum2-event_count*mean_energy^2)/(event_count-1))
    yields = Float64.(counts)./event_count
    yield_uncertainty = zeros(Float64,length(yields))
    for group in eachindex(yields)
        second_moment = counts_squared_per_event[group]/event_count
        event_variance = max(0.0,second_moment-yields[group]^2)
        yield_uncertainty[group] = sqrt(event_variance/event_count)
    end
    residual = isnothing(manifest.mean_recoil_energy_eV) ? nothing :
        manifest.q_value_eV-mean_energy-manifest.mean_recoil_energy_eV
    closure_pass = isnothing(residual) || abs(residual) <=
        Float64(closure_absolute_tolerance_eV)+
        Float64(closure_relative_tolerance)*manifest.q_value_eV
    requested_status = manifest.status
    status = if requested_status == :qualified && closure_pass &&
                parse_failures == 0 && actual_hash != "unbound"
        :qualified
    elseif requested_status in (:candidate,:qualified) && closure_pass
        :candidate
    else
        :downloaded
    end
    metadata = Dict(
        "parser" => "one-event-per-line-photon-energies",
        "maximum_events" => isnothing(maximum_events) ? "all" : string(maximum_events),
        "closure_pass" => string(closure_pass),
        "gamma_only" => "true",
        "conversion_and_auger_electrons_included" => "false",
        "source_path" => abspath(path),
    )
    return Gd_Cascade_Spectrum_Summary(
        manifest,actual_hash,Int64(event_count),Int64(photon_count),edges,yields,
        yield_uncertainty,mean_multiplicity,sqrt(multiplicity_variance),mean_energy,
        sqrt(energy_variance),residual,Int64(parse_failures),status,metadata,
    )
end

function build_gd_histogram_volume_source(
    summary::Gd_Cascade_Spectrum_Summary,
    rates::Capture_Rate_Field,
)
    summary.manifest.nuclide == rates.nuclide || error(
        "Gd cascade summary and capture-rate isotope do not match.",
    )
    voxel_count = length(rates.voxel_ids)
    group_count = length(summary.energy_edges_eV)-1
    values = zeros(Float64,voxel_count,group_count,1)
    variance = zeros(Float64,voxel_count,group_count,1)
    for voxel in 1:voxel_count, group in 1:group_count
        volume = rates.voxel_volumes_cm3[voxel]
        rate = rates.capture_rate_per_s[voxel]
        yield_value = summary.photon_yield_per_capture[group]
        values[voxel,group,1] = rate*yield_value/volume
        rate_variance = isnothing(rates.variance_per_s2) ? 0.0 :
            rates.variance_per_s2[voxel]
        yield_variance = summary.photon_yield_standard_uncertainty[group]^2
        variance[voxel,group,1] =
            (yield_value/volume)^2*rate_variance+(rate/volume)^2*yield_variance
    end
    normalization = Source_Normalization(
        basis=:per_second,source_rate_per_s=1.0,symmetry_factor=1.0,
        time_interval_s=(0.0,0.0),time_class=:prompt,
        source_hash=rates.transport_artifact_hash,
        provenance=Dict(
            "capture_cascade_sha256" => summary.source_file_sha256,
            "capture_rate_hash" => rates.transport_artifact_hash,
        ),
    )
    return Anisotropic_Volume_Source(
        Photon(),rates.voxel_ids,rates.voxel_volumes_cm3,summary.energy_edges_eV,
        :isotropic,values,normalization;
        variance=variance,parent_reaction=summary.manifest.nuclide*"(n,gamma)",
        provenance=Dict(
            "cascade_source_identifier" => summary.manifest.source_identifier,
            "cascade_status" => string(summary.status),
            "gamma_only" => "true",
            "capture_rates_are_self_shielded" => get(
                rates.provenance,"capture_rates_are_self_shielded","unbound",
            ),
        ),
    )
end

function gd_cascade_summary_is_production_ready(summary::Gd_Cascade_Spectrum_Summary)
    return summary.status == :qualified && summary.line_parse_failures == 0 &&
           summary.source_file_sha256 == summary.manifest.expected_file_sha256 &&
           get(summary.metadata,"closure_pass","false") == "true"
end

function maurina_zenodo_gd_manifest(
    nuclide::AbstractString;
    q_value_eV::Real,
    residual_nuclide::AbstractString,
    expected_file_sha256::AbstractString="unbound",
    skipped_numeric_columns::Integer=0,
    energy_multiplier_to_eV::Real=1.0e6,
    status::Symbol=:discovered,
)
    nuclide in ("Gd-155","Gd-157") || error(
        "Zenodo 7458654 helper supports Gd-155 or Gd-157.",
    )
    return Gd_Cascade_File_Manifest(
        nuclide=nuclide,residual_nuclide=residual_nuclide,q_value_eV=q_value_eV,
        source_title="Gamma cascades in gadolinium isotopes",
        source_identifier="doi:10.5281/zenodo.7458654",
        source_url="https://doi.org/10.5281/zenodo.7458654",
        expected_file_sha256=expected_file_sha256,
        energy_multiplier_to_eV=energy_multiplier_to_eV,
        skipped_numeric_columns=skipped_numeric_columns,status=status,
        metadata=Dict(
            "event_count_reported" => "10000000 per isotope",
            "model" => "Maurina Hauser-Feshbach cascade sample",
            "internal_conversion_in_source" => "false",
            "use" => "gamma-cascade spectral and multiplicity reference",
        ),
    )
end

function gd_cascade_data_sources()
    return Dict{String,Any}(
        "schema" => "radiant.hts.gd_cascade_sources/v1",
        "gamma_event_repository" => Dict(
            "identifier" => "doi:10.5281/zenodo.7458654",
            "description" => "10^7 gamma cascades for Gd-155 and Gd-157",
            "adapter" => "summarize_gd_gamma_cascade_file",
        ),
        "prompt_gamma_database" => Dict(
            "identifier" => "IAEA-PGAA",
            "url" => "https://nucleus.iaea.org/Pages/pgaa-iaea.aspx",
            "use" => "line energies, partial radiative cross sections, branching and checks",
        ),
        "measured_spectra" => [
            "doi:10.1093/ptep/ptz002",
            "doi:10.1093/ptep/ptaa015",
        ],
        "conversion_electron_option" => Dict(
            "identifier" => "Geant4-NuDEX",
            "use" => "event-level nuclear de-excitation with internal conversion",
            "qualification" => "compare gamma/electron yields to measured and evaluated data",
        ),
        "production_rule" =>
            "Do not label a gamma-only event repository as a complete capture heat source. " *
            "Internal-conversion electrons, atomic relaxation, and recoil remain separate ledgers.",
    )
end
