const PROCESS_RESPONSE_QUANTITIES = (
    "energy-deposition",
    "charge-deposition",
    "absorption",
    "total",
    "stopping-power",
    "momentum-transfer",
)

"""Spatial response fields resolved by the interaction/type/outgoing-particle kernel key."""
struct Process_Resolved_Score
    particle_tag::String
    quantity::String
    channels::Dict{String,Array{Float64,3}}
    total::Array{Float64,3}
    normalization_basis::Symbol
    units::String
    provenance::Dict{String,String}

    function Process_Resolved_Score(
        particle_tag::AbstractString,
        quantity::AbstractString,
        channels::AbstractDict,
        total::AbstractArray{<:Real,3};
        normalization_basis::Symbol=:transport_source_basis,
        units::AbstractString="unspecified",
        provenance::AbstractDict=Dict{String,String}(),
    )
        isempty(particle_tag) && error("Process-resolved score requires a particle tag.")
        quantity_string = String(quantity)
        quantity_string in PROCESS_RESPONSE_QUANTITIES || error(
            "Unsupported process response quantity: $(quantity_string).",
        )
        channel_dictionary = Dict{String,Array{Float64,3}}()
        total_array = Float64.(total)
        for (key,value) in channels
            field = Float64.(value)
            size(field) == size(total_array) || error(
                "Every process response field must match the total field shape.",
            )
            any(x -> !isfinite(x),field) && error("Process response fields must be finite.")
            channel_dictionary[string(key)] = field
        end
        isempty(channel_dictionary) && error("At least one process response channel is required.")
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            String(particle_tag),quantity_string,channel_dictionary,total_array,
            normalization_basis,String(units),provenance_string,
        )
    end
end

function get_process_score(this::Process_Resolved_Score,process_key::AbstractString)
    haskey(this.channels,String(process_key)) || error(
        "Process response channel is unavailable: $(process_key).",
    )
    return this.channels[String(process_key)]
end

get_process_keys(this::Process_Resolved_Score) = sort(collect(keys(this.channels)))

function _process_particle_index(cross_sections::Cross_Sections,particle::Particle)
    index = findfirst(
        candidate -> get_tag(candidate) == get_tag(particle),
        get_particles(cross_sections),
    )
    isnothing(index) && error("Cross sections do not contain particle $(get_tag(particle)).")
    return index
end

function _process_channel_keys(
    cross_sections::Cross_Sections,
    particle::Particle,
    quantity::String,
)
    ismissing(cross_sections.multigroup_cross_sections) && error(
        "Multigroup cross sections must be built before process scoring.",
    )
    particle_index = _process_particle_index(cross_sections,particle)
    keys_found = Set{String}()
    prefix = string(quantity,"|")
    for material_index in 1:get_number_of_materials(cross_sections)
        mcs = cross_sections.multigroup_cross_sections[particle_index,material_index]
        for key in get_response_channel_keys(mcs)
            startswith(key,prefix) && push!(keys_found,key)
        end
    end
    isempty(keys_found) && error(
        "No native $(quantity) channels are stored for particle $(get_tag(particle)).",
    )
    return sort(collect(keys_found))
end

function _score_output_units(quantity::String,mass_normalized::Bool)
    if quantity == "energy-deposition"
        return mass_normalized ? "MeV/g per transport source basis" :
                                 "MeV/cm3 per transport source basis"
    elseif quantity == "charge-deposition"
        return mass_normalized ? "elementary-charge/g per transport source basis" :
                                 "elementary-charge/cm3 per transport source basis"
    elseif quantity in ("absorption","total")
        return "response/cm3 per transport source basis"
    elseif quantity in ("stopping-power","momentum-transfer")
        return "model response per transport source basis"
    end
    return "unspecified"
end

"""
    score_process_responses(cross_sections, geometry, solvers, sources, flux, particle;
                            quantity="energy-deposition", mass_normalized=true,
                            physical_normalization=nothing)

Fold native per-interaction multigroup response channels with Radiant scalar flux. The raw kernel
key is retained, so photoelectric, Compton, pair-production, relaxation, stopping, and other
contributions are not collapsed into an assumed cryogenic heat fraction.

When `physical_normalization` is supplied, its source-rate and symmetry factor are applied after
the transport-source normalization. Prompt and delayed calculations therefore remain separate
unless the caller explicitly combines independently identified score objects.
"""
function score_process_responses(
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solvers::Solvers,
    sources::Fixed_Sources,
    flux::Flux,
    particle::Particle;
    quantity::AbstractString="energy-deposition",
    mass_normalized::Bool=true,
    physical_normalization::Union{Nothing,Source_Normalization}=nothing,
)
    quantity_string = String(quantity)
    quantity_string in PROCESS_RESPONSE_QUANTITIES || error(
        "Unsupported process response quantity: $(quantity_string).",
    )
    keys_full = _process_channel_keys(cross_sections,particle,quantity_string)
    particle_index = _process_particle_index(cross_sections,particle)
    method = get_method(solvers,particle)
    _,is_csd = get_solver_type(method)
    number_of_groups = get_number_of_groups(cross_sections,particle)
    voxel_counts = get_number_of_voxels(geometry)
    material_map = get_material_per_voxel(geometry)
    densities = get_densities(cross_sections)
    scalar_flux = get_flux(flux,particle)
    cutoff_flux = is_csd ? get_flux_cutoff(flux,particle) : nothing
    legacy_normalization = get_normalization_factor(sources)
    isfinite(legacy_normalization) && legacy_normalization > 0.0 || error(
        "Transport source normalization must be finite and positive.",
    )

    channels = Dict{String,Array{Float64,3}}()
    for full_key in keys_full
        process_key = split(full_key,"|";limit=2)[2]
        field = zeros(Float64,voxel_counts[1],voxel_counts[2],voxel_counts[3])
        for ix in 1:voxel_counts[1], iy in 1:voxel_counts[2], iz in 1:voxel_counts[3]
            material_index = material_map[ix,iy,iz]
            mcs = cross_sections.multigroup_cross_sections[particle_index,material_index]
            coefficients = has_response_channel(mcs,full_key) ?
                get_response_channel(mcs,full_key) : zeros(Float64,number_of_groups)
            length(coefficients) in (number_of_groups,number_of_groups+1) || error(
                "Native response channel $(full_key) has an invalid group dimension.",
            )
            value = 0.0
            for group in 1:number_of_groups
                value += coefficients[group]*scalar_flux[group,1,1,ix,iy,iz]
            end
            if length(coefficients) == number_of_groups+1
                is_csd || abs(coefficients[end]) ≤ 1.0e-14 || error(
                    "A non-CSD score contains a nonzero cutoff response coefficient.",
                )
                if is_csd
                    value += coefficients[end]*cutoff_flux[1,1,ix,iy,iz]
                end
            end
            value /= legacy_normalization
            if mass_normalized
                density = densities[material_index]
                isfinite(density) && density > 0.0 || error(
                    "Mass-normalized scoring requires finite positive density.",
                )
                value /= density
            end
            field[ix,iy,iz] = value
        end
        if !isnothing(physical_normalization)
            field .= apply_normalization(field,physical_normalization)
        end
        channels[process_key] = field
    end

    total = zeros(Float64,voxel_counts[1],voxel_counts[2],voxel_counts[3])
    for field in values(channels)
        total .+= field
    end
    basis = isnothing(physical_normalization) ? :transport_source_basis : :physical_rate
    provenance = Dict(
        "channel_schema" => "radiant.process_response_channel/v1",
        "particle" => get_tag(particle),
        "quantity" => quantity_string,
        "mass_normalized" => string(mass_normalized),
        "physical_normalization" => string(!isnothing(physical_normalization)),
    )
    return Process_Resolved_Score(
        get_tag(particle),quantity_string,channels,total;
        normalization_basis=basis,
        units=_score_output_units(quantity_string,mass_normalized),
        provenance=provenance,
    )
end

function assert_process_score_closure(
    score::Process_Resolved_Score,
    reference::AbstractArray{<:Real};
    rtol::Real=1.0e-9,
    atol::Real=1.0e-12,
)
    reference_array = Float64.(reference)
    expected_shape = size(score.total)
    if ndims(reference_array) == 1 && expected_shape[2] == 1 && expected_shape[3] == 1
        reference_array = reshape(reference_array,expected_shape)
    elseif ndims(reference_array) == 2 && expected_shape[3] == 1
        reference_array = reshape(reference_array,expected_shape)
    end
    size(reference_array) == expected_shape || error(
        "Process-score reference shape does not match the score field.",
    )
    for index in eachindex(score.total)
        isapprox(score.total[index],reference_array[index];rtol=rtol,atol=atol) || error(
            "Process-score closure failed at linear index $(index): " *
            "channels=$(score.total[index]), reference=$(reference_array[index]).",
        )
    end
    return true
end

"""Aggregate raw interaction/type keys by interaction class while preserving exact totals."""
function aggregate_process_scores(score::Process_Resolved_Score)
    aggregated = Dict{String,Array{Float64,3}}()
    for (key,field) in score.channels
        interaction = first(split(key,"/";limit=2))
        if haskey(aggregated,interaction)
            aggregated[interaction] .+= field
        else
            aggregated[interaction] = copy(field)
        end
    end
    total = zeros(Float64,size(score.total))
    for field in values(aggregated)
        total .+= field
    end
    return Process_Resolved_Score(
        score.particle_tag,score.quantity,aggregated,total;
        normalization_basis=score.normalization_basis,
        units=score.units,
        provenance=merge(
            copy(score.provenance),
            Dict("aggregation" => "interaction-class"),
        ),
    )
end
