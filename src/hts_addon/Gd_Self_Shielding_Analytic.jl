const GD_SELF_SHIELDING_DATA_STATUSES = (:verification,:candidate,:qualified)

"""Groupwise microscopic capture data and atomic density for one Gd isotope."""
struct Gd_Groupwise_Capture_Component
    nuclide::String
    atomic_density_atoms_cm3::Float64
    energy_edges_eV::Vector{Float64}
    capture_xs_barn::Vector{Float64}
    data_hash::String
    status::Symbol
    provenance::Dict{String,String}

    function Gd_Groupwise_Capture_Component(
        nuclide::AbstractString,
        atomic_density_atoms_cm3::Real,
        energy_edges_eV::AbstractVector{<:Real},
        capture_xs_barn::AbstractVector{<:Real};
        data_hash::AbstractString,
        status::Symbol=:candidate,
        provenance::AbstractDict=Dict{String,String}(),
    )
        startswith(String(nuclide),"Gd-") || error(
            "Gd self-shielding components require isotope-resolved Gd tags.",
        )
        density = Float64(atomic_density_atoms_cm3)
        isfinite(density) && density >= 0.0 || error(
            "Gd atomic density must be finite and nonnegative.",
        )
        edges = Float64.(energy_edges_eV)
        length(edges) >= 2 && all(isfinite,edges) && all(edges .>= 0.0) &&
            all(diff(edges) .> 0.0) || error(
            "Gd capture energy boundaries must be finite and strictly increasing.",
        )
        cross_sections = Float64.(capture_xs_barn)
        length(cross_sections) == length(edges)-1 || error(
            "Gd capture cross sections must contain one value per energy group.",
        )
        all(value -> isfinite(value) && value >= 0.0,cross_sections) || error(
            "Gd capture cross sections must be finite and nonnegative.",
        )
        isempty(data_hash) && error("Gd capture-data hash cannot be empty.")
        status in GD_SELF_SHIELDING_DATA_STATUSES || error(
            "Unknown Gd capture-data status: $(status).",
        )
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            String(nuclide),density,edges,cross_sections,String(data_hash),status,
            provenance_string,
        )
    end
end

macroscopic_capture_cm_inv(component::Gd_Groupwise_Capture_Component) =
    component.atomic_density_atoms_cm3 .* component.capture_xs_barn .* 1.0e-24

"""One physical Gd-bearing slab layer, optionally subdivided for capture-depth scoring."""
struct Gd_Self_Shielding_Layer
    layer_id::String
    material_tag::String
    thickness_cm::Float64
    subdivisions::Int64
    components::Vector{Gd_Groupwise_Capture_Component}
    background_absorption_cm_inv::Vector{Float64}
    metadata::Dict{String,String}

    function Gd_Self_Shielding_Layer(
        layer_id::AbstractString,
        material_tag::AbstractString,
        thickness_cm::Real,
        components::AbstractVector{Gd_Groupwise_Capture_Component};
        subdivisions::Integer=1,
        background_absorption_cm_inv::Union{Nothing,AbstractVector{<:Real}}=nothing,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(layer_id) && error("Gd self-shielding layer ID cannot be empty.")
        isempty(material_tag) && error("Gd self-shielding material tag cannot be empty.")
        thickness = Float64(thickness_cm)
        isfinite(thickness) && thickness > 0.0 || error(
            "Gd self-shielding layer thickness must be positive.",
        )
        subdivisions >= 1 || error("Gd self-shielding subdivisions must be positive.")
        component_vector = Gd_Groupwise_Capture_Component[components...]
        isempty(component_vector) && error(
            "A Gd self-shielding layer requires at least one isotope component.",
        )
        nuclides = getfield.(component_vector,:nuclide)
        length(unique(nuclides)) == length(nuclides) || error(
            "Gd isotope components must be unique within a layer.",
        )
        reference_edges = component_vector[1].energy_edges_eV
        all(component -> component.energy_edges_eV == reference_edges,component_vector) || error(
            "All Gd isotope components in a layer must use identical energy groups.",
        )
        background = isnothing(background_absorption_cm_inv) ?
            zeros(Float64,length(reference_edges)-1) :
            Float64.(background_absorption_cm_inv)
        length(background) == length(reference_edges)-1 || error(
            "Background absorption must contain one value per energy group.",
        )
        all(value -> isfinite(value) && value >= 0.0,background) || error(
            "Background absorption must be finite and nonnegative.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(layer_id),String(material_tag),thickness,Int64(subdivisions),
            component_vector,background,metadata_string,
        )
    end
end

struct Gd_Self_Shielding_Cell_Result
    layer_id::String
    material_tag::String
    cell_index_in_layer::Int64
    depth_start_cm::Float64
    depth_stop_cm::Float64
    entering_current_per_s::Vector{Float64}
    exiting_current_per_s::Vector{Float64}
    capture_rate_per_s::Dict{String,Vector{Float64}}
    background_absorption_rate_per_s::Vector{Float64}
end

"""
    Gd_Groupwise_Self_Shielding_Result

Uncollided, pure-removal slab solution used to screen Gd capture self-shielding and produce
capture-depth fixtures. It is not a replacement for continuous-energy OpenMC or a scattering
OpenSn solve. Physical promotion therefore remains false even when evaluated capture inputs are
used.
"""
struct Gd_Groupwise_Self_Shielding_Result
    schema::String
    energy_edges_eV::Vector{Float64}
    incident_current_per_s::Vector{Float64}
    direction_cosine::Float64
    cells::Vector{Gd_Self_Shielding_Cell_Result}
    transmitted_current_per_s::Vector{Float64}
    capture_rate_per_s::Dict{String,Vector{Float64}}
    background_absorption_rate_per_s::Vector{Float64}
    particle_balance_residual_per_s::Vector{Float64}
    geometry_hash::String
    transport_artifact_hash::String
    input_data_status::Symbol
    data_hashes::Vector{String}
    metadata::Dict{String,String}
end

function solve_gd_groupwise_self_shielding(
    layers::AbstractVector{Gd_Self_Shielding_Layer},
    incident_current_per_s::AbstractVector{<:Real};
    direction_cosine::Real=1.0,
    geometry_hash::AbstractString,
    transport_artifact_hash::AbstractString="analytic-gd-self-shielding",
    metadata::AbstractDict=Dict{String,String}(),
)
    layer_vector = Gd_Self_Shielding_Layer[layers...]
    isempty(layer_vector) && error("Gd self-shielding calculation requires at least one layer.")
    reference_edges = layer_vector[1].components[1].energy_edges_eV
    for layer in layer_vector
        layer.components[1].energy_edges_eV == reference_edges || error(
            "Every Gd self-shielding layer must use the same energy groups.",
        )
    end
    incident = Float64.(incident_current_per_s)
    length(incident) == length(reference_edges)-1 || error(
        "Incident current must contain one value per energy group.",
    )
    all(value -> isfinite(value) && value >= 0.0,incident) || error(
        "Incident group currents must be finite and nonnegative.",
    )
    mu = Float64(direction_cosine)
    isfinite(mu) && 0.0 < mu <= 1.0 || error(
        "Gd self-shielding direction cosine must lie in (0,1].",
    )
    isempty(geometry_hash) && error("Gd self-shielding geometry hash cannot be empty.")
    isempty(transport_artifact_hash) && error(
        "Gd self-shielding transport-artifact hash cannot be empty.",
    )

    all_components = [component for layer in layer_vector for component in layer.components]
    nuclides = sort(unique(getfield.(all_components,:nuclide)))
    total_capture = Dict(nuclide => zeros(Float64,length(incident)) for nuclide in nuclides)
    background_total = zeros(Float64,length(incident))
    current = copy(incident)
    cells = Gd_Self_Shielding_Cell_Result[]
    global_depth = 0.0

    for layer in layer_vector
        cell_thickness = layer.thickness_cm/layer.subdivisions
        macroscopic = Dict(
            component.nuclide => macroscopic_capture_cm_inv(component)
            for component in layer.components
        )
        for cell_index in 1:layer.subdivisions
            entering = copy(current)
            captures = Dict(nuclide => zeros(Float64,length(incident)) for nuclide in nuclides)
            background = zeros(Float64,length(incident))
            for group in eachindex(incident)
                sigma_components = sum(values[group] for values in values(macroscopic))
                sigma_total = sigma_components+layer.background_absorption_cm_inv[group]
                chord_cm = cell_thickness/mu
                if sigma_total > 0.0 && entering[group] > 0.0
                    removed = entering[group]*(-expm1(-sigma_total*chord_cm))
                    for (nuclide,sigma) in macroscopic
                        captures[nuclide][group] = removed*sigma[group]/sigma_total
                        total_capture[nuclide][group] += captures[nuclide][group]
                    end
                    background[group] = removed*layer.background_absorption_cm_inv[group]/
                                        sigma_total
                    background_total[group] += background[group]
                    current[group] = entering[group]-removed
                end
            end
            push!(cells,Gd_Self_Shielding_Cell_Result(
                layer.layer_id,layer.material_tag,Int64(cell_index),global_depth,
                global_depth+cell_thickness,entering,copy(current),captures,background,
            ))
            global_depth += cell_thickness
        end
    end

    balance = copy(incident)-current-background_total
    for rates in values(total_capture)
        balance .-= rates
    end
    input_status = all(component -> component.status == :qualified,all_components) ?
        :qualified :
        (all(component -> component.status != :verification,all_components) ? :candidate :
                                                                       :verification)
    metadata_string = Dict{String,String}(
        "transport_model" => "uncollided-pure-removal-groupwise-slab",
        "scattering_included" => "false",
        "screening_only" => "true",
        "physical_transport_qualification" => "false",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Gd_Groupwise_Self_Shielding_Result(
        "radiant.hts.gd_groupwise_self_shielding/v1",copy(reference_edges),incident,mu,
        cells,copy(current),total_capture,background_total,balance,String(geometry_hash),
        String(transport_artifact_hash),input_status,
        sort(unique(getfield.(all_components,:data_hash))),metadata_string,
    )
end

function gd_groupwise_self_shielding_factor(
    result::Gd_Groupwise_Self_Shielding_Result,
    layers::AbstractVector{Gd_Self_Shielding_Layer},
    nuclide::AbstractString,
)
    key = String(nuclide)
    haskey(result.capture_rate_per_s,key) || error(
        "Gd self-shielding result does not contain nuclide $(key).",
    )
    thin_limit = zeros(Float64,length(result.incident_current_per_s))
    for layer in layers
        component_index = findfirst(component -> component.nuclide == key,layer.components)
        isnothing(component_index) && continue
        sigma = macroscopic_capture_cm_inv(layer.components[component_index])
        thin_limit .+= result.incident_current_per_s .* sigma .* layer.thickness_cm ./
                      result.direction_cosine
    end
    factors = ones(Float64,length(thin_limit))
    for group in eachindex(factors)
        if thin_limit[group] > 0.0
            factors[group] = result.capture_rate_per_s[key][group]/thin_limit[group]
        elseif result.capture_rate_per_s[key][group] > 0.0
            error("Positive capture rate has a zero thin-limit denominator.")
        end
    end
    return factors
end

function capture_rate_field_from_gd_self_shielding(
    result::Gd_Groupwise_Self_Shielding_Result,
    nuclide::AbstractString,
    voxel_ids::AbstractVector{<:Integer},
    voxel_volumes_cm3::AbstractVector{<:Real};
    material_tags::AbstractVector{<:AbstractString}=String[],
    provenance::AbstractDict=Dict{String,String}(),
)
    key = String(nuclide)
    haskey(result.capture_rate_per_s,key) || error(
        "Gd self-shielding result does not contain nuclide $(key).",
    )
    length(voxel_ids) == length(result.cells) || error(
        "Capture-rate voxel identifiers must match the self-shielding cell count.",
    )
    cell_rates = [sum(cell.capture_rate_per_s[key]) for cell in result.cells]
    tags = isempty(material_tags) ? getfield.(result.cells,:material_tag) :
                                   String[string(value) for value in material_tags]
    extra = Dict{String,String}(
        "producer" => "solve_gd_groupwise_self_shielding",
        "screening_only" => "true",
        "capture_rates_are_self_shielded" => "true",
        "geometry_hash" => result.geometry_hash,
        "input_data_status" => string(result.input_data_status),
    )
    for (key_value,value) in provenance
        extra[string(key_value)] = string(value)
    end
    return Capture_Rate_Field(
        key,voxel_ids,voxel_volumes_cm3,cell_rates;
        material_tags=tags,transport_artifact_hash=result.transport_artifact_hash,
        provenance=extra,
    )
end

function gd_self_shielding_receipt(result::Gd_Groupwise_Self_Shielding_Result)
    maximum_balance = maximum(abs.(result.particle_balance_residual_per_s))
    return Dict{String,Any}(
        "schema" => result.schema,
        "energy_group_count" => length(result.energy_edges_eV)-1,
        "cell_count" => length(result.cells),
        "direction_cosine" => result.direction_cosine,
        "incident_current_per_s" => result.incident_current_per_s,
        "transmitted_current_per_s" => result.transmitted_current_per_s,
        "capture_rate_per_s" => Dict(
            key => copy(value) for (key,value) in result.capture_rate_per_s
        ),
        "background_absorption_rate_per_s" => result.background_absorption_rate_per_s,
        "maximum_particle_balance_residual_per_s" => maximum_balance,
        "particle_balance_pass" => maximum_balance <= 1.0e-10*
            max(maximum(result.incident_current_per_s),1.0),
        "geometry_hash" => result.geometry_hash,
        "transport_artifact_hash" => result.transport_artifact_hash,
        "input_data_status" => string(result.input_data_status),
        "data_hashes" => result.data_hashes,
        "screening_only" => true,
        "physical_transport_qualification" => false,
        "metadata" => copy(result.metadata),
    )
end
