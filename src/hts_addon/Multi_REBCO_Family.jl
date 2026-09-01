const REBCO_MATERIAL_FAMILY_PACKET_SCHEMA = "radiant.hts.rebco_material_family_packet/v1"
const REBCO_COMPARISON_CONTRACT_SCHEMA = "radiant.hts.rebco_comparison_contract/v1"
const REBCO_SELF_SHIELDING_SCHEMA = "radiant.hts.rebco_self_shielding_axis_sweep/v1"
const REBCO_CONSUMER_BINDING_SCHEMA = "radiant.hts.rebco_consumer_binding/v1"
const REBCO_TIER_A_SLAB_SCHEMA = "radiant.hts.rebco_tier_a_slab_input/v1"

const REBCO_FAMILY_IDENTITIES = Dict{String,Tuple{String,Symbol}}(
    "YBCO" => ("Y",:confirmed),
    "GdBCO" => ("Gd",:confirmed),
    "EuBCO" => ("Eu",:confirmed),
    # Preserve the user token without inventing a chemical symbol or element.
    "SaBCO" => ("",:unresolved),
    "SmBCO" => ("Sm",:confirmed),
)
const REBCO_PACKET_STATUSES = (:verification,:candidate,:qualified_input)
const REBCO_COMPARISON_MODES = (:controlled_substitution,:realistic_product)
const REBCO_CONSUMER_CLASSES = (:source,:deposition,:pka,:activation)
const REBCO_BINDING_STATUSES = (:verification,:candidate,:blocked_input,:qualified_input)

"""Reference-candidate slab input. It cannot represent a lot-qualified tape."""
struct REBCO_Tier_A_Slab_Input
    schema::String
    material_tag::String
    density_g_cm3::Float64
    molar_mass_g_mol::Float64
    isotope_capture_barn::Dict{String,Tuple{Float64,Float64}}
    film_thickness_um::Float64
    material_packet_sha256::String
    source_energy_eV::Float64
    reference_physical_candidate::Bool
    lot_specific::Bool
    production_qualified::Bool
end

function REBCO_Tier_A_Slab_Input(
    material_tag::AbstractString,density_g_cm3::Real,molar_mass_g_mol::Real,
    isotope_capture_barn::AbstractDict,film_thickness_um::Real,
    material_packet_sha256::AbstractString;source_energy_eV::Real=0.0253,
)
    family = String(material_tag)
    haskey(REBCO_FAMILY_IDENTITIES,family) || error("Unsupported REBCO material family.")
    REBCO_FAMILY_IDENTITIES[family][2] == :confirmed || error(
        "Tier-A slab transport requires a resolved rare-earth identity.",
    )
    density = Float64(density_g_cm3)
    molar_mass = Float64(molar_mass_g_mol)
    thickness = Float64(film_thickness_um)
    energy = Float64(source_energy_eV)
    all(value -> isfinite(value) && value > 0.0,(density,molar_mass,thickness,energy)) ||
        error("Tier-A slab physical inputs must be finite and positive.")
    digest = String(material_packet_sha256)
    length(digest) == 64 && all(character -> character in "0123456789abcdef",digest) ||
        error("Tier-A material packet requires a lowercase SHA-256 digest.")
    components = Dict{String,Tuple{Float64,Float64}}()
    symbol = REBCO_FAMILY_IDENTITIES[family][1]
    for (nuclide,raw_values) in isotope_capture_barn
        values = Tuple(raw_values)
        length(values) == 2 || error("Each isotope requires abundance and capture cross section.")
        abundance,cross_section = Float64(values[1]),Float64(values[2])
        startswith(string(nuclide),symbol*"-") || error(
            "Tier-A isotope identity does not match the REBCO family.",
        )
        isfinite(abundance) && 0.0 <= abundance <= 1.0 || error(
            "Tier-A isotope abundance must lie in [0,1].",
        )
        isfinite(cross_section) && cross_section >= 0.0 || error(
            "Tier-A capture cross section must be finite and nonnegative.",
        )
        components[string(nuclide)] = (abundance,cross_section)
    end
    isempty(components) && error("Tier-A slab input requires isotope capture components.")
    return REBCO_Tier_A_Slab_Input(
        REBCO_TIER_A_SLAB_SCHEMA,family,density,molar_mass,components,thickness,digest,
        energy,true,false,false,
    )
end

function rebco_tier_a_macroscopic_capture(input::REBCO_Tier_A_Slab_Input)
    input.reference_physical_candidate && !input.lot_specific && !input.production_qualified ||
        error("Tier-A input cannot be lot-specific or production-qualified.")
    avogadro = 6.02214076e23
    rare_earth_density = input.density_g_cm3/input.molar_mass_g_mol*avogadro
    return sum(
        rare_earth_density*abundance*cross_section*1.0e-24
        for (abundance,cross_section) in values(input.isotope_capture_barn)
    )
end

"""One-dimensional reference-candidate capture solve in the local film coordinate."""
function solve_rebco_tier_a_slab(
    input::REBCO_Tier_A_Slab_Input;angle_from_tape_plane_deg::Real,cells::Integer,
    cache_attenuation::Bool=false,
)
    angle = Float64(angle_from_tape_plane_deg)
    isfinite(angle) && 0.0 < angle <= 90.0 || error("Incidence angle must lie in (0,90].")
    cells >= 1 || error("Tier-A slab cell count must be positive.")
    mu = sind(angle)
    sigma = rebco_tier_a_macroscopic_capture(input)
    normal_cell_width_cm = input.film_thickness_um*1.0e-4/cells
    optical_depth_cell = sigma*normal_cell_width_cm/mu
    attenuation = cache_attenuation ? exp(-optical_depth_cell) : 0.0
    incident = 1.0
    captures = Vector{Float64}(undef,cells)
    for index in eachindex(captures)
        transmitted = incident*(cache_attenuation ? attenuation : exp(-optical_depth_cell))
        captures[index] = incident-transmitted
        incident = transmitted
    end
    exact_capture = -expm1(-optical_depth_cell*cells)
    integrated_capture = sum(captures)
    exact_surface_capture_density = sigma/mu
    numerical_peak_capture_density = captures[1]/normal_cell_width_cm
    peak_relative_error = exact_surface_capture_density == 0.0 ? 0.0 :
        abs(numerical_peak_capture_density/exact_surface_capture_density-1.0)
    return Dict{String,Any}(
        "schema" => "radiant.hts.rebco_tier_a_slab_result/v1",
        "material_tag" => input.material_tag,"angle_from_tape_plane_deg" => angle,
        "cells" => Int64(cells),"unknown_count" => Int64(cells),
        "macroscopic_capture_cm_inv" => sigma,"capture_fraction" => integrated_capture,
        "exact_capture_fraction" => exact_capture,
        "integrated_capture_relative_error" =>
            abs(integrated_capture-exact_capture)/max(exact_capture,eps(Float64)),
        "peak_sublayer_capture_relative_error" => peak_relative_error,
        "transmitted_fraction" => incident,"capture_profile" => captures,
        "cache_attenuation" => cache_attenuation,
        "reference_physical_candidate" => true,"lot_specific" => false,
        "production_qualified" => false,
    )
end

_rebco_string_dict(values::AbstractDict) = Dict{String,String}(
    string(key) => string(value) for (key,value) in values
)

"""
Versioned, hash-bound family identity packet. `qualified_input` describes the packet inputs only;
this adapter never promotes transport or a response to physical qualification.
"""
struct REBCO_Material_Family_Packet
    schema::String
    packet_id::String
    material_tag::String
    rare_earth_symbol::String
    identity_status::Symbol
    identity_evidence_hash::String
    composition_hash::String
    product_id::String
    product_construction_hash::String
    controlled_design_hash::String
    status::Symbol
    metadata::Dict{String,String}
end

function REBCO_Material_Family_Packet(
    packet_id::AbstractString,
    material_tag::AbstractString;
    rare_earth_symbol::AbstractString,
    identity_status::Symbol,
    identity_evidence_hash::AbstractString,
    composition_hash::AbstractString,
    product_id::AbstractString,
    product_construction_hash::AbstractString,
    controlled_design_hash::AbstractString,
    status::Symbol=:candidate,
    schema::AbstractString=REBCO_MATERIAL_FAMILY_PACKET_SCHEMA,
    metadata::AbstractDict=Dict{String,String}(),
)
    packet = REBCO_Material_Family_Packet(
        String(schema),String(packet_id),String(material_tag),String(rare_earth_symbol),
        identity_status,String(identity_evidence_hash),String(composition_hash),
        String(product_id),String(product_construction_hash),String(controlled_design_hash),
        status,_rebco_string_dict(metadata),
    )
    validate_rebco_material_packet(packet)
    return packet
end

function validate_rebco_material_packet(packet::REBCO_Material_Family_Packet)
    packet.schema == REBCO_MATERIAL_FAMILY_PACKET_SCHEMA || error(
        "Unsupported REBCO material-family packet schema: $(packet.schema).",
    )
    isempty(packet.packet_id) && error("REBCO packet ID cannot be empty.")
    haskey(REBCO_FAMILY_IDENTITIES,packet.material_tag) || error(
        "Unsupported REBCO material family: $(packet.material_tag).",
    )
    expected_symbol,expected_status = REBCO_FAMILY_IDENTITIES[packet.material_tag]
    packet.rare_earth_symbol == expected_symbol || error(
        "Material tag $(packet.material_tag) does not match rare-earth identity " *
        "$(packet.rare_earth_symbol).",
    )
    packet.identity_status == expected_status || error(
        "Material tag $(packet.material_tag) requires identity status $(expected_status).",
    )
    for (label,value) in (
        ("identity-evidence hash",packet.identity_evidence_hash),
        ("composition hash",packet.composition_hash),
        ("product ID",packet.product_id),
        ("product-construction hash",packet.product_construction_hash),
        ("controlled-design hash",packet.controlled_design_hash),
    )
        isempty(value) && error("REBCO $(label) cannot be empty.")
    end
    packet.status in REBCO_PACKET_STATUSES || error(
        "Unknown REBCO material packet status: $(packet.status).",
    )
    packet.identity_status == :unresolved && packet.status == :qualified_input && error(
        "An unresolved rare-earth identity cannot have qualified-input status.",
    )
    return true
end

function adapt_rebco_material_packet(values::AbstractDict)
    getvalue(name) = begin
        haskey(values,name) && return values[name]
        symbol = Symbol(name)
        haskey(values,symbol) && return values[symbol]
        error("REBCO material packet is missing required field $(name).")
    end
    metadata = haskey(values,"metadata") ? values["metadata"] :
               (haskey(values,:metadata) ? values[:metadata] : Dict{String,String}())
    status = haskey(values,"status") ? Symbol(string(values["status"])) :
             (haskey(values,:status) ? Symbol(string(values[:status])) : :candidate)
    schema = haskey(values,"schema") ? string(values["schema"]) :
             (haskey(values,:schema) ? string(values[:schema]) :
              REBCO_MATERIAL_FAMILY_PACKET_SCHEMA)
    return REBCO_Material_Family_Packet(
        string(getvalue("packet_id")),string(getvalue("material_tag"));
        rare_earth_symbol=string(getvalue("rare_earth_symbol")),
        identity_status=Symbol(string(getvalue("identity_status"))),
        identity_evidence_hash=string(getvalue("identity_evidence_hash")),
        composition_hash=string(getvalue("composition_hash")),
        product_id=string(getvalue("product_id")),
        product_construction_hash=string(getvalue("product_construction_hash")),
        controlled_design_hash=string(getvalue("controlled_design_hash")),
        status=status,schema=schema,metadata=metadata,
    )
end

rebco_packet_is_physically_qualified(::REBCO_Material_Family_Packet) = false

"""Comparison-only contract that distinguishes controlled substitution from real products."""
struct REBCO_Comparison_Contract
    schema::String
    comparison_id::String
    mode::Symbol
    reference_packet_id::String
    candidate_packet_id::String
    response_ids::Vector{String}
    paired_axes::Vector{String}
    limitations::Vector{String}
    comparison_only::Bool
    physical_qualification::Bool
end

function REBCO_Comparison_Contract(
    comparison_id::AbstractString,
    reference::REBCO_Material_Family_Packet,
    candidate::REBCO_Material_Family_Packet;
    mode::Symbol,
    response_ids::AbstractVector{<:AbstractString},
    paired_axes::AbstractVector{<:AbstractString},
    limitations::AbstractVector{<:AbstractString}=String[],
)
    validate_rebco_material_packet(reference)
    validate_rebco_material_packet(candidate)
    isempty(comparison_id) && error("REBCO comparison ID cannot be empty.")
    mode in REBCO_COMPARISON_MODES || error("Unknown REBCO comparison mode: $(mode).")
    reference.packet_id != candidate.packet_id || error(
        "A family comparison requires distinct material packets.",
    )
    reference.material_tag != candidate.material_tag || error(
        "A multi-family comparison requires distinct REBCO families.",
    )
    responses = String[string(value) for value in response_ids]
    axes = String[string(value) for value in paired_axes]
    isempty(responses) && error("A REBCO comparison requires at least one response.")
    isempty(axes) && error("A REBCO comparison requires explicitly paired axes.")
    length(unique(responses)) == length(responses) || error(
        "REBCO comparison response IDs must be unique.",
    )
    if mode == :controlled_substitution
        reference.controlled_design_hash == candidate.controlled_design_hash || error(
            "Controlled substitution requires the same controlled-design hash.",
        )
        reference.product_construction_hash == candidate.product_construction_hash || error(
            "Controlled substitution requires matched non-RE construction.",
        )
    else
        reference.product_id != candidate.product_id || error(
            "A realistic-product comparison requires distinct product IDs.",
        )
        reference.product_construction_hash != candidate.product_construction_hash || error(
            "A realistic-product comparison requires independently bound constructions.",
        )
    end
    limitation_vector = String[string(value) for value in limitations]
    if mode == :realistic_product
        required = "product differences cannot be attributed solely to rare-earth identity"
        required in limitation_vector || push!(limitation_vector,required)
    end
    return REBCO_Comparison_Contract(
        REBCO_COMPARISON_CONTRACT_SCHEMA,String(comparison_id),mode,
        reference.packet_id,candidate.packet_id,responses,axes,limitation_vector,true,false,
    )
end

"""Isotope-resolved rare-earth capture input for analytic screening."""
struct REBCO_Isotope_Capture_Component
    family_packet_id::String
    nuclide::String
    atomic_density_atoms_cm3::Float64
    energy_edges_eV::Vector{Float64}
    capture_xs_barn::Vector{Float64}
    data_hash::String
    status::Symbol
end

function REBCO_Isotope_Capture_Component(
    packet::REBCO_Material_Family_Packet,
    nuclide::AbstractString,
    atomic_density_atoms_cm3::Real,
    energy_edges_eV::AbstractVector{<:Real},
    capture_xs_barn::AbstractVector{<:Real};
    data_hash::AbstractString,
    status::Symbol=:candidate,
)
    validate_rebco_material_packet(packet)
    packet.identity_status == :confirmed || error(
        "Isotope-resolved capture data cannot bind to an unresolved RE identity.",
    )
    startswith(String(nuclide),packet.rare_earth_symbol*"-") || error(
        "Nuclide $(nuclide) does not preserve packet rare-earth identity " *
        "$(packet.rare_earth_symbol).",
    )
    density = Float64(atomic_density_atoms_cm3)
    isfinite(density) && density >= 0.0 || error(
        "REBCO isotope atomic density must be finite and nonnegative.",
    )
    edges = Float64.(energy_edges_eV)
    length(edges) >= 2 && all(isfinite,edges) && all(edges .>= 0.0) &&
        all(diff(edges) .> 0.0) || error(
        "REBCO capture energy boundaries must be finite and strictly increasing.",
    )
    cross_sections = Float64.(capture_xs_barn)
    length(cross_sections) == length(edges)-1 || error(
        "REBCO capture data require one value per energy group.",
    )
    all(value -> isfinite(value) && value >= 0.0,cross_sections) || error(
        "REBCO capture cross sections must be finite and nonnegative.",
    )
    isempty(data_hash) && error("REBCO capture-data hash cannot be empty.")
    status in REBCO_PACKET_STATUSES || error("Unknown REBCO capture-data status.")
    return REBCO_Isotope_Capture_Component(
        packet.packet_id,String(nuclide),density,edges,cross_sections,String(data_hash),status,
    )
end

struct REBCO_Self_Shielding_Case
    family_packet_id::String
    material_tag::String
    thickness_cm::Float64
    turn_index::Int64
    incidence_direction::NTuple{3,Float64}
    layer_normal::NTuple{3,Float64}
    absolute_direction_cosine::Float64
    path_length_cm::Float64
    incident_current_per_s::Vector{Float64}
    transmitted_current_per_s::Vector{Float64}
    capture_rate_per_s::Dict{String,Vector{Float64}}
    particle_balance_residual_per_s::Vector{Float64}
end

struct REBCO_Self_Shielding_Axis_Result
    schema::String
    family_packet_id::String
    material_tag::String
    thickness_axis_cm::Vector{Float64}
    turn_axis::Vector{Int64}
    incidence_axis::Vector{NTuple{3,Float64}}
    layer_normal::NTuple{3,Float64}
    cases::Vector{REBCO_Self_Shielding_Case}
    geometry_hash::String
    transport_artifact_hash::String
    data_hashes::Vector{String}
    material_packet_status::Symbol
    isotope_data_status::Symbol
    screening_only::Bool
    physical_transport_qualification::Bool
end

function _rebco_unit_vector(values::AbstractVector{<:Real},label::AbstractString)
    length(values) == 3 || error("$(label) must have three components.")
    vector = Float64.(values)
    all(isfinite,vector) || error("$(label) must be finite.")
    magnitude = sqrt(sum(abs2,vector))
    magnitude > 0.0 || error("$(label) cannot be zero.")
    normalized = vector./magnitude
    return (normalized[1],normalized[2],normalized[3])
end

function solve_rebco_self_shielding_axis_sweep(
    packet::REBCO_Material_Family_Packet,
    components::AbstractVector{REBCO_Isotope_Capture_Component},
    incident_current_per_s::AbstractVector{<:Real};
    thickness_axis_cm::AbstractVector{<:Real},
    turn_axis::AbstractVector{<:Integer},
    incidence_axis::AbstractVector{<:AbstractVector{<:Real}},
    layer_normal::AbstractVector{<:Real},
    geometry_hash::AbstractString,
    transport_artifact_hash::AbstractString="analytic-rebco-self-shielding",
)
    validate_rebco_material_packet(packet)
    packet.identity_status == :confirmed || error(
        "Analytic isotope self-shielding cannot run for an unresolved RE identity.",
    )
    component_vector = REBCO_Isotope_Capture_Component[components...]
    isempty(component_vector) && error("REBCO self-shielding requires isotope components.")
    all(component -> component.family_packet_id == packet.packet_id,component_vector) || error(
        "Every capture component must bind to the selected material packet.",
    )
    nuclides = getfield.(component_vector,:nuclide)
    length(unique(nuclides)) == length(nuclides) || error(
        "REBCO self-shielding isotope components must be unique.",
    )
    edges = component_vector[1].energy_edges_eV
    all(component -> component.energy_edges_eV == edges,component_vector) || error(
        "Every REBCO isotope component must use identical energy groups.",
    )
    incident = Float64.(incident_current_per_s)
    length(incident) == length(edges)-1 || error(
        "Incident current must contain one value per energy group.",
    )
    all(value -> isfinite(value) && value >= 0.0,incident) || error(
        "Incident current must be finite and nonnegative.",
    )
    thicknesses = Float64.(thickness_axis_cm)
    isempty(thicknesses) && error("Thickness axis cannot be empty.")
    all(value -> isfinite(value) && value > 0.0,thicknesses) || error(
        "Thickness axis must be finite and positive.",
    )
    turns = Int64.(turn_axis)
    isempty(turns) && error("Turn axis cannot be empty.")
    all(>=(1),turns) || error("Turn indices must be positive.")
    length(unique(turns)) == length(turns) || error("Turn indices must be unique.")
    normal = _rebco_unit_vector(layer_normal,"Layer normal")
    directions = NTuple{3,Float64}[
        _rebco_unit_vector(direction,"Incidence direction") for direction in incidence_axis
    ]
    isempty(directions) && error("Incidence axis cannot be empty.")
    isempty(geometry_hash) && error("REBCO geometry hash cannot be empty.")
    isempty(transport_artifact_hash) && error("REBCO transport-artifact hash cannot be empty.")
    macroscopic = Dict(
        component.nuclide => component.atomic_density_atoms_cm3 .* component.capture_xs_barn .* 1.0e-24
        for component in component_vector
    )
    cases = REBCO_Self_Shielding_Case[]
    for thickness in thicknesses, turn in turns, direction in directions
        mu = abs(sum(direction[index]*normal[index] for index in 1:3))
        mu > sqrt(eps(Float64)) || error(
            "Grazing incidence has an unbounded analytic path and must be resolved geometrically.",
        )
        path = thickness/mu
        transmitted = copy(incident)
        captures = Dict(nuclide => zeros(Float64,length(incident)) for nuclide in nuclides)
        for group in eachindex(incident)
            sigma_total = sum(values[group] for values in values(macroscopic))
            if sigma_total > 0.0 && incident[group] > 0.0
                removed = incident[group]*(-expm1(-sigma_total*path))
                transmitted[group] -= removed
                for (nuclide,sigma) in macroscopic
                    captures[nuclide][group] = removed*sigma[group]/sigma_total
                end
            end
        end
        balance = incident-transmitted
        for rates in values(captures)
            balance .-= rates
        end
        push!(cases,REBCO_Self_Shielding_Case(
            packet.packet_id,packet.material_tag,thickness,turn,direction,normal,mu,path,
            copy(incident),transmitted,captures,balance,
        ))
    end
    isotope_status = all(component -> component.status == :qualified_input,component_vector) ?
        :qualified_input :
        (all(component -> component.status != :verification,component_vector) ? :candidate :
                                                                            :verification)
    input_hashes = sort(unique(vcat(
        getfield.(component_vector,:data_hash),
        [packet.identity_evidence_hash,packet.composition_hash],
    )))
    return REBCO_Self_Shielding_Axis_Result(
        REBCO_SELF_SHIELDING_SCHEMA,packet.packet_id,packet.material_tag,thicknesses,turns,
        directions,normal,cases,String(geometry_hash),String(transport_artifact_hash),
        input_hashes,packet.status,isotope_status,true,false,
    )
end

function rebco_self_shielding_receipt(result::REBCO_Self_Shielding_Axis_Result)
    maximum_residual = maximum(
        maximum(abs.(case.particle_balance_residual_per_s)) for case in result.cases
    )
    return Dict{String,Any}(
        "schema" => result.schema,
        "family_packet_id" => result.family_packet_id,
        "material_tag" => result.material_tag,
        "thickness_axis_cm" => copy(result.thickness_axis_cm),
        "turn_axis" => copy(result.turn_axis),
        "incidence_axis" => [collect(direction) for direction in result.incidence_axis],
        "path_lengths_cm" => [case.path_length_cm for case in result.cases],
        "isotope_identities" => sort(collect(keys(result.cases[1].capture_rate_per_s))),
        "maximum_particle_balance_residual_per_s" => maximum_residual,
        "screening_only" => true,
        "physical_transport_qualification" => false,
        "geometry_hash" => result.geometry_hash,
        "transport_artifact_hash" => result.transport_artifact_hash,
        "data_hashes" => copy(result.data_hashes),
        "material_packet_status" => string(result.material_packet_status),
        "isotope_data_status" => string(result.isotope_data_status),
    )
end

rebco_self_shielding_is_physically_qualified(::REBCO_Self_Shielding_Axis_Result) = false

"""Hash-bound family consumer; physical qualification is deliberately unavailable here."""
struct REBCO_Family_Consumer_Binding
    schema::String
    binding_id::String
    family_packet_id::String
    material_tag::String
    rare_earth_symbol::String
    consumer_class::Symbol
    producer_artifact_hash::String
    consumer_artifact_hash::String
    data_hashes::Vector{String}
    status::Symbol
    comparison_only::Bool
    physical_qualification::Bool
    metadata::Dict{String,String}
end

function bind_rebco_family_consumer(
    packet::REBCO_Material_Family_Packet,
    consumer_class::Symbol;
    binding_id::AbstractString,
    producer_artifact_hash::AbstractString,
    consumer_artifact_hash::AbstractString,
    data_hashes::AbstractVector{<:AbstractString},
    status::Symbol=:candidate,
    comparison_only::Bool=false,
    metadata::AbstractDict=Dict{String,String}(),
)
    validate_rebco_material_packet(packet)
    consumer_class in REBCO_CONSUMER_CLASSES || error(
        "Unknown REBCO consumer class: $(consumer_class).",
    )
    status in REBCO_BINDING_STATUSES || error("Unknown REBCO consumer binding status.")
    for (label,value) in (
        ("binding ID",binding_id),("producer-artifact hash",producer_artifact_hash),
        ("consumer-artifact hash",consumer_artifact_hash),
    )
        isempty(value) && error("REBCO $(label) cannot be empty.")
    end
    hashes = String[string(value) for value in data_hashes]
    isempty(hashes) && error("REBCO consumer binding requires explicit data hashes.")
    any(isempty,hashes) && error("REBCO consumer data hashes cannot be empty.")
    if consumer_class in (:pka,:activation) && packet.identity_status != :confirmed
        error("$(consumer_class) binding requires a confirmed rare-earth identity.")
    end
    if packet.identity_status == :unresolved && status != :blocked_input
        error("An unresolved-family binding must remain :blocked_input.")
    end
    return REBCO_Family_Consumer_Binding(
        REBCO_CONSUMER_BINDING_SCHEMA,String(binding_id),packet.packet_id,
        packet.material_tag,packet.rare_earth_symbol,consumer_class,
        String(producer_artifact_hash),String(consumer_artifact_hash),sort(unique(hashes)),
        status,comparison_only,false,_rebco_string_dict(metadata),
    )
end
