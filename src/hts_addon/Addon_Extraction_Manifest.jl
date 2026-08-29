const HTS_ADDON_COMPONENT_STATUSES = (
    :temporary_core,
    :generic_core_hook,
    :ready_to_extract,
    :blocked,
)

"""One independently extractable component of the temporary HTS extension."""
struct HTS_Addon_Component
    component_id::Symbol
    files::Vector{String}
    core_dependencies::Vector{String}
    status::Symbol
    notes::String

    function HTS_Addon_Component(
        component_id::Symbol,
        files::AbstractVector{<:AbstractString};
        core_dependencies::AbstractVector{<:AbstractString}=String[],
        status::Symbol=:temporary_core,
        notes::AbstractString="",
    )
        component_id == :none && error("HTS add-on component identifier cannot be :none.")
        status in HTS_ADDON_COMPONENT_STATUSES || error(
            "Unknown HTS add-on component status: $(status).",
        )
        file_vector = String[string(value) for value in files]
        isempty(file_vector) && error("An HTS add-on component must own at least one file.")
        length(unique(file_vector)) == length(file_vector) || error(
            "HTS add-on component file paths must be unique.",
        )
        dependency_vector = String[string(value) for value in core_dependencies]
        return new(
            component_id,file_vector,dependency_vector,status,String(notes),
        )
    end
end

"""
    HTS_Addon_Extraction_Manifest

Machine-readable boundary for code temporarily hosted in `Radiant.jl` but intended to become the
separate `RadiantHTS.jl` package. New HTS-specific implementation should remain under
`src/hts_addon/`; changes outside that directory must be generic hooks or include/export wiring.
"""
struct HTS_Addon_Extraction_Manifest
    schema::String
    temporary_host_repository::String
    intended_package::String
    components::Vector{HTS_Addon_Component}
    generic_core_touchpoints::Vector{String}
    extraction_invariants::Vector{String}
    metadata::Dict{String,String}

    function HTS_Addon_Extraction_Manifest(
        components::AbstractVector{HTS_Addon_Component};
        schema::AbstractString="radiant.hts_addon_extraction/v1",
        temporary_host_repository::AbstractString="FusionSandwich/Radiant.jl",
        intended_package::AbstractString="RadiantHTS.jl",
        generic_core_touchpoints::AbstractVector{<:AbstractString}=String[],
        extraction_invariants::AbstractVector{<:AbstractString}=String[],
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(schema) && error("Extraction manifest schema cannot be empty.")
        isempty(temporary_host_repository) && error("Temporary host repository cannot be empty.")
        isempty(intended_package) && error("Intended add-on package cannot be empty.")
        component_vector = HTS_Addon_Component[components...]
        isempty(component_vector) && error("Extraction manifest cannot be empty.")
        identifiers = getfield.(component_vector,:component_id)
        length(unique(identifiers)) == length(identifiers) || error(
            "HTS add-on component identifiers must be unique.",
        )
        all_files = String[]
        for component in component_vector
            append!(all_files,component.files)
        end
        length(unique(all_files)) == length(all_files) || error(
            "A file may be owned by only one extraction component.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(schema),String(temporary_host_repository),String(intended_package),
            component_vector,String[string(value) for value in generic_core_touchpoints],
            String[string(value) for value in extraction_invariants],metadata_string,
        )
    end
end

function validate_addon_extraction_manifest(this::HTS_Addon_Extraction_Manifest)
    for component in this.components
        for path in component.files
            startswith(path,"src/hts_addon/") || error(
                "HTS-specific implementation must remain under src/hts_addon/: $(path).",
            )
        end
    end
    required_invariants = (
        "no-default-branch-merge",
        "no-production-physics-without-qualified-data",
        "one-owner-per-particle-domain-response",
        "source-and-result-provenance-preserved",
    )
    for invariant in required_invariants
        invariant in this.extraction_invariants || error(
            "Extraction manifest is missing invariant $(invariant).",
        )
    end
    return true
end

function get_addon_component(
    this::HTS_Addon_Extraction_Manifest,
    component_id::Symbol,
)
    index = findfirst(component -> component.component_id == component_id,this.components)
    isnothing(index) && error("Unknown HTS add-on component: $(component_id).")
    return this.components[index]
end

"""Return the current temporary-core extraction plan for all HTS-specific capabilities."""
function default_hts_addon_extraction_manifest()
    components = HTS_Addon_Component[
        HTS_Addon_Component(
            :process_resolved_scoring,
            ["src/hts_addon/Process_Resolved_Scoring.jl"];
            core_dependencies=[
                "Multigroup_Cross_Sections response-channel hook",
                "Cross_Sections",
                "Geometry",
                "Flux",
            ],
            status=:temporary_core,
        ),
        HTS_Addon_Component(
            :spatial_field_maps,
            ["src/hts_addon/Spatial_Magnetic_Field_Map.jl"];
            core_dependencies=["Electromagnetic_Field","Geometry"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :piecewise_flat_atlas,
            ["src/hts_addon/Piecewise_Flat_Tape_Atlas.jl"];
            core_dependencies=["Spatial_Magnetic_Field_Map"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :gd_capture_cascade,
            ["src/hts_addon/Gd_Prompt_Capture_Cascade.jl"];
            core_dependencies=["Anisotropic_Volume_Source","Source_Normalization","Particle"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :subkev_thermalization,
            ["src/hts_addon/SubkeV_Thermalization.jl"];
            core_dependencies=["Energy_Partition"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :microdosimetry,
            ["src/hts_addon/Statistical_Microdosimetry.jl"];
            core_dependencies=["Energy_Partition","Random"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :cryogenic_electrothermal,
            ["src/hts_addon/Cryogenic_Electrothermal.jl"];
            core_dependencies=["LinearAlgebra","HTS_Detector_Definition"],
            status=:ready_to_extract,
        ),
    ]
    manifest = HTS_Addon_Extraction_Manifest(
        components;
        generic_core_touchpoints=[
            "src/Radiant.jl include/export wiring",
            "src/structures/Multigroup_Cross_Sections.jl generic response channels",
            "src/cross_sections/generate_cross_sections.jl response-channel population",
            "Project.toml stdlib dependencies",
        ],
        extraction_invariants=[
            "no-default-branch-merge",
            "no-production-physics-without-qualified-data",
            "one-owner-per-particle-domain-response",
            "source-and-result-provenance-preserved",
            "prompt-and-delayed-ledgers-separated",
            "mean-transport-and-event-sampling-distinguished",
        ],
        metadata=Dict(
            "hosting_policy" => "temporary-core-only",
            "future_repository" => "RadiantHTS.jl",
            "branch_policy" => "branch-only-no-pr",
        ),
    )
    validate_addon_extraction_manifest(manifest)
    return manifest
end
