const HTS_ADDON_COMPONENT_STATUSES = (
    :temporary_core,
    :generic_core_hook,
    :ready_to_extract,
    :blocked,
)

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
        return new(component_id,file_vector,dependency_vector,status,String(notes))
    end
end

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
        schema::AbstractString="radiant.hts_addon_extraction/v3",
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

function _manifest_owned_files(this::HTS_Addon_Extraction_Manifest)
    return sort(unique(String[
        path for component in this.components for path in component.files
    ]))
end

function _present_hts_addon_julia_files()
    return sort(String[
        joinpath("src","hts_addon",name)
        for name in readdir(@__DIR__)
        if endswith(name,".jl")
    ])
end

function validate_addon_extraction_manifest(this::HTS_Addon_Extraction_Manifest)
    for component in this.components
        for path in component.files
            startswith(path,"src/hts_addon/") || error(
                "HTS-specific implementation must remain under src/hts_addon/: $(path).",
            )
            local_path = normpath(joinpath(@__DIR__,"..","..",path))
            isfile(local_path) || error(
                "Extraction manifest references a missing HTS source file: $(path).",
            )
        end
    end
    required_invariants = (
        "no-default-branch-merge",
        "no-production-physics-without-qualified-data",
        "one-owner-per-particle-domain-response",
        "source-and-result-provenance-preserved",
        "synthetic-evidence-cannot-pass-physical-gates",
        "future-addon-extraction-preserved",
        "comparison-responses-excluded-from-production-sums",
        "all-hts-source-files-have-extraction-owners",
        "response-specific-group-condensation-preserved",
        "activation-source-energy-separated-from-deposited-heat",
        "protected-response-uncertainty-and-multi-axis-convergence-preserved",
    )
    for invariant in required_invariants
        invariant in this.extraction_invariants || error(
            "Extraction manifest is missing invariant $(invariant).",
        )
    end
    owned = _manifest_owned_files(this)
    present = _present_hts_addon_julia_files()
    owned == present || error(
        "HTS add-on extraction manifest is incomplete. Missing owners=" *
        string(setdiff(present,owned))*"; stale entries="*string(setdiff(owned,present))*".",
    )
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

function default_hts_addon_extraction_manifest()
    components = HTS_Addon_Component[
        HTS_Addon_Component(
            :addon_infrastructure,
            [
                "src/hts_addon/Addon_Extraction_Manifest.jl",
                "src/hts_addon/CSV_Utilities.jl",
                "src/hts_addon/Completion_Extensions.jl",
            ];
            core_dependencies=["Radiant module include/export wiring"],
            status=:temporary_core,
            notes="The loader becomes RadiantHTS.jl package wiring; CSV utilities remain private add-on helpers.",
        ),
        HTS_Addon_Component(
            :process_resolved_scoring,
            ["src/hts_addon/Process_Resolved_Scoring.jl"];
            core_dependencies=[
                "Multigroup_Cross_Sections response-channel hook","Cross_Sections",
                "Geometry","Flux",
            ],status=:temporary_core,
        ),
        HTS_Addon_Component(
            :layer_heating_ownership,
            ["src/hts_addon/Layer_Heating_Ledger.jl"];
            core_dependencies=["Process_Resolved_Scoring","Source_Normalization"],
            status=:ready_to_extract,
            notes="Separates deposition, heat, storage, escape, recoil/cutoff handoff, prompt/delayed populations, and comparison-only responses.",
        ),
        HTS_Addon_Component(
            :activation_delayed_sources,
            ["src/hts_addon/Activation_Delayed_Source.jl"];
            core_dependencies=[
                "Anisotropic_Volume_Source","Source_Normalization","Layer_Heating_Ledger",
            ],
            status=:ready_to_extract,
            notes="Builds parent-resolved delayed photon/electron/positron sources and separate alpha, recoil, neutrino, and unresolved-energy handoffs from activity fields.",
        ),
        HTS_Addon_Component(
            :response_preserving_group_condensation,
            ["src/hts_addon/Response_Preserving_Group_Condensation.jl"];
            core_dependencies=["SHA","multigroup energy/source conventions"],
            status=:ready_to_extract,
            notes="Creates separate spectrum-weighted condensation receipts for heating, Gd capture, activation, PKA, and neutron-to-photon transfer responses.",
        ),
        HTS_Addon_Component(
            :response_uncertainty_and_convergence,
            ["src/hts_addon/Response_Uncertainty_And_Convergence.jl"];
            core_dependencies=["protected response results","physical reference contract"],
            status=:ready_to_extract,
            notes="Combines correlated and independent uncertainty sources and requires independent source, mesh, angle, energy, field, facet, material, coupling, and time convergence where requested.",
        ),
        HTS_Addon_Component(
            :spatial_field_maps,
            [
                "src/hts_addon/Spatial_Magnetic_Field_Map.jl",
                "src/hts_addon/Spatial_Field_Transport.jl",
            ];
            core_dependencies=[
                "Electromagnetic_Field cell-centred field hook",
                "sn_flux cell-local operator construction",
                "sn_inner_pass cell-local source application",
            ],status=:temporary_core,
        ),
        HTS_Addon_Component(
            :piecewise_flat_atlas_and_curvature,
            [
                "src/hts_addon/Piecewise_Flat_Tape_Atlas.jl",
                "src/hts_addon/Analytic_Curvature_Benchmark.jl",
                "src/hts_addon/Curvature_Convergence_Study.jl",
            ];
            core_dependencies=["Spatial_Magnetic_Field_Map","Physical reference contract"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :faceted_geometry,
            [
                "src/hts_addon/Faceted_Geometry.jl",
                "src/hts_addon/Faceted_Local_Tape_Domain.jl",
            ];
            core_dependencies=[
                "HDF5 interchange helpers","Geometry","Tape_Stack_Definition",
                "Boundary_Angular_Current_Source",
            ],
            status=:ready_to_extract,
            notes="DAGMC/ParaStell facets are validated and converted to conservative local Cartesian domains; direct unstructured SN transport is not claimed.",
        ),
        HTS_Addon_Component(
            :gd_capture_and_self_shielding,
            [
                "src/hts_addon/Gd_Prompt_Capture_Cascade.jl",
                "src/hts_addon/Gd_Cascade_Data_Adapter.jl",
                "src/hts_addon/Gd_Evaluated_Cascade_Adapter.jl",
                "src/hts_addon/Gd_Self_Shielding_Analytic.jl",
            ];
            core_dependencies=[
                "Anisotropic_Volume_Source","Source_Normalization","Particle",
                "Capture_Rate_Field",
            ],
            status=:ready_to_extract,
            notes="Groupwise slab self-shielding is a screening/reference fixture; physical production remains owned by OpenMC/OpenSn evaluated-data transport.",
        ),
        HTS_Addon_Component(
            :subkev_thermalization,
            ["src/hts_addon/SubkeV_Thermalization.jl"];
            core_dependencies=["Energy_Partition"],status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :microdosimetry,
            ["src/hts_addon/Statistical_Microdosimetry.jl"];
            core_dependencies=["Energy_Partition","Random"],status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :cryogenic_electrothermal,
            ["src/hts_addon/Cryogenic_Electrothermal.jl"];
            core_dependencies=["LinearAlgebra","HTS_Detector_Definition"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :material_response_registry,
            [
                "src/hts_addon/Material_Response_Registry.jl",
                "src/hts_addon/Response_Table_Completion.jl",
                "src/hts_addon/Material_Data_Bindings.jl",
            ];
            core_dependencies=["Cryogenic_Electrothermal","Atomistic response tables"],
            status=:ready_to_extract,
            notes="Material_Data_Bindings uses the current Atomistic_Response_Table and vector-registry APIs; YBCO-to-GdBCO substitution remains prohibited without explicit surrogate uncertainty.",
        ),
        HTS_Addon_Component(
            :atomistic_response_tables,
            ["src/hts_addon/Atomistic_Response_Table_Pipeline.jl"];
            core_dependencies=["HDF5","Material_Response_Registry"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :physical_reference_qualification,
            [
                "src/hts_addon/Physical_Reference_Artifact.jl",
                "src/hts_addon/Physical_Reference_Qualification.jl",
            ];
            core_dependencies=["Process_Resolved_Scoring","HDF5"],
            status=:ready_to_extract,
        ),
        HTS_Addon_Component(
            :closed_opensn_radiant_coupling,
            ["src/hts_addon/OpenSn_Radiant_Closed_Coupling.jl"];
            core_dependencies=[
                "Boundary_Angular_Current_Source","Transport_Ownership_Map",
                "outgoing SN boundary angular flux",
            ],status=:ready_to_extract,
        ),
    ]
    manifest = HTS_Addon_Extraction_Manifest(
        components;
        generic_core_touchpoints=[
            "src/Radiant.jl temporary include/export wiring",
            "src/structures/Multigroup_Cross_Sections.jl generic response channels",
            "src/cross_sections/generate_cross_sections.jl response-channel population",
            "src/structures/Electromagnetic_Field.jl generic cell-centred field storage",
            "src/structures/Flux_Per_Particle.jl generic outgoing-boundary storage",
            "src/particle_transport/sn_flux.jl generic cell-local field operator",
            "src/particle_transport/sn_inner_pass.jl generic cell-local operator application",
            "src/particle_transport/sn_boundary_flux.jl generic outgoing-current reconstruction",
            "src/interchange/source_hdf5.jl generic source interchange",
            "Project.toml standard-library dependencies",
            "scripts/geometry/export_dagmc_facets.py external adapter",
            "scripts/data/download_zenodo_record.py external-data acquisition",
        ],
        extraction_invariants=[
            "no-default-branch-merge",
            "no-production-physics-without-qualified-data",
            "one-owner-per-particle-domain-response",
            "source-and-result-provenance-preserved",
            "prompt-and-delayed-ledgers-separated",
            "mean-transport-and-event-sampling-distinguished",
            "synthetic-evidence-cannot-pass-physical-gates",
            "future-addon-extraction-preserved",
            "comparison-responses-excluded-from-production-sums",
            "all-hts-source-files-have-extraction-owners",
            "response-specific-group-condensation-preserved",
            "activation-source-energy-separated-from-deposited-heat",
            "protected-response-uncertainty-and-multi-axis-convergence-preserved",
        ],
        metadata=Dict(
            "hosting_policy" => "temporary-core-only",
            "future_repository" => "RadiantHTS.jl",
            "branch_policy" => "branch-only-no-pr",
            "direct_faceted_transport" => "not-claimed",
            "spatial_field_transport" => "cell-local-within-one-sn-sweep",
            "gd_self_shielding_transport" => "analytic-screening-only",
            "heating_ownership" => "one-additive-owner-per-population-response",
            "group_condensation" => "response-and-reference-spectrum-specific",
            "activation_delayed_sources" => "parent-resolved-and-energy-ledger-separated",
            "uncertainty_and_convergence" => "response-specific-and-multi-axis",
        ),
    )
    validate_addon_extraction_manifest(manifest)
    return manifest
end
