const PHYSICS_COVERAGE_STATUSES = (
    :implemented,
    :external_production,
    :comparison_only,
    :deferred,
    :unsupported,
)

"""Ownership and qualification state for one physical process in one model domain."""
struct Physics_Coverage_Record
    process_id::String
    domain::String
    incoming_particles::Vector{String}
    outgoing_particles::Vector{String}
    production_owner::String
    status::Symbol
    energy_min_eV::Float64
    energy_max_eV::Float64
    handoff_schema::String
    notes::String

    function Physics_Coverage_Record(;
        process_id::AbstractString,
        domain::AbstractString,
        incoming_particles::AbstractVector{<:AbstractString},
        outgoing_particles::AbstractVector{<:AbstractString}=String[],
        production_owner::AbstractString,
        status::Symbol,
        energy_min_eV::Real=0.0,
        energy_max_eV::Real=Inf,
        handoff_schema::AbstractString="none",
        notes::AbstractString="",
    )
        isempty(process_id) && error("Physics process identifier cannot be empty.")
        isempty(domain) && error("Physics domain cannot be empty.")
        isempty(incoming_particles) && error("At least one incoming particle is required.")
        isempty(production_owner) && error("Every physics record requires a production owner.")
        status in PHYSICS_COVERAGE_STATUSES || error("Unknown physics coverage status: $(status).")
        lower = Float64(energy_min_eV)
        upper = Float64(energy_max_eV)
        isfinite(lower) && lower ≥ 0.0 || error("Physics lower energy must be finite and nonnegative.")
        (isfinite(upper) || isinf(upper)) && upper > lower || error(
            "Physics upper energy must exceed its lower energy.",
        )
        return new(
            String(process_id),String(domain),String[string(value) for value in incoming_particles],
            String[string(value) for value in outgoing_particles],String(production_owner),status,
            lower,upper,String(handoff_schema),String(notes),
        )
    end
end

struct Physics_Coverage_Register
    schema::String
    records::Vector{Physics_Coverage_Record}

    function Physics_Coverage_Register(
        records::AbstractVector{Physics_Coverage_Record};
        schema::AbstractString="radiant.hts.physics_coverage/v1",
    )
        isempty(records) && error("Physics coverage register cannot be empty.")
        record_vector = Physics_Coverage_Record[records...]
        keys = [(record.domain,record.process_id) for record in record_vector]
        length(unique(keys)) == length(keys) || error(
            "Physics domain/process keys must be unique.",
        )
        return new(String(schema),record_vector)
    end
end

function get_physics_record(
    this::Physics_Coverage_Register;
    domain::AbstractString,
    process_id::AbstractString,
)
    index = findfirst(
        record -> record.domain == domain && record.process_id == process_id,
        this.records,
    )
    isnothing(index) && error("No physics coverage record matches the requested key.")
    return this.records[index]
end

function validate_physics_coverage(this::Physics_Coverage_Register)
    for record in this.records
        if record.status == :implemented && record.production_owner != "Radiant"
            error("Implemented Radiant physics must identify Radiant as production owner.")
        end
        if record.status == :unsupported && record.handoff_schema == "none"
            error("Unsupported physics must name a fail-closed handoff or unavailable reason.")
        end
    end
    return true
end

"""Default production ownership for the high-fidelity HTS tape, fibre, and detector workflow."""
function default_hts_physics_coverage()
    records = Physics_Coverage_Record[
        Physics_Coverage_Record(
            process_id="photoatomic-coupled-em",
            domain="tape-and-diagnostic-microdomain",
            incoming_particles=["photon","electron","positron"],
            outgoing_particles=["photon","electron","positron"],
            production_owner="Radiant",
            status=:implemented,
            energy_min_eV=1.0e3,
            energy_max_eV=9.0e8,
            handoff_schema="radiant.multigroup_source/v1",
            notes="Photoelectric, Compton, Rayleigh, pair production, bremsstrahlung, impact ionization, annihilation, and atomic relaxation within the qualified data envelope.",
        ),
        Physics_Coverage_Record(
            process_id="neutron-transport-and-capture",
            domain="local-neutral-domain",
            incoming_particles=["neutron"],
            outgoing_particles=["neutron","photon","charged-recoil"],
            production_owner="OpenSn/OpenMC",
            status=:external_production,
            handoff_schema="magnet_boundary_phase_space/v2 + anisotropic_volume_source/v1",
            notes="Includes cryogenic low-energy sensitivity, capture, resonance self-shielding, and future GdBCO qualification.",
        ),
        Physics_Coverage_Record(
            process_id="photonuclear-reactions",
            domain="tape-and-diagnostic-microdomain",
            incoming_particles=["photon"],
            outgoing_particles=["neutron","proton","alpha","recoil","residual-nuclide"],
            production_owner="OpenMC/Geant4",
            status=:external_production,
            handoff_schema="reaction_event_kernel/v1",
            notes="Radiant EM attenuation must not be interpreted as photonuclear activation or recoil production.",
        ),
        Physics_Coverage_Record(
            process_id="nuclear-recoil-and-light-ion-transport",
            domain="atomistic-handoff-microdomain",
            incoming_particles=["recoil","proton","deuteron","triton","alpha","Li-7"],
            outgoing_particles=["recoil-tree","electronic-loss","nuclear-loss"],
            production_owner="SPECTRA-PKA/Geant4/BCA",
            status=:external_production,
            handoff_schema="weighted_atomistic_event_kernel/v1",
            notes="Includes converter reaction products and PKA cascades; mean EM deposition is not a substitute.",
        ),
        Physics_Coverage_Record(
            process_id="sub-keV-electron-track-structure",
            domain="selected-interface-microdomain",
            incoming_particles=["electron","positron","Auger-electron"],
            outgoing_particles=["ionization","excitation","phonon","escaping-electron"],
            production_owner="Geant4-supported-material reference",
            status=:external_production,
            energy_min_eV=0.0,
            energy_max_eV=1.0e3,
            handoff_schema="track_structure_handoff/v1",
            notes="YBCO-specific track-structure fidelity remains unsupported until compound data are implemented and validated.",
        ),
        Physics_Coverage_Record(
            process_id="activation-and-delayed-particles",
            domain="activation-regions",
            incoming_particles=["neutron","photon"],
            outgoing_particles=["decay-photon","beta-electron","positron","alpha","decay-recoil"],
            production_owner="OpenMC depletion/R2S",
            status=:external_production,
            handoff_schema="delayed_particle_source/v1",
            notes="Prompt and delayed source ledgers remain separate before Radiant or recoil-code transport.",
        ),
        Physics_Coverage_Record(
            process_id="fibre-optical-degradation",
            domain="fibre-diagnostic",
            incoming_particles=["deposited-energy","defect-state","optical-photon"],
            outgoing_particles=["RIA","FBG-shift","radioluminescence","Cherenkov-source"],
            production_owner="calibrated optical response model",
            status=:external_production,
            handoff_schema="fibre_optical_response/v1",
            notes="Radiant supplies radiation fields and optical-source terms, not calibrated colour-centre kinetics or wavelength drift.",
        ),
        Physics_Coverage_Record(
            process_id="hts-electrothermal-detector-response",
            domain="hts-detector",
            incoming_particles=["energy-deposition-impulse","converter-products"],
            outgoing_particles=["temperature-pulse","resistance-pulse","current-pulse"],
            production_owner="electrothermal response solver",
            status=:external_production,
            handoff_schema="hts_detector_impulse/v1",
            notes="Requires heat capacities, interface conductances, bias, transition curve, magnetic field, and operating temperature.",
        ),
    ]
    register = Physics_Coverage_Register(records)
    validate_physics_coverage(register)
    return register
end
