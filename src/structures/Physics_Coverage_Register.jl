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
            String(process_id),String(domain),
            String[string(value) for value in incoming_particles],
            String[string(value) for value in outgoing_particles],
            String(production_owner),status,lower,upper,String(handoff_schema),String(notes),
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
        if occursin("photonuclear",record.process_id) && occursin("OpenMC",record.production_owner)
            error("Native OpenMC photoatomic transport must not be assigned photonuclear ownership.")
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
            notes="Photoelectric, Compton, Rayleigh, pair/triplet production, bremsstrahlung, impact ionization, annihilation, and atomic relaxation within the qualified data envelope.",
        ),
        Physics_Coverage_Record(
            process_id="neutron-transport-capture-and-self-shielding",
            domain="local-neutral-domain",
            incoming_particles=["neutron"],
            outgoing_particles=["neutron","capture-photon","conversion-electron","Auger-electron","charged-recoil"],
            production_owner="OpenSn/OpenMC",
            status=:external_production,
            handoff_schema="magnet_boundary_phase_space/v2 + anisotropic_volume_source/v1",
            notes="Includes low-energy and resonance transport, temperature-dependent nuclear data, capture cascades, and thickness-dependent GdBCO self-shielding qualification.",
        ),
        Physics_Coverage_Record(
            process_id="photonuclear-reactions",
            domain="tape-and-diagnostic-microdomain",
            incoming_particles=["photon"],
            outgoing_particles=["neutron","proton","deuteron","triton","alpha","recoil","residual-nuclide"],
            production_owner="Geant4/MCNP/PHITS",
            status=:external_production,
            handoff_schema="reaction_event_kernel/v1",
            notes="Radiant and OpenMC photoatomic attenuation must not be interpreted as photonuclear activation, recoil, or photoneutron production.",
        ),
        Physics_Coverage_Record(
            process_id="electronuclear-and-positronuclear-reactions",
            domain="selected-high-energy-microdomain",
            incoming_particles=["electron","positron"],
            outgoing_particles=["neutron","proton","meson","recoil","residual-nuclide"],
            production_owner="Geant4",
            status=:external_production,
            handoff_schema="reaction_event_kernel/v1",
            notes="Included only when the incident lepton spectrum reaches a material reaction threshold; it is not part of Radiant's condensed electromagnetic deposition score.",
        ),
        Physics_Coverage_Record(
            process_id="nuclear-recoil-and-light-ion-transport",
            domain="atomistic-handoff-microdomain",
            incoming_particles=["recoil","proton","deuteron","triton","alpha","Li-7"],
            outgoing_particles=["recoil-tree","electronic-loss","nuclear-loss"],
            production_owner="SPECTRA-PKA/Geant4/BCA",
            status=:external_production,
            handoff_schema="weighted_atomistic_event_kernel/v1",
            notes="Includes neutron/photon reaction products, converter products, decay recoils, and PKA cascades; mean electromagnetic deposition is not a substitute.",
        ),
        Physics_Coverage_Record(
            process_id="sub-keV-electron-track-structure",
            domain="selected-interface-microdomain",
            incoming_particles=["electron","positron","Auger-electron"],
            outgoing_particles=["ionization","excitation","phonon","escaping-electron"],
            production_owner="validated track-structure reference",
            status=:external_production,
            energy_min_eV=0.0,
            energy_max_eV=1.0e3,
            handoff_schema="track_structure_handoff/v1",
            notes="YBCO/GdBCO compound track-structure fidelity is not claimed until material-specific low-energy data and interface escape are validated.",
        ),
        Physics_Coverage_Record(
            process_id="neutron-activation-and-delayed-particles",
            domain="activation-regions",
            incoming_particles=["neutron"],
            outgoing_particles=["decay-photon","beta-electron","positron","alpha","decay-recoil"],
            production_owner="OpenMC depletion/R2S",
            status=:external_production,
            handoff_schema="delayed_particle_source/v1",
            notes="OpenMC supplies neutron reaction rates, inventory evolution, and decay-photon sources; non-photon delayed particles require an inventory-to-source adapter.",
        ),
        Physics_Coverage_Record(
            process_id="photonuclear-activation-and-delayed-particles",
            domain="activation-regions",
            incoming_particles=["photon","electron","positron"],
            outgoing_particles=["residual-nuclide","decay-photon","beta-electron","positron","alpha","decay-recoil"],
            production_owner="Geant4/MCNP/PHITS rates + FISPACT-II/ALARA or verified inventory solver",
            status=:external_production,
            handoff_schema="photonuclear_reaction_rates/v1 + delayed_particle_source/v1",
            notes="Kept separate from native OpenMC neutron depletion because photo/electro-nuclear rates and residual products require their own qualified producer.",
        ),
        Physics_Coverage_Record(
            process_id="energy-to-heat-and-damage-partition",
            domain="all-explicit-layers",
            incoming_particles=["deposited-energy","nuclear-recoil-energy"],
            outgoing_particles=["prompt-lattice-heat","ionization","electronic-excitation","defect-stored-energy","optical-emission","sub-cutoff-handoff"],
            production_owner="response-specific energy-partition model",
            status=:external_production,
            handoff_schema="energy_partition/v1",
            notes="Transport energy deposition is not automatically instantaneous cryogenic heat; the partition must close energy and preserve defect/optical/cutoff ledgers.",
        ),
        Physics_Coverage_Record(
            process_id="superconducting-nonequilibrium-response",
            domain="HTS-active-layer",
            incoming_particles=["ionization","electronic-excitation","phonon","local-temperature"],
            outgoing_particles=["quasiparticle-state","pair-breaking","vortex-response","Jc-shift","Tc-shift","resistivity-change"],
            production_owner="HTS electrothermal/superconductivity response model",
            status=:external_production,
            handoff_schema="hts_nonequilibrium_response/v1",
            notes="Required for detector pulses and transient magnet response; it is distinct from long-lived displacement damage and from bulk thermal diffusion.",
        ),
        Physics_Coverage_Record(
            process_id="spatial-magnetic-field-charged-transport",
            domain="curved-tape-and-detector-domain",
            incoming_particles=["electron","positron","charged-ion"],
            outgoing_particles=["field-deflected-transport"],
            production_owner="Geant4 field-map reference; Radiant after field-map qualification",
            status=:comparison_only,
            handoff_schema="spatial_field_map/v1",
            notes="Radiant currently exposes a constant magnetic field. A curved stellarator tape requires a spatially varying field and local-frame convergence study before production ownership changes.",
        ),
        Physics_Coverage_Record(
            process_id="cryogenic-material-state-feedback",
            domain="all-explicit-layers",
            incoming_particles=["temperature","magnetic-field","strain","oxygen-content","irradiation-state"],
            outgoing_particles=["density","cross-section-state","thermal-properties","electrical-properties"],
            production_owner="coupled material-state controller",
            status=:external_production,
            handoff_schema="material_state/v1",
            notes="Temperature, field, strain, oxygenation, and evolving inventory must remain explicit rather than being hidden in a material tag.",
        ),
        Physics_Coverage_Record(
            process_id="fibre-optical-source-generation",
            domain="fibre-diagnostic",
            incoming_particles=["charged-particle-track","electronic-excitation"],
            outgoing_particles=["Cherenkov-photon","radioluminescence-photon","scintillation-photon"],
            production_owner="Geant4 optical or calibrated optical-source model",
            status=:external_production,
            handoff_schema="fibre_optical_source/v1",
            notes="Radiant may supply charged-particle and excitation source terms, but optical-photon production and trapping require refractive-index, emission, and interface data.",
        ),
        Physics_Coverage_Record(
            process_id="fibre-optical-degradation",
            domain="fibre-diagnostic",
            incoming_particles=["dose","fluence","defect-state","temperature","dose-rate","annealing-history"],
            outgoing_particles=["RIA","FBG-shift","radioluminescence","recovery-state"],
            production_owner="calibrated optical response model",
            status=:external_production,
            handoff_schema="fibre_optical_response/v1",
            notes="Radiant supplies radiation fields and source terms, not calibrated colour-centre kinetics, RIA, or wavelength drift.",
        ),
        Physics_Coverage_Record(
            process_id="hts-electrothermal-detector-response",
            domain="hts-detector",
            incoming_particles=["energy-deposition-impulse","converter-products","bias-state"],
            outgoing_particles=["temperature-pulse","resistance-pulse","current-pulse","pulse-height","decay-time"],
            production_owner="electrothermal response solver",
            status=:external_production,
            handoff_schema="hts_detector_impulse/v1",
            notes="Requires active/substrate heat capacities, interface conductances, bias, transition curve, magnetic field, operating temperature, and readout response.",
        ),
    ]
    register = Physics_Coverage_Register(records)
    validate_physics_coverage(register)
    return register
end
