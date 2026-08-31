#!/usr/bin/env julia

using Radiant

length(ARGS) == 1 || error(
    "Usage: create_microdosimetry_event_bank_fixture.jl OUTPUT.h5",
)

fractions = Event_Energy_Partition_Fractions(Dict(
    :ionization => 0.20,
    :electronic_excitation => 0.10,
    :prompt_lattice_heat => 0.40,
    :nuclear_recoil_handoff => 0.10,
    :defect_storage => 0.05,
    :optical_emission => 0.05,
    :escaping_particle => 0.05,
    :cutoff_handoff => 0.05,
    :unresolved => 0.0,
))

secondary = Correlated_Secondary(
    "electron",25.0;
    direction_model=:parent_correlated,
    correlation_id="cross-language-event",
)
prototype = Microdosimetry_Event_Prototype(
    prototype_id="cross-language-fixed-parent",
    particle_tag="photon",
    process_id="compton",
    material_tag="YBCO",
    layer_id="REBCO",
    position_cm=[1.0,2.0,3.0],
    direction_model=:fixed,
    fixed_direction=[0.0,1.0,0.0],
    incident_energy_eV=1000.0,
    mean_deposited_energy_eV=100.0,
    event_rate_per_s=2.0,
    partition_fractions=fractions,
    correlated_secondaries=[secondary],
    correlation_id="cross-language-event",
    provenance=Dict("fixture" => "julia-to-python-hdf5-replay"),
)
kernel = Microdosimetry_Kernel(
    [prototype];
    source_artifact_hash="cross-language-source-hash",
    geometry_hash="cross-language-geometry-hash",
    material_state_hash="cross-language-material-hash",
)
events = sample_microdosimetry_events(kernel,4;time_window_s=(0.0,2.0),seed=17)
write_weighted_microdosimetry_event_bank_hdf5(
    ARGS[1],events;
    source_hash=repeat("a",64),
    kernel_hash=repeat("b",64),
)
