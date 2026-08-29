# Temporary HTS add-on implementation

This directory contains HTS-specific functionality that is **temporarily compiled as part of
Radiant.jl** on branch `codex/radianthts-core-foundation-20260828`.

It is intended to become a separate package, provisionally named `RadiantHTS.jl`. It is not part
of the proposed upstream Radiant core API, and it must not be merged into `main` merely because the
verification tests pass.

## Repository policy

```text
current host: FusionSandwich/Radiant.jl
current branch: codex/radianthts-core-foundation-20260828
future package: RadiantHTS.jl
pull request: none
main merge: prohibited until an explicit future extraction/integration decision
force push: prohibited
```

`Addon_Extraction_Manifest.jl` is the machine-readable ownership contract. Every HTS-specific
implementation file is owned by one extraction component. Changes outside this directory must be
limited to generic hooks, include/export wiring, standard-library dependencies, or qualification
infrastructure.

## Current components

| Component | File | Current scope |
|---|---|---|
| Process-resolved scoring | `Process_Resolved_Scoring.jl` | Fold native interaction response channels with Radiant flux |
| Spatial field maps | `Spatial_Magnetic_Field_Map.jl` | Constant and Cartesian mapped magnetic fields, local-frame sampling |
| Piecewise-flat atlas | `Piecewise_Flat_Tape_Atlas.jl` | Parallel-transport frames, curvature/torsion/sagitta screening, local field |
| Gd capture cascades | `Gd_Prompt_Capture_Cascade.jl` | Isotope-resolved external evaluated cascades, photons/electrons/recoils |
| Sub-keV thermalization | `SubkeV_Thermalization.jl` | Calibrated energy partition and non-equilibrium reservoirs |
| Statistical microdosimetry | `Statistical_Microdosimetry.jl` | Weighted event prototypes, straggling, correlated secondaries |
| Cryogenic electrothermal response | `Cryogenic_Electrothermal.jl` | Implicit thermal network and HTS transition-resistance response |

## Generic Radiant hooks

Only these changes are intended to remain candidates for Radiant core:

1. `Multigroup_Cross_Sections.response_channels`, a generic process-response store.
2. Population of those channels during physical-model cross-section generation.
3. `src/Radiant.jl` include/export wiring while the add-on is temporarily co-hosted.
4. `Random` as a standard-library dependency while microdosimetry remains co-hosted.

At extraction, `RadiantHTS.jl` should use public Radiant APIs or a deliberately small extension API.
No HTS material, detector, tape-atlas, Gd cascade, or cryogenic response policy should remain hidden
inside generic Radiant transport kernels.

## Qualification classification

The following distinctions are mandatory:

```text
implemented-and-unit-tested
software-verified-with-synthetic-data
physical-data-candidate
physical-data-qualified
production-ready
```

A synthetic Gd cascade, sub-keV kernel, field map, microdosimetry kernel, or electrothermal model is
never physical qualification. A process response score is not instantaneous heat until a qualified
energy-partition and thermalization model assigns it to the heat channel.

## Ownership boundaries

```text
OpenMC / OpenSn:
  neutron transport, isotope-resolved capture rates, self-shielding,
  neutron-to-photon production, activation input

Radiant:
  photon/electron/positron transport inside the selected EM microdomain

RadiantHTS temporary add-on:
  local geometry/frame mapping, process scoring, radiation-to-response handoffs,
  weighted microdosimetry, cryogenic electrothermal response

SPECTRA-PKA / BCA / Geant4 / MD:
  nuclear recoil and light-ion event physics

OpenMC depletion / activation tools:
  delayed inventories and delayed particle sources
```

Overlapping photon heating from OpenSn and Radiant may be compared but never added. Prompt and
delayed sources remain separate through scoring and electrothermal application.

## Remaining production gates

- Physical process-channel closure against established Radiant total responses.
- Spatial field-map transport comparison against selected fully resolved calculations.
- Piecewise-flat response convergence and a curved OpenSn/Geant4 reference segment.
- Evaluated Gd-155 and Gd-157 cascade data with self-shielded capture-rate inputs.
- Material-specific YBCO/GdBCO sub-keV and non-equilibrium calibration.
- Detector microdosimetry validation against event-level reference simulations or experiment.
- Cryogenic thermal-property, interface-conductance, superconducting-transition, and readout
  calibration.
- Final extraction to `RadiantHTS.jl` before any upstream or default-branch integration decision.
