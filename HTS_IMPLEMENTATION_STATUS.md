# Radiant HTS implementation and qualification status

Branch: `codex/radianthts-core-foundation-20260828`

Base: `e95eab60716a1a4b010b3c8669705809dc9065cb` (`v1.1.67`)

Draft CI pull request: `#1`

The pull request is a non-merge qualification trigger. The branch is not approved for `main`, and
physical production gates remain runtime-controlled.

## Implemented

### Cross-section and build-state repairs

- Corrected the `Cross_Sections` constructor so `particles::Vector{Particle}` is initialized with
  `Particle[]` rather than `Vector{String}()`.
- Initialized custom cross-section arrays and made setters replacement/idempotency safe.
- Added scalar and vector setters compatible with Julia 1.6.7 and 1.10.12.
- Corrected analytic custom absorption: `Σa` is now the supplied absorption and
  `Σt = Σa + Σs`.
- Added strict energy-grid, material, particle, dimension, and finite-value validation.
- Made `Source.set_particle` accept only a solver-compatible `Particle` object.

### Source and normalization contracts

- Explicit per-history, per-source-particle, or per-second normalization.
- Separate physical source rate and symmetry factor.
- Prompt/delayed time class and time interval.
- Source hash and provenance.
- Patch/group/direction boundary angular-current source.
- Isotropic, ordinate, and angular-moment volume sources.
- Optional source variances.
- Explicit rejection of negative source densities and unresolved nonzero grazing-current bins.

### Executable source adapters

- Exact source-to-Radiant energy-group interval mapping.
- Exact source-to-SN ordinate and quadrature-weight mapping.
- Cartesian boundary patch to boundary-cell-face mapping.
- Half-range boundary-moment projection.
- Volume-moment projection.
- Current/rate closure receipts.
- Negative reconstructed angular-source rejection without clipping.
- Boundary values installed as angular-flux densities; face measure is used only in current
  integration.

### Plain HDF5 interchange

- `radiant.boundary_angular_current/v1`.
- `radiant.anisotropic_volume_source/v1`.
- Hash-bound artifact reads.
- Explicit external row-major versus Julia-native storage declaration.
- Python/OpenMC/OpenSn-compatible dimension reversal.
- Boundary current and volume rate closure on read.
- Source normalization, time state, parent reaction, uncertainty, and provenance.

### Fixed-source lifecycle

- Explicit and legacy normalization contracts cannot be mixed.
- Compatible explicit sources may be combined.
- Source building is idempotent.
- Rebuilding after source addition starts from empty derived state.
- Projection receipts remain attached to the fixed-source collection.

### Executable HTS tape geometry

- `Tape_Geometry_Map` and `build_planar_tape_geometry`.
- Every physical layer is a separate region.
- Layer-requested transport subdivisions are preserved and can be uniformly refined.
- One-dimensional thickness, two-dimensional finite-width, and three-dimensional local-segment
  geometries.
- Activation and atomistic voxel masks.
- Layer sum, hotspot maximum, and volume-average response reductions.
- Canonical 4 mm wide, 163.2 micrometre eight-layer verification fixture.
- The straight planar model is explicit; curvature remains a separately qualified atlas/mapped-
  coordinate problem.

### High-fidelity diagnostic and detector contracts

- Fibre core, cladding, coating, adhesive, and capillary radial definition.
- Fibre radiation responses: dose, neutron capture, displacement source, charged-secondary
  spectrum, gas production, and optical-source/degradation handoffs.
- HTS detector active film, substrate, converter layers, bias, temperature, and electrothermal
  parameter contract.
- B-10-enriched B4C/YBCO transition-edge detector verification schema.
- Protected tape, fibre, and detector response registry.

### Deposition and physics ownership

- Mechanism-resolved `Energy_Partition`: ionization, excitation, prompt lattice heat, recoil
  handoff, defect-stored energy, optical emission, escape, cutoff, and unresolved terms.
- Separate owners for photoatomic EM, neutrons/capture/self-shielding, photonuclear and
  electronuclear reactions, recoils/light ions, sub-keV track structure, activation, optical
  response, and electrothermal/superconducting response.
- Photonuclear activation is not assigned to native OpenMC photoatomic transport.
- `parastell_damage.transport_ownership_map/v1` enforces one additive owner in production and
  non-additive comparison mode.

### Qualification infrastructure

- Pinned Julia 1.6.7 and 1.10.12 matrix.
- Per-version TOML receipts with git, Julia, project, and manifest hashes.
- Analytic one-cell photon/DD energy-deposition qualification.
- OpenSn-shaped boundary replay with current-conserving projection.
- HDF5 boundary and volume round-trip tests.
- Explicit eight-layer geometry and layer-response tests.
- Physical replay requirements and hash-bound case manifest requirements.

## Runtime gate status

| Gate | Status |
|---|---|
| Branch isolated from `main` | PASS |
| Draft PR is non-merge qualification trigger | PASS |
| Source/HDF5/tape/detector code authored | PASS (static) |
| Regression and qualification tests authored | PASS |
| Julia package compilation | NOT EXECUTED |
| Julia 1.6.7 qualification | NOT EXECUTED |
| Julia 1.10.12 qualification | NOT EXECUTED |
| Analytic EM fixture | NOT EXECUTED |
| OpenSn-shaped interface fixture | NOT EXECUTED |
| Physical OpenMC/OpenSn source replay | BLOCKED — artifact absent |
| `RADIANT_BUILD_PASS` | NOT ESTABLISHED |
| `RADIANT_EM_ANALYTIC_PASS` | NOT ESTABLISHED |
| `RADIANT_EM_PASS` | NOT ESTABLISHED |
| `OPENSN_RADIANT_INTERFACE_SOFTWARE_PASS` | NOT ESTABLISHED |
| `OPENSN_RADIANT_COUPLING_PASS` | NOT ESTABLISHED |

GitHub currently reports no workflow run or commit-status check for the draft PR head. Connector-
authored branch and PR events have therefore not produced executable CI evidence. No runtime pass is
claimed.

## Physical gate distinction

The analytic EM fixture can establish only `RADIANT_EM_ANALYTIC_PASS`. A physical
`RADIANT_EM_PASS` additionally requires compound material data, matched continuous-energy or
Geant4 references, layer/interface responses, and energy/particle/charge closure.

The OpenSn-shaped source fixture can establish only `OPENSN_RADIANT_INTERFACE_SOFTWARE_PASS`. A
physical `OPENSN_RADIANT_COUPLING_PASS` requires a hash-matched OpenMC/OpenSn bundle, HDF5 import,
projection, Radiant solve, outgoing-current retention, and return-current replay where photon
ownership is decomposed across domains.

## Next runtime actions

1. Enable/authorize GitHub Actions for the fork or run the pinned matrix locally.
2. Correct every Julia 1.6.7 and 1.10.12 compile/test failure from the receipts.
3. Add native outgoing boundary angular-current retention and balance extraction.
4. Replay a physical OpenMC boundary source and OpenSn neutron-induced photon volume source.
5. Compare straight YBCO tape and diagnostic/detector fixtures against matched OpenMC/Geant4
   references.
6. Add the piecewise-flat tape atlas and spatial magnetic-field qualification.
7. Qualify GdBCO capture/self-shielding against continuous-energy reference calculations.
