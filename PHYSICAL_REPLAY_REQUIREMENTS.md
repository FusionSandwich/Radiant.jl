# Physical OpenMC–OpenSn–Radiant replay requirements

A physical replay is accepted only when all inputs are hash-bound to the same geometry, material,
source-normalization, energy-group, angular, and ownership contracts. Analytic fixtures and public
examples remain software verification and cannot establish a physical production gate.

## Required boundary source artifact

Schema: `radiant.boundary_angular_current/v1`

Required HDF5 groups and datasets:

```text
/meta/schema
/meta/particle
/meta/source_representation
/meta/storage_order
/meta/normalization_basis
/meta/source_rate_per_s
/meta/symmetry_factor
/meta/time_interval_s
/meta/time_class
/meta/source_hash
/patch/id
/patch/centroid_cm
/patch/area_cm2
/patch/normal
/patch/tangent_1
/patch/tangent_2
/energy/edges_eV
/angle/direction_cosines
/angle/quadrature_weights
/source/angular_flux
    or /source/directional_current_density
/source/incoming_current
/provenance/keys
/provenance/values
```

The file must preserve patch-local position, outward normal, two tangent vectors, energy, direction,
quadrature/current weighting, source uncertainty, and normalization. A surface-averaged scalar
spectrum is not sufficient for the physical coupling gate.

## Required OpenSn volume-source artifact

Schema: `radiant.anisotropic_volume_source/v1`

This is required for neutron-induced photons and any other OpenSn-owned internal production entering
Radiant. It contains voxel IDs and volumes, energy groups, angular representation, source values,
optional variance, parent reaction, integrated-rate ledger, source normalization, and provenance.

## Required matched case bundle

The physical case manifest must bind:

- global OpenMC source-bank SHA-256;
- OpenMC statepoint/tally SHA-256;
- OpenSn input and output SHA-256;
- Radiant source HDF5 SHA-256;
- exact local geometry and material SHA-256;
- cross-section library and converter SHA-256;
- energy group ordering and units;
- angular quadrature and basis convention;
- source rate, sector/full-device symmetry, and irradiation interval;
- prompt/delayed classification;
- transport ownership map SHA-256;
- protected-response registry SHA-256.

## Minimum physical replay matrix

1. Straight one-dimensional eight-layer verification stack.
2. Straight finite-width three-dimensional tape segment.
3. Normal, oblique, and grazing source cases.
4. Prompt external photon source.
5. Neutron-induced photon volume source.
6. Coupled photon/electron/positron transport.
7. YBCO layer and interface heating/current responses.
8. Selected fibre core/cladding/coating response.
9. Selected HTS detector converter/active-layer deposition response.
10. Continuous-energy OpenMC and/or Geant4 reference using the same source and geometry.

## Required gate evidence

`RADIANT_BUILD_PASS` requires successful pinned Julia 1.6.7 and 1.10.12 qualification receipts.

`RADIANT_EM_PASS` requires physical material data and matched reference comparisons for particle,
energy, charge, layer deposition, spectra, and leakage. The analytic DD fixture alone is insufficient.

`OPENSN_RADIANT_COUPLING_PASS` requires current/rate closure after HDF5 import, source projection,
Radiant solve, outgoing-current extraction, and—when domains exchange photons—return-current replay.
A one-way source import is not a closed domain-decomposition qualification.

## Current block

No physical OpenMC/OpenSn phase-space or volume-source artifact is committed on this branch. Until a
matched bundle is supplied, only software-verification receipts can be produced. Missing artifacts
must be reported as unavailable; they must not be substituted by a public example or reconstructed
from a scalar spectrum.
