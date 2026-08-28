# HTS deterministic-transport foundation

This module couples reactor-scale neutral transport to layer-resolved Radiant calculations while
keeping source normalization, physics ownership, heat partitioning, diagnostics, and detector
response explicit.

## Source normalization

`Source_Normalization` stores the source basis, physical source rate, symmetry factor, prompt or
delayed classification, time interval, source hash, and provenance. Source tensors remain in their
declared basis; physical scaling is applied only through `apply_normalization`.

## Boundary angular current

`Boundary_Angular_Current_Source` stores patch-resolved geometry and local frames together with
energy boundaries, ordinates, quadrature weights, angular flux, uncertainty, and provenance. The
protected invariant is

```math
J^-_{k,g} = A_k \sum_{m:\Omega_m\cdot n_k<0}
                  w_m |\Omega_m\cdot n_k|\psi^-_{k,g,m}.
```

`boundary_source_from_directional_current` converts directional current-density contributions to
angular flux and immediately checks closure. It rejects unresolved nonzero grazing bins rather than
clipping them or dividing by an arbitrarily small direction cosine.

`project_boundary_source` maps exact source energy intervals and ordinates into Radiant's
half-range SN boundary basis. Each current implementation patch must match one Cartesian boundary-
cell face. Current-density and extensive patch-current closure are reported in a
`Boundary_Projection_Receipt`.

## Anisotropic volume sources

`Anisotropic_Volume_Source` supports isotropic, ordinate-collocated, and explicitly declared moment
representations for neutron-induced photons, decay photons, delayed charged particles, and other
coupled-particle sources. `project_volume_source` preserves integrated source rate and rejects
negative reconstructed angular-source lobes.

## Plain HDF5 interchange

The production interchange is plain HDF5 rather than Julia-only serialization:

- `radiant.boundary_angular_current/v1`;
- `radiant.anisotropic_volume_source/v1`.

The reader verifies an optional file SHA-256, schema, particle identity, source normalization,
storage order, array dimensions, current/rate ledger, uncertainty arrays, and provenance.
Multidimensional arrays written for external use are dimension-reversed so Python/h5py and other
row-major readers see canonical patch/group/direction and voxel/group/coefficient shapes.

## Executable tape geometry

`build_planar_tape_geometry` turns a `Tape_Stack_Definition` into a built Cartesian Radiant
geometry. The x-axis is the physical layer-normal direction; y is tape width; z is local segment
length. Every physical layer is a separate material region and retains its requested transport
subdivisions. The returned `Tape_Geometry_Map` provides layer voxel ranges, activation/atomistic
voxel IDs, and layer sum, maximum, and volume-average response reductions.

`verification_eight_layer_stack()` constructs the fixed 4 mm wide, 163.2 micrometre verification
fixture:

1. 20 micrometre front copper;
2. 2 micrometre silver cap;
3. 1 micrometre REBCO;
4. 0.2 micrometre buffer;
5. 50 micrometre Hastelloy;
6. 20 micrometre rear copper;
7. 20 micrometre solder;
8. 50 micrometre insulation.

The fixture is marked `verification-only`. The REBCO material tag can be changed to `GdBCO` for
schema testing, but no Gd self-shielding qualification is implied. The geometry is straight; curved
winding use requires a piecewise-flat or mapped-coordinate qualification.

## Heat and damage partition

`Energy_Partition` prevents a transport energy-deposition score from being silently interpreted as
instantaneous cryogenic heat. It separately accounts for ionization, electronic excitation, prompt
lattice heat, recoil-code handoff, metastable defect-stored energy, optical emission, escaping
energy, sub-cutoff handoff, and unresolved energy. The ledger must close before a response is
accepted.

## Fibre diagnostics

`Fibre_Diagnostic_Definition` represents concentric core, cladding, coating, adhesive, and capillary
layers. It protects radiation responses such as core dose, neutron capture, displacement source,
secondary charged-particle spectrum, gas production, and optical-source terms. Radiation-induced
attenuation, colour-centre kinetics, FBG wavelength shift, and recovery remain calibrated optical
response models rather than being inferred from dose alone.

## HTS detector use

`HTS_Detector_Definition` separates converter reactions and transport from the active-film
superconducting/electrothermal response. It stores active HTS dimensions, substrate, temperature,
bias, converter layers, and an optional thermal/electrical model. The B-10-enriched B4C/YBCO
verification detector is a schema fixture, not a performance prediction.

## Balance, ownership, and physics coverage

`Transport_Balance` records particle, kinetic-energy, rest-mass-exchange, cutoff, leakage,
deposition, and charge terms. `Transport_Ownership_Map` distinguishes production from non-additive
comparison calculations. `Physics_Coverage_Register` identifies the production owner and handoff
for photoatomic EM, neutrons, photonuclear/electronuclear reactions, nuclear recoils, activation,
sub-keV track structure, optical response, and HTS electrothermal/superconducting response.

## Qualification hierarchy

The pinned qualification runner targets Julia 1.6.7 and 1.10.12 and writes TOML receipts. It
contains:

- package instantiate, precompile, import, and tests;
- HDF5 boundary/volume round trips;
- explicit tape-geometry tests;
- analytic one-cell photon/DD flux and deposition;
- an OpenSn-shaped boundary-source replay.

The analytic and synthetic fixtures are software verification only. They do not establish physical
`RADIANT_EM_PASS` or `OPENSN_RADIANT_COUPLING_PASS`. Those gates require a hash-matched OpenMC/
OpenSn source/result bundle and matched physical-material reference calculations.
