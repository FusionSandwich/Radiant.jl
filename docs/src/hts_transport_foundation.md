# HTS deterministic-transport foundation

This branch begins the in-core HTS extension for coupling reactor-scale neutral transport to
layer-resolved Radiant photon/electron/positron calculations. The implementation is deliberately
conservative: source normalization, phase-space mapping, and code ownership are explicit and
validated before additional geometry or activation physics is added.

## Explicit normalization

`Source_Normalization` stores the source basis, physical source rate, symmetry factor, prompt or
delayed classification, time interval, producer hash, and provenance. Source tensors remain in
their declared basis; physical scaling is applied only through `apply_normalization`.

New explicit sources are not mixed with legacy `Surface_Source` or `Volume_Source` objects in one
`Fixed_Sources` collection. Their normalization semantics are different. For explicit sources the
legacy transport divisor is fixed at one so the solved response remains per history, per source
particle, or per second exactly as declared.

`Fixed_Sources.build` is now idempotent. Repeated builds do not duplicate source arrays or
normalization, and adding a new source invalidates the previous build.

## Boundary angular current

`Boundary_Angular_Current_Source` stores patch-resolved geometry and local frames together with
energy boundaries, ordinates, quadrature weights, angular flux, uncertainty, and provenance. The
protected invariant is

```math
J^-_{k,g} = A_k \sum_{m:\Omega_m\cdot n_k<0}
                  w_m |\Omega_m\cdot n_k|\psi^-_{k,g,m}.
```

`boundary_source_from_directional_current` converts directional-current contributions to angular
flux and immediately checks closure. It rejects unresolved nonzero grazing bins rather than
clipping them or dividing by an arbitrarily small direction cosine.

`project_boundary_source` now connects this object to the existing SN sweep. The initial qualified
mapping is intentionally narrow:

- exact source-to-Radiant energy-group interval matching;
- source directions and weights must match a subset of the selected Radiant quadrature;
- one source patch maps to one Cartesian boundary-cell face;
- patch position, outward normal, and measure must match that face;
- half-range reconstruction must remain nonnegative;
- patch/group current must close after angular projection.

The moments installed in the sweep remain angular-flux densities. Patch measure appears only in
the extensive-current closure receipt; it is not folded into the transport boundary value.

## Anisotropic volume sources

`Anisotropic_Volume_Source` supports isotropic, ordinate-collocated, and externally declared
moment representations. It is intended for neutron-induced photons, decay photons, delayed
charged particles, and other coupled-particle sources.

`project_volume_source` maps these sources into the existing SN volume-source tensor while
protecting voxel volume, scalar source rate, energy-group identity, and angular positivity. Signed
higher angular moments are allowed, but a native moment source must explicitly declare the
Radiant moment convention and its angle-integrated zeroth coefficient.

Both boundary and volume projections produce hash-bound receipts containing the source-to-target
maps and closure errors. `Fixed_Sources` retains those receipts for downstream audit.

## Balance and ownership

`Transport_Balance` records particle, kinetic-energy, rest-mass-exchange, cutoff, leakage,
deposition, and charge terms. It is currently a validated ledger type; native extraction of every
term from the sweep is still pending.

`Transport_Ownership_Map` implements
`parastell_damage.transport_ownership_map/v1` and distinguishes:

- production mode: exactly one additive owner for every domain/particle/source/response key;
- comparison mode: at least two non-additive comparison solvers and no production owner.

This is the core guard against adding overlapping OpenSn and Radiant photon heating.

## HTS tape definitions

`Layer_Definition` and `Tape_Stack_Definition` represent ordered physical layers and identify
activation and atomistic-output regions. `verification_eight_layer_stack()` constructs the fixed
4 mm wide, 163.2 micrometre verification fixture:

1. 20 micrometre front copper;
2. 2 micrometre silver cap;
3. 1 micrometre REBCO;
4. 0.2 micrometre buffer;
5. 50 micrometre Hastelloy;
6. 20 micrometre rear copper;
7. 20 micrometre solder;
8. 50 micrometre insulation.

The fixture is marked `verification-only`; it is not a measured production tape. The REBCO
material tag can be changed to `GdBCO` for schema testing, but this does not qualify Gd neutron
self-shielding or capture-gamma physics.

## Regression coverage added in this branch

The former placeholder test has been replaced with tests for:

- normalization bases, physical scale, time class, and invalid inputs;
- boundary-current integration and exact directional-current conversion;
- rejection of negative and unresolved grazing sources;
- isotropic and signed-moment volume-source contracts;
- one-voxel boundary and volume projection into Radiant arrays;
- projection receipts and closure;
- explicit-source build idempotency;
- rejection of mixed legacy/explicit normalization;
- the eight-layer HTS verification fixture;
- no-double-counting production and comparison maps;
- particle, energy, and charge ledger residuals.

## Remaining implementation sequence

1. execute and qualify the new test suite on pinned Julia versions;
2. add versioned cross-language HDF5 source interchange;
3. retain outgoing boundary angular moments and currents in `Flux` results;
4. compute native particle, energy, charge, and cutoff balances from the sweep;
5. add layer-, process-, and reaction-resolved response objects;
6. add OpenMC and OpenSn source/receipt adapters;
7. add piecewise-flat tape-atlas frames and curvature convergence;
8. add delayed-source, activation-rate, PKA, and atomistic event-kernel exports;
9. qualify GdBCO self-shielding against continuous-energy references.
