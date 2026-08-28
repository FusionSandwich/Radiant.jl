# HTS deterministic-transport foundation

This module is the first core implementation for coupling reactor-scale neutral transport to
layer-resolved Radiant calculations. It defines conservative source and ownership contracts before
those contracts are wired into the existing Cartesian sweep kernels.

## Explicit normalization

`Source_Normalization` stores the source basis, physical source rate, symmetry factor, prompt or
delayed classification, time interval, producer hash, and provenance. Source tensors remain in
their declared basis; physical scaling is applied only through `apply_normalization`.

## Boundary angular current

`Boundary_Angular_Current_Source` stores patch-resolved geometry and local frames together with
energy boundaries, ordinates, quadrature weights, angular flux, uncertainty, and provenance. The
protected invariant is

```math
J^-_{k,g} = A_k \sum_{m:\Omega_m\cdot n_k<0}
                  w_m |\Omega_m\cdot n_k|\psi^-_{k,g,m}.
```

`boundary_source_from_directional_current` converts directional current contributions to angular
flux and immediately checks closure. It rejects unresolved nonzero grazing bins rather than
clipping them or dividing by an arbitrarily small direction cosine.

## Anisotropic volume sources

`Anisotropic_Volume_Source` supports isotropic, ordinate-collocated, and externally declared
moment representations. It is intended for neutron-induced photons, decay photons, delayed
charged particles, and other coupled-particle sources. Moment sources cannot be silently reduced
to scalar rates without an explicit basis mapping.

## Balance and ownership

`Transport_Balance` records particle, kinetic-energy, rest-mass-exchange, cutoff, leakage,
deposition, and charge terms. `Transport_Ownership_Map` implements
`parastell_damage.transport_ownership_map/v1` and distinguishes:

- production mode: exactly one additive owner for every domain/particle/source/response key;
- comparison mode: at least two non-additive comparison solvers and no production owner.

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

The fixture is marked `verification-only`; it is not a production tape definition. The REBCO
material tag can be changed to `GdBCO` for schema testing, but no Gd self-shielding qualification is
implied.

## Current implementation boundary

The contracts and regression tests are implemented. They are not yet connected to the legacy
`Surface_Source`, `Volume_Source`, `Source`, or sweep result objects. The next implementation slice
will add:

1. HDF5 boundary and volume source interchange;
2. projection into Radiant half-range source moments;
3. outgoing boundary-current retention;
4. native particle/energy/charge balance extraction;
5. layer- and process-resolved response objects;
6. OpenSn and OpenMC adapters;
7. piecewise-flat tape-atlas frames.
