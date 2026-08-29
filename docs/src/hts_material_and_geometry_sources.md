# HTS material, cascade, CAD, and atomistic response sources

This document is part of the temporary `RadiantHTS` implementation on branch
`codex/radianthts-core-foundation-20260828`. It records the source and qualification boundary for
numerical material data. It is not an assertion that all listed literature is transferable to the
specific YBCO/GdBCO tape, oxygen state, texture, magnetic field, damage state, or temperature used
in a production calculation.

## Classification

Every numerical response must be classified as one of:

```text
literature-fit
literature-table-pending
atomistic-candidate
qualified
```

`literature-fit` means that an explicit published equation and coefficients have been transcribed
with their stated temperature range and uncertainty. `literature-table-pending` means that a
reference has been identified but no numerical curve has been silently invented. `qualified`
requires a hash-bound table, material state, uncertainty, and independent comparison.

## Cryogenic tape constituents

### OFHC copper

The numerical registry contains NIST Monograph 177 fits for:

- thermal conductivity at RRR 50, 100, 150, 300, and 500;
- specific heat.

Source:

- https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm

The selected RRR must match the stabilizer material. Irradiation, cold work, plating, soldering,
and magnetoresistance can change the effective conductivity; the NIST unirradiated fit is therefore
a candidate baseline rather than a universal tape value.

### Kapton and G-10 CR

The registry contains official NIST log-polynomial fits for Kapton thermal conductivity and
specific heat, and G-10 CR normal/warp conductivity plus specific heat.

Sources:

- https://trc.nist.gov/cryogenics/materials/Polyimide%20Kapton/PolyimideKapton_rev.htm
- https://trc.nist.gov/cryogenics/materials/G-10%20CR%20Fiberglass%20Epoxy/G10CRFiberglassEpoxy_rev.htm

### Coated-conductor heat capacity

The following source reports specific heat and thermal diffusivity for Ag/YBCO/Hastelloy and
Cu-reinforced coated conductors and shows that constituent-layer summation can reproduce the
measured conductor heat capacity:

- https://doi.org/10.1016/j.phpro.2012.06.316

The branch records this as `literature-table-pending`; the open curves still need digitization or a
layer-by-layer reconstruction with exact tape mass fractions.

### Hastelloy C-276

Cryogenic specific heat and thermal conductivity measurements are available in:

- https://doi.org/10.1063/1.2899058

The branch records the source and required digitization path, but does not substitute C-276 values
for every Hastelloy substrate grade without composition and specimen matching.

### REBCO tape longitudinal conductivity

A recent 4.2–200 K comparison of REBCO tapes with Cu, Ag, and Ag-Au stabilizers is available at:

- https://arxiv.org/abs/2410.09915

The stabilizer thickness and measured RRR must be retained with each digitized curve.

### YBCO and GdBCO tensors

A classic YBCO single-crystal anisotropy reference is:

- https://doi.org/10.1103/PhysRevB.40.9389

It is not automatically transferable to a textured coated-conductor film. No transferable GdBCO
thermal tensor has been selected. Both are therefore source-only records pending film-specific
literature extraction or atomistic generation.

## Gd capture cascades

The implementation distinguishes:

```text
capture gamma rays
characteristic x rays
internal-conversion electrons
Auger electrons
capture recoil
```

A gamma-only cascade repository is not a complete local heat source.

Primary sources and adapters:

- Zenodo event repository, DOI 10.5281/zenodo.7458654, containing reported large Gd-155 and Gd-157
  gamma-cascade samples. Use `scripts/data/download_zenodo_record.py` outside Git, then
  `summarize_gd_gamma_cascade_file` with the verified SHA-256.
- IAEA PGAA database: https://nucleus.iaea.org/Pages/pgaa-iaea.aspx
- Measured Gd-155/natural-Gd and Gd-157 spectra:
  https://doi.org/10.1093/ptep/ptz002 and https://doi.org/10.1093/ptep/ptaa015
- Geant4 NuDEX can provide an event-level nuclear-de-excitation route including internal
  conversion; it must be compared with measured/evaluated yields before production use.

The upstream neutron solve must supply isotope-resolved, spatially resolved, self-shielded capture
rates. The cascade adapter must not multiply an infinitely dilute microscopic cross section by an
unshielded scalar flux and call that a GdBCO capture source.

## CAD and faceted geometry

DAGMC/OpenMC CAD geometry uses MOAB H5M triangle facets and requires clean, watertight geometry.
The branch provides:

```text
scripts/geometry/export_dagmc_facets.py
src/hts_addon/Faceted_Geometry.jl
```

The Python adapter uses `pymoab` to preserve DAGMC geometric tags in a normalized HDF5 artifact.
Radiant independently checks degeneracy, duplicate triangles, boundary edges, nonmanifold edges,
and signed orientation. The current production route is:

```text
DAGMC/ParaStell H5M
  -> normalized faceted HDF5
  -> topology and overlap checks
  -> tape-surface local frames and source mapping
  -> piecewise-flat local domains or conservatively sampled Cartesian cells
  -> Radiant Cartesian SN sweep
```

This does not claim a native unstructured-mesh sweep. A direct tetrahedral/polyhedral transport
kernel would be a separate solver development and must not be hidden inside a stochastic
voxelization step.

## Generating YBCO/GdBCO response tables from DFT and MD

All generated tables use `radiant.hts.atomistic_response_table/v1` and must include exact material
state, structure, input, code, functional/potential, convergence settings, uncertainty, and source
hashes.

### Thermal-conductivity tensor and phonon lifetimes

Primary route:

```text
1. Enumerate YBCO6+x or GdBCO6+x oxygen ordering, texture, isotope, strain, and defect states.
2. Relax structures with converged DFT.
3. Generate second- and third-order force constants.
4. Run phono3py with q-mesh, supercell, interaction-cutoff, isotope-scattering, and boundary-length
   convergence.
5. Import `kappa-*.hdf5` with `read_phono3py_kappa_hdf5`.
6. Cross-check selected states with Green-Kubo heat-flux MD using a validated MLIP.
```

References:

- https://phonopy.github.io/phono3py/
- https://docs.lammps.org/compute_heat_flux.html

The MLIP route needs an ensemble or otherwise quantified model uncertainty. A potential that fits
bulk energies but is not stable for oxygen disorder, defects, and high-energy recoils is not
sufficient.

### Electron-phonon and quasiparticle relaxation

Primary route:

```text
1. Compute electronic and phonon states with converged DFT/DFPT.
2. Construct Wannier interpolation.
3. Run EPW for alpha2F, mode-resolved coupling, linewidths, and temperature dependence.
4. Import spectral output with `read_epw_spectral_function`.
5. Convert spectral functions to the detector thermalization model only through an explicit
   material-state-dependent reduction with uncertainty.
```

References:

- https://docs.epw-code.org/
- https://www.quantum-espresso.org/

Published YBCO ultrafast response values can be used as priors and cross-checks, not as universal
GdBCO/YBCO tape tables. The branch records DOI 10.1063/1.123388 as a source-only prior.

### Electronic stopping and sub-keV partition

Primary route:

```text
1. Run real-time TDDFT trajectories across crystallographic directions, impact parameters,
   projectile species, velocity, and charge states.
2. Converge cell size, time step, pseudopotential, exchange-correlation functional, and trajectory
   sampling.
3. Export stopping versus kinetic energy and direction.
4. Combine electronic excitation spectra with phonon/electron relaxation calculations to build
   ionization, excitation, prompt-phonon, athermal-phonon, quasiparticle, escape, and stored-energy
   partitions.
5. Import stopping tables with `read_electronic_stopping_table`; write complete response arrays with
   `write_atomistic_response_hdf5`.
```

Reference implementation option:

- https://octopus-code.org/

### Displacement thresholds, defect production, and stored energy

Primary route:

```text
1. Run directional DFT/MLIP recoil calculations for Y, Ba/Gd, Cu, and O sublattices.
2. Map threshold-displacement energy versus crystal direction and oxygen state.
3. Run cascade MD with an explicit electronic-stopping handoff and thermostat exclusion region.
4. Record prompt recombination, Frenkel pairs, antisites, clustering, amorphous volume, stored
   energy, strain, and escaping cascades.
5. Repeat with a potential ensemble and multiple seeds; preserve covariance among outputs.
```

Relevant YBCO radiation-damage modeling references include:

- https://doi.org/10.1088/1361-6668/ac47dc

The table pipeline records these as atomistic candidates until finite-size, time-window, potential,
and state uncertainties are quantified.

## Physical reference comparisons

`Physical_Reference_Qualification.jl` accepts hash-matched OpenMC, OpenSn, Geant4, or experimental
reference responses. A synthetic result cannot pass a physical gate. Required lineage includes:

```text
source artifact hash
result artifact hash
geometry hash
material-state hash
normalization hash
response units and uncertainty
```

For process-resolved scoring, the preferred references are matched continuous-energy or event-level
calculations with the same geometry/source and mutually exclusive score semantics. For curvature,
the piecewise-flat atlas must be compared to a genuinely curved OpenSn, Geant4, or explicit
Monte-Carlo segment.

## External data policy

Large H5M files, Gd event repositories, DFT/MD trajectories, force constants, and production HDF5
response tables remain outside Git. Only manifests, hashes, adapters, small analytic fixtures, and
qualification receipts are committed.
