# HTS material-response data

This directory is temporary add-on data for the branch-only `RadiantHTS` work. It is not a claim
that generic Radiant contains qualified YBCO or GdBCO material physics.

## Source hierarchy

1. Official evaluated fits or machine-readable primary datasets.
2. Digitized primary measurements with an independent digitization receipt.
3. Hash-bound first-principles or atomistic calculations with convergence and ensemble evidence.
4. Labelled surrogate data carrying an explicit model-discrepancy uncertainty.
5. Verification-only synthetic fixtures.

A lower class never silently replaces a higher class. Extrapolation outside the stated temperature,
field, energy, composition, oxygen-content, or irradiation-state domain is rejected.

## Existing numerical fits

`Material_Response_Registry.jl` contains the published NIST cryogenic polynomial forms for:

- OFHC copper specific heat and RRR-dependent thermal conductivity;
- Kapton specific heat and thermal conductivity;
- G-10 CR specific heat and normal/warp thermal conductivity.

The source records are mirrored in `material_sources.toml`. Copper purity and irradiation history
must be retained because RRR strongly changes low-temperature conductivity.

## YBCO and GdBCO tables to create

### Equilibrium heat capacity

Use Phonopy with converged force constants to generate `thermal_properties.yaml`. Convert with
`read_phonopy_thermal_properties_yaml`. Convert molar heat capacity to volumetric heat capacity only
with a separately hash-bound density and molar mass.

### Lattice thermal conductivity

Use phono3py with second- and third-order force constants, or independent Green--Kubo MD replicas.
The minimum table dimensions are temperature and tensor component. Grain texture, oxygen content,
defect state, and strain require separate material-state hashes.

### Electron--phonon and non-equilibrium relaxation

Use EPW or an equivalent converged electron--phonon workflow to produce mode- and temperature-
resolved scattering rates or Eliashberg spectral functions. Use ultrafast YBCO measurements as
state-specific calibration targets, not universal constants.

### Electronic stopping and sub-keV partition

Use real-time TDDFT projectile trajectories or a dielectric-response/track-structure calculation to
partition energy among ionization, electronic excitation, prompt and athermal phonons,
quasiparticles, defect storage, optical emission, escape, and unresolved channels. Export a closed
`energy_eV × temperature_K × channel` table using the template below.

### Cascade damage and stored energy

Use DFT/MLIP-qualified cascade MD over target species, recoil energy, crystal direction,
temperature, oxygen configuration, and independent velocity seeds. Export surviving defects,
antisites, clustering, stored energy, strain, cascade-interface crossings, and prompt recombination.

## Templates

- `subkev_partition_template.csv`
- `green_kubo_template.csv`
- `cascade_md_template.csv`

These are schema examples only. Their values are not physical data.

## Required lineage

Every production table must bind:

```text
material composition and oxygen content
crystal texture/orientation
temperature and magnetic field
strain and irradiation state
DFT/MD code and version
input and pseudopotential hashes
supercell, k/q meshes, cutoffs, and convergence evidence
MLIP hash and cascade qualification status where applicable
independent replica/trajectory identifiers
source and output SHA-256 hashes
units and normalization
uncertainty or ensemble definition
```
