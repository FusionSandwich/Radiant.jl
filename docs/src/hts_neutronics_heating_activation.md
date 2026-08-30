# HTS neutronics, heating, activation, and response control

This page documents the temporary HTS extension currently compiled inside `Radiant.jl`. The code is
intended to move to `RadiantHTS.jl`. The machine-readable reference registry is
`data/hts_physics/references.toml`.

The implementation distinguishes software verification from physical qualification. Analytic and
synthetic fixtures exercise normalization, dimensional, conservation, and ownership contracts;
they cannot set physical OpenMC, OpenSn, evaluated-Gd, activation, curved-geometry, or material-data
gates.

## Transport ownership

The intended production decomposition is:

```text
ParaStell/OpenMC
  global continuous-energy neutron/photon field and boundary-bank reference

OpenSn
  local neutron transport, isotope-resolved reaction rates, neutron-to-photon production,
  and optional external-domain photon transport

Radiant
  photon/electron/positron transport inside selected tape, fibre, or detector microdomains

OpenMC depletion / activation code
  inventories and activity-versus-time fields

SPECTRA-PKA / BCA / atomistic calculation
  nuclear recoil and light-ion handoff
```

Comparison calculations may overlap, but `Layer_Heating_Ledger` permits only one production-additive
owner for one population/domain/layer/response/time key. An OpenSn photon-heating result can be
stored beside a Radiant result for comparison without being added to it.

## Energy and heating ledger

Transport energy deposition is not automatically prompt lattice heat. The protected partition is

```math
E_{dep}=E_{prompt\ heat}+E_{delayed\ heat}+E_{non-eq}+E_{defect}
       +E_{chemical}+E_{escape}+E_{cutoff}+E_{recoil}+E_{unresolved}.
```

Only the prompt/delayed lattice-heat channels and explicitly calculated Joule heat enter the default
production thermal sum. Source-particle energy from activation or Gd capture is transported before
its deposited fraction is added to heat. Neutrino energy, escaped particles, defect storage, and
nuclear-recoil handoffs are never silently counted as heat.

## GdBCO self-shielding and prompt products

`Gd_Self_Shielding_Analytic.jl` provides a groupwise pure-removal slab solution for screening,
source-chain verification, and depth-profile tests. It resolves isotope-specific macroscopic capture
for Gd-155/Gd-157, background absorption, oblique chord length, transmission, subcell capture,
and particle balance.

The analytic solver is deliberately not a replacement for continuous-energy OpenMC or a scattering
OpenSn calculation. A production GdBCO result still requires isotope-resolved evaluated cross
sections, thermal/resonance resolution, explicit film/repeated-tape geometry, and comparison to
continuous-energy transport.

A self-shielded `Capture_Rate_Field` can feed the evaluated-cascade adapter. Prompt gamma rays,
characteristic x rays, conversion electrons, Auger electrons, and residual-nucleus recoil remain
separate source or handoff products. Every complete cascade must close its branch probability and
reaction Q value.

## Response-preserving multigroup condensation

`Response_Preserving_Group_Condensation.jl` creates an exact nested fine-to-coarse group mapping.
For a protected response coefficient `r_g` and reference group-integrated flux `phi_g`, the coarse
coefficient is

```math
\bar r_G=\frac{\sum_{g\in G}\phi_g r_g}{\sum_{g\in G}\phi_g}.
```

It therefore preserves

```math
\sum_g\phi_g r_g=\sum_G\Phi_G\bar r_G
```

for the declared reference spectrum. The result is spectrum-specific; it does not claim response
preservation after the spectrum changes.

Different protected quantities should use independent receipts even when they share boundaries:

```text
REBCO prompt heating
Gd-155/Gd-157 capture
activation reactions
oxygen/Y/Ba/Cu PKA production
neutron-to-photon production
interface photon current
```

The transfer-matrix routine assumes `(incoming group, outgoing group)` ordering, weights incoming
rows by the reference incoming flux, and sums outgoing fine groups. It checks that the aggregated
outgoing source is preserved. Full fine-group covariance can be propagated through the linear
weighting matrix.

## Activation and delayed-source reinjection

`Activation_Delayed_Source.jl` consumes parent-resolved activity in Bq on explicit voxels. It builds
separate delayed photon, electron, and positron volume sources and creates explicit handoffs for
alpha particles, daughter recoils, and neutrinos. Each source carries parent/daughter identity,
decay mode, cooling interval, inventory hash, activation-artifact hash, and decay-data hash.

The decay ledger requires

```math
Q=E_{\gamma}+E_{e^-}+E_{e^+}+E_{\alpha}+E_{recoil}+E_{\nu}+E_{unresolved}.
```

A candidate scheme may retain unresolved energy, but it cannot become production-ready. The
activity-to-source conversion uses `Bq = decays/s`; it does not multiply by an additional physical
source rate.

Boundary current must not be used directly for activation. The activation producer must supply
volume scalar flux or reaction rates, exact material/nuclide inventory, region volumes, irradiation
history, and microscopic cross sections or an equivalent depletion result.

## Material-response bindings

`Material_Data_Bindings.jl` converts hash-bound atomistic response tables into the property objects
used by the cryogenic electrothermal solver. Current supported routes include:

```text
Phonopy -> harmonic heat capacity
phono3py -> thermal-conductivity tensor and phonon lifetime
independent Green-Kubo MD replicas -> conductivity tensor and sampling uncertainty
EPW -> electron-phonon spectral inputs
real-time TDDFT / dielectric response -> electronic stopping and sub-keV partition inputs
cascade MD / DFT -> displacement thresholds, defects, interface crossings, and stored energy
```

Tensor averaging is not implicit. A conductivity component/direction must be selected explicitly.
YBCO data cannot be reused for GdBCO without a declared surrogate model and discrepancy
uncertainty.

## Protected-response uncertainty

`Response_Uncertainty_And_Convergence.jl` records separate standard-uncertainty fields for global
Monte Carlo statistics, source-bank sampling, projection, group condensation, angular and spatial
discretization, faceting/curvature, field maps, evaluated data, material properties, sub-keV and
atomistic models, coupling, electrothermal response, and time integration.

Components assigned to the same declared correlation group are added linearly within that group;
different groups are combined in quadrature. This is conservative for fully correlated components
and avoids pretending that source-generation and source-projection errors are independent when they
are not.

Each protected response can require independent convergence studies over selected axes:

```text
source-bank size
spatial mesh
angular quadrature
energy groups
scattering order
facet refinement
field-map resolution
material/sub-keV response-table resolution
OpenSn-Radiant coupling iteration
time step
```

Passing numerical convergence does not establish physical validity. Physical promotion requires a
matched, hash-bound physical reference under the same source, geometry, material state,
normalization, and score definition.

## Required physical evidence

The remaining external inputs are:

1. The actual patch/energy/direction-resolved OpenMC boundary bank and matching current tallies.
2. Physical neutron, photon, and neutron-to-photon multigroup data with processing provenance.
3. A real OpenSn local result and, for closed decomposition, forward and return angular-current
   artifacts.
4. Complete evaluated Gd-155/Gd-157 de-excitation data including atomic relaxation and recoil.
5. Generated or digitized YBCO/GdBCO thermal, sub-keV, and atomistic response tables.
6. An activation inventory/activity bundle at explicit irradiation and cooling times.
7. At least one genuinely curved neutral-transport reference for the faceted local-tape method.

Geant4 remains useful as an optional event-wise or de-excitation cross-check, but it is not the
current critical path. The first physical milestone is a matched straight-YBCO OpenMC/OpenSn/Radiant
calculation with layer heating, interface current, particle/energy balance, and independent
spatial/angular/energy convergence.
