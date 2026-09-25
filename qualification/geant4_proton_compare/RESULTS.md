# Matched Geant4/Radiant proton slab comparison (2026-09-25)

The [low-energy extension](LOW_ENERGY_RESULTS.md) adds 1.2, 2 and 5 MeV matched Geant4/Radiant cases and an independent NIST PSTAR electronic-CSDA reference. Its claim remains limited to electronic stopping and thin Al/Cu slab transport.

## What was compared

This is a controlled **electronic-stopping transport** comparison at 10, 100 and 499.99 MeV. Geant4 11.4.2 transports one proton per history through a 0.01 cm aluminum layer followed by a 0.01 cm copper layer. The source is uniform in the Al layer and in a 0.01 MeV-wide energy interval centred on the stated energy. Its two possible x-direction cosines are `±1/sqrt(3)` with equal probability. The Geant4 model registers proton `G4hIonisation`, transportation and a step limiter, sets a 1 m secondary production cut, disables energy-loss fluctuations, and has no multiple scattering or hadronic reactions. The final production run uses a 0.001 mm maximum step, one million histories per energy, and seed `731293`. Geant4 material densities, carried verbatim to Radiant, are Al 2.699 and Cu 8.96 g/cm³. Geant4's [application guide](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/TrackingAndPhysics/physicsProcess.html) documents the ionisation/stopping and production-cut interfaces; [its cut reference](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/energy_loss/setcuts.html) distinguishes production cuts from tracking cutoffs.

The Geant4 proton ionisation process exports a 301-point, 1–500 MeV electronic stopping table in MeV/cm. Its SHA-256 is `929e06363ccd6f702b95c09cb1354f096dfc150369fadd7ea0ec648929f0cfa9`. Radiant imports that same table as candidate material data and uses native 1D SN/CSD with S2 Gauss–Legendre directions, DD in space, DG1 in energy, no scattering or nonelastic removal, and the identical layer and source definitions. The source rate is checked as one proton per source. Radiant's `get_energy_deposition` result is mass-normalized per voxel (`MeV/g × cm` in 1D per source), so each score was multiplied by its material density and voxel width before comparing with Geant4 deposited MeV. The binding's `PHYSICAL_CANDIDATE` receipt describes table provenance, **not** a qualified full proton model; the nonelastic declaration remains synthetic zero.

## Independent deterministic reference

[`continuous_reference.py`](continuous_reference.py) integrates `dE/ds = -S(E)` along the two S2 characteristics using the exported Geant4 table. It averages source position and energy with Simpson quadrature and uses RK4 path steps no longer than 0.00001 cm. Doubling source-position quadrature from 100 to 200 intervals and halving the path step changed the 10 MeV Al/Cu means by less than `7×10⁻⁹` MeV; the earlier values are preserved in `continuous_reference_initial.csv`. This reference shares stopping physics with both runs but does not share their transport discretization. The minimum predicted exit energy at 10 MeV is 1.021 MeV, close to Radiant's 1 MeV lower boundary; the 100 and 499.99 MeV exits remain well above their lower energy boundaries.

## Results

Raw Geant4 deposited-energy means and standard errors are in [`geant4_scores.csv`](geant4_scores.csv). The Radiant grid study is in [`radiant_scores_mass_corrected.csv`](radiant_scores_mass_corrected.csv); [`deposition_comparison.csv`](deposition_comparison.csv) gives all pairwise differences. Values below are MeV deposited in each material per source proton. Errors shown for Geant4 are Monte Carlo standard errors only.

| Incident energy | Material | Geant4, 1 µm max step | Characteristic reference | Radiant grid | Radiant score |
|---:|---|---:|---:|---:|---:|
| 10 MeV | Al | 0.822799 ± 0.000486 | 0.823082 | 64 voxels × 800 lower groups | 0.823192 |
| 10 MeV | Cu | 2.987165 ± 0.003014 | 3.000524 | 64 voxels × 800 lower groups | 2.997493 |
| 100 MeV | Al | 0.132768 ± 0.000077 | 0.132773 | 32 voxels × 400 lower groups | 0.132774 |
| 100 MeV | Cu | 0.377041 ± 0.000376 | 0.376502 | 32 voxels × 400 lower groups | 0.376506 |
| 499.99 MeV | Al | 0.050696 ± 0.000029 | 0.050680 | 32 voxels × 400 lower groups | 0.050680 |
| 499.99 MeV | Cu | 0.145948 ± 0.000146 | 0.145952 | 32 voxels × 400 lower groups | 0.145952 |

At 100 and 499.99 MeV, the selected Radiant scores are within 0.15% of Geant4 and nearly equal to the independent characteristic calculation with shared stopping. At 10 MeV, both codes show numerical sensitivity: Radiant's 8-voxel/100-group copper score is 2.9451 MeV, while 64 voxels/400 groups gives 2.9989 MeV, 32 voxels/800 groups gives 2.9915 MeV, and the joint 64-voxel/800-group result is 2.9975 MeV. The latter is 0.10% below the characteristic reference. Geant4's 10 MeV Al/Cu scores change from 0.82148/2.97573 MeV at 5 µm steps to 0.82280/2.98716 MeV at 1 µm steps, moving toward the characteristic values 0.82308/3.00052 MeV. The joint-grid Radiant result is 0.05% above Geant4 in Al and 0.35% above it in Cu; the Cu difference is 3.43 Geant4 standard errors. Geant4 step convergence at 10 MeV is incomplete, so this residual cannot be assigned to Radiant physics alone.

## Independent stopping-model screen

[`compare_bethe.jl`](compare_bethe.jl) evaluates Radiant's inherited analytic `bethe` expression directly; this expression is **not** the stopping input used in the transport comparison above. Against Geant4's `hIoni` table for these elemental solids, Radiant differs by +0.32%/+0.74% in Al/Cu at 10 MeV, +0.11%/+0.53% at 100 MeV, and +0.16%/+0.18% at 499.99 MeV. At 1 MeV it is 0.95%/2.26% lower. The raw values are in [`bethe_vs_geant4_stopping.csv`](bethe_vs_geant4_stopping.csv). Agreement between two models is a screening result; it is not an experimental stopping benchmark or a test of HTS compounds.

## Limits and reproducibility

This case probes Radiant's native source, multigroup stopping, layer score, and unit conversion against a matched Geant4 ionisation-only model. Because Radiant imports Geant4's stopping table, the deposition comparison **cannot independently validate electronic stopping physics**. Neither run represents complete 1–500 MeV proton transport: nonelastic reactions, secondary neutrons/light ions/photons, nuclear recoil, multiple scattering, energy straggling, outgoing-current closure, magnetic fields and exact HTS materials are outside this matched case. Radiant's physical-mode guard still blocks a full model. The next physical accuracy gate needs exact-material stopping and nuclear data, a qualified secondary/event handoff, and measurements or independent full-physics reference cases.

The source, CMake build and run commands are in this directory. [`LOCAL_BUILD_INVENTORY.md`](LOCAL_BUILD_INVENTORY.md) records the local-only environment and acquisition gate. `build_local.sh` and `run_local.sh` use the existing WSL Geant4 installation; `run_radiant.jl` uses the existing Julia 1.10.12 depot/manifest. No software or data were downloaded. `geant4_scores_1m_0p005mm.csv` preserves the 5 µm comparison. The raw Geant4 and Radiant logs preserve the process settings and convergence lines; `ρ = NaN` in Radiant's zero-scattering fixture is an undefined rate diagnostic alongside finite converged scores.
