# HTS, diagnostic, and detector physics coverage

This document separates transport physics implemented in Radiant from external production owners
and downstream material-response models. A scalar deposited-energy score must not be treated as a
complete damage, cryogenic-heating, fibre-degradation, or detector-response model.

## Production decomposition

| Process | Production owner | Radiant role |
|---|---|---|
| Photon/electron/positron photoatomic transport above the qualified cutoff | Radiant | Full local transport and process-resolved source/deposition ledger |
| Neutron transport, capture, and Gd self-shielding | OpenSn/OpenMC | Consume boundary photons/electrons and neutron-induced volume sources |
| Photonuclear reactions and photoneutrons | Geant4/MCNP/PHITS | Consume reaction/event kernels; do not infer them from photoatomic attenuation |
| Electronuclear/positronuclear reactions | Geant4 where threshold-relevant | Consume reaction/event kernels |
| Nuclear recoils, p/d/t/alpha/Li ions | SPECTRA-PKA, Geant4, and BCA | Receive weighted recoil/ion energy, direction, position, species, and time |
| Neutron activation | OpenMC depletion/R2S | Transport delayed photons and selected delayed charged-particle sources |
| Photonuclear activation | Qualified photonuclear rate producer plus FISPACT-II/ALARA or verified inventory solver | Transport delayed emissions after a hash-bound inventory handoff |
| Fibre RIA, colour-centre kinetics, FBG shift, recovery | Calibrated optical-response model | Supply dose, spectra, reaction rates, defect/source terms, T, rate, and history |
| HTS detector pulse | Electrothermal/superconductivity response solver | Supply event-wise/layer-wise energy partitions and converter-product transport |

## Energy deposition is not automatically heat

`Energy_Partition` records a closed ledger for

- ionization;
- electronic excitation;
- prompt lattice heat;
- nuclear-recoil energy handed to a recoil model;
- metastable defect-stored energy;
- optical emission;
- escaping-particle energy;
- sub-cutoff handoff;
- unresolved energy.

Only the prompt lattice term is immediately available as heat by definition. Subsequent defect
annealing, quasiparticle recombination, optical absorption, delayed decay, and thermal diffusion are
separate time-dependent sources. Every production result must identify the model used to convert
transport deposition to heat.

## HTS detector use case

The current detector schema supports an explicit converter layer, active YBCO/GdBCO region,
substrate, operating temperature, bias, and a separate thermal/electrical model. The included
B-10-enriched B4C/YBCO example is a verification fixture only. Converter reaction rates and
alpha/Li-7 escape require a nuclear reaction/charged-particle owner. Radiant handles the coupled
electromagnetic deposition after those particles and photons enter its qualified domain.

Protected detector responses include active-layer deposition distributions, prompt heat,
ionization/excitation, converter-product escape, pulse height, decay time, and resistance/current
response. Pulse predictions require heat capacities, active/substrate and substrate/bath thermal
conductances, transition curve, magnetic field, bias, and readout parameters.

## Fibre diagnostic use case

`Fibre_Diagnostic_Definition` represents concentric core, cladding, coating, adhesive, and capillary
layers. Protected radiation responses include core/cladding dose, neutron capture, displacement
source, secondary electron spectrum, H/He production, and optical source terms.

Cherenkov, scintillation, and radioluminescence generation require optical material data and
interface transport. RIA and FBG wavelength drift additionally require dose/fluence, dose rate,
temperature, spectrum, impurity/dopant state, stress, and annealing history. These are not inferred
from dose alone.

## Physics additions that are intentionally not duplicated

Radiant should not duplicate the full neutron, photonuclear, activation, BCA, molecular-dynamics,
optical-photon, or electrothermal solvers. It should provide conservative interfaces, explicit
ownership, process-resolved responses, and closure ledgers. Neutrinos and cosmic muons are not added
merely for completeness; they become explicit only for a defined detector/background case with a
non-negligible response.
