# Proton transport requirements and architecture

## Boundary and native ownership

Proton support extends the existing multigroup Boltzmann Fokker-Planck method:
the shared SN sweeps, space/energy closure relations, angular operators and
cross-section/response interfaces remain the production framework. Proton
interaction models supply coefficients and response data to that framework.
Existing DD, DG and admitted AWD closures should be assessed before proposing
a new positivity treatment. Any generic numerical change also needs regression
evidence for the existing charged particles. See [Bienvenue and Hebert (2022)](https://doi.org/10.1016/j.anucene.2022.109032)
and [Bienvenue et al. (2025)](https://publications.polymtl.ca/61927/) for the
method lineage; their lepton validation does not validate proton physics.

The contract is 1–500 MeV **kinetic energy per proton**. Every source group, material-state table and production kernel must cover its requested interval without extrapolation. `Proton()` uses native Radiant particle identity. `Cross_Sections.build` owns binding to native multigroup coefficients; `SN` with CSD/BFP/FP owns the charged transport sweep, electronic stopping, optional nuclear-stopping energy transfer, field coupling, and Cartesian layer-resolved scores. `Proton_Native_Binding.jl` is the HTS adapter and preflight, not a second stepper. The standalone `ion_transport_step` remains a software fixture.

Nonzero angular variance requires the native BFP or FP angular Fokker-Planck
operator. BFP-EF does not apply that operator and accepts only zero angular
variance in this binding. Nonzero energy straggling rejects because a native
energy-straggling operator is not connected. A schema field or a standalone
step utility does not establish that the native transport applies its physics.
Angular-operator monotonicity also does not guarantee nonnegative space/energy
flux reconstruction.

The present executable slice uses an explicit synthetic zero nonelastic declaration and a native CSD solve. It does **not** claim a physical proton solve. The native adapter stores electronic deposition separately from nuclear-stopping recoil energy. Nuclear recoil energy is a handoff channel, not instant heat. The current binding rejects a nonzero nonelastic removal coefficient because no coupled production matrix and correlated event handoff are yet connected. A future producer must supply both before a solve can continue.

## Required model bindings

1. A source bank declares species/charge/mass, geometry and material hashes, spatial/angular/energy/time distribution, source probability/weight, physical rate and symmetry factor. `Source_Normalization` converts from history or source-particle basis once. The source group containing 500 MeV must terminate exactly at 500 MeV. Point energy at an edge requires a declared line-source convention; a wide multigroup bin is not a monoenergetic beam.
2. For every exact material and phase, stopping/scattering inputs declare density and temperature or other material state, units, interpolation law, energy bounds, source artifact SHA-256, uncertainty/covariance status, and evidence class. Candidate and qualified bindings verify the local artifact hash. Per-material provenance survives group generation and scoring. Electronic stopping, nuclear stopping, straggling and angular laws remain separate channels.
3. A nonelastic kernel declares target isotope, energy domain, removal rate, outgoing-species multigroup production matrix, correlated event-family law, residual inventory and all closure terms. A deterministic matrix supplies mean transport. A weighted event bank is a quadrature or importance-sampling representation of those rates: sample each family with declared proposal probability `q_f`, assign weight proportional to its deterministic rate divided by `q_f`, record covariance and normalization, and never label samples analog histories. Matrix and bank must reconcile by family, product and energy interval without adding both contributions twice.
4. Each emitted population has one production and one additive transport/response owner under `(population, domain, response, time class, source class)`. `Transport_Ownership_Map` enforces unique production records; comparison solvers are non-additive. Neutrons go to an owned neutral solver; photons and e± to the chosen Radiant/OpenSn domain; protons and light ions to native charged transport or an explicit BCA handoff; heavy recoils to directional PKA/BCA/MD; residuals to activation; significant pion/muon channels to qualified transport or an explicit decay/escape handoff; neutrino energy to escape. No OpenSn/Radiant heating sum is accepted without this map.
5. A native sweep must couple field data sampled on the same geometry, with source, material, field and geometry hashes. The present `Electromagnetic_Field` supports uniform and cell-centred magnetic input; the rest-mass parameter now follows `Particle`. `sn_flux.jl` still requires Cartesian geometry. Mapped/curved streaming has a manufactured metric kernel but no native collision/CSD solve: that integration is a separate runtime gate.
6. Every layer score must preserve process, particle, direction or explicitly state its angular integration, time class, and owner. Ledger closure includes source current and rate, particle current, kinetic energy, rest-mass exchange, charge, baryon number, available momentum, electronic deposition, recoil/cutoff handoff, leakage, activation precursor mass and delayed emission. Current `Flux` stores cell/energy moments and cutoff flux but not outgoing SN surface currents; therefore a full current/energy balance receipt cannot yet pass. Do not infer leakage from a residual and call it measured.

## Data and producer routes

| Route | Directly evidenced scope | Use and open issue |
|---|---|---|
| Evaluated data | [IAEA charged-particle database](https://nucleus.iaea.org/Pages/charged-particle-cross-section-database.aspx), static 2011 catalogue, covers selected medical reactions to at most 100 MeV. [TENDL-2025](https://tendl.ic.ac.uk/tendl_2025/tendl2025.html) lists evaluated particle libraries to 200 MeV and separate *test files* to 600 MeV. | Candidate below its verified isotope/channel/energy range only; no 200 MeV evaluation is extended to 500 MeV. Test files are not silently promoted. Processing, differential products, covariance and reuse terms require a per-artifact review. |
| Open-source event generator | [Geant4 Bertini reference](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/hadronic/BertiniCascade/Cascade.html) describes proton/pion induced reactions within its stated model range; [Geant4 license](https://geant4.web.cern.ch/license/license) permits source redistribution subject to notice and naming conditions. | A 200–500 MeV candidate producer, not a chosen or qualified physical kernel. Pin Geant4 release, physics list, datasets, model handoffs, seeds and event/residual output; compare alternatives and closure against independent data. No Geant4 execution is authorized here. |
| Independent references | [NIST PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) provides proton stopping/range calculations through 10 GeV for its stated materials; [IAEA EXFOR](https://nucleus.iaea.org/Pages/experimental-nuclear-reaction-data.aspx) compiles experimental nuclear data. | PSTAR can benchmark matching stopping materials and states, not arbitrary HTS phases or nonelastic event yields. Select isotope/energy/observable-matched EXFOR experiments and independent transport/model results for secondary spectra, residuals and pion channels. These are comparisons, not inputs that make a producer qualified automatically. |

For each imported data family, record official URL or DOI, version/date, license or reuse terms, exact isotope and material state, interval, interpolation, uncertainty, source file SHA-256, processing recipe, and whether a statement is direct evidence or inference. The present gate ledger contains no imported physical rows.

## Qualification sequence

Formula checks cover relativistic values, units, range integration and no extrapolation. Manufactured solver checks cover source normalization, zero and nonzero fields, two Cartesian layers, process-score closure, nonelastic fail-closed behavior and mapped-coordinate invariance. Cross-language tests must validate the DPA `ion_source_bank`, `spallation_event_family_ledger`, `secondary_routing_ledger`, activation and directional PKA/Beyond-DPA schemas with canonical hashes. Physical-reference tests need exact material, beam and model provenance with independent stopping/range, lateral-spread, nonelastic, secondaries/residuals and leakage comparisons on both Julia 1.6.7 and 1.10.12. No software fixture promotes a physical gate.

Future downstream “pollination modeling” is ambiguous and out of scope. A later application could consume a qualified source/fluence field, deposited-energy and activation ledgers, and uncertainty/ownership receipts without changing this proton physics contract.
