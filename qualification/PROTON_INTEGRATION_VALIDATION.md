# Proton repository integration validation (2026-09-25)

## Published scope and implementation plan

The [Radiant upstream README](https://github.com/CBienvenue/Radiant.jl) and [Getting Started guide](https://cbienvenue.github.io/Radiant.jl/user_guide_getting_started/) describe a solver specialized for coupled photons, electrons and positrons, with heavy charged particle building blocks. The [Particles guide](https://cbienvenue.github.io/Radiant.jl/user_guide_particles/) documents `Proton()` as a reserved particle with mass and charge. The public guide's complete calculation uses `Cross_Sections`, `Geometry`, `Solvers`, `Fixed_Sources` and `Computation_Unit.run`; it does not present a qualified proton calculation. The public upstream GitHub issue list and the FusionSandwich fork's issue list showed no proton roadmap at inspection. This is a finding from the pages inspected, not a claim about private plans.

The concrete 1–500 MeV plan for this branch is local: `docs/src/proton_transport_requirements.md`, `docs/src/proton_transport_capability_audit.md`, `PROTON_TRANSPORT_GATE_MATRIX.json` and `PROTON_DATA_COVERAGE_LEDGER.json`. The staged sequence is native source and stopping integration, exact material and source provenance, nonelastic production/event coupling, owned secondary handoffs and conserved outgoing currents, mapped geometry and field coupling, then independent physical benchmarks and uncertainty. The present slice is explicitly synthetic and cannot satisfy the physical gate.

## Integration exercised

`test/proton_native_tests.jl` now includes 18 additional assertions through the repository's public `Proton()`, `Cross_Sections.build`, `Geometry.build`, `Fixed_Sources`, `SN`, `Computation_Unit.run`, `get_flux` and `get_energy_deposition` APIs. A source in the upper of two energy groups maps to Radiant's descending group order, crosses two different material layers, and produces finite group-by-voxel flux and layer-dependent scores. Tripling source strength triples flux and deposition; changing the declared physical source rate scales only the physical-rate score. The process-resolved electronic channel agrees with native deposition.

The pass exposed a guard gap by code inspection: replacing an entry in a bound data vector or mutating a built multigroup array could bypass the original per-table snapshots. The adapter now hashes the bound input structure and generated native library and checks both before transport. Regression cases reject dropped nonelastic declarations, altered native total cross sections, changed group boundaries, changed material density and in-place table changes. These are programmatic integrity checks within Radiant's mutable object model; deliberately replacing the public preflight callback is outside this guarantee.

The first draft of one assertion used the internal six-dimensional flux shape for the public `get_flux(cu, proton)` result. The public method returns a two-dimensional group-by-voxel matrix in this 1D case. Correcting the assertion resolved that test error; the solver result itself was finite.

## Reproducible results and environment

- Julia 1.10.12: full `test/runtests.jl` **357/357 assertions**, exit 0.
- Julia 1.6.7: full `test/runtests.jl` **357/357 assertions**, exit 0.
- Count: 293 prior repository assertions, 46 original proton native assertions, 18 new integration assertions. Runs set `JULIA_PKG_OFFLINE=true`, `JULIA_PKG_PRECOMPILE_AUTO=0`, the version-matched local `JULIA_DEPOT_PATH`, and `--startup-file=no --compiled-modules=no --project=.`.
- Before reuse of local manifests, current RAM was 31.91 GiB total and 3.06–4.08 GiB free; C: free 79.71–80.47 GiB and D: free 470.61–470.68 GiB. Largest processes included WSL VM, ChatGPT, Chrome and Codex. Existing Julia 1.10.12/1.6.7, Python 3.12.10, Git 2.55.0, Julia depots, pip cache and source worktrees were observed. Python prefix equalled base prefix; no Conda executable was on PATH. Existing manifests were 11,601 and 9,678 bytes. Planned external acquisition was 0 bytes; only one manifest at a time was copied into this isolated worktree, and the ignored copy was removed afterward.

## Limits

Passing synthetic native transport and API integration does not qualify physical 1–500 MeV proton transport. Exact HTS stopping/scattering data, nonzero nonelastic event families, outgoing current/energy closure, mapped/curved native solves and independent reference comparisons remain blocked in the gate matrix. For zero-scattering manufactured groups, the solver prints an undefined convergence-rate diagnostic (`ρ = NaN`) alongside zero residual and finite converged flux; this did not cause a test failure and is not a physical-validation result.
