# Post-ASC proton implementation handoff

## Executable result

The isolated branch `JS/radianthts-full-proton-modeling-20260925` starts from verified Radiant `4f7f83937fcd23ffeff195bf75ceaa2c2e22c6d6`. `test/proton_native_tests.jl` constructs `Proton()` and an explicitly normalized `Anisotropic_Volume_Source`, binds two manufactured material tables through `Cross_Sections.build`, builds two native Cartesian layers, selects `SN` CSD or BFP, and runs `Computation_Unit.run` with zero and nonzero magnetic fields. The BFP case exercises a manufactured angular-variance moment. It scores electronic deposition and nuclear-stopping recoil energy separately; the native process sum agrees with the native total deposition score. The physical classification remains **software verification only**.

The generic core extension is the `provided` multigroup builder plus a transport preflight hook in `Cross_Sections`. `Proton_Native_Binding.jl` populates native group and response structures from exact-domain tables. `transport.jl` invokes the guard on every native solve, including a direct `Computation_Unit.run`, so absent or nonzero unconnected nonelastic physics cannot be silently zeroed. The guard also compares table snapshots to detect mutation after binding. Candidate/qualified tables require a matching local SHA-256 artifact. `Geometry.build` now records a successful build so explicit source projection can run. The magnetic operator receives the transported particle rest mass instead of always using the electron value.

## Exact breaks found and fixed

1. Existing `ion_transport_step` was disconnected from `Cross_Sections`, `SN` and `Computation_Unit`. Generic provided-library and preflight hooks now connect the manufactured table to native CSD.
2. Legacy `Volume_Source` produces a tensor incompatible with CSD `Source` (`Source.jl:add_source` shape check). The explicit `Anisotropic_Volume_Source` projects into the CSD representation and preserves a declared `Source_Normalization`.
3. `Geometry.build` did not set `is_build`, while explicit source projection requires it. The core method now sets it after successful construction.
4. `electromagnetic_scattering_matrix` used the electron rest energy for all charged particles. `sn_flux.jl` now passes `get_mass(part)`.

## Remaining critical path

1. Bind exact product/lot layer and material-state packets, physical source phase space/rate, and licensed, hashed 1–500 MeV stopping, straggling and screened Coulomb/nuclear elastic data. Keep all uncovered intervals blocked. DPA reports zero qualified material-energy cells and zero stopping/scattering physical benchmark rows.
2. Qualify a 200–500 MeV pion-capable nuclear event producer against independent isotope/energy/differential observables, with residual and correlated product closure. The IAEA selected database stops at 100 MeV; TENDL-2025 main evaluations stop at 200 MeV. Geant4 is only a producer candidate. No external transport/model code was run here.
3. Add a generic coupled production-matrix interface and a correlated weighted-event handoff to the native library; reconcile matrix means against event-family rates. Define one additive owner for every population/domain/response/time tuple and block unsupported products before transport.
4. Retain outgoing SN surface angular currents and cutoff energy flow in the native `Flux` result, then close source/current/energy/charge/baryon/momentum ledgers without inferred leakage. The current two-layer test scores process and layer deposition, but full balance cannot pass with the present `Flux` shape.
5. Couple mapped/curved metrics and cell-centred fields to the same native collision/CSD sweep. Existing mapped streaming is manufactured and separate; `sn_flux.jl` currently accepts Cartesian transport only.
6. Bind residual yields to activation and directional recoils to PKA/Beyond-DPA, with source and daughter energy owned once. Perform exact-material physical benchmarks, covariance/model-discrepancy and convergence studies on both Julia versions before any production claim.

## Validation receipts

- Julia 1.10.12, offline existing depot and local 11,601-byte manifest: `C:\rhts-sandbox\julia\julia-1.10.12\bin\julia.exe --startup-file=no --compiled-modules=no --project=. test\runtests.jl` **339/339 assertions passed** (293 inherited/HTS start, 46 new proton native tests).
- Julia 1.6.7, offline existing depot and local 9,678-byte manifest: `C:\rhts-sandbox\julia\julia-1.6.7\bin\julia.exe --startup-file=no --compiled-modules=no --project=. test\runtests.jl` **339/339 assertions passed**.
- JSON gate and coverage files parse with Python; `git diff --check` has no whitespace errors. Physical reference gates remain unrun and failed closed.

No runtime was downloaded. `Manifest.toml` in this isolated worktree was ignored and copied from local versioned manifests, one Julia version at a time; it was removed after validation. No SSH, job, PR, merge, push or ASC/poster write was performed. The authoritative DPA, Radiant and Para worktrees were used read-only and verified clean at their recorded HEADs. No ASC poster path was accessed.
