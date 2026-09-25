# Independent local assessment of the native proton workflow

This report was produced in a separate Codex task and copied into this branch after that task completed. Its raw logs, Julia scripts and TOML measurements remain under `C:\Users\joshu\Documents\Codex\2026-09-25\independently-assess-the-new-proton-workflow\work`.

**Commit:** `e595628759e53d19b98a576dc64ac6fa18a973d6` on `JS/radianthts-full-proton-modeling-20260925` (25 September 2026). Measured runs and tests used a detached worktree of that exact commit. The source checkout was read-only; execution used existing local Julia packages with offline mode enabled. No software or data were acquired. An initial unmanifested package-load probe failed because Julia could not resolve HDF5; the copied local manifest resolved this without installation.

## Finding

The native proton path runs and reproduces an exact synthetic streaming control. The two-group, two-material case also approaches an independently derived constant-stopping characteristic solution when space and energy are refined. Its default one-voxel-per-material grid has a **47.94% error in the small lower-energy flux component of the front layer**, even though total flux and electronic deposition are within 0.033% of the analytic values. This is a coarse-discretization limitation in the tested fixture, not evidence of physical proton accuracy. The code deliberately blocks a physical-mode run with these synthetic tables and blocks nonzero nonelastic removal without a production/event handoff.

## Fixture and independent reference

The baseline uses two 1 cm Cartesian layers, one energy group spanning 1–500 MeV, constant electronic stopping of 2 MeV/cm in both layers, a unit isotropic volume source in the front layer, vacuum ends, zero magnetic field, and zero nuclear stopping, straggling, angular variance, and nonelastic removal. The meaningful variant has two groups, 500–250 and 250–1 MeV, with stopping of 2 and 4 MeV/cm in the front and back layers. Its source occupies the **250–500 MeV band**, uniformly in the source group's first-order representation. It is not a monoenergetic 500 MeV beam. The source basis is `per_source_particle`, with a declared rate of 1 particle/s; the projected source rate is `[0, 1]` in ascending energy order and the native map is `[2, 1]`. Fluxes and scores below are per source particle.

The solver is native SN CSD, two-point Gauss–Legendre quadrature (`μ = ±1/√3`), Legendre order 0, `galerkin-d` angular treatment, order-1 diamond difference in x, and order-1 discontinuous Galerkin in energy. Each layer has one voxel unless a refinement count is shown. The default in-group tolerance is `1e-7`, with a 300-iteration cap. Every benchmark group reported convergence at iteration 2 with `ϵ = 0`; the printed `ρ = NaN` is an undefined acceleration ratio in this zero-scattering fixture, not a finite convergence measure.

For the zero-stopping control, a unit source over the front 1 cm and vacuum boundaries gives a cell-averaged scalar flux of `1/(2|μ|) = √3/2 = 0.8660254037844386` in each layer. The computed values were `0.8660254037844387` and `0.8660254037844389`; deposition was exactly zero.

For the two-group variant, the independent reference integrates the one-dimensional characteristics over source position and a uniform 250 MeV source band. The largest source-to-exit energy loss is `6/|μ| ≈ 10.39 MeV`, so no particle can cross the 1 MeV cutoff. Total group-integrated scalar flux is again `√3/2` in either layer. With `ΔE = 250 MeV`, `S₁ = 2 MeV/cm`, and `S₂ = 4 MeV/cm`, the lower-energy fluxes are `S₁/(6 μ² ΔE) = 0.004` (front) and `(S₁+S₂)/(4 μ² ΔE) = 0.018` (back). The high-group values are total minus low-group flux. Electronic deposition is `Sᵢ√3/2` per 1 cm layer. Independently, constant-stopping CSDA ranges from 500 to 1 MeV are `(500−1)/Sᵢ`: 249.5 and 124.75 cm; the binding returned both exact values.

| Two-group output, per source particle | Analytic | Computed on 1 voxel/layer | Relative error |
|---|---:|---:|---:|
| Front high-group flux | 0.862025404 | 0.860066687 | −0.2272% |
| Front low-group flux | 0.004000000 | 0.005917554 | +47.9389% |
| Back high-group flux | 0.848025404 | 0.848312129 | +0.0338% |
| Back low-group flux | 0.018000000 | 0.017429630 | −3.1687% |
| Front electronic deposition (MeV) | 1.732050808 | 1.731968482 | −0.00475% |
| Back electronic deposition (MeV) | 3.464101615 | 3.462967037 | −0.03275% |

Space and energy refinement used the same 1 cm layer widths and unit source rate. Source density was distributed uniformly over front-layer voxels; high-energy source-group weights summed to one. The table compares the region-averaged low-energy group-integrated flux with the same characteristic reference.

| Voxels per material | Energy groups | Front low-group error | Back low-group error |
|---:|---:|---:|---:|
| 1 | 2 | +47.94% | −3.17% |
| 8 | 2 | +0.068% | −2.45% |
| 16 | 2 | −0.501% | −2.44% |
| 16 | 4 | +0.188% | −0.0877% |
| 16 | 8 | +0.195% | −0.000242% |
| 32 | 8 | +0.0488% | −0.000239% |

Energy refinement alone, with one voxel per layer, left the front low-group error near +50%; space refinement was necessary there. The back residual on a two-group energy grid was reduced by energy refinement. These results establish convergence for this manufactured streaming/stopping case only.

## Local performance and tests

Host: CORSAIR VENGEANCE i7200, Intel i9-10850K (10 physical cores, 20 threads), 32 GB installed RAM. Runs used one Julia 1.10.12 thread on Windows, with HDF5 0.17.2, JLD2 0.5.15, and SpecialFunctions 2.5.1 from a copied local manifest. BenchmarkTools was not installed, so timing used Julia `@timed`. Package loading had already been precompiled by a probe before these fresh-process benchmarks. The cold solve is the first call in each fresh process and includes JIT compilation; three warm solves each used a newly constructed case, with `GC.gc()` outside timing. This is a very small fixture, so sub-millisecond figures should not be extrapolated to production geometry.

| Case | First construction | First solve | Median warm construction | Median warm solve | Warm solve allocation | Process high-water RSS | Whole-process wall |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 group, 2 layers | 1.728 s | 8.641 s | 0.356 ms | 0.765 ms | 127 kB | 1,041 MiB | 51.1 s |
| 2 groups, 2 materials | 1.823 s | 9.201 s | 0.374 ms | 1.201 ms | 169 kB | 1,023 MiB | 53.0 s |

`Sys.maxrss()` is a whole-process high-water value that includes the Julia runtime, package loading, and compilation; it does not isolate solver peak memory. First-call solve allocations were 355 MB and 371 MB, respectively, largely associated with compilation. The repository's targeted `proton_native_tests.jl` passed **46/46** manufactured-binding checks and **18/18** integration checks (exit code 0). No fast SN mode was present in this exact checkout, so none was timed.

## Interpretation and limits

The analytic reference validates one-dimensional, zero-reaction transport and constant stopping for the declared synthetic source. It does not validate measured stopping powers, nonelastic interactions, secondaries, scattering, magnetic-field accuracy, energy straggling, HTS material response, or a real beam. The default one-group case includes source energies down to the 1 MeV cutoff, so its flux is reported as a runtime/software result rather than compared with the two-group analytic reference. The physical-mode guard and nonelastic preflight failures in the targeted tests are appropriate blockers; a synthetic pass must not be promoted to physical qualification. No outgoing-current record was available for an independent event-by-event energy or particle-balance audit.

The copied manifest and scripts reside only in this assessment's detached worktree and output folder. `computer_inventory.json` records RAM, drives, processes, executables, environments, caches, local checkouts, and the zero-byte acquisition plan before the manifest copy and benchmark runs. The failed probe had already created a transient Julia compile cache in `work/julia_depot`; it did not acquire a package. Raw results are in `benchmark_*.toml`, `spatial_refinement.toml`, and `joint_refinement.toml`; exact test and solver output is in the corresponding `*_console.log` files. Reproduction scripts are `proton_benchmark.jl`, `proton_spatial_refinement.jl`, `proton_joint_refinement.jl`, `proton_range_check.jl`, and `run_native_tests.jl`.
