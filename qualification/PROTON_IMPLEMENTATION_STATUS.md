# Post-ASC proton implementation status

## M0 — isolated start (2026-09-25)

- Authoritative Radiant local branch, HEAD, and remote branch all equal `4f7f83937fcd23ffeff195bf75ceaa2c2e22c6d6`; `origin/main` equals `e95eab60716a1a4b010b3c8669705809dc9065cb`. The two histories have no merge base, so comparisons use `git diff origin/main HEAD`.
- First write was `git worktree add -b JS/radianthts-full-proton-modeling-20260925 C:\rhts-proton\radiant 4f7f83937fcd23ffeff195bf75ceaa2c2e22c6d6`.
- DPA and Para read-only HEADs match `b0cc63506a82c721035019979f3a47ca50baaf15` and `913532e4f99bdcc7687e7ec259feab9c0e51cec0`.
- No `AGENTS.md` was found inside the new worktree or its checked parent directories. User-wide rules supplied in the task apply.
- Existing authoritative worktrees and ASC poster locations are excluded from all writes. No SSH, job, push, PR, software acquisition, or physical-data acquisition is authorized.

## Local resource and software inventory before any possible acquisition

- Windows visible RAM 31.91 GiB, free RAM 2.60 GiB at check. C: (system and target) free 81.00 GiB; D: (data) free 470.98 GiB. `Temp` maps to D:.
- Largest working sets: WSL VM 2.06 GiB, Chrome 1.30/0.88/0.62/0.53 GiB, ChatGPT 1.28/0.52 GiB, Codex 0.65 GiB, LEDKeeper2 0.58 GiB, Explorer 0.54 GiB, OneDrive 0.51 GiB, Defender 0.48 GiB.
- Existing executables: portable Julia 1.10.12 and 1.6.7 paths under `C:\rhts-sandbox\julia`; Python 3.12 at `C:\Users\joshu\AppData\Local\Programs\Python\Python312\python.exe`; Git 2.55.0 at `C:\Program Files\Git\cmd\git.exe`. Python's prefix and base prefix both resolve to that Python312 installation. No `conda` command on PATH and no standard venv/Conda directory was found.
- Existing Julia depots: `C:\rhts-sandbox\julia-depots\1.10.12`, `C:\rhts-sandbox\julia-depots\1.6.7`, plus `C:\Users\joshu\.julia`; versioned manifests under `C:\rhts-sandbox\manifests`. Existing pip cache under `C:\Users\joshu\AppData\Local\pip\cache`. Local source checkouts include `C:\rhts-ta\radiant`, `C:\rhts-ta\dpa`, `C:\rhts-ta\para` and other registered Radiant worktrees.
- Planned external acquisition: **0 bytes**. Existing local runtimes and sources are sufficient for software implementation and tests. Missing physical nuclear data will remain explicit blockers.
- Local dependency-environment reuse required for tests: the isolated worktree initially lacks ignored `Manifest.toml`, and Julia cannot resolve HDF5 from the project without one. Two already-local manifests exist, 11,601 bytes for Julia 1.10.12 and 9,678 bytes for 1.6.7; the latter matches the existing authoritative checkout's manifest hash. The `Project.toml` files in those checkouts have identical content (Git reports only working-tree line-ending conversion). Plan: copy one manifest at a time into **this isolated worktree only**, set `JULIA_PKG_OFFLINE=true`, disable compiled modules and automatic precompile, run tests, then replace it with the other local manifest. Maximum copied bytes in the target at once: 11,601. Rollback: remove this worktree's ignored `Manifest.toml`. No package download, installation, or update is planned.

## Milestones

- M1 audit and requirements: complete in `docs/src/proton_transport_capability_audit.md`, `docs/src/proton_transport_requirements.md`, and JSON gate/data ledgers. Inherited upstream mechanisms, HTS additions and blockers are separated.
- M2 public API trace: complete. First break was missing binding from `Tabulated_Ion_Transport_Model` to `Cross_Sections`; legacy `Volume_Source` tensor also failed CSD projection. Exact evidence is in `PROTON_IMPLEMENTATION_HANDOFF.md`.
- M3 native manufactured slice: complete. Two-layer Proton `Cross_Sections.build` → explicit `Fixed_Sources` → native `SN` CSD/BFP → `Computation_Unit.run` → process/layer scores, with zero/nonzero field cases.
- M4 provenance and fail-closed ownership: complete for the manufactured binding. Tables declare domain, units, material state/density, interpolation, source hash, uncertainty and status. Candidate/qualified source artifacts are hash checked. Missing nonelastic input, nonzero removal without production, unsupported angular solver and physical-mode synthetic input block before transport.
- M5 tests: complete and extended by the repository integration pass. Julia 1.10.12 and 1.6.7 each passed the full `test/runtests.jl` suite: 357/357 assertions (293 existing, 46 native manufactured proton, and 18 new two-group integration assertions). Runs used local versioned manifests, offline depots and `--compiled-modules=no`.
- M5 follow-up: two native groups and two differing material layers, source-strength linearity, source-rate normalization, public flux/score agreement and post-build mutation guards pass. Findings and the public upstream documentation comparison are in `PROTON_INTEGRATION_VALIDATION.md`.
- M6 documentation and handoff: complete, with physical gates explicitly blocked. Final read-only non-interference check found the authoritative Radiant, DPA and Para worktrees clean at their recorded HEADs; no ASC poster path was accessed or written.
- M7 scoped local commit: complete on the isolated branch; the commit ID is this branch's HEAD. No push or PR was created.
