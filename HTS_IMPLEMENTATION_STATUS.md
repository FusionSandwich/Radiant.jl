# Radiant HTS implementation status

Branch: `codex/radianthts-core-foundation-20260828`

Base: `e95eab60716a1a4b010b3c8669705809dc9065cb` (`v1.1.67`)

This branch is an isolated implementation checkpoint. It has not been merged into `main`, and no
pull request has been opened.

## Implemented in this checkpoint

### Source and normalization contracts

- Explicit per-history, per-source-particle, or per-second normalization.
- Separate physical source rate and symmetry factor.
- Prompt/delayed time class and time interval.
- Source hash and provenance.
- Patch/group/direction boundary angular-current source.
- Isotropic, ordinate, and angular-moment volume sources.
- Optional source variances.
- Explicit rejection of negative source densities and unresolved nonzero grazing-current bins.

### Executable Radiant adapters

- Exact source-to-Radiant energy-group interval mapping.
- Exact source-to-SN ordinate and quadrature-weight mapping.
- Cartesian boundary patch to boundary-cell-face mapping.
- Half-range boundary-moment projection.
- Volume-moment projection.
- Current/rate closure receipts.
- Negative reconstructed angular-source rejection without clipping.
- Boundary values installed as angular-flux densities; face measure is used only in current
  integration.

### Fixed-source lifecycle

- Explicit and legacy normalization contracts cannot be mixed.
- Compatible explicit sources may be combined.
- Source building is idempotent.
- Rebuilding after source addition starts from empty derived state.
- Projection receipts remain attached to the fixed-source collection.

### HTS structures

- Layer definitions with transport subdivision, activation, and atomistic-output flags.
- Ordered tape-stack definition.
- Canonical 4 mm wide, 163.2 micrometre, eight-layer verification fixture.
- YBCO default and GdBCO material-tag schema test; no Gd self-shielding qualification is claimed.

### Accounting and ownership

- Particle, kinetic-energy, rest-mass, leakage, cutoff, deposition, and charge ledger type.
- `parastell_damage.transport_ownership_map/v1`.
- Production uniqueness and non-additive comparison modes.

### Tests authored

- Source normalization and validation.
- Boundary current integration and conservative conversion.
- Grazing and negative-source rejection.
- Volume-source scalar and moment conventions.
- One-voxel boundary and volume projection.
- Projection receipts.
- Fixed-source idempotency and mixed-contract rejection.
- Eight-layer fixture.
- Ownership uniqueness.
- Particle/energy/charge ledger residuals.

## Qualification status

| Gate | Status |
|---|---|
| Branch isolated from `main` | PASS |
| No PR/default-branch merge | PASS |
| Source schema implemented | PASS (static) |
| Boundary-current projection implemented | PASS (static) |
| Volume-source projection implemented | PASS (static) |
| No-double-counting map implemented | PASS (static) |
| Regression tests authored | PASS |
| Julia package compilation | NOT YET EXECUTED |
| Julia 1.6 test run | NOT YET EXECUTED |
| Julia 1.10 test run | NOT YET EXECUTED |
| Physical OpenMC/OpenSn source replay | NOT YET EXECUTED |
| RADIANT_BUILD_PASS | NOT YET ESTABLISHED |
| RADIANT_EM_PASS | NOT YET ESTABLISHED |
| OPENSN_RADIANT_COUPLING_PASS | NOT YET ESTABLISHED |

The code was written through the GitHub connector. Connector-authored commits did not start a
GitHub Actions run, and the execution environment used for this checkpoint did not contain a Julia
runtime. Therefore, no runtime gate is claimed from this checkpoint.

## Next implementation slice

1. Run `Pkg.test()` on pinned Julia 1.6 and 1.10 and correct all compile/runtime failures.
2. Add standard cross-language HDF5 boundary and volume source interchange.
3. Preserve outgoing boundary angular moments/current in transport results.
4. Compute native particle, energy, charge, leakage, and cutoff balances.
5. Add layer-, process-, and reaction-resolved responses.
6. Add OpenMC/OpenSn adapters and replay the analytic boundary fixtures.
7. Add the piecewise-flat tape atlas and curvature convergence metrics.
8. Add activation/delayed-source and PKA/atomistic event-kernel exports.
9. Qualify GdBCO self-shielding against continuous-energy reference calculations.
