# CSD energy moment scoring correction

The native SN/GN CSD operator represents stopping within each energy group by
the affine polynomial through its boundary stopping values. With the basis
`sqrt(2n+1) P_n`, its nonzero coefficients are

`S0 = (S_hi + S_lo)/2`, `S1 = (S_hi - S_lo)/(2 sqrt(3))`.

The previous deposition score used only `Sigma_e * Phi0`. Its cross-section
generator embeds a separately quadrature-averaged soft stopping `Saverage` in
`Sigma_e`. For native libraries with that decomposition available, scoring now
uses

`(Sigma_e - Saverage) Phi0 + S0 Phi0 + S1 Phi1`,

with the last term present only when the solver carries the first energy moment.
The pure energy moment has index 2 in full and reduced SN/GN spatial/energy
layouts. Nonconstant spatial and angular moments average to zero in this score.
Higher energy moments do not contribute to the represented affine soft stopping.
The existing cutoff score is retained separately.

FMAC-M and MATXS imports can provide authoritative deposition and boundary
stopping without a separately stored group-average soft term. Those materials
retain their existing deposition score. A boundary interpolant cannot identify
the missing decomposition. Higher-order scoring of such imports is unqualified.
The non-CSD scoring branch is unchanged.

## Independent numerical regression

`test/csd_energy_scoring.jl` supplies mathematical coefficients to native score
objects; it does not replace native cross-section generation as a physics gate.
For `S(x)=2+x` and `psi(x)=1+x`, direct integration on `[-1,1]` gives
`mean(S psi)=7/3`. The test deliberately stores `Saverage=7`, distinct from
`S0=2`, to expose both an unreplaced soft term and an omitted first moment.
Large nonconstant spatial moments detect contamination of the voxel average.
It covers SN/GN/DPN configuration, 1D/2D/3D, full/reduced coupling, energy
orders 1–3, the non-CSD limit and imported-style missing stopping data.

Root's unchanged-input suite `focused_implementation_06` passed 2489 assertions.
Independent Sol review accepted the narrow mathematical correction and the
import compatibility guard. The runtime receipt and source hashes are preserved
in the external DPA validation campaign.

## Units and remaining limits

After source normalization, the returned score has units
`MeV cm^(3-dimension)/g per source particle`. Multiplying by the actual material
density and the Cartesian cell measure and summing gives scored MeV per source
particle. No second division by source normalization is required.

This correction aligns soft scoring with the existing discrete CSD polynomial.
It does not establish higher-order hard-collision energy moments, positivity,
outgoing particle/energy currents, numerical convergence, or physical accuracy.
Energy conservation requires an independent source/deposition/escape/cutoff
ledger; deposition alone is insufficient.
