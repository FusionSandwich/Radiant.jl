# Hard electronic collision moment foundation

Owner: bounded scientific implementation by the alpha elastic agent, following
root suite03 acceptance of the preceding candidate. Root owns execution and
independent acceptance. This change owns only the heavy collision kernel source,
`test/alpha_collision_moments.jl`, and this receipt. No software acquisition,
Julia execution, SSH, commits or external submission was performed by this agent.

## Derivation and units

The existing bare-charge hard kernel is

`dσ/dW = A [W^-2 - b W^-1 + c]`,
`b = β²/Wmax`, `c = (1-β²)/(2 M²)`,
`A = 2π r_e² z² Zi/β²`.

Here M is the projectile/electron mass ratio, Zi the electron population, and W
the transfer in electron-rest-energy units. Let `d=h-l` and `L=log1p(d/l)` for
positive ordered bounds. Exact raw integrals before multiplying by A are

- n=0: `d/(l h) - b L + c d`;
- n=1: `L - b d + c d(h+l)/2`;
- n=2: `d [1 - b(h+l)/2 + c(h²+hl+l²)/3]`.

The source evaluates n=0's leading term as `(d/h)/l`, avoiding an unnecessary
product of the bounds. It removes the previous n=0 thin-interval `d/h²`
approximation and n=1 truncated logarithm. Shared logarithms use `log1p`, with
`log(h)-log(l)` only when the finite bounds' ratio overflows. A common interval
helper clips the upper endpoint to `min(Ei,Wmax)`. Negative, reversed and
nonfinite bounds error explicitly. Number and first moments on a nonempty
interval beginning at zero diverge and error; n=2 is finite there. An ordered
interval wholly above support is empty. Unsupported moment orders error even
when the physical interval would be empty.

Returned moment units are `cm² (m_e c²)^n`: n=0 is a cross section; n=1 is cross
section times electron rest energy; n=2 is cross section times electron rest
energy squared. The n=2 primitive is a **raw hard collision moment**, not a
central variance, soft straggling coefficient or released diffusion term. A
Fokker–Planck factor one half is not included in the moment coefficient.

## Shared behavior and tests

Per-subshell n=0/1 APIs and analytical feed number/first moments share the new
primitive. Singular and angular-correction feed integrals share its interval
validation and logarithm. The regular DCS partition remains signed and is
evaluated directly, avoiding subtraction from the singular term; negative
regular terms are not clipped. Native Proton/Alpha are the supported species;
other particles are rejected explicitly. Stable cache kinematics use
`tau=Ei/M`, `β²=tau(tau+2)/(1+tau)²`, and the corresponding Wmax expression.

Tests cover Alpha/Proton at 1, 5, 20 and 500 MeV, with an independent invariant
Wmax expression, 256-bit monomial antiderivatives, and 96-node Gauss–Legendre
integration after a log-W transformation. They test n=0/1/2, cutoff additivity,
relative thin widths 1e-4 through 1e-9, the former 1e-6 branch, endpoint clipping,
zero-transfer divergence, invalid values, unsupported orders/species, and signed
regular partitions. Feed comparisons use the actual representable transfer
bounds `Ei-Ef`; the existing outgoing-energy API cannot recover transfer widths
lost when its Float64 energy endpoints were initially rounded.

## Remaining qualification limits

This establishes exact finite moments of the specified existing hard kernel.
It does not validate physical hard/soft matching, soft energy-loss variance,
electron binding approximations, angular transport, or a diffusion release.
The hard kernel uses bare projectile charge squared, while the stopping model
can use effective charge; that mismatch remains unresolved. No claim of a
physically closed stopping/straggling model follows from these kernel tests.

At handoff, the implementation agent performed static diff review only. Root's
captured test receipt is authoritative for execution outcomes and final bytes.
