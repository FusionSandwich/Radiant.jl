# A1 electronic stopping implementation handoff

Date: 2026-10-08. Implementer: `/root/alpha_stopping`. Root owns runtime execution,
acceptance, repository includes, commits, and external publication.

## Source identity and attribution

Candidate branch: `JS/alpha-transport-20261008`, based on exact upstream
`2dbc7514f452b889c5708b85c44c8796c487e0da`.

Adapted selectively from zmeskma's
[PR 127](https://github.com/CBienvenue/Radiant.jl/pull/127), immutable author head
`7054f9ffc5ce107d29fde6a16249aebad1f5566f`. The integration does not establish the
author's reported numerical accuracy or reproduce ASTAR qualification.

Author source was read through GitHub UTF-8 file-fetch tools. Git blob SHA-1 identities:

| Author path | Git blob SHA-1 |
| --- | --- |
| `src/cross_sections/bethe.jl` | `1a087b5ae771da74542c3af536569d087d4e0725` |
| `src/cross_sections/fermi_density_effect.jl` | `e838bc72c342f16c11c8cfc9cb9ae6ef082ebf5f` |
| `src/cross_sections/low_energy_sp.jl` | `da38b372afcff199fe747be8124752fe00089ae7` |
| `src/structures/Inelastic_Collision.jl` | `7d642634e14c21041250773a8edddcf553d47a05` |

All four are text sources, approximately 81,000 characters altogether, below the
preflight's 2 MiB author-source budget. No binary coefficient dataset was acquired.
No package installation, environment resolution, SSH, push, PR, or commit was performed.
Existing source was unchanged until root's `BASELINE_CAPTURED` message. Root identified
the successful native baseline receipt as
`C:\rhts-alpha\dpa\validation\alphas\results\upstream_baseline_02\receipt.json`.

## Changed paths and integration dependency

- `src/cross_sections/bethe.jl`
- `src/cross_sections/fermi_density_effect.jl`
- `src/cross_sections/low_energy_sp.jl` (new)
- `src/cross_sections/plasma_energy.jl` (shared validation/composition helpers appended)
- `src/structures/Inelastic_Collision.jl`
- `test/alpha_stopping.jl` (new)
- `docs/src/api_inelastic_collision.md`
- `qualification/A1_IMPLEMENTATION.md` (this handoff)

Root must include `low_energy_sp.jl` after `bethe.jl` in `src/Radiant.jl`. No additional
exports are needed under the existing qualified/object-method convention. The actual
include loop loads `structures/` before `cross_sections/`; particle type identities
exist before `ECUT_DEFAULTS` is evaluated. Constructor references to those constants
are resolved when the constructor is called after module loading.

## Existing API and new semantics

The old mass-fraction `bethe(Z,w,rho,E,particle,...;atomic_weights=...,kwargs...)`
entry point remains. It derives atom fractions and nuclei density from one common set
of elemental/isotope masses and calls the new author-style overload
`bethe(Z,w,atom_fractions,rho,N_atoms,E,particle,...)`.

The old `fermi_density_effect(Z,w,rho,E,state,type,mass_ratio=1,I_eff=NaN;
atomic_weights=...)` convention remains. The author's proposed positional method had
the exact same Julia signature while changing `w,rho` to `atom_fractions,N_e`.
The canonical implementation is therefore explicitly named
`fermi_density_effect_from_electrons(Z,atom_fractions,N_e,E,state,type,mass_ratio,I)`.

The old `shell_correction(Z,w,gamma,particle)` remains a mass-fraction wrapper, with
an optional isotope-mass keyword. `shell_correction_from_atoms` is the explicit
atom-fraction implementation. The old five-positional-argument density/mass-fraction
`barkas_correction` remains a wrapper; the author's six-positional-argument method
accepts atom fractions, electron density, and an already resolved excitation energy.

`Inelastic_Collision` adds per-particle `set/get_ecut`, `set/get_emax` (including the
author's `"f"` option for protons/alphas), low-energy model/source selectors, and
optional Bragg additivity. Defaults remain `analytic`, `srim`, and Bragg disabled.
Alpha default cutoff/maximum are 5/0.8 MeV total kinetic energy. Numeric setter values
must be finite with `0 < Emax < Ecut`. Changing the cutoff can still reach a shell-table
coverage error, which is reported explicitly. When `"f"` is selected, the resolved
material-dependent maximum is checked during evaluation.

The optional Bragg mode evaluates each isolated elemental model at the supplied mass
density and sums by mass fraction. It uses each elemental electron density and I,
including overriding a compound I override; it does not preserve the composition
model's common electron density. This is an alternative approximation.

## Derivations and input domains

For mass fractions `w_j` and masses `A_j`, let `u_j = w_j/A_j`,
`a_j = u_j/sum(u)`, and `N_atoms = rho*N_A*sum(u)`. Then

`N_e = N_atoms*sum(a_j*Z_j) = rho*N_A*sum(w_j*Z_j/A_j)`.

Oscillator probabilities are `a_j*n_jk/Zbar`, where
`Zbar = sum(a_j*Z_j)` and `sum_k(n_jk)=Z_j`; their total is one. The shell correction
uses elemental electron fractions `a_j*Z_j/Zbar`. Barkas uses the single oscillator
sum with `CB=max(1,Zbar/10)` and the already resolved material I; it does not multiply
the entire compound oscillator sum by a second `sum(w/Z)` factor.

For projectile mass ratio `r=M/m_e` and kinetic energy `E` in electron rest energies,
`gamma=1+E/r` and `beta*gamma=sqrt(2*E*r+E^2)/r`. Sternheimer consequently uses
`X=log10(beta*gamma)`. The normalization matters for heavy particles and cancels
only for an electron/positron. The rounded factor `4.606` was replaced by `2*log(10)`.

The analytic extension is `exp(A-B*log(t)^(3/2))` above `t=1` and
`1.5*exp(A)*sqrt(t)*(1-t/3)` below it, with `t=E/Emax`. For cutoff value `Sc` and slope
`Sc_prime`, the value/slope join gives

`B = -2*Ecut*Sc_prime / (3*Sc*sqrt(log(Ecut/Emax)))`,

`A = log(Sc) + B*log(Ecut/Emax)^(3/2)`.

The lower piece has the same value and zero derivative at `t=1`. The cutoff derivative
uses the author's forward finite difference with `h=0.00005` in internal energy units.
A finite positive anchor and a nonpositive cutoff slope are required. An increasing
cutoff branch cannot match the prescribed earlier maximum and raises `DomainError`.

Composition arrays must be nonempty with matching lengths, Z in 1:100, finite
nonnegative mass fractions of positive sum, and finite positive masses and density.
Atom fractions must be finite, nonnegative, and normalized. The supplied atom fractions
and nuclei density must agree with the mass composition and isotope masses to relative
tolerance 1e-8. Energy must be finite and nonnegative; zero returns zero stopping.
I must be finite and positive, with NaN reserved as the existing default sentinel.
Canonical compound I values remain; generic mixture log-I weights use electron
fractions, including isotope-dependent atom fractions.

## Physical and data gates still open

The original effective-charge prescription remains in total electronic stopping.
Catastrophic heavy-particle collision integrals still use the incident bare charge.
No blanket charge correction was applied to differential kernels or nuclear scattering.
At sufficiently small cutoffs, `Sc` may exceed `Stot`; `checked_soft_stopping` rejects
negative or nonfinite components/residuals, including tiny negative residuals, without
clipping. This exposes an incompatible hard/soft partition; it does not resolve it.

The optional empirical `icru49` branch requires `data/low_energy_ion.jld2`. Absence
produces an explicit error. No data file was acquired or qualified. The branch rejects
energies below 10 keV total ion energy, validates coefficient entries, and retains the
explicitly rejects unsupported empirical species/Z instead of substituting analytic stopping. Nominal coefficient coverage,
phase conventions, fit quality, and SRIM/IAEA/ASTAR accuracy remain data/review gates.
The author's empirical maximum-energy fit and its compound extrapolation are unqualified.

No soft electronic second-moment API was added: no shared consistent second-moment
partition exists between the low-energy total stopping approximation and the current
catastrophic differential model. A standalone Bohr-like coefficient would conceal that
missing relationship. Electronic straggling and effective-charge hard/soft consistency
remain unsupported scientific work.

## Verification prepared for root

The implementer did not launch Julia. Static `git diff --check` on owned tracked paths
passed (Git only reported its ordinary LF/CRLF conversion notices).

Root command, using its existing depot, Julia 1.10.12, one Julia thread and one BLAS thread:

```powershell
$env:JULIA_DEPOT_PATH='C:\rhts-sandbox\julia-depots\1.10.12'
$env:OPENBLAS_NUM_THREADS='1'
& 'C:\rhts-sandbox\julia\julia-1.10.12\bin\julia.exe' --startup-file=no --threads=1 --project='C:\rhts-alpha\Radiant.jl' 'C:\rhts-alpha\Radiant.jl\test\alpha_stopping.jl'
```

Prepared tests include independent MeV relativistic maximum-transfer calculations;
equal-beta-gamma density-effect invariance across electron/proton/alpha masses with
independent asymptotic values; unscreened Bethe oracles with explicit isotope electron
densities and proton/electron regressions; hydrogen's closed-form Fano oscillator root;
mixed-element Barkas weighting and excitation-energy sensitivity; analytic value/slope
joins at Emax/Ecut; invalid NaN/Inf/negative controls; inconsistent compositions; and
negative/nonfinite soft-stopping rejection. Root runtime receipts and a separate
adversarial reviewer must determine acceptance; this handoff claims implementation
readiness, not test passage or physical qualification.
