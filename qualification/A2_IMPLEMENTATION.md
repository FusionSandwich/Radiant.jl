# A2 alpha elastic and recoil transport mechanics

Base: `2dbc7514f452b889c5708b85c44c8796c487e0da`, branch
`JS/alpha-transport-20261008`. Implemented by the A2 elastic agent under the
root's bounded assignment. Root captures runtime evidence and independently
accepts the candidate. This agent did not install software, run Julia, submit
remote work, or commit changes.

## Defined outputs and event convention

Defaults include Proton→Proton S/P, Alpha→Alpha S/P and Alpha→Proton P.
S produces one scattered projectile of the incoming species. P produces the
target recoil, with outgoing Proton selecting H1 and outgoing Alpha selecting
He4. A hydrogen recoil cannot enter an Alpha output channel, and a He4 recoil
cannot enter a Proton output channel. Configuration changes refresh incoming
and outgoing species lists.

For identical Alpha–He4 (and Proton–H1) events, use the existing LAW=5 event
hemisphere `0 ≤ mu_cm ≤ 1`. S produces the forward CM member once; P produces
the opposite CM member once. Inclusive same-species production is the sum S+P:
two particles per event. The total cross section counts each event once. S by
itself is never multiplied by two. LAW=5 identical-particle symmetry and spin
weight remain in the event DCS, independently of outgoing-particle multiplicity.
Native Alpha requires spin zero and the matching identical flag; inconsistent
LAW=5 metadata errors explicitly.

## Independent derivation and acceptance oracles

Let E be incident kinetic energy, m the projectile nuclear mass, M the stationary
target mass, and mu the outgoing projectile CM cosine. Lorentz transformation
and energy conservation give

`s = (m+M)^2 + 2 M E`,
`K = M E (E+2m)/s`,
`T_recoil = K(1-mu)`, and `E_scattered = E-T_recoil`.

The target lab longitudinal momentum is `gamma_cm p_star (1-mu)`, while its
transverse magnitude is `p_star sqrt(1-mu^2)`. Its lab angle therefore differs
from the projectile angle for unequal masses. P integrates the original event
DCS at mu, with the target energy and angle, rather than evaluating a projectile
map with a sign flip. Both number and outgoing energy moments use the same
event quadrature weights; group midpoint energy weighting is removed.

The tests independently reconstruct outgoing four-vectors for alpha/proton,
proton/alpha, and equal masses. Synthetic full-domain `q(mu)=1+mu/2` gives event
number `4 pi` and recoil energy moment `(10 pi/3) K` before unit conversion.
Synthetic constant identical-particle hemisphere DCS gives event number `2 pi`,
inclusive output number `4 pi`, scattered energy moment `3 E/4` per event, and
recoil energy moment `E/4` per event. S+P conserves incident kinetic energy.
Additional tests challenge high-energy H1 recoils, recoil lab angle, asymmetric
energy-bin partitions, finite/positive energy and cutoff inputs, wrong outgoing
species, wrong projectile mass, wrong identical flag, wrong spin, missing data,
and alpha energies outside supplied data coverage. Proton–H1 has its own
equal-mass regression loop.

## Qualification limits

`test/alpha_elastic.jl` uses deliberately synthetic angular distributions. These
tests can establish transport mechanics and catch the specified counterexamples;
they cannot qualify physical alpha scattering, nuclear interference, or a real
ENDF evaluation. No alpha ENDF library is present in the checkout's data directory.
Missing files, missing selected projectile data, inconsistent metadata and alpha
energy extrapolation remain explicit failures. Actual coefficient and tabulated
LAW=5 alpha data still need provenance, parser/normalization validation and physical
reference comparisons before scattering qualification.

H1 and He4 elastic targets now use the native Proton and Alpha nuclear masses,
respectively, consistently in data loading, kinematics and transport moments.
Other targets retain the pre-existing neutral atomic mass approximation, which
remains a limitation. Generic recoil species beyond H1/He4 are unsupported.
Fixed Gauss–Legendre quadrature remains inherited for group moments; sharply
forward-peaked physical angular laws still require convergence validation.

Root owns the runtime receipt and final acceptance. The A2 agent performed static
diff/whitespace review only; test outcomes must be read from the root's captured
execution evidence.
