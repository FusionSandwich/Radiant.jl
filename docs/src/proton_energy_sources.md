# Proton source energy moments and discrete balance

These APIs support software verification of the Cartesian SN CSD path. Physical
material qualification, nonnegative transport spectra and the strict energy
ledger remain separate validation requirements.

## Public volume source representation

`Energy_Moment_Volume_Source(base, moments)` wraps an
`Anisotropic_Volume_Source`. Its array has dimensions
`(voxel, ascending energy group, angular coefficient, energy mode)`.
The first energy mode must equal the base source's group-integrated `values`.
Use one to four explicitly supplied energy modes and set the solver's DG energy
order to the same count. Sources are constant in space within each source voxel.

For intragroup coordinate `u=(2E-Ehi-Elo)/(Ehi-Elo)`, the convention is

```math
Q_k=\int_{E_{lo}}^{E_{hi}}q(E)\sqrt{2k+1}P_k(u)\,dE,
\qquad
q(E)=\frac{1}{\Delta E}\sum_k Q_k\sqrt{2k+1}P_k(u).
```

Source boundaries are ascending and in **eV**. Native transport groups are
descending and in **MeV**. Projection reverses the group indices; the local
coordinate still increases with energy, so odd modes retain their signs.
Signed higher moments are valid. Positivity is checked on the reconstructed
polynomial at its endpoints and interior stationary points. The projection
receipt retains raw minima and a scale-dependent floating-point roundoff bound
for a numerically unresolved zero. Coefficients are preserved.

`source_energy_moments(edges, coefficients; energy_order)` integrates polynomials
analytically. The coefficient array uses the same first three dimensions and
ascending powers of `u`; differential input densities are per eV.

```julia
edges = [1e6, 5e6, 10e6]
coefficients = zeros(1, 2, 1, 3)
coefficients[:, :, :, 1] .= 2e-7
coefficients[:, :, :, 2] .= 3e-8
coefficients[:, :, :, 3] .= 4e-8
moments = source_energy_moments(edges, coefficients; energy_order=3)
base = Anisotropic_Volume_Source(
    proton, [1], [voxel_volume_cm3], edges, :isotropic,
    moments[:, :, :, 1],
    Source_Normalization(basis=:per_source_particle, source_hash="polynomial-example"),
)
source = Energy_Moment_Volume_Source(base, moments)
Radiant.set_scheme(solver, "E", "DG", 3)
Radiant.add_source(sources, source)
Radiant.build(sources)
```

The existing source normalization applies once. Projection supports Cartesian
1D, 2D and 3D SN CSD-family solvers, with both full and reduced coupling. Ordinate
data must be faithfully represented by the selected angular basis. Moment input
must declare the native angle-integrated zeroth-moment convention. Missing modes,
incompatible order and unsupported angular conventions raise errors.

`get_volume_source_energy_MeV(source; physical=false)` returns the represented
first energy moment per ascending group, summed over voxels and angle. It uses
`(mid_eV*Q0 + width_eV*Q1/(2sqrt(3)))/1e6`; a one-mode source has zero Q1.
With `physical=true`, the declared physical rate and symmetry scaling apply once.
This estimates declared source energy and does not score transported deposition.

Higher-order boundary sources and serialization through the legacy source
interchange are unsupported. The legacy boundary source remains groupwise
constant in energy. Mapped geometry utilities are not connected to this native
proton source and transport path.

## Discrete energy balance

`proton_discrete_energy_balance(binding, geometry, solvers, sources, flux;
injected_energy_MeV=missing)` requires a retained, converged, single-generation
native capture. Its bounded adapter supports 1D Cartesian pure CSD, DG space and
energy, void boundaries, volume sources and zero electromagnetic field.

For energy order at least two, the native energy test includes the first moment.
The transport kernel interpolates stopping from group endpoints. Its stopping
contraction is `S0*Phi0 + S1*Phi1`, where
`S0=(Shi+Slo)/2` and `S1=(Shi-Slo)/(2sqrt(3))`.
The diagnostic reports this contraction separately from the existing midpoint
process response. Energy DG1 has only the constant test: it enforces particle
balance, and its midpoint-weighted interface balance is a discretization
diagnostic rather than a physical deposition score.

The returned fields include group residuals, represented source energy,
independently supplied injection differences, signed cutoff observations and
stored versus reconstructed cutoff differences. A small discrete residual does
not establish nonnegative spectra or physical energy conservation. The ordinary
deposition scorer and strict accounting report retain their definitions and
thresholds. Physical validation, publication readiness and PR readiness remain
false in this diagnostic.
