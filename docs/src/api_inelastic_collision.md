# Inelastic_Collision

## Structure
```@docs
Radiant.Inelastic_Collision
```

## Setters
```@docs
Radiant.set_interaction_types(this::Radiant.Inelastic_Collision,interaction_types)
Radiant.set_density_correction(this::Radiant.Inelastic_Collision,density_correction::String)
Radiant.set_is_shell_correction(this::Radiant.Inelastic_Collision,is_shell_correction::Bool)
Radiant.set_is_subshells_dependant(this::Radiant.Inelastic_Collision,is_subshells_dependant::Bool)
Radiant.set_scattering_model(this::Radiant.Inelastic_Collision,scattering_model::String)
Radiant.set_ecut(this::Radiant.Inelastic_Collision,particle::Type,energy::Real)
Radiant.get_ecut(this::Radiant.Inelastic_Collision,particle::Type)
Radiant.set_emax(this::Radiant.Inelastic_Collision,particle::Type,energy::Real)
Radiant.set_emax(this::Radiant.Inelastic_Collision,particle::Type,mode::String)
Radiant.get_emax(this::Radiant.Inelastic_Collision,particle::Type)
Radiant.set_low_energy_model(this::Radiant.Inelastic_Collision,model::String)
Radiant.get_low_energy_model(this::Radiant.Inelastic_Collision)
Radiant.set_low_energy_source(this::Radiant.Inelastic_Collision,source::String)
Radiant.get_low_energy_source(this::Radiant.Inelastic_Collision)
Radiant.set_is_bragg_additivity(this::Radiant.Inelastic_Collision,is_bragg_additivity::Bool)
Radiant.get_is_bragg_additivity(this::Radiant.Inelastic_Collision)
```

## Electronic stopping controls

The default remains composition-based stopping with the analytic low-energy extension.
For alpha particles, the default Bethe cutoff is 5 MeV total kinetic energy and the
analytic maximum is 0.8 MeV. Setter energies are in MeV; internal stopping routines use
electron rest energy units. All numeric setter values must be finite and positive, with
`Emax < Ecut`. `set_emax(Alpha,"f")` selects the author's empirical maximum-energy fit;
its compound extrapolation is unqualified.

Atom fractions and nuclei density are derived using the same elemental or isotope masses.
Consequently `Ne = N_atoms * sum(atom_fraction .* Z) = rho * N_A * sum(w .* Z ./ A)`.
Shell, Fano oscillator and Barkas weights use electron fractions, and the resolved material
excitation energy is shared by those terms. Canonical compound excitation energies are
retained; other mixtures use the electron-weighted logarithmic rule. Sternheimer uses
`log10(beta*gamma)` with projectile mass normalization.

The optional `icru49` branch requires the separately qualified `data/low_energy_ion.jld2`
coefficient dataset. It is unavailable when that file is absent and rejects energies below
10 keV total ion energy. `srim` selects the author's base coefficients; `iaea` selects
per-element, per-phase overrides with SRIM fallback. Author accuracy and coverage claims
have not been independently reproduced, and these controls do not establish ASTAR
qualification. Selecting this empirical model for unsupported species or elements
raises an error; it does not substitute the analytic extension.

The optional Bragg mode sums isolated elemental stopping models by mass fraction. Each
element is evaluated at the supplied mass density with its own electron density and
excitation energy. This changes the compound approximation, including a user-set compound
excitation energy.

The electronic stopping correction retains the existing effective-charge prescription.
It is not applied to nuclear elastic scattering or to every differential collision
kernel. The soft stopping residual is `Stot - Sc`; a negative or nonfinite result raises
`DomainError` without clipping. Below the Bethe cutoff, total stopping and catastrophic
stopping arise from different models, so a consistent soft electronic second moment is
unsupported. No new electronic straggling coefficient is provided.

The optional controls and oscillator corrections were adapted selectively from
[zmeskma's PR #127](https://github.com/CBienvenue/Radiant.jl/pull/127), pinned at commit
`7054f9ffc5ce107d29fde6a16249aebad1f5566f`. Compatibility wrappers, finite-value guards,
partition rejection and independent formula tests were added during integration.
