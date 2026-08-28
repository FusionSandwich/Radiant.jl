"""
    Transport_Balance

Particle, energy, and charge ledger for one transported particle or for a coupled system.

The structure stores extensive quantities in one declared source-normalization basis. It does not
apply source-rate or symmetry factors internally. Production, rest-mass exchange, and cutoff terms
are explicit so coupled photon/electron/positron calculations can be balanced without hiding sinks.
"""
struct Transport_Balance
    label::String
    injected_particles::Float64
    produced_particles::Float64
    absorbed_particles::Float64
    leaked_particles::Float64
    cutoff_particles::Float64
    injected_kinetic_energy_MeV::Float64
    produced_kinetic_energy_MeV::Float64
    rest_mass_exchange_MeV::Float64
    deposited_energy_MeV::Float64
    leaked_energy_MeV::Float64
    cutoff_energy_MeV::Float64
    injected_charge::Float64
    deposited_charge::Float64
    leaked_charge::Float64

    function Transport_Balance(;
        label::AbstractString = "coupled",
        injected_particles::Real = 0.0,
        produced_particles::Real = 0.0,
        absorbed_particles::Real = 0.0,
        leaked_particles::Real = 0.0,
        cutoff_particles::Real = 0.0,
        injected_kinetic_energy_MeV::Real = 0.0,
        produced_kinetic_energy_MeV::Real = 0.0,
        rest_mass_exchange_MeV::Real = 0.0,
        deposited_energy_MeV::Real = 0.0,
        leaked_energy_MeV::Real = 0.0,
        cutoff_energy_MeV::Real = 0.0,
        injected_charge::Real = 0.0,
        deposited_charge::Real = 0.0,
        leaked_charge::Real = 0.0,
    )
        values = Float64[
            injected_particles,
            produced_particles,
            absorbed_particles,
            leaked_particles,
            cutoff_particles,
            injected_kinetic_energy_MeV,
            produced_kinetic_energy_MeV,
            rest_mass_exchange_MeV,
            deposited_energy_MeV,
            leaked_energy_MeV,
            cutoff_energy_MeV,
            injected_charge,
            deposited_charge,
            leaked_charge,
        ]
        if any(x -> !isfinite(x), values)
            error("Transport-balance entries must be finite.")
        end
        nonnegative_indices = (1,2,3,4,5,6,7,9,10,11)
        if any(values[index] < 0.0 for index in nonnegative_indices)
            error("Particle and kinetic-energy ledger entries must be nonnegative; rest_mass_exchange_MeV is the only signed energy term.")
        end

        return new(
            String(label),
            values[1], values[2], values[3], values[4], values[5],
            values[6], values[7], values[8], values[9], values[10], values[11],
            values[12], values[13], values[14],
        )
    end
end

"""
    get_particle_residual(this::Transport_Balance)

Return injected plus produced particles minus absorption, leakage, and cutoff.
"""
function get_particle_residual(this::Transport_Balance)
    return this.injected_particles + this.produced_particles -
           this.absorbed_particles - this.leaked_particles - this.cutoff_particles
end

"""
    get_energy_residual(this::Transport_Balance)

Return kinetic-energy input and production plus rest-mass exchange minus deposition, leakage, and
cutoff energy.
"""
function get_energy_residual(this::Transport_Balance)
    return this.injected_kinetic_energy_MeV + this.produced_kinetic_energy_MeV +
           this.rest_mass_exchange_MeV - this.deposited_energy_MeV -
           this.leaked_energy_MeV - this.cutoff_energy_MeV
end

"""
    get_charge_residual(this::Transport_Balance)

Return injected charge minus deposited and leaked charge.
"""
function get_charge_residual(this::Transport_Balance)
    return this.injected_charge - this.deposited_charge - this.leaked_charge
end

function _relative_balance_error(residual::Real, terms::AbstractVector{<:Real})
    scale = maximum(abs.(Float64.(terms)))
    return scale == 0.0 ? abs(Float64(residual)) : abs(Float64(residual)) / scale
end

"""
    get_relative_particle_residual(this::Transport_Balance)
"""
function get_relative_particle_residual(this::Transport_Balance)
    terms = [
        this.injected_particles + this.produced_particles,
        this.absorbed_particles + this.leaked_particles + this.cutoff_particles,
    ]
    return _relative_balance_error(get_particle_residual(this),terms)
end

"""
    get_relative_energy_residual(this::Transport_Balance)
"""
function get_relative_energy_residual(this::Transport_Balance)
    terms = [
        this.injected_kinetic_energy_MeV + this.produced_kinetic_energy_MeV + abs(this.rest_mass_exchange_MeV),
        this.deposited_energy_MeV + this.leaked_energy_MeV + this.cutoff_energy_MeV,
    ]
    return _relative_balance_error(get_energy_residual(this),terms)
end

"""
    is_balanced(this::Transport_Balance; particle_rtol=1e-8, energy_rtol=1e-8, charge_atol=1e-12)

Return `true` when the particle, energy, and charge ledgers satisfy the supplied tolerances.
"""
function is_balanced(
    this::Transport_Balance;
    particle_rtol::Real = 1.0e-8,
    energy_rtol::Real = 1.0e-8,
    charge_atol::Real = 1.0e-12,
)
    return get_relative_particle_residual(this) ≤ particle_rtol &&
           get_relative_energy_residual(this) ≤ energy_rtol &&
           abs(get_charge_residual(this)) ≤ charge_atol
end
