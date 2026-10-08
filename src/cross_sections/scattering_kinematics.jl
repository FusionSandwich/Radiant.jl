"""
    momentum_cm(E::Float64, m::Float64, M::Float64;
    kinematics::String="standard")

Gives the projectile center-of-mass momentum using either the ENDF
nonrelativistic convention (`"standard"`) or relativistic two-body kinematics
(`"relativistic"`). `E`, `m`, and `M` must use compatible energy units; the
returned momentum is in the same energy/c units.
"""
function momentum_cm(E::Float64, m::Float64, M::Float64;
                     kinematics::String="standard")::Float64
    kinematics_lc = lowercase(kinematics)
    if kinematics_lc == "standard"
        return (M / (m + M)) * sqrt(2.0 * m * E)
    elseif kinematics_lc == "relativistic"
        E1_total = m + E
        s = m^2 + M^2 + 2.0 * M * E1_total
        p2 = (s - (m + M)^2) * (s - (m - M)^2) / (4.0 * s)
        return sqrt(max(p2, 0.0))
    end
    error("Unknown momentum_cm kinematics $(kinematics) (should be standard or relativistic).")
end

"""
    momentum_lab(E::Float64, m::Float64; kinematics::String="standard")

Gives the initial projectile momentum in the laboratory frame using either
nonrelativistic (`"standard"`) or relativistic kinematics. `E` and `m` must use
compatible energy units; the returned momentum is in the same energy/c units.
"""
function momentum_lab(E::Float64, m::Float64; kinematics::String="standard")::Float64
    kinematics_lc = lowercase(kinematics)
    if kinematics_lc == "standard"
        return sqrt(2.0 * m * E)
    elseif kinematics_lc == "relativistic"
        return sqrt(max((E + m)^2 - m^2, 0.0))
    end
    error("Unknown momentum_lab kinematics $(kinematics) (should be standard or relativistic).")
end

"""
Relativistic scattered projectile energy for elastic scattering.
Assumes mu = cos(theta_cm) and target at rest.

Used as the forward map for mu inversion in
`mu_cm_from_Eout_relativistic` when computing mu bounds.
"""
function scattered_energy_relativistic(E_in_mev::Float64, mu::Float64,
                                       m_mev::Float64, M_mev::Float64)::Float64
    return max(E_in_mev - recoil_energy_relativistic(E_in_mev, mu, m_mev, M_mev), 0.0)
end

"""Energy-transfer slope K: the target recoil energy is K*(1-mu_cm).

This form follows from the Lorentz boost of the two outgoing CM momenta and
avoids subtracting nearly equal rest energies at small incident energy.
"""
function elastic_transfer_slope(E::Float64, m::Float64, M::Float64)
    all(isfinite, (E, m, M)) && E > 0.0 && m > 0.0 && M > 0.0 ||
        throw(ArgumentError("Elastic kinematics require finite positive incident energy and masses."))
    return M * E * (E + 2.0*m) / ((m + M)^2 + 2.0*M*E)
end

function recoil_energy_relativistic(E::Float64, mu::Float64, m::Float64, M::Float64)
    isfinite(mu) && -1.0 <= mu <= 1.0 || throw(ArgumentError("CM cosine must be in [-1,1]."))
    return elastic_transfer_slope(E, m, M) * (1.0 - mu)
end

"""Lab cosine of the outgoing target, expressed in the projectile event CM cosine.

The target CM momentum is opposite the projectile CM momentum. Its longitudinal
lab momentum is gamma_cm*p_star*(1-mu), with transverse magnitude
p_star*sqrt(1-mu^2). At zero recoil use the limiting cosine, zero.
"""
function recoil_mu_lab_from_mu_cm(mu::Float64, E::Float64, m::Float64, M::Float64)
    recoil_energy_relativistic(E, mu, m, M) # validate inputs
    gamma_cm = (E + m + M) / sqrt((m + M)^2 + 2.0*M*E)
    return gamma_cm * sqrt((1.0 - mu) / (1.0 + mu + gamma_cm^2 * (1.0 - mu)))
end

"""Use native nuclear masses for transported H1 and He4 recoils.

Other isotopes retain the existing neutral atomic mass approximation.
"""
function elastic_target_mass(Z::Int, A::Int)
    Z == 1 && A == 1 && return Float64(get_mass(Proton()))
    Z == 2 && A == 4 && return Float64(get_mass(Alpha()))
    return get_mass(Z, A)
end

"""
Invert the affine scattered-energy map and clamp to the physical domain [-1,1].

Used to convert Ef bounds into mu_cm bounds for elastic scattering
integration.
"""
function mu_cm_from_Eout_relativistic(E_in_mev::Float64, E_out_mev::Float64,
                                      m_mev::Float64, M_mev::Float64)

    isfinite(E_out_mev) || throw(ArgumentError("Outgoing energy must be finite."))
    K = elastic_transfer_slope(E_in_mev, m_mev, M_mev)
    return clamp(1.0 - (E_in_mev - E_out_mev) / K, -1.0, 1.0)
end

"""
Lab-frame mu (cos(theta_lab)) from CM mu for elastic scattering of projectile.
"""
function mu_lab_from_mu_cm(mu_cm::Float64, E_in_mev::Float64,
                           m_mev::Float64, M_mev::Float64)
    recoil_energy_relativistic(E_in_mev, mu_cm, m_mev, M_mev) # validate inputs
    E1 = E_in_mev + m_mev
    p1 = sqrt(E_in_mev * (E_in_mev + 2.0*m_mev))
    s = (m_mev + M_mev)^2 + 2.0 * M_mev * E_in_mev
    sqrt_s = sqrt(s)
    E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt_s)
    p_star = M_mev * p1 / sqrt_s
    beta_cm = p1 / (E1 + M_mev)
    gamma_cm = 1.0 / sqrt(1.0 - beta_cm^2)

    p_par = gamma_cm * (p_star * mu_cm + beta_cm * E1_star)
    p_perp = p_star * sqrt(max(1.0 - mu_cm^2, 0.0))
    p_lab = sqrt(p_par^2 + p_perp^2)
    return p_lab > 0.0 ? p_par / p_lab : 1.0
end
