# Adapted from zmeskma/Radiant.jl PR #127 at 7054f9ffc5ce107d29fde6a16249aebad1f5566f.
# Empirical accuracy claims and coefficient data have not been independently qualified.
# Reference(s)
# - Varelas, C. and Biersack, J. (1970), Reflection of energetic particles from
#   atomic or ionic chains in single crystals, Nuclear Instruments and Methods
#   79(2), 213-218. doi:10.1016/0029-554X(70)90141-2
#   Reciprocal interpolation formula: Sec. 4.2, Eq. (6), p. 216.
# - Andersen, H. H. and Ziegler, J. F. (1977), Hydrogen: Stopping Powers and
#   Ranges in All Elements, The Stopping and Ranges of Ions in Matter vol. 3,
#   Pergamon Press, New York.
# - Ziegler, J. F. (1977), Helium: Stopping Powers and Ranges in All Elemental
#   Matter, The Stopping and Ranges of Ions in Matter vol. 4, Pergamon Press,
#   New York.
# - ICRU Report 49 (1993), Stopping Powers and Ranges for Protons and Alpha
#   Particles, ICRU, Bethesda, MD.
# - Salvat, F. (2022), Bethe stopping-power formula and its corrections,
#   Physical Review A 106, 032809. doi:10.1103/PhysRevA.106.032809
#   Free-exponent low-energy parametrization: Eqs. (134)-(135); fitting window
#   and objective: Eq. (136) and the text above it. Source of the fitting
#   prescription, not of the functional form -- see icru49_stopping.

"""
    low_energy_stopping(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},
    ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,
    density_correction::String,state_of_matter::String,I_eff::Float64,Ecut::Float64,
    Emax::Float64)

Analytic continuation with a value and slope join at Ecut, or optional empirical stopping. Energies use electron rest energy units.
"""
function low_energy_stopping(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,density_correction::String,state_of_matter::String,I_eff::Float64,Ecut::Float64,Emax::Float64; atomic_weights::Union{Nothing,Vector{Float64}}=nothing,is_bethe_correction::Bool=true,is_density_correction::Bool=true,is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true,user_ecut::Float64=NaN,low_energy_model::String="analytic",low_energy_source::String="srim")

    mₑ = 0.51099895069
    (isfinite(Ecut) && isfinite(Emax) && 0 < Emax < Ecut) ||
        throw(ArgumentError("The analytic maximum must satisfy 0 < Emax < Ecut."))
    (isfinite(Ei) && Ei >= 0) || throw(ArgumentError("Energy must be finite and nonnegative."))
    Ei == 0 && return 0.0

    #----
    # Fitted ICRU 49 / Andersen-Ziegler branch
    #----
    # The unit change sits on this side of the dispatch: the fitted form works in MeV and
    # returns a mass stopping power, while the analytic extension already works in bethe's
    # units and returns a linear one.
    if low_energy_model != "analytic" && has_low_energy_data(Z,particle)
        return icru49_stopping(Z,ωz,Ei*mₑ,particle,state_of_matter,low_energy_source) * ρ/mₑ
    end

    #----
    # Estimate parameters A and B
    #----
    # The anchor energies are both ≥ Ecut, so these two calls take the corrected-Bethe
    # branch and never come back here: the recursion is one level deep by construction.
    # They pass low_energy_model = "analytic" whatever the caller selected, because under
    # the handover window the fitted form is what bethe returns at Ecut -- the blend weight
    # is zero there -- so an anchor taken with the caller's setting would anchor this
    # extension on the fitted form rather than on the Bethe branch. The correction switches
    # are forwarded unchanged, or the extension would be anchored on a different formula
    # from the one it is joined to.
    h = 0.00005
    Sc = bethe(Z,ωz,atz,ρ,N_density,Ecut,particle,density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_bethe_correction=is_bethe_correction,
        is_density_correction=is_density_correction,
        is_shell_correction=is_shell_correction,
        is_lindhard_sorensen=is_lindhard_sorensen,
        is_barkas=is_barkas,
        user_ecut=Ecut*mₑ,
        user_emax=Emax*mₑ,
        low_energy_model="analytic")
    Sc_plus_h = bethe(Z,ωz,atz,ρ,N_density,Ecut+h,particle,density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_bethe_correction=is_bethe_correction,
        is_density_correction=is_density_correction,
        is_shell_correction=is_shell_correction,
        is_lindhard_sorensen=is_lindhard_sorensen,
        is_barkas=is_barkas,
        user_ecut=Ecut*mₑ,
        user_emax=Emax*mₑ,
        low_energy_model="analytic")
    (isfinite(Sc) && Sc > 0 && isfinite(Sc_plus_h)) ||
        throw(DomainError(Sc,"The analytic stopping anchor must be finite and positive."))
    dSdE = (Sc_plus_h - Sc)/h
    dSdE <= 0 || throw(DomainError(dSdE,
        "The cutoff lies on an increasing stopping branch; no earlier analytic maximum can match this slope."))
    A = (3*log(Sc)*Sc - 2*dSdE*Ecut*log(Ecut/Emax))/(3*Sc)
    B = -2*Ecut*dSdE/(3*Sc*sqrt(log(Ecut/Emax)))

    #----
    # Compute total cross-sections using analytical formula
    #----
    t = Ei/Emax
    if t ≥ 1
        S = exp(A-B*(log(t))^1.5)
    else
        S = 1.5*exp(A)*sqrt(t)*(1-t/3)
    end
    return S
end

"""
    emax_stopping_maximum(Z::Vector{Int64},atz::Vector{Float64},particle::Particle)

Author-supplied empirical maximum-energy fit [MeV]. The fit and its compound extension are unqualified here.
"""
function emax_stopping_maximum(Z::Vector{Int64},atz::Vector{Float64},particle::Particle)
    Zbar = sum(atz .* Z)
    if is_proton(particle) || is_antiproton(particle)
        return 0.0481 * Zbar^0.190
    elseif is_alpha(particle)
        return 0.4829 * Zbar^0.114
    else
        error("No stopping-maximum fit for the given particle.")
    end
end

"""
    low_energy_tables(particle::Particle,source::String)

Load the author-supplied SRIM base and optional IAEA override coefficients. The external coefficient dataset is required and must be qualified separately.
"""
function low_energy_tables(particle::Particle,source::String)
    (is_proton(particle) || is_alpha(particle)) ||
        throw(ArgumentError("Empirical ion stopping supports protons and alphas only."))
    file = joinpath(find_package_root(),"data","low_energy_ion.jld2")
    isfile(file) || throw(ArgumentError(
        "Optional empirical stopping requires the separately qualified low_energy_ion.jld2 dataset."))
    data = fast_load("low_energy_ion.jld2")
    key  = is_alpha(particle) ? "alpha" : "proton"
    base = data["fitted_exponent"][key]
    if source == "srim"
        override = Dict{String,Any}()
    elseif source == "iaea"
        override = data["fitted_exponent_iaea"][key]
    else
        error("Unknown '$source' low-energy coefficient source.")
    end
    # gas_elements is a fact about the elements rather than a borrowing from either table,
    # and is the same eleven either way.
    return base, override, data["fitted_exponent"]["gas_elements"]
end

"""
    low_energy_slot(table::Dict,Z::Int64,is_gas::Bool,gas_elements::Vector{Int64})

Resolve the phase slot independently for each coefficient table.
"""
function low_energy_slot(table,Z::Int64,is_gas::Bool,gas_elements::Vector{Int64})
    if is_gas
        # The haskey on "gas" is what makes this reader work against a table written before
        # that slot existed, and the slot is not restricted to the elements that are gases
        # at STP.
        return (haskey(table,"gas") && haskey(table["gas"],Z)) ? "gas" : "coefficients"
    else
        return (Z in gas_elements && haskey(table,"condensed")) ? "condensed" : "coefficients"
    end
end

"""
    has_low_energy_data(Z::Vector{Int64},particle::Particle)

Check the author model's nominal Z/species domain. This does not establish that coefficient data are installed or qualified.
"""
function has_low_energy_data(Z::Vector{Int64},particle::Particle)
    if !(is_proton(particle) || is_alpha(particle)) return false end
    return all(1 .≤ Z .≤ 92)
end

"""
    icru49_stopping(Z::Vector{Int64},ωz::Vector{Float64},E_MeV::Float64,particle::Particle,
    state_of_matter::String,source::String)

Optional reciprocal interpolation of elemental mass stopping powers [MeV cm²/g]. The input energy is in MeV. Coefficient accuracy is unqualified here.
"""
function icru49_stopping(Z::Vector{Int64},ωz::Vector{Float64},E_MeV::Float64,particle::Particle,state_of_matter::String,source::String)

    (isfinite(E_MeV) && E_MeV >= 0.01) || throw(DomainError(E_MeV,
        "Optional empirical stopping is unqualified below 10 keV total ion energy."))
    state_of_matter in ("solid","liquid","gaz") || throw(ArgumentError("Unknown material state."))
    base, override, gas_elements = low_energy_tables(particle,source)
    T = E_MeV/base["T_unit_MeV"]
    is_gas = state_of_matter == "gaz"

    S = 0.0
    for i in eachindex(Z)

        # Resolve the phase first, then the source, and ask the coverage question of the
        # slot rather than of the element -- see low_energy_slot. Each table is asked with
        # its own phase convention. Selecting a source does not switch the whole material
        # over: it is resolved per element, at the point of use, and a compound routinely
        # ends up drawing from both tables at once.
        base_slot = low_energy_slot(base,Z[i],is_gas,gas_elements)
        over_slot = low_energy_slot(override,Z[i],is_gas,gas_elements)
        coef = haskey(override,over_slot) ? get(override[over_slot],Z[i],base[base_slot][Z[i]]) :
                                            base[base_slot][Z[i]]

        (length(coef) == 5 && all(isfinite,coef)) ||
            throw(ArgumentError("Invalid empirical stopping coefficient entry."))
        a1,a2,a3,a4,a5 = coef
        (a1 > 0 && 0.30 <= a2 <= 0.80 && a3 > 0 && a4 >= 0 && a5 >= 0) ||
            throw(ArgumentError("Empirical stopping coefficients are outside the supported domain."))
        s_l = a1 * T^a2
        s_h = (a3/T) * log(1 + a4/T + a5*T)
        S += ωz[i] * s_l*s_h/(s_l + s_h)
    end
    return S
end

"""
    blend_weight(E_MeV::Float64,Ecut_MeV::Float64,E_top_MeV::Float64)

Smoothstep Bethe weight over a finite increasing energy interval; its derivative vanishes at both ends.
"""
function blend_weight(E_MeV::Float64,Ecut_MeV::Float64,E_top_MeV::Float64)
    (all(isfinite,(E_MeV,Ecut_MeV,E_top_MeV)) && 0 < Ecut_MeV < E_top_MeV) ||
        throw(ArgumentError("Blend energies must be finite with 0 < cutoff < top."))
    t = (E_MeV - Ecut_MeV)/(E_top_MeV - Ecut_MeV)
    t ≤ 0 && return 0.0
    t ≥ 1 && return 1.0
    return t^2*(3 - 2*t)
end
