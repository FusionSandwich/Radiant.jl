"""
    custom_cross_sections(cross_sections::Cross_Sections)

Build the bounded one-group isotropic library used by analytic software-verification cases.
The custom path intentionally supports one particle only. Absorption and scattering are entered
as macroscopic cross sections per material; total is their sum. Energy and charge deposition are
zero unless a later response-specific custom interface supplies them explicitly.
"""
function custom_cross_sections(cross_sections::Cross_Sections)
    Nmat = get_number_of_materials(cross_sections)
    get_number_of_particles(cross_sections) == 1 || error(
        "Custom analytic cross sections support exactly one particle.",
    )
    length(cross_sections.custom_absorption) == Nmat || error(
        "Custom absorption requires one value per material.",
    )
    length(cross_sections.custom_scattering) == Nmat || error(
        "Custom scattering requires one value per material.",
    )

    library = Array{Multigroup_Cross_Sections}(undef,1,Nmat)
    for material in 1:Nmat
        absorption = cross_sections.custom_absorption[material]
        scattering = cross_sections.custom_scattering[material]
        isfinite(absorption) && absorption ≥ 0.0 || error(
            "Custom absorption must be finite and nonnegative.",
        )
        isfinite(scattering) && scattering ≥ 0.0 || error(
            "Custom scattering must be finite and nonnegative.",
        )

        mcs = Multigroup_Cross_Sections(1)
        Σt = [absorption + scattering]
        Σa = [absorption]
        Σsl = reshape([scattering],1,1,1)
        zeros_group = zeros(Float64,1)
        zeros_boundary = zeros(Float64,2)

        set_scattering(mcs,Σsl)
        set_total(mcs,Σt)
        set_absorption(mcs,Σa)
        set_boundary_stopping_powers(mcs,copy(zeros_boundary))
        set_stopping_powers(mcs,copy(zeros_group))
        set_momentum_transfer(mcs,copy(zeros_group))
        set_energy_deposition(mcs,copy(zeros_boundary))
        set_charge_deposition(mcs,copy(zeros_boundary))
        library[1,material] = mcs
    end

    set_number_of_groups(cross_sections,[1])
    if ismissing(cross_sections.energy_boundaries)
        set_energy_boundaries(cross_sections,[[1.0,0.0]])
    end
    boundaries = cross_sections.energy_boundaries[1]
    set_energy(cross_sections,[(boundaries[1]+boundaries[2])/2.0])
    set_cutoff(cross_sections,[boundaries[end]])
    set_multigroup_cross_sections(cross_sections,library)
    return cross_sections
end
