"""
    generate_cross_sections(cross_sections::Cross_Sections)

Generate coupled multigroup cross sections. In addition to the historical totals, every
interaction/type/incoming/outgoing contribution is stored in generic response channels using
`quantity|Interaction/type/in=tag/out=tag` keys. This lets the HTS add-on score process-resolved
energy and charge without reconstructing physics from a collapsed total coefficient.
"""
function generate_cross_sections(cross_sections::Cross_Sections)

    L = cross_sections.get_legendre_order()
    group_structure = cross_sections.get_group_structure()
    Nmat = cross_sections.get_number_of_materials()
    Npart = cross_sections.get_number_of_particles()
    materials = cross_sections.get_materials()
    particles = cross_sections.get_particles()
    interactions = cross_sections.get_interactions()
    ρ = zeros(Nmat)
    state_of_matter = Vector{String}()
    Z = Vector{Vector{Int64}}(undef,Nmat)
    ωz = Vector{Vector{Float64}}(undef,Nmat)
    I_eff = fill(NaN,Nmat)
    for n in range(1,Nmat)
        ρ[n] = materials[n].get_density()
        Zn = materials[n].get_atomic_numbers()
        ωn = materials[n].get_weight_fractions()
        keep = findall(>(0.0),ωn)
        Z[n] = Zn[keep]
        ωz[n] = ωn[keep]
        Iov = materials[n].get_mean_excitation_energy()
        if !ismissing(Iov)
            I_eff[n] = Iov/(1e6*0.510999)
        end
        push!(state_of_matter,materials[n].get_state_of_matter())
    end
    A,atpercentA = material_isotopic_composition(materials,Z)
    isotopes = unique_isotopes(Z,A)
    for interaction in interactions
        initialize_dispatch(interaction,particles,isotopes)
    end

    E₁ = Vector{Float64}(undef,Npart)
    Ec = Vector{Float64}(undef,Npart)
    Ng = Vector{Int64}(undef,Npart)
    Eᵇ = Vector{Vector{Float64}}(undef,Npart)
    for n in range(1,Npart)
        if haskey(group_structure,particles[n])
            Eᵇ[n] = group_structure[particles[n]]
        elseif haskey(group_structure,"default")
            Eᵇ[n] = group_structure["default"]
        else
            error("Group structure is not defined for $(particles[n].type).")
        end
        E₁[n] = (Eᵇ[n][1]+Eᵇ[n][2])/2
        Ec[n] = Eᵇ[n][end]
        Ng[n] = length(Eᵇ[n])-1
    end

    interaction_interdependances(interactions,particles)

    multigroup_cross_sections = Array{Multigroup_Cross_Sections}(undef,Npart,Nmat)
    for i in range(1,Npart), n in range(1,Nmat)
        mcs = Multigroup_Cross_Sections(Ng[i])
        Σt = zeros(Ng[i])
        Σa = zeros(Ng[i])
        Σe = zeros(Ng[i]+1)
        Σc = zeros(Ng[i]+1)
        Sb = zeros(Ng[i]+1)
        S = zeros(Ng[i])
        T = zeros(Ng[i])
        for j in range(1,Npart)
            Σsl = zeros(Ng[i],Ng[j],L+1)
            for interaction in interactions
                for pin in interaction.get_in_particles(), pout in interaction.get_out_particles()
                    if pin != get_type(particles[i]) || pout != get_type(particles[j])
                        continue
                    end
                    for type in interaction.get_types(pin,pout)
                        println(
                            "\n Interaction: $(typeof(interaction)) | Type: $(type) | " *
                            "Incoming particle: $(get_type(particles[i])) | " *
                            "Outgoing particle: $(get_type(particles[j]))",
                        )
                        @time Σsli,Σti,Σai,Σei,Σci,Sbi,Si,Ti = multigroup(
                            Z[n],ωz[n],ρ[n],state_of_matter[n],Eᵇ[i],Eᵇ[j],L,
                            interaction,type,particles[i],particles[j],particles,interactions,
                            I_eff[n],A[n],atpercentA[n],
                        )
                        Σsl .+= Σsli
                        Σt .+= Σti
                        Σa .+= Σai
                        Σe .+= Σei
                        Σc .+= Σci
                        Sb .+= Sbi
                        S .+= Si
                        T .+= Ti

                        process_key = string(
                            nameof(typeof(interaction)),"/",type,
                            "/in=",get_tag(particles[i]),
                            "/out=",get_tag(particles[j]),
                        )
                        add_response_channel!(mcs,string("total|",process_key),Σti)
                        add_response_channel!(mcs,string("absorption|",process_key),Σai)
                        add_response_channel!(
                            mcs,string("energy-deposition|",process_key),Σei,
                        )
                        add_response_channel!(
                            mcs,string("charge-deposition|",process_key),Σci,
                        )
                        add_response_channel!(
                            mcs,string("boundary-stopping-power|",process_key),Sbi,
                        )
                        add_response_channel!(
                            mcs,string("stopping-power|",process_key),Si,
                        )
                        add_response_channel!(
                            mcs,string("momentum-transfer|",process_key),Ti,
                        )
                    end
                end
            end
            mcs.set_scattering(Σsl)
        end
        mcs.set_total(Σt)
        mcs.set_absorption(Σa)
        mcs.set_boundary_stopping_powers(Sb)
        mcs.set_stopping_powers(S)
        mcs.set_momentum_transfer(T)
        mcs.set_energy_deposition(Σe)
        mcs.set_charge_deposition(Σc)

        for (quantity,reference) in (
            ("total",Σt),
            ("absorption",Σa),
            ("energy-deposition",Σe),
            ("charge-deposition",Σc),
            ("boundary-stopping-power",Sb),
            ("stopping-power",S),
            ("momentum-transfer",T),
        )
            channel_keys = get_response_channel_keys(mcs;quantity=quantity)
            if !isempty(channel_keys)
                reconstructed = sum_response_channels(mcs,quantity)
                isapprox(reconstructed,reference;rtol=1.0e-10,atol=1.0e-14) || error(
                    "Process-response channels do not reconstruct $(quantity) totals.",
                )
            end
        end
        multigroup_cross_sections[i,n] = mcs
    end

    cross_sections.set_energy(E₁)
    cross_sections.set_cutoff(Ec)
    cross_sections.set_number_of_groups(Ng)
    cross_sections.set_multigroup_cross_sections(multigroup_cross_sections)
    cross_sections.set_energy_boundaries(Eᵇ)
    return nothing
end
