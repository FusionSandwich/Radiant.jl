# Admission around the selected author features. These checks establish supported
# input identity and fail closed; they do not qualify the underlying approximations.
function author_nuclear_table(E, S, label::AbstractString)
    length(E) >= 2 && length(S) == length(E) || error("$label: missing or malformed energy table.")
    all(isfinite, E) && all(isfinite, S) || error("$label: nonfinite table.")
    all(diff(E) .> 0.0) && all(E .>= 0.0) && all(S .>= 0.0) ||
        error("$label: nonmonotone energies or negative coefficient.")
    return nothing
end

function author_nuclear_energy(E::Float64, grid, label::AbstractString; threshold::Bool=false)
    length(grid) >= 2 && all(isfinite, grid) && all(grid .>= 0.0) && all(diff(grid) .> 0.0) ||
        error("$label: missing or malformed support.")
    isfinite(E) && E >= 0.0 || error("$label: invalid requested energy $E.")
    E > grid[end] && error("$label: requested energy $E exceeds evaluated maximum $(grid[end]); endpoint holding is unsupported.")
    if E < grid[1]
        threshold && return false
        error("$label: requested energy $E is below verified support $(grid[1]).")
    end
    return true
end

function author_evaluation_identity(lines::Vector{String}, Z::Int, A::Int,
        particle::Particle, source::AbstractString)
    i = find_section_start(lines, 1, 451)
    isnothing(i) && error("$source: required MF1/MT451 projectile/target identity is absent.")
    (ZA, _, _, _, _, _, mat, mf, mt, _, i, _) = read_CONT(lines, i)
    validate_section_ids(mat, mf, mt, mat, 1, 451, "$source identity HEAD")
    ZA == 1000*Z + A || error("$source: target ZA mismatch.")
    (_, _, _, _, _, _, mat2, mf2, mt2, _, i, _) = read_CONT(lines, i)
    validate_section_ids(mat2, mf2, mt2, mat, 1, 451, "$source identity second CONT")
    (_, _, _, _, NSUB, _, mat3, mf3, mt3, _, _, _) = read_CONT(lines, i)
    validate_section_ids(mat3, mf3, mt3, mat, 1, 451, "$source identity third CONT")
    zap = ZAP(particle)
    isnothing(zap) && error("$source: unsupported incident particle identity.")
    NSUB == 10*zap || error("$source: projectile NSUB=$NSUB does not match incident ZAP=$zap.")
    # MF3/MF6 section target identity must agree with the evaluated target.
    seen = Set{Tuple{Int,Int}}()
    for j in eachindex(lines)
        (_,_,_,_,_,_,_,mfj,mtj,_) = parse_endf_line(lines[j])
        if mfj in (3,6) && mtj > 0 && !((mfj,mtj) in seen)
            push!(seen,(mfj,mtj))
            (ZAj,AWRj,_,_,_,_,_,_,_,_,_,_) = read_CONT(lines,j)
            ZAj == 1000*Z+A && isfinite(AWRj) && AWRj > 0.0 ||
                error("$source: MF$mfj MT$mtj target identity/mass mismatch.")
        end
    end
    return nothing
end

function author_range_stopping(S::Float64, E::Float64, Z, phase)
    isfinite(S) && S > 0.0 ||
        error("Equivalent-proton range: nonfinite/nonpositive stopping $S at energy $E m_e*c^2, elements $Z, phase $phase.")
    return S
end

function author_straggling_requested(cross_sections::Cross_Sections, particle::Particle)
    (is_proton(particle) || is_alpha(particle)) || return false # flag is a heavy-only feature; no lepton diffusion capability
    interactions = cross_sections.get_interactions()
    ismissing(interactions) && return false # legacy/custom library has no requested native model
    ptype = get_type(particle)
    return any(interactions) do interaction
        interaction isa Inelastic_Collision && interaction.is_straggling &&
            any(kv -> kv.first[1] == ptype && "S" in kv.second, interaction.interaction_types)
    end
end

function author_straggling_configuration(interaction::Interaction, particle::Particle, type::String)
    (is_proton(particle) || is_alpha(particle)) || return nothing # preserve author/main lepton stopping defaults
    if interaction isa Inelastic_Collision && interaction.is_straggling && startswith(type, "S")
        is_proton(particle) || is_alpha(particle) ||
            error("Author straggling is implemented only for Proton and Alpha; disable it explicitly for $(get_type(particle)).")
        interaction.is_subshells_dependant || error("Author straggling does not support unbound-electron configuration; disable it explicitly.")
        interaction.is_CSD && interaction.scattering_model != "BTE" ||
            error("Requested author straggling is unsupported without continuous energy transport.")
    end
    return nothing
end

function author_straggling_material(material::Material)
    any(a -> !isempty(a), material.mass_numbers) && error("Author straggling assumes natural elemental densities; explicit isotope composition is unsupported.")
    return nothing
end

function prepare_author_straggling(cross_sections::Cross_Sections, particle::Particle,
        route::Symbol, is_CSD::Bool, fast_requested::Bool=false)
    requested = author_straggling_requested(cross_sections, particle)
    ismissing(cross_sections.multigroup_cross_sections) &&
        error("Straggling admission: cross-section coefficients unavailable.")
    index = findfirst(==(particle), cross_sections.get_particles())
    isnothing(index) && error("Straggling admission: particle absent.")
    Ng = cross_sections.get_number_of_groups(particle)
    Nmat = cross_sections.get_number_of_materials()
    Tb = zeros(Ng+1, Nmat)
    for n in 1:Nmat
        tb = cross_sections.multigroup_cross_sections[index,n].get_boundary_straggling_coefficients()
        if ismissing(tb)
            requested && error("Requested straggling has missing coefficients in material $n.")
            continue # explicit no-model compatibility for legacy/custom data
        end
        length(tb) == Ng+1 && all(isfinite, tb) && all(tb .>= 0.0) ||
            error("Straggling admission: invalid boundary coefficients in material $n.")
        Tb[:,n] = tb
    end
    active = any(Tb .!= 0.0)
    if requested || active
        for material in cross_sections.get_materials()
            author_straggling_material(material)
        end
        is_proton(particle) || is_alpha(particle) ||
            error("Straggling admission: unsupported particle $(get_type(particle)).")
        is_CSD || error("Straggling admission: requested/nonzero coefficients on a non-CSD route.")
        route == :SN || error("Straggling is unsupported by $route (including GN/DPN and CP).")
    end
    active && fast_requested &&
        error("Nonzero author straggling is unsupported by the optimized SN chain; request the experimental reference path explicitly.")
    active && println(">>>Experimental author local SN straggling: UNQUALIFIED; energy-face/intergroup diffusion is incomplete.")
    return Tb
end
