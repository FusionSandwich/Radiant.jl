"""
    _prepare_electromagnetic_operator(cross_sections, geometry, solver, field, ...)

Build either the historical uniform angular operator or one angular operator per Cartesian voxel.
Spatial fields are sampled before the solve and remain fixed throughout source iteration. This is
a continuously varying field within a single Radiant sweep, not a sequence of independently solved
patches.
"""
function _prepare_electromagnetic_operator(
    electromagnetic_field::Electromagnetic_Field,
    part::Particle,
    geometry::Geometry,
    is_CSD::Bool,
    Ng::Int64,
    Np::Int64,
    Ns,
    Ω,
    w,
    Ndims::Int64,
    Mn,
    Dn,
    pl,
    pm,
    Eb,
    ΔE,
    Qdims::Int64;
    lambda_tolerance::Real=1.0e-14,
)
    is_EM,𝓔,𝓑 = get_electromagnetic_field(electromagnetic_field)
    if !is_EM
        return false,zeros(Float64,Ng,Np,Np),zeros(Float64,Ng)
    end
    any(𝓔 .!= 0.0) && error(
        "Electric-field transport is not implemented; only magnetic fields are supported.",
    )
    charge = get_charge(part)
    if charge == 0
        return false,zeros(Float64,Ng,Np,Np),zeros(Float64,Ng)
    end
    is_CSD || error(
        "Magnetic-field transport requires a CSD energy discretization for charged particles.",
    )

    if !has_spatial_magnetic_field(electromagnetic_field)
        matrix,lambda = electromagnetic_scattering_matrix(
            𝓔,𝓑,charge,Ω,w,Ndims,Mn,Dn,pl,pm,Np,Ng,Eb,ΔE,Qdims,
        )
        return true,matrix,lambda
    end

    field_values = get_spatial_magnetic_field(electromagnetic_field)
    expected_shape = (3,Ns[1],Ns[2],Ns[3])
    size(field_values) == expected_shape || error(
        "Spatial magnetic field shape " * string(size(field_values)) *
        " does not match transport geometry " * string(expected_shape) * ".",
    )
    matrix = zeros(Float64,Ng,Np,Np,Ns[1],Ns[2],Ns[3])
    lambda = zeros(Float64,Ng)
    cache = Dict{NTuple{3,Float64},Tuple{Array{Float64,3},Vector{Float64}}}()
    for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
        local_field = view(field_values,:,ix,iy,iz)
        key = (local_field[1],local_field[2],local_field[3])
        if !haskey(cache,key)
            local_matrix,local_lambda = electromagnetic_scattering_matrix(
                𝓔,collect(local_field),charge,Ω,w,Ndims,Mn,Dn,pl,pm,Np,Ng,Eb,ΔE,
                Qdims,
            )
            maximum(abs.(local_lambda)) <= Float64(lambda_tolerance) || error(
                "Spatial electromagnetic operator produced a cell-dependent removal term. " *
                "The current material-only total-cross-section representation cannot carry it.",
            )
            cache[key] = (local_matrix,local_lambda)
        end
        local_matrix,local_lambda = cache[key]
        matrix[:,:,:,ix,iy,iz] .= local_matrix
        lambda .= max.(lambda,abs.(local_lambda))
    end
    return true,matrix,lambda
end

"""
    compute_flux(cross_sections::Cross_Sections, geometry::Geometry,
                 solver::SN, source::Source, electromagnetic_field=Electromagnetic_Field())

Solve the discrete-ordinates transport equation for one particle. Uniform and cell-centred spatial
magnetic fields are supported for charged-particle CSD calculations. The spatial operator is
applied independently in every voxel during each source-iteration pass.
"""
function compute_flux(
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN,
    source::Source,
    electromagnetic_field::Electromagnetic_Field=Electromagnetic_Field(),
)
    Ndims = get_dimension(geometry)
    geo_type = get_type(geometry)
    geo_type == "cartesian" || error(
        "Direct Radiant transport currently requires Cartesian geometry; use the faceted " *
        "localization/voxelization adapter for CAD meshes.",
    )
    Ns = get_number_of_voxels(geometry)
    Δs = get_voxels_width(geometry)
    mat = get_material_per_voxel(geometry)
    boundary_conditions = get_boundary_conditions(geometry)

    L = get_legendre_order(solver)
    N = get_quadrature_order(solver)
    quadrature_type = get_quadrature_type(solver)
    SN_type = get_angular_boltzmann(solver)
    Qdims = get_quadrature_dimension(solver,Ndims)

    Ω,w = quadrature(N,quadrature_type,Ndims,Qdims)
    if typeof(Ω) == Vector{Float64}
        Ω = [Ω,0*Ω,0*Ω]
    end
    Nd = length(w)
    Np,Mn,Dn,pl,pm = angular_polynomial_basis(Ω,w,L,SN_type,Qdims)
    Np_surf,Mn_surf,Dn_surf,n⁺_to_n,n_to_n⁺,pl_surf,pm_surf =
        surface_angular_polynomial_basis(Ω,w,L,SN_type,Qdims,Ndims,geo_type)

    part = get_particle(solver)
    solver_type,is_CSD = get_solver_type(solver)
    Nmat = get_number_of_materials(cross_sections)
    Ng = get_number_of_groups(cross_sections,part)
    E = Float64[]
    Eb = Float64[]
    ΔE = Float64[]
    if is_CSD
        ΔE = get_energy_width(cross_sections,part)
        E = get_energies(cross_sections,part)
        Eb = get_energy_boundaries(cross_sections,part)
    end
    isFC = get_is_full_coupling(solver)
    schemes,𝒪,Nm = get_schemes(solver,geometry,isFC)
    ω,𝒞,is_adaptive,𝒲 = scheme_weights(𝒪,schemes,Ndims,is_CSD)

    println(">>>Particle: $(get_type(part)) <<<")

    Σtot = zeros(Ng,Nmat)
    if solver_type in [4,5]
        Σtot = get_absorption(cross_sections,part)
    else
        Σtot = get_total(cross_sections,part)
    end

    Ls = maximum(pl)
    Σs = zeros(Nmat,Ng,Ng,Ls+1)
    if !(solver_type in [4,5])
        Σs = get_scattering(cross_sections,part,part,Ls)
    end

    S⁻ = zeros(Ng,Nmat)
    S⁺ = zeros(Ng,Nmat)
    S = zeros(Ng,Nmat,𝒪[4])
    if is_CSD
        Sb = get_boundary_stopping_powers(cross_sections,part)
        for n in range(1,Nmat)
            S⁻[:,n] = Sb[1:Ng,n]
            S⁺[:,n] = Sb[2:Ng+1,n]
        end
        for n in range(1,Nmat), ig in range(1,Ng)
            S[ig,n,1] = (S⁻[ig,n]+S⁺[ig,n])/2
            if 𝒪[4] > 1
                S[ig,n,2] = (S⁻[ig,n]-S⁺[ig,n])/(2*sqrt(3))
            end
        end
    end

    T = zeros(Ng,Nmat)
    ℳ = Array{Float64}(undef)
    if solver_type in [2,4]
        T = get_momentum_transfer(cross_sections,part)
        fokker_planck_type = get_angular_fokker_planck(solver)
        ℳ,λ₀ = fokker_planck_scattering_matrix(
            N,Nd,quadrature_type,Ndims,fokker_planck_type,Mn,Dn,pl,Np,Qdims,
        )
        Σtot .+= T.*λ₀
    end

    if solver_type == 6
        for n in range(1,Nmat), ig in range(1,Ng)
            Σtot[ig,n] -= Σs[n,ig,ig,1]
        end
    end

    is_EM,ℳ_EM,λ₀_EM = _prepare_electromagnetic_operator(
        electromagnetic_field,part,geometry,is_CSD,Ng,Np,Ns,Ω,w,Ndims,Mn,Dn,pl,pm,
        Eb,ΔE,Qdims,
    )
    if ndims(ℳ_EM) == 3
        for ig in range(1,Ng)
            Σtot[ig,:] .+= λ₀_EM[ig]
        end
    end

    𝒜 = get_acceleration(solver)
    gmres_restart = get_gmres_restart(solver)
    anderson_depth = get_anderson_depth(solver)

    surface_sources = get_surface_sources(source)
    volume_sources = get_volume_sources(source)
    Np_source = Int64(min(Np_surf,length(surface_sources[1,:,1])))

    ϵ_max = get_convergence_criterion(solver)
    I_max = get_maximum_iteration(solver)
    𝚽l = zeros(Ng,Np,Nm[5],Ns[1],Ns[2],Ns[3])
    if is_CSD
        𝚽cutoff = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
    end

    i_out = 1
    is_outer_convergence = false
    ϵ_out = Inf
    is_outer_iteration = false
    Ntot = 0
    ρ_in = -ones(Ng)
    if is_outer_iteration
        𝚽l⁻ = zeros(Ng,Np,Nm[5],Ns[1],Ns[2],Ns[3])
    end

    while !is_outer_convergence
        ρ_in = -ones(Ng)
        if is_CSD
            𝚽E12 = zeros(Nd,Nm[4],Ns[1],Ns[2],Ns[3])
        else
            𝚽E12 = Array{Float64}(undef)
        end

        for ig in range(1,Ng)
            Qlout = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            if !(solver_type in [4,5])
                Qlout = scattering_source(
                    Qlout,𝚽l,Σs[:,:,ig,:],mat,Np,pl,Nm[5],Ns,Ng,ig,
                )
            end
            Qlout .+= volume_sources[ig,:,:,:,:,:]

            if is_CSD
                if ig != 1
                    𝚽E12 = 𝚽E12.*ΔE[ig]/ΔE[ig-1]
                end
                ΔEg = ΔE[ig]
                Sg⁻ = S⁻[ig,:]/ΔEg
                Sg⁺ = S⁺[ig,:]/ΔEg
                Sg = S[ig,:,:]/ΔEg
                if solver_type in [2,4]
                    Tg = T[ig,:]
                else
                    Tg = Vector{Float64}()
                    ℳ = Array{Float64}(undef)
                end
            else
                ΔEg = 0.0
                Sg⁻ = Vector{Float64}()
                Sg⁺ = Vector{Float64}()
                Sg = Vector{Float64}()
                Tg = Vector{Float64}()
                ℳ = Array{Float64}(undef)
            end
            electromagnetic_group = ndims(ℳ_EM) == 3 ?
                ℳ_EM[ig,:,:] : ℳ_EM[ig,:,:,:,:,:]
            𝚽l[ig,:,:,:,:,:],𝚽E12,ρ_in[ig],Ntot = sn_one_speed(
                𝚽l[ig,:,:,:,:,:],Qlout,Σtot[ig,:],Σs[:,ig,ig,:],mat,Ndims,Nd,
                ig,Ns,Δs,Ω,Mn,Dn,Np,pl,Mn_surf,Dn_surf,Np_surf,n_to_n⁺,𝒪,Nm,
                isFC,𝒞,ω,I_max,ϵ_max,surface_sources[ig,:,:],is_adaptive,is_CSD,
                solver_type,ΔEg,𝚽E12,Sg⁻,Sg⁺,Sg,Tg,ℳ,𝒜,Ntot,is_EM,
                electromagnetic_group,𝒲,boundary_conditions,Np_source,gmres_restart,
                anderson_depth,
            )
        end

        if is_outer_iteration
            ϵ_out = norm(𝚽l-𝚽l⁻)/max(norm(𝚽l),1.0e-16)
            𝚽l⁻ = 𝚽l
        end
        if (ϵ_out < ϵ_max || i_out >= I_max) || !is_outer_iteration
            is_outer_convergence = true
            if is_CSD
                for n in range(1,Nd), ix in range(1,Ns[1]), iy in range(1,Ns[2]),
                    iz in range(1,Ns[3]), is in range(1,Nm[4]), p in range(1,Np)
                    𝚽cutoff[p,is,ix,iy,iz] += Dn[p,n]*𝚽E12[n,is,ix,iy,iz]
                end
            end
        else
            i_out += 1
        end
    end

    particle_flux = Flux_Per_Particle(part)
    add_flux(particle_flux,𝚽l)
    if is_CSD
        add_flux_cutoff(particle_flux,𝚽cutoff)
    end
    add_spectral_radius(particle_flux,ρ_in)
    return particle_flux
end
