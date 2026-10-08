"""
    energy_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,
    sources::Fixed_Sources,flux::Flux,particles::Vector{<:Particle})

Calculate and extract the energy deposition. 

# Input Argument(s)
- `cross_sections::Cross_Sections`: cross section informations.
- `geometry::Geometry`: geometry informations.
- `solvers::Solvers`: solvers informations.
- `sources::Sources`: sources informations.
- `flux::Flux`: flux informations.
- `particles::Vector{<:Particle}`: particle list.

# Output Argument(s)
- `Dtot::Array{Float64}`: voxel-average energy deposition
  [MeV cm^(3-n)/g per source particle], where n is the geometry dimension.
  Multiplying by density and the n-dimensional cell measure and summing gives
  scored energy per source particle.
  Native CSD soft stopping is contracted with its represented energy moments.
  Imported libraries without a stored group-average soft stopping retain their
  authoritative legacy score; higher-order energy scoring is unqualified there.

# Reference(s)
- Morel et al. (1996), A Hybrid Multigroup/Continuous-Energy Monte Carlo Method for
  Solving the Boltzmann-Fokker-Planck Equation.

"""
function energy_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,sources::Fixed_Sources,flux::Flux,particles::Vector{<:Particle})

#----
# Extract geometry data
#----
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

#----
# Calculate the total energy deposition
#----
Dtot = zeros(Ns[1],Ns[2],Ns[3])
for part in particles

    # Extract discrete_ordinates data
    discrete_ordinates = solvers.get_method(part)
    _,isCSD = discrete_ordinates.get_solver_type()

    norm = sources.get_normalization_factor()
    D = zeros(Ns[1],Ns[2],Ns[3])

    #----
    # Extract flux data
    #----
    𝚽l = flux.get_flux(part)
    if isCSD 𝚽cutoff = flux.get_flux_cutoff(part) end

    #----
    # Extract cross sections data
    #----
    Nmat = cross_sections.get_number_of_materials()
    Ng = cross_sections.get_number_of_groups(part)
    ρ = cross_sections.get_densities()

    # Extract energy deposition cross sections
    Σe = cross_sections.get_energy_deposition(part)
    if isCSD
        # Match the same boundary-linear stopping polynomial as SN/GN transport.
        # The stored score contains a separately quadrature-averaged soft stopping
        # term. Replace that term, retaining the discrete hard-collision score.
        Saverage = zeros(Ng,Nmat)
        Sboundary = zeros(Ng+1,Nmat)
        correct_soft_score = falses(Nmat)
        index_particle = findfirst(==(part),cross_sections.get_particles())
        for m in 1:Nmat
            mcs = cross_sections.multigroup_cross_sections[index_particle,m]
            # FMAC-M/MATXS import stores an authoritative deposition coefficient
            # and boundary stopping, but not its separate group-average soft term.
            # Its decomposition cannot be reconstructed from boundary data alone.
            if !ismissing(mcs.stopping_powers)
                Saverage[:,m] = mcs.get_stopping_powers()
                Sboundary[:,m] = mcs.get_boundary_stopping_powers()
                correct_soft_score[m] = true
            end
        end
        has_energy_slope = false
        if any(correct_soft_score)
            _,orders,_ = discrete_ordinates.get_schemes(geometry,discrete_ordinates.get_is_full_coupling())
            has_energy_slope = orders[4] > 1
        end
    end

    #----
    # Energy deposition calculations
    #----
    for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])

        # In-group energy deposition
        for ig in range(1,Ng)
            D[ix,iy,iz] += Σe[ig,mat[ix,iy,iz]] * 𝚽l[ig,1,1,ix,iy,iz]
            if isCSD && correct_soft_score[mat[ix,iy,iz]]
                m = mat[ix,iy,iz]
                S0 = (Sboundary[ig,m]+Sboundary[ig+1,m])/2
                D[ix,iy,iz] += (S0-Saverage[ig,m]) * 𝚽l[ig,1,1,ix,iy,iz]
                if has_energy_slope
                    S1 = (Sboundary[ig,m]-Sboundary[ig+1,m])/(2*sqrt(3))
                    # Pure energy moments precede all nonconstant spatial moments
                    # in both full and reduced SN/GN moment layouts.
                    D[ix,iy,iz] += S1 * 𝚽l[ig,1,2,ix,iy,iz]
                end
            end
        end
        if isCSD D[ix,iy,iz] += Σe[Ng+1,mat[ix,iy,iz]] * 𝚽cutoff[1,1,ix,iy,iz] end

        # Normalization
        D[ix,iy,iz] /= ρ[mat[ix,iy,iz]] * norm

    end
    Dtot += D
end

if Ndims == 1
    return Dtot[:,1,1]
elseif Ndims == 2
    return Dtot[:,:,1]
elseif Ndims == 3
    return Dtot
end

end
