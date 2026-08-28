"""
    Anisotropic_Volume_Source

Voxel-, energy-, and angle-resolved volume source used for neutron-induced photons, delayed
emissions, and other coupled-particle sources.

Supported angular representations are:

- `:isotropic`: one angle-integrated coefficient per voxel and group.
- `:ordinates`: values collocated on the supplied directions and quadrature weights.
- `:moments`: coefficients in an externally declared angular basis; the basis convention must be
  recorded in provenance.

`values` has shape `(voxel, energy group, angular coefficient)` and is never implicitly
normalized by voxel volume or source rate.
"""
struct Anisotropic_Volume_Source <: Abstract_Radiant_Source
    particle::Particle
    voxel_ids::Vector{Int64}
    voxel_volumes_cm3::Vector{Float64}
    energy_edges_eV::Vector{Float64}
    angular_representation::Symbol
    directions::Matrix{Float64}
    quadrature_weights::Vector{Float64}
    values::Array{Float64,3}
    variance::Union{Nothing,Array{Float64,3}}
    parent_reaction::Union{Nothing,String}
    normalization::Source_Normalization
    provenance::Dict{String,String}

    function Anisotropic_Volume_Source(
        particle::Particle,
        voxel_ids::AbstractVector{<:Integer},
        voxel_volumes_cm3::AbstractVector{<:Real},
        energy_edges_eV::AbstractVector{<:Real},
        angular_representation::Symbol,
        values::AbstractArray{<:Real,3},
        normalization::Source_Normalization;
        directions::AbstractMatrix{<:Real} = zeros(Float64,0,3),
        quadrature_weights::AbstractVector{<:Real} = Float64[],
        variance::Union{Nothing,AbstractArray{<:Real,3}} = nothing,
        parent_reaction::Union{Nothing,AbstractString} = nothing,
        provenance::AbstractDict = Dict{String,String}(),
        direction_tolerance::Real = 1.0e-8,
        source_tolerance::Real = 1.0e-14,
    )
        if angular_representation ∉ (:isotropic,:ordinates,:moments)
            error("Angular representation must be :isotropic, :ordinates, or :moments.")
        end

        ids = Int64.(voxel_ids)
        volumes = Float64.(voxel_volumes_cm3)
        edges = Float64.(energy_edges_eV)
        ordinates = Float64.(directions)
        weights = Float64.(quadrature_weights)
        source_values = Float64.(values)
        variance_array = isnothing(variance) ? nothing : Float64.(variance)
        direction_tol = Float64(direction_tolerance)
        source_tol = Float64(source_tolerance)

        Nvoxel = length(ids)
        Ngroup = length(edges) - 1
        Ncoefficient = size(source_values,3)

        if Nvoxel == 0 || length(unique(ids)) != Nvoxel
            error("Volume-source voxel identifiers must be nonempty and unique.")
        end
        if length(volumes) != Nvoxel || any(x -> !isfinite(x) || x ≤ 0.0, volumes)
            error("Every source voxel must have a finite, positive volume.")
        end
        if length(edges) < 2 || any(x -> !isfinite(x) || x < 0.0, edges) || any(diff(edges) .≤ 0.0)
            error("Volume-source energy boundaries must be finite, nonnegative, and increasing.")
        end
        if size(source_values,1) != Nvoxel || size(source_values,2) != Ngroup || Ncoefficient == 0
            error("Volume-source values must have shape (voxel, energy group, angular coefficient).")
        end
        if !isfinite(source_tol) || source_tol < 0.0
            error("Source tolerance must be finite and nonnegative.")
        end
        if any(x -> !isfinite(x) || x < -source_tol, source_values)
            error("Volume-source values must be finite and nonnegative.")
        end
        source_values[abs.(source_values) .≤ source_tol] .= 0.0

        if !isnothing(variance_array)
            if size(variance_array) != size(source_values)
                error("Volume-source variance must have the same shape as source values.")
            end
            if any(x -> !isfinite(x) || x < 0.0, variance_array)
                error("Volume-source variance must be finite and nonnegative.")
            end
        end

        if angular_representation == :isotropic
            if Ncoefficient != 1 || size(ordinates,1) != 0 || length(weights) != 0
                error("An isotropic source must have one coefficient and no explicit quadrature.")
            end
        elseif angular_representation == :ordinates
            Ndir = length(weights)
            if Ndir == 0 || size(ordinates) != (Ndir,3) || Ncoefficient != Ndir
                error("An ordinate source requires matching directions, weights, and angular coefficients.")
            end
            if any(x -> !isfinite(x) || x ≤ 0.0, weights)
                error("Volume-source quadrature weights must be finite and positive.")
            end
            if !isfinite(direction_tol) || direction_tol ≤ 0.0
                error("Direction tolerance must be finite and positive.")
            end
            for idir in 1:Ndir
                if abs(norm(view(ordinates,idir,:)) - 1.0) > direction_tol
                    error("Every volume-source ordinate must be a unit vector.")
                end
            end
        else
            if size(ordinates,1) != 0 || length(weights) != 0
                error("Moment sources do not accept explicit quadrature arrays.")
            end
            if !haskey(provenance,"angular_basis") && !haskey(provenance,:angular_basis)
                error("Moment-source provenance must declare an angular_basis.")
            end
        end

        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end

        reaction = isnothing(parent_reaction) ? nothing : String(parent_reaction)
        return new(
            particle,
            ids,
            volumes,
            edges,
            angular_representation,
            ordinates,
            weights,
            source_values,
            variance_array,
            reaction,
            normalization,
            provenance_string,
        )
    end
end

"""
    get_volume_source_rate(this::Anisotropic_Volume_Source; physical=false)

Return the angle-integrated source rate per energy group, summed over source voxels. Isotropic
values are interpreted as angle-integrated source densities. Ordinate values are integrated using
the supplied quadrature. Moment sources require an explicit downstream basis mapping and are not
silently reduced to a scalar rate.
"""
function get_volume_source_rate(this::Anisotropic_Volume_Source; physical::Bool=false)
    Nvoxel, Ngroup, Ncoefficient = size(this.values)
    rate = zeros(Float64,Ngroup)

    if this.angular_representation == :isotropic
        for ivoxel in 1:Nvoxel, igroup in 1:Ngroup
            rate[igroup] += this.values[ivoxel,igroup,1] * this.voxel_volumes_cm3[ivoxel]
        end
    elseif this.angular_representation == :ordinates
        for ivoxel in 1:Nvoxel, igroup in 1:Ngroup, idir in 1:Ncoefficient
            rate[igroup] += this.values[ivoxel,igroup,idir] *
                            this.quadrature_weights[idir] *
                            this.voxel_volumes_cm3[ivoxel]
        end
    else
        error("A moment source requires an explicit basis-to-scalar mapping.")
    end

    return physical ? apply_normalization(rate,this.normalization) : rate
end
