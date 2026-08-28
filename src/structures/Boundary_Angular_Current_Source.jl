abstract type Abstract_Radiant_Source end

"""
    Boundary_Angular_Current_Source

Patch-, energy-, and direction-resolved incoming boundary source for deterministic coupling.

`angular_flux` is stored as a group-integrated angular flux on the supplied quadrature. For a
patch `k` and energy group `g`, the inward current density is

```math
j⁻[k,g] = Σ(m: Ω[m] ⋅ n[k] < 0) w[m] |Ω[m] ⋅ n[k]| ψ⁻[k,g,m].
```

The constructor rejects negative values, non-unit frames, and nonzero source values in outgoing
or tangential directions. It never clips or silently reorients a source.
"""
struct Boundary_Angular_Current_Source <: Abstract_Radiant_Source
    particle::Particle
    patch_ids::Vector{Int64}
    centroids_cm::Matrix{Float64}
    areas_cm2::Vector{Float64}
    normals::Matrix{Float64}
    tangent_1::Matrix{Float64}
    tangent_2::Matrix{Float64}
    energy_edges_eV::Vector{Float64}
    directions::Matrix{Float64}
    quadrature_weights::Vector{Float64}
    angular_flux::Array{Float64,3}
    variance::Union{Nothing,Array{Float64,3}}
    normalization::Source_Normalization
    provenance::Dict{String,String}

    function Boundary_Angular_Current_Source(
        particle::Particle,
        patch_ids::AbstractVector{<:Integer},
        centroids_cm::AbstractMatrix{<:Real},
        areas_cm2::AbstractVector{<:Real},
        normals::AbstractMatrix{<:Real},
        tangent_1::AbstractMatrix{<:Real},
        tangent_2::AbstractMatrix{<:Real},
        energy_edges_eV::AbstractVector{<:Real},
        directions::AbstractMatrix{<:Real},
        quadrature_weights::AbstractVector{<:Real},
        angular_flux::AbstractArray{<:Real,3},
        normalization::Source_Normalization;
        variance::Union{Nothing,AbstractArray{<:Real,3}} = nothing,
        provenance::AbstractDict = Dict{String,String}(),
        frame_tolerance::Real = 1.0e-8,
        source_tolerance::Real = 1.0e-14,
    )
        ids = Int64.(patch_ids)
        centroids = Float64.(centroids_cm)
        areas = Float64.(areas_cm2)
        normal_vectors = Float64.(normals)
        first_tangents = Float64.(tangent_1)
        second_tangents = Float64.(tangent_2)
        edges = Float64.(energy_edges_eV)
        ordinates = Float64.(directions)
        weights = Float64.(quadrature_weights)
        psi = Float64.(angular_flux)
        variance_array = isnothing(variance) ? nothing : Float64.(variance)
        frame_tol = Float64(frame_tolerance)
        source_tol = Float64(source_tolerance)
        if !isfinite(frame_tol) || frame_tol ≤ 0.0
            error("Frame tolerance must be finite and positive.")
        end
        if !isfinite(source_tol) || source_tol < 0.0
            error("Source tolerance must be finite and nonnegative.")
        end

        Npatch = length(ids)
        Ngroup = length(edges) - 1
        Ndir = length(weights)

        if Npatch == 0
            error("At least one boundary patch is required.")
        end
        if length(unique(ids)) != Npatch
            error("Boundary patch identifiers must be unique.")
        end
        if size(centroids) != (Npatch,3)
            error("Patch centroids must have shape (number of patches, 3).")
        end
        if length(areas) != Npatch || any(x -> !isfinite(x) || x ≤ 0.0, areas)
            error("Every patch must have a finite, positive area.")
        end
        if size(normal_vectors) != (Npatch,3) ||
           size(first_tangents) != (Npatch,3) ||
           size(second_tangents) != (Npatch,3)
            error("Patch normals and tangents must each have shape (number of patches, 3).")
        end
        if length(edges) < 2 || any(x -> !isfinite(x) || x < 0.0, edges)
            error("At least two finite, nonnegative energy boundaries are required.")
        end
        if any(diff(edges) .≤ 0.0)
            error("Energy boundaries must be strictly increasing.")
        end
        if Ndir == 0 || size(ordinates) != (Ndir,3)
            error("Directions must have shape (number of ordinates, 3).")
        end
        if any(x -> !isfinite(x) || x ≤ 0.0, weights)
            error("Quadrature weights must be finite and positive.")
        end
        if size(psi) != (Npatch,Ngroup,Ndir)
            error("Angular flux must have shape (patch, energy group, direction).")
        end
        if any(x -> !isfinite(x) || x < -source_tol, psi)
            error("Boundary angular flux must be finite and nonnegative.")
        end
        psi[abs.(psi) .≤ source_tol] .= 0.0

        if !isnothing(variance_array)
            if size(variance_array) != size(psi)
                error("Boundary-source variance must have the same shape as angular flux.")
            end
            if any(x -> !isfinite(x) || x < 0.0, variance_array)
                error("Boundary-source variance must be finite and nonnegative.")
            end
        end
        for idir in 1:Ndir
            if abs(norm(view(ordinates,idir,:)) - 1.0) > frame_tol
                error("Every source direction must be a unit vector.")
            end
        end

        for ipatch in 1:Npatch
            n = collect(view(normal_vectors,ipatch,:))
            t1 = collect(view(first_tangents,ipatch,:))
            t2 = collect(view(second_tangents,ipatch,:))

            if abs(norm(n) - 1.0) > frame_tol ||
               abs(norm(t1) - 1.0) > frame_tol ||
               abs(norm(t2) - 1.0) > frame_tol
                error("Patch normals and tangents must be unit vectors.")
            end
            if abs(dot(n,t1)) > frame_tol ||
               abs(dot(n,t2)) > frame_tol ||
               abs(dot(t1,t2)) > frame_tol
                error("Patch normal and tangent vectors must be mutually orthogonal.")
            end
            if dot(cross(t1,t2),n) < 1.0 - 10.0*frame_tol
                error("Patch frames must be right-handed: tangent_1 × tangent_2 = normal.")
            end

            for idir in 1:Ndir
                mu = dot(view(ordinates,idir,:),n)
                if mu ≥ -frame_tol && any(view(psi,ipatch,:,idir) .> source_tol)
                    error("Incoming boundary source has nonzero values in an outgoing or tangential direction.")
                end
            end
        end

        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end

        return new(
            particle,
            ids,
            centroids,
            areas,
            normal_vectors,
            first_tangents,
            second_tangents,
            edges,
            ordinates,
            weights,
            psi,
            variance_array,
            normalization,
            provenance_string,
        )
    end
end

"""
    get_incoming_current_density(this::Boundary_Angular_Current_Source)

Return inward particle current density for every patch and energy group. The returned array has
shape `(patch, group)` and retains the source normalization basis.
"""
function get_incoming_current_density(this::Boundary_Angular_Current_Source)
    Npatch, Ngroup, Ndir = size(this.angular_flux)
    current = zeros(Float64,Npatch,Ngroup)
    for ipatch in 1:Npatch
        n = view(this.normals,ipatch,:)
        for idir in 1:Ndir
            mu = dot(view(this.directions,idir,:),n)
            if mu < 0.0
                for igroup in 1:Ngroup
                    current[ipatch,igroup] += this.quadrature_weights[idir] * abs(mu) *
                                               this.angular_flux[ipatch,igroup,idir]
                end
            end
        end
    end
    return current
end

"""
    get_incoming_current(this::Boundary_Angular_Current_Source)

Return inward particle current for every patch and energy group after multiplication by patch
area. The returned values retain the source normalization basis.
"""
function get_incoming_current(this::Boundary_Angular_Current_Source)
    return get_incoming_current_density(this) .* reshape(this.areas_cm2,length(this.areas_cm2),1)
end

"""
    get_total_incoming_current(this::Boundary_Angular_Current_Source; physical=false)

Return total incoming current summed over all patches and groups. Set `physical=true` to apply the
source rate and symmetry factor.
"""
function get_total_incoming_current(this::Boundary_Angular_Current_Source; physical::Bool=false)
    total = sum(get_incoming_current(this))
    return physical ? apply_normalization(total,this.normalization) : total
end

"""
    assert_current_closure(this, expected_current; rtol=1e-10, atol=1e-12)

Verify that the integrated patch/group currents equal a reference array. An error is raised on
shape mismatch or failed closure; no rescaling is performed.
"""
function assert_current_closure(
    this::Boundary_Angular_Current_Source,
    expected_current::AbstractMatrix{<:Real};
    rtol::Real = 1.0e-10,
    atol::Real = 1.0e-12,
)
    calculated = get_incoming_current(this)
    reference = Float64.(expected_current)
    if size(reference) != size(calculated)
        error("Expected current must have shape (patch, energy group).")
    end
    for index in eachindex(calculated)
        difference = abs(calculated[index] - reference[index])
        tolerance = Float64(atol) + Float64(rtol) * abs(reference[index])
        if difference > tolerance
            error("Boundary-current closure failed at linear index $(index): calculated=$(calculated[index]), expected=$(reference[index]).")
        end
    end
    return true
end

"""
    boundary_source_from_directional_current(...)

Construct a boundary angular-flux source from directional current-density contributions. Each
input entry is interpreted as the contribution
`w_m * |Ω_m ⋅ n_k| * ψ[k,g,m]` to inward current density. The mapping is exactly conservative for
non-grazing incoming directions. A nonzero grazing contribution raises an error rather than being
clipped or divided by an arbitrarily small cosine. When supplied, `variance` is the variance of the
directional-current contribution and is transformed to angular-flux variance.
"""
function boundary_source_from_directional_current(
    particle::Particle,
    patch_ids::AbstractVector{<:Integer},
    centroids_cm::AbstractMatrix{<:Real},
    areas_cm2::AbstractVector{<:Real},
    normals::AbstractMatrix{<:Real},
    tangent_1::AbstractMatrix{<:Real},
    tangent_2::AbstractMatrix{<:Real},
    energy_edges_eV::AbstractVector{<:Real},
    directions::AbstractMatrix{<:Real},
    quadrature_weights::AbstractVector{<:Real},
    directional_current_density::AbstractArray{<:Real,3},
    normalization::Source_Normalization;
    variance::Union{Nothing,AbstractArray{<:Real,3}} = nothing,
    provenance::AbstractDict = Dict{String,String}(),
    minimum_transport_factor::Real = 1.0e-12,
    source_tolerance::Real = 1.0e-14,
)
    source_tol = Float64(source_tolerance)
    minimum_factor = Float64(minimum_transport_factor)
    if !isfinite(source_tol) || source_tol < 0.0
        error("Source tolerance must be finite and nonnegative.")
    end
    if !isfinite(minimum_factor) || minimum_factor ≤ 0.0
        error("Minimum transport factor must be finite and positive.")
    end

    current = Float64.(directional_current_density)
    current_variance = isnothing(variance) ? nothing : Float64.(variance)
    if any(x -> !isfinite(x) || x < -source_tol, current)
        error("Directional current-density contributions must be finite and nonnegative.")
    end
    current[abs.(current) .≤ source_tol] .= 0.0
    if !isnothing(current_variance)
        if size(current_variance) != size(current)
            error("Directional-current variance must have the same shape as directional current.")
        end
        if any(x -> !isfinite(x) || x < 0.0, current_variance)
            error("Directional-current variance must be finite and nonnegative.")
        end
    end

    Npatch, Ngroup, Ndir = size(current)
    if Npatch != length(patch_ids) ||
       Ngroup != length(energy_edges_eV)-1 ||
       Ndir != length(quadrature_weights)
        error("Directional current density must have shape (patch, energy group, direction).")
    end

    normal_vectors = Float64.(normals)
    ordinates = Float64.(directions)
    weights = Float64.(quadrature_weights)
    if size(normal_vectors) != (Npatch,3) || size(ordinates) != (Ndir,3)
        error("Normals and directions must have shapes (patch, 3) and (direction, 3).")
    end
    if any(x -> !isfinite(x) || x ≤ 0.0, weights)
        error("Quadrature weights must be finite and positive.")
    end
    psi = zeros(Float64,size(current))
    psi_variance = isnothing(current_variance) ? nothing : zeros(Float64,size(current))

    for ipatch in 1:Npatch, idir in 1:Ndir
        mu = dot(view(ordinates,idir,:),view(normal_vectors,ipatch,:))
        factor = weights[idir] * abs(mu)
        nonzero_mean = any(view(current,ipatch,:,idir) .> source_tol)
        nonzero_variance = !isnothing(current_variance) && any(view(current_variance,ipatch,:,idir) .> 0.0)
        nonzero = nonzero_mean || nonzero_variance
        if nonzero && mu ≥ 0.0
            error("Directional current was assigned to a non-incoming direction.")
        end
        if nonzero && factor ≤ minimum_factor
            error("Nonzero directional current lies in an unresolved grazing bin; use a conservative angular projection.")
        end
        if mu < 0.0 && factor > minimum_factor
            psi[ipatch,:,idir] .= view(current,ipatch,:,idir) ./ factor
            if !isnothing(psi_variance)
                psi_variance[ipatch,:,idir] .= view(current_variance,ipatch,:,idir) ./ factor^2
            end
        end
    end

    source = Boundary_Angular_Current_Source(
        particle,
        patch_ids,
        centroids_cm,
        areas_cm2,
        normals,
        tangent_1,
        tangent_2,
        energy_edges_eV,
        directions,
        quadrature_weights,
        psi,
        normalization;
        variance=psi_variance,
        provenance=provenance,
        source_tolerance=source_tolerance,
    )

    expected = dropdims(sum(current,dims=3),dims=3) .* reshape(Float64.(areas_cm2),length(areas_cm2),1)
    assert_current_closure(source,expected)
    return source
end
