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
normalized by voxel volume or source rate. Isotropic and ordinate values must be nonnegative.
Moment coefficients may be signed; their reconstructed angular source must be checked separately.
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
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end

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
        if any(x -> !isfinite(x),source_values)
            error("Volume-source values must be finite.")
        end
        if angular_representation != :moments && any(x -> x < -source_tol,source_values)
            error("Isotropic and ordinate volume-source values must be nonnegative.")
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
            if !haskey(provenance_string,"angular_basis")
                error("Moment-source provenance must declare an angular_basis.")
            end
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

get_particle(this::Anisotropic_Volume_Source) = this.particle
get_source_normalization(this::Anisotropic_Volume_Source) = this.normalization

# Normalized Legendre polynomials, in ascending powers of the physical energy
# coordinate u=(2E-Ehi-Elo)/(Ehi-Elo). The public first increment is DG1--DG4.
const _source_energy_power_basis = (
    [1.0], [0.0,sqrt(3.0)], [-sqrt(5.0)/2,0.0,3sqrt(5.0)/2],
    [0.0,-3sqrt(7.0)/2,0.0,5sqrt(7.0)/2],
)

function _source_energy_polynomial_diagnostic(moments)
    1 <= length(moments) <= 4 || error("Energy source supports one to four explicit modes.")
    powers = zeros(Float64,4)
    for k in eachindex(moments), j in eachindex(_source_energy_power_basis[k])
        powers[j] += moments[k] * _source_energy_power_basis[k][j]
    end
    all(isfinite,powers) || error("Represented energy polynomial overflowed.")
    scale = max(maximum(abs,powers),maximum(abs,moments))
    scale == 0.0 && return (minimum=0.0,roundoff_bound=0.0)
    # Normalize before solving the derivative: finite large coefficients must
    # not overflow the quadratic discriminant and hide an interior minimum.
    powers ./= scale
    candidates = [-1.0,1.0]
    a,b,c = 3powers[4],2powers[3],powers[2]
    if a == 0.0
        b != 0.0 && -1.0 <= -c/b <= 1.0 && push!(candidates,-c/b)
    else
        discriminant = b*b-4a*c
        isfinite(discriminant) || error("Energy-source derivative discriminant is nonfinite.")
        if discriminant >= 0.0
            q = -0.5*(b+copysign(sqrt(discriminant),b))
            roots = q == 0.0 ? [-b/(2a)] : [q/a,c/q]
            # A finite near-quadratic cubic can have an unrepresentably distant
            # root. +/-Inf is outside [-1,1]; NaN cannot classify a stationary point.
            any(isnan,roots) && error("Energy-source derivative roots are undefined.")
            for u in roots
                -1.0 <= u <= 1.0 && push!(candidates,u)
            end
        end
    end
    raw_minimum = scale*minimum(((powers[4]*u+powers[3])*u+powers[2])*u+powers[1] for u in candidates)
    # This bounds a numerically unresolved sign near a tangential zero. It is
    # reported, never used to clip moments or adjust transport positivity.
    roundoff_bound = scale*(64eps(Float64)*sum(abs(moments[k]/scale)*sum(abs,_source_energy_power_basis[k]) for k in eachindex(moments)))
    isfinite(raw_minimum) && isfinite(roundoff_bound) || error("Energy-source positivity diagnostic is nonfinite.")
    return (minimum=raw_minimum,roundoff_bound=roundoff_bound)
end

_source_energy_polynomial_minimum(moments) = _source_energy_polynomial_diagnostic(moments).minimum

"""
    source_energy_moments(energy_edges_eV, local_power_coefficients; energy_order)

Exactly integrate a declared intragroup polynomial into normalized energy Legendre moments.
Input has shape `(entity, ascending energy group, angular coefficient, power)` with ascending
powers of `u=(2E-Ehi-Elo)/(Ehi-Elo)`. Coefficients are differential densities per eV. Output has
the same first three dimensions and `energy_order` modes, with
`Qk = integral_bin q(E)*sqrt(2k+1)*Pk(u) dE`. Q0 is the exact group integral.
One to four modes are supported; every polynomial power must be represented. No quadrature,
clipping, normalization or inference of missing polynomial coefficients is performed.
"""
function source_energy_moments(edges_eV::AbstractVector{<:Real},
    coefficients::AbstractArray{<:Real,4}; energy_order::Integer=size(coefficients,4))
    edges = Float64.(edges_eV)
    1 <= energy_order <= 4 || error("Energy source order must lie between one and four.")
    1 <= size(coefficients,4) <= energy_order || error("Polynomial powers exceed the declared energy order.")
    length(edges) >= 2 && all(isfinite,edges) && all(edges .>= 0.0) && all(diff(edges) .> 0.0) ||
        error("Energy boundaries must be finite, nonnegative and strictly increasing.")
    size(coefficients,2) == length(edges)-1 || error("Polynomial energy groups do not match the boundaries.")
    size(coefficients,1) > 0 && size(coefficients,3) > 0 && all(isfinite,coefficients) ||
        error("Polynomial source coefficients must be finite and nonempty.")
    output = zeros(Float64,size(coefficients,1),size(coefficients,2),size(coefficients,3),energy_order)
    for entity in axes(output,1), group in axes(output,2), angular in axes(output,3)
        width = edges[group+1]-edges[group]
        a = [power <= size(coefficients,4) ? Float64(coefficients[entity,group,angular,power]) : 0.0 for power in 1:4]
        # Analytic integrals avoid cancellation in orthogonal modes of a lower
        # degree polynomial, which are exactly zero.
        output[entity,group,angular,1] = width*(a[1]+a[3]/3)
        energy_order >= 2 && (output[entity,group,angular,2] = width*(a[2]/sqrt(3.0)+sqrt(3.0)*a[4]/5))
        energy_order >= 3 && (output[entity,group,angular,3] = width*(2sqrt(5.0)*a[3]/15))
        energy_order >= 4 && (output[entity,group,angular,4] = width*(2sqrt(7.0)*a[4]/35))
    end
    all(isfinite,output) || error("Integrated source moments overflowed.")
    return output
end

"""
    Energy_Moment_Volume_Source(source, energy_moments)

Public energy-resolved wrapper for an existing explicitly normalized volume source. Moments have
shape `(voxel, ascending energy group, angular coefficient, energy mode)` in the normalized,
group-integrated convention documented by `source_energy_moments`. Q0 must equal `source.values`
exactly. Higher coefficients may be signed; the represented differential source must remain
nonnegative over each complete energy interval, allowing only a reported scale-dependent
floating-point roundoff bound at a numerically unresolved zero. No coefficient is clipped.
Angular moment data are checked after native
angular reconstruction during projection. One to four modes are supported, with no silent
truncation or zero-padding at transport. Sources remain spatially constant within each voxel.

Use `Fixed_Sources` or `project_volume_source` to project the wrapper; do not alter native arrays.
The selected SN CSD-family solver must use DG energy with exactly the supplied number of modes.
Serialization through the legacy three-dimensional source interchange is unsupported and must
not discard this wrapper's higher modes. Physical source-rate scaling retains the base contract.
"""
struct Energy_Moment_Volume_Source <: Abstract_Radiant_Source
    source::Anisotropic_Volume_Source
    energy_moments::Array{Float64,4}
    function Energy_Moment_Volume_Source(source::Anisotropic_Volume_Source,
        energy_moments::AbstractArray{<:Real}; positivity_tolerance::Real=0.0)
        ndims(energy_moments) == 4 || error("Energy source moments must have four dimensions.")
        size(energy_moments)[1:3] == size(source.values) || error("Energy source moment dimensions must match the base source.")
        1 <= size(energy_moments,4) <= 4 || error("Energy source requires one to four explicit modes.")
        positivity_tolerance == 0 || error("Energy source positivity is strict; a negative allowance is unsupported.")
        moments = Float64.(energy_moments)
        all(isfinite,moments) || error("Energy source moments must be finite.")
        moments[:,:,:,1] == source.values || error("Zeroth energy moments must equal the declared exact group integrals.")
        if source.angular_representation != :moments
            for voxel in axes(moments,1), group in axes(moments,2), angular in axes(moments,3)
                diagnostic = _source_energy_polynomial_diagnostic(view(moments,voxel,group,angular,:))
                diagnostic.minimum >= -diagnostic.roundoff_bound ||
                    error("Represented energy source has a negative differential value.")
            end
        end
        return new(deepcopy(source),moments)
    end
end

get_particle(this::Energy_Moment_Volume_Source) = get_particle(this.source)
get_source_normalization(this::Energy_Moment_Volume_Source) = get_source_normalization(this.source)
get_volume_source_rate(this::Energy_Moment_Volume_Source;physical::Bool=false) =
    get_volume_source_rate(this.source;physical=physical)

"""
    get_volume_source_rate(this::Anisotropic_Volume_Source; physical=false)

Return the angle-integrated source rate per energy group, summed over source voxels. Isotropic
values are interpreted as angle-integrated source densities. Ordinate values are integrated using
the supplied quadrature.

A moment source is reduced only when provenance declares both:

- `zeroth_moment_index`: one-based coefficient index;
- `zeroth_moment_is_angle_integrated = "true"`.

This prevents an implicit factor of `2`, `4π`, or a basis-normalization constant.
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
        if get(this.provenance,"zeroth_moment_is_angle_integrated","false") != "true"
            error("Moment-source scalar-rate extraction requires zeroth_moment_is_angle_integrated=\"true\".")
        end
        if !haskey(this.provenance,"zeroth_moment_index")
            error("Moment-source scalar-rate extraction requires zeroth_moment_index provenance.")
        end
        index = try
            parse(Int,this.provenance["zeroth_moment_index"])
        catch
            error("zeroth_moment_index must be a valid one-based integer.")
        end
        if index < 1 || index > Ncoefficient
            error("zeroth_moment_index lies outside the source coefficient dimension.")
        end
        for ivoxel in 1:Nvoxel, igroup in 1:Ngroup
            rate[igroup] += this.values[ivoxel,igroup,index] * this.voxel_volumes_cm3[ivoxel]
        end
    end

    return physical ? apply_normalization(rate,this.normalization) : rate
end
