"""
    Source_Normalization

Explicit normalization and provenance attached to a fixed source.

All source arrays in the HTS coupling interfaces are interpreted according to `basis`:

- `:per_history`: values are per sampled source history.
- `:per_source_particle`: values are per source particle.
- `:per_second`: values are already physical rates.

For per-history or per-source-particle data, `source_rate_per_s` converts the stored values to
physical rates. `symmetry_factor` is applied separately so sector and full-device normalization
cannot be hidden in the source tensor.
"""
struct Source_Normalization
    basis::Symbol
    source_rate_per_s::Float64
    symmetry_factor::Float64
    time_interval_s::Tuple{Float64,Float64}
    time_class::Symbol
    source_hash::String
    provenance::Dict{String,String}

    function Source_Normalization(;
        basis::Symbol = :per_history,
        source_rate_per_s::Real = 1.0,
        symmetry_factor::Real = 1.0,
        time_interval_s::Tuple{<:Real,<:Real} = (0.0,0.0),
        time_class::Symbol = :prompt,
        source_hash::AbstractString = "unbound",
        provenance::AbstractDict = Dict{String,String}(),
    )
        if basis ∉ (:per_history,:per_source_particle,:per_second)
            error("Unknown source normalization basis: $(basis).")
        end
        if time_class ∉ (:prompt,:delayed)
            error("Source time class must be :prompt or :delayed.")
        end

        rate = Float64(source_rate_per_s)
        symmetry = Float64(symmetry_factor)
        t0 = Float64(time_interval_s[1])
        t1 = Float64(time_interval_s[2])

        if !isfinite(rate) || rate ≤ 0.0
            error("Source rate must be finite and greater than zero.")
        end
        if !isfinite(symmetry) || symmetry ≤ 0.0
            error("Symmetry factor must be finite and greater than zero.")
        end
        if !isfinite(t0) || !isfinite(t1) || t1 < t0
            error("Source time interval must be finite and ordered as (start, stop).")
        end
        if isempty(source_hash)
            error("Source hash cannot be empty; use \"unbound\" for an analytic fixture.")
        end

        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end

        return new(
            basis,
            rate,
            symmetry,
            (t0,t1),
            time_class,
            String(source_hash),
            provenance_string,
        )
    end
end

"""
    get_physical_scale(this::Source_Normalization)

Return the multiplicative factor used to convert stored source-normalized values to full-device
physical rates. Values already expressed per second are not multiplied by `source_rate_per_s`.
"""
function get_physical_scale(this::Source_Normalization)
    if this.basis == :per_second
        return this.symmetry_factor
    end
    return this.source_rate_per_s * this.symmetry_factor
end

"""
    get_duration(this::Source_Normalization)

Return the duration of the source interval in seconds.
"""
function get_duration(this::Source_Normalization)
    return this.time_interval_s[2] - this.time_interval_s[1]
end

"""
    apply_normalization(value, this::Source_Normalization)

Convert a stored source-normalized value to a full-device physical rate.
"""
function apply_normalization(value, this::Source_Normalization)
    return value .* get_physical_scale(this)
end

# `get_normalization_factor` is retained for legacy Source/Fixed_Sources objects and is exported
# here so tests and downstream adapters can inspect the legacy divisor explicitly.
export get_normalization_factor
