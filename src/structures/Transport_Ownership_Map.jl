"""
    Transport_Ownership_Record

One particle/domain/source/response ownership declaration used to prevent double counting across
OpenSn, Radiant, OpenMC, Geant4, and downstream damage tools.
"""
struct Transport_Ownership_Record
    domain::String
    particle::String
    source_class::String
    response::String
    time_class::Symbol
    owner::String
    comparison_only::Bool
    artifact_hash::String

    function Transport_Ownership_Record(;
        domain::AbstractString,
        particle::AbstractString,
        source_class::AbstractString,
        response::AbstractString,
        time_class::Symbol = :prompt,
        owner::AbstractString,
        comparison_only::Bool = false,
        artifact_hash::AbstractString = "unbound",
    )
        fields = (domain,particle,source_class,response,owner,artifact_hash)
        if any(isempty,fields)
            error("Ownership-record strings cannot be empty.")
        end
        if time_class ∉ (:prompt,:delayed)
            error("Ownership time class must be :prompt or :delayed.")
        end
        return new(
            String(domain),
            String(particle),
            String(source_class),
            String(response),
            time_class,
            String(owner),
            comparison_only,
            String(artifact_hash),
        )
    end
end

function _ownership_key(record::Transport_Ownership_Record)
    return (
        record.domain,
        record.particle,
        record.source_class,
        record.response,
        record.time_class,
    )
end

"""
    Transport_Ownership_Map

Versioned no-double-counting map.

- `mode = :production`: every key must have exactly one non-comparison owner.
- `mode = :comparison`: every key must have at least two distinct comparison-only owners, and no
  result from that key may be added to another.
"""
struct Transport_Ownership_Map
    schema::String
    mode::Symbol
    records::Vector{Transport_Ownership_Record}

    function Transport_Ownership_Map(
        records::AbstractVector{Transport_Ownership_Record};
        mode::Symbol = :production,
        schema::AbstractString = "parastell_damage.transport_ownership_map/v1",
    )
        if mode ∉ (:production,:comparison)
            error("Ownership-map mode must be :production or :comparison.")
        end
        if isempty(records)
            error("Ownership map cannot be empty.")
        end
        this = new(String(schema),mode,collect(records))
        validate_ownership(this)
        return this
    end
end

"""
    validate_ownership(this::Transport_Ownership_Map)

Validate uniqueness in production mode and non-additive multiplicity in comparison mode.
"""
function validate_ownership(this::Transport_Ownership_Map)
    grouped = Dict{Tuple{String,String,String,String,Symbol},Vector{Transport_Ownership_Record}}()
    for record in this.records
        push!(get!(grouped,_ownership_key(record),Transport_Ownership_Record[]),record)
    end

    for (key,records) in grouped
        production_records = [record for record in records if !record.comparison_only]
        comparison_records = [record for record in records if record.comparison_only]

        if this.mode == :production
            if length(production_records) != 1
                error("Production ownership key $(key) must have exactly one production owner; found $(length(production_records)).")
            end
            production_owner_name = production_records[1].owner
            if any(record -> record.owner == production_owner_name, comparison_records)
                error("A production owner must not also be duplicated as a comparison-only record for key $(key).")
            end
        else
            if !isempty(production_records)
                error("Comparison ownership maps may contain only comparison-only records.")
            end
            owners = unique([record.owner for record in comparison_records])
            if length(owners) < 2
                error("Comparison ownership key $(key) must contain at least two distinct solvers.")
            end
        end
    end
    return true
end

"""
    get_production_owner(this::Transport_Ownership_Map; ...)

Return the unique production owner for a key. This function is not available for comparison maps.
"""
function get_production_owner(
    this::Transport_Ownership_Map;
    domain::AbstractString,
    particle::AbstractString,
    source_class::AbstractString,
    response::AbstractString,
    time_class::Symbol = :prompt,
)
    if this.mode != :production
        error("Comparison maps do not define an additive production owner.")
    end
    key = (String(domain),String(particle),String(source_class),String(response),time_class)
    matches = [record for record in this.records if _ownership_key(record) == key && !record.comparison_only]
    if length(matches) != 1
        error("Ownership key $(key) does not have exactly one production owner.")
    end
    return matches[1].owner
end
