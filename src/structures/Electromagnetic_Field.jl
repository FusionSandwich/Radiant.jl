"""
    Electromagnetic_Field

Generic Radiant electromagnetic-field input. A field can be either uniform or represented by a
cell-centred Cartesian magnetic-field array with shape `(3,Nx,Ny,Nz)`. Spatial fields are a small
generic core hook used by the temporary HTS add-on; geometry construction, interpolation, CAD
mapping, and field provenance remain outside the transport core.

Electric-field transport is still unavailable because it requires an energy-redistribution
operator. Any nonzero electric field is rejected.
"""
mutable struct Electromagnetic_Field
    electric_field              ::Vector{Float64}
    magnetic_field              ::Vector{Float64}
    spatial_magnetic_field      ::Union{Nothing,Array{Float64,4}}
    spatial_field_hash          ::String
    spatial_field_provenance    ::Dict{String,String}

    function Electromagnetic_Field()
        this = new()
        this.electric_field = [0.0,0.0,0.0]
        this.magnetic_field = [0.0,0.0,0.0]
        this.spatial_magnetic_field = nothing
        this.spatial_field_hash = "none"
        this.spatial_field_provenance = Dict{String,String}()
        return this
    end
end

function _validate_field_vector(field::AbstractVector{<:Real},label::AbstractString)
    length(field) == 3 || error("$(label) must be a three-vector.")
    value = Float64.(field)
    all(isfinite,value) || error("$(label) components must be finite.")
    return value
end

"""Set one uniform magnetic field [T] and clear any previously assigned spatial map."""
function set_magnetic_field(
    this::Electromagnetic_Field,
    magnetic_field::AbstractVector{<:Real},
)
    this.magnetic_field = _validate_field_vector(magnetic_field,"Magnetic field")
    this.spatial_magnetic_field = nothing
    this.spatial_field_hash = "none"
    empty!(this.spatial_field_provenance)
    return this
end

"""
    set_spatial_magnetic_field(this, values_T; field_hash, provenance=Dict())

Assign a cell-centred magnetic field with shape `(3,Nx,Ny,Nz)`. The dimensions are checked again
against the transport geometry when a solve begins. No interpolation occurs in the sweep: every
cell receives the field already sampled at its declared representative position.
"""
function set_spatial_magnetic_field(
    this::Electromagnetic_Field,
    values_T::AbstractArray{<:Real,4};
    field_hash::AbstractString,
    provenance::AbstractDict=Dict{String,String}(),
)
    size(values_T,1) == 3 || error(
        "Spatial magnetic field must have shape (3,Nx,Ny,Nz).",
    )
    all(size(values_T,dimension) >= 1 for dimension in 2:4) || error(
        "Spatial magnetic field dimensions must be nonempty.",
    )
    values = Float64.(values_T)
    all(isfinite,values) || error("Spatial magnetic-field values must be finite.")
    isempty(field_hash) && error("Spatial magnetic field requires a nonempty artifact hash.")
    provenance_string = Dict{String,String}()
    for (key,value) in provenance
        provenance_string[string(key)] = string(value)
    end
    this.spatial_magnetic_field = values
    this.spatial_field_hash = String(field_hash)
    this.spatial_field_provenance = provenance_string
    this.magnetic_field .= 0.0
    return this
end

function clear_spatial_magnetic_field(this::Electromagnetic_Field)
    this.spatial_magnetic_field = nothing
    this.spatial_field_hash = "none"
    empty!(this.spatial_field_provenance)
    return this
end

function set_electric_field(
    this::Electromagnetic_Field,
    electric_field::AbstractVector{<:Real},
)
    value = _validate_field_vector(electric_field,"Electric field")
    any(value .!= 0.0) && error(
        "Electric-field transport is not implemented; only magnetic fields are supported.",
    )
    this.electric_field = value
    return this
end

get_magnetic_field(this::Electromagnetic_Field) = copy(this.magnetic_field)
get_electric_field(this::Electromagnetic_Field) = copy(this.electric_field)
has_spatial_magnetic_field(this::Electromagnetic_Field) =
    !isnothing(this.spatial_magnetic_field)

function get_spatial_magnetic_field(this::Electromagnetic_Field)
    isnothing(this.spatial_magnetic_field) && error(
        "No spatial magnetic field is assigned.",
    )
    return this.spatial_magnetic_field
end

get_spatial_field_hash(this::Electromagnetic_Field) = this.spatial_field_hash
get_spatial_field_provenance(this::Electromagnetic_Field) =
    copy(this.spatial_field_provenance)

function get_field_mode(this::Electromagnetic_Field)
    has_spatial_magnetic_field(this) && return :spatial_cell_centered
    any(this.magnetic_field .!= 0.0) && return :uniform
    return :none
end

"""
Return whether a field is active together with the uniform electric and magnetic vectors. For a
spatial magnetic field, the returned uniform magnetic vector is zero; callers must inspect
`has_spatial_magnetic_field` and `get_spatial_magnetic_field`.
"""
function get_electromagnetic_field(this::Electromagnetic_Field)
    is_active = any(this.electric_field .!= 0.0) ||
                any(this.magnetic_field .!= 0.0) ||
                has_spatial_magnetic_field(this)
    return is_active,copy(this.electric_field),copy(this.magnetic_field)
end
