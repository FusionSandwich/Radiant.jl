abstract type Abstract_Spatial_Magnetic_Field end

"""Uniform magnetic field implementing the same spatial interface as mapped fields."""
struct Constant_Magnetic_Field <: Abstract_Spatial_Magnetic_Field
    magnetic_field_T::NTuple{3,Float64}
    provenance::Dict{String,String}

    function Constant_Magnetic_Field(
        magnetic_field_T::AbstractVector{<:Real};
        provenance::AbstractDict=Dict{String,String}(),
    )
        length(magnetic_field_T) == 3 || error("Magnetic field must be a three-vector.")
        field = Float64.(magnetic_field_T)
        all(isfinite,field) || error("Magnetic field components must be finite.")
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new((field[1],field[2],field[3]),provenance_string)
    end
end

"""
    Cartesian_Magnetic_Field_Map

Global Cartesian magnetic-field map. `values_T` has shape `(3,Nx,Ny,Nz)`. Coordinates are in cm,
field components are in tesla, and interpolation is either `:nearest` or `:trilinear`.

The map is an HTS add-on object. Radiant's current sweep still consumes one constant field per
local calculation; `electromagnetic_field_at` and atlas helpers provide conservative piecewise-
constant field realization until a spatially varying sweep operator is independently qualified.
"""
struct Cartesian_Magnetic_Field_Map <: Abstract_Spatial_Magnetic_Field
    x_cm::Vector{Float64}
    y_cm::Vector{Float64}
    z_cm::Vector{Float64}
    values_T::Array{Float64,4}
    interpolation::Symbol
    field_hash::String
    provenance::Dict{String,String}

    function Cartesian_Magnetic_Field_Map(
        x_cm::AbstractVector{<:Real},
        y_cm::AbstractVector{<:Real},
        z_cm::AbstractVector{<:Real},
        values_T::AbstractArray{<:Real,4};
        interpolation::Symbol=:trilinear,
        field_hash::AbstractString="unbound",
        provenance::AbstractDict=Dict{String,String}(),
    )
        interpolation in (:nearest,:trilinear) || error(
            "Magnetic-field interpolation must be :nearest or :trilinear.",
        )
        axes = (Float64.(x_cm),Float64.(y_cm),Float64.(z_cm))
        for (name,axis) in zip(("x","y","z"),axes)
            isempty(axis) && error("Magnetic-field $(name)-axis cannot be empty.")
            all(isfinite,axis) || error("Magnetic-field coordinates must be finite.")
            length(axis) == 1 || all(diff(axis) .> 0.0) || error(
                "Magnetic-field coordinates must be strictly increasing.",
            )
        end
        values = Float64.(values_T)
        size(values) == (3,length(axes[1]),length(axes[2]),length(axes[3])) || error(
            "Magnetic-field values must have shape (3,Nx,Ny,Nz).",
        )
        all(isfinite,values) || error("Magnetic-field map values must be finite.")
        isempty(field_hash) && error("Field-map hash cannot be empty; use unbound for fixtures.")
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            axes[1],axes[2],axes[3],values,interpolation,String(field_hash),
            provenance_string,
        )
    end
end

field_at(this::Constant_Magnetic_Field,position_cm::AbstractVector{<:Real}) = begin
    length(position_cm) == 3 || error("Field position must be a three-vector.")
    all(isfinite,position_cm) || error("Field position must be finite.")
    collect(this.magnetic_field_T)
end

function _field_axis_bracket(axis::Vector{Float64},coordinate::Float64)
    coordinate < axis[1] && error("Field query lies below the tabulated domain.")
    coordinate > axis[end] && error("Field query lies above the tabulated domain.")
    length(axis) == 1 && return 1,1,0.0
    upper = searchsortedfirst(axis,coordinate)
    if upper == 1
        return 1,1,0.0
    elseif upper > length(axis)
        return length(axis),length(axis),0.0
    elseif axis[upper] == coordinate
        return upper,upper,0.0
    end
    lower = upper-1
    fraction = (coordinate-axis[lower])/(axis[upper]-axis[lower])
    return lower,upper,fraction
end

function field_at(
    this::Cartesian_Magnetic_Field_Map,
    position_cm::AbstractVector{<:Real},
)
    length(position_cm) == 3 || error("Field position must be a three-vector.")
    position = Float64.(position_cm)
    all(isfinite,position) || error("Field position must be finite.")

    if this.interpolation == :nearest
        indices = ntuple(dimension -> begin
            axis = (this.x_cm,this.y_cm,this.z_cm)[dimension]
            argmin(abs.(axis .- position[dimension]))
        end,3)
        return collect(view(this.values_T,:,indices[1],indices[2],indices[3]))
    end

    x0,x1,fx = _field_axis_bracket(this.x_cm,position[1])
    y0,y1,fy = _field_axis_bracket(this.y_cm,position[2])
    z0,z1,fz = _field_axis_bracket(this.z_cm,position[3])
    result = zeros(Float64,3)
    for (ix,wx) in ((x0,1.0-fx),(x1,fx)),
        (iy,wy) in ((y0,1.0-fy),(y1,fy)),
        (iz,wz) in ((z0,1.0-fz),(z1,fz))
        weight = wx*wy*wz
        weight == 0.0 && continue
        result .+= weight .* view(this.values_T,:,ix,iy,iz)
    end
    return result
end

"""Rotate a global magnetic-field vector into a right-handed local frame matrix."""
function field_in_local_frame(
    field::Abstract_Spatial_Magnetic_Field,
    position_cm::AbstractVector{<:Real},
    local_to_global::AbstractMatrix{<:Real};
    orthogonality_tolerance::Real=1.0e-9,
)
    frame = Float64.(local_to_global)
    size(frame) == (3,3) || error("Local frame must be a 3×3 matrix.")
    isapprox(frame'*frame,Matrix{Float64}(I,3,3);rtol=orthogonality_tolerance,
             atol=orthogonality_tolerance) || error("Local frame must be orthonormal.")
    det(frame) > 0.0 || error("Local frame must be right-handed.")
    return frame' * field_at(field,position_cm)
end

"""Create the constant Radiant field used for one atlas patch or local microdomain."""
function electromagnetic_field_at(
    field::Abstract_Spatial_Magnetic_Field,
    position_cm::AbstractVector{<:Real};
    local_to_global::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
)
    magnetic_field = isnothing(local_to_global) ? field_at(field,position_cm) :
        field_in_local_frame(field,position_cm,local_to_global)
    output = Electromagnetic_Field()
    set_magnetic_field(output,magnetic_field)
    return output
end

"""Evaluate the mapped magnetic field at all Cartesian Radiant voxel centres."""
function field_on_geometry(
    field::Abstract_Spatial_Magnetic_Field,
    geometry::Geometry,
)
    get_type(geometry) == "cartesian" || error(
        "Field-map sampling currently requires Cartesian Radiant geometry.",
    )
    dimensions = get_dimension(geometry)
    x = geometry.voxels_position["x"]
    y = dimensions ≥ 2 ? geometry.voxels_position["y"] : [0.0]
    z = dimensions ≥ 3 ? geometry.voxels_position["z"] : [0.0]
    values = zeros(Float64,3,length(x),length(y),length(z))
    for ix in eachindex(x), iy in eachindex(y), iz in eachindex(z)
        values[:,ix,iy,iz] .= field_at(field,[x[ix],y[iy],z[iz]])
    end
    return values
end

"""
    field_variation_bound(field, points_cm)

Maximum Euclidean field difference from the field at the point-cloud centroid. This is a screening
metric for deciding whether one constant field is adequate for a local patch. It is not a transport
error estimate and must be combined with response convergence.
"""
function field_variation_bound(
    field::Abstract_Spatial_Magnetic_Field,
    points_cm::AbstractMatrix{<:Real},
)
    size(points_cm,2) == 3 || error("Field variation points must have shape (N,3).")
    size(points_cm,1) ≥ 1 || error("At least one field variation point is required.")
    points = Float64.(points_cm)
    all(isfinite,points) || error("Field variation points must be finite.")
    centroid = vec(sum(points,dims=1)./size(points,1))
    reference = field_at(field,centroid)
    maximum_variation = 0.0
    for index in axes(points,1)
        maximum_variation = max(
            maximum_variation,
            norm(field_at(field,view(points,index,:))-reference),
        )
    end
    return maximum_variation
end

function field_map_receipt(field::Cartesian_Magnetic_Field_Map)
    magnitudes = [norm(view(field.values_T,:,ix,iy,iz))
                  for ix in axes(field.values_T,2),
                      iy in axes(field.values_T,3),
                      iz in axes(field.values_T,4)]
    return Dict{String,Any}(
        "schema" => "radiant.hts.spatial_magnetic_field/v1",
        "field_hash" => field.field_hash,
        "interpolation" => string(field.interpolation),
        "shape" => collect(size(field.values_T)),
        "minimum_magnitude_T" => minimum(magnitudes),
        "maximum_magnitude_T" => maximum(magnitudes),
        "provenance" => copy(field.provenance),
    )
end
