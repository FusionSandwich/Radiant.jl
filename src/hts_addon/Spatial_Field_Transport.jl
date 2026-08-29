"""
    electromagnetic_field_on_geometry(field, geometry; field_hash=nothing)

Sample an add-on spatial field at every Radiant voxel centre and return a generic
`Electromagnetic_Field` carrying the full cell-centred map. This activates the cell-varying
operator in the ordinary `compute_flux` sweep; it does not decompose the solve into independent
constant-field patches.
"""
function electromagnetic_field_on_geometry(
    field::Abstract_Spatial_Magnetic_Field,
    geometry::Geometry;
    field_hash::Union{Nothing,AbstractString}=nothing,
)
    values = field_on_geometry(field,geometry)
    hash_value = if isnothing(field_hash)
        field isa Cartesian_Magnetic_Field_Map ? field.field_hash : "constant-spatial-sampling"
    else
        String(field_hash)
    end
    output = Electromagnetic_Field()
    set_spatial_magnetic_field(
        output,values;
        field_hash=hash_value,
        provenance=Dict(
            "schema" => "radiant.hts.spatial_field_transport/v1",
            "sampling" => "voxel-centre",
            "geometry_type" => get_type(geometry),
            "geometry_dimension" => string(get_dimension(geometry)),
        ),
    )
    return output
end

"""Return a fail-closed receipt for one cell-varying magnetic-field transport input."""
function spatial_field_transport_receipt(
    electromagnetic_field::Electromagnetic_Field,
    geometry::Geometry,
)
    has_spatial_magnetic_field(electromagnetic_field) || error(
        "A spatial-field transport receipt requires a cell-centred field map.",
    )
    values = get_spatial_magnetic_field(electromagnetic_field)
    counts = get_number_of_voxels(geometry)
    expected = (3,counts[1],counts[2],counts[3])
    size(values) == expected || error(
        "Spatial field and geometry dimensions do not agree.",
    )
    magnitudes = Float64[]
    for ix in axes(values,2), iy in axes(values,3), iz in axes(values,4)
        push!(magnitudes,norm(view(values,:,ix,iy,iz)))
    end
    return Dict{String,Any}(
        "schema" => "radiant.hts.spatial_field_transport_receipt/v1",
        "field_hash" => get_spatial_field_hash(electromagnetic_field),
        "field_mode" => string(get_field_mode(electromagnetic_field)),
        "shape" => collect(size(values)),
        "minimum_magnitude_T" => minimum(magnitudes),
        "maximum_magnitude_T" => maximum(magnitudes),
        "operator_application" => "cell-local-within-one-sn-sweep",
        "electric_field_supported" => false,
        "software_qualification" => true,
        "physical_qualification" => false,
        "physical_qualification_reason" =>
            "Requires comparison against a charged-particle reference in a nonuniform field.",
        "provenance" => get_spatial_field_provenance(electromagnetic_field),
    )
end
