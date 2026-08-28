const BOUNDARY_SOURCE_HDF5_SCHEMA = "radiant.boundary_angular_current/v1"
const VOLUME_SOURCE_HDF5_SCHEMA = "radiant.anisotropic_volume_source/v1"
const EXTERNAL_ROW_MAJOR = "external-row-major"
const JULIA_NATIVE_ORDER = "julia-native"

function _source_file_sha256(path::AbstractString)
    isfile(path) || error("Source artifact does not exist: $(path).")
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _hdf5_required(parent,name::AbstractString)
    haskey(parent,name) || error("Required HDF5 object is missing: $(name).")
    return parent[name]
end

function _hdf5_scalar(value,label::AbstractString)
    if value isa Number || value isa AbstractString || value isa Symbol
        return value
    elseif value isa AbstractVector{UInt8}
        return String(value)
    elseif value isa AbstractArray && length(value) == 1
        return _hdf5_scalar(first(value),label)
    end
    error("HDF5 value $(label) must be scalar.")
end

function _hdf5_string(value,label::AbstractString)
    scalar = _hdf5_scalar(value,label)
    return scalar isa AbstractString ? String(scalar) : string(scalar)
end

function _hdf5_float(value,label::AbstractString)
    scalar = _hdf5_scalar(value,label)
    result = try
        Float64(scalar)
    catch
        error("HDF5 value $(label) is not numeric.")
    end
    isfinite(result) || error("HDF5 value $(label) must be finite.")
    return result
end

function _hdf5_vector(value,::Type{T},label::AbstractString) where T
    result = vec(T.(value))
    isempty(result) && error("HDF5 vector $(label) cannot be empty.")
    return result
end

function _validate_storage_order(value::AbstractString)
    value in (EXTERNAL_ROW_MAJOR,JULIA_NATIVE_ORDER) || error(
        "Unsupported HDF5 storage order $(value).",
    )
    return String(value)
end

function _canonical_matrix(
    value,
    rows::Int64,
    columns::Int64,
    label::AbstractString,
    storage_order::AbstractString,
    ::Type{T}=Float64,
) where T
    data = T.(value)
    canonical = if storage_order == EXTERNAL_ROW_MAJOR
        ndims(data) == 2 || error("HDF5 matrix $(label) must have rank two.")
        permutedims(data,(2,1))
    elseif storage_order == JULIA_NATIVE_ORDER
        data
    else
        error("Unsupported HDF5 storage order $(storage_order).")
    end
    size(canonical) == (rows,columns) || error(
        "HDF5 matrix $(label) resolves to shape $(size(canonical)); expected $(rows,columns).",
    )
    return Array{T,2}(canonical)
end

function _canonical_tensor3(
    value,
    first_dimension::Int64,
    second_dimension::Int64,
    third_dimension::Int64,
    label::AbstractString,
    storage_order::AbstractString,
    ::Type{T}=Float64,
) where T
    data = T.(value)
    canonical = if storage_order == EXTERNAL_ROW_MAJOR
        ndims(data) == 3 || error("HDF5 tensor $(label) must have rank three.")
        permutedims(data,(3,2,1))
    elseif storage_order == JULIA_NATIVE_ORDER
        data
    else
        error("Unsupported HDF5 storage order $(storage_order).")
    end
    expected = (first_dimension,second_dimension,third_dimension)
    size(canonical) == expected || error(
        "HDF5 tensor $(label) resolves to shape $(size(canonical)); expected $(expected).",
    )
    return Array{T,3}(canonical)
end

_read_string_dataset(group,name::AbstractString) =
    _hdf5_string(read(_hdf5_required(group,name)),name)
_read_float_dataset(group,name::AbstractString) =
    _hdf5_float(read(_hdf5_required(group,name)),name)

function _read_provenance(file)
    result = Dict{String,String}()
    haskey(file,"provenance") || return result
    group = file["provenance"]
    if haskey(group,"keys") || haskey(group,"values")
        haskey(group,"keys") && haskey(group,"values") || error(
            "HDF5 provenance keys and values must be supplied together.",
        )
        keys_vector = vec(read(group["keys"]))
        values_vector = vec(read(group["values"]))
        length(keys_vector) == length(values_vector) || error(
            "HDF5 provenance keys and values have different lengths.",
        )
        for index in eachindex(keys_vector)
            result[_hdf5_string(keys_vector[index],"provenance key")] =
                _hdf5_string(values_vector[index],"provenance value")
        end
    end
    return result
end

function _write_provenance(file,provenance::AbstractDict)
    normalized = Dict{String,String}()
    for (key,value) in provenance
        normalized[string(key)] = string(value)
    end
    group = HDF5.create_group(file,"provenance")
    ordered_keys = sort!(collect(keys(normalized)))
    group["keys"] = ordered_keys
    group["values"] = String[normalized[key] for key in ordered_keys]
    return group
end

function _normalization_from_hdf5(meta)
    interval = Float64.(vec(read(_hdf5_required(meta,"time_interval_s"))))
    length(interval) == 2 || error("Source time interval must contain start and stop times.")
    return Source_Normalization(
        basis=Symbol(_read_string_dataset(meta,"normalization_basis")),
        source_rate_per_s=_read_float_dataset(meta,"source_rate_per_s"),
        symmetry_factor=_read_float_dataset(meta,"symmetry_factor"),
        time_interval_s=(interval[1],interval[2]),
        time_class=Symbol(_read_string_dataset(meta,"time_class")),
        source_hash=_read_string_dataset(meta,"source_hash"),
    )
end

function _write_normalization(meta,normalization::Source_Normalization)
    meta["normalization_basis"] = String(normalization.basis)
    meta["source_rate_per_s"] = normalization.source_rate_per_s
    meta["symmetry_factor"] = normalization.symmetry_factor
    meta["time_interval_s"] = collect(normalization.time_interval_s)
    meta["time_class"] = String(normalization.time_class)
    meta["source_hash"] = normalization.source_hash
    return meta
end

function _verify_file_hash(path::AbstractString,expected_file_sha256)
    actual = _source_file_sha256(path)
    if !isnothing(expected_file_sha256)
        expected = lowercase(String(expected_file_sha256))
        lowercase(actual) == expected || error(
            "Source artifact SHA-256 mismatch: expected $(expected), calculated $(actual).",
        )
    end
    return actual
end

_write_external_matrix(group,name::AbstractString,value::AbstractMatrix) =
    (group[name] = permutedims(value,(2,1)))
_write_external_tensor3(group,name::AbstractString,value::AbstractArray{<:Any,3}) =
    (group[name] = permutedims(value,(3,2,1)))

"""
    write_boundary_angular_current_hdf5(path, source; overwrite=false)

Write a plain HDF5 boundary source. Multidimensional arrays are reversed before HDF5.jl writes
them so Python/h5py and other row-major readers see canonical `(patch, group, direction)` and
`(patch, component)` shapes. The returned value is the completed file SHA-256.
"""
function write_boundary_angular_current_hdf5(
    path::AbstractString,
    source::Boundary_Angular_Current_Source;
    overwrite::Bool=false,
)
    isfile(path) && !overwrite && error("Refusing to overwrite source artifact $(path).")
    mkpath(dirname(abspath(path)))

    HDF5.h5open(path,"w") do file
        meta = HDF5.create_group(file,"meta")
        patch = HDF5.create_group(file,"patch")
        energy = HDF5.create_group(file,"energy")
        angle = HDF5.create_group(file,"angle")
        source_group = HDF5.create_group(file,"source")

        meta["schema"] = BOUNDARY_SOURCE_HDF5_SCHEMA
        meta["particle"] = get_tag(source.particle)
        meta["source_representation"] = "angular_flux"
        meta["storage_order"] = EXTERNAL_ROW_MAJOR
        _write_normalization(meta,source.normalization)

        patch["id"] = source.patch_ids
        _write_external_matrix(patch,"centroid_cm",source.centroids_cm)
        patch["area_cm2"] = source.areas_cm2
        _write_external_matrix(patch,"normal",source.normals)
        _write_external_matrix(patch,"tangent_1",source.tangent_1)
        _write_external_matrix(patch,"tangent_2",source.tangent_2)
        energy["edges_eV"] = source.energy_edges_eV
        _write_external_matrix(angle,"direction_cosines",source.directions)
        angle["quadrature_weights"] = source.quadrature_weights
        _write_external_tensor3(source_group,"angular_flux",source.angular_flux)
        _write_external_matrix(source_group,"incoming_current",get_incoming_current(source))
        if !isnothing(source.variance)
            _write_external_tensor3(source_group,"variance",source.variance)
        end
        _write_provenance(file,merge(
            source.provenance,
            Dict("writer" => "Radiant.jl","schema" => BOUNDARY_SOURCE_HDF5_SCHEMA),
        ))
    end
    return _source_file_sha256(path)
end

"""
    read_boundary_angular_current_hdf5(path, particle; expected_file_sha256=nothing)

Read an OpenMC/OpenSn/Radiant boundary artifact in the canonical HDF5 schema. The reader supports
`angular_flux` and `directional_current_density` representations, enforces declared storage order,
and validates an optional patch/group incoming-current ledger.
"""
function read_boundary_angular_current_hdf5(
    path::AbstractString,
    particle::Particle;
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
)
    artifact_hash = _verify_file_hash(path,expected_file_sha256)
    return HDF5.h5open(path,"r") do file
        meta = _hdf5_required(file,"meta")
        patch = _hdf5_required(file,"patch")
        energy = _hdf5_required(file,"energy")
        angle = _hdf5_required(file,"angle")
        source_group = _hdf5_required(file,"source")

        schema = _read_string_dataset(meta,"schema")
        schema == BOUNDARY_SOURCE_HDF5_SCHEMA || error(
            "Unsupported boundary-source schema $(schema).",
        )
        storage_order = _validate_storage_order(_read_string_dataset(meta,"storage_order"))
        stored_particle = _read_string_dataset(meta,"particle")
        stored_particle == get_tag(particle) || error(
            "Boundary-source particle $(stored_particle) does not match $(get_tag(particle)).",
        )
        representation = _read_string_dataset(meta,"source_representation")
        normalization = _normalization_from_hdf5(meta)

        patch_ids = _hdf5_vector(read(_hdf5_required(patch,"id")),Int64,"patch/id")
        Npatch = length(patch_ids)
        centroids = _canonical_matrix(
            read(_hdf5_required(patch,"centroid_cm")),Npatch,3,
            "patch/centroid_cm",storage_order,
        )
        areas = _hdf5_vector(read(_hdf5_required(patch,"area_cm2")),Float64,"patch/area_cm2")
        length(areas) == Npatch || error("Patch areas do not match patch identifiers.")
        normals = _canonical_matrix(
            read(_hdf5_required(patch,"normal")),Npatch,3,"patch/normal",storage_order,
        )
        tangent_1 = _canonical_matrix(
            read(_hdf5_required(patch,"tangent_1")),Npatch,3,
            "patch/tangent_1",storage_order,
        )
        tangent_2 = _canonical_matrix(
            read(_hdf5_required(patch,"tangent_2")),Npatch,3,
            "patch/tangent_2",storage_order,
        )

        energy_edges = _hdf5_vector(
            read(_hdf5_required(energy,"edges_eV")),Float64,"energy/edges_eV",
        )
        Ngroup = length(energy_edges)-1
        Ngroup > 0 || error("Boundary source requires at least one energy group.")
        weights = _hdf5_vector(
            read(_hdf5_required(angle,"quadrature_weights")),Float64,
            "angle/quadrature_weights",
        )
        Ndir = length(weights)
        directions = _canonical_matrix(
            read(_hdf5_required(angle,"direction_cosines")),Ndir,3,
            "angle/direction_cosines",storage_order,
        )
        provenance = _read_provenance(file)
        provenance["file_sha256"] = artifact_hash
        provenance["reader"] = "Radiant.jl"

        source = if representation == "angular_flux"
            values = _canonical_tensor3(
                read(_hdf5_required(source_group,"angular_flux")),
                Npatch,Ngroup,Ndir,"source/angular_flux",storage_order,
            )
            variance = haskey(source_group,"variance") ? _canonical_tensor3(
                read(source_group["variance"]),Npatch,Ngroup,Ndir,
                "source/variance",storage_order,
            ) : nothing
            Boundary_Angular_Current_Source(
                particle,patch_ids,centroids,areas,normals,tangent_1,tangent_2,
                energy_edges,directions,weights,values,normalization;
                variance=variance,provenance=provenance,
            )
        elseif representation == "directional_current_density"
            values = _canonical_tensor3(
                read(_hdf5_required(source_group,"directional_current_density")),
                Npatch,Ngroup,Ndir,"source/directional_current_density",storage_order,
            )
            variance = haskey(source_group,"variance") ? _canonical_tensor3(
                read(source_group["variance"]),Npatch,Ngroup,Ndir,
                "source/variance",storage_order,
            ) : nothing
            boundary_source_from_directional_current(
                particle,patch_ids,centroids,areas,normals,tangent_1,tangent_2,
                energy_edges,directions,weights,values,normalization;
                variance=variance,provenance=provenance,
            )
        else
            error("Unsupported boundary source representation $(representation).")
        end

        if haskey(source_group,"incoming_current")
            reference = _canonical_matrix(
                read(source_group["incoming_current"]),Npatch,Ngroup,
                "source/incoming_current",storage_order,
            )
            assert_current_closure(source,reference)
        end
        source
    end
end

"""Write an anisotropic volume source to the plain HDF5 interchange schema."""
function write_anisotropic_volume_source_hdf5(
    path::AbstractString,
    source::Anisotropic_Volume_Source;
    overwrite::Bool=false,
)
    isfile(path) && !overwrite && error("Refusing to overwrite source artifact $(path).")
    mkpath(dirname(abspath(path)))

    HDF5.h5open(path,"w") do file
        meta = HDF5.create_group(file,"meta")
        voxel = HDF5.create_group(file,"voxel")
        energy = HDF5.create_group(file,"energy")
        angle = HDF5.create_group(file,"angle")
        source_group = HDF5.create_group(file,"source")

        meta["schema"] = VOLUME_SOURCE_HDF5_SCHEMA
        meta["particle"] = get_tag(source.particle)
        meta["angular_representation"] = String(source.angular_representation)
        meta["storage_order"] = EXTERNAL_ROW_MAJOR
        meta["parent_reaction"] = isnothing(source.parent_reaction) ? "none" : source.parent_reaction
        _write_normalization(meta,source.normalization)

        voxel["id"] = source.voxel_ids
        voxel["volume_cm3"] = source.voxel_volumes_cm3
        energy["edges_eV"] = source.energy_edges_eV
        if source.angular_representation == :ordinates
            _write_external_matrix(angle,"direction_cosines",source.directions)
            angle["quadrature_weights"] = source.quadrature_weights
        end
        _write_external_tensor3(source_group,"values",source.values)
        source_group["integrated_rate"] = get_volume_source_rate(source)
        if !isnothing(source.variance)
            _write_external_tensor3(source_group,"variance",source.variance)
        end
        _write_provenance(file,merge(
            source.provenance,
            Dict("writer" => "Radiant.jl","schema" => VOLUME_SOURCE_HDF5_SCHEMA),
        ))
    end
    return _source_file_sha256(path)
end

"""Read a hash-bound anisotropic volume source from plain HDF5."""
function read_anisotropic_volume_source_hdf5(
    path::AbstractString,
    particle::Particle;
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
)
    artifact_hash = _verify_file_hash(path,expected_file_sha256)
    return HDF5.h5open(path,"r") do file
        meta = _hdf5_required(file,"meta")
        voxel = _hdf5_required(file,"voxel")
        energy = _hdf5_required(file,"energy")
        angle = _hdf5_required(file,"angle")
        source_group = _hdf5_required(file,"source")

        schema = _read_string_dataset(meta,"schema")
        schema == VOLUME_SOURCE_HDF5_SCHEMA || error(
            "Unsupported volume-source schema $(schema).",
        )
        storage_order = _validate_storage_order(_read_string_dataset(meta,"storage_order"))
        stored_particle = _read_string_dataset(meta,"particle")
        stored_particle == get_tag(particle) || error(
            "Volume-source particle $(stored_particle) does not match $(get_tag(particle)).",
        )
        representation = Symbol(_read_string_dataset(meta,"angular_representation"))
        normalization = _normalization_from_hdf5(meta)
        reaction_value = _read_string_dataset(meta,"parent_reaction")
        parent_reaction = reaction_value == "none" ? nothing : reaction_value

        voxel_ids = _hdf5_vector(read(_hdf5_required(voxel,"id")),Int64,"voxel/id")
        Nvoxel = length(voxel_ids)
        volumes = _hdf5_vector(
            read(_hdf5_required(voxel,"volume_cm3")),Float64,"voxel/volume_cm3",
        )
        length(volumes) == Nvoxel || error("Volume-source volumes do not match voxel identifiers.")
        energy_edges = _hdf5_vector(
            read(_hdf5_required(energy,"edges_eV")),Float64,"energy/edges_eV",
        )
        Ngroup = length(energy_edges)-1
        Ngroup > 0 || error("Volume source requires at least one energy group.")

        directions = zeros(Float64,0,3)
        weights = Float64[]
        Ncoefficient = 1
        if representation == :ordinates
            weights = _hdf5_vector(
                read(_hdf5_required(angle,"quadrature_weights")),Float64,
                "angle/quadrature_weights",
            )
            Ncoefficient = length(weights)
            directions = _canonical_matrix(
                read(_hdf5_required(angle,"direction_cosines")),Ncoefficient,3,
                "angle/direction_cosines",storage_order,
            )
        elseif representation == :moments
            raw_shape = size(read(_hdf5_required(source_group,"values")))
            length(raw_shape) == 3 || error("Moment source values must be rank three.")
            Ncoefficient = storage_order == EXTERNAL_ROW_MAJOR ? raw_shape[1] : raw_shape[3]
        elseif representation != :isotropic
            error("Unsupported volume-source angular representation $(representation).")
        end

        raw_values = read(_hdf5_required(source_group,"values"))
        values = _canonical_tensor3(
            raw_values,Nvoxel,Ngroup,Ncoefficient,"source/values",storage_order,
        )
        variance = haskey(source_group,"variance") ? _canonical_tensor3(
            read(source_group["variance"]),Nvoxel,Ngroup,Ncoefficient,
            "source/variance",storage_order,
        ) : nothing
        provenance = _read_provenance(file)
        provenance["file_sha256"] = artifact_hash
        provenance["reader"] = "Radiant.jl"

        source = Anisotropic_Volume_Source(
            particle,voxel_ids,volumes,energy_edges,representation,values,normalization;
            directions=directions,quadrature_weights=weights,variance=variance,
            parent_reaction=parent_reaction,provenance=provenance,
        )
        if haskey(source_group,"integrated_rate")
            reference = _hdf5_vector(
                read(source_group["integrated_rate"]),Float64,"source/integrated_rate",
            )
            calculated = get_volume_source_rate(source)
            length(reference) == length(calculated) || error(
                "Stored and calculated volume-source rate vectors have different lengths.",
            )
            for index in eachindex(reference)
                isapprox(calculated[index],reference[index];rtol=1.0e-10,atol=1.0e-12) || error(
                    "Volume-source integrated-rate closure failed in group $(index).",
                )
            end
        end
        source
    end
end
