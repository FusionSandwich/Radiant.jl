const FACETED_GEOMETRY_HDF5_SCHEMA = "radiant.faceted_surface_mesh/v1"

"""
    Faceted_Surface_Mesh

Triangle surface representation for DAGMC/ParaStell/CAD-derived local domains. Coordinates are
centimetres. Each triangle carries a surface ID, adjacent volume ID, and material tag. Direct
unstructured transport is not implied: the mesh can be converted to a local tape atlas or a
conservatively sampled Cartesian material map for the existing Radiant sweep.
"""
struct Faceted_Surface_Mesh
    vertices_cm::Matrix{Float64}
    triangles::Matrix{Int64}
    surface_ids::Vector{Int64}
    volume_ids::Vector{Int64}
    material_tags::Vector{String}
    geometry_hash::String
    provenance::Dict{String,String}

    function Faceted_Surface_Mesh(
        vertices_cm::AbstractMatrix{<:Real},
        triangles::AbstractMatrix{<:Integer};
        surface_ids::AbstractVector{<:Integer}=collect(1:size(triangles,1)),
        volume_ids::AbstractVector{<:Integer}=ones(Int64,size(triangles,1)),
        material_tags::AbstractVector{<:AbstractString}=fill("unbound",size(triangles,1)),
        geometry_hash::AbstractString="unbound",
        provenance::AbstractDict=Dict{String,String}(),
        minimum_area_cm2::Real=1.0e-20,
    )
        vertices = Float64.(vertices_cm)
        connectivity = Int64.(triangles)
        size(vertices,2) == 3 || error("Faceted vertices must have shape (N,3).")
        size(connectivity,2) == 3 || error("Triangle connectivity must have shape (M,3).")
        size(vertices,1) >= 3 || error("Faceted geometry requires at least three vertices.")
        size(connectivity,1) >= 1 || error("Faceted geometry requires at least one triangle.")
        all(isfinite,vertices) || error("Faceted coordinates must be finite.")
        minimum(connectivity) >= 1 && maximum(connectivity) <= size(vertices,1) || error(
            "Triangle connectivity contains an out-of-range vertex index.",
        )
        count = size(connectivity,1)
        surfaces = Int64.(surface_ids)
        volumes = Int64.(volume_ids)
        materials = String[string(value) for value in material_tags]
        length(surfaces) == count && length(volumes) == count && length(materials) == count ||
            error("Facet tag arrays must match the triangle count.")
        all(value -> value >= 1,surfaces) || error("Surface identifiers must be positive.")
        all(value -> value >= 1,volumes) || error("Volume identifiers must be positive.")
        all(value -> !isempty(value),materials) || error("Facet material tags cannot be empty.")
        isempty(geometry_hash) && error("Faceted geometry hash cannot be empty.")
        threshold = Float64(minimum_area_cm2)
        isfinite(threshold) && threshold > 0.0 || error(
            "Minimum facet area must be finite and positive.",
        )
        for triangle_index in 1:count
            indices = vec(connectivity[triangle_index,:])
            length(unique(indices)) == 3 || error(
                "Triangle $(triangle_index) repeats a vertex index.",
            )
            p1 = view(vertices,indices[1],:)
            p2 = view(vertices,indices[2],:)
            p3 = view(vertices,indices[3],:)
            area = 0.5*norm(cross(p2-p1,p3-p1))
            isfinite(area) && area > threshold || error(
                "Triangle $(triangle_index) is degenerate or below the minimum area.",
            )
        end
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            vertices,connectivity,surfaces,volumes,materials,String(geometry_hash),
            provenance_string,
        )
    end
end

function facet_geometry(mesh::Faceted_Surface_Mesh)
    count = size(mesh.triangles,1)
    centroids = zeros(Float64,count,3)
    normals = zeros(Float64,count,3)
    areas = zeros(Float64,count)
    for index in 1:count
        ids = vec(mesh.triangles[index,:])
        p1 = view(mesh.vertices_cm,ids[1],:)
        p2 = view(mesh.vertices_cm,ids[2],:)
        p3 = view(mesh.vertices_cm,ids[3],:)
        cross_product = cross(p2-p1,p3-p1)
        norm_value = norm(cross_product)
        areas[index] = 0.5*norm_value
        normals[index,:] .= cross_product/norm_value
        centroids[index,:] .= (p1+p2+p3)/3.0
    end
    return centroids,normals,areas
end

function _facet_edges(triangle::AbstractVector{<:Integer})
    return (
        Tuple(sort(Int64[triangle[1],triangle[2]])),
        Tuple(sort(Int64[triangle[2],triangle[3]])),
        Tuple(sort(Int64[triangle[3],triangle[1]])),
    )
end

function faceted_topology_report(mesh::Faceted_Surface_Mesh)
    edge_incidence = Dict{Tuple{Int64,Int64},Vector{Int64}}()
    duplicate_lookup = Dict{NTuple{3,Int64},Vector{Int64}}()
    for triangle_index in axes(mesh.triangles,1)
        triangle = vec(mesh.triangles[triangle_index,:])
        for edge in _facet_edges(triangle)
            push!(get!(edge_incidence,edge,Int64[]),triangle_index)
        end
        canonical = Tuple(sort(Int64.(triangle)))
        push!(get!(duplicate_lookup,canonical,Int64[]),triangle_index)
    end
    boundary_edges = [edge for (edge,owners) in edge_incidence if length(owners) == 1]
    nonmanifold_edges = [edge for (edge,owners) in edge_incidence if length(owners) > 2]
    duplicate_groups = [owners for owners in values(duplicate_lookup) if length(owners) > 1]
    centroids,normals,areas = facet_geometry(mesh)
    signed_volume_cm3 = 0.0
    for triangle_index in axes(mesh.triangles,1)
        ids = vec(mesh.triangles[triangle_index,:])
        p1 = view(mesh.vertices_cm,ids[1],:)
        p2 = view(mesh.vertices_cm,ids[2],:)
        p3 = view(mesh.vertices_cm,ids[3],:)
        signed_volume_cm3 += dot(p1,cross(p2,p3))/6.0
    end
    return Dict{String,Any}(
        "schema" => "radiant.faceted_topology_report/v1",
        "geometry_hash" => mesh.geometry_hash,
        "vertex_count" => size(mesh.vertices_cm,1),
        "triangle_count" => size(mesh.triangles,1),
        "surface_count" => length(unique(mesh.surface_ids)),
        "volume_count" => length(unique(mesh.volume_ids)),
        "boundary_edge_count" => length(boundary_edges),
        "nonmanifold_edge_count" => length(nonmanifold_edges),
        "duplicate_triangle_group_count" => length(duplicate_groups),
        "watertight" => isempty(boundary_edges) && isempty(nonmanifold_edges),
        "manifold" => isempty(nonmanifold_edges),
        "minimum_area_cm2" => minimum(areas),
        "maximum_area_cm2" => maximum(areas),
        "signed_enclosed_volume_cm3" => signed_volume_cm3,
        "orientation_is_positive" => signed_volume_cm3 > 0.0,
    )
end

function assert_transport_ready_facets(
    mesh::Faceted_Surface_Mesh;
    require_positive_orientation::Bool=true,
)
    report = faceted_topology_report(mesh)
    report["watertight"] || error(
        "Faceted geometry is not watertight; open or nonmanifold edges remain.",
    )
    report["duplicate_triangle_group_count"] == 0 || error(
        "Faceted geometry contains duplicate triangles.",
    )
    if require_positive_orientation
        report["orientation_is_positive"] || error(
            "Faceted geometry has non-positive signed volume; orient triangles consistently.",
        )
    end
    return report
end

function write_faceted_geometry_hdf5(
    path::AbstractString,
    mesh::Faceted_Surface_Mesh;
    overwrite::Bool=false,
)
    isfile(path) && !overwrite && error("Refusing to overwrite faceted geometry $(path).")
    mkpath(dirname(abspath(path)))
    HDF5.h5open(path,"w") do file
        meta = HDF5.create_group(file,"meta")
        geometry = HDF5.create_group(file,"geometry")
        meta["schema"] = FACETED_GEOMETRY_HDF5_SCHEMA
        meta["storage_order"] = EXTERNAL_ROW_MAJOR
        meta["length_unit"] = "cm"
        meta["geometry_hash"] = mesh.geometry_hash
        _write_external_matrix(geometry,"vertices_cm",mesh.vertices_cm)
        _write_external_matrix(geometry,"triangles",mesh.triangles)
        geometry["surface_ids"] = mesh.surface_ids
        geometry["volume_ids"] = mesh.volume_ids
        geometry["material_tags"] = mesh.material_tags
        _write_provenance(file,merge(
            mesh.provenance,
            Dict("writer" => "Radiant.jl","schema" => FACETED_GEOMETRY_HDF5_SCHEMA),
        ))
    end
    return _source_file_sha256(path)
end

function read_faceted_geometry_hdf5(
    path::AbstractString;
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
)
    artifact_hash = _verify_file_hash(path,expected_file_sha256)
    return HDF5.h5open(path,"r") do file
        meta = _hdf5_required(file,"meta")
        geometry = _hdf5_required(file,"geometry")
        _read_string_dataset(meta,"schema") == FACETED_GEOMETRY_HDF5_SCHEMA || error(
            "Unsupported faceted-geometry HDF5 schema.",
        )
        _read_string_dataset(meta,"length_unit") == "cm" || error(
            "Faceted-geometry coordinates must be in centimetres.",
        )
        storage_order = _validate_storage_order(_read_string_dataset(meta,"storage_order"))
        surface_ids = _hdf5_vector(
            read(_hdf5_required(geometry,"surface_ids")),Int64,"surface_ids",
        )
        triangle_count = length(surface_ids)
        raw_vertices = read(_hdf5_required(geometry,"vertices_cm"))
        vertex_count = storage_order == EXTERNAL_ROW_MAJOR ? size(raw_vertices,2) :
                                                             size(raw_vertices,1)
        vertices = _canonical_matrix(
            raw_vertices,vertex_count,3,"geometry/vertices_cm",storage_order,
        )
        triangles = _canonical_matrix(
            read(_hdf5_required(geometry,"triangles")),triangle_count,3,
            "geometry/triangles",storage_order,Int64,
        )
        volume_ids = _hdf5_vector(
            read(_hdf5_required(geometry,"volume_ids")),Int64,"volume_ids",
        )
        material_raw = vec(read(_hdf5_required(geometry,"material_tags")))
        material_tags = String[_hdf5_string(value,"material tag") for value in material_raw]
        provenance = _read_provenance(file)
        provenance["file_sha256"] = artifact_hash
        provenance["reader"] = "Radiant.jl"
        Faceted_Surface_Mesh(
            vertices,triangles;
            surface_ids=surface_ids,volume_ids=volume_ids,material_tags=material_tags,
            geometry_hash=_read_string_dataset(meta,"geometry_hash"),
            provenance=provenance,
        )
    end
end

function read_obj_facets(
    path::AbstractString;
    length_multiplier_to_cm::Real=1.0,
    geometry_hash::AbstractString="unbound",
)
    isfile(path) || error("OBJ file does not exist: $(path).")
    multiplier = Float64(length_multiplier_to_cm)
    isfinite(multiplier) && multiplier > 0.0 || error("OBJ length multiplier must be positive.")
    vertices = Vector{NTuple{3,Float64}}()
    triangles = Vector{NTuple{3,Int64}}()
    surface_ids = Int64[]
    current_surface = 1
    surface_names = Dict{Int64,String}(1 => "default")
    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line,"#") && continue
        fields = split(line)
        isempty(fields) && continue
        if fields[1] == "v"
            length(fields) >= 4 || error("Malformed OBJ vertex line.")
            push!(vertices,(
                multiplier*parse(Float64,fields[2]),
                multiplier*parse(Float64,fields[3]),
                multiplier*parse(Float64,fields[4]),
            ))
        elseif fields[1] in ("g","o")
            current_surface += 1
            surface_names[current_surface] = length(fields) >= 2 ? fields[2] :
                                                               "surface-$(current_surface)"
        elseif fields[1] == "f"
            length(fields) >= 4 || error("OBJ face requires at least three vertices.")
            indices = Int64[]
            for token in fields[2:end]
                vertex_token = first(split(token,"/"))
                value = parse(Int64,vertex_token)
                value < 0 && (value = length(vertices)+value+1)
                push!(indices,value)
            end
            for local_index in 2:length(indices)-1
                push!(triangles,(indices[1],indices[local_index],indices[local_index+1]))
                push!(surface_ids,current_surface)
            end
        end
    end
    isempty(vertices) && error("OBJ file contains no vertices.")
    isempty(triangles) && error("OBJ file contains no faces.")
    vertex_matrix = reduce(vcat,[reshape(collect(value),1,3) for value in vertices])
    triangle_matrix = reduce(vcat,[reshape(collect(value),1,3) for value in triangles])
    materials = [get(surface_names,id,"surface-$(id)") for id in surface_ids]
    return Faceted_Surface_Mesh(
        vertex_matrix,triangle_matrix;
        surface_ids=surface_ids,volume_ids=ones(Int64,length(surface_ids)),
        material_tags=materials,geometry_hash=geometry_hash,
        provenance=Dict(
            "source_format" => "OBJ",
            "source_path" => abspath(path),
            "length_multiplier_to_cm" => string(multiplier),
        ),
    )
end

function read_ascii_stl_facets(
    path::AbstractString;
    length_multiplier_to_cm::Real=1.0,
    geometry_hash::AbstractString="unbound",
    material_tag::AbstractString="unbound",
)
    isfile(path) || error("STL file does not exist: $(path).")
    multiplier = Float64(length_multiplier_to_cm)
    isfinite(multiplier) && multiplier > 0.0 || error("STL length multiplier must be positive.")
    vertices = Vector{NTuple{3,Float64}}()
    for raw_line in eachline(path)
        fields = split(strip(raw_line))
        length(fields) == 4 && lowercase(fields[1]) == "vertex" || continue
        push!(vertices,(
            multiplier*parse(Float64,fields[2]),
            multiplier*parse(Float64,fields[3]),
            multiplier*parse(Float64,fields[4]),
        ))
    end
    length(vertices) >= 3 && length(vertices)%3 == 0 || error(
        "ASCII STL must contain a multiple of three vertex records.",
    )
    unique_vertices = Vector{NTuple{3,Float64}}()
    vertex_lookup = Dict{NTuple{3,Float64},Int64}()
    connectivity = Matrix{Int64}(undef,length(vertices)÷3,3)
    for (index,vertex) in enumerate(vertices)
        if !haskey(vertex_lookup,vertex)
            push!(unique_vertices,vertex)
            vertex_lookup[vertex] = length(unique_vertices)
        end
        connectivity[(index-1)÷3+1,(index-1)%3+1] = vertex_lookup[vertex]
    end
    vertex_matrix = reduce(vcat,[reshape(collect(value),1,3) for value in unique_vertices])
    triangle_count = size(connectivity,1)
    return Faceted_Surface_Mesh(
        vertex_matrix,connectivity;
        surface_ids=ones(Int64,triangle_count),volume_ids=ones(Int64,triangle_count),
        material_tags=fill(String(material_tag),triangle_count),geometry_hash=geometry_hash,
        provenance=Dict(
            "source_format" => "ASCII-STL",
            "source_path" => abspath(path),
            "length_multiplier_to_cm" => string(multiplier),
        ),
    )
end

function _ray_triangle_intersection(
    origin::Vector{Float64},
    direction::Vector{Float64},
    p1,
    p2,
    p3;
    tolerance::Real=1.0e-12,
)
    edge1 = p2-p1
    edge2 = p3-p1
    h = cross(direction,edge2)
    determinant = dot(edge1,h)
    abs(determinant) <= tolerance && return nothing
    inverse = 1.0/determinant
    s = origin-p1
    u = inverse*dot(s,h)
    (u < -tolerance || u > 1.0+tolerance) && return nothing
    q = cross(s,edge1)
    v = inverse*dot(direction,q)
    (v < -tolerance || u+v > 1.0+tolerance) && return nothing
    distance = inverse*dot(edge2,q)
    distance > tolerance || return nothing
    return distance
end

function point_in_faceted_volume(
    mesh::Faceted_Surface_Mesh,
    point_cm::AbstractVector{<:Real};
    volume_id::Integer=first(sort(unique(mesh.volume_ids))),
)
    length(point_cm) == 3 || error("Point-in-mesh query requires a three-vector.")
    point = Float64.(point_cm)
    all(isfinite,point) || error("Point-in-mesh query must be finite.")
    selected = findall(==(Int64(volume_id)),mesh.volume_ids)
    isempty(selected) && error("Requested volume ID is absent from the faceted mesh.")
    submesh = Faceted_Surface_Mesh(
        mesh.vertices_cm,mesh.triangles[selected,:];
        surface_ids=mesh.surface_ids[selected],volume_ids=mesh.volume_ids[selected],
        material_tags=mesh.material_tags[selected],geometry_hash=mesh.geometry_hash,
        provenance=mesh.provenance,
    )
    assert_transport_ready_facets(submesh;require_positive_orientation=false)
    ray_directions = (
        normalize([1.0,0.17320508075688773,0.071]),
        normalize([0.137,1.0,0.22360679774997896]),
        normalize([0.193,0.113,1.0]),
    )
    votes = 0
    for direction in ray_directions
        distances = Float64[]
        for triangle_index in axes(submesh.triangles,1)
            ids = vec(submesh.triangles[triangle_index,:])
            distance = _ray_triangle_intersection(
                point,direction,view(submesh.vertices_cm,ids[1],:),
                view(submesh.vertices_cm,ids[2],:),view(submesh.vertices_cm,ids[3],:),
            )
            isnothing(distance) || push!(distances,distance)
        end
        sort!(distances)
        unique_distances = Float64[]
        for distance in distances
            if isempty(unique_distances) || abs(distance-unique_distances[end]) > 1.0e-9
                push!(unique_distances,distance)
            end
        end
        isodd(length(unique_distances)) && (votes += 1)
    end
    return votes >= 2
end

struct Faceted_Voxelization
    x_boundaries_cm::Vector{Float64}
    y_boundaries_cm::Vector{Float64}
    z_boundaries_cm::Vector{Float64}
    volume_ids::Vector{Int64}
    volume_fractions::Array{Float64,4}
    dominant_volume_id::Array{Int64,3}
    unresolved_fraction::Array{Float64,3}
    samples_per_axis::Int64
    geometry_hash::String
    provenance::Dict{String,String}
end

function voxelize_faceted_mesh(
    mesh::Faceted_Surface_Mesh,
    x_boundaries_cm::AbstractVector{<:Real},
    y_boundaries_cm::AbstractVector{<:Real},
    z_boundaries_cm::AbstractVector{<:Real};
    samples_per_axis::Integer=3,
)
    assert_transport_ready_facets(mesh;require_positive_orientation=false)
    boundaries = (
        Float64.(x_boundaries_cm),Float64.(y_boundaries_cm),Float64.(z_boundaries_cm),
    )
    for axis in boundaries
        length(axis) >= 2 && all(isfinite,axis) && all(diff(axis) .> 0.0) || error(
            "Voxelization boundaries must be finite and strictly increasing.",
        )
    end
    samples_per_axis >= 1 || error("Voxelization sampling order must be positive.")
    volume_ids = sort(unique(mesh.volume_ids))
    counts = (length(boundaries[1])-1,length(boundaries[2])-1,length(boundaries[3])-1)
    fractions = zeros(Float64,length(volume_ids),counts[1],counts[2],counts[3])
    dominant = zeros(Int64,counts...)
    unresolved = zeros(Float64,counts...)
    offsets = [(index-0.5)/samples_per_axis for index in 1:samples_per_axis]
    total_samples = samples_per_axis^3
    for ix in 1:counts[1], iy in 1:counts[2], iz in 1:counts[3]
        assigned = 0
        for fx in offsets, fy in offsets, fz in offsets
            point = [
                boundaries[1][ix]+fx*(boundaries[1][ix+1]-boundaries[1][ix]),
                boundaries[2][iy]+fy*(boundaries[2][iy+1]-boundaries[2][iy]),
                boundaries[3][iz]+fz*(boundaries[3][iz+1]-boundaries[3][iz]),
            ]
            inside = Int64[]
            for (volume_index,volume_id) in enumerate(volume_ids)
                if point_in_faceted_volume(mesh,point;volume_id=volume_id)
                    push!(inside,volume_index)
                end
            end
            length(inside) <= 1 || error(
                "Faceted volumes overlap at a voxelization sample point.",
            )
            if length(inside) == 1
                fractions[inside[1],ix,iy,iz] += 1.0/total_samples
                assigned += 1
            end
        end
        unresolved[ix,iy,iz] = 1.0-assigned/total_samples
        local_fractions = view(fractions,:,ix,iy,iz)
        if maximum(local_fractions) > 0.0
            dominant[ix,iy,iz] = volume_ids[argmax(local_fractions)]
        end
    end
    return Faceted_Voxelization(
        boundaries[1],boundaries[2],boundaries[3],volume_ids,fractions,dominant,
        unresolved,Int64(samples_per_axis),mesh.geometry_hash,Dict(
            "schema" => "radiant.faceted_voxelization/v1",
            "sampling" => "deterministic-subcell-centres",
            "direct_unstructured_transport" => "false",
        ),
    )
end

struct Faceted_Tape_Patch
    facet_index::Int64
    surface_id::Int64
    centroid_cm::NTuple{3,Float64}
    area_cm2::Float64
    local_to_global::Matrix{Float64}
    curvature_proxy_cm_inv::Float64
    material_tag::String
end

function build_faceted_tape_patches(
    mesh::Faceted_Surface_Mesh;
    selected_surface_ids::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    reference_tangent::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    centroids,normals,areas = facet_geometry(mesh)
    selected = isnothing(selected_surface_ids) ? collect(axes(mesh.triangles,1)) :
        findall(id -> id in Int64.(selected_surface_ids),mesh.surface_ids)
    isempty(selected) && error("No facets match the requested tape surfaces.")
    edge_to_facets = Dict{Tuple{Int64,Int64},Vector{Int64}}()
    for triangle_index in selected
        for edge in _facet_edges(vec(mesh.triangles[triangle_index,:]))
            push!(get!(edge_to_facets,edge,Int64[]),triangle_index)
        end
    end
    patches = Faceted_Tape_Patch[]
    for triangle_index in selected
        normal = vec(normals[triangle_index,:])
        triangle = vec(mesh.triangles[triangle_index,:])
        tangent = if isnothing(reference_tangent)
            edges = [
                vec(mesh.vertices_cm[triangle[2],:]-mesh.vertices_cm[triangle[1],:]),
                vec(mesh.vertices_cm[triangle[3],:]-mesh.vertices_cm[triangle[2],:]),
                vec(mesh.vertices_cm[triangle[1],:]-mesh.vertices_cm[triangle[3],:]),
            ]
            normalize(edges[argmax(norm.(edges))])
        else
            candidate = Float64.(reference_tangent)
            candidate = candidate-dot(candidate,normal)*normal
            norm(candidate) > 1.0e-12 || error(
                "Reference tangent is parallel to a selected facet normal.",
            )
            normalize(candidate)
        end
        width_axis = normalize(cross(normal,tangent))
        tangent = normalize(cross(width_axis,normal))
        frame = hcat(tangent,width_axis,normal)
        curvature = 0.0
        for edge in _facet_edges(triangle)
            owners = get(edge_to_facets,edge,Int64[])
            for neighbour in owners
                neighbour == triangle_index && continue
                distance = norm(vec(centroids[neighbour,:]-centroids[triangle_index,:]))
                distance <= 1.0e-14 && continue
                angle = acos(clamp(dot(normal,vec(normals[neighbour,:])),-1.0,1.0))
                curvature = max(curvature,angle/distance)
            end
        end
        centroid = vec(centroids[triangle_index,:])
        push!(patches,Faceted_Tape_Patch(
            triangle_index,mesh.surface_ids[triangle_index],
            (centroid[1],centroid[2],centroid[3]),areas[triangle_index],frame,curvature,
            mesh.material_tags[triangle_index],
        ))
    end
    return patches
end

function map_boundary_patches_to_facets(
    source::Boundary_Angular_Current_Source,
    mesh::Faceted_Surface_Mesh;
    maximum_centroid_distance_cm::Real,
    minimum_normal_dot::Real=0.5,
)
    centroids,normals,areas = facet_geometry(mesh)
    maximum_distance = Float64(maximum_centroid_distance_cm)
    isfinite(maximum_distance) && maximum_distance > 0.0 || error(
        "Maximum source-to-facet distance must be positive.",
    )
    mappings = Vector{Dict{String,Any}}()
    used_facets = Set{Int64}()
    for patch_index in eachindex(source.patch_ids)
        distances = [
            norm(vec(centroids[facet_index,:])-vec(source.centroids_cm[patch_index,:]))
            for facet_index in axes(mesh.triangles,1)
        ]
        facet_index = argmin(distances)
        distances[facet_index] <= maximum_distance || error(
            "Boundary source patch $(source.patch_ids[patch_index]) has no nearby CAD facet.",
        )
        normal_dot = dot(
            vec(source.normals[patch_index,:]),vec(normals[facet_index,:]),
        )
        abs(normal_dot) >= minimum_normal_dot || error(
            "Boundary source and CAD facet normals are inconsistent.",
        )
        facet_index in used_facets && error(
            "Multiple source patches map to one facet; refine or provide explicit IDs.",
        )
        push!(used_facets,facet_index)
        push!(mappings,Dict{String,Any}(
            "source_patch_id" => source.patch_ids[patch_index],
            "facet_index" => facet_index,
            "surface_id" => mesh.surface_ids[facet_index],
            "centroid_distance_cm" => distances[facet_index],
            "normal_dot" => normal_dot,
            "facet_area_cm2" => areas[facet_index],
            "source_area_cm2" => source.areas_cm2[patch_index],
            "area_ratio" => source.areas_cm2[patch_index]/areas[facet_index],
        ))
    end
    return Dict{String,Any}(
        "schema" => "radiant.faceted_source_mapping/v1",
        "source_hash" => source.normalization.source_hash,
        "geometry_hash" => mesh.geometry_hash,
        "mapping_count" => length(mappings),
        "mappings" => mappings,
        "current_conservation_required_after_projection" => true,
    )
end
