"""Precomputed metrics for a logically Cartesian, physically mapped hexahedral grid."""
struct Mapped_Structured_Grid
    schema::String
    coordinate_system::Symbol
    logical_edges::NTuple{3,Vector{Float64}}
    nodes_cm::Array{Float64,4}
    cell_centers_cm::Array{Float64,4}
    jacobians::Array{Float64,5}
    inverse_jacobians::Array{Float64,5}
    jacobian_determinants::Array{Float64,3}
    cell_volumes_cm3::Array{Float64,3}
    face_area_vectors_cm2::Array{Float64,6}
    face_areas_cm2::Array{Float64,5}
    face_normals::Array{Float64,6}
    geometry_hash::String
end

function _mapped_edges(values,label)
    edges = Float64.(values)
    length(edges) >= 2 || error("$(label) requires at least two logical edges.")
    all(isfinite,edges) || error("$(label) logical edges must be finite.")
    all(diff(edges) .> 0.0) || error("$(label) logical edges must be strictly increasing.")
    return edges
end

function _mapped_node(nodes,i,j,k)
    return collect(view(nodes,:,i,j,k))
end

function _mapped_face_geometry(nodes,i,j,k,axis,side)
    offset = side == 1 ? 0 : 1
    if axis == 1
        points = (
            _mapped_node(nodes,i+offset,j,k),_mapped_node(nodes,i+offset,j+1,k),
            _mapped_node(nodes,i+offset,j,k+1),_mapped_node(nodes,i+offset,j+1,k+1),
        )
        tangent_a = 0.5*((points[2]-points[1])+(points[4]-points[3]))
        tangent_b = 0.5*((points[3]-points[1])+(points[4]-points[2]))
        positive = cross(tangent_a,tangent_b)
    elseif axis == 2
        points = (
            _mapped_node(nodes,i,j+offset,k),_mapped_node(nodes,i,j+offset,k+1),
            _mapped_node(nodes,i+1,j+offset,k),_mapped_node(nodes,i+1,j+offset,k+1),
        )
        tangent_a = 0.5*((points[2]-points[1])+(points[4]-points[3]))
        tangent_b = 0.5*((points[3]-points[1])+(points[4]-points[2]))
        positive = cross(tangent_a,tangent_b)
    else
        points = (
            _mapped_node(nodes,i,j,k+offset),_mapped_node(nodes,i+1,j,k+offset),
            _mapped_node(nodes,i,j+1,k+offset),_mapped_node(nodes,i+1,j+1,k+offset),
        )
        tangent_a = 0.5*((points[2]-points[1])+(points[4]-points[3]))
        tangent_b = 0.5*((points[3]-points[1])+(points[4]-points[2]))
        positive = cross(tangent_a,tangent_b)
    end
    area_vector = side == 1 ? -positive : positive
    return sum(points)/4.0,area_vector
end

"""
    build_mapped_structured_grid(edges, mapping; coordinate_system, geometry_hash)

Sample a native structured mapping at logical nodes and cache cell-centred Jacobians,
contravariant inverse metrics, outward integrated face metrics, areas, normals, and volumes.
The represented cell is the straight-sided hexahedron defined by the sampled nodes.
"""
function build_mapped_structured_grid(
    logical_edges::NTuple{3,Any},
    mapping::Function;
    coordinate_system::Symbol,
    geometry_hash::AbstractString,
)
    isempty(geometry_hash) && error("Mapped geometry hash cannot be empty.")
    edges = ntuple(index -> _mapped_edges(logical_edges[index],"axis $(index)"),3)
    counts = ntuple(index -> length(edges[index])-1,3)
    nodes = zeros(Float64,3,counts[1]+1,counts[2]+1,counts[3]+1)
    for i in 1:counts[1]+1, j in 1:counts[2]+1, k in 1:counts[3]+1
        point = Float64.(mapping(edges[1][i],edges[2][j],edges[3][k]))
        length(point) == 3 && all(isfinite,point) || error(
            "Mapped coordinates must be finite three-vectors.",
        )
        nodes[:,i,j,k] .= point
    end

    centers = zeros(Float64,3,counts...)
    jacobians = zeros(Float64,3,3,counts...)
    inverse_jacobians = similar(jacobians)
    determinants = zeros(Float64,counts...)
    volumes = zeros(Float64,counts...)
    face_vectors = zeros(Float64,3,2,3,counts...)
    face_areas = zeros(Float64,2,3,counts...)
    face_normals = zeros(Float64,3,2,3,counts...)

    for i in 1:counts[1], j in 1:counts[2], k in 1:counts[3]
        vertices = [
            _mapped_node(nodes,i+di,j+dj,k+dk)
            for di in 0:1, dj in 0:1, dk in 0:1
        ]
        center = sum(vertices)/8.0
        centers[:,i,j,k] .= center
        delta = (edges[1][i+1]-edges[1][i],edges[2][j+1]-edges[2][j],edges[3][k+1]-edges[3][k])
        derivative_1 = sum(vertices[2,jj,kk]-vertices[1,jj,kk] for jj in 1:2, kk in 1:2)/(4.0*delta[1])
        derivative_2 = sum(vertices[ii,2,kk]-vertices[ii,1,kk] for ii in 1:2, kk in 1:2)/(4.0*delta[2])
        derivative_3 = sum(vertices[ii,jj,2]-vertices[ii,jj,1] for ii in 1:2, jj in 1:2)/(4.0*delta[3])
        jacobian = hcat(derivative_1,derivative_2,derivative_3)
        determinant = det(jacobian)
        isfinite(determinant) && determinant > 0.0 || error(
            "Mapped cell ($(i),$(j),$(k)) has a non-positive Jacobian.",
        )
        jacobians[:,:,i,j,k] .= jacobian
        inverse_jacobians[:,:,i,j,k] .= inv(jacobian)
        determinants[i,j,k] = determinant
        volume = 0.0
        for axis in 1:3, side in 1:2
            face_center,area_vector = _mapped_face_geometry(nodes,i,j,k,axis,side)
            area = norm(area_vector)
            isfinite(area) && area > 0.0 || error("Mapped cell has a degenerate face.")
            dot(area_vector,face_center-center) > 0.0 || error(
                "Mapped face orientation is not outward.",
            )
            face_vectors[:,side,axis,i,j,k] .= area_vector
            face_areas[side,axis,i,j,k] = area
            face_normals[:,side,axis,i,j,k] .= area_vector/area
            volume += dot(face_center,area_vector)/3.0
        end
        isfinite(volume) && volume > 0.0 || error("Mapped cell volume must be positive.")
        volumes[i,j,k] = volume
    end
    grid = Mapped_Structured_Grid(
        "radiant.mapped_structured_grid/v1",coordinate_system,edges,nodes,centers,
        jacobians,inverse_jacobians,determinants,volumes,face_vectors,face_areas,
        face_normals,String(geometry_hash),
    )
    assert_mapped_grid_invariants(grid)
    return grid
end

function mapped_grid_invariants(grid::Mapped_Structured_Grid; tolerance::Real=1.0e-9)
    tol = Float64(tolerance)
    isfinite(tol) && tol > 0.0 || error("Mapped-grid tolerance must be positive.")
    counts = size(grid.cell_volumes_cm3)
    maximum_metric_error = 0.0
    maximum_normal_error = 0.0
    maximum_face_closure = 0.0
    for i in 1:counts[1], j in 1:counts[2], k in 1:counts[3]
        product = grid.jacobians[:,:,i,j,k]*grid.inverse_jacobians[:,:,i,j,k]
        maximum_metric_error = max(maximum_metric_error,maximum(abs.(product-Matrix{Float64}(I,3,3))))
        closure = zeros(Float64,3)
        total_area = 0.0
        for axis in 1:3, side in 1:2
            normal = view(grid.face_normals,:,side,axis,i,j,k)
            maximum_normal_error = max(maximum_normal_error,abs(norm(normal)-1.0))
            closure .+= view(grid.face_area_vectors_cm2,:,side,axis,i,j,k)
            total_area += grid.face_areas_cm2[side,axis,i,j,k]
        end
        maximum_face_closure = max(maximum_face_closure,norm(closure)/total_area)
    end
    passed = minimum(grid.jacobian_determinants) > 0.0 &&
        minimum(grid.cell_volumes_cm3) > 0.0 && maximum_metric_error <= tol &&
        maximum_normal_error <= tol && maximum_face_closure <= 10.0*tol
    return Dict{String,Any}(
        "schema" => "radiant.mapped_grid_invariants/v1",
        "coordinate_system" => string(grid.coordinate_system),
        "cell_count" => prod(counts),
        "minimum_jacobian_determinant" => minimum(grid.jacobian_determinants),
        "minimum_cell_volume_cm3" => minimum(grid.cell_volumes_cm3),
        "maximum_metric_inverse_error" => maximum_metric_error,
        "maximum_face_normal_error" => maximum_normal_error,
        "maximum_relative_face_closure" => maximum_face_closure,
        "passed" => passed,
    )
end

function assert_mapped_grid_invariants(grid::Mapped_Structured_Grid; tolerance::Real=1.0e-9)
    result = mapped_grid_invariants(grid;tolerance=tolerance)
    result["passed"] || error("Mapped structured-grid invariants do not close.")
    return true
end

function mapped_cartesian_grid(x,y,z; geometry_hash::AbstractString="analytic-cartesian")
    return build_mapped_structured_grid((x,y,z),(a,b,c) -> (a,b,c);
        coordinate_system=:cartesian,geometry_hash=geometry_hash)
end

function mapped_cylindrical_grid(r,theta,z; geometry_hash::AbstractString="analytic-cylindrical")
    return build_mapped_structured_grid((r,theta,z),(radius,angle,axial) ->
        (radius*cos(angle),radius*sin(angle),axial);
        coordinate_system=:cylindrical,geometry_hash=geometry_hash)
end

function mapped_toroidal_grid(r,theta,phi; major_radius_cm::Real,
    geometry_hash::AbstractString="analytic-toroidal")
    major = Float64(major_radius_cm)
    major > maximum(Float64.(r)) || error("Toroidal major radius must exceed minor radius.")
    return build_mapped_structured_grid((r,theta,phi),(radius,poloidal,toroidal) -> begin
        radial = major+radius*cos(poloidal)
        # `(r,poloidal,toroidal)` is kept right-handed for the transport metric cache.
        (radial*cos(toroidal),radial*sin(toroidal),-radius*sin(poloidal))
    end;coordinate_system=:toroidal,geometry_hash=geometry_hash)
end

function mapped_helical_grid(u,v,angle; radius_cm::Real,pitch_per_turn_cm::Real,
    coordinate_system::Symbol=:helical,geometry_hash::AbstractString="analytic-helical")
    radius = Float64(radius_cm)
    pitch = Float64(pitch_per_turn_cm)
    radius > 0.0 && isfinite(pitch) || error("Helical dimensions are invalid.")
    return build_mapped_structured_grid((u,v,angle),(a,b,t) -> begin
        radial = [cos(t),sin(t),0.0]
        tangent = [-radius*sin(t),radius*cos(t),pitch/(2.0*pi)]
        tangent ./= norm(tangent)
        binormal = cross(tangent,radial)
        center = [radius*cos(t),radius*sin(t),pitch*t/(2.0*pi)]
        center+a*radial+b*binormal
    end;coordinate_system=coordinate_system,geometry_hash=geometry_hash)
end

function mapped_frenet_grid(u,v,s,centerline::Function,normal::Function,binormal::Function;
    coordinate_system::Symbol=:varying_curvature_torsion,
    geometry_hash::AbstractString="analytic-varying-curvature-torsion")
    return build_mapped_structured_grid((u,v,s),(a,b,arc) -> begin
        center = Float64.(centerline(arc))
        n = _atlas_normalize(normal(arc),"mapped Frenet normal")
        bin = _atlas_normalize(binormal(arc),"mapped Frenet binormal")
        length(center) == 3 && all(isfinite,center) || error("Mapped centerline is invalid.")
        abs(dot(n,bin)) <= 1.0e-9 || error("Mapped Frenet frame is not orthogonal.")
        center+a*n+b*bin
    end;coordinate_system=coordinate_system,geometry_hash=geometry_hash)
end

function mapped_fibre_grid(r,theta,z; geometry_hash::AbstractString="analytic-fibre")
    return build_mapped_structured_grid((r,theta,z),(radius,angle,axial) ->
        (radius*cos(angle),radius*sin(angle),axial);
        coordinate_system=:fibre,geometry_hash=geometry_hash)
end

"""Cache immutable atlas frames for repeated point/vector remapping."""
struct Atlas_Mapping_Cache
    atlas_hash::String
    centroids_cm::Matrix{Float64}
    frames::Array{Float64,3}
    inverse_frames::Array{Float64,3}
end

function build_atlas_mapping_cache(atlas::Piecewise_Flat_Tape_Atlas)
    count = length(atlas.patches)
    centroids = zeros(Float64,count,3)
    frames = zeros(Float64,3,3,count)
    inverse_frames = similar(frames)
    for patch in atlas.patches
        centroids[patch.patch_id,:] .= collect(patch.centroid_cm)
        frame = patch_frame(patch)
        frames[:,:,patch.patch_id] .= frame
        inverse_frames[:,:,patch.patch_id] .= frame'
    end
    return Atlas_Mapping_Cache(atlas.geometry_hash,centroids,frames,inverse_frames)
end

function atlas_global_to_local(cache::Atlas_Mapping_Cache,patch_id::Integer,point)
    1 <= patch_id <= size(cache.frames,3) || error("Atlas cache patch is out of range.")
    return cache.inverse_frames[:,:,patch_id]*(Float64.(point)-cache.centroids_cm[patch_id,:])
end

function atlas_local_to_global(cache::Atlas_Mapping_Cache,patch_id::Integer,point)
    1 <= patch_id <= size(cache.frames,3) || error("Atlas cache patch is out of range.")
    return cache.centroids_cm[patch_id,:]+cache.frames[:,:,patch_id]*Float64.(point)
end

struct Conservative_Interface_Transfer
    source_edges::Vector{Float64}
    target_edges::Vector{Float64}
    matrix::Matrix{Float64}
end

function conservative_interface_transfer(source_edges,target_edges; tolerance::Real=1.0e-12)
    source = _mapped_edges(source_edges,"source interface")
    target = _mapped_edges(target_edges,"target interface")
    scale = max(abs(source[1]),abs(source[end]),abs(target[1]),abs(target[end]),1.0)
    isapprox(source[1],target[1];rtol=0.0,atol=tolerance*scale) &&
        isapprox(source[end],target[end];rtol=0.0,atol=tolerance*scale) || error(
        "Conservative interfaces must cover the same physical interval.",
    )
    matrix = zeros(Float64,length(target)-1,length(source)-1)
    for j in axes(matrix,1), i in axes(matrix,2)
        overlap = max(0.0,min(target[j+1],source[i+1])-max(target[j],source[i]))
        matrix[j,i] = overlap/(source[i+1]-source[i])
    end
    maximum(abs.(vec(sum(matrix,dims=1)) .- 1.0)) <= 100.0*eps(Float64) || error(
        "Conservative interface weights do not close.",
    )
    return Conservative_Interface_Transfer(source,target,matrix)
end

function apply_interface_transfer(transfer::Conservative_Interface_Transfer,source_integrals)
    values = Float64.(source_integrals)
    length(values) == size(transfer.matrix,2) || error("Source interface value count mismatch.")
    all(isfinite,values) || error("Source interface values must be finite.")
    target = transfer.matrix*values
    isapprox(sum(target),sum(values);rtol=1.0e-12,atol=1.0e-12*max(sum(abs,values),1.0)) || error(
        "Conservative interface transfer failed integral closure.",
    )
    return target
end

struct Multiscale_Local_Layer_Fixture
    scale_ratio::Float64
    layer_ids::Vector{String}
    layer_thickness_cm::Vector{Float64}
    block_local_origins_cm::Vector{Float64}
    block_length_scales_cm::Vector{Float64}
    nondimensional_edges::Vector{Float64}
    length_scale_cm::Float64
    grid::Mapped_Structured_Grid
    jacobian_condition_number::Float64
end

function multiscale_local_layer_fixture(scale_ratio::Real)
    ratio = Float64(scale_ratio)
    ratio in (1.0e3,1.0e5,1.0e7) || error("Multiscale fixture ratio must be 1e3, 1e5, or 1e7.")
    length_scale = 1.0e-5 # 0.1 micrometre in cm
    layer_ids = [
        "film-0.1um","film-0.25um","film-0.5um","film-1.0um",
        "buffer-0.3um","metallic-cap-10um","substrate-50um",
        "cable-scale","macro-block",
    ]
    # At 1e7 this ordered analytic hierarchy reaches from a 0.1 um film through
    # a 1 mm cable to a 1 m macro block. Smaller fixtures keep the same ordered
    # roles while bounding the cable so that max(thickness)/min(thickness) is
    # exactly the requested 1e3 or 1e5 ratio.
    macro_thickness = ratio*length_scale
    cable_thickness = min(0.1,0.5*macro_thickness) # 1 mm when the scale permits it
    thicknesses = length_scale.*[1.0,2.5,5.0,10.0,3.0,100.0,500.0]
    append!(thicknesses,[cable_thickness,macro_thickness])
    isapprox(maximum(thicknesses)/minimum(thicknesses),ratio;rtol=8.0*eps(Float64),atol=0.0) ||
        error("Multiscale fixture failed its requested thickness-ratio invariant.")
    physical_edges = vcat(0.0,cumsum(thicknesses))
    block_origins = physical_edges[1:end-1]
    block_scales = copy(thicknesses)
    nondimensional = physical_edges/length_scale
    grid = build_mapped_structured_grid(
        (nondimensional,[0.0,1.0],[0.0,1.0]),
        (xi,eta,zeta) -> (length_scale*xi,eta,zeta);
        coordinate_system=:local_nondimensional_layers,
        geometry_hash="analytic-multiscale-$(Int(round(ratio)))",
    )
    condition = maximum([
        cond(grid.jacobians[:,:,i,j,k])
        for i in axes(grid.cell_volumes_cm3,1), j in axes(grid.cell_volumes_cm3,2),
            k in axes(grid.cell_volumes_cm3,3)
    ])
    return Multiscale_Local_Layer_Fixture(
        ratio,layer_ids,thicknesses,block_origins,block_scales,
        nondimensional,length_scale,grid,condition,
    )
end

local_layer_coordinate(fixture::Multiscale_Local_Layer_Fixture,physical_cm::Real) =
    Float64(physical_cm)/fixture.length_scale_cm
physical_layer_coordinate(fixture::Multiscale_Local_Layer_Fixture,local_coordinate::Real) =
    Float64(local_coordinate)*fixture.length_scale_cm

function local_layer_coordinate(fixture::Multiscale_Local_Layer_Fixture,
    layer_index::Integer,physical_cm::Real)
    checkbounds(fixture.layer_ids,layer_index)
    return (Float64(physical_cm)-fixture.block_local_origins_cm[layer_index])/
        fixture.block_length_scales_cm[layer_index]
end

function physical_layer_coordinate(fixture::Multiscale_Local_Layer_Fixture,
    layer_index::Integer,local_coordinate::Real)
    checkbounds(fixture.layer_ids,layer_index)
    return fixture.block_local_origins_cm[layer_index]+
        Float64(local_coordinate)*fixture.block_length_scales_cm[layer_index]
end

struct Mapped_Streaming_Result
    angular_divergence::Array{Float64,3}
    direction::NTuple{3,Float64}
    signed_boundary_current::Float64
    outward_boundary_current::Float64
    inward_boundary_current::Float64
    integrated_divergence::Float64
    particle_balance_residual::Float64
    unique_face_count::Int64
    expected_unique_face_count::Int64
    maximum_internal_face_metric_mismatch::Float64
end

function _mapped_neighbor(index::NTuple{3,Int},axis::Int,side::Int,counts)
    candidate = collect(index)
    candidate[axis] += side == 1 ? -1 : 1
    inside = all(1 <= candidate[value] <= counts[value] for value in 1:3)
    return inside ? (candidate[1],candidate[2],candidate[3]) : nothing
end

"""
    mapped_streaming_divergence(grid, direction, angular_flux; boundary_inflow)

Conservative first-order upwind streaming divergence using cached integrated outward face-area
vectors and cell volumes. `boundary_inflow(axis,side,i,j,k)` supplies incoming angular flux only
where the directed face current is inward. This is a manufactured software kernel, not a full
multigroup collision/source solve.
"""
function mapped_streaming_divergence(grid::Mapped_Structured_Grid,direction,angular_flux;
    boundary_inflow::Function=(axis,side,i,j,k) -> 0.0)
    omega = _atlas_normalize(collect(direction),"mapped streaming direction")
    values = Float64.(angular_flux)
    counts = size(grid.cell_volumes_cm3)
    size(values) == counts || error("Mapped angular flux shape does not match the grid.")
    all(value -> isfinite(value) && value >= 0.0,values) || error(
        "Mapped angular flux must be finite and nonnegative.",
    )
    divergence = zeros(Float64,counts)
    signed_boundary = 0.0
    outward_boundary = 0.0
    inward_boundary = 0.0
    maximum_mismatch = 0.0
    face_incidence = Dict{NTuple{4,Int64},Int64}()
    for i in 1:counts[1], j in 1:counts[2], k in 1:counts[3]
        index = (i,j,k)
        numerator = 0.0
        for axis in 1:3, side in 1:2
            face_coordinate = axis == 1 ? (side == 1 ? i-1 : i) :
                axis == 2 ? (side == 1 ? j-1 : j) : (side == 1 ? k-1 : k)
            face_key = axis == 1 ? (Int64(axis),Int64(face_coordinate),Int64(j),Int64(k)) :
                axis == 2 ? (Int64(axis),Int64(i),Int64(face_coordinate),Int64(k)) :
                (Int64(axis),Int64(i),Int64(j),Int64(face_coordinate))
            face_incidence[face_key] = get(face_incidence,face_key,Int64(0))+1
            area_vector = view(grid.face_area_vectors_cm2,:,side,axis,i,j,k)
            directed_area = dot(omega,area_vector)
            neighbor = _mapped_neighbor(index,axis,side,counts)
            upwind = if directed_area >= 0.0
                values[i,j,k]
            elseif isnothing(neighbor)
                incoming = Float64(boundary_inflow(axis,side,i,j,k))
                isfinite(incoming) && incoming >= 0.0 || error(
                    "Mapped boundary inflow must be finite and nonnegative.",
                )
                incoming
            else
                values[neighbor...]
            end
            current = directed_area*upwind
            numerator += current
            if isnothing(neighbor)
                signed_boundary += current
                if current >= 0.0
                    outward_boundary += current
                else
                    inward_boundary -= current
                end
            elseif side == 2
                opposite = view(
                    grid.face_area_vectors_cm2,:,1,axis,neighbor[1],neighbor[2],neighbor[3],
                )
                scale = max(norm(area_vector),norm(opposite),1.0)
                maximum_mismatch = max(maximum_mismatch,norm(area_vector+opposite)/scale)
            end
        end
        divergence[i,j,k] = numerator/grid.cell_volumes_cm3[i,j,k]
    end
    integrated = sum(divergence.*grid.cell_volumes_cm3)
    expected_faces = (counts[1]+1)*counts[2]*counts[3] +
        counts[1]*(counts[2]+1)*counts[3] + counts[1]*counts[2]*(counts[3]+1)
    length(face_incidence) == expected_faces || error(
        "Mapped topology has a missing or duplicate logical face.",
    )
    for (key,incidence) in face_incidence
        axis = key[1]
        face_coordinate = key[axis+1]
        axis_count = counts[axis]
        expected_incidence = face_coordinate in (0,axis_count) ? 1 : 2
        incidence == expected_incidence || error("Mapped logical face incidence is invalid.")
    end
    result = Mapped_Streaming_Result(
        divergence,(omega[1],omega[2],omega[3]),signed_boundary,outward_boundary,
        inward_boundary,integrated,integrated-signed_boundary,Int64(length(face_incidence)),
        Int64(expected_faces),maximum_mismatch,
    )
    abs(result.particle_balance_residual) <=
        1.0e-10*max(abs(integrated),abs(signed_boundary),1.0) || error(
        "Mapped streaming particle balance does not close.",
    )
    maximum_mismatch <= 1.0e-9 || error("Mapped internal face metrics do not match.")
    return result
end

function mapped_streaming_closure(result::Mapped_Streaming_Result;particle_energy_MeV::Real)
    energy = Float64(particle_energy_MeV)
    isfinite(energy) && energy >= 0.0 || error("Particle energy must be nonnegative.")
    return Dict{String,Any}(
        "schema" => "radiant.mapped_streaming_closure/v1",
        "particle_balance_residual" => result.particle_balance_residual,
        "energy_current_balance_residual_MeV" => energy*result.particle_balance_residual,
        "signed_boundary_current" => result.signed_boundary_current,
        "outward_boundary_current" => result.outward_boundary_current,
        "inward_boundary_current" => result.inward_boundary_current,
        "unique_face_count" => result.unique_face_count,
        "expected_unique_face_count" => result.expected_unique_face_count,
        "maximum_internal_face_metric_mismatch" =>
            result.maximum_internal_face_metric_mismatch,
        "software_manufactured_only" => true,
    )
end

function transform_polar_vector(transform,vector)
    matrix = Float64.(transform)
    size(matrix) == (3,3) && isapprox(matrix'*matrix,Matrix{Float64}(I,3,3);atol=1.0e-10) ||
        error("Polar-vector transform must be orthogonal.")
    return matrix*Float64.(vector)
end


function transform_axial_vector(transform,vector)
    matrix = Float64.(transform)
    size(matrix) == (3,3) && isapprox(matrix'*matrix,Matrix{Float64}(I,3,3);atol=1.0e-10) ||
        error("Axial-vector transform must be orthogonal.")
    return det(matrix)*matrix*Float64.(vector)
end

struct Curved_Performance_Target
    target_id::String
    maximum_sweep_per_unknown_ratio::Float64
    maximum_equal_error_solve_ratio::Float64
    maximum_memory_ratio::Float64
    maximum_cached_atlas_ratio::Float64
    maximum_cells::Int64
end

function Curved_Performance_Target(;target_id::AbstractString,
    maximum_sweep_per_unknown_ratio::Real=1.25,
    maximum_equal_error_solve_ratio::Real=1.35,
    maximum_memory_ratio::Real=1.25,
    maximum_cached_atlas_ratio::Real=1.50,
    maximum_cells::Integer)
    isempty(target_id) && error("Performance target ID cannot be empty.")
    limits = Float64[
        maximum_sweep_per_unknown_ratio,maximum_equal_error_solve_ratio,
        maximum_memory_ratio,maximum_cached_atlas_ratio,
    ]
    all(value -> isfinite(value) && value >= 1.0,limits) || error(
        "Performance ratios must be finite and at least one.",
    )
    maximum_cells >= 1 || error("Performance target cell capacity must be positive.")
    return Curved_Performance_Target(String(target_id),limits...,Int64(maximum_cells))
end

function benchmark_curved_pipeline(target::Curved_Performance_Target,grid_builder::Function;
    baseline_sweep::Function,improved_sweep::Function,scoring::Function,
    baseline_remap::Function,cached_atlas_remap::Function,samples::Integer=7)
    samples >= 3 || error("Curved timing requires at least three samples.")
    grid_ref = Ref{Any}()
    preprocessing_s = @elapsed grid_ref[] = grid_builder()
    built = grid_ref[]
    baseline_grid = built isa NamedTuple ? built.baseline : built
    improved_grid = built isa NamedTuple ? built.improved : built
    prod(size(improved_grid.cell_volumes_cm3)) <= target.maximum_cells || error(
        "Benchmark grid exceeds the preregistered cell capacity.",
    )
    warmup_s = @elapsed begin
        baseline_warm = baseline_sweep(baseline_grid)
        improved_warm = improved_sweep(improved_grid)
        scoring(improved_warm,improved_grid)
        baseline_remap(baseline_warm,baseline_grid)
        cached_atlas_remap(improved_warm,improved_grid)
    end
    baseline_times = Float64[]
    improved_times = Float64[]
    baseline_remap_times = Float64[]
    cached_remap_times = Float64[]
    baseline_output = baseline_sweep(baseline_grid)
    improved_output = improved_sweep(improved_grid)
    for _ in 1:samples
        push!(baseline_times,@elapsed baseline_output = baseline_sweep(baseline_grid))
        push!(improved_times,@elapsed improved_output = improved_sweep(improved_grid))
        push!(baseline_remap_times,@elapsed baseline_remap(baseline_output,baseline_grid))
        push!(cached_remap_times,@elapsed cached_atlas_remap(improved_output,improved_grid))
    end
    median_value(values) = sort(values)[cld(length(values),2)]
    baseline_sweep_s = median_value(baseline_times)
    improved_sweep_s = median_value(improved_times)
    baseline_remap_s = median_value(baseline_remap_times)
    cached_remap_s = median_value(cached_remap_times)
    scoring_ref = Ref{Any}()
    scoring_s = @elapsed scoring_ref[] = scoring(improved_output,improved_grid)
    unknowns = improved_output isa Mapped_Streaming_Result ?
        length(improved_output.angular_divergence) : length(improved_output)
    baseline_sweep_per_unknown = baseline_sweep_s/unknowns
    improved_sweep_per_unknown = improved_sweep_s/unknowns
    sweep_ratio = improved_sweep_per_unknown/max(baseline_sweep_per_unknown,eps(Float64))
    matched_unknown_kernel_time_ratio =
        improved_sweep_s/max(baseline_sweep_s,eps(Float64))
    baseline_bytes = Base.summarysize(baseline_grid)+Base.summarysize(baseline_output)
    improved_bytes = Base.summarysize(improved_grid)+Base.summarysize(improved_output)
    memory_ratio = improved_bytes/max(baseline_bytes,1)
    atlas_ratio = cached_remap_s/max(baseline_remap_s,eps(Float64))
    passed = sweep_ratio <= target.maximum_sweep_per_unknown_ratio &&
        memory_ratio <= target.maximum_memory_ratio &&
        atlas_ratio <= target.maximum_cached_atlas_ratio
    return Dict{String,Any}(
        "schema" => "radiant.curved_performance_benchmark/v2",
        "target_id" => target.target_id,
        "target_preregistered_before_timing" => true,
        "cell_count" => prod(size(improved_grid.cell_volumes_cm3)),
        "preprocessing_s" => preprocessing_s,"warmup_s" => warmup_s,"scoring_s" => scoring_s,
        "baseline" => Dict{String,Any}(
            "sweep_s" => baseline_sweep_s,"sweep_per_unknown_s" => baseline_sweep_per_unknown,
            "remap_s" => baseline_remap_s,"memory_bytes" => baseline_bytes,
        ),
        "improved" => Dict{String,Any}(
            "sweep_s" => improved_sweep_s,"sweep_per_unknown_s" => improved_sweep_per_unknown,
            "cached_atlas_remap_s" => cached_remap_s,"memory_bytes" => improved_bytes,
        ),
        "ratios" => Dict{String,Any}(
            "sweep_per_unknown" => sweep_ratio,
            "matched_unknown_kernel_time" => matched_unknown_kernel_time_ratio,
            "equal_error_solve" => nothing,
            "memory" => memory_ratio,"cached_atlas" => atlas_ratio,
        ),
        "thresholds" => Dict{String,Any}(
            "sweep_per_unknown" => target.maximum_sweep_per_unknown_ratio,
            "equal_error_solve" => target.maximum_equal_error_solve_ratio,
            "memory" => target.maximum_memory_ratio,
            "cached_atlas" => target.maximum_cached_atlas_ratio,
        ),
        "measurement_status" => Dict{String,Any}(
            "sweep_per_unknown" => "MEASURED",
            "matched_unknown_kernel_time" => "MEASURED_DIAGNOSTIC_ONLY",
            "equal_error_solve" => "NOT_MEASURED_DIAGNOSTIC_ONLY",
            "memory" => "MEASURED","cached_atlas" => "MEASURED",
        ),
        "passed" => passed,"physical_qualification" => false,
    )
end
