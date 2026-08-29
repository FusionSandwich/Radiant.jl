function _atlas_normalize(vector::AbstractVector{<:Real},label::AbstractString)
    length(vector) == 3 || error("$(label) must be a three-vector.")
    value = Float64.(vector)
    all(isfinite,value) || error("$(label) must be finite.")
    magnitude = norm(value)
    magnitude > 0.0 || error("$(label) cannot be zero.")
    return value/magnitude
end

function _atlas_choose_initial_normal(tangent::Vector{Float64})
    candidates = (
        [1.0,0.0,0.0],
        [0.0,1.0,0.0],
        [0.0,0.0,1.0],
    )
    candidate = candidates[argmin([abs(dot(tangent,value)) for value in candidates])]
    normal = candidate-dot(candidate,tangent)*tangent
    return _atlas_normalize(normal,"initial tape normal")
end

function _atlas_rodrigues(
    vector::Vector{Float64},
    axis::Vector{Float64},
    angle::Float64,
)
    return vector*cos(angle) + cross(axis,vector)*sin(angle) +
           axis*dot(axis,vector)*(1.0-cos(angle))
end

function _atlas_parallel_transport(
    normal::Vector{Float64},
    tangent_old::Vector{Float64},
    tangent_new::Vector{Float64};
    tolerance::Real=1.0e-12,
)
    cosine = clamp(dot(tangent_old,tangent_new),-1.0,1.0)
    rotation_axis = cross(tangent_old,tangent_new)
    sine = norm(rotation_axis)
    if sine ≤ tolerance
        cosine < 0.0 && error(
            "Parallel transport is undefined across an unresolved 180-degree tangent reversal.",
        )
        transported = normal
    else
        axis = rotation_axis/sine
        transported = _atlas_rodrigues(normal,axis,atan(sine,cosine))
    end
    transported .-= dot(transported,tangent_new)*tangent_new
    return _atlas_normalize(transported,"parallel-transported tape normal")
end

"""One locally flat tape patch with a right-handed `(tangent,width,normal)` frame."""
struct Tape_Atlas_Patch
    patch_id::Int64
    s_start_cm::Float64
    s_stop_cm::Float64
    centroid_cm::NTuple{3,Float64}
    tangent::NTuple{3,Float64}
    width_axis::NTuple{3,Float64}
    normal::NTuple{3,Float64}
    width_cm::Float64
    thickness_cm::Float64
    turn_angle_rad::Float64
    curvature_cm_inv::Float64
    torsion_cm_inv::Float64
    sagitta_bound_cm::Float64
    magnetic_field_T::NTuple{3,Float64}
    field_variation_bound_T::Float64
    source_patch_ids::Vector{Int64}
    metadata::Dict{String,String}

    function Tape_Atlas_Patch(;
        patch_id::Integer,
        s_start_cm::Real,
        s_stop_cm::Real,
        centroid_cm::AbstractVector{<:Real},
        tangent::AbstractVector{<:Real},
        width_axis::AbstractVector{<:Real},
        normal::AbstractVector{<:Real},
        width_cm::Real,
        thickness_cm::Real,
        turn_angle_rad::Real=0.0,
        curvature_cm_inv::Real=0.0,
        torsion_cm_inv::Real=0.0,
        sagitta_bound_cm::Real=0.0,
        magnetic_field_T::AbstractVector{<:Real}=[0.0,0.0,0.0],
        field_variation_bound_T::Real=0.0,
        source_patch_ids::AbstractVector{<:Integer}=Int64[],
        metadata::AbstractDict=Dict{String,String}(),
        frame_tolerance::Real=1.0e-9,
    )
        patch_id ≥ 1 || error("Tape atlas patch identifiers must be positive.")
        start_value = Float64(s_start_cm)
        stop_value = Float64(s_stop_cm)
        isfinite(start_value) && isfinite(stop_value) && stop_value > start_value || error(
            "Tape patch arc-length interval must be finite and increasing.",
        )
        centroid = Float64.(centroid_cm)
        length(centroid) == 3 && all(isfinite,centroid) || error(
            "Tape patch centroid must be a finite three-vector.",
        )
        t = _atlas_normalize(tangent,"tape tangent")
        w = _atlas_normalize(width_axis,"tape width axis")
        n = _atlas_normalize(normal,"tape normal")
        abs(dot(t,w)) ≤ frame_tolerance && abs(dot(t,n)) ≤ frame_tolerance &&
            abs(dot(w,n)) ≤ frame_tolerance || error("Tape patch frame must be orthogonal.")
        dot(cross(t,w),n) ≥ 1.0-10.0*frame_tolerance || error(
            "Tape patch frame must be right-handed: tangent × width = normal.",
        )
        width_value = Float64(width_cm)
        thickness_value = Float64(thickness_cm)
        isfinite(width_value) && width_value > 0.0 || error("Tape width must be positive.")
        isfinite(thickness_value) && thickness_value > 0.0 || error(
            "Tape thickness must be positive.",
        )
        geometry_values = Float64[
            turn_angle_rad,curvature_cm_inv,torsion_cm_inv,sagitta_bound_cm,
            field_variation_bound_T,
        ]
        all(isfinite,geometry_values) || error("Tape patch geometry metrics must be finite.")
        geometry_values[1] ≥ 0.0 || error("Tape turn angle must be nonnegative.")
        geometry_values[2] ≥ 0.0 || error("Tape curvature magnitude must be nonnegative.")
        geometry_values[4] ≥ 0.0 || error("Tape sagitta bound must be nonnegative.")
        geometry_values[5] ≥ 0.0 || error("Field variation bound must be nonnegative.")
        field = Float64.(magnetic_field_T)
        length(field) == 3 && all(isfinite,field) || error(
            "Tape patch field must be a finite three-vector.",
        )
        identifiers = Int64.(source_patch_ids)
        length(unique(identifiers)) == length(identifiers) || error(
            "Source patch identifiers must be unique within a tape patch.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            Int64(patch_id),start_value,stop_value,(centroid[1],centroid[2],centroid[3]),
            (t[1],t[2],t[3]),(w[1],w[2],w[3]),(n[1],n[2],n[3]),width_value,
            thickness_value,geometry_values[1],geometry_values[2],geometry_values[3],
            geometry_values[4],(field[1],field[2],field[3]),geometry_values[5],identifiers,
            metadata_string,
        )
    end
end

struct Piecewise_Flat_Tape_Atlas
    schema::String
    patches::Vector{Tape_Atlas_Patch}
    total_length_cm::Float64
    periodic::Bool
    frame_convention::String
    geometry_hash::String
    metadata::Dict{String,String}

    function Piecewise_Flat_Tape_Atlas(
        patches::AbstractVector{Tape_Atlas_Patch};
        periodic::Bool=false,
        frame_convention::AbstractString="parallel-transport",
        geometry_hash::AbstractString="unbound",
        metadata::AbstractDict=Dict{String,String}(),
        continuity_tolerance_cm::Real=1.0e-10,
    )
        patch_vector = Tape_Atlas_Patch[patches...]
        isempty(patch_vector) && error("Tape atlas requires at least one patch.")
        identifiers = getfield.(patch_vector,:patch_id)
        identifiers == collect(1:length(patch_vector)) || error(
            "Tape atlas patch identifiers must be contiguous and one-based.",
        )
        abs(patch_vector[1].s_start_cm) ≤ continuity_tolerance_cm || error(
            "Tape atlas must begin at arc length zero.",
        )
        for index in 2:length(patch_vector)
            isapprox(
                patch_vector[index-1].s_stop_cm,patch_vector[index].s_start_cm;
                rtol=0.0,atol=continuity_tolerance_cm,
            ) || error("Tape atlas arc-length intervals must be contiguous.")
        end
        isempty(frame_convention) && error("Tape atlas frame convention cannot be empty.")
        isempty(geometry_hash) && error("Tape atlas geometry hash cannot be empty.")
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            "radiant.hts.piecewise_flat_atlas/v1",patch_vector,
            patch_vector[end].s_stop_cm,periodic,String(frame_convention),
            String(geometry_hash),metadata_string,
        )
    end
end

function patch_frame(this::Tape_Atlas_Patch)
    return hcat(collect(this.tangent),collect(this.width_axis),collect(this.normal))
end

function global_to_patch_coordinates(
    this::Tape_Atlas_Patch,
    position_cm::AbstractVector{<:Real},
)
    length(position_cm) == 3 || error("Global position must be a three-vector.")
    offset = Float64.(position_cm)-collect(this.centroid_cm)
    return patch_frame(this)'*offset
end

function patch_to_global_coordinates(
    this::Tape_Atlas_Patch,
    local_position_cm::AbstractVector{<:Real},
)
    length(local_position_cm) == 3 || error("Local position must be a three-vector.")
    return collect(this.centroid_cm)+patch_frame(this)*Float64.(local_position_cm)
end

function global_direction_to_patch(
    this::Tape_Atlas_Patch,
    direction::AbstractVector{<:Real},
)
    value = _atlas_normalize(direction,"global direction")
    return patch_frame(this)'*value
end

function patch_direction_to_global(
    this::Tape_Atlas_Patch,
    direction::AbstractVector{<:Real},
)
    value = _atlas_normalize(direction,"local direction")
    return patch_frame(this)*value
end

function electromagnetic_field_for_patch(this::Tape_Atlas_Patch)
    field = Electromagnetic_Field()
    set_magnetic_field(field,collect(this.magnetic_field_T))
    return field
end

function _atlas_discrete_torsion(
    tangents::Vector{Vector{Float64}},
    lengths::Vector{Float64},
    index::Int,
)
    if index ≤ 1 || index ≥ length(tangents)
        return 0.0
    end
    cross_before = cross(tangents[index-1],tangents[index])
    cross_after = cross(tangents[index],tangents[index+1])
    norm(cross_before) ≤ 1.0e-12 && return 0.0
    norm(cross_after) ≤ 1.0e-12 && return 0.0
    binormal_before = cross_before/norm(cross_before)
    binormal_after = cross_after/norm(cross_after)
    signed_angle = atan(
        dot(cross(binormal_before,binormal_after),tangents[index]),
        clamp(dot(binormal_before,binormal_after),-1.0,1.0),
    )
    length_scale = 0.5*(lengths[index-1]+lengths[index])
    return signed_angle/length_scale
end

"""
    build_piecewise_flat_tape_atlas(centerline_cm; ...)

Construct one flat patch per centerline segment using a minimal-twist parallel-transport frame.
Curvature and torsion are discrete screening estimates. `sagitta_bound_cm = κL²/8` is recorded per
patch; production acceptance still requires protected-response convergence against a curved
reference segment.
"""
function build_piecewise_flat_tape_atlas(
    centerline_cm::AbstractMatrix{<:Real};
    width_cm::Real,
    thickness_cm::Real,
    initial_normal::Union{Nothing,AbstractVector{<:Real}}=nothing,
    magnetic_field::Union{Nothing,Abstract_Spatial_Magnetic_Field}=nothing,
    source_patch_ids::Union{Nothing,AbstractVector}=nothing,
    periodic::Bool=false,
    geometry_hash::AbstractString="unbound",
    metadata::AbstractDict=Dict{String,String}(),
)
    points = Float64.(centerline_cm)
    size(points,2) == 3 || error("Tape centerline must have shape (N,3).")
    size(points,1) ≥ 2 || error("Tape centerline requires at least two points.")
    all(isfinite,points) || error("Tape centerline coordinates must be finite.")
    segment_count = size(points,1)-1
    lengths = zeros(Float64,segment_count)
    tangents = Vector{Vector{Float64}}(undef,segment_count)
    for index in 1:segment_count
        delta = vec(points[index+1,:]-points[index,:])
        lengths[index] = norm(delta)
        lengths[index] > 0.0 || error("Tape centerline contains a zero-length segment.")
        tangents[index] = delta/lengths[index]
    end

    normal = isnothing(initial_normal) ? _atlas_choose_initial_normal(tangents[1]) : begin
        candidate = _atlas_normalize(initial_normal,"initial tape normal")
        candidate .-= dot(candidate,tangents[1])*tangents[1]
        _atlas_normalize(candidate,"orthogonalized initial tape normal")
    end
    normals = Vector{Vector{Float64}}(undef,segment_count)
    normals[1] = normal
    for index in 2:segment_count
        normals[index] = _atlas_parallel_transport(
            normals[index-1],tangents[index-1],tangents[index],
        )
    end

    turns = zeros(Float64,segment_count)
    curvatures = zeros(Float64,segment_count)
    for index in 1:segment_count
        adjacent_angles = Float64[]
        adjacent_curvatures = Float64[]
        if index > 1
            angle = acos(clamp(dot(tangents[index-1],tangents[index]),-1.0,1.0))
            push!(adjacent_angles,angle)
            push!(adjacent_curvatures,angle/(0.5*(lengths[index-1]+lengths[index])))
        end
        if index < segment_count
            angle = acos(clamp(dot(tangents[index],tangents[index+1]),-1.0,1.0))
            push!(adjacent_angles,angle)
            push!(adjacent_curvatures,angle/(0.5*(lengths[index]+lengths[index+1])))
        end
        turns[index] = isempty(adjacent_angles) ? 0.0 : maximum(adjacent_angles)
        curvatures[index] = isempty(adjacent_curvatures) ? 0.0 :
            sum(adjacent_curvatures)/length(adjacent_curvatures)
    end

    s_value = 0.0
    patches = Tape_Atlas_Patch[]
    for index in 1:segment_count
        tangent = tangents[index]
        normal_vector = normals[index]
        width_axis = _atlas_normalize(cross(normal_vector,tangent),"tape width axis")
        centroid = vec(0.5.*(points[index,:]+points[index+1,:]))
        field_global = isnothing(magnetic_field) ? zeros(Float64,3) :
            field_at(magnetic_field,centroid)
        frame = hcat(tangent,width_axis,normal_vector)
        field_local = frame'*field_global
        field_variation = isnothing(magnetic_field) ? 0.0 : field_variation_bound(
            magnetic_field,
            vcat(
                reshape(vec(points[index,:]),1,3),
                reshape(centroid,1,3),
                reshape(vec(points[index+1,:]),1,3),
            ),
        )
        identifiers = isnothing(source_patch_ids) ? Int64[] :
            Int64.(source_patch_ids[index])
        patch = Tape_Atlas_Patch(
            patch_id=index,
            s_start_cm=s_value,
            s_stop_cm=s_value+lengths[index],
            centroid_cm=centroid,
            tangent=tangent,
            width_axis=width_axis,
            normal=normal_vector,
            width_cm=width_cm,
            thickness_cm=thickness_cm,
            turn_angle_rad=turns[index],
            curvature_cm_inv=curvatures[index],
            torsion_cm_inv=_atlas_discrete_torsion(tangents,lengths,index),
            sagitta_bound_cm=curvatures[index]*lengths[index]^2/8.0,
            magnetic_field_T=field_local,
            field_variation_bound_T=field_variation,
            source_patch_ids=identifiers,
            metadata=Dict(
                "field_frame" => "patch-local",
                "centerline_segment" => string(index),
            ),
        )
        push!(patches,patch)
        s_value += lengths[index]
    end
    return Piecewise_Flat_Tape_Atlas(
        patches;
        periodic=periodic,
        geometry_hash=geometry_hash,
        metadata=metadata,
    )
end

function get_atlas_patch(this::Piecewise_Flat_Tape_Atlas,patch_id::Integer)
    1 ≤ patch_id ≤ length(this.patches) || error("Tape atlas patch identifier is out of range.")
    return this.patches[Int(patch_id)]
end

function get_atlas_patch_at_arc_length(
    this::Piecewise_Flat_Tape_Atlas,
    arc_length_cm::Real;
    tolerance_cm::Real=1.0e-12,
)
    value = Float64(arc_length_cm)
    value ≥ -tolerance_cm && value ≤ this.total_length_cm+tolerance_cm || error(
        "Arc length lies outside the tape atlas.",
    )
    for patch in this.patches
        if value ≥ patch.s_start_cm-tolerance_cm && value ≤ patch.s_stop_cm+tolerance_cm
            return patch
        end
    end
    error("Arc length could not be assigned to a tape patch.")
end

function atlas_refinement_report(
    this::Piecewise_Flat_Tape_Atlas;
    maximum_turn_angle_rad::Real,
    maximum_sagitta_cm::Real,
    maximum_field_variation_T::Real,
)
    limits = Float64[
        maximum_turn_angle_rad,maximum_sagitta_cm,maximum_field_variation_T,
    ]
    all(isfinite,limits) && all(limits .≥ 0.0) || error(
        "Tape atlas refinement limits must be finite and nonnegative.",
    )
    failed_turn = Int64[]
    failed_sagitta = Int64[]
    failed_field = Int64[]
    for patch in this.patches
        patch.turn_angle_rad > limits[1] && push!(failed_turn,patch.patch_id)
        patch.sagitta_bound_cm > limits[2] && push!(failed_sagitta,patch.patch_id)
        patch.field_variation_bound_T > limits[3] && push!(failed_field,patch.patch_id)
    end
    passing = isempty(failed_turn) && isempty(failed_sagitta) && isempty(failed_field)
    return Dict{String,Any}(
        "schema" => "radiant.hts.atlas_refinement/v1",
        "geometry_hash" => this.geometry_hash,
        "patch_count" => length(this.patches),
        "total_length_cm" => this.total_length_cm,
        "maximum_observed_turn_angle_rad" => maximum(getfield.(this.patches,:turn_angle_rad)),
        "maximum_observed_sagitta_cm" => maximum(getfield.(this.patches,:sagitta_bound_cm)),
        "maximum_observed_field_variation_T" => maximum(
            getfield.(this.patches,:field_variation_bound_T),
        ),
        "failed_turn_patch_ids" => failed_turn,
        "failed_sagitta_patch_ids" => failed_sagitta,
        "failed_field_patch_ids" => failed_field,
        "screening_pass" => passing,
        "production_pass" => false,
        "production_pass_reason" =>
            "A curved OpenSn or Geant4 response comparison is still required.",
    )
end
