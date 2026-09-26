"""
    SN_Boundary_Flux

Actual outgoing discrete-ordinate face moments, captured during the final physical sweep.
`faces[f]` has axes `(group, ordinate, face moment, tangent cell 1, tangent cell 2)`;
face order is x-, x+, y-, y+, z-, z+. Only active faces exist. Moment 1 is the
face average, integrated over the energy group. Higher moments are retained, not clipped.
1D current is per unit transverse area; 2D current is per unit extruded length.
Reflective/periodic face crossings are retained and are not escaped particles.
Scores retain the raw input-source normalization; divide by the applicable Fixed_Sources
normalization before claiming per-history values. Energy contraction requires explicit
native normalized-Legendre metadata; legacy captures support particle scoring only.
"""
struct SN_Boundary_Flux
    dimension::Int64
    directions::Vector{Vector{Float64}}
    weights::Vector{Float64}
    widths::Vector{Vector{Float64}}
    energy_boundaries::Vector{Float64}
    boundary_conditions::Vector{Int64}
    faces::Vector{Array{Float64,5}}
    convergence::Vector{NamedTuple}
    orders::Union{Nothing,Vector{Int64}}
    fully_coupled::Union{Nothing,Bool}
    is_csd::Union{Nothing,Bool}
    basis::Symbol
    outer_convergence::Base.RefValue{NamedTuple}
end

function SN_Boundary_Flux(dimension,Ω,w,Δs,Eb,boundary_conditions,Ng,Nm;
    orders=nothing,fully_coupled=nothing,is_csd=nothing)
    dimension in 1:3 || error("SN boundary scoring supports Cartesian dimensions 1, 2 and 3.")
    length(Ω) == 3 && all(length(v) == length(w) for v in Ω) || error("Invalid boundary quadrature shape.")
    all(isfinite, w) && all(w .> 0) || error("Boundary quadrature weights must be finite and positive.")
    all(all(isfinite,v) for v in Ω) || error("Boundary directions must be finite.")
    length(boundary_conditions) == 2dimension || error("Boundary labels must match active faces.")
    Ng > 0 && length(Nm) >= dimension && all(Nm[a] > 0 for a in 1:dimension) ||
        error("Boundary group and face moment counts must be positive.")
    length(Eb) == Ng+1 && all(isfinite,Eb) && all(diff(Eb) .< 0) || error("Boundary energy edges must be strictly decreasing.")
    for a in 1:dimension
        all(isfinite,Δs[a]) && all(Δs[a] .> 0) || error("Boundary cell widths must be finite and positive.")
    end
    metadata = (orders,fully_coupled,is_csd)
    all(isnothing,metadata) || all(!isnothing(v) for v in metadata) ||
        error("Boundary basis orders, coupling and CSD flag must be supplied together.")
    native_orders = nothing
    if orders !== nothing
        length(orders) == 4 && all(v isa Integer && v >= 1 for v in orders) ||
            error("Native boundary orders must contain four positive integers.")
        fully_coupled isa Bool && is_csd isa Bool || error("Boundary scheme flags must be Bool.")
        native_orders = Int64.(orders)
        all(native_orders[a] == 1 for a in dimension+1:3) || error("Inactive axis orders must be one.")
        is_csd || native_orders[4] == 1 || error("Non-CSD capture cannot carry energy slopes.")
        for a in 1:dimension
            face_orders = native_orders[[b for b in 1:4 if b != a]]
            expected = fully_coupled ? prod(big.(face_orders)) : 1+sum(big.(face_orders).-1)
            expected == Nm[a] || error("Boundary face shape does not match native basis orders.")
        end
    end
    faces = Array{Float64,5}[]
    for a in 1:dimension
        tangents = [b for b in 1:dimension if b != a]
        shape = [length(Δs[b]) for b in tangents]
        while length(shape) < 2; push!(shape,1); end
        for side in 1:2
            push!(faces,zeros(Ng,length(w),Nm[a],shape...))
        end
    end
    status = NamedTuple[(converged=false,terminal=:not_solved,iterations=0,
        solver_residual=NaN,reconstruction_residual=NaN,tolerance=NaN,iteration_cap=0) for g in 1:Ng]
    # Native lower-dimensional geometry leaves inactive Δs slots undefined. Retain
    # only physical axes and canonical empty placeholders for inactive axes, so
    # generation compatibility checks never traverse uninitialized references.
    widths = [a <= dimension ? copy(Δs[a]) : Float64[] for a in 1:3]
    outer = Ref{NamedTuple}((converged=false,required=missing,residual=NaN,tolerance=NaN,
        iterations=0,iteration_cap=0))
    return SN_Boundary_Flux(Int64(dimension),deepcopy(Ω),copy(w),widths,copy(Eb),copy(boundary_conditions),faces,status,
        native_orders,fully_coupled,is_csd,orders === nothing ? :unknown : :native_group_integrated_normalized_legendre,outer)
end

# The optional capture tuple contains this generation, group and ordinate. Capture the
# raw sweep value before any half-range angular projection used for boundary iteration.
function _record_sn_face!(capture,axis,cosine,values,tangent1=1,tangent2=1)
    capture === nothing && return nothing
    cosine == 0 && return nothing
    data,group,ordinate = capture
    face = 2*axis - (cosine < 0 ? 1 : 0)
    data.faces[face][group,ordinate,:,tangent1,tangent2] .= values
    return nothing
end

"""
    get_outgoing_current(data::SN_Boundary_Flux)

Return outward crossing particle current by `(face, energy group)`, using quadrature
weights, outward normal cosine and the actual tangential cell measures. Negative
face values remain negative. This score is not net current or leakage on closed faces.
It retains raw input-source normalization, not an automatic per-history normalization.
"""
function get_outgoing_current(data::SN_Boundary_Flux)
    Ng = size(data.faces[1],1)
    current = zeros(2*data.dimension,Ng)
    for face in 1:2*data.dimension
        axis = (face+1) ÷ 2
        sign = isodd(face) ? -1.0 : 1.0
        tangents = [a for a in 1:data.dimension if a != axis]
        for g in 1:Ng, n in eachindex(data.weights), j in axes(data.faces[face],4), k in axes(data.faces[face],5)
            cosine = sign * data.directions[axis][n]
            cosine > 0 || continue
            measure = isempty(tangents) ? 1.0 : data.widths[tangents[1]][j]
            if length(tangents) == 2; measure *= data.widths[tangents[2]][k]; end
            current[face,g] += data.weights[n]*cosine*measure*data.faces[face][g,n,1,j,k]
        end
    end
    return current
end

function _validate_sn_energy_basis(data::SN_Boundary_Flux)
    data.orders !== nothing && data.fully_coupled !== nothing && data.is_csd !== nothing &&
        data.basis == :native_group_integrated_normalized_legendre ||
        error("Energy current requires authoritative native boundary basis metadata.")
    length(data.orders) == 4 && all(data.orders .>= 1) || error("Invalid captured basis orders.")
    all(data.orders[a] == 1 for a in data.dimension+1:3) || error("Inactive captured axis order differs.")
    data.is_csd || data.orders[4] == 1 || error("Non-CSD capture has an energy slope.")
    length(data.faces) == 2*data.dimension && length(data.boundary_conditions) == length(data.faces) ||
        error("Captured boundary face count differs.")
    Ng = length(data.energy_boundaries)-1
    Ng > 0 && all(isfinite,data.energy_boundaries) && all(diff(data.energy_boundaries) .< 0) ||
        error("Captured energy boundaries must be finite and decreasing.")
    length(data.directions) == 3 && all(length(v) == length(data.weights) && all(isfinite,v) for v in data.directions) &&
        all(v -> isfinite(v) && v > 0,data.weights) || error("Invalid captured quadrature.")
    for a in 1:data.dimension
        all(v -> isfinite(v) && v > 0,data.widths[a]) || error("Invalid captured cell widths.")
        face_orders = data.orders[[b for b in 1:4 if b != a]]
        expected = data.fully_coupled ? prod(big.(face_orders)) : 1+sum(big.(face_orders).-1)
        tangents = [b for b in 1:data.dimension if b != a]
        shape = [length(data.widths[b]) for b in tangents]
        while length(shape) < 2; push!(shape,1); end
        for f in (2a-1,2a)
            size(data.faces[f]) == (Ng,length(data.weights),expected,shape...) ||
                error("Captured face shape differs from basis or geometry.")
            all(isfinite,data.faces[f]) || error("Captured boundary moments must be finite.")
        end
    end
    return nothing
end

"""
Integrate represented crossing kinetic energy by `(face,group)`. Native face moments
are group integrated and use normalized Legendre polynomials with upper E at +1.
Orthogonality leaves only E midpoint times F0 plus ΔE/(2sqrt(3)) times FE1.
The FE1 mode is index 2 only when the authoritative energy order is at least two.
Signed diagnostics are retained; reflective and periodic crossings are included.
"""
function get_outgoing_energy_current(data::SN_Boundary_Flux)
    _validate_sn_energy_basis(data)
    Ng = length(data.energy_boundaries)-1
    current = zeros(2*data.dimension,Ng)
    for face in eachindex(data.faces)
        axis = (face+1) ÷ 2
        outward_sign = isodd(face) ? -1.0 : 1.0
        tangents = [a for a in 1:data.dimension if a != axis]
        for g in 1:Ng
            upper,lower = data.energy_boundaries[g:g+1]
            midpoint = upper/2+lower/2
            slope_weight = (upper-lower)/(2sqrt(3.0))
            isfinite(midpoint) && isfinite(slope_weight) || error("Energy integration weights overflowed.")
            for n in eachindex(data.weights), j in axes(data.faces[face],4), k in axes(data.faces[face],5)
                cosine = outward_sign*data.directions[axis][n]
                cosine > 0 || continue
                measure = isempty(tangents) ? 1.0 : data.widths[tangents[1]][j]
                if length(tangents) == 2; measure *= data.widths[tangents[2]][k]; end
                energy = midpoint*data.faces[face][g,n,1,j,k]
                if data.orders[4] >= 2; energy += slope_weight*data.faces[face][g,n,2,j,k]; end
                contribution = data.weights[n]*cosine*measure*energy
                isfinite(measure) && isfinite(energy) && isfinite(contribution) || error("Boundary energy score overflowed.")
                current[face,g] += contribution
                isfinite(current[face,g]) || error("Accumulated boundary energy score overflowed.")
            end
        end
    end
    return current
end

"""Void-only escaped kinetic energy. Closed faces are zero; unknown labels fail closed."""
function get_escaped_energy_current(data::SN_Boundary_Flux)
    all(code in (0,1,2) for code in data.boundary_conditions) || error("Unknown boundary condition code in escape score.")
    score = get_outgoing_energy_current(data)
    for face in eachindex(data.boundary_conditions)
        data.boundary_conditions[face] == 0 || (score[face,:] .= 0.0)
    end
    return score
end
