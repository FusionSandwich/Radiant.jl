"""
    SN_Boundary_Flux

Actual outgoing discrete-ordinate face moments, captured during the final physical sweep.
`faces[f]` has axes `(group, ordinate, face moment, tangent cell 1, tangent cell 2)`;
face order is x-, x+, y-, y+, z-, z+. Only active faces exist. Moment 1 is the
face average, integrated over the energy group. Higher moments are retained, not clipped.
1D current is per unit transverse area; 2D current is per unit extruded length.
Reflective/periodic face crossings are retained and are not escaped particles.
Scores retain the raw input-source normalization; divide by the applicable Fixed_Sources
normalization before claiming per-history values. Exact escaped energy is unsupported:
group midpoint times current is generally insufficient, and the energy/face moment
contraction and moment mapping have not been qualified by this first increment.
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
end

function SN_Boundary_Flux(dimension,Ω,w,Δs,Eb,boundary_conditions,Ng,Nm)
    dimension in 1:3 || error("SN boundary scoring supports Cartesian dimensions 1, 2 and 3.")
    length(Ω) == 3 && all(length(v) == length(w) for v in Ω) || error("Invalid boundary quadrature shape.")
    all(isfinite, w) && all(w .> 0) || error("Boundary quadrature weights must be finite and positive.")
    all(all(isfinite,v) for v in Ω) || error("Boundary directions must be finite.")
    length(Eb) == Ng+1 && all(isfinite,Eb) && all(diff(Eb) .< 0) || error("Boundary energy edges must be strictly decreasing.")
    for a in 1:dimension
        all(isfinite,Δs[a]) && all(Δs[a] .> 0) || error("Boundary cell widths must be finite and positive.")
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
    return SN_Boundary_Flux(Int64(dimension),deepcopy(Ω),copy(w),widths,copy(Eb),copy(boundary_conditions),faces,status)
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
