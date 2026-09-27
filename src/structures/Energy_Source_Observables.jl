"""
    get_volume_source_energy_MeV(source::Energy_Moment_Volume_Source; physical=false)

Return the represented first energy moment per ascending source energy group, summed over
source voxels and integrated over angle. This is a declared-source diagnostic; it does not
represent transported, deposited, escaped, or physically validated energy.
"""
function get_volume_source_energy_MeV(
    wrapper::Energy_Moment_Volume_Source;
    physical::Bool=false,
)
    source = wrapper.source
    moments = wrapper.energy_moments
    size(moments,4) in 1:4 || error("Energy source requires one to four energy modes.")
    all(isfinite,moments) || error("Energy source moments must be finite.")
    values = source.values
    all(isfinite,values) || error("Volume source values must be finite.")
    size(moments)[1:3] == size(values) || error("Energy moments and base source shapes differ.")
    moments[:,:,:,1] == values || error("Zeroth energy moments differ from the declared base source.")

    nvoxel, ngroup, nangular = size(values)
    length(source.voxel_volumes_cm3) == nvoxel || error("Source voxel volumes do not match source shape.")
    all(v -> isfinite(v) && v > 0.0,source.voxel_volumes_cm3) ||
        error("Source voxel volumes must be finite and positive.")
    length(source.energy_edges_eV) == ngroup+1 || error("Source energy edges do not match source groups.")
    edges = source.energy_edges_eV
    all(e -> isfinite(e) && e >= 0.0,edges) && all(diff(edges) .> 0.0) ||
        error("Source energy edges must be finite, nonnegative, and strictly ascending.")

    angular_indices = Int[]
    angular_weights = Float64[]
    if source.angular_representation == :isotropic
        nangular == 1 || error("Isotropic source must have one angle-integrated coefficient.")
        angular_indices = [1]
        angular_weights = [1.0]
    elseif source.angular_representation == :ordinates
        length(source.quadrature_weights) == nangular ||
            error("Ordinate source weights do not match angular coefficients.")
        size(source.directions) == (nangular,3) || error("Ordinate directions do not match angular coefficients.")
        all(w -> isfinite(w) && w > 0.0,source.quadrature_weights) ||
            error("Ordinate weights must be finite and positive.")
        all(isfinite,source.directions) || error("Ordinate directions must be finite.")
        angular_indices = collect(1:nangular)
        angular_weights = source.quadrature_weights
    elseif source.angular_representation == :moments
        get(source.provenance,"zeroth_moment_is_angle_integrated","false") == "true" ||
            error("Moment-source energy extraction requires zeroth_moment_is_angle_integrated=\"true\".")
        haskey(source.provenance,"zeroth_moment_index") ||
            error("Moment-source energy extraction requires zeroth_moment_index provenance.")
        index = try
            parse(Int,source.provenance["zeroth_moment_index"])
        catch
            error("zeroth_moment_index must be a valid one-based integer.")
        end
        1 <= index <= nangular || error("zeroth_moment_index lies outside the angular coefficient dimension.")
        angular_indices = [index]
        angular_weights = [1.0]
    else
        error("Unsupported volume-source angular representation.")
    end

    result = zeros(Float64,ngroup)
    for group in 1:ngroup
        lower = edges[group]
        upper = edges[group+1]
        midpoint = lower/2 + upper/2
        width = upper-lower
        isfinite(midpoint) && isfinite(width) || error("Source energy integration weights overflowed.")
        for voxel in 1:nvoxel, slot in eachindex(angular_indices)
            angular = angular_indices[slot]
            q0 = moments[voxel,group,angular,1]
            q1 = size(moments,4) >= 2 ? moments[voxel,group,angular,2] : 0.0
            first_moment_eV = midpoint*q0 + width*q1/(2sqrt(3.0))
            contribution = source.voxel_volumes_cm3[voxel] * angular_weights[slot] * first_moment_eV
            isfinite(first_moment_eV) && isfinite(contribution) ||
                error("Volume source energy moment overflowed.")
            result[group] += contribution
            isfinite(result[group]) || error("Accumulated volume source energy moment overflowed.")
        end
    end
    result .*= 1.0e-6 # eV to MeV
    all(isfinite,result) || error("Converted volume source energy moment overflowed.")
    if physical
        result = apply_normalization(result,get_source_normalization(wrapper))
        all(isfinite,result) || error("Physically scaled volume source energy moment overflowed.")
    end
    return result
end
