"""
    Cylindrical_Shell_Path_Result

Exact two-dimensional ray path through a cylindrical shell, used to screen when a locally planar
HTS-tape approximation becomes unreliable at grazing incidence. The particle starts on the outer
surface and points inward. `normal_cosine` is the cosine with the local inward radial normal.

This is an analytic geometry/heating benchmark, not a replacement for a curved OpenSn or Monte
Carlo transport calculation.
"""
struct Cylindrical_Shell_Path_Result
    midplane_radius_cm::Float64
    thickness_cm::Float64
    normal_cosine::Float64
    critical_through_thickness_cosine::Float64
    reaches_inner_surface::Bool
    exact_path_cm::Float64
    planar_path_cm::Float64
    path_relative_difference::Float64
    endpoint_radius_cm::Float64
    polar_displacement_rad::Float64
    midplane_arc_displacement_cm::Float64
end

function cylindrical_shell_critical_cosine(
    midplane_radius_cm::Real,
    thickness_cm::Real,
)
    radius = Float64(midplane_radius_cm)
    thickness = Float64(thickness_cm)
    isfinite(radius) && radius > 0.0 || error("Midplane radius must be finite and positive.")
    isfinite(thickness) && thickness > 0.0 || error("Shell thickness must be finite and positive.")
    radius > 0.5*thickness || error(
        "Shell midplane radius must exceed half the thickness.",
    )
    outer = radius+0.5*thickness
    inner = radius-0.5*thickness
    return sqrt(max(0.0,1.0-(inner/outer)^2))
end

"""
    cylindrical_shell_path(midplane_radius_cm, thickness_cm, normal_cosine)

Return the exact distance to the next shell boundary. If the ray reaches the inner cylindrical
surface, the smaller quadratic root is used. At more grazing incidence the ray misses the inner
surface and leaves again through the outer surface; the exact shell path is then the outer-circle
chord. The planar comparator is always `thickness/normal_cosine` and therefore exposes the topology
change rather than hiding it.
"""
function cylindrical_shell_path(
    midplane_radius_cm::Real,
    thickness_cm::Real,
    normal_cosine::Real;
    discriminant_tolerance::Real=1.0e-14,
)
    radius = Float64(midplane_radius_cm)
    thickness = Float64(thickness_cm)
    mu = Float64(normal_cosine)
    isfinite(radius) && radius > 0.0 || error("Midplane radius must be finite and positive.")
    isfinite(thickness) && thickness > 0.0 || error("Shell thickness must be finite and positive.")
    radius > 0.5*thickness || error(
        "Shell midplane radius must exceed half the thickness.",
    )
    isfinite(mu) && 0.0 < mu <= 1.0 || error(
        "Inward normal cosine must lie in (0,1].",
    )
    tolerance = Float64(discriminant_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 || error(
        "Discriminant tolerance must be finite and nonnegative.",
    )

    outer = radius+0.5*thickness
    inner = radius-0.5*thickness
    tangential_cosine = sqrt(max(0.0,1.0-mu^2))
    discriminant = inner^2-outer^2*(1.0-mu^2)
    scale = max(inner^2,outer^2,1.0)
    reaches_inner = discriminant >= -tolerance*scale
    exact_path = if reaches_inner
        outer*mu-sqrt(max(0.0,discriminant))
    else
        2.0*outer*mu
    end
    exact_path >= 0.0 || error("Calculated cylindrical-shell path is negative.")
    planar_path = thickness/mu
    relative_difference = abs(planar_path-exact_path)/max(exact_path,eps(Float64))

    endpoint_x = outer-mu*exact_path
    endpoint_y = tangential_cosine*exact_path
    endpoint_radius = hypot(endpoint_x,endpoint_y)
    polar_displacement = abs(atan(endpoint_y,endpoint_x))
    arc_displacement = radius*polar_displacement
    critical_cosine = cylindrical_shell_critical_cosine(radius,thickness)

    expected_radius = reaches_inner ? inner : outer
    isapprox(endpoint_radius,expected_radius;rtol=1.0e-11,atol=1.0e-13) || error(
        "Cylindrical-shell ray endpoint does not lie on the expected boundary.",
    )
    return Cylindrical_Shell_Path_Result(
        radius,thickness,mu,critical_cosine,reaches_inner,exact_path,planar_path,
        relative_difference,endpoint_radius,polar_displacement,arc_displacement,
    )
end

pure_absorption_fraction(macroscopic_absorption_cm_inv::Real,path_cm::Real) = begin
    sigma = Float64(macroscopic_absorption_cm_inv)
    path = Float64(path_cm)
    isfinite(sigma) && sigma >= 0.0 || error(
        "Macroscopic absorption must be finite and nonnegative.",
    )
    isfinite(path) && path >= 0.0 || error("Transport path must be finite and nonnegative.")
    -expm1(-sigma*path)
end

struct Cylindrical_Shell_Heating_Result
    path::Cylindrical_Shell_Path_Result
    macroscopic_absorption_cm_inv::Float64
    exact_absorbed_fraction::Float64
    planar_absorbed_fraction::Float64
    absorbed_fraction_relative_difference::Float64
    independent_patch_local::Union{Nothing,Bool}
    patch_length_cm::Union{Nothing,Float64}
end

function cylindrical_shell_heating(
    midplane_radius_cm::Real,
    thickness_cm::Real,
    normal_cosine::Real,
    macroscopic_absorption_cm_inv::Real;
    patch_length_cm::Union{Nothing,Real}=nothing,
)
    path = cylindrical_shell_path(midplane_radius_cm,thickness_cm,normal_cosine)
    sigma = Float64(macroscopic_absorption_cm_inv)
    exact = pure_absorption_fraction(sigma,path.exact_path_cm)
    planar = pure_absorption_fraction(sigma,path.planar_path_cm)
    relative = abs(planar-exact)/max(abs(exact),eps(Float64))
    patch_length = isnothing(patch_length_cm) ? nothing : Float64(patch_length_cm)
    locality = if isnothing(patch_length)
        nothing
    else
        isfinite(patch_length) && patch_length > 0.0 || error(
            "Patch length must be finite and positive.",
        )
        # A source entering at a patch centre remains within that independently solved patch only
        # when its projected curved-surface displacement is no greater than half the patch length.
        path.midplane_arc_displacement_cm <= 0.5*patch_length
    end
    return Cylindrical_Shell_Heating_Result(
        path,sigma,exact,planar,relative,locality,patch_length,
    )
end

struct Curvature_Heating_Screen
    results::Vector{Cylindrical_Shell_Heating_Result}
    maximum_path_relative_difference::Float64
    maximum_absorbed_fraction_relative_difference::Float64
    all_rays_reach_inner_surface::Bool
    all_rays_patch_local::Union{Nothing,Bool}
    planar_response_pass::Bool
    independent_patch_pass::Union{Nothing,Bool}
    response_relative_tolerance::Float64
    metadata::Dict{String,String}
end

"""
    screen_cylindrical_shell_heating(radius, thickness, sigma, normal_cosines; ...)

Screen a set of angular bins against the exact cylindrical-shell path and pure-absorption response.
A planar response passes only when every ray reaches the inner surface and every absorbed-fraction
error is within tolerance. If `patch_length_cm` is supplied, independent-patch locality is checked
separately. Neither result is a physical curved-transport qualification.
"""
function screen_cylindrical_shell_heating(
    midplane_radius_cm::Real,
    thickness_cm::Real,
    macroscopic_absorption_cm_inv::Real,
    normal_cosines::AbstractVector{<:Real};
    patch_length_cm::Union{Nothing,Real}=nothing,
    response_relative_tolerance::Real=0.01,
    metadata::AbstractDict=Dict{String,String}(),
)
    isempty(normal_cosines) && error("Curvature heating screen requires angular samples.")
    tolerance = Float64(response_relative_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 || error(
        "Response tolerance must be finite and nonnegative.",
    )
    results = Cylindrical_Shell_Heating_Result[
        cylindrical_shell_heating(
            midplane_radius_cm,thickness_cm,mu,macroscopic_absorption_cm_inv;
            patch_length_cm=patch_length_cm,
        ) for mu in normal_cosines
    ]
    all_inner = all(result.path.reaches_inner_surface for result in results)
    all_local = isnothing(patch_length_cm) ? nothing :
        all(result.independent_patch_local === true for result in results)
    max_path = maximum(result.path.path_relative_difference for result in results)
    max_heating = maximum(
        result.absorbed_fraction_relative_difference for result in results
    )
    planar_pass = all_inner && max_heating <= tolerance
    patch_pass = isnothing(all_local) ? nothing : all_local
    metadata_string = Dict{String,String}(
        "schema" => "radiant.hts.analytic_curvature_heating/v1",
        "classification" => "analytic-screening-only",
        "geometry" => "cylindrical-shell-cross-section",
        "transport" => "uncollided-pure-absorption",
        "physical_curved_reference_present" => "false",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Curvature_Heating_Screen(
        results,max_path,max_heating,all_inner,all_local,planar_pass,patch_pass,
        tolerance,metadata_string,
    )
end

"""Exact chord length corresponding to a maximum circle sagitta."""
function maximum_facet_chord_for_sagitta(
    radius_cm::Real,
    maximum_sagitta_cm::Real,
)
    radius = Float64(radius_cm)
    sagitta = Float64(maximum_sagitta_cm)
    isfinite(radius) && radius > 0.0 || error("Curvature radius must be positive.")
    isfinite(sagitta) && 0.0 <= sagitta <= radius || error(
        "Sagitta must lie between zero and the radius.",
    )
    return 2.0*sqrt(max(0.0,2.0*radius*sagitta-sagitta^2))
end

function curvature_heating_receipt(screen::Curvature_Heating_Screen)
    return Dict{String,Any}(
        "schema" => screen.metadata["schema"],
        "classification" => screen.metadata["classification"],
        "normal_cosines" => [result.path.normal_cosine for result in screen.results],
        "critical_through_thickness_cosine" =>
            screen.results[1].path.critical_through_thickness_cosine,
        "exact_path_cm" => [result.path.exact_path_cm for result in screen.results],
        "planar_path_cm" => [result.path.planar_path_cm for result in screen.results],
        "exact_absorbed_fraction" =>
            [result.exact_absorbed_fraction for result in screen.results],
        "planar_absorbed_fraction" =>
            [result.planar_absorbed_fraction for result in screen.results],
        "maximum_path_relative_difference" => screen.maximum_path_relative_difference,
        "maximum_absorbed_fraction_relative_difference" =>
            screen.maximum_absorbed_fraction_relative_difference,
        "all_rays_reach_inner_surface" => screen.all_rays_reach_inner_surface,
        "all_rays_patch_local" => isnothing(screen.all_rays_patch_local) ?
            "not-evaluated" : screen.all_rays_patch_local,
        "planar_response_pass" => screen.planar_response_pass,
        "independent_patch_pass" => isnothing(screen.independent_patch_pass) ?
            "not-evaluated" : screen.independent_patch_pass,
        "physical_curved_transport_pass" => false,
        "metadata" => copy(screen.metadata),
    )
end
