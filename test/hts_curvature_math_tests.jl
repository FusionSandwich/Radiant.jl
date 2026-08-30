using Radiant
using Test

@testset "Analytic cylindrical-shell curvature and grazing heating" begin
    radius = 10.0
    thickness = 0.01

    normal = cylindrical_shell_path(radius,thickness,1.0)
    @test normal.reaches_inner_surface
    @test normal.exact_path_cm ≈ thickness rtol=1.0e-12 atol=1.0e-14
    @test normal.planar_path_cm ≈ thickness
    @test normal.path_relative_difference ≤ 1.0e-12
    @test normal.midplane_arc_displacement_cm ≈ 0.0 atol=1.0e-14

    critical = cylindrical_shell_critical_cosine(radius,thickness)
    @test 0.0 < critical < 1.0
    through = cylindrical_shell_path(radius,thickness,min(1.0,1.01*critical))
    miss = cylindrical_shell_path(radius,thickness,0.99*critical)
    @test through.reaches_inner_surface
    @test !miss.reaches_inner_surface
    @test miss.endpoint_radius_cm ≈ radius+0.5*thickness rtol=1.0e-11
    # Below the tangency threshold the ray remains in the shell until it re-crosses the outer
    # surface; the curved chord is therefore longer than the local planar thickness/mu path.
    @test miss.exact_path_cm > miss.planar_path_cm

    sigma = 0.3
    normal_heating = cylindrical_shell_heating(radius,thickness,1.0,sigma)
    @test normal_heating.exact_absorbed_fraction ≈
          pure_absorption_fraction(sigma,thickness) rtol=1.0e-12
    @test normal_heating.planar_absorbed_fraction ≈
          normal_heating.exact_absorbed_fraction rtol=1.0e-12

    high_mu_screen = screen_cylindrical_shell_heating(
        radius,thickness,sigma,[1.0,0.95,0.90];
        patch_length_cm=1.0,response_relative_tolerance=1.0e-3,
    )
    @test high_mu_screen.all_rays_reach_inner_surface
    @test high_mu_screen.planar_response_pass
    @test high_mu_screen.independent_patch_pass === true
    @test curvature_heating_receipt(high_mu_screen)["physical_curved_transport_pass"] == false

    grazing_screen = screen_cylindrical_shell_heating(
        radius,thickness,sigma,[1.0,0.99*critical];
        patch_length_cm=1.0,response_relative_tolerance=0.05,
    )
    @test !grazing_screen.all_rays_reach_inner_surface
    @test !grazing_screen.planar_response_pass
    @test grazing_screen.maximum_path_relative_difference > 0.0
    @test grazing_screen.maximum_absorbed_fraction_relative_difference > 0.0

    @test maximum_facet_chord_for_sagitta(10.0,0.0) == 0.0
    chord = maximum_facet_chord_for_sagitta(10.0,0.01)
    recovered_sagitta = 10.0-sqrt(10.0^2-(0.5*chord)^2)
    @test recovered_sagitta ≈ 0.01 rtol=1.0e-12 atol=1.0e-14

    @test_throws ErrorException cylindrical_shell_path(radius,thickness,0.0)
    @test_throws ErrorException cylindrical_shell_path(0.001,0.01,1.0)
    @test_throws ErrorException pure_absorption_fraction(-1.0,1.0)
end
