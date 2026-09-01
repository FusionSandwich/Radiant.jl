using Radiant
using Test
using LinearAlgebra

@testset "Mapped structured curvilinear metrics" begin
    cartesian = mapped_cartesian_grid([0.0,1.0,2.0],[0.0,1.0],[0.0,2.0])
    @test assert_mapped_grid_invariants(cartesian)
    @test cartesian.cell_volumes_cm3 == fill(2.0,2,1,1)
    @test cartesian.jacobian_determinants == ones(2,1,1)
    @test all(cartesian.face_areas_cm2 .> 0.0)

    cylindrical = mapped_cylindrical_grid(
        [1.0,1.5,2.0],collect(range(0.0,pi/2.0,length=5)),[0.0,1.0,2.0],
    )
    @test mapped_grid_invariants(cylindrical)["passed"]
    @test cylindrical.coordinate_system == :cylindrical

    toroidal = mapped_toroidal_grid(
        [0.2,0.3],collect(range(0.0,2.0*pi,length=9)),
        collect(range(0.0,pi/2.0,length=5));major_radius_cm=5.0,
    )
    @test assert_mapped_grid_invariants(toroidal)

    helical = mapped_helical_grid(
        [-0.05,0.0,0.05],[-0.02,0.02],collect(range(0.0,2.0*pi,length=17));
        radius_cm=2.0,pitch_per_turn_cm=0.5,
    )
    @test assert_mapped_grid_invariants(helical)
    repeated_turns = mapped_helical_grid(
        [-0.05,0.05],[-0.02,0.02],collect(range(0.0,4.0*pi,length=33));
        radius_cm=2.0,pitch_per_turn_cm=0.5,coordinate_system=:repeated_turns,
    )
    @test mapped_grid_invariants(repeated_turns)["passed"]
    @test sum(repeated_turns.cell_volumes_cm3) ≈
        2.0*sum(mapped_helical_grid(
            [-0.05,0.05],[-0.02,0.02],collect(range(0.0,2.0*pi,length=17));
            radius_cm=2.0,pitch_per_turn_cm=0.5,
        ).cell_volumes_cm3) rtol=1.0e-10

    radius = 3.0
    pitch = 0.4
    tangent(angle) = normalize([-radius*sin(angle),radius*cos(angle),pitch])
    normal(angle) = [cos(angle),sin(angle),0.0]
    binormal(angle) = cross(tangent(angle),normal(angle))
    centerline(angle) = [
        radius*cos(angle),radius*sin(angle),pitch*angle+0.02*sin(2.0*angle),
    ]
    varying = mapped_frenet_grid(
        [-0.02,0.02],[-0.01,0.01],collect(range(0.0,pi,length=17)),
        centerline,normal,binormal,
    )
    @test assert_mapped_grid_invariants(varying)
    @test varying.coordinate_system == :varying_curvature_torsion

    fibre = mapped_fibre_grid(
        [0.1,0.2],collect(range(0.0,2.0*pi,length=17)),[0.0,1.0],
    )
    @test assert_mapped_grid_invariants(fibre)

    @test_throws ErrorException mapped_cartesian_grid([0.0,0.0],[0.0,1.0],[0.0,1.0])
    @test_throws ErrorException mapped_toroidal_grid(
        [1.0,2.0],[0.0,1.0],[0.0,1.0];major_radius_cm=1.5,
    )
end

@testset "Cached atlas and multiscale conservative transfer" begin
    atlas = build_piecewise_flat_tape_atlas(
        [0.0 0.0 0.0; 1.0 0.0 0.0; 2.0 0.1 0.0];
        width_cm=0.4,thickness_cm=0.01,geometry_hash="cache-fixture",
    )
    cache = build_atlas_mapping_cache(atlas)
    local_point = [0.1,0.02,0.001]
    global_point = atlas_local_to_global(cache,1,local_point)
    @test atlas_global_to_local(cache,1,global_point) ≈ local_point atol=1.0e-12
    @test cache.atlas_hash == atlas.geometry_hash

    source_values = [1.0,2.0,3.0]
    for scale in (1.0e3,1.0e5,1.0e7)
        fixture = multiscale_local_layer_fixture(scale)
        @test fixture.layer_ids == [
            "film-0.1um","film-0.25um","film-0.5um","film-1.0um",
            "buffer-0.3um","metallic-cap-10um","substrate-50um",
            "cable-scale","macro-block",
        ]
        @test fixture.layer_thickness_cm[1:7] ==
            1.0e-5.*[1.0,2.5,5.0,10.0,3.0,100.0,500.0]
        @test all(fixture.layer_thickness_cm .> 0.0)
        @test maximum(fixture.layer_thickness_cm)/minimum(fixture.layer_thickness_cm) ≈
            scale rtol=8.0*eps(Float64) atol=0.0
        @test fixture.nondimensional_edges == sort(fixture.nondimensional_edges)
        @test fixture.grid.cell_volumes_cm3[:,1,1] ≈ fixture.layer_thickness_cm
        @test fixture.jacobian_condition_number <= 1.0e5+1.0
        for physical in cumsum(fixture.layer_thickness_cm)
            @test physical_layer_coordinate(
                fixture,local_layer_coordinate(fixture,physical),
            ) ≈ physical rtol=0.0 atol=eps(physical)*2.0
        end
        for layer_index in eachindex(fixture.layer_ids)
            origin = fixture.block_local_origins_cm[layer_index]
            scale_cm = fixture.block_length_scales_cm[layer_index]
            @test local_layer_coordinate(fixture,layer_index,origin) == 0.0
            @test physical_layer_coordinate(fixture,layer_index,1.0) ≈ origin+scale_cm
        end
        if scale == 1.0e7
            @test fixture.layer_thickness_cm[end-1] == 0.1 # 1 mm cable
            @test fixture.layer_thickness_cm[end] ≈ 100.0 rtol=8.0*eps(Float64) # 1 m macro block
        end
        transfer = conservative_interface_transfer(
            scale.*[0.0,0.2,0.6,1.0],scale.*[0.0,0.1,0.5,0.8,1.0],
        )
        target = apply_interface_transfer(transfer,source_values)
        @test sum(target) ≈ sum(source_values) atol=1.0e-12
        @test all(target .>= 0.0)
    end
    @test_throws ErrorException conservative_interface_transfer(
        [0.0,1.0],[0.0,0.5],
    )
end

@testset "Conservative mapped-coordinate streaming" begin
    grid = mapped_cartesian_grid([0.0,1.0,2.0,3.0],[0.0,1.0],[0.0,1.0])
    flux = reshape([1.0,2.0,4.0],3,1,1)
    result = mapped_streaming_divergence(grid,[1.0,0.0,0.0],flux)
    @test result.angular_divergence[:,1,1] ≈ [1.0,1.0,2.0]
    @test result.integrated_divergence ≈ 4.0
    @test result.signed_boundary_current ≈ 4.0
    @test result.particle_balance_residual ≈ 0.0 atol=1.0e-12
    @test result.unique_face_count == result.expected_unique_face_count == 16
    @test result.maximum_internal_face_metric_mismatch <= 1.0e-12
    closure = mapped_streaming_closure(result;particle_energy_MeV=500.0)
    @test closure["particle_balance_residual"] ≈ 0.0 atol=1.0e-12
    @test closure["energy_current_balance_residual_MeV"] ≈ 0.0 atol=1.0e-12

    constant_result = mapped_streaming_divergence(
        grid,[1.0,0.0,0.0],fill(2.0,3,1,1);
        boundary_inflow=(axis,side,i,j,k) -> 2.0,
    )
    @test maximum(abs.(constant_result.angular_divergence)) <= 1.0e-12
    @test constant_result.signed_boundary_current ≈ 0.0 atol=1.0e-12

    rotation = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
    rotated = build_mapped_structured_grid(
        ([0.0,1.0,2.0,3.0],[0.0,1.0],[0.0,1.0]),
        (x,y,z) -> rotation*[x,y,z];coordinate_system=:rotated_cartesian,
        geometry_hash="rigid-rotation-fixture",
    )
    rotated_result = mapped_streaming_divergence(
        rotated,transform_polar_vector(rotation,[1.0,0.0,0.0]),flux,
    )
    @test rotated_result.angular_divergence ≈ result.angular_divergence atol=1.0e-12

    reflection = [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test transform_polar_vector(reflection,[1.0,2.0,3.0]) == [-1.0,2.0,3.0]
    @test transform_axial_vector(reflection,[0.0,0.0,2.0]) == [0.0,0.0,-2.0]
    species = proton_species()
    original_direction,_ = magnetic_direction_step(
        species,500.0,[1.0,0.0,0.0],[0.0,0.0,2.0],1.0,
    )
    reflected_direction,_ = magnetic_direction_step(
        species,500.0,transform_polar_vector(reflection,[1.0,0.0,0.0]),
        transform_axial_vector(reflection,[0.0,0.0,2.0]),1.0,
    )
    @test collect(reflected_direction) ≈
        transform_polar_vector(reflection,collect(original_direction)) atol=1.0e-12
end

@testset "Preregistered curved performance phases" begin
    target = Curved_Performance_Target(
        target_id="dpa-b653f509-curved-performance-target",maximum_cells=10_000,
    )
    result = benchmark_curved_pipeline(
        target,
        () -> begin
            baseline = mapped_cartesian_grid(
                [-0.02,0.0,0.02],[-0.01,0.01],collect(range(0.0,pi,length=33)),
            )
            improved = mapped_helical_grid(
                [-0.02,0.0,0.02],[-0.01,0.01],collect(range(0.0,pi,length=33));
                radius_cm=2.0,pitch_per_turn_cm=0.2,
            )
            (baseline=baseline,improved=improved)
        end;
        baseline_sweep=grid -> mapped_streaming_divergence(
            grid,[0.0,0.0,1.0],ones(size(grid.cell_volumes_cm3));
            boundary_inflow=(axis,side,i,j,k) -> 1.0,
        ),
        improved_sweep=grid -> mapped_streaming_divergence(
            grid,[0.0,0.0,1.0],ones(size(grid.cell_volumes_cm3));
            boundary_inflow=(axis,side,i,j,k) -> 1.0,
        ),
        scoring=(result,grid) -> sum(abs,result.angular_divergence),
        baseline_remap=(result,grid) ->
            reverse(reverse(result.angular_divergence,dims=3),dims=3),
        cached_atlas_remap=(result,grid) -> reverse(result.angular_divergence,dims=3),
    )
    @test result["target_preregistered_before_timing"]
    @test result["passed"]
    @test result["cell_count"] == 64
    @test all(result[key] >= 0.0 for key in ("preprocessing_s","warmup_s","scoring_s"))
    @test result["ratios"]["equal_error_solve"] === nothing
    @test result["ratios"]["matched_unknown_kernel_time"] >= 0.0
    @test result["measurement_status"]["equal_error_solve"] ==
        "NOT_MEASURED_DIAGNOSTIC_ONLY"
    @test result["measurement_status"]["matched_unknown_kernel_time"] ==
        "MEASURED_DIAGNOSTIC_ONLY"
    @test result["thresholds"] == Dict{String,Any}(
        "sweep_per_unknown" => 1.25,"equal_error_solve" => 1.35,
        "memory" => 1.25,"cached_atlas" => 1.50,
    )
end

@testset "Species-generic proton and ion transport primitives" begin
    species,source,model = synthetic_proton_transport_fixture()
    zero = relativistic_ion_kinematics(species,0.0)
    high = relativistic_ion_kinematics(species,500.0)
    @test zero.beta == 0.0
    @test high.total_energy_MeV == species.rest_mass_MeV+500.0
    @test high.momentum_MeV_c^2 ≈
        high.total_energy_MeV^2-species.rest_mass_MeV^2 rtol=1.0e-12
    @test 0.0 < high.beta < 1.0
    @test high.speed_m_s < 299_792_458.0
    @test_throws ErrorException relativistic_ion_kinematics(species,500.1)

    particle = particle_from_ion(species)
    @test Radiant.get_tag(particle) == "proton"
    @test Radiant.get_mass(particle) == species.rest_mass_MeV
    @test Radiant.get_charge(particle) == 1.0

    coefficients = ion_transport_coefficients(model,species,250.0)
    @test coefficients.electronic_stopping_MeV_cm == 1.0
    @test coefficients.nuclear_stopping_MeV_cm == 0.1
    route = Nonelastic_Secondary_Route(
        channel_id="p-nonelastic",production_owner="evaluated-ion-producer",
        handoff_schema="weighted_nonelastic_secondary_bank/v1",
        differential_data_hash="synthetic-differential-fixture",status=:candidate,
    )
    step = ion_transport_step(
        source,model,2.0;magnetic_field_T=[0.0,0.0,2.0],nonelastic_route=route,
    )
    @test step.final_energy_MeV ≈ 497.8
    @test step.electronic_loss_MeV == 2.0
    @test step.nuclear_loss_MeV == 0.2
    @test step.energy_variance_MeV2 == 0.02
    @test step.angular_variance_rad2 == 2.0e-4
    @test norm(collect(step.final_direction)) ≈ 1.0 atol=1.0e-12
    @test step.final_direction != step.initial_direction
    @test step.nonelastic_route === route

    @test_throws ErrorException ion_transport_step(source,model,500.0)
    @test_throws ErrorException ion_transport_coefficients(model,species,501.0)
    @test_throws ErrorException Nonelastic_Secondary_Route(
        channel_id="p-nonelastic",production_owner="",handoff_schema="bank/v1",
        differential_data_hash="missing",status=:candidate,
    )
    @test_throws ErrorException Tabulated_Ion_Transport_Model(
        species_id="proton",material_id="synthetic",energy_MeV=[0.0,600.0],
        electronic_stopping_MeV_cm=[1.0,1.0],nuclear_stopping_MeV_cm=[0.0,0.0],
        energy_straggling_variance_MeV2_cm=[0.0,0.0],angular_variance_rad2_cm=[0.0,0.0],
        data_hash="bad-range",
    )
end
