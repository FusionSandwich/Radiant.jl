using Radiant
using Test

@testset "Radiant.jl" begin

    @testset "Explicit source normalization" begin
        normalization = Source_Normalization(
            basis=:per_history,
            source_rate_per_s=2.0,
            symmetry_factor=4.0,
            time_interval_s=(1.0,6.0),
            time_class=:prompt,
            source_hash="analytic-fixture",
            provenance=Dict("producer" => "unit-test"),
        )
        @test get_physical_scale(normalization) == 8.0
        @test get_duration(normalization) == 5.0
        @test apply_normalization(3.0,normalization) == 24.0
        @test_throws ErrorException Source_Normalization(basis=:unknown)
        @test_throws ErrorException Source_Normalization(source_rate_per_s=0.0)
        @test_throws ErrorException Source_Normalization(time_interval_s=(2.0,1.0))
    end

    @testset "Current-conserving boundary source" begin
        particle = Photon()
        normalization = Source_Normalization(
            source_rate_per_s=2.0,
            symmetry_factor=4.0,
            source_hash="boundary-fixture",
        )
        patch_ids = [11]
        centroids = [0.0 0.0 0.0]
        areas = [2.0]
        normals = [1.0 0.0 0.0]
        tangent_1 = [0.0 1.0 0.0]
        tangent_2 = [0.0 0.0 1.0]
        energy_edges = [1.0,10.0,100.0]
        directions = [
            -1.0 0.0 0.0
            -0.5 sqrt(3.0)/2.0 0.0
        ]
        weights = [0.25,0.75]
        angular_flux = zeros(Float64,1,2,2)
        angular_flux[1,1,:] .= [4.0,8.0]
        angular_flux[1,2,:] .= [2.0,4.0]

        source = Boundary_Angular_Current_Source(
            particle,
            patch_ids,
            centroids,
            areas,
            normals,
            tangent_1,
            tangent_2,
            energy_edges,
            directions,
            weights,
            angular_flux,
            normalization,
        )

        @test get_incoming_current_density(source) ≈ [4.0 2.0]
        @test get_incoming_current(source) ≈ [8.0 4.0]
        @test get_total_incoming_current(source) ≈ 12.0
        @test get_total_incoming_current(source;physical=true) ≈ 96.0
        @test assert_current_closure(source,[8.0 4.0])
        @test_throws ErrorException assert_current_closure(source,[8.1 4.0])

        directional_current = zeros(Float64,1,2,2)
        directional_current[1,1,:] .= [1.0,3.0]
        directional_current[1,2,:] .= [0.5,1.5]
        projected = boundary_source_from_directional_current(
            particle,
            patch_ids,
            centroids,
            areas,
            normals,
            tangent_1,
            tangent_2,
            energy_edges,
            directions,
            weights,
            directional_current,
            normalization,
        )
        @test get_incoming_current(projected) ≈ [8.0 4.0]

        grazing_directions = [-1.0e-15 sqrt(1.0-1.0e-30) 0.0]
        grazing_current = ones(Float64,1,2,1)
        @test_throws ErrorException boundary_source_from_directional_current(
            particle,
            patch_ids,
            centroids,
            areas,
            normals,
            tangent_1,
            tangent_2,
            energy_edges,
            grazing_directions,
            [1.0],
            grazing_current,
            normalization,
        )

        bad_flux = copy(angular_flux)
        bad_flux[1,1,1] = -1.0
        @test_throws ErrorException Boundary_Angular_Current_Source(
            particle,
            patch_ids,
            centroids,
            areas,
            normals,
            tangent_1,
            tangent_2,
            energy_edges,
            directions,
            weights,
            bad_flux,
            normalization,
        )
    end

    @testset "Anisotropic volume-source schema" begin
        normalization = Source_Normalization(
            source_rate_per_s=2.0,
            symmetry_factor=4.0,
            source_hash="volume-fixture",
        )
        values = zeros(Float64,2,2,1)
        values[1,:,1] .= [1.0,2.0]
        values[2,:,1] .= [3.0,4.0]
        source = Anisotropic_Volume_Source(
            Photon(),
            [1,2],
            [1.0,2.0],
            [1.0,10.0,100.0],
            :isotropic,
            values,
            normalization;
            parent_reaction="delayed-photon",
        )
        @test get_volume_source_rate(source) ≈ [7.0,10.0]
        @test get_volume_source_rate(source;physical=true) ≈ [56.0,80.0]

        moment_values = ones(Float64,1,1,2)
        @test_throws ErrorException Anisotropic_Volume_Source(
            Photon(),
            [1],
            [1.0],
            [1.0,2.0],
            :moments,
            moment_values,
            normalization,
        )
    end

    @testset "HTS tape-stack definitions" begin
        stack = verification_eight_layer_stack()
        @test stack.width_cm == 0.4
        @test get_total_thickness(stack) ≈ 0.01632 atol=1.0e-14
        @test length(stack.layers) == 8
        @test stack.layers[3].material_tag == "YBCO"
        @test stack.layers[3].activation_enabled
        @test stack.layers[3].atomistic_enabled
        @test stack.metadata["classification"] == "verification-only"

        boundaries = get_layer_boundaries(stack)
        @test length(boundaries) == 9
        @test boundaries[end] ≈ 0.01632 atol=1.0e-14
        @test get_layer_index(stack,(boundaries[3]+boundaries[4])/2.0) == 3
        @test get_layer_index(stack,boundaries[end]) == 8
        @test_throws ErrorException get_layer_index(stack,-1.0e-6)

        gd_stack = verification_eight_layer_stack(rebco_material_tag="GdBCO")
        @test gd_stack.layers[3].material_tag == "GdBCO"
    end

    @testset "Transport ownership and no double counting" begin
        production_record = Transport_Ownership_Record(
            domain="tape-microdomain",
            particle="photon",
            source_class="prompt-external",
            response="heating",
            owner="Radiant",
            artifact_hash="fixture-a",
        )
        comparison_record = Transport_Ownership_Record(
            domain="tape-microdomain",
            particle="photon",
            source_class="prompt-external",
            response="heating",
            owner="OpenSn",
            comparison_only=true,
            artifact_hash="fixture-b",
        )
        production_map = Transport_Ownership_Map(
            [production_record,comparison_record];
            mode=:production,
        )
        @test validate_ownership(production_map)
        @test get_production_owner(
            production_map;
            domain="tape-microdomain",
            particle="photon",
            source_class="prompt-external",
            response="heating",
        ) == "Radiant"

        duplicate_record = Transport_Ownership_Record(
            domain="tape-microdomain",
            particle="photon",
            source_class="prompt-external",
            response="heating",
            owner="OpenSn",
            artifact_hash="fixture-c",
        )
        @test_throws ErrorException Transport_Ownership_Map(
            [production_record,duplicate_record];
            mode=:production,
        )

        comparison_map = Transport_Ownership_Map(
            [
                Transport_Ownership_Record(
                    domain="verification-stack",
                    particle="photon",
                    source_class="analytic",
                    response="heating",
                    owner="Radiant",
                    comparison_only=true,
                ),
                Transport_Ownership_Record(
                    domain="verification-stack",
                    particle="photon",
                    source_class="analytic",
                    response="heating",
                    owner="OpenSn",
                    comparison_only=true,
                ),
            ];
            mode=:comparison,
        )
        @test validate_ownership(comparison_map)
    end

    @testset "Explicit transport ledgers" begin
        balance = Transport_Balance(
            label="photon-electron",
            injected_particles=10.0,
            produced_particles=2.0,
            absorbed_particles=5.0,
            leaked_particles=6.0,
            cutoff_particles=1.0,
            injected_kinetic_energy_MeV=100.0,
            produced_kinetic_energy_MeV=10.0,
            deposited_energy_MeV=70.0,
            leaked_energy_MeV=30.0,
            cutoff_energy_MeV=10.0,
            injected_charge=1.0,
            deposited_charge=0.4,
            leaked_charge=0.6,
        )
        @test get_particle_residual(balance) == 0.0
        @test get_energy_residual(balance) == 0.0
        @test abs(get_charge_residual(balance)) < 1.0e-14
        @test is_balanced(balance)

        unbalanced = Transport_Balance(injected_particles=1.0)
        @test !is_balanced(unbalanced)
    end

end
nothing
