using Radiant
using Test

@testset "HTS extended qualification contracts" begin
    @testset "Cross_Sections constructor and analytic custom library" begin
        photon = Photon()
        material = Material("analytic-medium")
        Radiant.set_density(material,1.0)

        cross_sections = Cross_Sections()
        @test eltype(cross_sections.particles) == Particle
        @test isempty(cross_sections.custom_absorption)
        @test isempty(cross_sections.custom_scattering)

        Radiant.set_source(cross_sections,"custom")
        Radiant.set_materials(cross_sections,material)
        Radiant.set_particles(cross_sections,photon)
        Radiant.set_absorption(cross_sections,[0.2])
        Radiant.set_scattering(cross_sections,[0.3])
        Radiant.set_energy_boundaries(cross_sections,[[1.0,0.001]])
        Radiant.build(cross_sections)

        @test cross_sections.is_build
        @test Radiant.get_number_of_particles(cross_sections) == 1
        @test Radiant.get_number_of_materials(cross_sections) == 1
        @test Radiant.get_number_of_groups(cross_sections,photon) == 1
        @test Radiant.get_total(cross_sections,photon)[1,1] ≈ 0.5
        @test Radiant.get_absorption(cross_sections,photon)[1,1] ≈ 0.2
        @test Radiant.get_scattering(cross_sections,photon,photon,0)[1,1,1,1] ≈ 0.3
        @test Radiant.get_energy_width(cross_sections,photon)[1] ≈ 0.999

        first_library = cross_sections.multigroup_cross_sections
        Radiant.build(cross_sections)
        @test cross_sections.multigroup_cross_sections === first_library

        Radiant.set_energy(cross_sections,2.0)
        Radiant.set_cutoff(cross_sections,0.01)
        Radiant.set_number_of_groups(cross_sections,1)
        @test cross_sections.energy == [2.0]
        @test cross_sections.cutoff == [0.01]
        @test cross_sections.number_of_groups == [1]
    end

    @testset "Mechanism-resolved energy partition" begin
        partition = Energy_Partition(
            label="YBCO-event",
            available_energy_eV=1000.0,
            ionization_eV=300.0,
            electronic_excitation_eV=150.0,
            prompt_lattice_heat_eV=250.0,
            nuclear_recoil_handoff_eV=100.0,
            defect_stored_eV=50.0,
            optical_emission_eV=25.0,
            escaping_particle_eV=75.0,
            cutoff_handoff_eV=50.0,
        )
        @test get_accounted_energy(partition) ≈ 1000.0
        @test get_energy_partition_residual(partition) ≈ 0.0
        @test is_energy_partition_closed(partition)
        @test is_energy_partition_resolved(partition)
        @test get_prompt_heat_fraction(partition) ≈ 0.25

        unresolved = Energy_Partition(
            available_energy_eV=10.0,
            prompt_lattice_heat_eV=8.0,
            unresolved_eV=2.0,
        )
        @test is_energy_partition_closed(unresolved)
        @test !is_energy_partition_resolved(unresolved)
        @test_throws ErrorException Energy_Partition(
            available_energy_eV=1.0,
            ionization_eV=-1.0,
        )
    end

    @testset "Fibre diagnostic radial geometry" begin
        layers = Fibre_Radial_Layer[
            Fibre_Radial_Layer("core","silica-core",0.0,0.0025;optical_role=:core),
            Fibre_Radial_Layer("cladding","silica-cladding",0.0025,0.00625;optical_role=:cladding),
            Fibre_Radial_Layer("coating","polyimide",0.00625,0.00775;optical_role=:coating),
            Fibre_Radial_Layer("adhesive","cryogenic-adhesive",0.00775,0.009;optical_role=:adhesive),
            Fibre_Radial_Layer("capillary","steel-capillary",0.009,0.012;optical_role=:capillary),
        ]
        fibre = Fibre_Diagnostic_Definition(
            "HTS winding FBG",:fbg,layers;
            operating_temperature_K=77.0,
            interrogation_wavelength_nm=1550.0,
            responses=Symbol[
                :ionizing_dose,
                :capture_rate,
                :secondary_electron_spectrum,
                :radiation_induced_attenuation_state,
                :fbg_wavelength_shift,
            ],
        )
        @test get_outer_radius(fibre) ≈ 0.012
        @test get_cross_sectional_area(fibre) ≈ π*0.012^2
        @test get_fibre_layer_index(fibre,0.001) == 1
        @test get_fibre_layer_index(fibre,0.010) == 5
        @test requires_optical_response_model(fibre)

        bad_layers = Fibre_Radial_Layer[
            Fibre_Radial_Layer("core","silica",0.0,0.002),
            Fibre_Radial_Layer("gap","polymer",0.003,0.004),
        ]
        @test_throws ErrorException Fibre_Diagnostic_Definition(
            "nonconformal",:dosimetry,bad_layers,
        )
    end

    @testset "HTS detector transport and electrothermal separation" begin
        fixture = verification_b10_ybco_detector()
        @test fixture.mode == :transition_edge
        @test fixture.active_material_tag == "YBCO"
        @test fixture.operating_temperature_K == 77.0
        @test length(fixture.converter_layers) == 1
        @test fixture.converter_layers[1].reaction_channel == "B-10(n,alpha)Li-7"
        @test fixture.metadata["classification"] == "verification-only"
        @test get_active_volume(fixture) ≈ 1.0e-7
        @test !is_ready_for_transient_response(fixture)

        thermal = Detector_Thermal_Model(
            active_heat_capacity_J_K=1.0e-9,
            substrate_heat_capacity_J_K=1.0e-6,
            active_substrate_conductance_W_K=1.0e-5,
            substrate_bath_conductance_W_K=1.0e-4,
            transition_temperature_K=90.0,
            transition_width_K=1.0,
            normal_resistance_ohm=10.0,
            inductance_H=1.0e-6,
        )
        detector = HTS_Detector_Definition(
            "qualified-schema-detector",:transition_edge,"YBCO";
            operating_temperature_K=80.0,
            bias_current_A=1.0e-3,
            active_width_cm=0.1,
            active_length_cm=0.1,
            active_thickness_cm=1.0e-5,
            substrate_material_tag="sapphire",
            converter_layers=fixture.converter_layers,
            thermal_model=thermal,
        )
        @test is_ready_for_transient_response(detector)
    end

    @testset "Physics ownership and protected responses" begin
        coverage = default_hts_physics_coverage()
        @test validate_physics_coverage(coverage)
        @test get_physics_record(
            coverage;
            domain="tape-and-diagnostic-microdomain",
            process_id="photoatomic-coupled-em",
        ).production_owner == "Radiant"
        @test get_physics_record(
            coverage;
            domain="tape-and-diagnostic-microdomain",
            process_id="photonuclear-reactions",
        ).status == :external_production
        @test get_physics_record(
            coverage;
            domain="selected-interface-microdomain",
            process_id="sub-keV-electron-track-structure",
        ).energy_max_eV == 1.0e3

        registry = default_hts_protected_responses()
        ybco = get_protected_response(registry,"ybco_em_deposition")
        @test ybco.production_owner == "Radiant"
        @test response_is_converged(ybco,100.0,104.0)
        @test !response_is_converged(ybco,100.0,106.0)
        @test get_protected_response(
            registry,"fibre_core_ionizing_dose",
        ).reference_solver == "Geant4 matched-physics"
    end
end
