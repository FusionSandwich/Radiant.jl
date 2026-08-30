using Radiant
using Test

@testset "Parent-resolved activation delayed-source bundle" begin
    scheme = synthetic_activation_decay_fixture()
    activity = Delayed_Activity_Field(
        "X-100",[11,12],[1.0,2.0],[2.0,4.0];
        activity_standard_uncertainty_Bq=[0.2,0.4],
        material_tags=["YBCO","YBCO"],
        cooling_interval_s=(100.0,200.0),
        inventory_hash="synthetic-inventory-hash",
        activation_artifact_hash="synthetic-activation-artifact",
        provenance=Dict(
            "classification" => "software-verification",
            "scalar_flux_not_surface_current" => "true",
        ),
    )
    grids = Dict{Symbol,Vector{Float64}}(
        :photon => [0.0,1.0e5,1.0e6,2.0e6],
        :electron => [0.0,1.0e5,1.0e6,2.0e6],
        :positron => [0.0,1.0e5,1.0e6,2.0e6],
    )
    bundle = build_activation_delayed_source_bundle(
        scheme,activity;energy_grids=grids,
    )
    @test bundle.schema == "radiant.hts.activation_delayed_source_bundle/v1"
    @test bundle.parent_nuclide == "X-100"
    @test bundle.daughter_nuclide == "Y-100"
    @test Set(keys(bundle.em_sources)) == Set([:photon,:electron])
    @test bundle.total_activity_Bq == 6.0
    @test sum(get_volume_source_rate(bundle.em_sources[:photon])) ≈ 6.0
    @test sum(get_volume_source_rate(bundle.em_sources[:electron])) ≈ 6.0
    @test bundle.em_sources[:photon].normalization.time_class == :delayed
    @test bundle.em_sources[:photon].normalization.basis == :per_second
    @test bundle.em_sources[:photon].normalization.time_interval_s == (100.0,200.0)
    @test bundle.em_sources[:photon].provenance["source_energy_is_not_deposited_heat"] ==
          "true"

    recoil = only([handoff for handoff in bundle.non_em_handoffs
                   if handoff.species == :recoil])
    neutrino = only([handoff for handoff in bundle.non_em_handoffs
                     if handoff.species == :neutrino])
    @test sum(recoil.event_rate_per_s) == 6.0
    @test sum(neutrino.event_rate_per_s) == 6.0
    @test recoil.transport_owner == "SPECTRA-PKA/BCA/MD"
    @test neutrino.transport_owner == "escape-energy-ledger"

    @test bundle.q_value_energy_rate_eV_per_s == 18.0e6
    @test bundle.em_source_energy_rate_eV_per_s == 9.0e6
    @test bundle.non_em_energy_rate_eV_per_s == 0.6e6
    @test bundle.neutrino_energy_rate_eV_per_s == 8.4e6
    @test bundle.unresolved_energy_rate_eV_per_s == 0.0
    @test abs(bundle.energy_balance_residual_eV_per_s) <= 1.0e-8
    @test !activation_delayed_bundle_is_production_ready(bundle)

    receipt = activation_delayed_source_receipt(bundle)
    @test receipt["source_energy_is_not_deposited_heat"]
    @test !receipt["neutrino_energy_is_local_heat"]
    @test !receipt["production_ready"]
    @test receipt["transported_em_species"] == ["electron","photon"]
    @test receipt["non_em_handoff_species"] == ["neutrino","recoil"]
end

@testset "Delayed source uncertainty and group placement" begin
    scheme = synthetic_activation_decay_fixture()
    activity = Delayed_Activity_Field(
        "X-100",[1],[2.0],[4.0];
        activity_standard_uncertainty_Bq=[0.4],
        material_tags=["REBCO"],
        cooling_interval_s=(3600.0,3600.0),
        inventory_hash="inventory-fixture",
        activation_artifact_hash="activation-fixture",
    )
    bundle = build_activation_delayed_source_bundle(
        scheme,activity;
        energy_grids=Dict(
            :photon => [0.0,5.0e5,1.5e6],
            :electron => [0.0,7.5e5,1.5e6],
        ),
    )
    photon = bundle.em_sources[:photon]
    electron = bundle.em_sources[:electron]
    @test photon.values[1,2,1] == 2.0 # 4 Bq / 2 cm3
    @test photon.values[1,1,1] == 0.0
    @test photon.variance[1,2,1] ≈ 0.04 # (0.4 Bq / 2 cm3)^2
    @test electron.values[1,1,1] == 2.0
    @test electron.values[1,2,1] == 0.0
    @test electron.parent_reaction == "X-100 beta-minus"
end

@testset "Decay energy cannot disappear into heat" begin
    incomplete = Delayed_Emission_Bin[
        Delayed_Emission_Bin(:photon,1.0e6,1.0),
        Delayed_Emission_Bin(:electron,5.0e5,1.0),
    ]
    @test_throws ErrorException Delayed_Decay_Scheme(
        "X-100","Y-100","beta-minus",3.0e6,incomplete;
        evaluation_id="incomplete",evaluation_hash="incomplete",
        status=:verification,
    )

    candidate = Delayed_Decay_Scheme(
        "X-100","Y-100","beta-minus",3.0e6,incomplete;
        unresolved_energy_eV_per_decay=1.5e6,
        evaluation_id="candidate-unresolved",evaluation_hash="candidate-unresolved",
        status=:candidate,
    )
    activity = Delayed_Activity_Field(
        "X-100",[1],[1.0],[1.0];
        cooling_interval_s=(0.0,1.0),
        inventory_hash="candidate-inventory",
        activation_artifact_hash="candidate-activation",
    )
    bundle = build_activation_delayed_source_bundle(
        candidate,activity;
        energy_grids=Dict(
            :photon => [0.0,2.0e6],
            :electron => [0.0,2.0e6],
        ),
    )
    @test bundle.unresolved_energy_rate_eV_per_s == 1.5e6
    @test !activation_delayed_bundle_is_production_ready(bundle)
    @test bundle.provenance["source_energy_is_not_deposited_heat"] == "true"
end

@testset "Delayed transported heat enters the ownership ledger only after scoring" begin
    ledger = Layer_Heating_Ledger()
    delayed_total = Layer_Heating_Contribution(
        contribution_id="delayed-photon-deposition",
        population_id="X-100-decay-photons",
        owner="Radiant",
        domain_id="tape-microdomain",
        layer_id="REBCO",
        particle_tag="photon-electron",
        source_class=:activation_delayed,
        time_class=:delayed,
        response_channel=:total_deposition,
        values=reshape([2.0],1,1,1),
        units="W/cm3",
        source_hash="delayed-source",
        geometry_hash="tape-geometry",
        normalization_hash="activity-normalization",
    )
    delayed_heat = Layer_Heating_Contribution(
        contribution_id="delayed-lattice-heat",
        population_id="X-100-decay-photons",
        owner="Radiant",
        domain_id="tape-microdomain",
        layer_id="REBCO",
        particle_tag="photon-electron",
        source_class=:activation_delayed,
        time_class=:delayed,
        response_channel=:delayed_lattice_heat,
        values=reshape([1.5],1,1,1),
        units="W/cm3",
        source_hash="delayed-source",
        geometry_hash="tape-geometry",
        normalization_hash="activity-normalization",
    )
    delayed_escape = Layer_Heating_Contribution(
        contribution_id="delayed-escape",
        population_id="X-100-decay-photons",
        owner="Radiant",
        domain_id="tape-microdomain",
        layer_id="REBCO",
        particle_tag="photon-electron",
        source_class=:activation_delayed,
        time_class=:delayed,
        response_channel=:escaped,
        values=reshape([0.5],1,1,1),
        units="W/cm3",
        source_hash="delayed-source",
        geometry_hash="tape-geometry",
        normalization_hash="activity-normalization",
    )
    add_heating_contribution!(ledger,delayed_total)
    add_heating_contribution!(ledger,delayed_heat)
    add_heating_contribution!(ledger,delayed_escape)
    @test validate_population_energy_closure(ledger,"X-100-decay-photons").closed
    result = production_heating_total(ledger;time_class=:delayed)
    @test only(result.values) == 1.5
    @test result.contribution_ids == ["delayed-lattice-heat"]
end
