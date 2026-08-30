#!/usr/bin/env julia

using Radiant
using TOML
using Dates

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT = joinpath(
    ROOT,"qualification-results","julia-$(VERSION)",
    "gdbco-capture-source-chain-analytic.toml",
)

function main()
    neutron_edges_eV = [0.0,1.0e3,1.0e6]
    gd157 = Gd_Groupwise_Capture_Component(
        "Gd-157",2.0e22,neutron_edges_eV,[20.0,2.0];
        data_hash="synthetic-gd157-capture-xs-v1",
        status=:verification,
        provenance=Dict(
            "classification" => "analytic-verification",
            "evaluated_nuclear_data" => "false",
        ),
    )
    layer = Gd_Self_Shielding_Layer(
        "gdbco-active-layer","GdBCO-analytic",0.2,[gd157];
        subdivisions=4,
        background_absorption_cm_inv=[0.02,0.0],
        metadata=Dict(
            "cross_sectional_area_cm2" => "1.0",
            "classification" => "analytic-verification",
        ),
    )
    incident_current_per_s = [1000.0,2000.0]
    shielding = solve_gd_groupwise_self_shielding(
        [layer],incident_current_per_s;
        direction_cosine=0.6,
        geometry_hash="analytic-gdbco-four-cell-layer",
        transport_artifact_hash="analytic-gdbco-capture-field-v1",
        metadata=Dict("source" => "software-verification-neutron-current"),
    )
    shielding_receipt = gd_self_shielding_receipt(shielding)
    shielding_receipt["particle_balance_pass"] || error(
        "Analytic Gd self-shielding particle balance failed.",
    )

    cell_volume_cm3 = layer.thickness_cm/layer.subdivisions
    capture_field = capture_rate_field_from_gd_self_shielding(
        shielding,"Gd-157",collect(1:layer.subdivisions),
        fill(cell_volume_cm3,layer.subdivisions);
        material_tags=fill("GdBCO-analytic",layer.subdivisions),
        provenance=Dict(
            "classification" => "software-verification",
            "cross_sectional_area_cm2" => "1.0",
        ),
    )
    total_capture_rate_per_s = sum(capture_field.capture_rate_per_s)
    total_capture_rate_per_s > 0.0 || error("Analytic Gd capture rate is zero.")
    isapprox(
        total_capture_rate_per_s,
        sum(shielding.capture_rate_per_s["Gd-157"]);
        rtol=1.0e-13,atol=1.0e-13,
    ) || error("Capture-rate field does not preserve the self-shielding capture rate.")

    cascade = synthetic_gd157_capture_fixture()
    source_edges_eV = [0.0,1.0e4,1.0e5,1.0e6,1.0e7]
    bundle = build_gd_capture_source_bundle(
        cascade,capture_field;
        photon_energy_edges_eV=source_edges_eV,
        electron_energy_edges_eV=source_edges_eV,
    )
    expected_species_rates = Dict{String,Float64}()
    observed_species_rates = Dict{String,Float64}()
    for species in (:photon,:xray,:conversion_electron,:auger_electron)
        lines = [line for line in cascade.emissions if line.species == species]
        expected_rate = total_capture_rate_per_s*
                        sum(line.yield_per_capture for line in lines)
        haskey(bundle.sources,species) || error(
            "Gd capture-source bundle is missing species $(species).",
        )
        observed_rate = sum(get_volume_source_rate(bundle.sources[species]))
        isapprox(observed_rate,expected_rate;rtol=1.0e-12,atol=1.0e-12) || error(
            "Gd capture-source rate closure failed for $(species).",
        )
        bundle.sources[species].normalization.time_class == :prompt || error(
            "Gd prompt-capture source was not classified as prompt.",
        )
        expected_species_rates[string(species)] = expected_rate
        observed_species_rates[string(species)] = observed_rate
    end

    isnothing(bundle.recoil_source) && error("Gd capture recoil source is missing.")
    isapprox(
        sum(bundle.recoil_source.event_rate_per_s),total_capture_rate_per_s;
        rtol=1.0e-12,atol=1.0e-12,
    ) || error("Gd recoil event rate does not equal the capture rate.")
    q_rate_eV_per_s = total_capture_rate_per_s*cascade.q_value_eV
    emitted_rate_eV_per_s = total_capture_rate_per_s*
                            bundle.prompt_emitted_energy_eV_per_capture
    recoil_rate_eV_per_s = total_capture_rate_per_s*
                           bundle.recoil_energy_eV_per_capture
    isapprox(
        emitted_rate_eV_per_s+recoil_rate_eV_per_s,q_rate_eV_per_s;
        rtol=1.0e-12,atol=1.0e-8,
    ) || error("Gd capture-source energy rate does not close the reaction Q value.")

    # The source bundle represents particles created by capture, not energy deposited locally.
    # Heating ownership begins only after those particles are transported or handed to the recoil
    # model. Keep this explicit in the qualification receipt.
    receipt = Dict{String,Any}(
        "schema" => "radiant.hts.gdbco_capture_source_chain_analytic/v1",
        "classification" => "software-verification",
        "status" => "GDBCO_CAPTURE_SOURCE_CHAIN_ANALYTIC_PASS",
        "completed_at" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "geometry_hash" => shielding.geometry_hash,
        "transport_artifact_hash" => shielding.transport_artifact_hash,
        "capture_data_hashes" => shielding.data_hashes,
        "capture_rate_per_s" => total_capture_rate_per_s,
        "expected_species_rates_per_s" => expected_species_rates,
        "observed_species_rates_per_s" => observed_species_rates,
        "q_value_eV_per_capture" => cascade.q_value_eV,
        "emitted_energy_rate_eV_per_s" => emitted_rate_eV_per_s,
        "recoil_energy_rate_eV_per_s" => recoil_rate_eV_per_s,
        "q_value_energy_rate_eV_per_s" => q_rate_eV_per_s,
        "self_shielding_particle_balance_pass" => true,
        "capture_to_secondary_rate_closure_pass" => true,
        "capture_energy_closure_pass" => true,
        "capture_depth_resolved" => true,
        "prompt_and_recoil_sources_separate" => true,
        "source_energy_is_not_deposited_heat" => true,
        "physical_neutron_transport" => false,
        "evaluated_gd_data" => false,
        "physical_gate_promoted" => false,
        "shielding_receipt" => shielding_receipt,
    )
    mkpath(dirname(OUTPUT))
    open(OUTPUT,"w") do io
        TOML.print(io,receipt)
    end
    println("GDBCO_CAPTURE_SOURCE_CHAIN_ANALYTIC_PASS")
    println("GdBCO capture-source receipt: $(OUTPUT)")
end

try
    main()
catch exception
    println(stderr,"GDBCO_CAPTURE_SOURCE_CHAIN_ANALYTIC_PASS not established.")
    showerror(stderr,exception,catch_backtrace())
    println(stderr)
    exit(1)
end
