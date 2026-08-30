using Radiant
using Test
using SHA

function _test_file_sha256(path::AbstractString)
    return open(path,"r") do io
        bytes2hex(SHA.sha256(io))
    end
end

@testset "Complete evaluated Gd cascade adapter" begin
    directory = mktempdir()
    path = joinpath(directory,"gd-complete-verification.csv")
    open(path,"w") do io
        println(io,
            "target_nuclide,residual_nuclide,branch_id,branch_probability,q_value_eV," *
            "species,energy_eV,yield_per_branch,delay_s,correlation_id,transition_id," *
            "recoil_energy_eV",
        )
        for row in (
            ("Gd-155","Gd-156","b1",1.0,8.0e6,"photon",7.0e6,1.0,0.0,"c1","g1",100.0),
            ("Gd-155","Gd-156","b1",1.0,8.0e6,"xray",5.0e4,1.0,0.0,"c1","x1",100.0),
            ("Gd-155","Gd-156","b1",1.0,8.0e6,"conversion_electron",9.0e5,1.0,0.0,"c1","ce1",100.0),
            ("Gd-155","Gd-156","b1",1.0,8.0e6,"auger_electron",4.99e4,1.0,0.0,"c1","a1",100.0),
            ("Gd-157","Gd-158","b1",1.0,7.9e6,"photon",7.0e6,1.0,0.0,"c2","g1",100.0),
            ("Gd-157","Gd-158","b1",1.0,7.9e6,"xray",4.0e4,1.0,0.0,"c2","x1",100.0),
            ("Gd-157","Gd-158","b1",1.0,7.9e6,"conversion_electron",8.0e5,1.0,0.0,"c2","ce1",100.0),
            ("Gd-157","Gd-158","b1",1.0,7.9e6,"auger_electron",5.99e4,1.0,0.0,"c2","a1",100.0),
        )
            println(io,join(row,','))
        end
    end
    digest = _test_file_sha256(path)
    dataset = read_evaluated_gd_cascade_csv(
        path;
        evaluation_id="synthetic-complete-cascade",
        evaluation_version="verification-v1",
        expected_file_sha256=digest,
        qualification_status=:verification,
        energy_rtol=1.0e-12,
        energy_atol_eV=1.0e-8,
        provenance=Dict(
            "source_classification" => "synthetic-verification",
            "physical_qualification" => "false",
        ),
    )

    @test Set(keys(dataset.cascades)) == Set(["Gd-155","Gd-157"])
    @test dataset.input_file_sha256 == digest
    @test !evaluated_gd_dataset_is_production_ready(dataset)
    receipt = gd_dataset_receipt(dataset)
    @test receipt["production_ready"] == false
    @test receipt["isotopes"]["Gd-155"]["branch_probabilities"]["b1"] == 1.0

    expected_species = Set([
        :photon,:xray,:conversion_electron,:auger_electron,
    ])
    energy_edges_eV = [0.0,1.0e4,1.0e5,1.0e6,1.0e7]
    for target in ("Gd-155","Gd-157")
        cascade = dataset.cascades[target]
        @test Set(getfield.(cascade.emissions,:species)) == expected_species
        emitted = sum(
            line.energy_eV*line.yield_per_capture for line in cascade.emissions
        )
        @test emitted+cascade.mean_recoil_energy_eV ≈ cascade.q_value_eV atol=1.0e-8
        rates = Capture_Rate_Field(
            target,[1],[2.0],[4.0];
            variance_per_s2=[0.04],
            material_tags=["GdBCO"],
            transport_artifact_hash="self-shielded-$(target)-capture-rate",
            provenance=Dict("capture_rates_are_self_shielded" => "true"),
        )
        bundle = build_gd_capture_source_bundle(
            cascade,rates;
            photon_energy_edges_eV=energy_edges_eV,
            electron_energy_edges_eV=energy_edges_eV,
            energy_closure_rtol=1.0e-12,
            energy_closure_atol_eV=1.0e-8,
        )
        @test Set(keys(bundle.sources)) == expected_species
        @test bundle.recoil_source !== nothing
        @test abs(bundle.energy_residual_eV_per_capture) ≤ 1.0e-8
        for source in values(bundle.sources)
            @test source.normalization.time_class == :prompt
            @test source.provenance["capture_rates_are_self_shielded"] == "true"
            @test sum(get_volume_source_rate(source)) ≈ 4.0
        end
    end

    incomplete_path = joinpath(directory,"gd-incomplete.csv")
    open(incomplete_path,"w") do io
        println(io,
            "target_nuclide,residual_nuclide,branch_id,branch_probability,q_value_eV," *
            "species,energy_eV,yield_per_branch,delay_s,correlation_id,transition_id," *
            "recoil_energy_eV",
        )
        println(io,"Gd-155,Gd-156,b1,1.0,8000000,photon,7000000,1.0,0,c1,g1,100")
    end
    @test_throws ErrorException read_evaluated_gd_cascade_csv(
        incomplete_path;
        evaluation_id="incomplete",
        evaluation_version="verification-v1",
        qualification_status=:verification,
        energy_rtol=1.0e-12,
        energy_atol_eV=1.0e-8,
    )

    wrong_hash = repeat("0",64)
    @test_throws ErrorException read_evaluated_gd_cascade_csv(
        path;
        evaluation_id="wrong-hash",
        evaluation_version="verification-v1",
        expected_file_sha256=wrong_hash,
        qualification_status=:verification,
    )
end
