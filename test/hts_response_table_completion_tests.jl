using Radiant
using Test
using SHA

function _response_calculation(method::Symbol; id::String, material::String="YBCO")
    return Atomistic_Calculation_Manifest(
        calculation_id=id,
        material_tag=material,
        material_state_hash="$(material)-synthetic-state",
        method=method,
        code="synthetic-test-producer",
        code_version="1",
        input_hash="$(id)-input-hash",
        structure_hash="$(id)-structure-hash",
        potential_or_functional="verification-only",
        convergence_parameters=Dict(
            "replicas" => "2",
            "classification" => "software-verification",
        ),
        references=["verification-only-no-physical-claim"],
        status=:candidate,
        metadata=Dict("physical_qualification" => "false"),
    )
end

@testset "Current-format HTS material response completion" begin
    directory = mktempdir()

    green_kubo_path = joinpath(directory,"green-kubo.csv")
    open(green_kubo_path,"w") do io
        println(io,"temperature_K,replica_id,kxx_W_mK,kyy_W_mK,kzz_W_mK")
        println(io,"10,r1,10,8,2")
        println(io,"10,r2,14,10,4")
        println(io,"20,r1,8,6,1")
        println(io,"20,r2,12,8,3")
    end
    green_kubo = read_green_kubo_replicas_csv(
        green_kubo_path,
        _response_calculation(:green_kubo_md;id="gk");
        status=:candidate,
    )
    @test green_kubo.quantity == :thermal_conductivity_tensor
    @test green_kubo.axis_order == ["temperature_K","component_index"]
    @test green_kubo.values ≈ [12.0 9.0 3.0; 10.0 7.0 2.0]
    @test green_kubo.standard_uncertainty ≈ [2.0 1.0 1.0; 2.0 1.0 1.0]
    @test !atomistic_table_is_production_ready(green_kubo)
    conductance = cryogenic_conductance_from_conductivity_table(
        green_kubo;
        component_index=1,
        area_cm2=0.4,
        length_cm=1.0,
    )
    @test conductance.units == "W/K"
    @test conductance.values ≈ [0.048,0.040]

    subkev_path = joinpath(directory,"subkev.csv")
    channels = SUBKEV_RESPONSE_CHANNEL_NAMES
    sigma_columns = [string(channel,"_sigma") for channel in channels]
    open(subkev_path,"w") do io
        println(io,join(vcat(["energy_eV","temperature_K"],channels,sigma_columns),','))
        for energy in (10.0,100.0), temperature in (10.0,20.0)
            fractions = temperature == 10.0 ?
                [0.10,0.10,0.20,0.20,0.20,0.05,0.05,0.10] :
                [0.15,0.10,0.20,0.15,0.15,0.05,0.10,0.10]
            sigmas = fill(0.01,length(channels))
            println(io,join(vcat([energy,temperature],fractions,sigmas),','))
        end
    end
    subkev = read_subkev_partition_csv(
        subkev_path,
        _response_calculation(:rt_tddft;id="subkev");
        status=:candidate,
    )
    @test subkev.quantity == :energy_partition
    @test size(subkev.values) == (2,2,8)
    @test subkev.standard_uncertainty !== nothing
    binding = subkev_kernel_from_response_table(subkev,15.0)
    @test binding.interpolation == "linear"
    @test binding.kernel.temperature_K == 15.0
    @test binding.channel_standard_uncertainty !== nothing
    for energy_index in eachindex(binding.kernel.energy_grid_eV)
        @test sum(
            binding.kernel.fractions[channel][energy_index]
            for channel in SUBKEV_PARTITION_CHANNELS
        ) ≈ 1.0 atol=1.0e-12
    end
    thermalized = thermalize_subkev_event(binding.kernel,100.0)
    @test is_subkev_partition_closed(thermalized)
    @test binding.kernel.qualification_status == :candidate

    phonopy_path = joinpath(directory,"thermal_properties.yaml")
    open(phonopy_path,"w") do io
        println(io,"thermal_properties:")
        println(io,"- temperature: 10.0")
        println(io,"  free_energy: -1.0")
        println(io,"  entropy: 0.1")
        println(io,"  heat_capacity: 2.0")
        println(io,"- temperature: 20.0")
        println(io,"  free_energy: -0.9")
        println(io,"  entropy: 0.2")
        println(io,"  heat_capacity: 4.0")
    end
    heat_capacity = read_phonopy_heat_capacity_yaml(
        phonopy_path,
        _response_calculation(:dft_static;id="phonopy");
        density_g_cm3=6.0,
        molar_mass_g_mol=600.0,
        standard_uncertainty=[0.1,0.2],
        status=:candidate,
    )
    @test heat_capacity.units == "J/(cm^3*K)"
    @test heat_capacity.values ≈ [0.02,0.04]
    @test heat_capacity.standard_uncertainty ≈ [0.001,0.002]
    node_capacity = cryogenic_heat_capacity_from_response_table(heat_capacity,0.5)
    @test node_capacity.units == "J/K"
    @test node_capacity.values ≈ [0.01,0.02]

    bad_partition_path = joinpath(directory,"bad-subkev.csv")
    open(bad_partition_path,"w") do io
        println(io,join(vcat(["energy_eV","temperature_K"],channels),','))
        println(io,join(vcat([10.0,10.0],fill(0.2,8)),','))
    end
    @test_throws ErrorException read_subkev_partition_csv(
        bad_partition_path,
        _response_calculation(:rt_tddft;id="bad-subkev");
        status=:candidate,
    )

    duplicate_path = joinpath(directory,"duplicate-green-kubo.csv")
    open(duplicate_path,"w") do io
        println(io,"temperature_K,replica_id,kxx_W_mK,kyy_W_mK,kzz_W_mK")
        println(io,"10,r1,10,8,2")
        println(io,"10,r1,11,9,3")
    end
    @test_throws ErrorException read_green_kubo_replicas_csv(
        duplicate_path,
        _response_calculation(:green_kubo_md;id="duplicate");
        status=:candidate,
    )
end
