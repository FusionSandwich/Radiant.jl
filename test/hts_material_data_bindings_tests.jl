using Radiant
using Test

function _binding_calculation(
    method::Symbol,
    id::String;
    material::String="YBCO7",
    state_hash::String="ybco7-orthorhombic-verification",
)
    return Atomistic_Calculation_Manifest(
        calculation_id=id,
        material_tag=material,
        material_state_hash=state_hash,
        method=method,
        code="synthetic-response-producer",
        code_version="1",
        input_hash="$(id)-input-sha256",
        structure_hash="$(id)-structure-sha256",
        potential_or_functional="verification-only",
        convergence_parameters=Dict(
            "classification" => "software-verification",
            "replicas" => "2",
        ),
        references=["verification-only"],
        status=:candidate,
        metadata=Dict("physical_qualification" => "false"),
    )
end

@testset "Current atomistic heat-capacity binding" begin
    calculation = _binding_calculation(:dft_static,"heat-capacity")
    table = Atomistic_Response_Table(
        table_id="ybco-harmonic-cv",
        material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:heat_capacity,
        axis_order=["temperature_K"],
        axes=Dict("temperature_K" => [10.0,20.0]),
        values=[0.02,0.04],
        standard_uncertainty=[0.001,0.002],
        units="J/(cm^3*K)",
        calculation=calculation,
        table_hash="heat-capacity-table-sha256",
        status=:candidate,
        metadata=Dict("source_url" => "urn:test:heat-capacity"),
    )
    property = cryogenic_property_from_atomistic_table(table)
    @test property_value(property,10.0) == 0.02
    @test property_value(property,15.0) ≈ 0.03
    @test property.property_hash == table.table_hash
    @test property.qualification_status == :candidate

    record = material_record_from_atomistic_table(
        table;
        property=:volumetric_heat_capacity,
        direction=:isotropic,
        state_description="orthorhombic YBCO7 verification state",
    )
    @test record.record_id == table.table_id
    @test record.status == :atomistic_candidate
    @test record.material_tag == "YBCO7"
    @test record.nominal_relative_uncertainty ≈ 0.05
    @test material_property_value(record,15.0) ≈ 0.03
    @test material_property_standard_uncertainty(record.model,15.0) ≈ 0.0015
    @test record.metadata["calculation_input_hash"] == calculation.input_hash
end

@testset "Anisotropic conductivity binding and registry insertion" begin
    calculation = _binding_calculation(:green_kubo_md,"thermal-conductivity")
    table = Atomistic_Response_Table(
        table_id="ybco-k-tensor",
        material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:thermal_conductivity_tensor,
        axis_order=["temperature_K","component_index"],
        axes=Dict(
            "temperature_K" => [10.0,20.0],
            "component_index" => [1.0,2.0,3.0],
        ),
        values=[10.0 8.0 2.0; 12.0 9.0 3.0],
        standard_uncertainty=[1.0 0.8 0.2; 1.2 0.9 0.3],
        units="W/(m*K)",
        calculation=calculation,
        table_hash="conductivity-table-sha256",
        status=:candidate,
        metadata=Dict(
            "component_convention" => "xx,yy,zz",
            "source_url" => "urn:test:conductivity",
        ),
    )
    @test_throws ErrorException cryogenic_property_from_atomistic_table(table)
    normal_property = cryogenic_property_from_atomistic_table(table;component=3)
    @test normal_property.values == [2.0,3.0]
    @test property_value(normal_property,15.0) ≈ 2.5

    normal_record = material_record_from_atomistic_table(
        table;
        record_id="ybco-k-normal",
        property=:thermal_conductivity,
        direction=:normal,
        component=3,
    )
    in_plane_record = material_record_from_atomistic_table(
        table;
        record_id="ybco-k-in-plane",
        property=:thermal_conductivity,
        direction=:warp,
        component=1,
    )
    registry = Material_Response_Registry([normal_record])
    register_material_property!(registry,in_plane_record)
    @test length(registry.records) == 2
    @test get_material_property(
        registry;material_tag="YBCO7",property=:thermal_conductivity,
        direction=:normal,
    ).record_id == "ybco-k-normal"
    @test_throws ErrorException register_material_property!(registry,in_plane_record)
    replacement = material_record_from_atomistic_table(
        table;
        record_id="ybco-k-in-plane",
        property=:thermal_conductivity,
        direction=:warp,
        component=2,
        notes="Replacement used only to test vector-registry update semantics.",
    )
    register_material_property!(registry,replacement;replace=true)
    @test length(registry.records) == 2
    @test get_material_property(
        registry;material_tag="YBCO7",property=:thermal_conductivity,
        direction=:warp,
    ).model.values == [8.0,9.0]
end

@testset "Sub-keV response-table compatibility binding" begin
    calculation = _binding_calculation(:rt_tddft,"subkev")
    channels = String[
        "ionization","electronic_excitation","prompt_phonon","athermal_phonon",
        "quasiparticle","defect_storage","optical_emission","electron_escape",
    ]
    values = zeros(Float64,2,2,8)
    first_partition = [0.10,0.10,0.20,0.20,0.20,0.05,0.05,0.10]
    second_partition = [0.15,0.10,0.20,0.15,0.15,0.05,0.10,0.10]
    for energy_index in 1:2
        values[energy_index,1,:] .= first_partition
        values[energy_index,2,:] .= second_partition
    end
    table = Atomistic_Response_Table(
        table_id="ybco-subkev-partition",
        material_tag=calculation.material_tag,
        material_state_hash=calculation.material_state_hash,
        quantity=:energy_partition,
        axis_order=["energy_eV","temperature_K","channel_index"],
        axes=Dict(
            "energy_eV" => [10.0,100.0],
            "temperature_K" => [10.0,20.0],
            "channel_index" => collect(1.0:8.0),
        ),
        values=values,
        standard_uncertainty=fill(0.01,size(values)),
        units="fraction",
        calculation=calculation,
        table_hash="subkev-table-sha256",
        status=:candidate,
        metadata=Dict("channel_names" => join(channels,',')),
    )
    current_binding = subkev_kernel_from_response_table(table,15.0)
    compatibility_kernel = subkev_kernel_from_atomistic_table(
        table,15.0;model_id="custom-subkev-binding",
    )
    @test compatibility_kernel.model_id == "custom-subkev-binding"
    @test compatibility_kernel.model_hash == table.table_hash
    @test compatibility_kernel.material_state_hash == table.material_state_hash
    @test compatibility_kernel.energy_grid_eV == current_binding.kernel.energy_grid_eV
    for channel in SUBKEV_PARTITION_CHANNELS
        @test compatibility_kernel.fractions[channel] ≈
              current_binding.kernel.fractions[channel]
    end
    @test is_subkev_partition_closed(
        thermalize_subkev_event(compatibility_kernel,100.0),
    )
end

@testset "Material response source map rejects implicit YBCO-to-GdBCO substitution" begin
    source_map = hts_material_response_source_map()
    @test haskey(source_map,"YBCO")
    @test haskey(source_map,"GdBCO")
    @test occursin("Do not reuse YBCO",source_map["GdBCO"]["prohibited_shortcut"])
    @test "electronic stopping from real-time TDDFT" in
          source_map["YBCO"]["generated_properties"]
end
