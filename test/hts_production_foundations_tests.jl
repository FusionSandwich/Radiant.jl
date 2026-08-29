using Radiant
using Test
using LinearAlgebra
using HDF5

function _two_voxel_geometry()
    geometry = Geometry()
    geometry.type = "cartesian"
    geometry.dimension = 1
    geometry.axis = ["x"]
    geometry.number_of_voxels["x"] = 2
    geometry.voxels_width["x"] = [0.5,0.5]
    geometry.voxels_position["x"] = [0.25,0.75]
    geometry.voxels_boundaries["x"] = [0.0,0.5,1.0]
    geometry.volume_per_voxel = [0.5,0.5]
    geometry.material_per_voxel = ones(Int64,2,1,1)
    geometry.is_build = true
    return geometry
end

function _unit_boundary_source(;source_hash="unit-boundary",scale=1.0)
    normalization = Source_Normalization(
        basis=:per_history,source_rate_per_s=1.0,symmetry_factor=1.0,
        time_interval_s=(0.0,1.0),time_class=:prompt,source_hash=source_hash,
    )
    return Boundary_Angular_Current_Source(
        Photon(),[1],[0.0 0.0 0.0],[1.0],[-1.0 0.0 0.0],
        [0.0 1.0 0.0],[0.0 0.0 -1.0],[1.0e3,1.0e6],
        [1.0 0.0 0.0],[1.0],reshape([Float64(scale)],1,1,1),normalization,
    )
end

function _outward_cube_mesh()
    vertices = [
        0.0 0.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
        1.0 0.0 1.0
        1.0 1.0 1.0
        0.0 1.0 1.0
    ]
    triangles = Int64[
        1 3 2
        1 4 3
        5 6 7
        5 7 8
        1 2 6
        1 6 5
        4 8 7
        4 7 3
        1 5 8
        1 8 4
        2 3 7
        2 7 6
    ]
    return Faceted_Surface_Mesh(
        vertices,triangles;
        surface_ids=repeat(collect(1:6),inner=2),volume_ids=ones(Int64,12),
        material_tags=fill("cube-material",12),geometry_hash="cube-fixture",
        provenance=Dict("classification" => "verification-only"),
    )
end

@testset "Cell-varying magnetic field transport hook" begin
    field_values = zeros(Float64,3,2,1,1)
    field_values[3,1,1,1] = 10.0
    field_values[3,2,1,1] = 14.0
    field = Electromagnetic_Field()
    set_spatial_magnetic_field(
        field,field_values;field_hash="two-cell-field",
        provenance=Dict("producer" => "analytic-test"),
    )
    @test has_spatial_magnetic_field(field)
    @test get_field_mode(field) == :spatial_cell_centered
    @test get_spatial_field_hash(field) == "two-cell-field"
    @test get_spatial_magnetic_field(field) == field_values
    receipt = spatial_field_transport_receipt(field,_two_voxel_geometry())
    @test receipt["operator_application"] == "cell-local-within-one-sn-sweep"
    @test receipt["maximum_magnitude_T"] == 14.0
    @test !receipt["physical_qualification"]

    flux = reshape([2.0,3.0],1,1,2,1,1)
    source = zeros(Float64,1,1,2,1,1)
    operator = zeros(Float64,1,1,2,1,1)
    operator[1,1,1,1,1] = 4.0
    operator[1,1,2,1,1] = 5.0
    Radiant._add_electromagnetic_source!(source,flux,operator,[2,1,1],1,1)
    @test source[1,1,1,1,1] == 8.0
    @test source[1,1,2,1,1] == 15.0

    map_values = zeros(Float64,3,2,1,1)
    map_values[3,1,1,1] = 10.0
    map_values[3,2,1,1] = 14.0
    mapped = Cartesian_Magnetic_Field_Map(
        [0.25,0.75],[0.0],[0.0],map_values;field_hash="mapped-two-cell",
    )
    sampled = electromagnetic_field_on_geometry(mapped,_two_voxel_geometry())
    @test get_spatial_magnetic_field(sampled) == map_values
end

@testset "Faceted CAD and DAGMC interchange foundation" begin
    mesh = _outward_cube_mesh()
    report = faceted_topology_report(mesh)
    @test report["watertight"]
    @test report["manifold"]
    @test report["orientation_is_positive"]
    @test report["signed_enclosed_volume_cm3"] ≈ 1.0 atol=1.0e-12
    @test assert_transport_ready_facets(mesh)["watertight"]
    @test point_in_faceted_volume(mesh,[0.5,0.5,0.5])
    @test !point_in_faceted_volume(mesh,[1.5,0.5,0.5])

    temporary = mktempdir()
    path = joinpath(temporary,"cube_facets.h5")
    digest = write_faceted_geometry_hdf5(path,mesh)
    replay = read_faceted_geometry_hdf5(path;expected_file_sha256=digest)
    @test replay.vertices_cm == mesh.vertices_cm
    @test replay.triangles == mesh.triangles
    @test faceted_topology_report(replay)["watertight"]

    voxelization = voxelize_faceted_mesh(
        replay,[0.0,1.0],[0.0,1.0],[0.0,1.0];samples_per_axis=2,
    )
    @test voxelization.volume_fractions[1,1,1,1] ≈ 1.0
    @test voxelization.unresolved_fraction[1,1,1] ≈ 0.0
    @test voxelization.dominant_volume_id[1,1,1] == 1

    patches = build_faceted_tape_patches(replay;selected_surface_ids=[5])
    @test length(patches) == 2
    @test all(patch.area_cm2 ≈ 0.5 for patch in patches)
    @test all(
        patch.local_to_global'*patch.local_to_global ≈ Matrix{Float64}(I,3,3)
        for patch in patches
    )

    centroids,normals,areas = facet_geometry(replay)
    facet_index = findfirst(==(5),replay.surface_ids)
    source = Boundary_Angular_Current_Source(
        Photon(),[101],reshape(centroids[facet_index,:],1,3),[areas[facet_index]],
        reshape(normals[facet_index,:],1,3),
        reshape(patches[1].local_to_global[:,1],1,3),
        reshape(patches[1].local_to_global[:,2],1,3),
        [1.0e3,1.0e6],reshape(-normals[facet_index,:],1,3),[1.0],
        reshape([1.0],1,1,1),
        Source_Normalization(basis=:per_history,source_hash="facet-source"),
    )
    mapping = map_boundary_patches_to_facets(
        source,replay;maximum_centroid_distance_cm=1.0e-10,
    )
    @test mapping["mapping_count"] == 1
    @test mapping["mappings"][1]["facet_index"] == facet_index
end

@testset "Source-backed cryogenic material registry" begin
    registry = default_material_response_registry()
    copper_k = get_material_property(
        registry;material_tag="Cu-OFHC",property=:thermal_conductivity,
        state_contains="RRR=100",
    )
    copper_cp = get_material_property(
        registry;material_tag="Cu-OFHC",property=:specific_heat,
        state_contains="OFHC copper",
    )
    kapton_k = get_material_property(
        registry;material_tag="Kapton",property=:thermal_conductivity,
    )
    g10_normal = get_material_property(
        registry;material_tag="G10-CR",property=:thermal_conductivity,
        direction=:normal,
    )
    @test material_property_value(copper_k,20.0) > 0.0
    @test material_property_value(copper_cp,20.0) > 0.0
    @test material_property_value(kapton_k,20.0) > 0.0
    @test material_property_value(g10_normal,20.0) > 0.0
    table = tabulate_material_property(copper_k,[4.0,20.0,77.0,300.0])
    @test length(table["value"]) == 4
    @test all(table["standard_uncertainty"] .>= 0.0)
    cryogenic = cryogenic_property_from_record(copper_k,[4.0,20.0,77.0,300.0])
    @test property_value(cryogenic,20.0) ≈ material_property_value(copper_k,20.0)

    gdbco = get_material_property(
        registry;material_tag="GdBCO",property=:thermal_conductivity,
        direction=:crystal_tensor,
    )
    @test gdbco.status == :literature_table_pending
    @test_throws ErrorException material_property_value(gdbco,20.0)
end

@testset "Atomistic response-table pipeline" begin
    calculation = Atomistic_Calculation_Manifest(
        calculation_id="synthetic-phono3py",material_tag="YBCO",
        material_state_hash="ybco-state",method=:phono3py,code="phono3py",
        code_version="fixture",input_hash="input-hash",structure_hash="structure-hash",
        potential_or_functional="fixture-functional",
        convergence_parameters=Dict(
            "supercell" => "2x2x1","q_mesh" => "8x8x4","cutoff" => "fixture",
        ),status=:candidate,
    )
    table = Atomistic_Response_Table(
        table_id="fixture-kappa",material_tag="YBCO",material_state_hash="ybco-state",
        quantity=:thermal_conductivity_tensor,
        axis_order=["temperature_K","component_index"],
        axes=Dict("temperature_K" => [20.0,77.0],"component_index" => [1.0,2.0]),
        values=[10.0 2.0; 6.0 1.0],
        standard_uncertainty=[1.0 0.2; 0.6 0.1],units="W/(m*K)",
        calculation=calculation,table_hash="fixture-table",status=:candidate,
    )
    @test !atomistic_table_is_production_ready(table)
    temporary = mktempdir()
    path = joinpath(temporary,"atomistic_response.h5")
    digest = write_atomistic_response_hdf5(path,table)
    replay = read_atomistic_response_hdf5(path;expected_file_sha256=digest)
    @test replay.values == table.values
    @test replay.axes == table.axes
    @test replay.metadata["artifact_sha256"] == digest

    phono3py_path = joinpath(temporary,"kappa.hdf5")
    HDF5.h5open(phono3py_path,"w") do file
        file["temperature"] = [20.0,77.0]
        file["kappa"] = [10.0 8.0 6.0 0.0 0.0 0.0; 5.0 4.0 3.0 0.0 0.0 0.0]
    end
    phono_table = read_phono3py_kappa_hdf5(phono3py_path,calculation)
    @test size(phono_table.values) == (2,6)
    @test phono_table.metadata["physical_qualification"] == "false"

    epw_path = joinpath(temporary,"a2f.dat")
    write(epw_path,"# omega alpha2F\n1.0 0.1\n2.0 0.2\n3.0 0.15\n")
    epw_calculation = Atomistic_Calculation_Manifest(
        calculation_id="synthetic-epw",material_tag="YBCO",
        material_state_hash="ybco-state",method=:epw,code="EPW",code_version="fixture",
        input_hash="epw-input",structure_hash="structure-hash",
        potential_or_functional="fixture-functional",
        convergence_parameters=Dict("k_mesh" => "8x8x4","q_mesh" => "4x4x2"),
    )
    epw_table = read_epw_spectral_function(epw_path,epw_calculation)
    @test length(epw_table.values) == 3
    @test epw_table.metadata["interpretation"] == "spectral-input-not-direct-lifetime"
    @test length(default_ybco_gdbco_atomistic_plan()["calculations"]) == 4
end

@testset "Matched physical-reference qualification" begin
    channels = Dict(
        "Photoelectric/test" => reshape([10.0,20.0],2,1,1),
        "Compton/test" => reshape([5.0,8.0],2,1,1),
    )
    total = channels["Photoelectric/test"]+channels["Compton/test"]
    score = Process_Resolved_Score(
        "photon","energy-deposition",channels,total;
        units="MeV/g per transport source basis",
    )
    references = Dict{String,Physical_Reference_Response}()
    for key in keys(channels)
        references[key] = Physical_Reference_Response(
            response_id=key,process_key=key,particle_tag="photon",
            units=score.units,values=1.001.*channels[key],
            standard_uncertainty=0.01.*channels[key],
            classification=:continuous_energy_openmc,producer="OpenMC 0.16 fixture",
            source_artifact_hash="openmc-source",result_artifact_hash="openmc-result",
            geometry_hash="matched-geometry",material_state_hash="matched-material",
            normalization_hash="matched-normalization",
        )
    end
    qualification = qualify_process_resolved_score(
        score,references;relative_tolerance=0.01,absolute_tolerance=1.0e-12,
        physical_reference_only=true,
    )
    @test qualification.passed
    @test physical_qualification_receipt(qualification)["passed"]

    synthetic_reference = Physical_Reference_Response(
        response_id="synthetic",process_key="Photoelectric/test",particle_tag="photon",
        units=score.units,values=channels["Photoelectric/test"],
        standard_uncertainty=zeros(2,1,1),classification=:synthetic,
        producer="synthetic",source_artifact_hash="s",result_artifact_hash="r",
        geometry_hash="g",material_state_hash="m",normalization_hash="n",
    )
    @test !reference_is_physical(synthetic_reference)
end

@testset "Streaming Gd cascade data adapter" begin
    temporary = mktempdir()
    event_path = joinpath(temporary,"gd157_events.dat")
    write(event_path,join(fill("6.0 1.8 0.15",4),"\n")*"\n")
    digest = Radiant._source_file_sha256(event_path)
    q_value = 7.95e6+100.0
    manifest = Gd_Cascade_File_Manifest(
        nuclide="Gd-157",residual_nuclide="Gd-158",q_value_eV=q_value,
        mean_recoil_energy_eV=100.0,
        source_title="synthetic cascade fixture",source_identifier="fixture",
        source_url="https://example.invalid/fixture",expected_file_sha256=digest,
        energy_multiplier_to_eV=1.0e6,status=:qualified,
    )
    summary = summarize_gd_gamma_cascade_file(
        event_path,manifest,[0.0,2.0e5,2.0e6,7.0e6,8.0e6],
    )
    @test summary.event_count == 4
    @test summary.mean_multiplicity == 3.0
    @test abs(summary.energy_residual_eV) <= 1.0e-8
    @test gd_cascade_summary_is_production_ready(summary)
    rates = Capture_Rate_Field(
        "Gd-157",[1],[1.0],[2.0];transport_artifact_hash="self-shielded-rates",
        provenance=Dict("capture_rates_are_self_shielded" => "true"),
    )
    source = build_gd_histogram_volume_source(summary,rates)
    @test sum(get_volume_source_rate(source)) ≈ 6.0
    @test gd_cascade_data_sources()["gamma_event_repository"]["identifier"] ==
          "doi:10.5281/zenodo.7458654"
end

@testset "Two-way OpenSn-Radiant closed-coupling algorithm" begin
    initial = _unit_boundary_source(source_hash="forward-initial",scale=1.0)
    return_source = _unit_boundary_source(source_hash="return-fixed",scale=0.25)
    target_forward = _unit_boundary_source(source_hash="forward-fixed",scale=1.0)
    radiant_callback = forward -> Coupling_Solver_Output(
        return_source;protected_responses=[get_total_incoming_current(forward)],
        particle_balance_residual=0.0,energy_balance_residual=0.0,
        solver_id="Radiant-fixture",result_artifact_hash="radiant-result",
    )
    opensn_callback = returned -> Coupling_Solver_Output(
        target_forward;protected_responses=[get_total_incoming_current(returned)],
        particle_balance_residual=0.0,energy_balance_residual=0.0,
        solver_id="OpenSn-fixture",result_artifact_hash="opensn-result",
    )
    result = solve_closed_opensn_radiant_coupling(
        initial,radiant_callback,opensn_callback;
        settings=Closed_Coupling_Settings(
            maximum_iterations=5,minimum_iterations=2,relaxation=0.5,
            current_relative_tolerance=1.0e-12,
            energy_current_relative_tolerance=1.0e-12,
            response_relative_tolerance=1.0e-12,
        ),ownership_map_hash="ownership-map-fixture",
    )
    @test result.converged
    @test length(result.iterations) == 2
    receipt = closed_coupling_receipt(result)
    @test receipt["forward_current_closure_pass"]
    @test receipt["return_current_closure_pass"]
    @test receipt["energy_current_closure_pass"]
    @test !receipt["clipping_used"]
end
