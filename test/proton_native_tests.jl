using SHA
using TOML

function _manufactured_proton_case(;nuclear_stopping=0.0,angular_variance=0.0,
    nonelastic_removal=0.0,field_T=0.0,solver_type="CSD",
    qualification_mode=:software,electronic_stopping_by_layer=(2.0,2.0),
    group_boundaries_MeV=[500.0,1.0],source_values=[1.0],
    source_rate_per_s=1.0)
    proton = Proton()
    materials = Material[]
    tables = Proton_Material_Data[]
    nonelastic = Proton_Nonelastic_Data[]
    for (layer,id) in enumerate(("manufactured-front","manufactured-back"))
        material = Material(id)
        Radiant.set_density(material,1.0)
        push!(materials,material)
        hash = bytes2hex(sha256("$id/stopping/linear/1-500"))
        model = Tabulated_Ion_Transport_Model(
            species_id="proton",material_id=id,energy_MeV=[1.0,500.0],
            electronic_stopping_MeV_cm=fill(electronic_stopping_by_layer[layer],2),
            nuclear_stopping_MeV_cm=[nuclear_stopping,nuclear_stopping],
            energy_straggling_variance_MeV2_cm=[0.0,0.0],
            angular_variance_rad2_cm=[angular_variance,angular_variance],
            data_hash=hash,qualification_status=:synthetic,
        )
        push!(tables,Proton_Material_Data(material,model;
            material_state="solid",density_g_cm3=1.0,source_sha256=hash,
            uncertainty_status=:unquantified))
        event_hash = bytes2hex(sha256("$id/nonelastic/manufactured"))
        push!(nonelastic,Proton_Nonelastic_Data(
            material_id=id,energy_MeV=[1.0,500.0],
            removal_cm_inv=[nonelastic_removal,nonelastic_removal],
            event_family_ids=nonelastic_removal == 0.0 ? String[] : ["manufactured-event"],
            source_sha256=event_hash,qualification_status=:synthetic,
        ))
    end
    binding = bind_proton_native(proton,tables,nonelastic,group_boundaries_MeV;
        qualification_mode=qualification_mode)
    cs = binding.cross_sections
    Radiant.build(cs)

    geometry = Geometry()
    Radiant.set_dimension(geometry,1)
    Radiant.set_number_of_regions(geometry,"x",2)
    Radiant.set_region_boundaries(geometry,"x",[0.0,1.0,2.0])
    Radiant.set_voxels_per_region(geometry,"x",[1,1])
    Radiant.set_boundary_conditions(geometry,"x-","void")
    Radiant.set_boundary_conditions(geometry,"x+","void")
    Radiant.set_material_per_region(geometry,materials)
    Radiant.build(geometry,cs)

    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,solver_type)
    if field_T == 0.0
        Radiant.set_quadrature(solver,"gauss-legendre",2)
        Radiant.set_legendre_order(solver,0)
    else
        Radiant.set_quadrature(solver,"gauss-legendre-chebychev",2,3)
        Radiant.set_legendre_order(solver,1)
    end
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DD",1)
    Radiant.set_scheme(solver,"E","DG",1)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)

    normalization = Source_Normalization(
        basis=:per_source_particle,source_rate_per_s=source_rate_per_s,
        source_hash=bytes2hex(sha256("manufactured-proton-source")),
        provenance=Dict("classification" => "SOFTWARE_VERIFIED_SYNTHETIC"),
    )
    source = Anisotropic_Volume_Source(
        proton,[1],[1.0],reverse(group_boundaries_MeV).*1.0e6,:isotropic,
        reshape(Float64.(source_values),1,length(source_values),1),normalization,
    )
    sources = Fixed_Sources(cs,geometry,solvers)
    Radiant.add_source(sources,source)

    cu = Computation_Unit()
    Radiant.set_cross_sections(cu,cs)
    Radiant.set_geometry(cu,geometry)
    Radiant.set_solvers(cu,solvers)
    Radiant.set_sources(cu,sources)
    field = Electromagnetic_Field()
    Radiant.set_magnetic_field(field,[0.0,0.0,field_T])
    Radiant.set_electromagnetic_field(cu,field)
    return cu,binding,proton
end

@testset "Native manufactured proton binding" begin
    plus = proton_species()
    states = [relativistic_ion_kinematics(plus,E) for E in (1.0,100.0,500.0)]
    for (E,state) in zip((1.0,100.0,500.0),states)
        @test state.gamma ≈ 1.0+E/plus.rest_mass_MeV
        @test state.momentum_MeV_c^2 ≈ E*(E+2.0*plus.rest_mass_MeV)
    end
    @test states[1].beta < states[2].beta < states[3].beta < 1.0
    minus = Ion_Species(species_id="manufactured-negative-proton-charge",
        atomic_number=1,mass_number=1,rest_mass_MeV=plus.rest_mass_MeV,
        charge_e=-1.0,provenance_hash="synthetic-charge-sign")
    straight,angle_zero = magnetic_direction_step(plus,100.0,[1.0,0.0,0.0],
        [0.0,0.0,0.0],1.0)
    @test straight == (1.0,0.0,0.0) && angle_zero == 0.0
    _,angle_zero_path = magnetic_direction_step(plus,100.0,[1.0,0.0,0.0],
        [0.0,0.0,1.0],0.0)
    @test angle_zero_path == 0.0
    direction_plus,angle_plus = magnetic_direction_step(plus,100.0,
        [1.0,0.0,0.0],[0.0,0.0,1.0],1.0)
    direction_minus,angle_minus = magnetic_direction_step(minus,100.0,
        [1.0,0.0,0.0],[0.0,0.0,1.0],1.0)
    @test angle_plus ≈ -angle_minus
    @test direction_plus[2] ≈ -direction_minus[2]
    @test norm(collect(direction_plus)) ≈ norm(collect(direction_minus)) ≈ 1.0

    cu,binding,proton = _manufactured_proton_case()
    @test Radiant.get_number_of_groups(binding.cross_sections,proton) == 1
    @test Radiant.get_energy_boundaries(binding.cross_sections,proton) == [500.0,1.0]
    @test proton_binding_receipt(binding)["classification"] == "SOFTWARE_VERIFIED_SYNTHETIC"
    receipt_buffer = IOBuffer()
    TOML.print(receipt_buffer,proton_binding_receipt(binding))
    @test TOML.parse(String(take!(receipt_buffer)))["materials"][1]["id"] ==
        "manufactured-front"
    Radiant.run(cu)
    deposited = Radiant.get_energy_deposition(cu,proton)
    score = score_process_responses(binding.cross_sections,cu.geometry,cu.solvers,
        cu.sources,cu.flux,proton;quantity="energy-deposition")
    @test vec(score.total) ≈ deposited
    @test all(isfinite,score.total)
    @test all(score.total .>= 0.0)
    @test sum(score.total) > 0.0
    @test Radiant.get_normalization_factor(cu.sources) > 0.0
    @test all(iszero,Radiant.get_absorption(binding.cross_sections,proton))
    @test all(==(2.0),Radiant.get_stopping_powers(binding.cross_sections,proton))
    @test ion_csda_range_cm(binding.material_data[1].model,proton_species(),500.0,1.0) ≈ 249.5
    @test_throws ErrorException ion_csda_range_cm(binding.material_data[1].model,
        proton_species(),500.1,1.0)
    @test_throws ErrorException ion_csda_range_cm(binding.material_data[1].model,
        proton_species(),500.0,0.5)
    @test get_volume_source_rate(cu.sources.source_collection[1]) == [1.0]
    @test get_volume_source_rate(cu.sources.source_collection[1];physical=true) == [1.0]

    cu_n,bind_n,p_n = _manufactured_proton_case(nuclear_stopping=0.5)
    Radiant.run(cu_n)
    recoil = score_process_responses(bind_n.cross_sections,cu_n.geometry,cu_n.solvers,
        cu_n.sources,cu_n.flux,p_n;quantity="recoil-handoff")
    electronic = score_process_responses(bind_n.cross_sections,cu_n.geometry,cu_n.solvers,
        cu_n.sources,cu_n.flux,p_n;quantity="energy-deposition")
    @test vec(recoil.total) ≈ 0.25 .* vec(electronic.total)
    @test Radiant.get_stopping_powers(bind_n.cross_sections,p_n) == fill(2.5,1,2)

    cu_field,binding_field,proton_field = _manufactured_proton_case(field_T=0.01)
    Radiant.run(cu_field)
    @test all(isfinite,Radiant.get_flux(cu_field,proton_field))
    @test sum(Radiant.get_energy_deposition(cu_field,proton_field)) > 0.0

    cu_bfp,binding_bfp,proton_bfp = _manufactured_proton_case(
        angular_variance=0.001,solver_type="BFP")
    Radiant.run(cu_bfp)
    @test all(isfinite,Radiant.get_flux(cu_bfp,proton_bfp))
    @test Radiant.get_momentum_transfer(binding_bfp.cross_sections,proton_bfp) ==
        fill(0.0005,1,2)

    cu_bad,_,_ = _manufactured_proton_case(nonelastic_removal=0.01)
    @test_throws ErrorException Radiant.run(cu_bad)
    cu_ang,_,_ = _manufactured_proton_case(angular_variance=0.001)
    @test_throws ErrorException Radiant.run(cu_ang)
    cu_physical,binding_physical,_ = _manufactured_proton_case(
        qualification_mode=:physical)
    @test proton_binding_receipt(binding_physical)["classification"] == "BLOCKED_INPUT"
    @test_throws ErrorException Radiant.run(cu_physical)
    @test_throws ErrorException bind_proton_native(proton,binding.material_data,
        Proton_Nonelastic_Data[],[500.0,1.0])
    @test_throws ErrorException bind_proton_native(proton,binding.material_data,
        binding.nonelastic_data,[500.1,1.0])
    @test_throws ErrorException bind_proton_native(proton,binding.material_data,
        binding.nonelastic_data,[500.0,0.5])
    @test_throws ErrorException Proton_Material_Data(
        binding.material_data[1].material,binding.material_data[1].model;
        material_state="liquid",density_g_cm3=1.0,
        source_sha256=binding.material_data[1].source_sha256)
    @test_throws ErrorException Proton_Nonelastic_Data(
        material_id="manufactured-front",energy_MeV=[1.0,500.0],
        removal_cm_inv=[0.0,0.0],source_sha256=bytes2hex(sha256("zero")),
        qualification_status=:qualified)
    @test_throws ErrorException Proton_Nonelastic_Data(
        material_id="manufactured-front",energy_MeV=[1.0,500.0],
        removal_cm_inv=[0.1,0.1],source_sha256=bytes2hex(sha256("nonzero")))
    prior = binding.material_data[1]
    candidate = Tabulated_Ion_Transport_Model(
        species_id="proton",material_id="manufactured-front",energy_MeV=[1.0,500.0],
        electronic_stopping_MeV_cm=[2.0,2.0],nuclear_stopping_MeV_cm=[0.0,0.0],
        energy_straggling_variance_MeV2_cm=[0.0,0.0],
        angular_variance_rad2_cm=[0.0,0.0],data_hash=prior.source_sha256,
        qualification_status=:candidate)
    @test_throws ErrorException Proton_Material_Data(prior.material,candidate;
        material_state="solid",density_g_cm3=1.0,
        source_sha256=prior.source_sha256)
    cu_mutated,binding_mutated,_ = _manufactured_proton_case()
    binding_mutated.material_data[1].model.electronic_stopping_MeV_cm[1] = 3.0
    @test_throws ErrorException Radiant.run(cu_mutated)
end

@testset "Proton native repository integration" begin
    boundaries = [500.0,250.0,1.0]
    cu,binding,proton = _manufactured_proton_case(
        electronic_stopping_by_layer=(2.0,4.0),
        group_boundaries_MeV=boundaries,source_values=[0.0,1.0],
    )
    @test Radiant.get_number_of_groups(binding.cross_sections,proton) == 2
    @test Radiant.get_stopping_powers(binding.cross_sections,proton) ==
        [2.0 4.0; 2.0 4.0]
    Radiant.run(cu)
    @test get_projection_receipts(cu.sources)[1].energy_group_map == [2,1]
    flux = Radiant.get_flux(cu,proton)
    @test size(flux,1) == 2 && all(isfinite,flux)
    score = score_process_responses(binding.cross_sections,cu.geometry,cu.solvers,
        cu.sources,cu.flux,proton;quantity="energy-deposition")
    @test score.total[1,1,1] ≈ 2.0*sum(flux[:,1])
    @test score.total[2,1,1] ≈ 4.0*sum(flux[:,2])
    @test vec(score.total) ≈ Radiant.get_energy_deposition(cu,proton)
    @test all(score.total .> 0.0)

    cu_scaled,binding_scaled,proton_scaled = _manufactured_proton_case(
        electronic_stopping_by_layer=(2.0,4.0),
        group_boundaries_MeV=boundaries,source_values=[0.0,3.0],
    )
    Radiant.run(cu_scaled)
    @test Radiant.get_flux(cu_scaled,proton_scaled) ≈ 3.0 .* flux
    @test Radiant.get_energy_deposition(cu_scaled,proton_scaled) ≈
        3.0 .* Radiant.get_energy_deposition(cu,proton)

    cu_rate,binding_rate,proton_rate = _manufactured_proton_case(
        electronic_stopping_by_layer=(2.0,4.0),
        group_boundaries_MeV=boundaries,source_values=[0.0,1.0],
        source_rate_per_s=7.0,
    )
    Radiant.run(cu_rate)
    @test Radiant.get_flux(cu_rate,proton_rate) ≈ flux
    rate_score = score_process_responses(binding_rate.cross_sections,cu_rate.geometry,
        cu_rate.solvers,cu_rate.sources,cu_rate.flux,proton_rate;
        quantity="energy-deposition",
        physical_normalization=get_source_normalization(cu_rate.sources))
    @test rate_score.total ≈ 7.0 .* score.total
    @test rate_score.normalization_basis == :physical_rate

    cu_nonelastic,binding_nonelastic,_ = _manufactured_proton_case()
    binding_nonelastic.nonelastic_data[1].removal_cm_inv[1] = 0.1
    @test_throws ErrorException Radiant.run(cu_nonelastic)
    cu_density,binding_density,_ = _manufactured_proton_case()
    Radiant.set_density(binding_density.material_data[1].material,2.0)
    @test_throws ErrorException Radiant.run(cu_density)
    cu_boundaries,binding_boundaries,_ = _manufactured_proton_case()
    binding_boundaries.cross_sections.energy_boundaries[1][1] = 499.0
    @test_throws ErrorException Radiant.run(cu_boundaries)
    cu_dropped,binding_dropped,_ = _manufactured_proton_case()
    empty!(binding_dropped.nonelastic_data)
    @test_throws ErrorException Radiant.run(cu_dropped)
    cu_library,binding_library,_ = _manufactured_proton_case()
    binding_library.cross_sections.multigroup_cross_sections[1,1].total[1] = 0.1
    @test_throws ErrorException Radiant.run(cu_library)
end
