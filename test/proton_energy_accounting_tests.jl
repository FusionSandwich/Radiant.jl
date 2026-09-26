using Test
using SHA

function _cutoff_accounting_fixture(;kwargs...)
    defaults = (scalar_cutoff_flux=[8.0,4.0],last_group_width_MeV=4.0,
        cutoff_MeV=1.0,total_stopping_MeV_cm=[3.0,6.0],
        cell_measure_cm_d=[2.0,5.0],density_g_cm3=[2.0,4.0],
        geometry_dimension=1,source_normalization_divisor=2.0,
        source_basis=:per_source_particle)
    return proton_cutoff_handoff(;merge(defaults,(;kwargs...))...)
end

@testset "Explicit proton cutoff handoff algebra" begin
    # Independent pencil-and-paper result: 3*8/(4*2)=3, 6*4/(4*2)=3.
    score = _cutoff_accounting_fixture()
    @test score.cutoff_particle_flow_density == [3.0,3.0]
    @test score.cutoff_particles == 21.0
    @test score.cutoff_kinetic_handoff_MeV == 21.0
    @test score.cutoff_kinetic_handoff_per_mass == [1.5,0.75]
    @test score.energy_per_mass_units == "MeV*cm^2/g per per_source_particle"
    @test !score.physical_validation && !score.full_stopping_qualified
    @test _cutoff_accounting_fixture(scalar_cutoff_flux=zeros(2)).cutoff_particles == 0.0
    @test _cutoff_accounting_fixture(total_stopping_MeV_cm=zeros(2)).cutoff_particles == 0.0
    # Stored boundary coefficient doubles when last group width doubles;
    # physical φ(Ecut), boundary particle flow and energy must stay unchanged.
    wide = _cutoff_accounting_fixture(last_group_width_MeV=8.0,
        scalar_cutoff_flux=[16.0,8.0])
    @test wide.cutoff_particles == score.cutoff_particles
    @test wide.cutoff_kinetic_handoff_MeV == score.cutoff_kinetic_handoff_MeV
    @test _cutoff_accounting_fixture(cutoff_MeV=2.0).cutoff_kinetic_handoff_MeV == 42.0
    dense = _cutoff_accounting_fixture(density_g_cm3=[4.0,8.0])
    @test dense.cutoff_particles == score.cutoff_particles
    @test dense.cutoff_kinetic_handoff_per_mass == score.cutoff_kinetic_handoff_per_mass ./ 2
    @test _cutoff_accounting_fixture(cell_measure_cm_d=[4.0,10.0]).cutoff_particles == 42.0
    @test _cutoff_accounting_fixture(source_normalization_divisor=4.0).cutoff_particles == 10.5
    @test _cutoff_accounting_fixture(scalar_cutoff_flux=[16.0,8.0],
        source_normalization_divisor=4.0).cutoff_particles == 21.0
    @test _cutoff_accounting_fixture(geometry_dimension=3).energy_per_mass_units ==
        "MeV*cm^0/g per per_source_particle"
    normalization = Source_Normalization(basis=:per_source_particle,
        source_rate_per_s=2.0,symmetry_factor=3.0)
    @test _cutoff_accounting_fixture(physical_normalization=normalization).cutoff_particles == 126.0
    per_second = Source_Normalization(basis=:per_second,source_rate_per_s=100.0,symmetry_factor=3.0)
    @test _cutoff_accounting_fixture(source_basis=:per_second,
        physical_normalization=per_second).cutoff_particles == 63.0
    @test_throws ErrorException _cutoff_accounting_fixture(physical_normalization=per_second)
    for invalid in (-1.0,NaN,Inf)
        @test_throws ErrorException _cutoff_accounting_fixture(scalar_cutoff_flux=[invalid,0.0])
        @test_throws ErrorException _cutoff_accounting_fixture(total_stopping_MeV_cm=[invalid,0.0])
        @test_throws ErrorException _cutoff_accounting_fixture(last_group_width_MeV=invalid)
        @test_throws ErrorException _cutoff_accounting_fixture(source_normalization_divisor=invalid)
        @test_throws ErrorException _cutoff_accounting_fixture(density_g_cm3=[invalid,1.0])
        @test_throws ErrorException _cutoff_accounting_fixture(cell_measure_cm_d=[invalid,1.0])
    end
    @test_throws ErrorException _cutoff_accounting_fixture(last_group_width_MeV=0.0)
    @test_throws ErrorException _cutoff_accounting_fixture(source_normalization_divisor=0.0)
    @test_throws ErrorException _cutoff_accounting_fixture(density_g_cm3=[0.0,1.0])
    @test_throws ErrorException _cutoff_accounting_fixture(cell_measure_cm_d=[0.0,1.0])
    @test_throws ErrorException _cutoff_accounting_fixture(scalar_cutoff_flux=[1.0])
    @test_throws ErrorException _cutoff_accounting_fixture(geometry_dimension=0)
    @test_throws ErrorException _cutoff_accounting_fixture(source_basis=:unknown)
    @test_throws ErrorException _cutoff_accounting_fixture(cutoff_MeV=0.999)
    @test_throws ErrorException _cutoff_accounting_fixture(cutoff_MeV=500.001)
    @test_throws ErrorException _cutoff_accounting_fixture(cell_measure_cm_d=fill(floatmax(Float64),2))
    @test_throws ErrorException _cutoff_accounting_fixture(density_g_cm3=fill(nextfloat(0.0),2))
    # A finite denominator product can exceed representability. Never turn a
    # physically nonzero .25 flow into zero through an infinite denominator.
    @test_throws ErrorException _cutoff_accounting_fixture(
        total_stopping_MeV_cm=fill(floatmax(Float64),2),scalar_cutoff_flux=ones(2),
        source_normalization_divisor=floatmax(Float64))
end

@testset "Missing energy observations remain missing" begin
    score = merge(_cutoff_accounting_fixture(),
        (electronic_deposition_MeV=4.0,recoil_handoff_MeV=2.0))
    partial = proton_energy_accounting_report(score;injected_energy_MeV=27.0)
    @test partial["classification"] == "BLOCKED_INCOMPLETE_LEDGER"
    @test ismissing(partial["escaped_energy_MeV"])
    @test ismissing(partial["arithmetic_closure"])
    @test "verified_convergence" in partial["unobserved_or_unverified"]
    complete = (injected_energy_MeV=27.0,escaped_energy_MeV=0.0,
        secondary_transfer_MeV=0.0,other_transfer_MeV=0.0,convergence_verified=true)
    @test proton_energy_accounting_report(score;complete...)["classification"] ==
        "SOFTWARE_ARITHMETIC_CLOSED"
    # Double-counting the residual cutoff as electronic deposition must fail.
    false_heat = merge(score,(electronic_deposition_MeV=25.0,))
    @test proton_energy_accounting_report(false_heat;complete...)["classification"] ==
        "FAIL_ENERGY_CLOSURE"
    unverified = merge(complete,(convergence_verified=false,))
    @test proton_energy_accounting_report(score;unverified...)["classification"] ==
        "BLOCKED_INCOMPLETE_LEDGER"
    @test_throws ErrorException proton_energy_accounting_report(score;escaped_energy_MeV=-1.0)
    @test_throws ErrorException proton_energy_accounting_report(score;secondary_transfer_MeV=NaN)
    @test_throws ErrorException proton_energy_accounting_report(score;closure_atol_MeV=-1.0)
    @test_throws ErrorException proton_energy_accounting_report(score;convergence_verified="yes")
    huge = floatmax(Float64)
    huge_score = merge(score,(electronic_deposition_MeV=huge,recoil_handoff_MeV=huge,))
    @test_throws ErrorException proton_energy_accounting_report(huge_score;complete...)
    @test_throws ErrorException proton_energy_accounting_report(score;
        merge(complete,(injected_energy_MeV=huge,))...,closure_atol_MeV=huge,closure_rtol=1.0)
    @test_throws ErrorException proton_energy_accounting_report(score;
        merge(complete,(escaped_energy_MeV=huge,secondary_transfer_MeV=huge))...)
end

@testset "Native proton prescribed moments and ownership" begin
    proton = Proton()
    material = Material("accounting-manufactured")
    Radiant.set_density(material,2.0)
    hash = bytes2hex(sha256("accounting-manufactured-linear-table"))
    model = Tabulated_Ion_Transport_Model(species_id="proton",material_id="accounting-manufactured",
        energy_MeV=[1.0,10.0],electronic_stopping_MeV_cm=[2.0,2.0],
        nuclear_stopping_MeV_cm=[1.0,1.0],energy_straggling_variance_MeV2_cm=[0.0,0.0],
        angular_variance_rad2_cm=[0.0,0.0],data_hash=hash,qualification_status=:synthetic)
    data = Proton_Material_Data(material,model;material_state="solid",density_g_cm3=2.0,
        source_sha256=hash)
    nonelastic = Proton_Nonelastic_Data(material_id="accounting-manufactured",
        energy_MeV=[1.0,10.0],removal_cm_inv=[0.0,0.0],source_sha256=hash)
    binding = bind_proton_native(proton,[data],[nonelastic],[10.0,5.0,1.0])
    cs = binding.cross_sections
    Radiant.build(cs)
    geometry = Geometry()
    Radiant.set_dimension(geometry,1)
    Radiant.set_number_of_regions(geometry,"x",1)
    Radiant.set_region_boundaries(geometry,"x",[0.0,2.0])
    Radiant.set_voxels_per_region(geometry,"x",[1])
    Radiant.set_boundary_conditions(geometry,"x-","void")
    Radiant.set_boundary_conditions(geometry,"x+","void")
    Radiant.set_material_per_region(geometry,[material])
    Radiant.build(geometry,cs)
    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,"CSD")
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
    # A prescribed postsolve observation fixture, not a numerical transport run.
    sources = Fixed_Sources(cs,geometry,solvers)
    sources.is_build = true
    sources.particles = [proton]
    sources.number_of_particles = 1
    sources.normalization_factor = 2.0
    moments = zeros(2,2,1,1,1,1)
    moments[:,1,1,1,1,1] .= [3.0,5.0]
    moments[:,2,1,1,1,1] .= -100.0
    cutoff = zeros(2,1,1,1,1)
    cutoff[1,1,1,1,1] = 8.0
    cutoff[2,1,1,1,1] = -100.0
    fpp = Radiant.Flux_Per_Particle(proton)
    Radiant.add_flux(fpp,moments)
    Radiant.add_flux_cutoff(fpp,cutoff)
    flux = Radiant.Flux()
    Radiant.add_flux(flux,fpp)
    score = proton_energy_accounting(binding,geometry,solvers,sources,flux)
    @test score.cutoff_particles == 6.0 # 3*8/(4*2)*2cm
    @test score.cutoff_kinetic_handoff_MeV == 6.0
    @test score.electronic_deposition_MeV == 16.0 # 2*(3+5)/2*2cm
    @test score.recoil_handoff_MeV == 8.0
    @test vec(Radiant.energy_deposition(cs,geometry,solvers,sources,flux,[proton])) == [4.0]
    @test score.electronic_deposition_per_mass[1] == 4.0
    @test !score.physical_validation
    fpp.flux_cutoff[1][1,1,1,1,1] = -1.0
    @test_throws ErrorException proton_energy_accounting(binding,geometry,solvers,sources,flux)
    fpp.flux_cutoff[1][1,1,1,1,1] = 8.0
    nonelastic_binding = bind_proton_native(proton,[data],[Proton_Nonelastic_Data(
        material_id="accounting-manufactured",energy_MeV=[1.0,10.0],
        removal_cm_inv=[1.0,1.0],event_family_ids=["synthetic-event"],source_sha256=hash)],
        [10.0,5.0,1.0])
    Radiant.build(nonelastic_binding.cross_sections)
    @test_throws ErrorException proton_energy_accounting(nonelastic_binding,geometry,solvers,sources,flux)
    physical_binding = bind_proton_native(proton,[data],[nonelastic],[10.0,5.0,1.0];
        qualification_mode=:physical)
    Radiant.build(physical_binding.cross_sections)
    @test_throws ErrorException proton_energy_accounting(physical_binding,geometry,solvers,sources,flux)

    # Native runtime convention check: conservative constant-loss CSD, two
    # unequal energy groups and void faces. This is not an energy leakage score.
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,0)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DG",1)
    Radiant.set_scheme(solver,"E","DG",1)
    Radiant.set_maximum_iteration(solver,100)
    Radiant.set_convergence_criterion(solver,1.0e-10)
    live_source = Anisotropic_Volume_Source(proton,[1],[2.0],
        [1.0e6,5.0e6,10.0e6],:isotropic,reshape([0.0,0.5],1,2,1),
        Source_Normalization(basis=:per_source_particle,source_hash="manufactured-csd-accounting"))
    live_sources = Fixed_Sources(cs,geometry,solvers)
    Radiant.add_source(live_sources,live_source)
    Radiant.build(live_sources)
    live_flux = Radiant.transport(cs,geometry,solvers,live_sources;retain_boundary_flux=true)
    live_score = proton_energy_accounting(binding,geometry,solvers,live_sources,live_flux)
    boundary_receipts = Radiant.get_boundary_flux(live_flux,proton)
    @test length(boundary_receipts) == 1
    @test all(row.converged && row.iterations <= row.iteration_cap &&
        isfinite(row.reconstruction_residual) && row.reconstruction_residual <= row.tolerance
        for row in boundary_receipts[1].convergence)
    outgoing = Radiant.get_outgoing_current(live_flux,proton) ./
        Radiant.get_normalization_factor(live_sources)
    @test all(outgoing .>= -1.0e-12)
    injected_particles = sum(get_volume_source_rate(live_source)) # 0.5/cm * 2cm = 1
    @test injected_particles == 1.0
    @test live_score.cutoff_particles > 0.0
    @test isapprox(sum(outgoing)+live_score.cutoff_particles,injected_particles;atol=1.0e-8,rtol=0.0)
    @test proton_energy_accounting_report(live_score)["classification"] == "BLOCKED_INCOMPLETE_LEDGER"
    println("CSD_ACCOUNTING injected=",injected_particles," outgoing=",sum(outgoing),
        " cutoff=",live_score.cutoff_particles," residual=",
        injected_particles-sum(outgoing)-live_score.cutoff_particles)
    println("CSD_ENERGY electronic=",live_score.electronic_deposition_MeV,
        " recoil_handoff=",live_score.recoil_handoff_MeV,
        " cutoff_handoff=",live_score.cutoff_kinetic_handoff_MeV,
        " escaped_energy=MISSING secondary_energy=MISSING physical_validation=false")
    for (group,status) in enumerate(boundary_receipts[1].convergence)
        println("CSD_CONVERGENCE group=",group," receipt=",status)
    end
    sources.number_of_particles = 2
    @test_throws ErrorException proton_energy_accounting(binding,geometry,solvers,sources,flux)
end
