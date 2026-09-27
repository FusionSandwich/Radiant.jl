using Test
using SHA

# Independent manufactured differential scalar flux phi(E)=10-E in [1,10].
# One spatial DG1 cell, length 2, GL2 directions: h=|mu|/length. The prescribed
# scalar source q=h*phi-d(S_h*phi)/dE satisfies the native weak equations with
# zero incoming x and upper-energy values. This is a discretization fixture,
# not a continuum-space solution or a physical material.
function _discrete_energy_fixture(;oe=2,coupled=true,stopping=:constant,bounds=[10.0,1.0])
    println("CSD_WEAK_CASE_BEGIN oe=",oe," coupled=",coupled," stopping=",stopping," bounds_MeV=",bounds)
    upper,lower,length_cm = 10.0,1.0,2.0
    h = inv(sqrt(3.0))/length_cm
    slope = stopping == :affine ? 0.1 : 0.0
    native_s(E) = 2.0+slope*(E-lower)
    phi(E) = upper-E
    q(E) = h*phi(E)-slope*phi(E)+native_s(E)
    simpson(f) = (upper-lower)/6*(f(lower)+4*f((upper+lower)/2)+f(upper))
    group_integral(f,a,b) = (a-b)/6*(f(b)+4*f((a+b)/2)+f(a))
    q0_groups = [group_integral(q,bounds[g],bounds[g+1]) for g in 1:length(bounds)-1]
    q1_groups = [group_integral(E -> q(E)*sqrt(3.0)*(2*E-bounds[g]-bounds[g+1])/
        (bounds[g]-bounds[g+1]),bounds[g],bounds[g+1]) for g in 1:length(bounds)-1]
    q0 = sum(q0_groups)
    injection = length_cm*simpson(E -> E*q(E))
    expected_loss = length_cm*simpson(E -> native_s(E)*phi(E))
    expected_escape = inv(sqrt(3.0))*simpson(E -> E*phi(E))
    expected_cutoff = length_cm*lower*native_s(lower)*phi(lower)
    proton = Proton()
    material = Material("weak-balance-manufactured")
    Radiant.set_density(material,2.0)
    hash = bytes2hex(sha256("weak-balance-manufactured-$(stopping)-no-production"))
    energies = stopping == :midpoint_bump ? [1.0,5.5,10.0] : [1.0,10.0]
    values = stopping == :midpoint_bump ? [2.0,4.0,2.0] : native_s.(energies)
    model = Tabulated_Ion_Transport_Model(species_id="proton",material_id="weak-balance-manufactured",
        energy_MeV=energies,electronic_stopping_MeV_cm=values,
        nuclear_stopping_MeV_cm=zeros(length(energies)),
        energy_straggling_variance_MeV2_cm=zeros(length(energies)),
        angular_variance_rad2_cm=zeros(length(energies)),
        data_hash=hash,qualification_status=:synthetic)
    data = Proton_Material_Data(material,model;material_state="solid",density_g_cm3=2.0,source_sha256=hash)
    nonelastic = Proton_Nonelastic_Data(material_id="weak-balance-manufactured",
        energy_MeV=[1.0,10.0],removal_cm_inv=[0.0,0.0],source_sha256=hash)
    binding = bind_proton_native(proton,[data],[nonelastic],bounds)
    cs = binding.cross_sections
    Radiant.build(cs)
    geometry = Geometry()
    Radiant.set_dimension(geometry,1)
    Radiant.set_number_of_regions(geometry,"x",1)
    Radiant.set_region_boundaries(geometry,"x",[-0.7,1.3])
    Radiant.set_voxels_per_region(geometry,"x",[1])
    Radiant.set_boundary_conditions(geometry,"x-","void")
    Radiant.set_boundary_conditions(geometry,"x+","void")
    Radiant.set_material_per_region(geometry,[material])
    Radiant.build(geometry,cs)
    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,"CSD")
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,0)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DG",1)
    Radiant.set_scheme(solver,"E","DG",oe)
    Radiant.set_is_full_coupling(solver,coupled)
    Radiant.set_maximum_iteration(solver,100)
    Radiant.set_convergence_criterion(solver,1e-11)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
    source = Anisotropic_Volume_Source(proton,[1],[length_cm],reverse(bounds).*1e6,
        :isotropic,reshape(reverse(q0_groups),1,length(q0_groups),1),
        Source_Normalization(basis=:per_source_particle,source_hash=hash))
    sources = Fixed_Sources(cs,geometry,solvers)
    Radiant.add_source(sources,source)
    Radiant.build(sources)
    # Explicit fixture-only source energy slope. q is affine, so coefficients
    # >=2 vanish exactly. This does not qualify public higher-energy source support.
    if oe >= 2
        for g in eachindex(q1_groups)
            Radiant.get_source(sources,proton).volume_sources[g,1,2,1,1,1] = q1_groups[g]
        end
    end
    flux = Radiant.transport(cs,geometry,solvers,sources;retain_boundary_flux=true)
    diagnostic = Radiant.proton_discrete_energy_balance(binding,geometry,solvers,sources,flux;
        injected_energy_MeV=injection)
    println("CSD_WEAK_CASE_RESULT oe=",oe," coupled=",coupled," stopping=",stopping,
        " bounds_MeV=",bounds," source_energy_MeV=",diagnostic.represented_source_energy_MeV,
        " stopping_volume_MeV=",diagnostic.stopping_volume_MeV,
        " expected_stopping_MeV=",expected_loss," escape_MeV=",diagnostic.outgoing_energy_MeV,
        " expected_escape_MeV=",expected_escape," cutoff_MeV=",diagnostic.cutoff_kinetic_energy_MeV,
        " expected_cutoff_MeV=",expected_cutoff," weak_residual_MeV=",diagnostic.weak_balance_residual_MeV,
        " native_scalar_energy_coefficients=",Radiant.get_flux(flux,proton)[:,1,1:oe,1,1,1])
    return (;binding,geometry,solvers,sources,flux,proton,solver,diagnostic,
        injection,expected_loss,expected_escape,expected_cutoff,q0,q1_groups)
end

@testset "Native CSD first-energy-moment weak balance" begin
    for coupled in (false,true), oe in (2,4), stopping in (:constant,:affine,:midpoint_bump)
        f = _discrete_energy_fixture(;oe,coupled,stopping)
        d = f.diagnostic
        @test d.convergence_verified && d.energy_test_represented
        @test d.classification == :signed_native_weak_balance_diagnostic
        @test d.orders == [1,1,1,oe] && d.fully_coupled == coupled
        @test d.represented_source_energy_MeV ≈ f.injection atol=1e-9 rtol=1e-11
        @test d.stopping_volume_MeV ≈ f.expected_loss atol=1e-9 rtol=1e-11
        @test d.outgoing_energy_MeV ≈ f.expected_escape atol=1e-9 rtol=1e-11
        @test d.cutoff_kinetic_energy_MeV ≈ f.expected_cutoff atol=1e-9 rtol=1e-11
        @test maximum(abs,d.particle_balance_residual_by_group) < 1e-9
        @test maximum(abs,d.tested_balance_residual_by_group_MeV) < 1e-9
        @test maximum(abs,d.stored_minus_reconstructed_cutoff) < 1e-9
        @test abs(d.weak_balance_residual_MeV) < 1e-9
        native_coefficients = Radiant.get_flux(f.flux,f.proton)[1,1,1:oe,1,1,1]
        expected_coefficients = vcat(40.5,-81.0/(2*sqrt(3.0)),zeros(oe-2))
        @test native_coefficients ≈ expected_coefficients atol=1e-9 rtol=1e-11
        @test !d.response_score_replaced && !d.strict_ledger_accepted
        @test !d.physical_validation && !d.publication_ready && !d.pr_ready
        # Independently integrated physical source must not be replaced by an
        # endpoint energy label times source particle count.
        wrong = Radiant.proton_discrete_energy_balance(f.binding,f.geometry,f.solvers,f.sources,f.flux;
            injected_energy_MeV=10.0*2.0*f.q0)
        @test abs(wrong.injection_projection_difference_MeV) > 1.0
        @test wrong.weak_balance_residual_MeV == d.weak_balance_residual_MeV
        # Default process response remains midpoint*Phi0, even when it differs
        # from the native weak stopping integral. No closure-forcing correction.
        score = proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
        @test score.electronic_deposition_MeV ≈ d.midpoint_response_MeV atol=1e-9
        report = proton_energy_accounting_report(score;injected_energy_MeV=f.injection,
            secondary_transfer_MeV=0.0,other_transfer_MeV=0.0,
            closure_atol_MeV=1e-9,closure_rtol=1e-11)
        if stopping == :constant
            @test abs(d.stopping_response_difference_MeV) < 1e-9
            @test report["classification"] == "SOFTWARE_ARITHMETIC_CLOSED"
        else
            @test abs(d.stopping_response_difference_MeV) > 1.0
            @test report["classification"] == "FAIL_ENERGY_CLOSURE"
            @test report["energy_residual_MeV"] ≈ d.stopping_response_difference_MeV atol=1e-9
        end
        if oe == 2 && coupled && stopping == :affine
            # Negative cutoff stays signed in the diagnostic and is rejected
            # by ordinary physical ownership scoring. A modified trace also
            # appears explicitly as a mismatch; convergence is not a waiver.
            cut = only(f.flux.flux_per_particle).flux_cutoff[1]
            saved = cut[1,1,1,1,1]
            cut[1,1,1,1,1] = -1.0
            adverse = Radiant.proton_discrete_energy_balance(f.binding,f.geometry,f.solvers,f.sources,f.flux)
            @test adverse.minimum_stored_scalar_cutoff == -1.0
            @test adverse.cutoff_kinetic_energy_MeV < 0.0
            @test abs(only(adverse.stored_minus_reconstructed_cutoff)) > 1.0
            @test abs(adverse.weak_balance_residual_MeV) > 1.0
            @test_throws ErrorException proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
            cut[1,1,1,1,1] = NaN
            @test_throws ErrorException Radiant.proton_discrete_energy_balance(f.binding,f.geometry,f.solvers,f.sources,f.flux)
            cut[1,1,1,1,1] = saved
            source = Radiant.get_source(f.sources,f.proton)
            saved_surface = source.surface_sources[1,1,1]
            source.surface_sources[1,1,1] = 1.0
            @test_throws ErrorException Radiant.proton_discrete_energy_balance(f.binding,f.geometry,f.solvers,f.sources,f.flux)
            source.surface_sources[1,1,1] = saved_surface
            Radiant.set_scheme(f.solver,"E","DD",2)
            @test_throws ErrorException Radiant.proton_discrete_energy_balance(f.binding,f.geometry,f.solvers,f.sources,f.flux)
            Radiant.set_scheme(f.solver,"E","DG",2)
        end
    end
    for coupled in (false,true)
        f = _discrete_energy_fixture(;oe=1,coupled,stopping=:constant)
        d = f.diagnostic
        @test !d.energy_test_represented
        @test maximum(abs,d.particle_balance_residual_by_group) < 1e-9
        @test abs(d.weak_balance_residual_MeV) < 1e-9
        # For one group with zero upper input, midpoint-weighted particle
        # balance yields (mid-Ecut)*S(Ecut)*Phi0/width = S*Phi0/2.
        @test d.dg1_midpoint_interface_loss_MeV ≈ d.stopping_volume_MeV/2 atol=1e-9
        @test abs(d.represented_midpoint_response_residual_MeV) > 1.0
        @test abs(d.injection_projection_difference_MeV) > 1.0
        @test !d.strict_ledger_accepted
    end
    # Unequal energy groups require native incoming trace width rescaling.
    # Interface flow is continuous and cancels at the physical boundary energy,
    # whereas a DG1 midpoint test uses the gap between group midpoints.
    for oe in (1,2,4)
        f = _discrete_energy_fixture(;oe,stopping=:affine,bounds=[10.0,6.0,1.0])
        d = f.diagnostic
        @test maximum(abs,d.particle_balance_residual_by_group) < 1e-9
        @test abs(d.weak_balance_residual_MeV) < 1e-9
        @test maximum(abs,d.stored_minus_reconstructed_cutoff) < 1e-9
        if oe >= 2
            @test d.lower_energy_particle_flow_by_group ≈ [2.0*2.5*4.0,2.0*2.0*9.0] atol=1e-9
            @test d.stopping_volume_MeV ≈ f.expected_loss atol=1e-9
            @test d.represented_source_energy_MeV ≈ f.injection atol=1e-9
        end
    end
end
