using Test
using SHA

# Declared software fixture: only constant electronic/recoil stopping exists.
# Secondary and nonelastic energy transfers are exactly zero by construction.
# Injection is the represented DG source polynomial, not a monoenergy label.
function _complete_proton_ledger_fixture(;coupled=false,energy_order=2,source_slope=0.1,incident_group=1)
    proton = Proton()
    material = Material("complete-ledger-manufactured")
    Radiant.set_density(material,2.0)
    hash = bytes2hex(sha256("complete-ledger-constant-stopping-zero-production"))
    model = Tabulated_Ion_Transport_Model(species_id="proton",material_id="complete-ledger-manufactured",
        energy_MeV=[1.0,10.0],electronic_stopping_MeV_cm=[2.0,2.0],
        nuclear_stopping_MeV_cm=[1.0,1.0],energy_straggling_variance_MeV2_cm=[0.0,0.0],
        angular_variance_rad2_cm=[0.0,0.0],data_hash=hash,qualification_status=:synthetic)
    data = Proton_Material_Data(material,model;material_state="solid",density_g_cm3=2.0,source_sha256=hash)
    nonelastic = Proton_Nonelastic_Data(material_id="complete-ledger-manufactured",
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
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,0)
    Radiant.set_angular_boltzmann(solver,"standard")
    Radiant.set_scheme(solver,"x","DG",1)
    Radiant.set_scheme(solver,"E","DG",energy_order)
    Radiant.set_is_full_coupling(solver,coupled)
    Radiant.set_maximum_iteration(solver,100)
    Radiant.set_convergence_criterion(solver,1.0e-10)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
    source_values = incident_group == 1 ? [0.0,0.5] : [0.5,0.0]
    source = Anisotropic_Volume_Source(proton,[1],[2.0],[1.0e6,5.0e6,10.0e6],
        :isotropic,reshape(source_values,1,2,1),
        Source_Normalization(basis=:per_source_particle,source_hash="manufactured-represented-dg-energy"))
    sources = Fixed_Sources(cs,geometry,solvers)
    Radiant.add_source(sources,source)
    Radiant.build(sources)
    # Explicit native polynomial fixture augmentation, after the public constant
    # group source is projected. .5 + sqrt(3)*.1*u is nonnegative on [-1,1].
    if energy_order >= 2
        Radiant.get_source(sources,proton).volume_sources[incident_group,1,2,1,1,1] = source_slope
    end
    flux = Radiant.transport(cs,geometry,solvers,sources;retain_boundary_flux=true)
    score = proton_energy_accounting(binding,geometry,solvers,sources,flux)
    # Independent Simpson integration of E*q(E) for the known source polynomial.
    upper,lower = incident_group == 1 ? (10.0,5.0) : (5.0,1.0)
    slope = energy_order >= 2 ? source_slope : 0.0
    source_at(E) = (0.5+sqrt(3)*slope*(2*E-upper-lower)/(upper-lower))/(upper-lower)
    middle = (upper+lower)/2
    injection = 2.0*(upper-lower)/6*(lower*source_at(lower)+4middle*source_at(middle)+upper*source_at(upper))
    return (;binding,geometry,solvers,sources,flux,proton,score,injection,source_slope=slope)
end

@testset "Native DG represented energy closes the manufactured ledger" begin
    for coupled in (false,true)
        f = _complete_proton_ledger_fixture(;coupled)
        captured = only(Radiant.get_boundary_flux(f.flux,f.proton))
        @test captured.orders == [1,1,1,2]
        @test captured.fully_coupled == coupled && captured.is_csd
        @test captured.outer_convergence[].converged
        @test f.score.convergence_verified === true
        @test f.score.escape_status == :verified_native_capture
        @test f.score.escaped_energy_MeV > 0
        # Independent DG weak equations in one cell, using angular q=scalar q/2.
        # μ=±1/sqrt(3), S=3, incoming x flux=0, upper energy incoming=0.
        # Lower energy trace is a-sqrt(3)*b; its stored coefficient scales with
        # the next group width before being used as an incoming boundary value.
        incoming = 0.0
        h = inv(sqrt(3.0))/2.0
        for (group,width) in enumerate((5.0,4.0))
            s = 3.0/width
            q0 = group == 1 ? 0.25 : 0.0
            q1 = group == 1 ? f.source_slope/2 : 0.0
            # Solve the independent 2x2 system by its determinant, not native helpers.
            A,B,C,D = h+s,-sqrt(3)*s,sqrt(3)*s,h+3s
            r0,r1 = q0+s*incoming,q1+sqrt(3)*s*incoming
            a = (D*r0-B*r1)/(A*D-B*C)
            b = (A*r1-C*r0)/(A*D-B*C)
            for face in 1:2
                n = findfirst(c -> (isodd(face) ? -c : c) > 0,captured.directions[1])
                @test captured.faces[face][group,n,1,1,1] ≈ a atol=1e-11
                @test captured.faces[face][group,n,2,1,1] ≈ b atol=1e-11
                # Independently integrate the native reconstructed differential trace.
                upper,lower = captured.energy_boundaries[group:group+1]
                middle = (upper+lower)/2
                trace(E) = (a+sqrt(3)*b*(2*E-upper-lower)/width)/width
                integrated = width/6*(lower*trace(lower)+4middle*trace(middle)+upper*trace(upper))
                expected = captured.weights[n]*abs(captured.directions[1][n])*integrated
                @test Radiant.get_outgoing_energy_current(captured)[face,group] ≈ expected atol=1e-11
            end
            incoming = (a-sqrt(3)*b)*(group == 1 ? 4.0/5.0 : 1.0)
        end
        midpoint_only = sum(Radiant.get_outgoing_current(captured) .* reshape([7.5,3.0],1,2))
        @test abs(midpoint_only-f.score.escaped_energy_MeV) > 1e-5
        report = proton_energy_accounting_report(f.score;injected_energy_MeV=f.injection,
            secondary_transfer_MeV=0.0,other_transfer_MeV=0.0,closure_atol_MeV=1e-9)
        @test report["classification"] == "SOFTWARE_ARITHMETIC_CLOSED"
        @test abs(report["energy_residual_MeV"]) < 1e-9
        @test !report["physical_validation"] && !report["full_stopping_qualified"]
        @test proton_energy_accounting_report(f.score;injected_energy_MeV=f.injection)["classification"] ==
            "BLOCKED_INCOMPLETE_LEDGER" # undeclared production remains unknown
        scaled = proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux;
            physical_normalization=Source_Normalization(basis=:per_source_particle,source_rate_per_s=2.0,symmetry_factor=3.0))
        for key in (:escaped_energy_MeV,:electronic_deposition_MeV,:recoil_handoff_MeV,:cutoff_kinetic_handoff_MeV)
            @test getproperty(scaled,key) ≈ 6*getproperty(f.score,key)
        end
        @test scaled.source_basis == :physical_rate
        # A failed or nonfinite receipt cannot become a normalized ledger escape.
        original = captured.convergence[1]
        for replacement in ((converged=false,), (solver_residual=NaN,),
            (reconstruction_residual=2*original.tolerance,))
            captured.convergence[1] = merge(original,replacement)
            blocked = proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
            @test ismissing(blocked.escaped_energy_MeV)
            @test blocked.convergence_verified === false
        end
        captured.convergence[1] = original
        outer = captured.outer_convergence[]
        captured.outer_convergence[] = merge(outer,(converged=false,))
        @test ismissing(proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux).escaped_energy_MeV)
        captured.outer_convergence[] = outer
        captured.orders[4] = 1
        @test_throws ErrorException proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
        captured.orders[4] = 2
        saved_faces = deepcopy(captured.faces)
        # A negative void face/group must remain blocked even when the other
        # face offsets it enough to make the all-face total positive.
        for (face,value) in ((1,-1.0),(2,3.0))
            n = findfirst(c -> (isodd(face) ? -c : c) > 0,captured.directions[1])
            captured.faces[face][1,n,1,1,1] = value
            captured.faces[face][1,n,2,1,1] = 0.0
        end
        @test sum(Radiant.get_escaped_energy_current(captured)) > 0
        adverse = proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
        @test ismissing(adverse.escaped_energy_MeV)
        @test adverse.escape_status == :negative_integrated_void_energy
        @test minimum(adverse.raw_escaped_energy_by_face_group_MeV) < 0
        for face in eachindex(captured.faces); captured.faces[face] .= saved_faces[face]; end
        # Compatible-generation cancellation must also fail: generation 1 has
        # negative escape on x-, generation 2 offsets that same face/group.
        fpp = only(f.flux.flux_per_particle)
        second_capture = deepcopy(captured)
        for (data,value) in ((captured,-1.0),(second_capture,3.0))
            n = findfirst(c -> c < 0,data.directions[1])
            data.faces[1][1,n,1,1,1] = value
            data.faces[1][1,n,2,1,1] = 0.0
        end
        push!(fpp.flux,deepcopy(fpp.flux[1]))
        push!(fpp.flux_cutoff,deepcopy(fpp.flux_cutoff[1]))
        push!(fpp.boundary_flux,second_capture)
        @test minimum(Radiant.get_escaped_energy_current(fpp)) >= 0
        cancelled = proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
        @test ismissing(cancelled.escaped_energy_MeV)
        @test cancelled.escape_status == :negative_integrated_void_energy
        @test minimum(cancelled.raw_escaped_energy_by_generation_MeV[1]) < 0
        @test minimum(cancelled.raw_escaped_energy_by_face_group_MeV) >= 0
        pop!(fpp.flux); pop!(fpp.flux_cutoff); pop!(fpp.boundary_flux)
        for face in eachindex(captured.faces); captured.faces[face] .= saved_faces[face]; end
        println("MANUFACTURED_COMPLETE_LEDGER coupled=",coupled," injected=",f.injection,
            " electronic=",f.score.electronic_deposition_MeV," recoil=",f.score.recoil_handoff_MeV,
            " cutoff=",f.score.cutoff_kinetic_handoff_MeV," void_escape=",f.score.escaped_energy_MeV,
            " residual=",report["energy_residual_MeV"]," physical_validation=false")
    end
    empty = _complete_proton_ledger_fixture(;incident_group=2)
    capture = only(Radiant.get_boundary_flux(empty.flux,empty.proton))
    @test capture.convergence[1].terminal == :zero_source
    @test capture.convergence[1].iterations == 0
    @test empty.score.convergence_verified === true
    @test empty.score.escaped_energy_by_face_group_MeV[:,1] == zeros(2)
    @test proton_energy_accounting_report(empty.score;injected_energy_MeV=empty.injection,
        secondary_transfer_MeV=0.0,other_transfer_MeV=0.0,closure_atol_MeV=1e-9)["classification"] ==
        "SOFTWARE_ARITHMETIC_CLOSED"
end
