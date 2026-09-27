# Public source counterpart of the existing independently integrated DG2 fixture.
# Original private-source tests remain unchanged. This file executes actual SN
# transport in shifted unequal 2D/3D boxes; reference integration is Simpson's
# rule and an independent two-by-two weak-system solution. Synthetic only.
using Test, SHA

# Additive native normalization controls. Every spatial order is one, so the
# reduced and fully coupled spaces contain the same two pure energy modes.
# No mixed spatial/energy approximation or physical material accuracy is claimed.
function _public_multidimensional_proton_ledger(dimension,coupled,upper_energy;coarse_negative=false,account=true)
    proton = Proton()
    material = Material("multidimensional-ledger-manufactured")
    Radiant.set_density(material,2.0)
    hash = bytes2hex(sha256("multidimensional-ledger-Sel2-Srecoil1-zero-production-domain1-500"))
    model = Tabulated_Ion_Transport_Model(species_id="proton",material_id="multidimensional-ledger-manufactured",
        energy_MeV=[1.0,500.0],electronic_stopping_MeV_cm=[2.0,2.0],
        nuclear_stopping_MeV_cm=[1.0,1.0],energy_straggling_variance_MeV2_cm=[0.0,0.0],
        angular_variance_rad2_cm=[0.0,0.0],data_hash=hash,qualification_status=:synthetic)
    data = Proton_Material_Data(material,model;material_state="solid",density_g_cm3=2.0,source_sha256=hash)
    nonelastic = Proton_Nonelastic_Data(material_id="multidimensional-ledger-manufactured",
        energy_MeV=[1.0,500.0],removal_cm_inv=[0.0,0.0],source_sha256=hash)
    lower_width = coarse_negative ? 0.35*(upper_energy-1.0) : min(0.35*(upper_energy-1.0),1.0)
    edges = [Float64(upper_energy),1.0+lower_width,1.0]
    binding = bind_proton_native(proton,[data],[nonelastic],edges)
    cs = binding.cross_sections
    Radiant.build(cs)
    widths = [0.7,1.3,2.1]
    origins = [-0.4,1.7,-2.3]
    geometry = Geometry()
    Radiant.set_dimension(geometry,dimension)
    for a in 1:dimension
        axis = ("x","y","z")[a]
        Radiant.set_number_of_regions(geometry,axis,1)
        Radiant.set_region_boundaries(geometry,axis,[origins[a],origins[a]+widths[a]])
        Radiant.set_voxels_per_region(geometry,axis,[1])
        Radiant.set_boundary_conditions(geometry,axis*"-","void")
        Radiant.set_boundary_conditions(geometry,axis*"+","void")
    end
    Radiant.set_material_per_region(geometry,reshape([material],ntuple(_ -> 1,dimension)))
    Radiant.build(geometry,cs)
    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,"CSD")
    Radiant.set_quadrature(solver,"gauss-legendre-chebychev",2,3)
    Radiant.set_legendre_order(solver,0)
    Radiant.set_angular_boltzmann(solver,"standard")
    for axis in ("x","y","z")[1:dimension]
        Radiant.set_scheme(solver,axis,"DG",1)
    end
    Radiant.set_scheme(solver,"E","DG",2)
    Radiant.set_is_full_coupling(solver,coupled)
    Radiant.set_maximum_iteration(solver,100)
    Radiant.set_convergence_criterion(solver,1.0e-10)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
    volume = prod(widths[1:dimension])
    scalar_source = 1.0/volume
    scalar_slope = 0.2/volume
    source = Anisotropic_Volume_Source(proton,[1],[volume],reverse(edges).*1.0e6,
        :isotropic,reshape([0.0,scalar_source],1,2,1),
        Source_Normalization(basis=:per_source_particle,source_hash="multidimensional-represented-polynomial"))
    # Public API: ascending source group 2 is the upper native group.
    # Keep signed Q1; its differential polynomial is nonnegative on [-1,1].
    moments = zeros(Float64,1,2,1,2)
    moments[:,:,:,1] .= source.values
    moments[1,2,1,2] = scalar_slope
    source = Energy_Moment_Volume_Source(source,moments)
    sources = Fixed_Sources(cs,geometry,solvers)
    Radiant.add_source(sources,source)
    Radiant.build(sources)
    flux = Radiant.transport(cs,geometry,solvers,sources;retain_boundary_flux=true)
    capture = only(Radiant.get_boundary_flux(flux,proton))
    raw_cutoff = Radiant.get_flux_cutoff(flux,proton)[1,1,1,1,1]
    println("PUBLIC_MULTIDIMENSIONAL_RAW_CUTOFF dimension=",dimension," coupled=",coupled,
        " spectrum_upper_MeV=",upper_energy," coarse_negative=",coarse_negative,
        " edges_MeV=",edges," scalar_cutoff_coefficient=",raw_cutoff,
        " groups=",capture.convergence," outer=",capture.outer_convergence[])
    score = account ? proton_energy_accounting(binding,geometry,solvers,sources,flux) : nothing
    return (;proton,binding,geometry,solvers,sources,source,flux,capture,score,
        volume,widths,edges,scalar_source,scalar_slope)
end

# Independent diagnostic for the retained original negative-cutoff control.
# With no source in the lower group, its outgoing energy trace is exactly
# 2*s*incoming*(3*s-h)/determinant: convergence does not imply positivity.
function _public_multidimensional_cutoff_reference(f)
    data = f.capture
    traces = Float64[]
    analytic_traces = Float64[]
    for n in eachindex(data.weights)
        h = sum(abs(data.directions[a][n])/f.widths[a] for a in 1:data.dimension)
        incoming = 0.0
        for g in 1:2
            width = f.edges[g]-f.edges[g+1]
            s = 3.0/width
            q0 = g == 1 ? f.scalar_source/sum(data.weights) : 0.0
            q1 = g == 1 ? f.scalar_slope/sum(data.weights) : 0.0
            A,B,C,D = h+s,-sqrt(3.0)*s,sqrt(3.0)*s,h+3.0*s
            determinant = A*D-B*C
            r0,r1 = q0+s*incoming,q1+sqrt(3.0)*s*incoming
            a = (D*r0-B*r1)/determinant
            b = (A*r1-C*r0)/determinant
            edge = a-sqrt(3.0)*b
            if g == 1
                incoming = edge*(f.edges[2]-f.edges[3])/width
            else
                push!(traces,edge)
                push!(analytic_traces,2.0*s*incoming*(3.0*s-h)/determinant)
            end
        end
    end
    return (;traces,analytic_traces,scalar=sum(data.weights .* traces))
end

function _check_public_multidimensional_proton_ledger(dimension,coupled,upper_energy)
    f = _public_multidimensional_proton_ledger(dimension,coupled,upper_energy)
    data = f.capture
    @test data.dimension == dimension
    @test data.orders == [1,1,1,2]
    @test data.fully_coupled == coupled && data.is_csd
    @test data.basis == :native_group_integrated_normalized_legendre
    @test all(size(face,3) == 2 for face in data.faces)
    # Cauchy-Schwarz bounds the directional streaming rate for all unit vectors.
    h_bound = sqrt(sum(inv(f.widths[a])^2 for a in 1:dimension))
    @test h_bound < 3.0*(3.0/(f.edges[2]-f.edges[3]))
    @test all(sum(data.directions[a][n]^2 for a in 1:3) ≈ 1.0 for n in eachindex(data.weights))
    # Native GL-cheby weights integrate the full unit sphere. Angular source
    # density is scalar/sum(weights), independent of any implicit 4pi factor.
    angular_measure = sum(data.weights)
    @test angular_measure ≈ 4pi atol=1e-12
    @test Radiant.get_normalization_factor(f.sources) == 1.0
    @test sum(get_volume_source_rate(f.source)) ≈ 1.0 atol=1e-12
    @test f.score.convergence_verified === true
    @test f.score.escape_status == :verified_native_capture
    @test all(r.converged && isfinite(r.solver_residual) && r.solver_residual <= r.tolerance &&
        isfinite(r.reconstruction_residual) && r.reconstruction_residual <= r.tolerance for r in data.convergence)

    exact_faces = zeros(2dimension,2)
    cutoff_coefficient = 0.0
    for n in eachindex(data.weights)
        incoming = 0.0
        h = sum(abs(data.directions[a][n])/f.widths[a] for a in 1:dimension)
        for g in 1:2
            upper,lower = f.edges[g:g+1]
            width = upper-lower
            s = 3.0/width
            q0 = g == 1 ? f.scalar_source/angular_measure : 0.0
            q1 = g == 1 ? f.scalar_slope/angular_measure : 0.0
            # Independent weak DG moment equations. Constant stopping and
            # no incoming spatial faces make streaming h times each mode.
            A,B,C,D = h+s,-sqrt(3.0)*s,sqrt(3.0)*s,h+3.0*s
            r0,r1 = q0+s*incoming,q1+sqrt(3.0)*s*incoming
            determinant = A*D-B*C
            a = (D*r0-B*r1)/determinant
            b = (A*r1-C*r0)/determinant
            midpoint = (upper+lower)/2
            trace_at(E) = (a+sqrt(3.0)*b*(2*E-upper-lower)/width)/width
            # Polynomial projection and kinetic-energy integral by independent
            # Simpson rule; no scorer moment helper or energy contraction used.
            projected0 = width/6*(trace_at(lower)+4*trace_at(midpoint)+trace_at(upper))
            projected1 = width/6*(-sqrt(3.0)*trace_at(lower)+sqrt(3.0)*trace_at(upper))
            energy_integral = width/6*(lower*trace_at(lower)+4*midpoint*trace_at(midpoint)+upper*trace_at(upper))
            @test projected0 ≈ a atol=1e-12 rtol=1e-11
            @test projected1 ≈ b atol=1e-12 rtol=1e-11
            for axis in 1:dimension
                cosine = data.directions[axis][n]
                cosine == 0 && continue
                face = 2axis-(cosine < 0 ? 1 : 0)
                @test data.faces[face][g,n,1,1,1] ≈ projected0 atol=1e-11 rtol=1e-11
                @test data.faces[face][g,n,2,1,1] ≈ projected1 atol=1e-11 rtol=1e-11
                measure = f.volume/f.widths[axis]
                exact_faces[face,g] += data.weights[n]*abs(cosine)*measure*energy_integral
            end
            outgoing_energy_edge = a-sqrt(3.0)*b
            if g == 1
                incoming = outgoing_energy_edge*(f.edges[2]-f.edges[3])/width
            else
                cutoff_coefficient += data.weights[n]*outgoing_energy_edge
            end
        end
    end
    @test Radiant.get_outgoing_energy_current(data) ≈ exact_faces atol=1e-10 rtol=1e-11
    @test Radiant.get_escaped_energy_current(data) == Radiant.get_outgoing_energy_current(data)
    @test f.score.escaped_energy_by_face_group_MeV ≈ exact_faces atol=1e-10 rtol=1e-11
    native_cutoff = Radiant.get_flux_cutoff(f.flux,f.proton)[1,1,1,1,1]
    @test native_cutoff ≈ cutoff_coefficient atol=1e-11 rtol=1e-11
    cutoff_particles = 3.0*cutoff_coefficient/(f.edges[2]-f.edges[3])*f.volume
    @test f.score.cutoff_particles ≈ cutoff_particles atol=1e-10 rtol=1e-11
    @test f.score.cutoff_kinetic_handoff_MeV ≈ f.edges[3]*cutoff_particles atol=1e-10 rtol=1e-11

    # Injection independently integrates the declared nonnegative source
    # polynomial and cell volume. Eupper is a spectrum endpoint, not monoenergy.
    upper,lower = f.edges[1:2]
    width = upper-lower
    midpoint = (upper+lower)/2
    source_at(E) = (f.scalar_source+sqrt(3.0)*f.scalar_slope*(2*E-upper-lower)/width)/width
    injection = f.volume*width/6*(lower*source_at(lower)+4*midpoint*source_at(midpoint)+upper*source_at(upper))
    declared_atol,declared_rtol = 1e-9,1e-11
    report = proton_energy_accounting_report(f.score;injected_energy_MeV=injection,
        secondary_transfer_MeV=0.0,other_transfer_MeV=0.0,
        closure_atol_MeV=declared_atol,closure_rtol=declared_rtol)
    println("PUBLIC_MULTIDIMENSIONAL_LEDGER dimension=",dimension," coupled=",coupled,
        " spectrum_upper_MeV=",upper_energy," represented_injection_MeV=",injection,
        " electronic_MeV=",f.score.electronic_deposition_MeV," recoil_MeV=",f.score.recoil_handoff_MeV,
        " cutoff_MeV=",f.score.cutoff_kinetic_handoff_MeV," void_escape_MeV=",f.score.escaped_energy_MeV,
        " residual_MeV=",report["energy_residual_MeV"]," physical_validation=false")
    println("PUBLIC_MULTIDIMENSIONAL_CONVERGENCE dimension=",dimension," coupled=",coupled,
        " spectrum_upper_MeV=",upper_energy," groups=",data.convergence," outer=",data.outer_convergence[])
    @test report["classification"] == "SOFTWARE_ARITHMETIC_CLOSED"
    @test abs(report["energy_residual_MeV"]) <= declared_atol+declared_rtol*injection
    @test !report["physical_validation"] && !report["full_stopping_qualified"]
    @test proton_energy_accounting_report(f.score;injected_energy_MeV=injection)["classification"] ==
        "BLOCKED_INCOMPLETE_LEDGER"
end

# Base normalization controls must pass before the wider endpoint sweep runs.
@testset "Public source native 2D and 3D DG2 ledger base controls" begin
    for dimension in (2,3), coupled in (false,true)
        _check_public_multidimensional_proton_ledger(dimension,coupled,20.0)
    end
end

@testset "Public source coarse multidimensional DG2 negative cutoff is rejected" begin
    f = _public_multidimensional_proton_ledger(3,false,20.0;coarse_negative=true,account=false)
    reference = _public_multidimensional_cutoff_reference(f)
    raw_cutoff = Radiant.get_flux_cutoff(f.flux,f.proton)[1,1,1,1,1]
    println("PUBLIC_MULTIDIMENSIONAL_NEGATIVE_CUTOFF independent_ordinate_traces=",reference.traces,
        " independent_scalar=",reference.scalar," native_scalar=",raw_cutoff,
        " physical_validation=false expected_accounting_rejection=true")
    @test all(row.converged for row in f.capture.convergence)
    @test reference.traces ≈ reference.analytic_traces atol=1e-12 rtol=1e-11
    @test minimum(reference.traces) < 0.0
    @test reference.scalar < 0.0
    @test raw_cutoff ≈ reference.scalar atol=1e-11 rtol=1e-11
    @test raw_cutoff < 0.0
    @test_throws ErrorException proton_energy_accounting(f.binding,f.geometry,f.solvers,f.sources,f.flux)
end

@testset "Public source native low-energy and upper-endpoint geometry controls" begin
    for upper_energy in (1.001,1.01,1.1,1.5,3.0,5.0,100.0,200.0,499.99,500.0), dimension in (2,3), coupled in (false,true)
        _check_public_multidimensional_proton_ledger(dimension,coupled,upper_energy)
    end
end
