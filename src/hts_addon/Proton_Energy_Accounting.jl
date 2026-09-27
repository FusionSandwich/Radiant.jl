"""Finite, nonnegative scalar observables are required; negative moments are not clipped."""
function _proton_accounting_array(value, name; positive=false)
    data = Float64.(value)
    isempty(data) && error("$(name) cannot be empty.")
    all(x -> isfinite(x) && (positive ? x > 0.0 : x >= 0.0),data) ||
        error("$(name) must be finite and $(positive ? "positive" : "nonnegative").")
    return data
end

"""Finite signed diagnostics; a Legendre coefficient is not a positive observable."""
function _proton_balance_signed(value,name)
    data = Float64.(value)
    isempty(data) && error("$(name) cannot be empty.")
    all(isfinite,data) || error("$(name) must be finite.")
    return data
end

"""
Scheme-consistent arithmetic for a pure-CSD, DG, volume-source problem with zero
incoming energy flow at the upper domain boundary. Arrays `phi` and `q` contain
the spatial-average, scalar-angular, group-integrated normalized Legendre energy
coefficients, shaped (group, energy order, cell). All signed observations survive.
This internal kernel does not establish that its inputs came from a native solve.
"""
function _proton_discrete_csd_terms(;energy_boundaries_MeV,phi,q,
    stopping_at_boundaries_MeV_cm,midpoint_response_MeV_cm,cell_measure_cm_d,
    outgoing_particles_by_face_group,outgoing_energy_by_face_group_MeV,
    scalar_cutoff_flux,source_normalization_divisor::Real,
    injected_energy_MeV=missing)
    eb = _proton_balance_signed(energy_boundaries_MeV,"Energy boundaries")
    eb isa AbstractVector && length(eb) >= 2 && all(diff(eb) .< 0) &&
        1.0 <= eb[end] < eb[1] <= 500.0 || error("Require descending proton energy boundaries in [1,500].")
    p = _proton_balance_signed(phi,"Scalar flux moments")
    z = _proton_balance_signed(q,"Scalar source moments")
    ndims(p) == 3 && size(p) == size(z) || error("Flux/source must have identical (group,energy order,cell) shapes.")
    ng,oe,nx = size(p)
    ng == length(eb)-1 && oe >= 1 || error("Energy moment shape differs from group metadata.")
    sb = _proton_accounting_array(stopping_at_boundaries_MeV_cm,"Boundary stopping")
    sm = _proton_accounting_array(midpoint_response_MeV_cm,"Midpoint response")
    v = _proton_accounting_array(cell_measure_cm_d,"Cell measure";positive=true)
    size(sb) == (ng+1,nx) && size(sm) == (ng,nx) && size(v) == (nx,) ||
        error("Stopping and cell measure shapes differ from moment metadata.")
    fp = _proton_balance_signed(outgoing_particles_by_face_group,"Outward particles")
    fe = _proton_balance_signed(outgoing_energy_by_face_group_MeV,"Outward kinetic energy")
    size(fp) == size(fe) && ndims(fp) == 2 && size(fp,2) == ng ||
        error("Face/group scores must have identical shapes with the declared group count.")
    cut = _proton_balance_signed(scalar_cutoff_flux,"Stored scalar cutoff")
    size(cut) == (nx,) || error("Cutoff shape differs from cell count.")
    norm = _proton_accounting_scalar(source_normalization_divisor,"Source divisor";positive=true)
    injection = ismissing(injected_energy_MeV) ? missing :
        _proton_accounting_scalar(injected_energy_MeV,"Independent injection")
    widths = eb[1:end-1] .- eb[2:end]
    mids = (eb[1:end-1] .+ eb[2:end])./2
    lower_trace = zeros(ng,nx)
    source_energy = zeros(ng)
    source_particles = zeros(ng)
    volume_loss = zeros(ng)
    midpoint_loss = zeros(ng)
    energy_flow = zeros(ng)
    for g in 1:ng, x in 1:nx
        # Pure energy coefficients have spatial degree zero in both coupling modes.
        for k in 1:oe
            lower_trace[g,x] += (-1)^(k-1)*sqrt(2*k-1)*p[g,k,x]
        end
        factor = v[x]/norm
        source_particles[g] += factor*z[g,1,x]
        source_energy[g] += factor*(mids[g]*z[g,1,x] +
            (oe >= 2 ? widths[g]*z[g,2,x]/(2*sqrt(3.0)) : 0.0))
        s0 = (sb[g,x]+sb[g+1,x])/2
        s1 = (sb[g,x]-sb[g+1,x])/(2*sqrt(3.0))
        volume_loss[g] += factor*(s0*p[g,1,x] + (oe >= 2 ? s1*p[g,2,x] : 0.0))
        midpoint_loss[g] += factor*sm[g,x]*p[g,1,x]
        energy_flow[g] += factor*sb[g+1,x]*lower_trace[g,x]/widths[g]
    end
    stored_cutoff_particles = sum(v .* sb[end,:] .* cut)/(norm*widths[end])
    cutoff_energy = eb[end]*stored_cutoff_particles
    outgoing_p = vec(sum(fp;dims=1))./norm
    outgoing_e = vec(sum(fe;dims=1))./norm
    incoming_flow = vcat(0.0,energy_flow[1:end-1])
    particle_residual = source_particles-outgoing_p-energy_flow+incoming_flow
    # DG1 only has the constant test. Its midpoint-weighted group equations
    # telescope to this loss; it is not integral S*phi and is never called heat.
    midpoint_transfer = mids .* (energy_flow-incoming_flow)
    dg1_loss = sum(midpoint_transfer)-eb[end]*energy_flow[end]
    loss = oe == 1 ? dg1_loss : sum(volume_loss)
    tested_group_residual = oe == 1 ?
        mids .* source_particles-outgoing_e-midpoint_transfer :
        source_energy-outgoing_e-eb[2:end].*energy_flow+
            eb[1:end-1].*incoming_flow-volume_loss
    represented_injection = sum(source_energy)
    escape = sum(outgoing_e)
    residual = represented_injection-escape-cutoff_energy-loss
    legacy_residual = represented_injection-escape-cutoff_energy-sum(midpoint_loss)
    cutoff_difference = cut-lower_trace[end,:]
    all(isfinite,lower_trace) && all(isfinite,source_energy) &&
        all(isfinite,source_particles) && all(isfinite,volume_loss) &&
        all(isfinite,midpoint_loss) && all(isfinite,energy_flow) &&
        all(isfinite,particle_residual) && all(isfinite,tested_group_residual) &&
        all(isfinite,cutoff_difference) && all(isfinite,
            (stored_cutoff_particles,cutoff_energy,dg1_loss,loss,represented_injection,
             escape,residual,legacy_residual)) || error("Discrete energy balance overflowed.")
    independent_residual = ismissing(injection) ? missing :
        injection-escape-cutoff_energy-sum(midpoint_loss)
    projection_delta = ismissing(injection) ? missing : injection-represented_injection
    (ismissing(independent_residual) || isfinite(independent_residual)) &&
        (ismissing(projection_delta) || isfinite(projection_delta)) || error("Injection comparison overflowed.")
    return (energy_order=oe,energy_test_represented=oe >= 2,
        represented_source_energy_MeV=represented_injection,
        independently_injected_energy_MeV=injection,
        injection_projection_difference_MeV=projection_delta,
        source_particles_by_group=source_particles,
        particle_balance_residual_by_group=particle_residual,
        tested_balance_residual_by_group_MeV=tested_group_residual,
        stopping_volume_by_group_MeV=volume_loss,
        stopping_volume_MeV=sum(volume_loss),
        midpoint_response_by_group_MeV=midpoint_loss,
        midpoint_response_MeV=sum(midpoint_loss),
        stopping_response_difference_MeV=sum(volume_loss)-sum(midpoint_loss),
        dg1_midpoint_interface_loss_MeV=dg1_loss,
        scheme_balance_loss_MeV=loss,
        lower_energy_particle_flow_by_group=energy_flow,
        reconstructed_lower_energy_scalar_trace=lower_trace,
        stored_minus_reconstructed_cutoff=cutoff_difference,
        stored_cutoff_particles=stored_cutoff_particles,
        cutoff_kinetic_energy_MeV=cutoff_energy,
        outgoing_energy_MeV=escape,
        weak_balance_residual_MeV=residual,
        represented_midpoint_response_residual_MeV=legacy_residual,
        independent_midpoint_response_residual_MeV=independent_residual,
        minimum_stored_scalar_cutoff=minimum(cut),
        minimum_outgoing_face_group_energy_MeV=minimum(fe)/norm,
        response_score_replaced=false,strict_ledger_accepted=false,
        physical_validation=false,publication_ready=false,pr_ready=false)
end

"""
    proton_discrete_energy_balance(binding, geometry, solvers, sources, flux;
                                   injected_energy_MeV=missing)

Opt-in signed weak-balance diagnostic for one-generation, one-dimensional
Cartesian pure CSD with DG space/energy schemes, void faces and volume sources.
No physical rate scaling is applied; every term uses the declared transport
source basis and its common divisor. Compatible native capture and receipts are
required. A converged signed solution may still fail positivity and the strict
ledger. Neither operator loss nor DG1 interface loss replaces process deposition.
Boundary sources, fields, scattering, removal and other schemes are unsupported.
"""
function proton_discrete_energy_balance(binding::Proton_Native_Binding,geometry::Geometry,
    solvers::Solvers,sources::Fixed_Sources,flux::Flux;
    injected_energy_MeV=missing,field::Electromagnetic_Field=Electromagnetic_Field())
    cs,particle = binding.cross_sections,binding.particle
    proton_transport_preflight(binding,cs,geometry,solvers,sources,field)
    sources.is_build && sources.number_of_particles == 1 || error("Require one built source species.")
    get_dimension(geometry) == 1 && all(iszero,get_boundary_conditions(geometry)) ||
        error("Discrete CSD diagnostic supports only one-dimensional void Cartesian geometry.")
    all(iszero,field.electric_field) && all(iszero,field.magnetic_field) &&
        isnothing(field.spatial_magnetic_field) || error("Discrete CSD diagnostic excludes fields.")
    solver = get_method(solvers,particle)
    kind,is_csd = get_solver_type(solver)
    kind == 5 && is_csd || error("Discrete balance supports pure CSD only.")
    get_quadrature_dimension(solver,1) == 1 ||
        error("Discrete balance requires the one-dimensional angular scalar convention.")
    coupling = get_is_full_coupling(solver)
    schemes,orders,nm = get_schemes(solver,geometry,coupling)
    schemes[1] == "DG" && schemes[4] == "DG" || error("Discrete balance requires native DG space and energy.")
    source = get_source(sources,particle)
    surface = get_surface_sources(source)
    all(value -> value isa Number ? iszero(value) : all(iszero,value),surface) ||
        error("Boundary-source injection is not supported by this diagnostic.")
    index = findfirst(x -> get_tag(x) == get_tag(particle),flux.particles)
    index === nothing && error("No proton flux for discrete balance.")
    fpp = flux.flux_per_particle[index]
    length(fpp.flux) == length(fpp.flux_cutoff) == length(fpp.boundary_flux) == 1 ||
        error("Discrete balance requires exactly one retained native generation.")
    basis = isnothing(sources.explicit_normalization) ? :legacy_source_integral :
        get_source_normalization(sources).basis
    capture = _proton_verified_escape(binding,geometry,solvers,sources,flux,particle,basis,nothing)
    capture.convergence_verified === true || error("Discrete balance requires verified native convergence.")
    ng = length(binding.energy_boundaries_MeV)-1
    nx = get_number_of_voxels(geometry)[1]
    phi = get_flux(flux,particle)
    q = get_volume_sources(source)
    ndims(phi) == 6 && ndims(q) == 6 && size(phi,1) == size(q,1) == ng &&
        size(phi,3) == size(q,3) == nm[5] && size(phi)[4:6] == size(q)[4:6] == (nx,1,1) ||
        error("Native volume moments differ from declared scheme and geometry.")
    # In both native coupling layouts, spatial-average pure energy modes are 1:OE.
    scalar_phi = copy(phi[:,1,1:orders[4],:,1,1])
    scalar_q = copy(q[:,1,1:orders[4],:,1,1])
    material = vec(get_material_per_voxel(geometry))
    eb = binding.energy_boundaries_MeV
    sb = get_boundary_stopping_powers(cs,particle)[:,material]
    midpoint = zeros(ng,nx)
    for g in 1:ng, x in 1:nx
        model = binding.material_data[material[x]].model
        value = ion_transport_coefficients(model,proton_species(),(eb[g]+eb[g+1])/2)
        midpoint[g,x] = value.electronic_stopping_MeV_cm+value.nuclear_stopping_MeV_cm
    end
    cutoff = get_flux_cutoff(flux,particle)
    ndims(cutoff) == 5 && size(cutoff)[3:5] == (nx,1,1) || error("Native cutoff shape differs from geometry.")
    result = _proton_discrete_csd_terms(energy_boundaries_MeV=eb,phi=scalar_phi,q=scalar_q,
        stopping_at_boundaries_MeV_cm=sb,midpoint_response_MeV_cm=midpoint,
        cell_measure_cm_d=vec(geometry.volume_per_voxel),
        outgoing_particles_by_face_group=get_outgoing_current(only(fpp.boundary_flux)),
        outgoing_energy_by_face_group_MeV=capture.raw_escaped_energy_by_face_group_MeV,
        scalar_cutoff_flux=vec(cutoff[1,1,:,1,1]),
        source_normalization_divisor=get_normalization_factor(sources),
        injected_energy_MeV=injected_energy_MeV)
    return merge(result,(schema="radianthts.proton_discrete_csd_balance/v1",
        classification=:signed_native_weak_balance_diagnostic,
        convergence_verified=true,source_basis=basis,
        capture_escape_status=capture.escape_status,
        orders=copy(orders),fully_coupled=coupling,
        stopping_representation=:native_endpoint_linear,
        incoming_upper_energy_particle_flow=0.0,
        binding_input_sha256=binding.input_snapshot_sha256,
        binding_library_sha256=binding.library_snapshot_sha256[]))
end

function _proton_accounting_scalar(value, name; positive=false)
    data = Float64(value)
    isfinite(data) && (positive ? data > 0.0 : data >= 0.0) ||
        error("$(name) must be finite and $(positive ? "positive" : "nonnegative").")
    return data
end

"""
    proton_cutoff_handoff(; scalar_cutoff_flux, last_group_width_MeV, ...)

Score the energy boundary of the native CSD equation, separately from deposition.
`scalar_cutoff_flux` is Radiant's stored `[1,1,...]` moment: it is the last group
width times differential scalar flux at Ecut. Thus particle flow density is
`S(Ecut)*scalar_cutoff_flux/(last_group_width_MeV*source_normalization_divisor)`.
S is TOTAL linear stopping (electronic plus nuclear), in MeV/cm. No 4π factor
is applied: the scalar angular moment already integrates the angular flux.

`cell_measure_cm_d` is length/area/volume in active Cartesian dimension d. No
unobserved transverse extent is invented. Integrated results are particles or
MeV per declared source basis; energy per mass fields have MeV*cm^(3-d)/g
per basis. With explicit 3D source densities these are MeV/g per basis.
Residual kinetic energy belongs to an unresolved below-cutoff handoff, never
automatically to electronic deposition or heat. This does not model full stopping.
"""
function proton_cutoff_handoff(;scalar_cutoff_flux, last_group_width_MeV::Real,
    cutoff_MeV::Real, total_stopping_MeV_cm, cell_measure_cm_d,
    density_g_cm3, geometry_dimension::Integer,
    source_normalization_divisor::Real, source_basis::Symbol,
    physical_normalization::Union{Nothing,Source_Normalization}=nothing)
    geometry_dimension in (1,2,3) || error("Accounting dimension must be 1, 2, or 3.")
    source_basis in (:legacy_source_integral,:per_history,:per_source_particle,:per_second) ||
        error("Declare the transport source basis explicitly.")
    width = _proton_accounting_scalar(last_group_width_MeV,"Last group width";positive=true)
    energy = _proton_accounting_scalar(cutoff_MeV,"Cutoff energy";positive=true)
    1.0 <= energy <= 500.0 || error("Proton cutoff must lie within [1,500] MeV.")
    divisor = _proton_accounting_scalar(source_normalization_divisor,"Source divisor";positive=true)
    phi = _proton_accounting_array(scalar_cutoff_flux,"Scalar cutoff flux")
    stopping = _proton_accounting_array(total_stopping_MeV_cm,"Total stopping")
    measure = _proton_accounting_array(cell_measure_cm_d,"Cell measure";positive=true)
    density = _proton_accounting_array(density_g_cm3,"Density";positive=true)
    size(phi) == size(stopping) == size(measure) == size(density) ||
        error("Cutoff fields, stopping, cell measure and density must have identical shapes.")
    scale = 1.0
    if !isnothing(physical_normalization)
        physical_normalization.basis == source_basis || error("Physical and transport source bases differ.")
        scale = get_physical_scale(physical_normalization)
    end
    denominator = _proton_accounting_scalar(width*divisor,"Cutoff normalization denominator";positive=true)
    multiplier = _proton_accounting_scalar(scale/denominator,"Cutoff normalization multiplier";positive=true)
    particles = stopping .* phi .* multiplier
    handoff = energy .* particles
    per_mass = handoff ./ density
    particle_total = sum(particles .* measure)
    energy_total = sum(handoff .* measure)
    all(isfinite,particles) && all(isfinite,handoff) && all(isfinite,per_mass) &&
        isfinite(particle_total) && isfinite(energy_total) || error("Cutoff score overflowed.")
    basis = isnothing(physical_normalization) ? source_basis : :physical_rate
    return (cutoff_particle_flow_density=particles,
        cutoff_kinetic_handoff_density=handoff,
        cutoff_kinetic_handoff_per_mass=per_mass,
        cutoff_particles=particle_total,
        cutoff_kinetic_handoff_MeV=energy_total,
        cell_measure_cm_d=measure, density_g_cm3=density,
        geometry_dimension=Int(geometry_dimension),source_basis=basis,
        source_normalization_divisor=divisor,
        cutoff_MeV=energy,last_group_width_MeV=width,
        particle_density_units="particles/cm^$(geometry_dimension) per $(basis)",
        energy_density_units="MeV/cm^$(geometry_dimension) per $(basis)",
        energy_per_mass_units="MeV*cm^$(3-geometry_dimension)/g per $(basis)",
        physical_validation=false,full_stopping_qualified=false)
end

"""
Extract separate electronic, recoil and residual-cutoff ownership from a built
single-species native proton solve. The binding's fail-closed preflight is reused;
nonelastic removal and unqualified physical mode cannot bypass it through scoring.
Void escape is attached only with compatible captured basis and finite native
solver/reconstruction/outer receipts. Injection and secondary ownership remain explicit.
"""
function proton_energy_accounting(binding::Proton_Native_Binding,geometry::Geometry,
    solvers::Solvers,sources::Fixed_Sources,flux::Flux;
    field::Electromagnetic_Field=Electromagnetic_Field(),
    physical_normalization::Union{Nothing,Source_Normalization}=nothing)
    cs,particle = binding.cross_sections,binding.particle
    proton_transport_preflight(binding,cs,geometry,solvers,sources,field)
    sources.is_build && sources.number_of_particles == 1 ||
        error("Accounting requires a built source collection with one species.")
    counts = Tuple(get_number_of_voxels(geometry))
    material = reshape(get_material_per_voxel(geometry),counts)
    measure = reshape(geometry.volume_per_voxel,counts)
    cutoff_raw = get_flux_cutoff(flux,particle)
    ndims(cutoff_raw) == 5 && size(cutoff_raw)[3:5] == counts ||
        error("Native cutoff flux shape does not match geometry.")
    cutoff = cutoff_raw[1,1,:,:,:]
    stopping_boundary = get_boundary_stopping_powers(cs,particle)
    stopping = reshape([stopping_boundary[end,m] for m in material],counts)
    density = reshape([get_densities(cs)[m] for m in material],counts)
    basis = isnothing(sources.explicit_normalization) ? :legacy_source_integral :
        get_source_normalization(sources).basis
    boundary = binding.energy_boundaries_MeV
    handoff = proton_cutoff_handoff(scalar_cutoff_flux=cutoff,
        last_group_width_MeV=boundary[end-1]-boundary[end],cutoff_MeV=boundary[end],
        total_stopping_MeV_cm=stopping,cell_measure_cm_d=measure,
        density_g_cm3=density,geometry_dimension=get_dimension(geometry),
        source_normalization_divisor=get_normalization_factor(sources),source_basis=basis,
        physical_normalization=physical_normalization)
    electronic = score_process_responses(cs,geometry,solvers,sources,flux,particle;
        quantity="energy-deposition",mass_normalized=false,
        physical_normalization=physical_normalization).total
    recoil = score_process_responses(cs,geometry,solvers,sources,flux,particle;
        quantity="recoil-handoff",mass_normalized=false,
        physical_normalization=physical_normalization).total
    electronic = _proton_accounting_array(electronic,"Electronic deposition")
    recoil = _proton_accounting_array(recoil,"Recoil handoff")
    electronic_mass = _proton_accounting_array(electronic ./ density,"Electronic mass score")
    recoil_mass = _proton_accounting_array(recoil ./ density,"Recoil mass score")
    electronic_total = _proton_accounting_scalar(sum(electronic .* measure),"Electronic total")
    recoil_total = _proton_accounting_scalar(sum(recoil .* measure),"Recoil total")
    escape = _proton_verified_escape(binding,geometry,solvers,sources,flux,particle,
        basis,physical_normalization)
    return merge(handoff,(electronic_deposition_density=electronic,
        electronic_deposition_per_mass=electronic_mass,
        recoil_handoff_density=recoil,recoil_handoff_per_mass=recoil_mass,
        electronic_deposition_MeV=electronic_total,
        recoil_handoff_MeV=recoil_total,
        binding_input_sha256=binding.input_snapshot_sha256,
        binding_library_sha256=binding.library_snapshot_sha256[]),escape)
end

function _proton_verified_escape(binding,geometry,solvers,sources,flux,particle,basis,physical_normalization)
    absent = (escaped_energy_MeV=missing,escaped_energy_by_face_group_MeV=missing,
        raw_crossing_energy_by_face_group_MeV=missing,raw_escaped_energy_by_face_group_MeV=missing,
        raw_escaped_energy_by_generation_MeV=missing,
        convergence_verified=missing,
        escape_status=:capture_absent,escape_units="MeV per declared source basis")
    index = findfirst(x -> get_tag(x) == get_tag(particle),flux.particles)
    index === nothing && error("No proton flux for energy accounting.")
    fpp = flux.flux_per_particle[index]
    isempty(fpp.boundary_flux) && return absent
    generations = _compatible_boundary_generations(fpp)
    solver = get_method(solvers,particle)
    coupling = get_is_full_coupling(solver)
    _,orders,_ = get_schemes(solver,geometry,coupling)
    _,is_csd = get_solver_type(solver)
    dimension = get_dimension(geometry)
    native_widths = get_voxels_width(geometry)
    widths = [a <= dimension ? native_widths[a] : Float64[] for a in 1:3]
    directions,weights = quadrature(get_quadrature_order(solver),get_quadrature_type(solver),
        dimension,get_quadrature_dimension(solver,dimension))
    directions isa Vector{Float64} && (directions = [directions,0*directions,0*directions])
    for data in generations
        data.energy_boundaries == binding.energy_boundaries_MeV &&
            data.dimension == dimension && data.widths == widths &&
            data.boundary_conditions == get_boundary_conditions(geometry) &&
            data.orders == orders && data.fully_coupled == coupling && data.is_csd == is_csd &&
            data.directions == directions && data.weights == weights ||
            error("Captured energy metadata differs from native binding, geometry or solver.")
    end
    crossing = get_outgoing_energy_current(fpp)
    raw_escape = get_escaped_energy_current(fpp)
    generation_escape = [get_escaped_energy_current(data) for data in generations]
    verified = true
    for data in generations
        length(data.convergence) == length(binding.energy_boundaries_MeV)-1 ||
            error("Captured convergence receipts do not match groups.")
        for (group,row) in enumerate(data.convergence)
            terminal_ok = row.terminal == :solver_tolerance ||
                (row.terminal == :zero_source && row.iterations == 0 &&
                    row.solver_residual == 0.0 && row.reconstruction_residual == 0.0 &&
                    all(all(iszero,view(face,group,:,:,:,:)) for face in data.faces))
            verified &= row.converged === true && terminal_ok &&
                0 <= row.iterations <= row.iteration_cap &&
                isfinite(row.tolerance) && row.tolerance > 0 &&
                isfinite(row.solver_residual) && 0 <= row.solver_residual <= row.tolerance &&
                isfinite(row.reconstruction_residual) && 0 <= row.reconstruction_residual <= row.tolerance
        end
        outer = data.outer_convergence[]
        verified &= outer.converged === true && outer.required isa Bool &&
            0 < outer.iterations <= outer.iteration_cap &&
            isfinite(outer.tolerance) && outer.tolerance > 0 &&
            isfinite(outer.residual) && 0 <= outer.residual <= outer.tolerance
    end
    if !verified
        return merge(absent,(raw_crossing_energy_by_face_group_MeV=crossing,
            raw_escaped_energy_by_face_group_MeV=raw_escape,
            raw_escaped_energy_by_generation_MeV=generation_escape,
            convergence_verified=false,escape_status=:unverified_native_convergence))
    end
    # Diagnose each generation before accumulation: a negative represented void
    # observable cannot be legitimized by a positive later generation.
    if any(any(<(0.0),score) for score in generation_escape)
        return merge(absent,(raw_crossing_energy_by_face_group_MeV=crossing,
            raw_escaped_energy_by_face_group_MeV=raw_escape,
            raw_escaped_energy_by_generation_MeV=generation_escape,
            convergence_verified=true,escape_status=:negative_integrated_void_energy))
    end
    divisor = _proton_accounting_scalar(get_normalization_factor(sources),"Escape source divisor";positive=true)
    scale = 1.0
    if physical_normalization !== nothing
        physical_normalization.basis == basis || error("Escape physical and transport source bases differ.")
        scale = get_physical_scale(physical_normalization)
    end
    multiplier = _proton_accounting_scalar(scale/divisor,"Escape normalization multiplier";positive=true)
    escaped = raw_escape .* multiplier
    all(isfinite,escaped) || error("Normalized escaped energy overflowed.")
    total = _proton_accounting_scalar(sum(escaped),"Integrated escaped energy")
    return (escaped_energy_MeV=total,escaped_energy_by_face_group_MeV=escaped,
        raw_crossing_energy_by_face_group_MeV=crossing,raw_escaped_energy_by_face_group_MeV=raw_escape,
        raw_escaped_energy_by_generation_MeV=generation_escape,
        convergence_verified=true,
        escape_status=:verified_native_capture,
        escape_units="MeV per $(physical_normalization === nothing ? basis : :physical_rate)")
end

"""
Report an incomplete energy ledger honestly. Omitted observations stay `missing`,
including secondaries and leakage; zero must be an explicitly supplied observation.
The optional numeric closure is software arithmetic only, with tolerance chosen
by the caller. It never certifies physical validation or a stopping endpoint.
"""
function proton_energy_accounting_report(score;injected_energy_MeV=missing,
    escaped_energy_MeV=hasproperty(score,:escaped_energy_MeV) ? score.escaped_energy_MeV : missing,
    secondary_transfer_MeV=missing,other_transfer_MeV=missing,
    convergence_verified=hasproperty(score,:convergence_verified) ? score.convergence_verified : missing,
    closure_atol_MeV::Real=0.0,closure_rtol::Real=0.0)
    atol = _proton_accounting_scalar(closure_atol_MeV,"Closure absolute tolerance")
    rtol = _proton_accounting_scalar(closure_rtol,"Closure relative tolerance")
    inputs = Dict("injected_energy_MeV"=>injected_energy_MeV,
        "escaped_energy_MeV"=>escaped_energy_MeV,
        "secondary_transfer_MeV"=>secondary_transfer_MeV,
        "other_transfer_MeV"=>other_transfer_MeV)
    values_checked = Dict{String,Union{Missing,Float64}}()
    for (key,value) in inputs
        values_checked[key] = ismissing(value) ? missing : _proton_accounting_scalar(value,key)
    end
    (ismissing(convergence_verified) || convergence_verified isa Bool) ||
        error("Convergence evidence must be Bool or missing.")
    accounted = _proton_accounting_scalar(sum(_proton_accounting_scalar(getproperty(score,key),string(key)) for
        key in (:electronic_deposition_MeV,:recoil_handoff_MeV,:cutoff_kinetic_handoff_MeV)),
        "Accounted energy sum")
    unknown = sort([key for (key,value) in values_checked if ismissing(value)])
    convergence_verified === true || push!(unknown,"verified_convergence")
    residual = missing
    closure = missing
    if all(!ismissing(value) for value in values(values_checked))
        injected = values_checked["injected_energy_MeV"]
        residual = injected-accounted-values_checked["escaped_energy_MeV"]-
            values_checked["secondary_transfer_MeV"]-values_checked["other_transfer_MeV"]
        isfinite(residual) || error("Derived energy residual is not finite.")
        threshold = _proton_accounting_scalar(atol+rtol*injected,"Derived closure threshold")
        closure = abs(residual) <= threshold
    end
    classification = !isempty(unknown) ? "BLOCKED_INCOMPLETE_LEDGER" :
        closure === true ? "SOFTWARE_ARITHMETIC_CLOSED" : "FAIL_ENERGY_CLOSURE"
    return merge(Dict{String,Any}(values_checked),Dict{String,Any}(
        "schema"=>"radianthts.proton_energy_accounting/v1",
        "classification"=>classification,"unobserved_or_unverified"=>unknown,
        "electronic_deposition_MeV"=>score.electronic_deposition_MeV,
        "recoil_handoff_MeV"=>score.recoil_handoff_MeV,
        "cutoff_kinetic_handoff_MeV"=>score.cutoff_kinetic_handoff_MeV,
        "energy_residual_MeV"=>residual,"arithmetic_closure"=>closure,
        "source_basis"=>string(score.source_basis),
        "physical_validation"=>false,"full_stopping_qualified"=>false))
end
