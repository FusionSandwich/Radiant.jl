"""Finite, nonnegative scalar observables are required; negative moments are not clipped."""
function _proton_accounting_array(value, name; positive=false)
    data = Float64.(value)
    isempty(data) && error("$(name) cannot be empty.")
    all(x -> isfinite(x) && (positive ? x > 0.0 : x >= 0.0),data) ||
        error("$(name) must be finite and $(positive ? "positive" : "nonnegative").")
    return data
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
This partial ledger has no inferred leakage, secondaries or convergence evidence.
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
    return merge(handoff,(electronic_deposition_density=electronic,
        electronic_deposition_per_mass=electronic_mass,
        recoil_handoff_density=recoil,recoil_handoff_per_mass=recoil_mass,
        electronic_deposition_MeV=electronic_total,
        recoil_handoff_MeV=recoil_total,
        binding_input_sha256=binding.input_snapshot_sha256,
        binding_library_sha256=binding.library_snapshot_sha256[]))
end

"""
Report an incomplete energy ledger honestly. Omitted observations stay `missing`,
including secondaries and leakage; zero must be an explicitly supplied observation.
The optional numeric closure is software arithmetic only, with tolerance chosen
by the caller. It never certifies physical validation or a stopping endpoint.
"""
function proton_energy_accounting_report(score;injected_energy_MeV=missing,
    escaped_energy_MeV=missing,secondary_transfer_MeV=missing,
    other_transfer_MeV=missing,convergence_verified=missing,
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
