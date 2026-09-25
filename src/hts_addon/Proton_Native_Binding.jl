"""Provenance-bound material state for a proton stopping/moment table.

The coefficient units are fixed by `Tabulated_Ion_Transport_Model`: MeV/cm, MeV²/cm,
and rad²/cm. Linear interpolation is the only implemented law. A binding never promotes
the table's evidence class; synthetic rows remain software evidence.
"""
_proton_model_snapshot(model::Tabulated_Ion_Transport_Model) = bytes2hex(sha256(repr((
    model.species_id,model.material_id,model.energy_MeV,
    model.electronic_stopping_MeV_cm,model.nuclear_stopping_MeV_cm,
    model.energy_straggling_variance_MeV2_cm,model.angular_variance_rad2_cm,
    model.data_hash,model.qualification_status,
))))

struct Proton_Material_Data
    material::Material
    model::Tabulated_Ion_Transport_Model
    material_state::String
    density_g_cm3::Float64
    interpolation::Symbol
    source_sha256::String
    source_path::Union{Nothing,String}
    uncertainty_status::Symbol
    snapshot_sha256::String

    function Proton_Material_Data(material::Material,model::Tabulated_Ion_Transport_Model;
        material_state::AbstractString,
        density_g_cm3::Real,
        interpolation::Symbol=:linear,
        source_sha256::AbstractString,
        source_path::Union{Nothing,AbstractString}=nothing,
        uncertainty_status::Symbol=:unquantified)
        get_tag(material) == model.material_id || error("Proton table material ID does not match Radiant material.")
        String(material_state) == get_state_of_matter(material) || error("Proton table material state does not match Radiant material.")
        density = Float64(density_g_cm3)
        isfinite(density) && density > 0.0 || error("Proton table density must be positive and finite.")
        isapprox(density,get_density(material);rtol=1.0e-12,atol=0.0) || error("Proton table density does not match Radiant material.")
        interpolation == :linear || error("Only declared linear proton-table interpolation is implemented.")
        hash = String(source_sha256)
        occursin(r"^[0-9a-f]{64}$",hash) || error("Proton table source SHA-256 must be 64 lowercase hex characters.")
        model.data_hash == hash || error("Proton table hash does not match source SHA-256.")
        path = isnothing(source_path) ? nothing : String(source_path)
        if model.qualification_status != :synthetic
            isnothing(path) && error("Candidate/qualified proton tables require a local source artifact.")
            isfile(path) || error("Proton source artifact is missing.")
            open(path,"r") do io
                bytes2hex(sha256(io)) == hash || error("Proton source artifact SHA-256 mismatch.")
            end
        end
        uncertainty_status in (:unquantified,:bounded,:covariance) || error("Unknown proton uncertainty status.")
        return new(material,model,String(material_state),density,interpolation,hash,path,
            uncertainty_status,_proton_model_snapshot(model))
    end
end

"""Explicit nonelastic domain declaration; zero removal is allowed only for a manufactured fixture.

Nonzero removal is represented but cannot enter the current native solver until a coupled
production matrix and correlated residual/event handoff are both bound. This is a fail-closed
extension point, not a zero-filled physical kernel.
"""
_proton_nonelastic_snapshot(data) = bytes2hex(sha256(repr((
    data.material_id,data.energy_MeV,data.removal_cm_inv,data.event_family_ids,
    data.source_sha256,data.qualification_status,
))))

struct Proton_Nonelastic_Data
    material_id::String
    energy_MeV::Vector{Float64}
    removal_cm_inv::Vector{Float64}
    event_family_ids::Vector{String}
    source_sha256::String
    source_path::Union{Nothing,String}
    qualification_status::Symbol
    snapshot_sha256::String

    function Proton_Nonelastic_Data(;material_id::AbstractString,energy_MeV,
        removal_cm_inv,event_family_ids=String[],source_sha256::AbstractString,
        source_path::Union{Nothing,AbstractString}=nothing,
        qualification_status::Symbol=:synthetic)
        isempty(material_id) && error("Nonelastic material ID is required.")
        energy = Float64.(energy_MeV)
        length(energy) >= 2 && all(diff(energy) .> 0.0) &&
            energy[1] >= 0.0 && energy[end] <= RADIANT_MAX_ION_KINETIC_ENERGY_MEV ||
            error("Nonelastic energy grid must increase within [0,500] MeV.")
        removal = Float64.(removal_cm_inv)
        length(removal) == length(energy) &&
            all(x -> isfinite(x) && x >= 0.0,removal) ||
            error("Nonelastic removal must be finite, nonnegative, and match the grid.")
        ids = String.(event_family_ids)
        length(unique(ids)) == length(ids) && all(id -> !isempty(id),ids) ||
            error("Nonelastic event-family IDs must be unique and nonempty.")
        hash = String(source_sha256)
        occursin(r"^[0-9a-f]{64}$",hash) || error("Nonelastic source SHA-256 must be 64 lowercase hex characters.")
        qualification_status in (:synthetic,:candidate,:qualified) || error("Unknown nonelastic status.")
        path = isnothing(source_path) ? nothing : String(source_path)
        if qualification_status != :synthetic
            isnothing(path) && error("Candidate/qualified nonelastic data require a local source artifact.")
            isfile(path) || error("Nonelastic source artifact is missing.")
            open(path,"r") do io
                bytes2hex(sha256(io)) == hash || error("Nonelastic source artifact SHA-256 mismatch.")
            end
        end
        if all(iszero,removal)
            qualification_status == :synthetic || error("Zero nonelastic removal is allowed only as a manufactured fixture.")
            isempty(ids) || error("Zero nonelastic removal cannot declare event families.")
        else
            isempty(ids) && error("Nonzero nonelastic removal requires correlated event-family IDs.")
        end
        snapshot = _proton_nonelastic_snapshot((;
            material_id=String(material_id),energy_MeV=energy,removal_cm_inv=removal,
            event_family_ids=ids,source_sha256=hash,qualification_status=qualification_status,
        ))
        return new(String(material_id),energy,removal,ids,hash,path,qualification_status,snapshot)
    end
end

struct Proton_Native_Binding
    cross_sections::Cross_Sections
    particle::Particle
    material_data::Vector{Proton_Material_Data}
    nonelastic_data::Vector{Proton_Nonelastic_Data}
    energy_boundaries_MeV::Vector{Float64}
    qualification_mode::Symbol
end

function _proton_nonelastic_value(data::Proton_Nonelastic_Data,energy_MeV::Float64)
    data.energy_MeV[1] <= energy_MeV <= data.energy_MeV[end] || error(
        "Nonelastic energy lies outside its declared table; extrapolation is forbidden.",
    )
    index = searchsortedlast(data.energy_MeV,energy_MeV)
    index == length(data.energy_MeV) && return data.removal_cm_inv[end]
    low = data.energy_MeV[index]
    weight = (energy_MeV-low)/(data.energy_MeV[index+1]-low)
    return (1-weight)*data.removal_cm_inv[index]+weight*data.removal_cm_inv[index+1]
end

function _proton_native_group_library(binding::Proton_Native_Binding)
    cs = binding.cross_sections
    boundaries = binding.energy_boundaries_MeV
    ng = length(boundaries)-1
    library = Array{Multigroup_Cross_Sections}(undef,1,length(binding.material_data))
    for (index,data) in enumerate(binding.material_data)
        nonelastic = binding.nonelastic_data[index]
        mcs = Multigroup_Cross_Sections(Int64(ng))
        centers = [(boundaries[g]+boundaries[g+1])/2 for g in 1:ng]
        coefficients = [ion_transport_coefficients(data.model,proton_species(),E) for E in centers]
        boundary_coefficients = [ion_transport_coefficients(data.model,proton_species(),E) for E in boundaries]
        electronic = [c.electronic_stopping_MeV_cm for c in coefficients]
        nuclear = [c.nuclear_stopping_MeV_cm for c in coefficients]
        angular = [c.angular_variance_rad2_cm/2 for c in coefficients]
        total_stopping = electronic .+ nuclear
        boundary_stopping = [c.electronic_stopping_MeV_cm+c.nuclear_stopping_MeV_cm
            for c in boundary_coefficients]
        removal = [_proton_nonelastic_value(nonelastic,E) for E in centers]
        set_total(mcs,removal)
        set_absorption(mcs,copy(removal))
        set_scattering(mcs,zeros(Float64,ng,ng,1))
        set_boundary_stopping_powers(mcs,boundary_stopping)
        set_stopping_powers(mcs,total_stopping)
        set_momentum_transfer(mcs,angular)
        set_energy_deposition(mcs,vcat(electronic,0.0))
        set_charge_deposition(mcs,zeros(Float64,ng+1))
        key = "/in=$(get_tag(binding.particle))/out=local"
        add_response_channel!(mcs,"energy-deposition|Proton/electronic-stopping"*key,vcat(electronic,0.0))
        add_response_channel!(mcs,"stopping-power|Proton/electronic-stopping"*key,electronic)
        add_response_channel!(mcs,"stopping-power|Proton/nuclear-stopping"*key,nuclear)
        add_response_channel!(mcs,"recoil-handoff|Proton/nuclear-stopping"*key,nuclear)
        add_response_channel!(mcs,"momentum-transfer|Proton/angular-variance"*key,angular)
        add_response_channel!(mcs,"absorption|Proton/nonelastic"*key,removal)
        add_response_channel!(mcs,"total|Proton/nonelastic"*key,removal)
        library[1,index] = mcs
    end
    set_number_of_groups(cs,[ng])
    set_energy_boundaries(cs,[boundaries])
    set_energy(cs,[(boundaries[1]+boundaries[2])/2])
    set_cutoff(cs,[boundaries[end]])
    set_multigroup_cross_sections(cs,library)
    return cs
end

"""Bind proton tables to the normal `Cross_Sections.build` and native `transport` path."""
function bind_proton_native(particle::Particle,material_data::Vector{Proton_Material_Data},
    nonelastic_data::Vector{Proton_Nonelastic_Data},energy_boundaries_MeV;
    qualification_mode::Symbol=:software)
    get_type(particle) == Proton || error("Native proton binding requires Proton().")
    !isempty(material_data) || error("At least one proton material table is required.")
    length(nonelastic_data) == length(material_data) || error("Every proton material needs an explicit nonelastic declaration.")
    boundaries = Float64.(energy_boundaries_MeV)
    length(boundaries) >= 2 && all(diff(boundaries) .< 0.0) &&
        boundaries[1] <= RADIANT_MAX_ION_KINETIC_ENERGY_MEV && boundaries[end] >= 1.0 ||
        error("Proton group boundaries must descend within [1,500] MeV.")
    qualification_mode in (:software,:physical) || error("Unknown proton qualification mode.")
    tags = get_tag.([data.material for data in material_data])
    length(unique(tags)) == length(tags) || error("Proton material IDs must be unique.")
    for (data,nonelastic) in zip(material_data,nonelastic_data)
        data.model.species_id == "proton" || error("Bound stopping table is not for protons.")
        nonelastic.material_id == get_tag(data.material) || error("Nonelastic material ID differs from stopping table.")
        data.model.energy_MeV[1] <= boundaries[end] <= boundaries[1] <= data.model.energy_MeV[end] ||
            error("Proton group boundaries exceed the stopping-table domain; extrapolation is forbidden.")
        nonelastic.energy_MeV[1] <= boundaries[end] <= boundaries[1] <= nonelastic.energy_MeV[end] ||
            error("Proton group boundaries exceed the nonelastic domain; extrapolation is forbidden.")
    end
    cs = Cross_Sections()
    set_source(cs,"provided")
    set_particles(cs,particle)
    set_materials(cs,[data.material for data in material_data])
    binding = Proton_Native_Binding(cs,particle,copy(material_data),copy(nonelastic_data),
        boundaries,qualification_mode)
    set_provided_builder(cs,_ -> _proton_native_group_library(binding))
    set_transport_preflight(cs,(xs,geo,solvers,sources,field) ->
        proton_transport_preflight(binding,xs,geo,solvers,sources,field))
    return binding
end

"""Check blockers before any native solve; direct `Computation_Unit.run` also invokes this guard."""
function proton_transport_preflight(binding::Proton_Native_Binding,cs::Cross_Sections,
    geometry::Geometry,solvers::Solvers,sources::Fixed_Sources,
    field::Electromagnetic_Field)
    cs === binding.cross_sections || error("Proton binding was detached from its cross sections.")
    get_energy_boundaries(cs,binding.particle) == binding.energy_boundaries_MeV ||
        error("Proton group boundaries changed after binding.")
    for (data,nonelastic) in zip(binding.material_data,binding.nonelastic_data)
        _proton_model_snapshot(data.model) == data.snapshot_sha256 ||
            error("Proton stopping table changed after binding.")
        _proton_nonelastic_snapshot(nonelastic) == nonelastic.snapshot_sha256 ||
            error("Proton nonelastic table changed after binding.")
        get_state_of_matter(data.material) == data.material_state &&
            isapprox(get_density(data.material),data.density_g_cm3;rtol=1.0e-12,atol=0.0) ||
            error("Proton material state or density changed after binding.")
        for (path,hash) in ((data.source_path,data.source_sha256),
            (nonelastic.source_path,nonelastic.source_sha256))
            if !isnothing(path)
                isfile(path) || error("Proton bound source artifact disappeared.")
                open(path,"r") do io
                    bytes2hex(sha256(io)) == hash ||
                        error("Proton bound source artifact changed after binding.")
                end
            end
        end
    end
    get_type(geometry) == "cartesian" || error("Mapped/curved proton transport is not yet native.")
    get_number_of_particles(solvers) == 1 || error("This proton binding supports one transported species.")
    solver = get_method(solvers,binding.particle)
    solver isa SN || error("This proton binding requires the native SN solver.")
    solver_type,is_csd = get_solver_type(solver)
    is_csd || error("Proton stopping requires a native CSD-capable solver.")
    angular = any(any(>(0.0),data.model.angular_variance_rad2_cm) for data in binding.material_data)
    angular && !(solver_type in (2,4,6)) && error("Nonzero angular variance requires a BFP/FP solver.")
    any(any(>(0.0),data.removal_cm_inv) for data in binding.nonelastic_data) && error(
        "Nonzero proton nonelastic removal has no bound production matrix and correlated event handoff.",
    )
    if binding.qualification_mode == :physical
        all(data.model.qualification_status == :qualified &&
            data.uncertainty_status != :unquantified for data in binding.material_data) ||
            error("Physical proton transport requires qualified stopping tables and uncertainties.")
        all(data.qualification_status == :qualified for data in binding.nonelastic_data) ||
            error("Physical proton transport requires qualified nonelastic event kernels.")
    end
    get_tag(binding.particle) in get_tag.(get_particles(sources)) || error("Proton fixed source is absent.")
    return true
end

function proton_binding_receipt(binding::Proton_Native_Binding)
    classification = if binding.qualification_mode == :physical
        "BLOCKED_INPUT"
    elseif all(data.model.qualification_status == :synthetic for data in binding.material_data) &&
        all(data.qualification_status == :synthetic for data in binding.nonelastic_data)
        "SOFTWARE_VERIFIED_SYNTHETIC"
    else
        "PHYSICAL_CANDIDATE"
    end
    return Dict(
        "schema" => "radianthts.proton_native_binding/v1",
        "classification" => classification,
        "particle" => get_tag(binding.particle),
        "energy_boundaries_MeV" => copy(binding.energy_boundaries_MeV),
        "materials" => [Dict(
            "id" => get_tag(data.material),"state" => data.material_state,
            "density_g_cm3" => data.density_g_cm3,
            "stopping_sha256" => data.source_sha256,
            "stopping_status" => string(data.model.qualification_status),
            "nonelastic_sha256" => binding.nonelastic_data[i].source_sha256,
            "nonelastic_status" => string(binding.nonelastic_data[i].qualification_status),
        ) for (i,data) in enumerate(binding.material_data)],
        "nonelastic_connected" => false,
    )
end
