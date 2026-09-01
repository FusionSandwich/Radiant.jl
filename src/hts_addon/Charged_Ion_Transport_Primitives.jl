const RADIANT_MAX_ION_KINETIC_ENERGY_MEV = 500.0
const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0

"""Species-generic charged ion definition; mass is rest energy in MeV/c^2."""
struct Ion_Species
    species_id::String
    atomic_number::Int64
    mass_number::Int64
    rest_mass_MeV::Float64
    charge_e::Float64
    provenance_hash::String
end

function Ion_Species(;species_id::AbstractString,atomic_number::Integer,mass_number::Integer,
    rest_mass_MeV::Real,charge_e::Real,provenance_hash::AbstractString)
    isempty(species_id) && error("Ion species ID cannot be empty.")
    atomic_number >= 1 && mass_number >= atomic_number || error("Ion Z/A values are invalid.")
    mass = Float64(rest_mass_MeV)
    charge = Float64(charge_e)
    isfinite(mass) && mass > 0.0 || error("Ion rest mass must be positive.")
    isfinite(charge) && charge != 0.0 || error("Ion charge state must be finite and nonzero.")
    abs(charge) <= atomic_number || error("Ion charge magnitude cannot exceed atomic number.")
    isempty(provenance_hash) && error("Ion mass/charge provenance hash cannot be empty.")
    return Ion_Species(String(species_id),Int64(atomic_number),Int64(mass_number),mass,
        charge,String(provenance_hash))
end

proton_species(;provenance_hash::AbstractString="codata-proton-rest-mass") = Ion_Species(
    species_id="proton",atomic_number=1,mass_number=1,rest_mass_MeV=938.2720894,
    charge_e=1.0,provenance_hash=provenance_hash,
)

function particle_from_ion(species::Ion_Species)
    return Particle(species.species_id,species.rest_mass_MeV,species.charge_e)
end

struct Relativistic_Ion_Kinematics
    kinetic_energy_MeV::Float64
    total_energy_MeV::Float64
    momentum_MeV_c::Float64
    gamma::Float64
    beta::Float64
    speed_m_s::Float64
end

function relativistic_ion_kinematics(species::Ion_Species,kinetic_energy_MeV::Real)
    kinetic = Float64(kinetic_energy_MeV)
    isfinite(kinetic) && 0.0 <= kinetic <= RADIANT_MAX_ION_KINETIC_ENERGY_MEV || error(
        "Ion kinetic energy must lie in [0,$(RADIANT_MAX_ION_KINETIC_ENERGY_MEV)] MeV.",
    )
    total = species.rest_mass_MeV+kinetic
    momentum = sqrt(max(0.0,total^2-species.rest_mass_MeV^2))
    gamma = total/species.rest_mass_MeV
    beta = momentum == 0.0 ? 0.0 : momentum/total
    return Relativistic_Ion_Kinematics(
        kinetic,total,momentum,gamma,beta,beta*SPEED_OF_LIGHT_M_PER_S,
    )
end

"""One weighted ion source state; synthetic and physical-candidate status are explicit."""
struct Charged_Ion_Source_State
    source_id::String
    species::Ion_Species
    kinetic_energy_MeV::Float64
    position_cm::NTuple{3,Float64}
    direction::NTuple{3,Float64}
    statistical_weight::Float64
    time_s::Float64
    source_hash::String
    qualification_status::Symbol
end

function Charged_Ion_Source_State(;source_id::AbstractString,species::Ion_Species,
    kinetic_energy_MeV::Real,position_cm,direction,statistical_weight::Real=1.0,
    time_s::Real=0.0,source_hash::AbstractString,qualification_status::Symbol=:synthetic)
    isempty(source_id) && error("Ion source ID cannot be empty.")
    relativistic_ion_kinematics(species,kinetic_energy_MeV)
    position = Float64.(position_cm)
    length(position) == 3 && all(isfinite,position) || error("Ion source position is invalid.")
    unit_direction = _atlas_normalize(direction,"ion source direction")
    weight = Float64(statistical_weight)
    time = Float64(time_s)
    isfinite(weight) && weight > 0.0 || error("Ion source weight must be positive.")
    isfinite(time) && time >= 0.0 || error("Ion source time must be nonnegative.")
    isempty(source_hash) && error("Ion source hash cannot be empty.")
    qualification_status in (:synthetic,:candidate,:qualified) || error(
        "Ion source qualification status is invalid.",
    )
    return Charged_Ion_Source_State(
        String(source_id),species,Float64(kinetic_energy_MeV),
        (position[1],position[2],position[3]),
        (unit_direction[1],unit_direction[2],unit_direction[3]),weight,time,
        String(source_hash),qualification_status,
    )
end

"""Hash-bound tabulation of electronic/nuclear stopping and transport scattering moments."""
struct Tabulated_Ion_Transport_Model
    species_id::String
    material_id::String
    energy_MeV::Vector{Float64}
    electronic_stopping_MeV_cm::Vector{Float64}
    nuclear_stopping_MeV_cm::Vector{Float64}
    energy_straggling_variance_MeV2_cm::Vector{Float64}
    angular_variance_rad2_cm::Vector{Float64}
    data_hash::String
    qualification_status::Symbol
end

function Tabulated_Ion_Transport_Model(;species_id::AbstractString,material_id::AbstractString,
    energy_MeV,electronic_stopping_MeV_cm,nuclear_stopping_MeV_cm,
    energy_straggling_variance_MeV2_cm,angular_variance_rad2_cm,
    data_hash::AbstractString,qualification_status::Symbol=:synthetic)
    isempty(species_id) && error("Ion transport species ID cannot be empty.")
    isempty(material_id) && error("Ion transport material ID cannot be empty.")
    energy = Float64.(energy_MeV)
    length(energy) >= 2 && all(diff(energy) .> 0.0) || error(
        "Ion transport energy grid must be increasing.",
    )
    energy[1] >= 0.0 && energy[end] <= RADIANT_MAX_ION_KINETIC_ENERGY_MEV || error(
        "Ion transport table lies outside the supported energy interval.",
    )
    arrays = (
        Float64.(electronic_stopping_MeV_cm),Float64.(nuclear_stopping_MeV_cm),
        Float64.(energy_straggling_variance_MeV2_cm),Float64.(angular_variance_rad2_cm),
    )
    all(values -> length(values) == length(energy) &&
        all(value -> isfinite(value) && value >= 0.0,values),arrays) || error(
        "Ion transport tables must be finite, nonnegative, and share the energy grid.",
    )
    isempty(data_hash) && error("Ion transport data hash cannot be empty.")
    qualification_status in (:synthetic,:candidate,:qualified) || error(
        "Ion transport model status is invalid.",
    )
    return Tabulated_Ion_Transport_Model(
        String(species_id),String(material_id),energy,arrays...,String(data_hash),
        qualification_status,
    )
end

function _ion_table_value(model::Tabulated_Ion_Transport_Model,values,energy_MeV)
    energy = Float64(energy_MeV)
    model.energy_MeV[1] <= energy <= model.energy_MeV[end] || error(
        "Ion energy lies outside the tabulated model; extrapolation is forbidden.",
    )
    index = searchsortedlast(model.energy_MeV,energy)
    index == length(model.energy_MeV) && return values[end]
    lower = model.energy_MeV[index]
    fraction = (energy-lower)/(model.energy_MeV[index+1]-lower)
    return (1.0-fraction)*values[index]+fraction*values[index+1]
end

function ion_transport_coefficients(model::Tabulated_Ion_Transport_Model,
    species::Ion_Species,energy_MeV::Real)
    model.species_id == species.species_id || error("Ion transport species does not match source.")
    relativistic_ion_kinematics(species,energy_MeV)
    return (
        electronic_stopping_MeV_cm=_ion_table_value(model,model.electronic_stopping_MeV_cm,energy_MeV),
        nuclear_stopping_MeV_cm=_ion_table_value(model,model.nuclear_stopping_MeV_cm,energy_MeV),
        energy_straggling_variance_MeV2_cm=_ion_table_value(
            model,model.energy_straggling_variance_MeV2_cm,energy_MeV),
        angular_variance_rad2_cm=_ion_table_value(model,model.angular_variance_rad2_cm,energy_MeV),
    )
end

struct Nonelastic_Secondary_Route
    channel_id::String
    production_owner::String
    handoff_schema::String
    differential_data_hash::String
    status::Symbol
end

function Nonelastic_Secondary_Route(;channel_id::AbstractString,production_owner::AbstractString,
    handoff_schema::AbstractString,differential_data_hash::AbstractString,status::Symbol)
    all(value -> !isempty(value),(channel_id,production_owner,handoff_schema,differential_data_hash)) || error(
        "Nonelastic routing requires complete ownership and data lineage.",
    )
    status in (:candidate,:qualified) || error(
        "Nonelastic routing requires candidate or qualified differential data.",
    )
    return Nonelastic_Secondary_Route(String(channel_id),String(production_owner),
        String(handoff_schema),String(differential_data_hash),status)
end

struct Ion_Transport_Step_Result
    initial_energy_MeV::Float64
    final_energy_MeV::Float64
    electronic_loss_MeV::Float64
    nuclear_loss_MeV::Float64
    energy_variance_MeV2::Float64
    angular_variance_rad2::Float64
    initial_direction::NTuple{3,Float64}
    final_direction::NTuple{3,Float64}
    path_length_cm::Float64
    magnetic_rotation_rad::Float64
    nonelastic_route::Union{Nothing,Nonelastic_Secondary_Route}
    model_hash::String
end

function magnetic_direction_step(species::Ion_Species,kinetic_energy_MeV::Real,
    direction,magnetic_field_T,path_length_cm::Real)
    unit_direction = _atlas_normalize(collect(direction),"ion transport direction")
    field = Float64.(magnetic_field_T)
    length(field) == 3 && all(isfinite,field) || error("Magnetic field must be finite.")
    length_cm = Float64(path_length_cm)
    isfinite(length_cm) && length_cm >= 0.0 || error("Ion path length must be nonnegative.")
    field_strength = norm(field)
    kinematics = relativistic_ion_kinematics(species,kinetic_energy_MeV)
    if field_strength == 0.0 || length_cm == 0.0 || kinematics.momentum_MeV_c == 0.0
        return (unit_direction[1],unit_direction[2],unit_direction[3]),0.0
    end
    momentum_GeV_c = kinematics.momentum_MeV_c/1000.0
    radius_m = momentum_GeV_c/(0.299792458*abs(species.charge_e)*field_strength)
    angle = sign(species.charge_e)*(length_cm/100.0)/radius_m
    axis = field/field_strength
    rotated = _atlas_rodrigues(unit_direction,axis,angle)
    rotated = _atlas_normalize(rotated,"magnetically rotated ion direction")
    return (rotated[1],rotated[2],rotated[3]),angle
end

"""
    ion_transport_step(source, model, path_length_cm; magnetic_field_T, nonelastic_route)

Advance deterministic energy-loss and angular-variance moments. No random angular law or
nonelastic yield is invented. A step that would exhaust the particle energy is rejected instead
of clipping the result to zero.
"""
function ion_transport_step(source::Charged_Ion_Source_State,
    model::Tabulated_Ion_Transport_Model,path_length_cm::Real;
    magnetic_field_T=[0.0,0.0,0.0],
    nonelastic_route::Union{Nothing,Nonelastic_Secondary_Route}=nothing)
    length_value = Float64(path_length_cm)
    isfinite(length_value) && length_value >= 0.0 || error("Ion step length must be nonnegative.")
    coefficients = ion_transport_coefficients(model,source.species,source.kinetic_energy_MeV)
    electronic = coefficients.electronic_stopping_MeV_cm*length_value
    nuclear = coefficients.nuclear_stopping_MeV_cm*length_value
    loss = electronic+nuclear
    loss <= source.kinetic_energy_MeV || error(
        "Ion step would exhaust kinetic energy; reduce the step instead of clipping.",
    )
    direction,angle = magnetic_direction_step(
        source.species,source.kinetic_energy_MeV,source.direction,magnetic_field_T,length_value,
    )
    return Ion_Transport_Step_Result(
        source.kinetic_energy_MeV,source.kinetic_energy_MeV-loss,electronic,nuclear,
        coefficients.energy_straggling_variance_MeV2_cm*length_value,
        coefficients.angular_variance_rad2_cm*length_value,source.direction,direction,
        length_value,angle,nonelastic_route,model.data_hash,
    )
end

function synthetic_proton_transport_fixture()
    species = proton_species(provenance_hash="synthetic-codata-binding")
    source = Charged_Ion_Source_State(
        source_id="synthetic-proton",species=species,kinetic_energy_MeV=500.0,
        position_cm=[0.0,0.0,0.0],direction=[1.0,0.0,0.0],source_hash="synthetic-source",
        qualification_status=:synthetic,
    )
    model = Tabulated_Ion_Transport_Model(
        species_id="proton",material_id="synthetic-material",
        energy_MeV=[0.0,100.0,500.0],electronic_stopping_MeV_cm=[1.0,1.0,1.0],
        nuclear_stopping_MeV_cm=[0.1,0.1,0.1],
        energy_straggling_variance_MeV2_cm=[0.01,0.01,0.01],
        angular_variance_rad2_cm=[1.0e-4,1.0e-4,1.0e-4],
        data_hash="synthetic-transport-table",qualification_status=:synthetic,
    )
    return species,source,model
end
