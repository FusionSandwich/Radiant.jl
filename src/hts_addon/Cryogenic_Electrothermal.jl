"""Positive temperature-dependent scalar property with linear interpolation."""
struct Tabulated_Cryogenic_Property
    temperature_K::Vector{Float64}
    values::Vector{Float64}
    units::String
    property_hash::String
    qualification_status::Symbol

    function Tabulated_Cryogenic_Property(
        temperature_K::AbstractVector{<:Real},
        values::AbstractVector{<:Real};
        units::AbstractString,
        property_hash::AbstractString,
        qualification_status::Symbol=:candidate,
        allow_zero::Bool=false,
    )
        temperatures = Float64.(temperature_K)
        data = Float64.(values)
        length(temperatures) ≥ 1 && length(data) == length(temperatures) || error(
            "Cryogenic property temperature and value arrays are inconsistent.",
        )
        all(isfinite,temperatures) && all(temperatures .> 0.0) || error(
            "Cryogenic property temperatures must be finite and positive.",
        )
        length(temperatures) == 1 || all(diff(temperatures) .> 0.0) || error(
            "Cryogenic property temperatures must be strictly increasing.",
        )
        predicate = allow_zero ? (value -> isfinite(value) && value ≥ 0.0) :
                                 (value -> isfinite(value) && value > 0.0)
        all(predicate,data) || error("Cryogenic property values are outside their allowed range.")
        isempty(units) && error("Cryogenic property units cannot be empty.")
        isempty(property_hash) && error("Cryogenic property hash cannot be empty.")
        qualification_status in (:verification,:candidate,:qualified) || error(
            "Cryogenic property qualification status is invalid.",
        )
        return new(
            temperatures,data,String(units),String(property_hash),qualification_status,
        )
    end
end

function constant_cryogenic_property(
    value::Real;
    units::AbstractString,
    property_hash::AbstractString="constant-fixture",
    qualification_status::Symbol=:verification,
    allow_zero::Bool=false,
)
    return Tabulated_Cryogenic_Property(
        [1.0],[Float64(value)];
        units=units,property_hash=property_hash,
        qualification_status=qualification_status,allow_zero=allow_zero,
    )
end

function property_value(this::Tabulated_Cryogenic_Property,temperature_K::Real)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error(
        "Cryogenic property query temperature must be finite and positive.",
    )
    length(this.temperature_K) == 1 && return this.values[1]
    temperature < this.temperature_K[1] && error(
        "Cryogenic property query lies below the qualified temperature range.",
    )
    temperature > this.temperature_K[end] && error(
        "Cryogenic property query lies above the qualified temperature range.",
    )
    temperature == this.temperature_K[end] && return this.values[end]
    upper = searchsortedfirst(this.temperature_K,temperature)
    upper == 1 && return this.values[1]
    this.temperature_K[upper] == temperature && return this.values[upper]
    lower = upper-1
    fraction = (temperature-this.temperature_K[lower])/
               (this.temperature_K[upper]-this.temperature_K[lower])
    return (1.0-fraction)*this.values[lower]+fraction*this.values[upper]
end

struct Cryogenic_Thermal_Node
    node_id::String
    material_tag::String
    heat_capacity_J_K::Tabulated_Cryogenic_Property
    bath_conductance_W_K::Tabulated_Cryogenic_Property
    initial_temperature_K::Float64
    metadata::Dict{String,String}

    function Cryogenic_Thermal_Node(
        node_id::AbstractString,
        material_tag::AbstractString,
        heat_capacity_J_K::Tabulated_Cryogenic_Property,
        bath_conductance_W_K::Tabulated_Cryogenic_Property;
        initial_temperature_K::Real,
        metadata::AbstractDict=Dict{String,String}(),
    )
        isempty(node_id) && error("Thermal node identifier cannot be empty.")
        isempty(material_tag) && error("Thermal node material tag cannot be empty.")
        temperature = Float64(initial_temperature_K)
        isfinite(temperature) && temperature > 0.0 || error(
            "Thermal node initial temperature must be finite and positive.",
        )
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(node_id),String(material_tag),heat_capacity_J_K,bath_conductance_W_K,
            temperature,metadata_string,
        )
    end
end

struct Cryogenic_Thermal_Link
    node_a::String
    node_b::String
    conductance_W_K::Tabulated_Cryogenic_Property
    interface_hash::String

    function Cryogenic_Thermal_Link(
        node_a::AbstractString,
        node_b::AbstractString,
        conductance_W_K::Tabulated_Cryogenic_Property;
        interface_hash::AbstractString,
    )
        (isempty(node_a) || isempty(node_b)) && error(
            "Thermal link node identifiers cannot be empty.",
        )
        node_a == node_b && error("Thermal link must connect two different nodes.")
        isempty(interface_hash) && error("Thermal interface hash cannot be empty.")
        return new(String(node_a),String(node_b),conductance_W_K,String(interface_hash))
    end
end

"""Smooth superconducting-to-normal resistance model with field and damage Tc shifts."""
struct HTS_Transition_Resistance_Model
    normal_resistance_ohm::Float64
    residual_resistance_ohm::Float64
    critical_temperature_K::Float64
    transition_width_K::Float64
    field_tc_suppression_K_per_T::Float64
    damage_tc_shift_K::Float64
    model_hash::String
    qualification_status::Symbol

    function HTS_Transition_Resistance_Model(;
        normal_resistance_ohm::Real,
        residual_resistance_ohm::Real=0.0,
        critical_temperature_K::Real,
        transition_width_K::Real,
        field_tc_suppression_K_per_T::Real=0.0,
        damage_tc_shift_K::Real=0.0,
        model_hash::AbstractString,
        qualification_status::Symbol=:candidate,
    )
        values = Float64[
            normal_resistance_ohm,residual_resistance_ohm,critical_temperature_K,
            transition_width_K,field_tc_suppression_K_per_T,damage_tc_shift_K,
        ]
        all(isfinite,values) || error("HTS resistance-model parameters must be finite.")
        values[1] > 0.0 || error("Normal resistance must be positive.")
        0.0 ≤ values[2] ≤ values[1] || error(
            "Residual resistance must lie between zero and normal resistance.",
        )
        values[3] > 0.0 && values[4] > 0.0 || error(
            "Critical temperature and transition width must be positive.",
        )
        values[5] ≥ 0.0 || error("Field Tc suppression coefficient must be nonnegative.")
        isempty(model_hash) && error("HTS resistance-model hash cannot be empty.")
        qualification_status in (:verification,:candidate,:qualified) || error(
            "HTS resistance-model qualification status is invalid.",
        )
        return new(values...,String(model_hash),qualification_status)
    end
end

function effective_critical_temperature(
    this::HTS_Transition_Resistance_Model,
    magnetic_field_T::Real,
)
    field = Float64(magnetic_field_T)
    isfinite(field) && field ≥ 0.0 || error("Magnetic-field magnitude must be nonnegative.")
    return this.critical_temperature_K-
           this.field_tc_suppression_K_per_T*field-
           this.damage_tc_shift_K
end

function hts_resistance(
    this::HTS_Transition_Resistance_Model,
    temperature_K::Real,
    magnetic_field_T::Real=0.0,
)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0.0 || error("HTS temperature must be positive.")
    transition_temperature = effective_critical_temperature(this,magnetic_field_T)
    logistic = 1.0/(1.0+exp(-(temperature-transition_temperature)/this.transition_width_K))
    return this.residual_resistance_ohm+
           (this.normal_resistance_ohm-this.residual_resistance_ohm)*logistic
end

struct Electrothermal_Circuit
    mode::Symbol
    bias_current_A::Float64
    bias_voltage_V::Float64
    series_resistance_ohm::Float64
    inductance_H::Float64

    function Electrothermal_Circuit(;
        mode::Symbol,
        bias_current_A::Real=0.0,
        bias_voltage_V::Real=0.0,
        series_resistance_ohm::Real=0.0,
        inductance_H::Real=0.0,
    )
        mode in (:current_bias,:voltage_bias) || error(
            "Electrothermal circuit mode must be current_bias or voltage_bias.",
        )
        values = Float64[
            bias_current_A,bias_voltage_V,series_resistance_ohm,inductance_H,
        ]
        all(isfinite,values) && all(values .≥ 0.0) || error(
            "Electrothermal circuit values must be finite and nonnegative.",
        )
        mode == :current_bias && values[1] == 0.0 && error(
            "Current-biased circuit requires a positive bias current.",
        )
        mode == :voltage_bias && values[2] == 0.0 && error(
            "Voltage-biased circuit requires a positive bias voltage.",
        )
        return new(mode,values...)
    end
end

struct Cryogenic_Electrothermal_Model
    nodes::Vector{Cryogenic_Thermal_Node}
    links::Vector{Cryogenic_Thermal_Link}
    active_node_id::String
    bath_temperature_K::Float64
    magnetic_field_T::Float64
    resistance_model::HTS_Transition_Resistance_Model
    circuit::Electrothermal_Circuit
    geometry_hash::String
    material_state_hash::String
    model_hash::String
    qualification_status::Symbol
end

function Cryogenic_Electrothermal_Model(
    nodes::AbstractVector{Cryogenic_Thermal_Node},
    links::AbstractVector{Cryogenic_Thermal_Link},
    active_node_id::AbstractString,
    resistance_model::HTS_Transition_Resistance_Model,
    circuit::Electrothermal_Circuit;
    bath_temperature_K::Real,
    magnetic_field_T::Real=0.0,
    geometry_hash::AbstractString,
    material_state_hash::AbstractString,
    model_hash::AbstractString,
    qualification_status::Symbol=:candidate,
)
    node_vector = Cryogenic_Thermal_Node[nodes...]
    isempty(node_vector) && error("Electrothermal model requires at least one node.")
    identifiers = getfield.(node_vector,:node_id)
    length(unique(identifiers)) == length(identifiers) || error("Thermal node identifiers must be unique.")
    String(active_node_id) in identifiers || error("Active detector node is absent.")
    link_vector = Cryogenic_Thermal_Link[links...]
    for link in link_vector
        link.node_a in identifiers && link.node_b in identifiers || error(
            "Thermal link references an unknown node.",
        )
    end
    bath = Float64(bath_temperature_K)
    field = Float64(magnetic_field_T)
    isfinite(bath) && bath > 0.0 || error("Bath temperature must be positive.")
    isfinite(field) && field ≥ 0.0 || error("Magnetic field must be nonnegative.")
    for value in (geometry_hash,material_state_hash,model_hash)
        isempty(value) && error("Electrothermal lineage hashes cannot be empty.")
    end
    qualification_status in (:verification,:candidate,:qualified) || error(
        "Electrothermal qualification status is invalid.",
    )
    return Cryogenic_Electrothermal_Model(
        node_vector,link_vector,String(active_node_id),bath,field,resistance_model,
        circuit,String(geometry_hash),String(material_state_hash),String(model_hash),
        qualification_status,
    )
end

struct Electrothermal_Energy_Impulse
    time_s::Float64
    node_id::String
    energy_J::Float64
    source_hash::String
    channel::Symbol

    function Electrothermal_Energy_Impulse(
        time_s::Real,
        node_id::AbstractString,
        energy_J::Real;
        source_hash::AbstractString,
        channel::Symbol=:prompt_lattice_heat,
    )
        time = Float64(time_s)
        energy = Float64(energy_J)
        isfinite(time) && time ≥ 0.0 || error("Impulse time must be nonnegative.")
        isfinite(energy) && energy ≥ 0.0 || error("Impulse energy must be nonnegative.")
        isempty(node_id) && error("Impulse node identifier cannot be empty.")
        isempty(source_hash) && error("Impulse source hash cannot be empty.")
        channel in (:prompt_lattice_heat,:delayed_lattice_heat,:joule,:external) || error(
            "Unknown electrothermal impulse channel.",
        )
        return new(time,String(node_id),energy,String(source_hash),channel)
    end
end

struct Cryogenic_Electrothermal_Result
    time_s::Vector{Float64}
    node_ids::Vector{String}
    temperature_K::Matrix{Float64}
    current_A::Vector{Float64}
    active_resistance_ohm::Vector{Float64}
    active_voltage_V::Vector{Float64}
    step_energy_residual_J::Vector{Float64}
    model_hash::String
    source_hashes::Vector{String}
    provenance::Dict{String,String}
end

function _thermal_node_index(model::Cryogenic_Electrothermal_Model)
    return Dict(node.node_id => index for (index,node) in enumerate(model.nodes))
end

function _circuit_current(
    circuit::Electrothermal_Circuit,
    previous_current_A::Float64,
    active_resistance_ohm::Float64,
    dt_s::Float64,
)
    circuit.mode == :current_bias && return circuit.bias_current_A
    total_resistance = circuit.series_resistance_ohm+active_resistance_ohm
    if circuit.inductance_H == 0.0
        total_resistance > 0.0 || error("Voltage-biased circuit has zero total resistance.")
        return circuit.bias_voltage_V/total_resistance
    end
    coefficient = circuit.inductance_H/dt_s
    return (coefficient*previous_current_A+circuit.bias_voltage_V)/
           (coefficient+total_resistance)
end

"""
    simulate_cryogenic_electrothermal(model, time_s; impulses=[])

Semi-implicit lumped thermal-network solve. Heat capacities and conductances are evaluated at the
previous state, while conductive coupling is solved implicitly. Impulse energy, Joule heating,
bath removal, and internal-link cancellation are recorded in a step energy residual. This solver
turns process-partitioned heat into temperature and voltage; it must not consume defect-stored,
recoil-handoff, escaping, or unresolved energy as heat.
"""
function simulate_cryogenic_electrothermal(
    model::Cryogenic_Electrothermal_Model,
    time_s::AbstractVector{<:Real};
    impulses::AbstractVector{Electrothermal_Energy_Impulse}=Electrothermal_Energy_Impulse[],
    energy_balance_rtol::Real=5.0e-3,
    energy_balance_atol_J::Real=1.0e-12,
)
    times = Float64.(time_s)
    length(times) ≥ 2 && all(isfinite,times) && times[1] ≥ 0.0 &&
        all(diff(times) .> 0.0) || error(
        "Electrothermal time grid must be finite and strictly increasing.",
    )
    node_index = _thermal_node_index(model)
    for impulse in impulses
        haskey(node_index,impulse.node_id) || error("Impulse references an unknown thermal node.")
        impulse.time_s ≥ times[1] && impulse.time_s ≤ times[end] || error(
            "Impulse time lies outside the simulation interval.",
        )
    end
    node_count = length(model.nodes)
    active_index = node_index[model.active_node_id]
    temperatures = zeros(Float64,node_count,length(times))
    temperatures[:,1] .= getfield.(model.nodes,:initial_temperature_K)
    currents = zeros(Float64,length(times))
    resistances = zeros(Float64,length(times))
    voltages = zeros(Float64,length(times))
    residuals = zeros(Float64,length(times)-1)
    resistances[1] = hts_resistance(
        model.resistance_model,temperatures[active_index,1],model.magnetic_field_T,
    )
    currents[1] = model.circuit.mode == :current_bias ? model.circuit.bias_current_A :
        _circuit_current(model.circuit,0.0,resistances[1],times[2]-times[1])
    voltages[1] = currents[1]*resistances[1]

    source_hashes = unique(String[impulse.source_hash for impulse in impulses])
    for step in 2:length(times)
        dt = times[step]-times[step-1]
        previous = view(temperatures,:,step-1)
        capacities = [
            property_value(model.nodes[index].heat_capacity_J_K,previous[index])
            for index in 1:node_count
        ]
        bath_conductances = [
            property_value(model.nodes[index].bath_conductance_W_K,previous[index])
            for index in 1:node_count
        ]
        matrix = zeros(Float64,node_count,node_count)
        right_hand_side = zeros(Float64,node_count)
        for index in 1:node_count
            matrix[index,index] += capacities[index]/dt+bath_conductances[index]
            right_hand_side[index] += capacities[index]/dt*previous[index]+
                                      bath_conductances[index]*model.bath_temperature_K
        end
        for link in model.links
            a = node_index[link.node_a]
            b = node_index[link.node_b]
            conductance = property_value(
                link.conductance_W_K,0.5*(previous[a]+previous[b]),
            )
            matrix[a,a] += conductance
            matrix[b,b] += conductance
            matrix[a,b] -= conductance
            matrix[b,a] -= conductance
        end

        active_resistance = hts_resistance(
            model.resistance_model,previous[active_index],model.magnetic_field_T,
        )
        current = _circuit_current(model.circuit,currents[step-1],active_resistance,dt)
        joule_power = current^2*active_resistance
        right_hand_side[active_index] += joule_power

        impulse_energy = zeros(Float64,node_count)
        for impulse in impulses
            if impulse.time_s > times[step-1] && impulse.time_s ≤ times[step]
                impulse_energy[node_index[impulse.node_id]] += impulse.energy_J
            elseif step == 2 && impulse.time_s == times[1]
                impulse_energy[node_index[impulse.node_id]] += impulse.energy_J
            end
        end
        right_hand_side .+= impulse_energy./dt
        next_temperature = matrix\right_hand_side
        all(isfinite,next_temperature) && all(next_temperature .> 0.0) || error(
            "Electrothermal solve produced a nonphysical temperature.",
        )
        temperatures[:,step] .= next_temperature
        currents[step] = current
        resistances[step] = hts_resistance(
            model.resistance_model,next_temperature[active_index],model.magnetic_field_T,
        )
        voltages[step] = current*resistances[step]

        stored_change = sum(
            capacities[index]*(next_temperature[index]-previous[index])
            for index in 1:node_count
        )
        input_energy = sum(impulse_energy)+joule_power*dt
        bath_energy = sum(
            bath_conductances[index]*(next_temperature[index]-model.bath_temperature_K)*dt
            for index in 1:node_count
        )
        residuals[step-1] = input_energy-bath_energy-stored_change
        scale = max(abs(input_energy),abs(bath_energy),abs(stored_change),energy_balance_atol_J)
        abs(residuals[step-1]) ≤ energy_balance_atol_J+energy_balance_rtol*scale || error(
            "Electrothermal step energy balance failed at step $(step-1).",
        )
    end
    return Cryogenic_Electrothermal_Result(
        times,getfield.(model.nodes,:node_id),temperatures,currents,resistances,voltages,
        residuals,model.model_hash,source_hashes,Dict(
            "schema" => "radiant.hts.cryogenic_electrothermal_result/v1",
            "geometry_hash" => model.geometry_hash,
            "material_state_hash" => model.material_state_hash,
            "qualification_status" => string(model.qualification_status),
        ),
    )
end

function peak_active_temperature(
    result::Cryogenic_Electrothermal_Result,
    active_node_id::AbstractString,
)
    index = findfirst(==(String(active_node_id)),result.node_ids)
    isnothing(index) && error("Active node is absent from the electrothermal result.")
    return maximum(view(result.temperature_K,index,:))
end

peak_voltage(result::Cryogenic_Electrothermal_Result) = maximum(result.active_voltage_V)

function recovery_time_s(
    result::Cryogenic_Electrothermal_Result,
    active_node_id::AbstractString;
    temperature_tolerance_K::Real,
)
    index = findfirst(==(String(active_node_id)),result.node_ids)
    isnothing(index) && error("Active node is absent from the electrothermal result.")
    tolerance = Float64(temperature_tolerance_K)
    tolerance ≥ 0.0 || error("Recovery temperature tolerance must be nonnegative.")
    series = view(result.temperature_K,index,:)
    peak_index = argmax(series)
    baseline = series[1]
    for time_index in peak_index:length(series)
        if abs(series[time_index]-baseline) ≤ tolerance
            return result.time_s[time_index]-result.time_s[peak_index]
        end
    end
    return Inf
end

"""Construct a two-node verification model from a fully specified detector definition."""
function cryogenic_model_from_detector(
    detector::HTS_Detector_Definition;
    bath_temperature_K::Real,
    magnetic_field_T::Real=0.0,
    geometry_hash::AbstractString="unbound",
    material_state_hash::AbstractString="unbound",
    model_hash::AbstractString="detector-lumped-model",
)
    is_ready_for_transient_response(detector) || error(
        "Detector definition lacks required electrothermal parameters.",
    )
    thermal = detector.thermal_model
    active = Cryogenic_Thermal_Node(
        "active",detector.active_material_tag,
        constant_cryogenic_property(
            thermal.active_heat_capacity_J_K;
            units="J/K",property_hash="detector-active-C",
        ),
        constant_cryogenic_property(
            0.0;units="W/K",property_hash="active-no-direct-bath",allow_zero=true,
        );
        initial_temperature_K=detector.operating_temperature_K,
    )
    substrate = Cryogenic_Thermal_Node(
        "substrate",detector.substrate_material_tag,
        constant_cryogenic_property(
            thermal.substrate_heat_capacity_J_K;
            units="J/K",property_hash="detector-substrate-C",
        ),
        constant_cryogenic_property(
            thermal.substrate_bath_conductance_W_K;
            units="W/K",property_hash="detector-substrate-bath-G",
        );
        initial_temperature_K=detector.operating_temperature_K,
    )
    link = Cryogenic_Thermal_Link(
        "active","substrate",
        constant_cryogenic_property(
            thermal.active_substrate_conductance_W_K;
            units="W/K",property_hash="detector-interface-G",
        );
        interface_hash="detector-interface",
    )
    resistance = HTS_Transition_Resistance_Model(
        normal_resistance_ohm=thermal.normal_resistance_ohm,
        residual_resistance_ohm=0.0,
        critical_temperature_K=thermal.transition_temperature_K,
        transition_width_K=thermal.transition_width_K,
        model_hash="detector-transition-model",
        qualification_status=:verification,
    )
    circuit = Electrothermal_Circuit(
        mode=:current_bias,
        bias_current_A=max(detector.bias_current_A,1.0e-9),
        inductance_H=isnothing(thermal.inductance_H) ? 0.0 : thermal.inductance_H,
    )
    return Cryogenic_Electrothermal_Model(
        [active,substrate],[link],"active",resistance,circuit;
        bath_temperature_K=bath_temperature_K,magnetic_field_T=magnetic_field_T,
        geometry_hash=geometry_hash,material_state_hash=material_state_hash,
        model_hash=model_hash,qualification_status=:verification,
    )
end
