"""
    Multigroup_Cross_Sections

Multigroup transport data for one incoming particle and material. `response_channels` is a generic
core hook used by the temporary HTS add-on to retain interaction-resolved deposition, charge,
absorption, stopping-power, and momentum-transfer contributions. Channel keys and values are
copied on insertion so downstream scoring cannot mutate the transport library silently.
"""
mutable struct Multigroup_Cross_Sections
    number_of_groups               ::Int64
    total                          ::Union{Missing,Vector{Float64}}
    absorption                     ::Union{Missing,Vector{Float64}}
    boundary_stopping_powers       ::Union{Missing,Vector{Float64}}
    stopping_powers                ::Union{Missing,Vector{Float64}}
    momentum_transfer              ::Union{Missing,Vector{Float64}}
    energy_deposition              ::Union{Missing,Vector{Float64}}
    charge_deposition              ::Union{Missing,Vector{Float64}}
    scattering                     ::Vector{Array{Float64,3}}
    response_channels              ::Dict{String,Vector{Float64}}

    function Multigroup_Cross_Sections(number_of_groups::Int64)
        number_of_groups ≥ 1 || error("Multigroup data require at least one energy group.")
        this = new()
        this.number_of_groups = number_of_groups
        this.total = missing
        this.absorption = missing
        this.boundary_stopping_powers = missing
        this.stopping_powers = missing
        this.momentum_transfer = missing
        this.energy_deposition = missing
        this.charge_deposition = missing
        this.scattering = Vector{Array{Float64,3}}()
        this.response_channels = Dict{String,Vector{Float64}}()
        return this
    end
end

function _validate_group_vector(
    this::Multigroup_Cross_Sections,
    values::Vector{Float64},
    label::AbstractString;
    allow_cutoff::Bool=false,
)
    expected = allow_cutoff ? (this.number_of_groups,this.number_of_groups+1) :
                              (this.number_of_groups,)
    length(values) in expected || error(
        "The $(label) vector length does not match the number of groups.",
    )
    all(isfinite,values) || error("The $(label) vector contains non-finite values.")
    return values
end

function set_total(this::Multigroup_Cross_Sections,total::Vector{Float64})
    this.total = copy(_validate_group_vector(this,total,"total cross section"))
end

function set_absorption(this::Multigroup_Cross_Sections,absorption::Vector{Float64})
    this.absorption = copy(_validate_group_vector(this,absorption,"absorption cross section"))
end

function set_boundary_stopping_powers(
    this::Multigroup_Cross_Sections,
    boundary_stopping_powers::Vector{Float64},
)
    length(boundary_stopping_powers) == this.number_of_groups+1 || error(
        "Boundary stopping powers require number_of_groups+1 values.",
    )
    all(isfinite,boundary_stopping_powers) || error(
        "Boundary stopping powers contain non-finite values.",
    )
    this.boundary_stopping_powers = copy(boundary_stopping_powers)
end

function set_stopping_powers(this::Multigroup_Cross_Sections,stopping_powers::Vector{Float64})
    this.stopping_powers = copy(
        _validate_group_vector(this,stopping_powers,"stopping power"),
    )
end

function set_momentum_transfer(
    this::Multigroup_Cross_Sections,
    momentum_transfer::Vector{Float64},
)
    this.momentum_transfer = copy(
        _validate_group_vector(this,momentum_transfer,"momentum transfer"),
    )
end

function set_energy_deposition(
    this::Multigroup_Cross_Sections,
    energy_deposition::Vector{Float64},
)
    length(energy_deposition) == this.number_of_groups+1 || error(
        "Energy deposition requires number_of_groups+1 values.",
    )
    all(isfinite,energy_deposition) || error(
        "Energy-deposition coefficients contain non-finite values.",
    )
    this.energy_deposition = copy(energy_deposition)
end

function set_charge_deposition(
    this::Multigroup_Cross_Sections,
    charge_deposition::Vector{Float64},
)
    length(charge_deposition) == this.number_of_groups+1 || error(
        "Charge deposition requires number_of_groups+1 values.",
    )
    all(isfinite,charge_deposition) || error(
        "Charge-deposition coefficients contain non-finite values.",
    )
    this.charge_deposition = copy(charge_deposition)
end

function set_scattering(
    this::Multigroup_Cross_Sections,
    scattering::Array{Float64,3},
)
    size(scattering,1) == this.number_of_groups || error(
        "Scattering incoming-group dimension is inconsistent.",
    )
    all(isfinite,scattering) || error("Scattering data contain non-finite values.")
    push!(this.scattering,copy(scattering))
end

"""
    add_response_channel!(this, key, values; accumulate=true)

Store one generic response vector. Values may contain one coefficient per energy group or one
additional cutoff coefficient. Repeated kernel contributions accumulate by default; shape mismatch
is rejected rather than padded or truncated.
"""
function add_response_channel!(
    this::Multigroup_Cross_Sections,
    key::AbstractString,
    values::AbstractVector{<:Real};
    accumulate::Bool=true,
)
    key_string = String(key)
    isempty(key_string) && error("Response-channel key cannot be empty.")
    occursin("|",key_string) || error(
        "Response-channel key must use quantity|process schema.",
    )
    vector = Float64.(values)
    length(vector) in (this.number_of_groups,this.number_of_groups+1) || error(
        "Response-channel length must be number_of_groups or number_of_groups+1.",
    )
    all(isfinite,vector) || error("Response-channel coefficients must be finite.")
    if haskey(this.response_channels,key_string)
        accumulate || error("Response channel already exists: $(key_string).")
        length(this.response_channels[key_string]) == length(vector) || error(
            "Cannot accumulate response channels with different group dimensions.",
        )
        this.response_channels[key_string] .+= vector
    else
        this.response_channels[key_string] = copy(vector)
    end
    return this
end

function set_response_channel!(
    this::Multigroup_Cross_Sections,
    key::AbstractString,
    values::AbstractVector{<:Real},
)
    key_string = String(key)
    haskey(this.response_channels,key_string) && delete!(this.response_channels,key_string)
    return add_response_channel!(this,key_string,values;accumulate=false)
end

has_response_channel(this::Multigroup_Cross_Sections,key::AbstractString) =
    haskey(this.response_channels,String(key))

function get_response_channel(this::Multigroup_Cross_Sections,key::AbstractString)
    key_string = String(key)
    haskey(this.response_channels,key_string) || error(
        "Multigroup response channel is unavailable: $(key_string).",
    )
    return copy(this.response_channels[key_string])
end

function get_response_channel_keys(
    this::Multigroup_Cross_Sections;
    quantity::Union{Nothing,AbstractString}=nothing,
)
    keys_all = sort(collect(keys(this.response_channels)))
    isnothing(quantity) && return keys_all
    prefix = string(quantity,"|")
    return [key for key in keys_all if startswith(key,prefix)]
end

function get_response_channels(this::Multigroup_Cross_Sections)
    return Dict(key => copy(value) for (key,value) in this.response_channels)
end

function sum_response_channels(
    this::Multigroup_Cross_Sections,
    quantity::AbstractString,
)
    keys_selected = get_response_channel_keys(this;quantity=quantity)
    isempty(keys_selected) && error("No response channels exist for quantity $(quantity).")
    lengths = unique(length(this.response_channels[key]) for key in keys_selected)
    length(lengths) == 1 || error("Response channels have inconsistent group dimensions.")
    total = zeros(Float64,first(lengths))
    for key in keys_selected
        total .+= this.response_channels[key]
    end
    return total
end

function get_total(this::Multigroup_Cross_Sections)
    ismissing(this.total) && error("Unable to get multigroup total cross sections. Missing data.")
    return this.total
end

function get_absorption(this::Multigroup_Cross_Sections)
    ismissing(this.absorption) && error(
        "Unable to get multigroup absorption cross sections. Missing data.",
    )
    return this.absorption
end

function get_scattering(this::Multigroup_Cross_Sections,index_particle_out::Int64)
    1 ≤ index_particle_out ≤ length(this.scattering) || error(
        "Outgoing-particle scattering index is out of range.",
    )
    return this.scattering[index_particle_out]
end

function get_boundary_stopping_powers(this::Multigroup_Cross_Sections)
    ismissing(this.boundary_stopping_powers) && error("Boundary stopping powers are missing.")
    return this.boundary_stopping_powers
end

function get_stopping_powers(this::Multigroup_Cross_Sections)
    ismissing(this.stopping_powers) && error("Stopping powers are missing.")
    return this.stopping_powers
end

function get_momentum_transfer(this::Multigroup_Cross_Sections)
    ismissing(this.momentum_transfer) && error("Momentum-transfer data are missing.")
    return this.momentum_transfer
end

function get_energy_deposition(this::Multigroup_Cross_Sections)
    ismissing(this.energy_deposition) && error("Energy-deposition data are missing.")
    return this.energy_deposition
end

function get_charge_deposition(this::Multigroup_Cross_Sections)
    ismissing(this.charge_deposition) && error("Charge-deposition data are missing.")
    return this.charge_deposition
end
