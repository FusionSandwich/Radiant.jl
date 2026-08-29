const SUBKEV_PARTITION_CHANNELS = (
    :ionization,
    :electronic_excitation,
    :prompt_phonon,
    :athermal_phonon,
    :quasiparticle,
    :defect_storage,
    :optical_emission,
    :electron_escape,
)

"""
    SubkeV_Thermalization_Kernel

Material- and state-specific partition model below Radiant's transport cutoff. Fractions are
tabulated versus electron kinetic energy and must sum to one at every grid point. The kernel is a
handoff model, not a universal YBCO track-structure law; `qualification_status=:qualified` requires
independent material-specific evidence and a bound state hash.
"""
struct SubkeV_Thermalization_Kernel
    material_tag::String
    temperature_K::Float64
    energy_grid_eV::Vector{Float64}
    fractions::Dict{Symbol,Vector{Float64}}
    model_id::String
    model_hash::String
    material_state_hash::String
    qualification_status::Symbol
    metadata::Dict{String,String}

    function SubkeV_Thermalization_Kernel(
        material_tag::AbstractString,
        temperature_K::Real,
        energy_grid_eV::AbstractVector{<:Real},
        fractions::AbstractDict;
        model_id::AbstractString,
        model_hash::AbstractString,
        material_state_hash::AbstractString,
        qualification_status::Symbol=:candidate,
        metadata::AbstractDict=Dict{String,String}(),
        closure_tolerance::Real=1.0e-10,
    )
        isempty(material_tag) && error("Sub-keV material tag cannot be empty.")
        temperature = Float64(temperature_K)
        isfinite(temperature) && temperature > 0.0 || error(
            "Sub-keV model temperature must be finite and positive.",
        )
        energies = Float64.(energy_grid_eV)
        length(energies) ≥ 2 && all(isfinite,energies) && all(energies .> 0.0) &&
            all(diff(energies) .> 0.0) || error(
            "Sub-keV energy grid must be finite, positive, and strictly increasing.",
        )
        qualification_status in (:verification,:candidate,:qualified) || error(
            "Sub-keV qualification status is invalid.",
        )
        isempty(model_id) && error("Sub-keV model identifier cannot be empty.")
        isempty(model_hash) && error("Sub-keV model hash cannot be empty.")
        isempty(material_state_hash) && error("Sub-keV material-state hash cannot be empty.")

        fraction_dictionary = Dict{Symbol,Vector{Float64}}()
        for channel in SUBKEV_PARTITION_CHANNELS
            haskey(fractions,channel) || error(
                "Sub-keV partition is missing channel $(channel).",
            )
            values = Float64.(fractions[channel])
            length(values) == length(energies) || error(
                "Sub-keV channel $(channel) does not match the energy grid.",
            )
            all(value -> isfinite(value) && value ≥ 0.0 && value ≤ 1.0,values) || error(
                "Sub-keV channel fractions must lie in [0,1].",
            )
            fraction_dictionary[channel] = values
        end
        for index in eachindex(energies)
            total = sum(fraction_dictionary[channel][index] for channel in SUBKEV_PARTITION_CHANNELS)
            isapprox(total,1.0;rtol=closure_tolerance,atol=closure_tolerance) || error(
                "Sub-keV channel fractions do not close at energy index $(index): $(total).",
            )
        end
        metadata_string = Dict{String,String}()
        for (key,value) in metadata
            metadata_string[string(key)] = string(value)
        end
        return new(
            String(material_tag),temperature,energies,fraction_dictionary,String(model_id),
            String(model_hash),String(material_state_hash),qualification_status,metadata_string,
        )
    end
end

struct SubkeV_Thermalization_Result
    incident_energy_eV::Float64
    base_partition::Energy_Partition
    athermal_phonon_eV::Float64
    quasiparticle_eV::Float64
    model_hash::String
    material_state_hash::String
    qualification_status::Symbol
    metadata::Dict{String,String}
end

function _subkev_interpolate(
    energy_grid::Vector{Float64},
    values::Vector{Float64},
    energy_eV::Float64,
)
    energy_eV < energy_grid[1] && error("Sub-keV event lies below the model energy range.")
    energy_eV > energy_grid[end] && error("Sub-keV event lies above the model energy range.")
    energy_eV == energy_grid[end] && return values[end]
    upper = searchsortedfirst(energy_grid,energy_eV)
    upper == 1 && return values[1]
    energy_grid[upper] == energy_eV && return values[upper]
    lower = upper-1
    fraction = (energy_eV-energy_grid[lower])/(energy_grid[upper]-energy_grid[lower])
    return (1.0-fraction)*values[lower]+fraction*values[upper]
end

function subkev_partition_fractions(
    kernel::SubkeV_Thermalization_Kernel,
    energy_eV::Real,
)
    energy = Float64(energy_eV)
    isfinite(energy) && energy > 0.0 || error("Sub-keV event energy must be positive.")
    output = Dict{Symbol,Float64}()
    for channel in SUBKEV_PARTITION_CHANNELS
        output[channel] = _subkev_interpolate(
            kernel.energy_grid_eV,kernel.fractions[channel],energy,
        )
    end
    total = sum(values(output))
    isapprox(total,1.0;rtol=1.0e-9,atol=1.0e-10) || error(
        "Interpolated sub-keV partition does not close: $(total).",
    )
    return output
end

"""Partition one sub-keV electron event without silently converting non-equilibrium energy to heat."""
function thermalize_subkev_event(
    kernel::SubkeV_Thermalization_Kernel,
    energy_eV::Real;
    label::AbstractString="sub-keV-electron",
    metadata::AbstractDict=Dict{String,String}(),
)
    energy = Float64(energy_eV)
    fractions = subkev_partition_fractions(kernel,energy)
    athermal = energy*fractions[:athermal_phonon]
    quasiparticle = energy*fractions[:quasiparticle]
    unresolved_non_equilibrium = athermal+quasiparticle
    combined_metadata = Dict{String,String}(
        "subkev_model_id" => kernel.model_id,
        "subkev_model_hash" => kernel.model_hash,
        "material_state_hash" => kernel.material_state_hash,
        "athermal_and_quasiparticle_bookkeeping" => "stored-in-result",
    )
    for (key,value) in metadata
        combined_metadata[string(key)] = string(value)
    end
    base = Energy_Partition(
        label=label,
        available_energy_eV=energy,
        ionization_eV=energy*fractions[:ionization],
        electronic_excitation_eV=energy*fractions[:electronic_excitation],
        prompt_lattice_heat_eV=energy*fractions[:prompt_phonon],
        nuclear_recoil_handoff_eV=0.0,
        defect_stored_eV=energy*fractions[:defect_storage],
        optical_emission_eV=energy*fractions[:optical_emission],
        escaping_particle_eV=energy*fractions[:electron_escape],
        cutoff_handoff_eV=0.0,
        unresolved_eV=unresolved_non_equilibrium,
        metadata=combined_metadata,
    )
    is_energy_partition_closed(base) || error("Sub-keV energy partition failed closure.")
    return SubkeV_Thermalization_Result(
        energy,base,athermal,quasiparticle,kernel.model_hash,kernel.material_state_hash,
        kernel.qualification_status,combined_metadata,
    )
end

function get_non_equilibrium_energy(this::SubkeV_Thermalization_Result)
    return this.athermal_phonon_eV+this.quasiparticle_eV
end

function get_resolved_subkev_energy(this::SubkeV_Thermalization_Result)
    return get_accounted_energy(this.base_partition)-this.base_partition.unresolved_eV+
           get_non_equilibrium_energy(this)
end

function is_subkev_partition_closed(
    this::SubkeV_Thermalization_Result;
    rtol::Real=1.0e-10,
    atol_eV::Real=1.0e-8,
)
    return isapprox(
        get_resolved_subkev_energy(this),this.incident_energy_eV;
        rtol=rtol,atol=atol_eV,
    )
end

"""One exponential non-equilibrium reservoir and its eventual branching."""
struct NonEquilibrium_Decay_Channel
    source::Symbol
    lifetime_s::Float64
    heat_fraction::Float64
    escape_fraction::Float64
    optical_fraction::Float64

    function NonEquilibrium_Decay_Channel(
        source::Symbol,
        lifetime_s::Real;
        heat_fraction::Real,
        escape_fraction::Real=0.0,
        optical_fraction::Real=0.0,
        closure_tolerance::Real=1.0e-10,
    )
        source in (:athermal_phonon,:quasiparticle) || error(
            "Non-equilibrium source must be :athermal_phonon or :quasiparticle.",
        )
        lifetime = Float64(lifetime_s)
        fractions = Float64[heat_fraction,escape_fraction,optical_fraction]
        isfinite(lifetime) && lifetime > 0.0 || error(
            "Non-equilibrium lifetime must be positive.",
        )
        all(value -> isfinite(value) && value ≥ 0.0, fractions) || error(
            "Non-equilibrium branching fractions must be nonnegative.",
        )
        isapprox(sum(fractions),1.0;rtol=closure_tolerance,atol=closure_tolerance) || error(
            "Non-equilibrium branching fractions must sum to one.",
        )
        return new(source,lifetime,fractions[1],fractions[2],fractions[3])
    end
end

struct NonEquilibrium_Thermalization_Model
    athermal_phonon::NonEquilibrium_Decay_Channel
    quasiparticle::NonEquilibrium_Decay_Channel
    model_hash::String
    qualification_status::Symbol

    function NonEquilibrium_Thermalization_Model(
        athermal_phonon::NonEquilibrium_Decay_Channel,
        quasiparticle::NonEquilibrium_Decay_Channel;
        model_hash::AbstractString,
        qualification_status::Symbol=:candidate,
    )
        athermal_phonon.source == :athermal_phonon || error(
            "Athermal channel has the wrong source identifier.",
        )
        quasiparticle.source == :quasiparticle || error(
            "Quasiparticle channel has the wrong source identifier.",
        )
        isempty(model_hash) && error("Non-equilibrium model hash cannot be empty.")
        qualification_status in (:verification,:candidate,:qualified) || error(
            "Non-equilibrium qualification status is invalid.",
        )
        return new(athermal_phonon,quasiparticle,String(model_hash),qualification_status)
    end
end

function non_equilibrium_release(
    result::SubkeV_Thermalization_Result,
    model::NonEquilibrium_Thermalization_Model,
    time_s::Real,
)
    time = Float64(time_s)
    isfinite(time) && time ≥ 0.0 || error("Thermalization time must be nonnegative.")
    heat = 0.0
    escaped = 0.0
    optical = 0.0
    remaining = 0.0
    for (energy,channel) in (
        (result.athermal_phonon_eV,model.athermal_phonon),
        (result.quasiparticle_eV,model.quasiparticle),
    )
        released = energy*(1.0-exp(-time/channel.lifetime_s))
        remaining += energy-released
        heat += released*channel.heat_fraction
        escaped += released*channel.escape_fraction
        optical += released*channel.optical_fraction
    end
    return (
        heat_eV=heat,
        escaped_eV=escaped,
        optical_eV=optical,
        remaining_non_equilibrium_eV=remaining,
    )
end

function non_equilibrium_power(
    result::SubkeV_Thermalization_Result,
    model::NonEquilibrium_Thermalization_Model,
    time_s::Real,
)
    time = Float64(time_s)
    isfinite(time) && time ≥ 0.0 || error("Thermalization time must be nonnegative.")
    heat_power = 0.0
    escape_power = 0.0
    optical_power = 0.0
    for (energy,channel) in (
        (result.athermal_phonon_eV,model.athermal_phonon),
        (result.quasiparticle_eV,model.quasiparticle),
    )
        release_rate = energy/channel.lifetime_s*exp(-time/channel.lifetime_s)
        heat_power += release_rate*channel.heat_fraction
        escape_power += release_rate*channel.escape_fraction
        optical_power += release_rate*channel.optical_fraction
    end
    return (
        heat_eV_per_s=heat_power,
        escape_eV_per_s=escape_power,
        optical_eV_per_s=optical_power,
    )
end

function subkev_model_is_production_ready(kernel::SubkeV_Thermalization_Kernel)
    return kernel.qualification_status == :qualified &&
           kernel.model_hash != "unbound" &&
           kernel.material_state_hash != "unbound"
end

"""Synthetic kernel for code verification; it is not YBCO material data."""
function synthetic_subkev_kernel_fixture()
    energies = [10.0,100.0,1000.0]
    fractions = Dict{Symbol,Vector{Float64}}(
        :ionization => [0.10,0.20,0.30],
        :electronic_excitation => [0.10,0.10,0.10],
        :prompt_phonon => [0.20,0.20,0.20],
        :athermal_phonon => [0.20,0.15,0.10],
        :quasiparticle => [0.20,0.15,0.10],
        :defect_storage => [0.05,0.05,0.05],
        :optical_emission => [0.05,0.05,0.05],
        :electron_escape => [0.10,0.10,0.10],
    )
    return SubkeV_Thermalization_Kernel(
        "synthetic-material",20.0,energies,fractions;
        model_id="synthetic-subkev-verification",
        model_hash="synthetic-subkev-kernel",
        material_state_hash="synthetic-state",
        qualification_status=:verification,
        metadata=Dict("classification" => "verification-only"),
    )
end
