"""
    Energy_Partition

Mechanism-resolved energy accounting for one event, voxel, layer, or response bin.

Radiant's conventional energy-deposition score is not automatically identical to instantaneous
cryogenic heat. This ledger separates ionization, electronic excitation, prompt lattice heat,
nuclear-recoil energy handed to a recoil code, energy stored in metastable defects, optical
emission, escaping-particle energy, sub-cutoff handoff, and an explicitly unresolved remainder.
All values are in eV and refer to the same normalization basis.
"""
struct Energy_Partition
    label::String
    available_energy_eV::Float64
    ionization_eV::Float64
    electronic_excitation_eV::Float64
    prompt_lattice_heat_eV::Float64
    nuclear_recoil_handoff_eV::Float64
    defect_stored_eV::Float64
    optical_emission_eV::Float64
    escaping_particle_eV::Float64
    cutoff_handoff_eV::Float64
    unresolved_eV::Float64
    metadata::Dict{String,String}
end

function Energy_Partition(;
    label::AbstractString="energy-partition",
    available_energy_eV::Real=0.0,
    ionization_eV::Real=0.0,
    electronic_excitation_eV::Real=0.0,
    prompt_lattice_heat_eV::Real=0.0,
    nuclear_recoil_handoff_eV::Real=0.0,
    defect_stored_eV::Real=0.0,
    optical_emission_eV::Real=0.0,
    escaping_particle_eV::Real=0.0,
    cutoff_handoff_eV::Real=0.0,
    unresolved_eV::Real=0.0,
    metadata::AbstractDict=Dict{String,String}(),
)
    values = Float64[
        available_energy_eV,
        ionization_eV,
        electronic_excitation_eV,
        prompt_lattice_heat_eV,
        nuclear_recoil_handoff_eV,
        defect_stored_eV,
        optical_emission_eV,
        escaping_particle_eV,
        cutoff_handoff_eV,
        unresolved_eV,
    ]
    any(x -> !isfinite(x) || x < 0.0,values) && error(
        "Energy-partition terms must be finite and nonnegative.",
    )
    metadata_string = Dict{String,String}()
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    return Energy_Partition(
        String(label),values[1],values[2],values[3],values[4],values[5],values[6],
        values[7],values[8],values[9],values[10],metadata_string,
    )
end

function get_accounted_energy(this::Energy_Partition)
    return this.ionization_eV +
           this.electronic_excitation_eV +
           this.prompt_lattice_heat_eV +
           this.nuclear_recoil_handoff_eV +
           this.defect_stored_eV +
           this.optical_emission_eV +
           this.escaping_particle_eV +
           this.cutoff_handoff_eV +
           this.unresolved_eV
end

get_energy_partition_residual(this::Energy_Partition) =
    this.available_energy_eV - get_accounted_energy(this)

function get_relative_energy_partition_residual(this::Energy_Partition;atol_eV::Real=1.0e-12)
    denominator = max(abs(this.available_energy_eV),Float64(atol_eV))
    return abs(get_energy_partition_residual(this))/denominator
end

function is_energy_partition_closed(
    this::Energy_Partition;
    rtol::Real=1.0e-10,
    atol_eV::Real=1.0e-8,
)
    return isapprox(
        get_accounted_energy(this),this.available_energy_eV;
        rtol=rtol,atol=atol_eV,
    )
end

function is_energy_partition_resolved(this::Energy_Partition;atol_eV::Real=1.0e-8)
    return is_energy_partition_closed(this;atol_eV=atol_eV) && this.unresolved_eV ≤ atol_eV
end

function get_prompt_heat_fraction(this::Energy_Partition)
    this.available_energy_eV == 0.0 && return 0.0
    return this.prompt_lattice_heat_eV/this.available_energy_eV
end
