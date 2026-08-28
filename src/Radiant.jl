module Radiant

    #----
    # External packages
    #----
    import Base.println
    using Printf: @sprintf
    using LinearAlgebra
    using HDF5
    using JLD2
    using SHA
    using SpecialFunctions

    # Load version-gated compatibility methods before any source files are evaluated.
    include("tools/julia_compat.jl")

    #----
    # Include files
    #----
    radiant_src = Dict{String,Vector{String}}()
    radiant_src["cross_sections/"] = [
        "read_fmac_m.jl",
        "write_fmac_m.jl",
        "read_matxs.jl",
        "write_matxs.jl",
        "generate_cross_sections.jl",
        "custom_cross_sections.jl",
        "energy_group_structure.jl",
        "atomic_electron_cascades.jl",
        "multigroup.jl",
        "feed.jl",
        "mean_excitation_energy.jl",
        "atomic_weight.jl",
        "isotopic_composition.jl",
        "density.jl",
        "nuclei_density.jl",
        "plasma_energy.jl",
        "electron_subshells.jl",
        "orbital_compton_profiles.jl",
        "resonance_energy.jl",
        "fermi_density_effect.jl",
        "atomic_number.jl",
        "transport_correction.jl",
        "angular_fokker_planck_decomposition.jl",
        "bethe.jl",
        "moller.jl",
        "bhabha.jl",
        "kawrakow_correction.jl",
        "moliere_screening.jl",
        "mott.jl",
        "sommerfield.jl",
        "poskus.jl",
        "ratio_positron_electron_bremsstrahlung.jl",
        "seltzer_berger.jl",
        "klein_nishina.jl",
        "waller_hartree.jl",
        "impulse_approximation.jl",
        "biggs_lighthill.jl",
        "photoelectric_per_subshell.jl",
        "sauter.jl",
        "rayleigh.jl",
        "high_energy_coulomb_correction.jl",
        "baro.jl",
        "heitler.jl",
        "relaxation.jl",
        "interaction_interdependances.jl",
        "inelastic_collision_heavy_particle.jl",
        "soft_catastrophic_cutoff.jl",
        "endf_reading.jl",
        "scattering_kinematics.jl",
        "elastic_scattering_endf.jl"
    ]
    radiant_src["particle_transport/"] = [
        "geometry.jl",
        "volume_source.jl",
        "surface_source.jl",
        "scattering_source.jl",
        "particle_source.jl",
        "fokker_planck_source.jl",
        "fokker_planck_scattering_matrix.jl",
        "fokker_planck_finite_difference.jl",
        "fokker_planck_differential_quadrature.jl",
        "fokker_planck_galerkin.jl",
        "fokker_planck_finite_difference_gn.jl",
        "electromagnetic_scattering_matrix.jl",
        "transport.jl",
        "sn_inner_pass.jl",
        "sn_flux.jl",
        "sn_one_speed.jl",
        "sn_sweep_1D.jl",
        "sn_sweep_2D.jl",
        "sn_sweep_3D.jl",
        "sn_1D_bte.jl",
        "sn_1D_bfp.jl",
        "sn_2D_bte.jl",
        "sn_2D_bfp.jl",
        "sn_3D_bte.jl",
        "sn_3D_bfp.jl",
        "adaptive.jl",
        "scheme_weights.jl",
        "map_moments.jl",
        "constant_linear.jl",
        "angular_polynomial_basis.jl",
        "hts_source_projection.jl",
        "energy_deposition.jl",
        "charge_deposition.jl",
        "flux.jl",
        "restricted_to_full_domain_matrix.jl",
        "gn_patch_geometry.jl",
        "gn_flux.jl",
        "gn_inner_pass.jl",
        "gn_one_speed.jl",
        "gn_sweep_1D.jl",
        "gn_sweep_2D.jl",
        "gn_sweep_3D.jl",
        "gn_1D_bte.jl",
        "gn_1D_bfp.jl",
        "gn_2D_bte.jl",
        "gn_2D_bfp.jl",
        "gn_3D_bte.jl",
        "gn_3D_bfp.jl",
        "gn_weights.jl",
        "cp_collision_probabilities.jl",
        "cp_surface_probabilities.jl",
        "cp_sweep_1D.jl",
        "cp_flux.jl"
    ]
    radiant_src["structures/"] = [
        "Particle.jl",
        "Source_Normalization.jl",
        "Boundary_Angular_Current_Source.jl",
        "Anisotropic_Volume_Source.jl",
        "Transport_Balance.jl",
        "Transport_Ownership_Map.jl",
        "Layer_Definition.jl",
        "Tape_Stack_Definition.jl",
        "Energy_Partition.jl",
        "Fibre_Diagnostic_Definition.jl",
        "HTS_Detector_Definition.jl",
        "Physics_Coverage_Register.jl",
        "Protected_Response_Registry.jl",
        "Interaction.jl",
        "Elastic_Collision.jl",
        "Elastic_Scattering.jl",
        "Inelastic_Collision.jl",
        "Bremsstrahlung.jl",
        "Compton.jl",
        "Pair_Production.jl",
        "Photoelectric.jl",
        "Annihilation.jl",
        "Rayleigh.jl",
        "Relaxation.jl",
        "Multigroup_Cross_Sections.jl",
        "Material.jl",
        "Material_List.jl",
        "Cross_Sections.jl",
        "Geometry.jl",
        "SN.jl",
        "DPN.jl",
        "GN.jl",
        "CP.jl",
        "Solvers.jl",
        "Surface_Source.jl",
        "Volume_Source.jl",
        "Source.jl",
        "Sources.jl",
        "Fixed_Sources.jl",
        "Flux_Per_Particle.jl",
        "Flux.jl",
        "Electromagnetic_Field.jl",
        "Computation_Unit.jl"
    ]
    radiant_src["tools/"] = [
        "quadrature.jl",
        "gauss_legendre.jl",
        "gauss_lobatto.jl",
        "gauss_legendre_chebychev.jl",
        "lebedev.jl",
        "carlson.jl",
        "legendre_polynomials.jl",
        "jacobi_polynomials.jl",
        "ferrer_associated_legendre.jl",
        "real_spherical_harmonics.jl",
        "spherical_triangle.jl",
        "heaviside.jl",
        "interpolation.jl",
        "integrals.jl",
        "spline.jl",
        "voronoi.jl",
        "newton_bisection.jl",
        "tanh_sinh_integral.jl",
        "one_space.jl",
        "cache.jl",
        "find_package_root.jl",
        "python_method_notation.jl",
        "factorial_factor.jl",
        "double_factorial.jl",
        "krylov_state.jl",
        "gmres.jl",
        "bicgstab.jl",
        "anderson.jl",
        "livolant.jl"
    ]
    radiant_src["interchange/"] = [
        "source_hdf5.jl"
    ]
    for folder in ["structures/","tools/","cross_sections/","particle_transport/","interchange/"]
        for file in radiant_src[folder]
            include(string(folder,file))
        end
    end

    # Export all generated material constructors (from structures/Material_list.jl).
    if isdefined(@__MODULE__, :MATERIAL_CONSTRUCTORS)
        eval(Expr(:export, getfield(@__MODULE__, :MATERIAL_CONSTRUCTORS)...))
    end

    #----
    # Export objects
    #----
    export Particle,Photon,Electron,Positron,Proton,Antiproton,Alpha,Muon,Antimuon
    export Elastic_Collision,Elastic_Scattering,Inelastic_Collision,Bremsstrahlung,Compton
    export Pair_Production,Photoelectric,Annihilation,Rayleigh,Relaxation,Fluorescence,Auger
    export Material,Cross_Sections,Geometry,SN,Solvers,Surface_Source,Volume_Source
    export Fixed_Sources,Computation_Unit,DPN,GN,CP,Electromagnetic_Field
    export Discrete_Ordinates

    # HTS deterministic-coupling foundations.
    export Source_Normalization,get_physical_scale,get_duration,apply_normalization
    export Abstract_Radiant_Source,Boundary_Angular_Current_Source
    export get_incoming_current_density,get_incoming_current,get_total_incoming_current
    export assert_current_closure,boundary_source_from_directional_current
    export Anisotropic_Volume_Source,get_volume_source_rate,get_source_normalization
    export Boundary_Projection_Receipt,Volume_Projection_Receipt
    export project_boundary_source,project_volume_source,get_projection_receipts
    export Transport_Balance,get_particle_residual,get_energy_residual,get_charge_residual
    export get_relative_particle_residual,get_relative_energy_residual,is_balanced
    export Transport_Ownership_Record,Transport_Ownership_Map,validate_ownership
    export get_production_owner
    export Layer_Definition,Tape_Stack_Definition,get_total_thickness,get_layer_boundaries
    export get_layer_index,verification_eight_layer_stack

    # High-fidelity deposition, diagnostic, detector, and qualification contracts.
    export Energy_Partition,get_accounted_energy,get_energy_partition_residual
    export get_relative_energy_partition_residual,is_energy_partition_closed
    export is_energy_partition_resolved,get_prompt_heat_fraction
    export Fibre_Radial_Layer,Fibre_Diagnostic_Definition,get_outer_radius
    export get_cross_sectional_area,get_fibre_layer_index,requires_optical_response_model
    export Detector_Converter_Layer,Detector_Thermal_Model,HTS_Detector_Definition
    export is_ready_for_transient_response,get_active_volume,verification_b10_ybco_detector
    export Physics_Coverage_Record,Physics_Coverage_Register,get_physics_record
    export validate_physics_coverage,default_hts_physics_coverage
    export Protected_Response,Protected_Response_Registry,get_protected_response
    export response_is_converged,default_hts_protected_responses

    # Hash-bound OpenMC/OpenSn/Radiant source interchange.
    export BOUNDARY_SOURCE_HDF5_SCHEMA,VOLUME_SOURCE_HDF5_SCHEMA
    export write_boundary_angular_current_hdf5,read_boundary_angular_current_hdf5
    export write_anisotropic_volume_source_hdf5,read_anisotropic_volume_source_hdf5

    #----
    # Execution of Radiant script files
    #----
    export @radiant_input
    macro radiant_input()
        return quote
            if abspath(PROGRAM_FILE) == @__FILE__
                using Radiant
                Radiant.run_script(@__FILE__)
                exit()
            end
        end
    end

    function run_script(script::AbstractString)
        isfile(script) || error("Input script not found: $(script)")
        include(script)
    end
end
