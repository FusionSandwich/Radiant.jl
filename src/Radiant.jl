module Radiant

    import Base.println
    using Printf: @sprintf
    using LinearAlgebra
    using HDF5
    using JLD2
    using Random
    using SHA
    using SpecialFunctions

    include("tools/julia_compat.jl")

    radiant_src = Dict{String,Vector{String}}()
    radiant_src["cross_sections/"] = [
        "read_fmac_m.jl","write_fmac_m.jl","read_matxs.jl","write_matxs.jl",
        "generate_cross_sections.jl","custom_cross_sections.jl",
        "energy_group_structure.jl","atomic_electron_cascades.jl","multigroup.jl",
        "feed.jl","mean_excitation_energy.jl","atomic_weight.jl",
        "isotopic_composition.jl","density.jl","nuclei_density.jl","plasma_energy.jl",
        "electron_subshells.jl","orbital_compton_profiles.jl","resonance_energy.jl",
        "fermi_density_effect.jl","atomic_number.jl","transport_correction.jl",
        "angular_fokker_planck_decomposition.jl","bethe.jl","moller.jl","bhabha.jl",
        "kawrakow_correction.jl","moliere_screening.jl","mott.jl","sommerfield.jl",
        "poskus.jl","ratio_positron_electron_bremsstrahlung.jl","seltzer_berger.jl",
        "klein_nishina.jl","waller_hartree.jl","impulse_approximation.jl",
        "biggs_lighthill.jl","photoelectric_per_subshell.jl","sauter.jl",
        "rayleigh.jl","high_energy_coulomb_correction.jl","baro.jl","heitler.jl",
        "relaxation.jl","interaction_interdependances.jl",
        "inelastic_collision_heavy_particle.jl","soft_catastrophic_cutoff.jl",
        "endf_reading.jl","scattering_kinematics.jl","elastic_scattering_endf.jl"
    ]
    radiant_src["particle_transport/"] = [
        "geometry.jl","volume_source.jl","surface_source.jl","scattering_source.jl",
        "particle_source.jl","fokker_planck_source.jl",
        "fokker_planck_scattering_matrix.jl","fokker_planck_finite_difference.jl",
        "fokker_planck_differential_quadrature.jl","fokker_planck_galerkin.jl",
        "fokker_planck_finite_difference_gn.jl","electromagnetic_scattering_matrix.jl",
        "transport.jl","sn_inner_pass.jl","sn_flux.jl","sn_one_speed.jl",
        "sn_sweep_1D.jl","sn_sweep_2D.jl","sn_sweep_3D.jl","sn_1D_bte.jl",
        "sn_1D_bfp.jl","sn_2D_bte.jl","sn_2D_bfp.jl","sn_3D_bte.jl",
        "sn_3D_bfp.jl","adaptive.jl","scheme_weights.jl","map_moments.jl",
        "constant_linear.jl","angular_polynomial_basis.jl","hts_source_projection.jl",
        "energy_deposition.jl","charge_deposition.jl","flux.jl",
        "restricted_to_full_domain_matrix.jl","gn_patch_geometry.jl","gn_flux.jl",
        "gn_inner_pass.jl","gn_one_speed.jl","gn_sweep_1D.jl","gn_sweep_2D.jl",
        "gn_sweep_3D.jl","gn_1D_bte.jl","gn_1D_bfp.jl","gn_2D_bte.jl",
        "gn_2D_bfp.jl","gn_3D_bte.jl","gn_3D_bfp.jl","gn_weights.jl",
        "cp_collision_probabilities.jl","cp_surface_probabilities.jl",
        "cp_sweep_1D.jl","cp_flux.jl"
    ]
    radiant_src["structures/"] = [
        "Particle.jl","Source_Normalization.jl","Boundary_Angular_Current_Source.jl",
        "Anisotropic_Volume_Source.jl","Transport_Balance.jl",
        "Transport_Ownership_Map.jl","Layer_Definition.jl","Tape_Stack_Definition.jl",
        "Energy_Partition.jl","Fibre_Diagnostic_Definition.jl",
        "HTS_Detector_Definition.jl","Physics_Coverage_Register.jl",
        "Protected_Response_Registry.jl","Interaction.jl","Elastic_Collision.jl",
        "Elastic_Scattering.jl","Inelastic_Collision.jl","Bremsstrahlung.jl",
        "Compton.jl","Pair_Production.jl","Photoelectric.jl","Annihilation.jl",
        "Rayleigh.jl","Relaxation.jl","Multigroup_Cross_Sections.jl","Material.jl",
        "Material_List.jl","Cross_Sections.jl","Geometry.jl","SN.jl","DPN.jl",
        "GN.jl","CP.jl","Solvers.jl","Surface_Source.jl","Volume_Source.jl",
        "Source.jl","Sources.jl","Fixed_Sources.jl","Flux_Per_Particle.jl",
        "Flux.jl","Electromagnetic_Field.jl","Computation_Unit.jl"
    ]
    radiant_src["tools/"] = [
        "quadrature.jl","gauss_legendre.jl","gauss_lobatto.jl",
        "gauss_legendre_chebychev.jl","lebedev.jl","carlson.jl",
        "legendre_polynomials.jl","jacobi_polynomials.jl",
        "ferrer_associated_legendre.jl","real_spherical_harmonics.jl",
        "spherical_triangle.jl","heaviside.jl","interpolation.jl","integrals.jl",
        "spline.jl","voronoi.jl","newton_bisection.jl","tanh_sinh_integral.jl",
        "one_space.jl","cache.jl","find_package_root.jl","python_method_notation.jl",
        "factorial_factor.jl","double_factorial.jl","krylov_state.jl","gmres.jl",
        "bicgstab.jl","anderson.jl","livolant.jl"
    ]
    radiant_src["interchange/"] = ["source_hdf5.jl"]
    radiant_src["hts_addon/"] = [
        "Addon_Extraction_Manifest.jl",
        "Process_Resolved_Scoring.jl",
        "Spatial_Magnetic_Field_Map.jl",
        "Spatial_Field_Transport.jl",
        "Piecewise_Flat_Tape_Atlas.jl",
        "Analytic_Curvature_Benchmark.jl",
        "Faceted_Geometry.jl",
        "Faceted_Local_Tape_Domain.jl",
        "Gd_Prompt_Capture_Cascade.jl",
        "Gd_Cascade_Data_Adapter.jl",
        "SubkeV_Thermalization.jl",
        "Statistical_Microdosimetry.jl",
        "Cryogenic_Electrothermal.jl",
        "Material_Response_Registry.jl",
        "Atomistic_Response_Table_Pipeline.jl",
        "Physical_Reference_Qualification.jl",
        "OpenSn_Radiant_Closed_Coupling.jl",
        "Completion_Extensions.jl"
    ]
    for folder in [
        "structures/","tools/","cross_sections/","particle_transport/","interchange/",
        "hts_addon/",
    ]
        for file in radiant_src[folder]
            include(string(folder,file))
        end
    end

    if isdefined(@__MODULE__, :MATERIAL_CONSTRUCTORS)
        eval(Expr(:export, getfield(@__MODULE__, :MATERIAL_CONSTRUCTORS)...))
    end

    export Particle,Photon,Electron,Positron,Proton,Antiproton,Alpha,Muon,Antimuon
    export Elastic_Collision,Elastic_Scattering,Inelastic_Collision,Bremsstrahlung,Compton
    export Pair_Production,Photoelectric,Annihilation,Rayleigh,Relaxation,Fluorescence,Auger
    export Material,Cross_Sections,Geometry,SN,Solvers,Surface_Source,Volume_Source
    export Fixed_Sources,Computation_Unit,DPN,GN,CP,Electromagnetic_Field
    export Discrete_Ordinates

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

    export BOUNDARY_SOURCE_HDF5_SCHEMA,VOLUME_SOURCE_HDF5_SCHEMA
    export write_boundary_angular_current_hdf5,read_boundary_angular_current_hdf5
    export write_anisotropic_volume_source_hdf5,read_anisotropic_volume_source_hdf5

    export add_response_channel!,set_response_channel!,has_response_channel
    export get_response_channel,get_response_channel_keys,get_response_channels
    export sum_response_channels

    export set_spatial_magnetic_field,clear_spatial_magnetic_field
    export has_spatial_magnetic_field,get_spatial_magnetic_field,get_spatial_field_hash
    export get_spatial_field_provenance,get_field_mode

    export HTS_Addon_Component,HTS_Addon_Extraction_Manifest
    export default_hts_addon_extraction_manifest,validate_addon_extraction_manifest
    export get_addon_component

    export Process_Resolved_Score,score_process_responses,get_process_score,get_process_keys
    export assert_process_score_closure,aggregate_process_scores

    export Abstract_Spatial_Magnetic_Field,Constant_Magnetic_Field
    export Cartesian_Magnetic_Field_Map,field_at,field_in_local_frame
    export electromagnetic_field_at,field_on_geometry,field_variation_bound,field_map_receipt
    export electromagnetic_field_on_geometry,spatial_field_transport_receipt

    export Tape_Atlas_Patch,Piecewise_Flat_Tape_Atlas,patch_frame
    export build_piecewise_flat_tape_atlas,global_to_patch_coordinates
    export patch_to_global_coordinates,global_direction_to_patch,patch_direction_to_global
    export electromagnetic_field_for_patch,get_atlas_patch,get_atlas_patch_at_arc_length
    export atlas_refinement_report

    export Cylindrical_Shell_Path_Result,cylindrical_shell_critical_cosine
    export cylindrical_shell_path,pure_absorption_fraction
    export Cylindrical_Shell_Heating_Result,cylindrical_shell_heating
    export Curvature_Heating_Screen,screen_cylindrical_shell_heating
    export maximum_facet_chord_for_sagitta,curvature_heating_receipt

    export FACETED_GEOMETRY_HDF5_SCHEMA,Faceted_Surface_Mesh,facet_geometry
    export faceted_topology_report,assert_transport_ready_facets
    export write_faceted_geometry_hdf5,read_faceted_geometry_hdf5
    export read_obj_facets,read_ascii_stl_facets,point_in_faceted_volume
    export Faceted_Voxelization,voxelize_faceted_mesh
    export Faceted_Tape_Patch,build_faceted_tape_patches,map_boundary_patches_to_facets
    export Facet_Boundary_Mapping_Entry,Facet_Boundary_Map,facet_boundary_map
    export Faceted_Local_Tape_Domain,Faceted_Local_Source_Receipt
    export build_faceted_local_tape_domain,local_domain_frame
    export global_to_local,local_to_global,local_layer_for_voxel
    export localize_boundary_source_to_facet,faceted_local_domain_receipt

    export Capture_Emission_Line,Gd_Prompt_Capture_Cascade,Capture_Rate_Field
    export Gd_Capture_Recoil_Source,Gd_Capture_Source_Bundle
    export build_gd_capture_source_bundle,synthetic_gd157_capture_fixture
    export Gd_Cascade_File_Manifest,Gd_Cascade_Spectrum_Summary
    export summarize_gd_gamma_cascade_file,build_gd_histogram_volume_source
    export gd_cascade_summary_is_production_ready,maurina_zenodo_gd_manifest
    export gd_cascade_data_sources

    export SubkeV_Thermalization_Kernel,SubkeV_Thermalization_Result
    export NonEquilibrium_Decay_Channel,NonEquilibrium_Thermalization_Model
    export subkev_partition_fractions,thermalize_subkev_event,get_non_equilibrium_energy
    export get_resolved_subkev_energy,is_subkev_partition_closed
    export non_equilibrium_release,non_equilibrium_power,subkev_model_is_production_ready
    export synthetic_subkev_kernel_fixture

    export Event_Energy_Partition_Fractions,Correlated_Secondary,Energy_Straggling_Model
    export Microdosimetry_Event_Prototype,Microdosimetry_Kernel
    export Weighted_Microdosimetry_Event,energy_partition_from_fractions
    export sample_microdosimetry_events,microdosimetry_effective_sample_size
    export weighted_deposited_energy_mean,weighted_deposited_energy_variance
    export detector_trigger_probability,expected_specific_energy_Gy
    export synthetic_microdosimetry_kernel_fixture

    export Tabulated_Cryogenic_Property,constant_cryogenic_property,property_value
    export Cryogenic_Thermal_Node,Cryogenic_Thermal_Link
    export HTS_Transition_Resistance_Model,effective_critical_temperature,hts_resistance
    export Electrothermal_Circuit,Cryogenic_Electrothermal_Model
    export Electrothermal_Energy_Impulse,Cryogenic_Electrothermal_Result
    export simulate_cryogenic_electrothermal,peak_active_temperature,peak_voltage
    export recovery_time_s,cryogenic_model_from_detector

    export Abstract_Material_Property_Model,NIST_Log10_Polynomial
    export NIST_Copper_Conductivity,Source_Only_Material_Property
    export Material_Property_Record,Material_Response_Registry
    export material_property_value,get_material_property,tabulate_material_property
    export cryogenic_property_from_record,default_material_response_registry

    export Atomistic_Calculation_Manifest,Atomistic_Response_Table
    export write_atomistic_response_hdf5,read_atomistic_response_hdf5
    export atomistic_table_is_production_ready,read_phono3py_kappa_hdf5
    export read_epw_spectral_function,read_electronic_stopping_table
    export default_ybco_gdbco_atomistic_plan

    export Physical_Reference_Response,Matched_Response_Comparison
    export compare_matched_response,reference_is_physical,Process_Score_Qualification
    export qualify_process_resolved_score,Curved_Reference_Qualification
    export qualify_atlas_against_curved_reference,physical_qualification_receipt
    export curved_qualification_receipt

    export Coupling_Solver_Output,Closed_Coupling_Settings,Coupling_Iteration_Record
    export Closed_Coupling_Result,scale_boundary_source,blend_boundary_sources
    export boundary_energy_current,boundary_current_relative_change,response_relative_change
    export solve_closed_opensn_radiant_coupling,closed_coupling_receipt

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
