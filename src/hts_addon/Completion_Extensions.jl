# This is the explicit temporary loader for HTS completion modules that will move together into
# RadiantHTS.jl. Keeping this list in one file avoids scattering nested includes through generic
# Radiant transport kernels. Include order is dependency order.

include(joinpath(@__DIR__,"CSV_Utilities.jl"))
include(joinpath(@__DIR__,"Gd_Evaluated_Cascade_Adapter.jl"))
include(joinpath(@__DIR__,"Response_Table_Completion.jl"))
include(joinpath(@__DIR__,"Curvature_Convergence_Study.jl"))
include(joinpath(@__DIR__,"Physical_Reference_Artifact.jl"))
include(joinpath(@__DIR__,"Layer_Heating_Ledger.jl"))
include(joinpath(@__DIR__,"Gd_Self_Shielding_Analytic.jl"))

export Evaluated_Gd_Cascade_Dataset,read_evaluated_gd_cascade_csv
export evaluated_gd_dataset_is_production_ready,gd_dataset_receipt
export SUBKEV_PARTITION_CHANNELS

export Layer_Heating_Contribution,Layer_Heating_Ledger,heating_ownership_key
export add_heating_contribution!,get_heating_contribution,production_heating_total
export heating_comparison_set,validate_population_energy_closure
export validate_heating_ledger,heating_ledger_receipt

export Gd_Groupwise_Capture_Component,Gd_Self_Shielding_Layer
export Gd_Self_Shielding_Cell_Result,Gd_Groupwise_Self_Shielding_Result
export macroscopic_capture_cm_inv,solve_gd_groupwise_self_shielding
export gd_groupwise_self_shielding_factor,capture_rate_field_from_gd_self_shielding
export gd_self_shielding_receipt
