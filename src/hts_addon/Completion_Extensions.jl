# This is the explicit temporary loader for HTS completion modules that will move together into
# RadiantHTS.jl. Keeping this list in one file avoids scattering nested includes through generic
# Radiant transport kernels. Include order is dependency order.

include(joinpath(@__DIR__,"CSV_Utilities.jl"))
include(joinpath(@__DIR__,"Gd_Evaluated_Cascade_Adapter.jl"))
include(joinpath(@__DIR__,"Response_Table_Completion.jl"))
include(joinpath(@__DIR__,"Curvature_Convergence_Study.jl"))
include(joinpath(@__DIR__,"Physical_Reference_Artifact.jl"))

export Evaluated_Gd_Cascade_Dataset,read_evaluated_gd_cascade_csv
export evaluated_gd_dataset_is_production_ready,gd_dataset_receipt
export SUBKEV_PARTITION_CHANNELS
