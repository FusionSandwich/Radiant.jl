using Radiant
using Test
import Radiant

@testset "Radiant.jl" begin
    
    include("alpha_stopping.jl")
    include("alpha_collision_moments.jl")
    include("alpha_elastic.jl")
    include("csd_kernel_normalization.jl")
    include("gn_csd_kernel_normalization.jl")
    include("geometry_origins.jl")

end
nothing


