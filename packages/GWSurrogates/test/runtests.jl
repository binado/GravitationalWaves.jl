using Test
using GWSurrogates

@testset "GWSurrogates" begin
    include("test_spline.jl")
    include("test_wignerd.jl")
    include("test_norms.jl")
    include("test_coorbital.jl")
    include("test_dynamics.jl")
    include("test_quality.jl")
end
