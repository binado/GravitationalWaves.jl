using Test
using GWSeries
using GWDetectors
using ForwardDiff
using LinearAlgebra
using StaticArrays

@testset "GWDetectors" begin
    include("test_geometry.jl")
    include("test_sensitivity.jl")
    include("test_response.jl")
    include("test_projection.jl")
    include("test_reference.jl")
    include("test_ad.jl")
    include("test_domains.jl")
    include("test_quality.jl")
end
