using Test
using GWDetectors
using ForwardDiff

@testset "GWDetectors" begin
    include("test_geometry.jl")
    include("test_response.jl")
    include("test_projection.jl")
    include("test_reference.jl")
    include("test_ad.jl")
    include("test_quality.jl")
end
