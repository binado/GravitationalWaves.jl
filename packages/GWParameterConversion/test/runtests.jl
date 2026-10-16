using Test
using GWParameterConversion

@testset "GWParameterConversion" begin
    include("test_masses.jl")
    include("test_reference.jl")
    include("test_quality.jl")
end
