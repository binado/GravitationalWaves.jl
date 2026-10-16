using Test
using GWDistributions

@testset "GWDistributions" begin
    include("test_utils.jl")
    include("test_mass.jl")
    include("test_redshift.jl")
    include("test_quality.jl")
end
