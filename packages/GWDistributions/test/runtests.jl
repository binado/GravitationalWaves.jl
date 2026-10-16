using Test
using GWDistributions

@testset "GWDistributions" begin
    include("test_utils.jl")
    include("test_interpolated.jl")
    include("test_mass.jl")
    include("test_delay.jl")
    include("test_orientation.jl")
    include("test_redshift.jl")
    include("test_time_delayed.jl")
    include("test_astrogwb_reference.jl")
    include("test_quality.jl")
end
