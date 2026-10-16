using Test
using GWWaveforms

@testset "GWWaveforms" begin
    include("test_settings.jl")
    include("test_internals.jl")
    include("test_lal_comparison.jl")
    include("test_ad.jl")
    include("test_quality.jl")
end
