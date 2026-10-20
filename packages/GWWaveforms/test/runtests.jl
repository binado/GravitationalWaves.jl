using Test
using GWWaveforms

include("../../AbstractGWWaveforms/test/contract.jl")

@testset "GWWaveforms" begin
    include("test_generators.jl")
    include("test_internals.jl")
    include("test_phenomxas.jl")
    include("test_contract.jl")
    include("test_lal_comparison.jl")
    include("test_interface.jl")
    include("test_ad.jl")
    include("test_quality.jl")
end
