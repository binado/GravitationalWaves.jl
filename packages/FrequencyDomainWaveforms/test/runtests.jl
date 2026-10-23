using GWSeries
using Test
using FrequencyDomainWaveforms

using AbstractWaveformsTests

@testset "FrequencyDomainWaveforms" begin
    include("test_generators.jl")
    include("test_internals.jl")
    include("test_phenomxas.jl")
    include("test_contract.jl")
    include("test_amplitude_phase.jl")
    include("test_interface.jl")
    include("test_ad.jl")
    include("test_quality.jl")
end
