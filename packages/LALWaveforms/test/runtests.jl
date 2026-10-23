using GWSeries
using Test
using LALWaveforms

using AbstractWaveformsTests

@testset "LALWaveforms" begin
    include("test_library.jl")
    include("test_conversion.jl")
    include("test_fd.jl")
    include("test_reference.jl")
    include("test_contract.jl")
    include("test_quality.jl")
end
