using Test
using LALWaveforms

@testset "LALWaveforms" begin
    include("test_library.jl")
    include("test_conversion.jl")
    include("test_fd.jl")
    include("test_reference.jl")
    include("test_quality.jl")
end
