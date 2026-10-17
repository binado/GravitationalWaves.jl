using Test
using LALWaveforms

include("../../AbstractGWWaveforms/test/contract.jl")

@testset "LALWaveforms" begin
    include("test_library.jl")
    include("test_conversion.jl")
    include("test_fd.jl")
    include("test_reference.jl")
    include("test_contract.jl")
    include("test_quality.jl")
end
