using Test
using GWSeries
using AbstractFFTs
using FFTW
using LinearAlgebra
using Random

include("test_series.jl")
include("test_fourier.jl")
include("test_windows.jl")
using ArrayInterface
@testset "ArrayInterface storage traits" begin
    x = TimeSeries(ones(8); epoch = 0, sampling_frequency = 8, channels = (:H1,))
    @test ArrayInterface.parent_type(typeof(x)) === typeof(parent(x))
    @test ArrayInterface.can_setindex(x)
    @test ArrayInterface.fast_scalar_indexing(x)
    @test !ArrayInterface.defines_strides(x)
end
include("test_quality.jl")
