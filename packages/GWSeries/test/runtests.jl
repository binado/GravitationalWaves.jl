using Test
using GWSeries
using DimensionalData
using AbstractFFTs
using FFTW
using LinearAlgebra
using Random

include("test_series.jl")
include("test_fourier.jl")
include("test_windows.jl")
# Storage traits and static arrays come from DimensionalData's own extensions.
using ArrayInterface
@testset "ArrayInterface storage traits" begin
    x = TimeSeries(ones(8); epoch = 0, sampling_frequency = 8, channels = (:H1,))
    @test ArrayInterface.parent_type(typeof(x)) === typeof(parent(x))
    @test ArrayInterface.can_setindex(x)
    @test ArrayInterface.fast_scalar_indexing(x)
    @test ArrayInterface.defines_strides(typeof(x)) == ArrayInterface.defines_strides(typeof(parent(x)))
end
using StaticArrays
@testset "Static-array interoperability" begin
    s = TimeSeries(ones(2, 2); epoch = 0, sampling_interval = 1, channels = (:H1, :L1))
    a = @SMatrix [1.0 0.0; 0.0 1.0]
    # DimensionalData's StaticArrays extension decides these: `+` and `-` return unlabelled
    # static arrays, and broadcasting a static matrix with any DimArray is unsupported upstream.
    for result in (a + s, s + a, a - s, s - a)
        @test result isa StaticArray
    end
    @test a + s == a + parent(s)
    @test_broken (a .+ s) isa TimeSeries
    out = copy(s)
    out .= 2 .* a
    @test parent(out) == 2 .* a
end
include("test_interface.jl")
include("test_quality.jl")
