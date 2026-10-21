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
using StaticArrays
@testset "Static-array interoperability" begin
    s = TimeSeries(ones(2, 2); epoch = 0, sampling_interval = 1, channels = (:H1, :L1))
    a = @SMatrix [1.0 0.0; 0.0 1.0]
    for result in (a + s, s + a, a - s, s - a, a * s, s * a, map(+, a, s), map(+, s, a), a .+ s, s .+ a)
        @test result isa TimeSeries
        @test metadata(result) == metadata(s)
    end
    out = copy(s)
    out .= 2 .* a
    @test parent(out) == 2 .* a
    @test isempty(detect_ambiguities(GWSeries, Base.get_extension(GWSeries, :GWSeriesStaticArraysCoreExt); recursive = true))
end
include("test_quality.jl")
