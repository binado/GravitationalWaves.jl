using Test
using Random
using GWNoise

@testset "GWNoise" begin
    @testset "Spectral density tables" begin
        curve = read_psd(IOBuffer("# frequency PSD\n\n0 2\n2 6 # comment\n4 10\n"))
        @test psd.([0, 1, 4], curve) == [2, 4, 10]
        @test psd.(0:2:4, curve) == [2, 6, 10]
        @test psd.(1, [curve, curve]) == [4, 4]
        @test asd(1, curve) == 2
        @test psd.([-1, 5], curve; outside = :zero) == [0, 0]
        @test_throws DomainError psd(5, curve)
        @test_throws ArgumentError TabulatedPSD([1, 1], [2, 3])
        @test_throws ArgumentError TabulatedPSD([1, 2], [2, -3])
        mktemp() do path, io
            write(io, "# ASD\n1,2\n3,4\n")
            flush(io)
            curve = read_asd(path)
            @test psd.([1, 2, 3], curve) == [4, 10, 16]
            @test asd.([1, 3], curve) == [2, 4]
        end
        powerlaw = TabulatedPSD([1, 10, 100], [1, 100, 10000]; interpolation = :loglog)
        @test psd.([2, 20, 50], powerlaw) ≈ [4, 400, 2500]
        @test_throws ArgumentError TabulatedPSD([0, 1], [1, 2]; interpolation = :loglog)
    end

    @testset "Reading files" begin
        loglog = read_psd(IOBuffer("1 1\n10 100\n100 10000\n"); interpolation = :loglog)
        @test psd(2, loglog) ≈ 4
        csv = read_psd(IOBuffer("# frequency, PSD\n1, 2\n3 ,4\n"))
        @test psd.([1, 3], csv) == [2, 4]
        tabs = read_psd(IOBuffer("1\t2\n3\t4\n"))
        @test psd.([1, 3], tabs) == [2, 4]
        @test_throws ArgumentError read_psd(IOBuffer("1 2 3\n4 5 6\n"))
        @test_throws ArgumentError read_psd(IOBuffer(""))
        @test_throws ArgumentError read_psd(IOBuffer("# only a comment\n\n"))
        @test_throws ErrorException read_psd(IOBuffer("1 2\n3 abc\n"))
        @test_throws ArgumentError read_psd(IOBuffer("1 2\n3 NaN\n"))
        @test_throws ArgumentError read_psd(IOBuffer("1 2\n3 -4\n"))
        @test_throws ArgumentError read_asd(IOBuffer("1 Inf\n3 4\n"))
    end

    include("test_bundled.jl")
    include("test_arrays.jl")
    include("test_quality.jl")
end
