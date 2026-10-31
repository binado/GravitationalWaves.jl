using TestItems

@testitem "Spectral density tables" setup = [Common] begin
    curve = TabulatedPSD([0 2; 2 6; 4 10])
    @test psd.([0, 1, 4], curve) == [2, 4, 10]
    @test psd.(0:2:4, curve) == [2, 6, 10]
    @test psd.(1, [curve, curve]) == [4, 4]
    @test asd(1, curve) == 2
    @test psd.([-1, 5], curve; outside = :zero) == [0, 0]
    @test_throws DomainError psd(5, curve)
    @test_throws ArgumentError TabulatedPSD([1, 1], [2, 3])
    @test_throws ArgumentError TabulatedPSD([1, 2], [2, -3])
    asd_curve = TabulatedPSD(joinpath(@__DIR__, "data", "asd.csv"), ','; kind = :asd)
    @test psd.([1, 2, 3], asd_curve) == [4, 10, 16]
    @test asd.([1, 3], asd_curve) == [2, 4]
    table = [1 2; 3 4]
    copied = TabulatedPSD(view(table, :, :))
    table .= 0
    @test copied.frequencies == [1, 3]
    @test copied.values == [2, 4]
    @test_throws ArgumentError TabulatedPSD(zeros(2, 1))
    @test_throws ArgumentError TabulatedPSD(zeros(2, 3))
    @test_throws ArgumentError TabulatedPSD(zeros(0, 2))
    @test_throws ArgumentError TabulatedPSD([1 2])
    powerlaw = TabulatedPSD([1, 10, 100], [1, 100, 10000]; interpolation = :loglog)
    @test psd.([2, 20, 50], powerlaw) ≈ [4, 400, 2500]
    @test_throws ArgumentError TabulatedPSD([0, 1], [1, 2]; interpolation = :loglog)
end

@testitem "Reading files" setup = [Common] begin
    path = joinpath(@__DIR__, "data", "psd.csv")
    curve = TabulatedPSD(path, ',')
    @test curve.frequencies == [0, 2, 4]
    @test curve.values == [2, 6, 10]
    open(path) do io
        @test TabulatedPSD(io, ',').values == curve.values
    end
    loglog = TabulatedPSD(IOBuffer("1,1\n10,100\n100,10000\n"), ','; interpolation = :loglog)
    @test psd(2, loglog) ≈ 4
    tabs = TabulatedPSD(IOBuffer("1\t2\n3\t4\n"), '\t')
    @test psd.([1, 3], tabs) == [2, 4]
    spaces = TabulatedPSD(IOBuffer("1 2\n3 4\n"), ' ')
    @test spaces.values == [2, 4]
    @test_throws ErrorException TabulatedPSD(path, ';')
    @test_throws ErrorException TabulatedPSD(IOBuffer("1  2\n3  4\n"), ' ')
    @test_throws ArgumentError TabulatedPSD(IOBuffer("1,2,3\n4,5,6\n"), ',')
    @test_throws ErrorException TabulatedPSD(IOBuffer(""), ',')
    @test_throws ErrorException TabulatedPSD(IOBuffer("# only a comment\n\n"), ',')
    @test_throws ErrorException TabulatedPSD(IOBuffer("1,2\n3,abc\n"), ',')
    @test_throws ArgumentError TabulatedPSD(IOBuffer("1,2\n3,NaN\n"), ',')
    @test_throws ArgumentError TabulatedPSD(IOBuffer("1,2\n3,-4\n"), ',')
    @test_throws ArgumentError TabulatedPSD(IOBuffer("1,Inf\n3,4\n"), ','; kind = :asd)
end
