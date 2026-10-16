@testset "library" begin
    LW = LALWaveforms

    @testset "discovery" begin
        # macOS wheel: lal/.dylibs
        mktempdir() do dir
            for n in ("liblal.20.dylib", "liblalsimulation.37.dylib", "liblalburst.8.dylib", "libgsl.28.dylib")
                touch(joinpath(dir, n))
            end
            lal, sim = LW.find_libraries(dir)
            @test basename(lal) == "liblal.20.dylib"
            @test basename(sim) == "liblalsimulation.37.dylib"
        end
        # Linux auditwheel: lalsuite.libs, hash-suffixed names
        mktempdir() do dir
            for n in (
                    "liblal-9eb973fe.so.20.5.1", "liblalsimulation-23824854.so.37.2.1",
                    "liblalburst-6637ed4a.so.8.0.1", "liblalframe-d713deed.so.14.0.4",
                )
                touch(joinpath(dir, n))
            end
            lal, sim = LW.find_libraries(dir)
            @test basename(lal) == "liblal-9eb973fe.so.20.5.1"
            @test basename(sim) == "liblalsimulation-23824854.so.37.2.1"
        end
        # Plain unix layout
        mktempdir() do dir
            touch(joinpath(dir, "liblal.so"))
            touch(joinpath(dir, "liblalsimulation.so.37"))
            @test LW.find_libraries(dir) !== nothing
        end
        # Incomplete / missing directories
        mktempdir() do dir
            touch(joinpath(dir, "liblal.20.dylib"))
            @test LW.find_libraries(dir) === nothing
            @test LW.find_libraries(joinpath(dir, "nope")) === nothing
        end
        @test_throws ArgumentError set_libdir!(tempdir())
    end

    @testset "load" begin
        fns = LW.lal()
        approx = ccall(
            fns.XLALSimInspiralGetApproximantFromString, Cint, (Cstring,), "IMRPhenomD"
        )
        @test approx >= 0
        @test ccall(
            fns.XLALSimInspiralImplementedFDApproximants, Cint, (Cint,), approx
        ) != 0
    end
end
