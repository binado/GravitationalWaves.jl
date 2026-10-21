# Run in a temporary environment containing the generated JLL and LALWaveforms.
# No Python, preferences, or downloads of waveform datasets are needed.
using Test, Libdl
using LALWaveforms

function lal_libraries()
    return filter(dllist()) do path
        occursin(r"^liblal(support|simulation)?([.-].*)?\.(dylib|so)(\.\d+)*$", basename(path))
    end
end

@testset "native artifact smoke" begin
    @test LALWaveforms.STATE[] === nothing
    @test isempty(lal_libraries())
    # Deliberately load the candidate JLL after checking package import laziness.
    jll = Base.require(Main, :LALSuite_jll)
    @test Base.invokelatest(jll.is_available)
    @test isempty(lal_libraries())
    @test isfile(jll.liblal_path)
    @test isfile(jll.liblalsupport_path)
    @test isfile(jll.liblalsimulation_path)
    @test dirname(jll.liblal_path) == dirname(jll.liblalsimulation_path)

    # Until registration, exercise these artifacts through the existing explicit
    # provider. The adoption PR must additionally test default provider selection.
    withenv("LALWAVEFORMS_LIBDIR" => dirname(jll.liblal_path)) do
        fns = LALWaveforms.lal()
        @test LALWaveforms.lal() === fns
        for symbol in (LALWaveforms.LIBLAL_SYMBOLS..., LALWaveforms.LIBSIM_SYMBOLS...)
            @test getproperty(fns, symbol) != C_NULL
        end
        f = range(0, 1024; step = 0.25)
        params = (
            mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0,
            theta_jn = 0.4, phase = 1.3,
        )
        precessing = merge(
            params, (
                a_1 = 0.5, a_2 = 0.3, tilt_1 = 0.6, tilt_2 = 1.0,
                phi_12 = 0.7, phi_jl = 0.2,
            )
        )
        for (approximant, p) in (
                ("TaylorF2", params), ("IMRPhenomD", params),
                ("IMRPhenomXPHM", precessing), ("TaylorT4", params),
            )
            @testset "$approximant" begin
                generator = LALWaveformGenerator(
                    approximant; minimum_frequency = 20.0, maximum_frequency = 1024.0
                )
                h = hphc_fd(f, p, generator)
                @test length(h.plus) == length(h.cross) == length(f)
                @test all(isfinite, h.plus) && all(isfinite, h.cross)
                @test any(!iszero, h.plus) && any(!iszero, h.cross)
                @test all(iszero, h.plus[f .< 20.0])
            end
        end
    end
    loaded = realpath.(lal_libraries())
    @test realpath(jll.liblal_path) in loaded
    @test realpath(jll.liblalsimulation_path) in loaded
    @test !any(path -> occursin("libpython", lowercase(basename(path))), dllist())
    @info "Native smoke passed" platform = Sys.MACHINE julia = VERSION libraries = loaded
end
