using LALWaveforms: LALWaveforms

# The LAL-backed `LALWaveforms.hphc_fd` is the oracle. It needs the LAL shared libraries
# (`LALWAVEFORMS_LIBDIR` or a saved preference); without them the comparison is skipped, or
# fails on CI.
lal_available() = try
    LALWaveforms.hphc_fd(
        range(0, 64; step = 1.0),
        (mass_1 = 30.0, mass_2 = 30.0, luminosity_distance = 100.0, theta_jn = 0.0, phase = 0.0),
        LALWaveforms.WaveformSettings("TaylorF2"; minimum_frequency = 20.0),
    )
    true
catch err
    err isa InterruptException && rethrow()
    false
end

# Largest deviation relative to the largest LAL sample.
relative_error(h, h_lal) = maximum(abs, h .- h_lal) / maximum(abs, h_lal)

const BBH = (
    mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
)
const ALIGNED = (a_1 = 0.5, a_2 = 0.3, tilt_1 = 0.0, tilt_2 = 0.0)
const ANTIALIGNED = (a_1 = 0.5, a_2 = 0.3, tilt_1 = Float64(π), tilt_2 = 0.0)
const BNS = (
    mass_1 = 1.4, mass_2 = 1.35, luminosity_distance = 100.0, theta_jn = 0.9, phase = 0.2,
)

# (name, approximant, parameters, WaveformSettings keywords, rtol)
const CASES = [
    ("F2 BNS", "TaylorF2", BNS, (;), 1.0e-10),
    ("F2 NSBH", "TaylorF2", merge(BNS, (; mass_1 = 8.0, mass_2 = 1.4)), (;), 1.0e-10),
    ("F2 BBH", "TaylorF2", BBH, (;), 1.0e-10),
    ("F2 aligned spins", "TaylorF2", merge(BNS, (; a_1 = 0.05, a_2 = 0.02, tilt_1 = 0.0, tilt_2 = 0.0)), (;), 1.0e-10),
    ("F2 anti-aligned spins", "TaylorF2", merge(BBH, ANTIALIGNED), (;), 1.0e-10),
    ("F2 other reference frequency", "TaylorF2", merge(BNS, ALIGNED), (; reference_frequency = 25.0), 1.0e-10),
    ("F2 zero reference frequency", "TaylorF2", merge(BNS, ALIGNED), (; reference_frequency = 0.0), 1.0e-10),
    ("F2 maximum frequency", "TaylorF2", BNS, (; maximum_frequency = 400.0), 1.0e-10),
    ("F2 phase order -1", "TaylorF2", merge(BNS, ALIGNED), (; pn_phase_order = -1), 1.0e-10),
    ("F2 phase order 4", "TaylorF2", merge(BNS, ALIGNED), (; pn_phase_order = 4), 1.0e-10),
    ("F2 phase order 6", "TaylorF2", merge(BNS, ALIGNED), (; pn_phase_order = 6), 1.0e-10),
    ("F2 phase order 0", "TaylorF2", merge(BNS, ALIGNED), (; pn_phase_order = 0), 1.0e-10),
    ("F2 spin order 3", "TaylorF2", merge(BNS, ALIGNED), (; pn_spin_order = 3), 1.0e-10),
    ("F2 spin order 5", "TaylorF2", merge(BNS, ALIGNED), (; pn_spin_order = 5), 1.0e-10),
    ("F2 spin order 0", "TaylorF2", merge(BNS, ALIGNED), (; pn_spin_order = 0), 1.0e-10),
    ("F2 amplitude order 3", "TaylorF2", merge(BNS, ALIGNED), (; pn_amplitude_order = 3), 1.0e-10),
    ("F2 amplitude order -1", "TaylorF2", merge(BNS, ALIGNED), (; pn_amplitude_order = -1), 1.0e-10),
    ("F2 amplitude order 7", "TaylorF2", merge(BNS, ALIGNED), (; pn_amplitude_order = 7), 1.0e-10),
]

@testset "comparison with LAL" begin
    if !lal_available()
        msg = "LAL is not configured (set LALWAVEFORMS_LIBDIR); skipping the comparison with LALWaveforms"
        get(ENV, "CI", "") == "true" ? error(msg) : @warn msg
    else
        for grid in (range(0, 1024; step = 0.25),), (name, approximant, params, kw, rtol) in CASES
            @testset "$name" begin
                h = hphc_fd(grid, params, WaveformSettings(approximant; kw...))
                h_lal = LALWaveforms.hphc_fd(grid, params, LALWaveforms.WaveformSettings(approximant; kw...))
                err_plus = relative_error(h.plus, h_lal.plus)
                err_cross = relative_error(h.cross, h_lal.cross)
                @info "$name: relative error" plus = err_plus cross = err_cross
                @test err_plus <= rtol
                @test err_cross <= rtol
                # The support must coincide exactly.
                @test iszero.(h.plus) == iszero.(h_lal.plus)
            end
        end

        @testset "other grids" begin
            for approximant in ("TaylorF2",), df in (0.125, 1.0 / 3.0, 1.0)
                grid = range(0, 512; step = df)
                params = merge(BBH, ALIGNED)
                h = hphc_fd(grid, params, WaveformSettings(approximant))
                h_lal = LALWaveforms.hphc_fd(grid, params, LALWaveforms.WaveformSettings(approximant))
                @test relative_error(h.plus, h_lal.plus) <= 1.0e-10
            end
        end
    end
end
