using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using Test
    using Aqua
    using FFTW
    using ForwardDiff
    using ADTypes: AutoForwardDiff
    using AbstractWaveforms
    using AbstractWaveforms: check_frequency_grid, frequency_bounds, AlignedSpinQuadrupolarWaveformGenerator
    using AbstractWaveformsTests
    using GWSeries

    # Smooth band-limited spectrum and its exact time/phase-shifted copies.
    spectrum(f) = exp(-((f - 100.0) / 30.0)^2 / 2) * cis(0.01 * f)
    shifted(f, τ, φ) = spectrum(f) * cis(-2π * f * τ + φ)

    struct Shifted <: AbstractWaveformGenerator
        time::Float64
        phase::Float64
    end
    AbstractWaveforms.approximant(g::Shifted) = "Shifted"
    AbstractWaveforms.reference_frequency(::Shifted) = 50.0
    AbstractWaveforms.minimum_frequency(::Shifted) = 20.0
    AbstractWaveforms.maximum_frequency(::Shifted) = nothing
    function AbstractWaveforms.hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, g::Shifted)
        check_frequency_grid(f)
        lo, hi = frequency_bounds(g, f)
        h = [lo <= fi <= hi ? shifted(fi, g.time, g.phase) : 0.0im for fi in f]
        return hcat(h, 0.5 .* h)
    end
    AbstractWaveforms.validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, ::Shifted) = nothing

    struct Pulse <: AbstractWaveformGenerator
        shift::Float64
        width::Float64
    end
    AbstractWaveforms.approximant(::Pulse) = "Pulse"
    AbstractWaveforms.reference_frequency(::Pulse) = 0.0
    AbstractWaveforms.minimum_frequency(::Pulse) = 0.0
    AbstractWaveforms.maximum_frequency(::Pulse) = nothing
    function AbstractWaveforms.hphc_td(prototype::TimeSeries, params::NamedTuple, g::Pulse)
        w = [exp(-((t - g.shift) / g.width)^2) * cos(2π * 100 * (t - g.shift)) for t in times(prototype)]
        return TimeSeries(
            hcat(w, 0.5 .* w); epoch = epoch(prototype), sampling_interval = sampling_interval(prototype),
            channels = (:plus, :cross),
        )
    end

    struct Dummy <: AbstractWaveformGenerator
        maximum::Union{Nothing, Float64}
    end
    AbstractWaveforms.approximant(::Dummy) = "Dummy"
    AbstractWaveforms.reference_frequency(::Dummy) = 50.0
    AbstractWaveforms.minimum_frequency(::Dummy) = 20.0
    AbstractWaveforms.maximum_frequency(g::Dummy) = g.maximum
    function AbstractWaveforms.hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, g::Dummy)
        check_frequency_grid(f)
        f_min, f_max = frequency_bounds(g, f)
        h = [f_min <= fi <= f_max ? complex(1.0, fi) : 0.0im for fi in f]
        return hcat(h, h)
    end
    function AbstractWaveforms.validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, ::Dummy)
        foreach(k -> k in (:mass_1, :mass_2) || throw(ArgumentError("unknown parameter $k")), keys(params))
        return nothing
    end

    # Implements only amplitude_phase_fd; polarization assembly comes from the interface.
    struct DummyQuadrupole <: AlignedSpinQuadrupolarWaveformGenerator end
    AbstractWaveforms.approximant(::DummyQuadrupole) = "DummyQuadrupole"
    AbstractWaveforms.reference_frequency(::DummyQuadrupole) = 50.0
    AbstractWaveforms.minimum_frequency(::DummyQuadrupole) = 20.0
    AbstractWaveforms.maximum_frequency(::DummyQuadrupole) = 60.0
    function AbstractWaveforms.amplitude_phase_fd(f::AbstractRange{Float64}, params::NamedTuple, g::DummyQuadrupole)
        check_frequency_grid(f)
        lo, hi = frequency_bounds(g, f)
        amplitude = [lo <= fi <= hi ? 2.0 : 0.0 for fi in f]
        phase = [lo <= fi <= hi ? 0.1 * fi - 2get(params, :phase, 0.0) : 0.0 for fi in f]
        return (; amplitude, phase)
    end
    function AbstractWaveforms.validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, ::DummyQuadrupole)
        check_frequency_grid(f)
        foreach(k -> k == :theta_jn || throw(ArgumentError("unknown parameter $k")), keys(params))
        return nothing
    end

    const F = 0.0:0.5:1024.0
end
