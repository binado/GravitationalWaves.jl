using TestItems

# Evaluated afresh in each test item, including mutable fixtures and dummy types.
@testsnippet Common begin
    using GWSeries
    using Test
    using Aqua
    using JET
    using AbstractWaveforms
    using AbstractWaveforms: check_frequency_grid, frequency_bounds
    using AbstractWaveforms: AlignedSpinQuadrupolarWaveformGenerator, amplitude_phase_fd
    using AbstractWaveformsTests

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

    # Records the coordinates a kernel receives, to check that series grids are unwrapped.
    struct GridProbe <: AbstractWaveformGenerator
        seen::Base.RefValue{Any}
    end
    function AbstractWaveforms.hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, g::GridProbe)
        g.seen[] = f
        return zeros(ComplexF64, length(f), 2)
    end
    function AbstractWaveforms.hphc_td(t::AbstractVector{<:Real}, params::NamedTuple, g::GridProbe)
        g.seen[] = t
        return zeros(length(t), 2)
    end
end
