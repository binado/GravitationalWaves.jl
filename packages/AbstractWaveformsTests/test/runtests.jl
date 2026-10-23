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

@testset "AbstractWaveformsTests" begin
    @testset "Aqua" begin
        Aqua.test_all(AbstractWaveformsTests)
    end

    @testset "contract on dummy generators" begin
        test_generator_contract(DummyQuadrupole(); params = (; theta_jn = 0.4), f = 0:0.5:100)
        test_generator_contract(Dummy(nothing); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
        test_generator_contract(Dummy(60.0); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
    end

    @testset "cases" begin
        cases = default_cases()
        @test haskey(cases, "bbh_precessing")
        @test_throws ArgumentError fd_parameters(cases["bbh_precessing"])
        ids(g) = first.(applicable_cases(g))
        @test !("bbh_precessing" in ids(Shifted(0.0, 0.0)))
        @test "bbh_aligned" in ids(Shifted(0.0, 0.0))
        @test first.(applicable_cases(Shifted(0.0, 0.0); ids = ("bbh_unequal", "bbh_precessing", "bbh_aligned"))) ==
            ["bbh_aligned", "bbh_unequal"]
        @test_throws ArgumentError applicable_cases(Shifted(0.0, 0.0); ids = ("nonexistent",))
        p = cases["bbh_unequal"]
        x, scales = scaled_coordinates(p)
        @test fd_adapter(p, scales)(x).mass_1 ≈ p.mass_1_msun
        @test fd_adapter(p, scales)(x).luminosity_distance ≈ p.distance_mpc
        @test surrogate_adapter(first(scaled_coordinates(p; precessing = true)) .* 0 .+ 1)(ones(11)).total_mass > 0
        @test_throws ErrorException load_cases(joinpath(@__DIR__, "runtests.jl"))
    end

    @testset "mismatch" begin
        metric = flat_metric(F; minimum_frequency = 20.0)
        h(τ, φ) = hcat([shifted(f, τ, φ) for f in F], [0.5 * shifted(f, τ, φ) for f in F])
        a = h(0.0, 0.0)
        @test mismatch(a, a, metric; maximize = :none) < 1.0e-12
        @test mismatch(a, 3 .* a, metric; maximize = :none) < 1.0e-12
        @test mismatch(a, h(0.0, 0.7), metric; maximize = :none) > 0.1
        @test mismatch(a, h(0.0, 0.7), metric; maximize = :phase) < 1.0e-12
        @test mismatch(a, h(0.003, 0.0), metric; maximize = :phase) > 1.0e-2
        @test mismatch(a, h(0.003, 0.7), metric; maximize = :time_phase) < 1.0e-9
        @test mismatch(a, h(-0.0123, -2.0), metric; maximize = :time_phase) < 1.0e-9
        # A different shape is not absorbed by any shift.
        other = hcat([exp(-((f - 150.0) / 30.0)^2 / 2) + 0im for f in F], zeros(ComplexF64, length(F)))
        @test mismatch(a, other, metric) > 1.0e-2
        @test_throws DimensionMismatch mismatch(a, a[1:(end - 1), :], metric)
        @test_throws ArgumentError mismatch(a, a, metric; maximize = :nonsense)
        @test_throws DomainError mismatch(a, zero(a), metric)
        series = FrequencySeries(a; epoch = 0.0, frequency_spacing = step(F), channels = (:plus, :cross))
        @test mismatch(series, series, metric) < 1.0e-12
    end

    @testset "test_reference" begin
        case = [("shift", (;))]
        test_reference(Shifted(0.002, 0.3), Shifted(0.0, 0.0), case; grid = F, parameters = identity, max_mismatch = 1.0e-9)
        test_reference(
            Shifted(0.0, 0.0), Shifted(0.0, 0.0), case; grid = F, parameters = identity,
            max_mismatch = 1.0e-12, rtol = 1.0e-12, exact_support = true,
        )
        @test_throws ArgumentError test_reference(Shifted(0.0, 0.0), Shifted(0.0, 0.0), []; grid = F, max_mismatch = 1.0)
        # Extension point: a user-supplied metric function receives the frequency grid.
        test_reference(
            Shifted(0.001, 0.0), Shifted(0.0, 0.0), case; grid = F, parameters = identity,
            metric = f -> flat_metric(f; minimum_frequency = 50.0, maximum_frequency = 200.0), max_mismatch = 1.0e-9,
        )
        # Time domain: the same pulse shifted by whole samples agrees after the physical rfft.
        t = range(-0.5, 0.5; step = 1 / 2048)
        test_reference(
            Pulse(10 / 2048, 0.02), Pulse(0.0, 0.02), case; grid = t, parameters = identity,
            target_domain = :td, max_mismatch = 1.0e-8,
        )
        @testset "different pulses are told apart" begin
            evaluate = AbstractWaveformsTests.evaluate_polarizations
            wide, f = evaluate(Pulse(0.0, 0.03), (;), t, :td, identity)
            narrow, _ = evaluate(Pulse(0.0, 0.02), (;), t, :td, identity)
            @test mismatch(wide, narrow, flat_metric(f)) > 1.0e-3
        end
    end

    @testset "gradients" begin
        backends = (AutoForwardDiff(),)
        x = [1.3, 0.4, 1.0, 0.2, 0.5]
        objective(z) = z[1]^2 * sin(z[2]) / z[3] + z[4] * z[5]^2
        test_gradients(objective, x; backends, rtol = 1.0e-8)
        test_gradients(objective, x; backends, rtol = 1.0e-8, reuse_index = 3)

        # The last three coordinates are distance (scaled to one), inclination and phase.
        fgrid = 0.0:1.0:200.0
        strain(z) = reshape(
            [(z[1] / z[3]) * (1 + 0.3 * z[4]) * exp(-((f - 100) / 30)^2 / 2) * cis(z[2] * f / 100 + z[5]) for f in fgrid], :, 1,
        )
        metric = flat_metric(fgrid; channels = (:strain,), minimum_frequency = 20.0)
        y = [2.0, 0.1, 1.0, 0.4, 0.3]
        test_snr_gradients(strain, metric, y; backends, rtol = 1.0e-8)
        test_phase_gradients(strain, y; backends, rtol = 1.0e-8)
    end

    @testset "waveform Jacobian" begin
        struct Parametrized <: AbstractWaveformGenerator end
        AbstractWaveforms.approximant(::Parametrized) = "Parametrized"
        AbstractWaveforms.reference_frequency(::Parametrized) = 50.0
        AbstractWaveforms.minimum_frequency(::Parametrized) = 20.0
        AbstractWaveforms.maximum_frequency(::Parametrized) = nothing
        function AbstractWaveforms.hphc_fd(f::AbstractRange{Float64}, p::NamedTuple, g::Parametrized)
            v = p.amplitude .* cis.(p.slope .* f)
            h = ifelse.(f .>= minimum_frequency(g), v, zero(eltype(v)))
            return hcat(h, 0.5 .* h)
        end
        grid = 0.0:2.0:200.0
        test_waveform_jacobian(
            Parametrized(), [1.5, 0.02], z -> (; amplitude = z[1], slope = z[2]);
            grid, backends = (AutoForwardDiff(),), rtol = 1.0e-8,
        )
    end
end
