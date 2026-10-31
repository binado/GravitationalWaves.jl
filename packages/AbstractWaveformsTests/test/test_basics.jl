using TestItems

@testitem "Aqua" setup = [Common] begin
    Aqua.test_all(AbstractWaveformsTests)
end

@testitem "contract on dummy generators" setup = [Common] begin
    test_generator_contract(DummyQuadrupole(); params = (; theta_jn = 0.4), f = 0:0.5:100)
    test_generator_contract(Dummy(nothing); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
    test_generator_contract(Dummy(60.0); params = (mass_1 = 30.0, mass_2 = 30.0), f = 0:0.5:100)
end

@testitem "cases" setup = [Common] begin
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

@testitem "mismatch" setup = [Common] begin
    W = flat_whitening(F; minimum_frequency = 20.0)
    h(τ, φ) = hcat([shifted(f, τ, φ) for f in F], [0.5 * shifted(f, τ, φ) for f in F])
    a = h(0.0, 0.0)
    @test mismatch(a, a, W, F; maximize = :none) < 1.0e-12
    @test mismatch(a, 3 .* a, W, F; maximize = :none) < 1.0e-12
    @test mismatch(a, h(0.0, 0.7), W, F; maximize = :none) > 0.1
    @test mismatch(a, h(0.0, 0.7), W, F; maximize = :phase) < 1.0e-12
    @test mismatch(a, h(0.003, 0.0), W, F; maximize = :phase) > 1.0e-2
    @test mismatch(a, h(0.003, 0.7), W, F; maximize = :time_phase) < 1.0e-9
    @test mismatch(a, h(-0.0123, -2.0), W, F; maximize = :time_phase) < 1.0e-9
    # A different shape is not absorbed by any shift.
    other = hcat([exp(-((f - 150.0) / 30.0)^2 / 2) + 0im for f in F], zeros(ComplexF64, length(F)))
    @test mismatch(a, other, W, F) > 1.0e-2
    @test_throws DimensionMismatch mismatch(a, a[1:(end - 1), :], W, F)
    @test_throws ArgumentError mismatch(a, a, W, F; maximize = :nonsense)
    @test_throws DomainError mismatch(a, zero(a), W, F)
    series = FrequencySeries(a; epoch = 0.0, frequency_spacing = step(F), channels = (:plus, :cross))
    @test mismatch(series, series, W) < 1.0e-12
end

@testitem "test_reference" setup = [Common] begin
    case = [("shift", (;))]
    test_reference(Shifted(0.002, 0.3), Shifted(0.0, 0.0), case; grid = F, parameters = identity, max_mismatch = 1.0e-9)
    test_reference(
        Shifted(0.0, 0.0), Shifted(0.0, 0.0), case; grid = F, parameters = identity,
        max_mismatch = 1.0e-12, rtol = 1.0e-12, exact_support = true,
    )
    @test_throws ArgumentError test_reference(Shifted(0.0, 0.0), Shifted(0.0, 0.0), []; grid = F, max_mismatch = 1.0)
    # Extension point: a user-supplied whitening function receives the frequency grid.
    test_reference(
        Shifted(0.001, 0.0), Shifted(0.0, 0.0), case; grid = F, parameters = identity,
        whitening = f -> flat_whitening(f; minimum_frequency = 50.0, maximum_frequency = 200.0), max_mismatch = 1.0e-9,
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
        @test mismatch(wide, narrow, flat_whitening(f), f) > 1.0e-3
    end
end

@testitem "gradients" setup = [Common] begin
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
    W = flat_whitening(fgrid; nchannels = 1, minimum_frequency = 20.0)
    y = [2.0, 0.1, 1.0, 0.4, 0.3]
    test_snr_gradients(strain, W, step(fgrid), y; backends, rtol = 1.0e-8)
    test_phase_gradients(strain, y; backends, rtol = 1.0e-8)
end

@testitem "waveform Jacobian" setup = [Common] begin
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
