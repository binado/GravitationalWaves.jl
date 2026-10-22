using Test, GWNoise, GWMatchedFiltering, DifferentiationInterface, ForwardDiff, Mooncake
using FiniteDifferences: central_fdm, grad

# Captures contain only fixed configuration, never aliases of differentiated x.
# The prepared metric is fixed configuration: declare it non-differentiable, since
# Mooncake cannot copy the TwicePrecision frequency metadata stored inside it.
struct FixedObjective{F}
    f::F
end
(f::FixedObjective)(x) = f.f(x)
Mooncake.tangent_type(::Type{<:FixedObjective}) = Mooncake.NoTangent
Mooncake.tangent_type(::Type{<:FrequencyDomainMetric}) = Mooncake.NoTangent

const SNR_BACKENDS = (AutoForwardDiff(), AutoMooncake(; config = nothing))
const PSD_NAMES = (:ET_000A_18_ETD_hlf, :cosmic_explorer_strain)
# A fixed response with both polarizations observable; geometry is not differentiated.
project_strain(h) = 0.7 .* h[:, 1] .+ 0.3 .* h[:, 2]

function snr_metric(f, name; nsamples = nothing)
    curve = read_psd(psd_path(name); interpolation = :loglog)
    low = max(20.0, first(curve.frequencies))
    high = min(512.0, last(curve.frequencies), last(f))
    spectrum = [low <= x <= high ? psd(x, curve) : 0.0 for x in f]
    return FrequencyDomainMetric(
        spectrum; frequencies = f, channels = (:strain,),
        minimum_frequency = low, maximum_frequency = high, nsamples
    )
end

function check_snr_ad(objective, x; rtol, finite_difference = false)
    loss = FixedObjective(objective)
    reference = loss(x)
    @test isfinite(reference) && reference > 0
    gradients = map(SNR_BACKENDS) do backend
        prep = prepare_gradient(loss, backend, x)
        value, g = value_and_gradient(loss, prep, backend, x)
        @test value ≈ reference rtol = 1.0e-10
        @test all(isfinite, g)
        # Reuse at a nearby distance (last three coordinates are distance/iota/phase).
        nearby = copy(x)
        nearby[end - 2] *= 1.001
        reused_value, reused = value_and_gradient(loss, prep, backend, nearby)
        expected_value, expected = value_and_gradient(loss, backend, nearby)
        @test reused_value ≈ expected_value rtol = 1.0e-10
        @test all(isfinite, reused)
        @test all(isapprox.(reused ./ reference, expected ./ reference; rtol, atol = 1.0e-8))
        g
    end
    @test all(isapprox.(gradients[1] ./ reference, gradients[2] ./ reference; rtol, atol = 1.0e-8))
    # The distance coordinate is dimensionless and starts at one: rho ∝ 1/d.
    @test gradients[1][end - 2] / reference ≈ -1 rtol = rtol
    if finite_difference
        fd = grad(central_fdm(5, 1), loss, x)[1]
        @test all(isfinite, fd)
        @test all(isapprox.(fd ./ reference, gradients[1] ./ reference; rtol = 1.0e-4, atol = 1.0e-7))
    end
    return nothing
end

function check_phase_ad(strain, x)
    shifted = copy(x)
    shifted[end] += 0.2
    template = conj.(strain(shifted))
    scale = sum(abs2, template)
    loss = FixedObjective(z -> sum(real, strain(z) .* template) / scale)
    gradients = map(backend -> gradient(loss, backend, x), SNR_BACKENDS)
    @test all(isfinite, gradients[1]) && all(isfinite, gradients[2])
    @test all(isapprox.(gradients[1], gradients[2]; rtol = 1.0e-6, atol = 1.0e-8))
    fd = grad(central_fdm(5, 1), loss, x)[1]
    @test gradients[1][end] ≈ fd[end] rtol = 1.0e-5 atol = 1.0e-7
    return @test abs(gradients[1][end]) > 1.0e-6
end
