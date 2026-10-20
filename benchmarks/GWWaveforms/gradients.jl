# Instantiate benchmarks/GWWaveforms and integration/GWWaveforms first, then run:
# julia --threads=1 --project=integration/GWWaveforms benchmarks/GWWaveforms/gradients.jl
push!(LOAD_PATH, @__DIR__)
using ForwardDiff
using Mooncake: Mooncake
using LinearAlgebra: norm
using FiniteDifferences: central_fdm, grad
using Test
include("scalar_kernels.jl")

struct WaveformLoss{F, D, G, T}
    implementation::F
    frequencies::D
    generator::G
    template::T
end

function gradient_params(x)
    return (;
        mass_1 = x[1], mass_2 = x[2], a_1 = x[3], a_2 = x[4],
        luminosity_distance = x[5], theta_jn = x[6], phase = x[7],
    )
end

# Phase-sensitive scalar functional of both polarizations, with strains rescaled to
# order one. All implementations use this identical downstream computation.
function (loss::WaveformLoss)(x)
    p = gradient_params(x)
    a = loss.implementation(loss.frequencies, p, loss.generator)
    h = (1.0e22 .* a.amplitude) .* cis.(a.phase)
    c = cos(p.theta_jn)
    plus = ((1 + c * c) / 2) .* h
    cross = (-im * c) .* h
    return sum(abs2, plus) + real(sum(cross .* loss.template))
end

relative_gradient_error(actual, expected) = norm(actual - expected) / max(norm(expected), eps(Float64))

function benchmark_gradient(loss, x, reference, sample_points; samples)
    output = similar(x)
    config = ForwardDiff.GradientConfig(loss, x, ForwardDiff.Chunk{7}())
    ForwardDiff.gradient!(output, loss, x, config)
    @test output ≈ reference rtol = 1.0e-10
    forward_error = relative_gradient_error(output, reference)

    cache = Mooncake.prepare_gradient_cache(loss, x)
    _, tangents = Mooncake.value_and_gradient!!(cache, loss, x)
    @test tangents[2] ≈ reference rtol = 1.0e-10
    reverse_error = relative_gradient_error(tangents[2], reference)

    # Reuse the same prepared reverse cache while crossing component ordering,
    # equal masses and many intrinsic cutoff bins. These checks are not timed.
    for point in sample_points
        ForwardDiff.gradient!(output, loss, point, config)
        value, tangents = Mooncake.value_and_gradient!!(cache, loss, point)
        @test isfinite(value) && all(isfinite, output) && all(isfinite, tangents[2])
        @test tangents[2] ≈ output rtol = 1.0e-9
    end
    # Return to the point used by every timing measurement.
    ForwardDiff.gradient!(output, loss, x, config)
    Mooncake.value_and_gradient!!(cache, loss, x)
    for _ in 1:10
        ForwardDiff.gradient!(output, loss, x, config)
        Mooncake.value_and_gradient!!(cache, loss, x)
    end
    forward = @benchmark ForwardDiff.gradient!($output, $loss, $x, $config) samples = samples evals = 1 seconds = 30
    reverse = @benchmark Mooncake.value_and_gradient!!($cache, $loss, $x) samples = samples evals = 1 seconds = 30
    return (;
        forward = BenchmarkTools.median(forward), reverse = BenchmarkTools.median(reverse),
        forward_error, reverse_error,
    )
end

function run_gradient_benchmarks(;
        samples = 100, spacings = (2.0, 0.25),
        implementations = (("scalar_loop", scalar_loop), ("tuple_broadcast", broadcasted), ("separate_broadcast", separate_broadcast)),
    )
    x = [36.0, 29.0, 0.5, 0.3, 500.0, 0.4, 0.2]
    points = (
        [180.0, 120.0, 0.3, 0.2, 500.0, 0.4, 0.2],
        [30.0, 30.0, 0.3, 0.2, 500.0, 0.4, 0.2],
        [29.0, 36.0, 0.3, 0.5, 500.0, 0.4, 0.2],
    )
    println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads())
    println("Prepared gradients of a phase-sensitive polarization loss with respect to seven parameters.")
    println("Compilation/cache preparation excluded; warmed median over ", samples, " single evaluations.")
    println("model,bins,implementation,backend,median_us,bytes,allocations,relative_gradient_error")
    flush(stdout)
    return @testset "gradient implementations" begin
        for GT in (IMRPhenomD, IMRPhenomXAS), df in spacings
            # These binary spacings are represented exactly, so this has the same
            # samples as 20.0:df:1024.0 without TwicePrecision range metadata. Mooncake
            # 0.5 cannot copy that metadata when reusing some compiled broadcast rules.
            f = StepRangeLen(20.0, df, floor(Int, (1024.0 - 20.0) / df) + 1)
            @test collect(f) == collect(20.0:df:1024.0)
            g = GT(; minimum_frequency = 20.0, maximum_frequency = 1000.0)
            template_params = gradient_params(x .* [1.01, 1, 1, 1, 1, 1, 1])
            a = amplitude_phase(f, template_params, g)
            template = conj.((-im * cos(template_params.theta_jn) * 1.0e22) .* a.amplitude .* cis.(a.phase))
            baseline = WaveformLoss(scalar_loop, f, g, template)
            reference = ForwardDiff.gradient(baseline, x)
            @test all(isfinite, reference)
            if df == first(spacings)
                finite_difference = grad(central_fdm(5, 1), baseline, x)[1]
                @test reference ≈ finite_difference rtol = 1.0e-5
            end
            for (name, implementation) in implementations
                loss = WaveformLoss(implementation, f, g, template)
                @test loss(x) == baseline(x)
                for point in points
                    @test ForwardDiff.gradient(loss, point) ≈ ForwardDiff.gradient(baseline, point) rtol = 1.0e-10
                end
                GC.gc()
                stats = benchmark_gradient(loss, x, reference, points; samples)
                for (backend, estimate, error) in (
                        ("ForwardDiff", stats.forward, stats.forward_error),
                        ("Mooncake", stats.reverse, stats.reverse_error),
                    )
                    @printf(
                        "%s,%d,%s,%s,%.3f,%d,%d,%.3e\n", nameof(GT), length(f), name,
                        backend, estimate.time / 1000, estimate.memory, estimate.allocs, error
                    )
                end
                flush(stdout)
            end
        end
    end
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    run_gradient_benchmarks()
end
