# Run with --project=integration/GWWaveforms after instantiating that environment and
# benchmarks/GWWaveforms. Both primal and AD timings exclude compilation/preparation.
include("gradients.jl")

function run_ref_benchmarks()
    params = (;
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0,
        theta_jn = 0.4, phase = 0.2, a_1 = 0.5, a_2 = 0.3,
    )
    implementations = (("anonymous_broadcast", anonymous_broadcast), ("ref_broadcast", amplitude_phase))
    println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads())
    println("Primal warmed medians, 200 single evaluations; ordinary frequency ranges.")
    println("model,bins,implementation,median_us,bytes,allocations")
    for GT in (IMRPhenomD, IMRPhenomXAS), df in (0.25, 0.03125)
        f = 20.0:df:1024.0
        g = GT(; minimum_frequency = 20.0, maximum_frequency = 1000.0)
        expected = anonymous_broadcast(f, params, g)
        for (name, implementation) in implementations
            actual = implementation(f, params, g)
            @assert actual.amplitude == expected.amplitude
            @assert actual.phase == expected.phase
            for _ in 1:10
                implementation(f, params, g)
            end
            trial = @benchmark $implementation($f, $params, $g) samples = 200 evals = 1 seconds = 10
            estimate = BenchmarkTools.median(trial)
            @printf(
                "%s,%d,%s,%.3f,%d,%d\n", nameof(GT), length(f), name,
                estimate.time / 1000, estimate.memory, estimate.allocs
            )
            flush(stdout)
        end
    end
    return run_gradient_benchmarks(; spacings = (0.25,), implementations)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    run_ref_benchmarks()
end
