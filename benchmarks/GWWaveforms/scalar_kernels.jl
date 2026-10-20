# julia --project=benchmarks/GWWaveforms -e 'using Pkg; Pkg.instantiate()'
# julia --project=benchmarks/GWWaveforms benchmarks/GWWaveforms/scalar_kernels.jl
using BenchmarkTools
using GWWaveforms
using Printf

include("previous_loops.jl")
const G = GWWaveforms

function previous(f, params, g)
    T, p = G.promote_parameters(G.aligned_parameters(params))
    f_min, f_max = G.frequency_bounds(g, f)
    f_start = g.pn_amplitude_order != 0 ? G.f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min
    return previous_fd(g, T, f, p, f_start, f_min, f_max)
end

function previous_fd(g::IMRPhenomD, T, f, p, f_start, f_min, f_max)
    lo, hi = G.phenom_frequency_slice(f, f_start, f_min, f_max)
    return PreviousLoops.baseline_phenomd(T, f, lo, hi, p, g, f_start, f_max)
end
previous_fd(g::IMRPhenomXAS, T, f, p, f_start, f_min, f_max) =
    PreviousLoops.baseline_phenomxas(g, T, f, p, f_start, f_min, f_max)

function broadcast_samples(kernel::K, ::Type{T}, f, lo, hi, model, scale, cutoff) where {K, T}
    amplitude = zeros(T, length(f))
    phase = zeros(T, length(f))
    values = kernel.(Float64.(view(f, lo:hi)), Ref(model), scale, cutoff)
    amplitude[lo:hi] .= first.(values)
    phase[lo:hi] .= last.(values)
    return (; amplitude, phase)
end

function broadcasted(f, params, g, assemble::A = broadcast_samples) where {A}
    T, p = G.promote_parameters(G.aligned_parameters(params))
    f_min, f_max = G.frequency_bounds(g, f)
    f_start = g.pn_amplitude_order != 0 ? G.f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min
    lo, hi = G.phenom_frequency_slice(f, f_start, f_min, f_max)
    lo > hi && return (; amplitude = zeros(T, length(f)), phase = zeros(T, length(f)))
    return broadcasted_fd(g, T, f, lo, hi, p, f_start, f_max, assemble)
end

function broadcasted_fd(g::IMRPhenomD, T, f, lo, hi, p, f_start, f_max, assemble)
    fref = G.phenom_reference_frequency(g, f_start)
    model = G.phenomd_model(p.mass_1, p.mass_2, p.chi_1, p.chi_2, fref, p.phase)
    scale = G.phenom_amplitude_scale(p.mass_1, p.mass_2, p.luminosity_distance)
    cutoff = G.F_CUT / ((p.mass_1 + p.mass_2) * G.MTSUN_SI)
    df = Float64(step(f))
    stop = floor(min(cutoff, f_max) / df) * df
    return assemble(G.phenomd_amplitude_phase, T, f, lo, hi, model, scale, stop)
end

function broadcasted_fd(g::IMRPhenomXAS, T, f, lo, hi, p, f_start, f_max, assemble)
    fref = G.phenom_reference_frequency(g, f_start)
    model = G.phenomxas_model(p.mass_1, p.mass_2, p.chi_1, p.chi_2, fref, p.phase)
    amp0 = G.phenom_amplitude_scale(p.mass_1, p.mass_2, p.luminosity_distance)
    scale = amp0 * sqrt(2model.r.eta / 3) * Float64(π)^(-1 / 6)
    return assemble(G.phenomxas_amplitude_phase, T, f, lo, hi, model, scale, 0.3 / model.M_sec)
end

function scalar_loop_samples(kernel::K, ::Type{T}, f, lo, hi, model, scale, cutoff) where {K, T}
    amplitude = zeros(T, length(f))
    phase = zeros(T, length(f))
    for i in lo:hi
        amplitude[i], phase[i] = kernel(Float64(f[i]), model, scale, cutoff)
    end
    return (; amplitude, phase)
end
scalar_loop(f, params, g) = broadcasted(f, params, g, scalar_loop_samples)
const separate_broadcast = amplitude_phase

function run_benchmarks()
    params = (;
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0,
        theta_jn = 0.4, phase = 0.2, a_1 = 0.5, a_2 = 0.3,
    )
    println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads())
    println("End-to-end amplitude_phase, including coefficient preparation and output allocation.")
    println("Trials run after warm-up; times are median and minimum over 200 single evaluations.")
    println("model,bins,implementation,median_us,min_us,bytes,allocations")
    for GType in (IMRPhenomD, IMRPhenomXAS), df in (2.0, 0.25, 0.03125)
        f = 20.0:df:1024.0
        g = GType(; minimum_frequency = 20.0, maximum_frequency = 1000.0)
        expected = previous(f, params, g)
        for (name, implementation) in (("previous_loop", previous), ("scalar_loop", scalar_loop), ("tuple_broadcast", broadcasted), ("separate_broadcast", separate_broadcast))
            actual = implementation(f, params, g)
            @assert actual.amplitude == expected.amplitude
            @assert actual.phase == expected.phase
            # Warm every specialization before collecting timing or allocation results.
            for _ in 1:10
                implementation(f, params, g)
            end
            trial = @benchmark $implementation($f, $params, $g) samples = 200 evals = 1 seconds = 10
            middle = BenchmarkTools.median(trial)
            fastest = minimum(trial)
            @printf(
                "%s,%d,%s,%.3f,%.3f,%d,%d\n", nameof(GType), length(f), name,
                middle.time / 1000, fastest.time / 1000, middle.memory, middle.allocs
            )
        end
    end
    return
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_benchmarks()
