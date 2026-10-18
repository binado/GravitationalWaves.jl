using BenchmarkTools
using GWMatchedFiltering
using GWSeries
using GWDetectors
using GWWaveforms
using LinearAlgebra
using Random
using StaticArrays
using Statistics: median

BLAS.set_num_threads(1)

# Each loop applies exactly the same whitening operator. Layout is a compile-time
# choice, and conversion is excluded from contraction timings but included below
# in the waveform/projection path.
@inline sample(a, k, i, ::Val{:frequency_first}) = a[k, i]
@inline sample(a, k, i, ::Val{:channel_first}) = a[i, k]
@inline sample(a, k, i, ::Val{:block_contiguous}) = a[k, i]
@inline factor(w, k, i, j, ::Val{:frequency_first}) = w[k, i, j]
@inline factor(w, k, i, j, ::Val{:channel_first}) = w[i, j, k]
@inline factor(w, k, i, j, ::Val{:block_contiguous}) = w[i, j, k]

function scalar_blocks(a, b, w, layout, ::Val{D}, nf) where {D}
    result = zero(promote_type(eltype(a), eltype(b), eltype(w)))
    @inbounds for k in 1:nf
        for i in 1:D
            wa = zero(result)
            wb = zero(result)
            for j in 1:D
                wa += factor(w, k, i, j, layout) * sample(a, k, j, layout)
                wb += factor(w, k, i, j, layout) * sample(b, k, j, layout)
            end
            result += conj(wa) * wb
        end
    end
    return result
end

function static_blocks(a, b, w, layout, ::Val{D}, nf) where {D}
    result = zero(promote_type(eltype(a), eltype(b), eltype(w)))
    @inbounds for k in 1:nf
        aa = SVector{D}(ntuple(i -> sample(a, k, i, layout), Val(D)))
        bb = SVector{D}(ntuple(i -> sample(b, k, i, layout), Val(D)))
        W = SMatrix{D, D}(ntuple(n -> factor(w, k, mod1(n, D), cld(n, D), layout), Val(D * D)))
        result += dot(W * aa, W * bb)
    end
    return result
end

diagonal(a, b, w) = sum(conj.(a) .* w .* b)

function record(label, benchmark)
    estimate = median(run(benchmark; seconds = 0.3))
    println(label, ",", estimate.time, ",", estimate.memory, ",", estimate.allocs)
    flush(stdout)
    return nothing
end

function benchmark_layouts()
    println("# Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads(), "; BLAS threads 1")
    println("case,time_ns,bytes,allocations")
    rng = Xoshiro(81)
    for D in (1, 3, 6), nf in (1025, 8193, 32769)
        f = range(0.0; step = 0.25, length = nf + 1)
        names = Tuple(Symbol("D", i) for i in 1:D)
        a = randn(rng, ComplexF64, nf, D)
        b = randn(rng, ComplexF64, nf, D)
        s0 = Matrix{ComplexF64}(I, D, D) + fill(0.1, D, D)
        S = repeat(reshape(s0, 1, D, D), nf + 1, 1, 1)
        m = FrequencyDomainMetric(f, names, S)
        w = permutedims(m.metric.whiteners, (3, 1, 2))
        psd = fill(2.0, nf + 1, D)
        md = FrequencyDomainMetric(f, names, psd)
        full_a = vcat(zeros(ComplexF64, 1, D), a)
        full_b = vcat(zeros(ComplexF64, 1, D), b)
        reference = complex_inner_product(full_a, full_b, m) / m.normalization
        for layout in (Val(:frequency_first), Val(:channel_first), Val(:block_contiguous))
            aa = layout isa Val{:channel_first} ? permutedims(a) : a
            bb = layout isa Val{:channel_first} ? permutedims(b) : b
            ww = layout isa Val{:frequency_first} ? w : permutedims(w, (2, 3, 1))
            weights = fill(0.5, size(aa))
            @assert scalar_blocks(aa, bb, ww, layout, Val(D), nf) ≈ reference
            @assert static_blocks(aa, bb, ww, layout, Val(D), nf) ≈ reference
            label = "D=$D/nf=$nf/$(typeof(layout))"
            record("$label/scalar", @benchmarkable scalar_blocks($aa, $bb, $ww, $layout, $(Val(D)), $nf))
            record("$label/static", @benchmarkable static_blocks($aa, $bb, $ww, $layout, $(Val(D)), $nf))
            record("$label/diagonal", @benchmarkable diagonal($aa, $bb, $weights))
        end
        series_a = FrequencySeries(full_a; frequencies = f, epoch = 1126259460.0, channels = names)
        series_b = FrequencySeries(full_b; frequencies = f, epoch = 1126259460.0, channels = names)
        @assert complex_inner_product(series_a, series_b, m) ≈ m.normalization * reference
        record("D=$D/nf=$nf/series/blocks", @benchmarkable complex_inner_product($series_a, $series_b, $m))
        record("D=$D/nf=$nf/series/diagonal", @benchmarkable complex_inner_product($series_a, $series_b, $md))
        record("D=$D/nf=$nf/production/blocks", @benchmarkable complex_inner_product($full_a, $full_b, $m))
        record("D=$D/nf=$nf/production/diagonal", @benchmarkable complex_inner_product($full_a, $full_b, $md))
        record("D=$D/nf=$nf/production/norm", @benchmarkable optimal_snr_squared($full_a, $m))
        record("D=$D/nf=$nf/production/matched", @benchmarkable matched_filter_snr($full_a, $full_b, $m))
        record("D=$D/nf=$nf/production/residual", @benchmarkable residual_squared($full_a, $full_b, $m))
        record("D=$D/nf=$nf/setup/blocks", @benchmarkable FrequencyDomainMetric($f, $names, $S))
        record("D=$D/nf=$nf/setup/diagonal", @benchmarkable FrequencyDomainMetric($f, $names, $psd))
    end
    return nothing
end

function waveform_path(f, generator, params, sky, w, layout)
    pol = hphc_fd(f, params, generator)
    strains = hcat([project_fd(f, pol, sky, d; start_time = sky.geocent_time - 2) for d in (detectors.H1, detectors.L1, detectors.V1)]...)
    selected = strains[2:end, :]
    a = layout isa Val{:channel_first} ? permutedims(selected) : selected
    return scalar_blocks(a, a, w, layout, Val(3), length(f) - 1)
end

function benchmark_pipeline()
    f = 0.0:0.25:1024.0
    generator = TaylorF2(; maximum_frequency = last(f))
    params = (; mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3)
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    W = repeat(reshape(Matrix{Float64}(I, 3, 3), 1, 3, 3), length(f) - 1, 1, 1)
    for layout in (Val(:frequency_first), Val(:channel_first), Val(:block_contiguous))
        w = layout isa Val{:frequency_first} ? W : permutedims(W, (2, 3, 1))
        record("pipeline/$(typeof(layout))", @benchmarkable waveform_path($f, $generator, $params, $sky, $w, $layout))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_layouts()
    benchmark_pipeline()
end
