# Focused wrapper comparison; run without competing compilation for useful timings.
include("layouts.jl")

function benchmark_series()
    println("# Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads(), "; BLAS threads 1")
    println("case,time_ns,bytes,allocations")
    f = range(0.0; step = 0.25, length = 8194)
    names = (:H1, :L1, :V1)
    rng = Xoshiro(31)
    a = randn(rng, ComplexF64, length(f), 3)
    b = randn(rng, ComplexF64, length(f), 3)
    psds = fill(2.0, length(f), 3)
    blocks = repeat(reshape([2.0 0.2 0.1; 0.2 2.0 0.3; 0.1 0.3 2.0], 1, 3, 3), length(f), 1, 1)
    for (kind, spectrum) in (("diagonal", psds), ("blocks", blocks)), (layout, storage) in (("matrix", identity), ("view", x -> view(x, :, :)))
        metric = FrequencyDomainMetric(f, names, spectrum)
        aa, bb = storage(a), storage(b)
        h = FrequencySeries(aa; frequencies = f, epoch = 1126259460.0, channels = names)
        d = FrequencySeries(bb; frequencies = f, epoch = epoch(h), channels = names)
        @assert complex_inner_product(h, d, metric) == complex_inner_product(aa, bb, metric)
        record("$kind/$layout/array", @benchmarkable complex_inner_product($aa, $bb, $metric))
        record("$kind/$layout/series", @benchmarkable complex_inner_product($h, $d, $metric))
    end
    return nothing
end

benchmark_series()
