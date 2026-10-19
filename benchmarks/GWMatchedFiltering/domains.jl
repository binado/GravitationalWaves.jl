include("layouts.jl")

function raw_inner_product(a, b, metric, df)
    bins = metric.active_bins
    return 4df * GWMatchedFiltering.contraction(GWMatchedFiltering.active_values(a, bins), GWMatchedFiltering.active_values(b, bins), metric.metric)
end

raw_inner_product(a::AbstractVector, b::AbstractVector, metric, df) =
    4df * GWMatchedFiltering.contraction(a, b, metric.metric, metric.active_bins)

function benchmark_domains()
    println("# Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads(), "; BLAS threads 1")
    println("case,time_ns,bytes,allocations")
    rng = Xoshiro(31)
    for D in (1, 3)
        names = D == 1 ? (:H1,) : (:H1, :L1, :V1)
        fd = FrequencyDomain(; frequencies = range(0.0; step = 0.25, length = 8194), epoch = 1126259460.0, channels = names)
        a = randn(rng, ComplexF64, nsamples(fd), D)
        b = randn(rng, ComplexF64, nsamples(fd), D)
        if D == 1
            a, b = vec(a), vec(b)
        end
        psds = fill(2.0, nsamples(fd), D)
        block = [2.0 0.2 0.1;0.2 2.0 0.3;0.1 0.3 2.0][1:D, 1:D]
        blocks = repeat(reshape(block, 1, D, D), nsamples(fd), 1, 1)
        df = frequency_spacing(fd)
        layouts = D == 1 ? (("vector", identity), ("view", x -> view(x, :))) : (("matrix", identity), ("view", x -> view(x, :, :)))
        for (kind, spectrum) in (("diagonal", psds), ("blocks", blocks)), (layout, storage) in layouts
            metric = FrequencyDomainMetric(fd, spectrum)
            aa, bb = storage(a), storage(b)
            check_data(fd, aa);check_data(fd, bb);check_compatible(fd, metric)
            @assert complex_inner_product(fd, aa, bb, metric) == raw_inner_product(aa, bb, metric, df)
            record("D=$D/$kind/$layout/kernel", @benchmarkable raw_inner_product($aa, $bb, $metric, $df))
            record("D=$D/$kind/$layout/domain", @benchmarkable complex_inner_product($fd, $aa, $bb, $metric))
        end
    end
    return nothing
end
benchmark_domains()
