using BenchmarkTools
using GWSeries
using GWMatchedFiltering
using Random

println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads())
println("ComplexF64 matrices, (frequency, channel) storage; 3 ordered channels")
rng = Xoshiro(37)
for n in (16, 1024, 8192, 131072)
    a = randn(rng, ComplexF64, n, 3)
    b = randn(rng, ComplexF64, n, 3)
    sa = FrequencySeries(a; epoch = 1126259460.125, frequency_start = 0, frequency_spacing = 0.25, channels = (:H1, :L1, :V1))
    sb = FrequencySeries(b; metadata(sa)...)
    metric = FrequencyDomainMetric(sa, ones(n, 3))
    println("\n", n, " bins")
    for (label, trial) in (
            "grid comparison" => @benchmark(check_same_grid($sa, $sb); seconds = 0.3),
            "raw contraction" => @benchmark(complex_inner_product($a, $b, $metric); seconds = 0.3),
            "checked series contraction" => @benchmark(complex_inner_product($sa, $sb, $metric); seconds = 0.3),
        )
        best = minimum(trial)
        println(rpad(label, 30), round(best.time; digits = 1), " ns; ", best.allocs, " allocations; ", best.memory, " bytes")
    end
end
