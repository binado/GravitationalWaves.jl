using BenchmarkTools
using GWMatchedFiltering
using LinearAlgebra
using Random

println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads())
println("ComplexF64 matrices, (frequency, channel) data; 3 ordered channels")
println("diagonal W: (frequency, channel); block W: (channel, channel, frequency)")
rng = Xoshiro(37)
df = 0.25
S0 = ComplexF64[2 0.2 + 0.1im 0.1; 0.2 - 0.1im 2 0.3im; 0.1 -0.3im 2]
L0 = inv(cholesky(Hermitian(S0)).L)
for n in (16, 1024, 8192, 131072)
    a = randn(rng, ComplexF64, n, 3)
    b = randn(rng, ComplexF64, n, 3)
    W_diagonal = fill(1 / sqrt(2.0), n, 3)
    W_blocks = repeat(Matrix(L0), 1, 1, n)
    println("\n", n, " bins")
    for (label, trial) in (
            "diagonal contraction" => @benchmark(complex_inner_product($a, $b, $W_diagonal, $df); seconds = 0.3),
            "block contraction" => @benchmark(complex_inner_product($a, $b, $W_blocks, $df); seconds = 0.3),
            "block whitening" => @benchmark(whiten($a, $W_blocks); seconds = 0.3),
        )
        best = minimum(trial)
        println(rpad(label, 30), round(best.time; digits = 1), " ns; ", best.allocs, " allocations; ", best.memory, " bytes")
    end
end
