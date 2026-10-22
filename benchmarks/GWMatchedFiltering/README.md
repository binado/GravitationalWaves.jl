# Inference storage and metadata costs

Run `just benchmark-inference` from the repository root. The harness reports
Julia/CPU configuration, grid-check costs, and raw versus checked-series inner
products for three channels over several frequency counts. Timings use warmed
BenchmarkTools samples and include allocation counts. The kernels use
`(frequency, channel)` ComplexF64 matrices and a fixed diagonal PSD.

Interpret very small absolute differences as measurement noise. The important
regression checks are allocation-free native contractions, fixed-cost metadata
comparisons, and unchanged raw numerical paths after flattening.
