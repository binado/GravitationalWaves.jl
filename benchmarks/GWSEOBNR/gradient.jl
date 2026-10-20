# Primal and ForwardDiff gradient cost of the SEOBNRv5HM spectrum.
# Reports minimum and median wall time, allocations, Julia version and CPU.
using ForwardDiff, GWSEOBNR
const S = GWSEOBNR

const GRID = 0:0.125:1024
const GENERATOR = SEOBNRv5HM()
const SAMPLES = 5

median(v) = sort(v)[cld(length(v), 2)]

# (mass_1, mass_2, chi_1, chi_2, luminosity_distance, inclination, phase)
const POINTS = (
    "M=40 q=1" => [20.0, 20.0, 0.0, 0.0, 500.0, 0.7, 0.3],
    "M=60 q=2 spinning" => [40.0, 20.0, 0.4, -0.3, 500.0, 0.7, 0.3],
)

function loss(x)
    h = S.evaluated_spectrum(
        GRID, (;
            mass_1 = x[1], mass_2 = x[2], chi_1 = x[3], chi_2 = x[4],
            luminosity_distance = x[5], inclination = x[6], phase = x[7],
        ), GENERATOR,
    )
    return 1.0e22 * sum(real(h.plus[k]) + 0.3imag(h.plus[k]) + 0.7real(h.cross[k]) for k in 161:37:4001)
end

function measure(f)
    f() # compile
    times = Float64[]
    bytes = Int[]
    for _ in 1:SAMPLES
        stats = @timed f()
        push!(times, stats.time)
        push!(bytes, stats.bytes)
    end
    return (; min = minimum(times), median = median(times), alloc_mib = median(bytes) / 2^20)
end

function report(label, m)
    return println(
        rpad(label, 34), lpad(string(round(m.min; digits = 3)), 9), " s min  ",
        lpad(string(round(m.median; digits = 3)), 9), " s median  ",
        lpad(string(round(m.alloc_mib; digits = 1)), 10), " MiB"
    )
end

println("Julia ", VERSION, "  threads=", Threads.nthreads(), "  CPU: ", Sys.cpu_info()[1].model)
println("grid length ", length(GRID), ", samples ", SAMPLES, "\n")

for (name, x) in POINTS
    println("== ", name)
    report("primal", measure(() -> loss(x)))
    # Chunk{7} is the full-gradient single pass; smaller chunks need ceil(7 / N) passes.
    for n in (1, 2, 4, 7)
        cfg = ForwardDiff.GradientConfig(loss, x, ForwardDiff.Chunk{n}())
        report("gradient, chunk $n", measure(() -> ForwardDiff.gradient(loss, x, cfg)))
    end
    cfg = ForwardDiff.GradientConfig(loss, x)
    println("default chunk size: ", ForwardDiff.chunksize(cfg))
    println()
end
