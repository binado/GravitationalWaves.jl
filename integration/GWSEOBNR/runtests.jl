using GWSEOBNR, GWDetectors, GWMatchedFiltering, GWTimeFrequencyDomains
using ForwardDiff, LinearAlgebra, NPZ, TOML, Test, FFTW
const S = GWSEOBNR
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const F = 0:0.0625:2048
const G = SEOBNRv5HM()
const STEPS = 10.0 .^ range(-2, -5; length = 5)
const PSD = [f >= 20 ? 1 + (40 / f)^4 + (f / 200)^2 : 1.0 for f in F]
const WEIGHTS = vcat([f >= 20 ? 1 / s : 0.0 for (f, s) in zip(F, PSD)], [f >= 20 ? 1 / s : 0.0 for (f, s) in zip(F, PSD)])
const DIAGNOSTICS = Dict{String, Any}("julia" => string(VERSION), "cases" => [])

function polarizations(x, g = G)
    p = (;
        mass_1 = x[1], mass_2 = x[2], chi_1 = x[3], chi_2 = x[4],
        luminosity_distance = x[5], inclination = x[6], phase = x[7],
    )
    return S.evaluated_spectrum(F, p, g)
end
function waveform(x, g = G)
    h = polarizations(x, g)
    return vcat(h.plus, h.cross)
end
relative(a, b) = norm(a - b) / norm(b)
fivepoint(fun, x, v, eps) = (-fun(x + 2eps * v) + 8fun(x + eps * v) - 8fun(x - eps * v) + fun(x - 2eps * v)) / (12eps)
adjacent_pass(errors, tol) = any(errors[i] < tol && errors[i + 1] < tol for i in 1:(length(errors) - 1))

function mismatches(a, b)
    product = conj.(b) .* WEIGHTS .* a
    normalizer = sqrt(real(dot(a, WEIGHTS .* a) * dot(b, WEIGHTS .* b)))
    unoptimized = 1 - real(sum(product)) / normalizer
    # Optimized time/phase mismatch is diagnostic only. A grid locates the
    # maximum correlation, followed by a continuous golden-section search.
    frequency = vcat(collect(F), collect(F))
    overlap(t) = abs(sum(product .* cis.(2π .* frequency .* t))) / normalizer
    times = range(-0.002, 0.002; length = 81)
    i = argmax(overlap.(times))
    lo, hi = times[max(i - 1, 1)], times[min(i + 1, length(times))]
    for _ in 1:50
        t1, t2 = lo + (hi - lo) * 0.3819660112501051, lo + (hi - lo) * 0.6180339887498949
        if overlap(t1) > overlap(t2)
            hi = t2
        else
            lo = t1
        end
    end
    return unoptimized, 1 - overlap((lo + hi) / 2)
end

function case_acceptance(M, q, s1, s2)
    println(stderr, "Acceptance case M=$M q=$q"); flush(stderr)
    record = Dict{String, Any}("mass" => M, "q" => q, "chi1" => s1, "chi2" => s2)
    push!(DIAGNOSTICS["cases"], record)
    x = [M * q / (1 + q), M / (1 + q), s1, s2, 500.0, 0.7, 0.3]
    # Warm each scalar type before measuring.
    h = waveform(x)
    timed = @timed waveform(x)
    record["seconds"], record["bytes"] = timed.time, timed.bytes
    oracle = npzread(joinpath(ROOT, "reference/seobnr/M$(Int(M))_q$(Int(q))_fd.npz"))
    ref = vcat(oracle["plus"], oracle["cross"])
    record["spectral_l2"] = relative(h, ref)
    record["mismatch"], record["optimized_mismatch"] = mismatches(h, ref)
    @test record["spectral_l2"] < 1.0e-3
    @test record["mismatch"] < 1.0e-5
    e = S.evolve(x[1], x[2], s1, s2, G)
    nqc = S.nqc_coefficients(e)
    record["t_isco"], record["t_attach"], record["t_ref"] = S.primal(e.t_isco), S.primal(e.t_attach), S.primal(e.t_ref)
    record["termination"] = string(e.sol.retcode)
    record["steps"] = e.sol.stats.naccept
    record["rejected_steps"] = e.sol.stats.nreject
    record["nqc_a"], record["nqc_b"] = collect(nqc.a), collect(nqc.b)
    mr = S.merger_coefficients(e, nqc)
    peak = S.peak_time(e, nqc, mr)
    record["peak_time"] = S.primal(peak)
    record["peak_at_attachment"] = abs(S.primal(peak - e.t_attach)) < 1.0e-10
    record["positive_delta_t"] = e.dt_nr > 0
    record["qnm_positive_spin"] = e.af > 0
    record["ringdown_c1f_limited"] = mr.c1f != S.mr_c1f(e.c)
    reference_nqc = oracle["nqc"][1:5]
    record["nqc_relative"] = relative(vcat(nqc.a, nqc.b), reference_nqc)
    @test record["nqc_relative"] < 1.0e-9
    # Symmetry-preserving mass derivative at equal mass; independent masses only
    # in the interior of the ordered domain. Spins are signed internally.
    directions = q == 1 ? [[0.4, 0.4, 0, 0, 0, 0, 0]] :
        [[1.0, 0, 0, 0, 0, 0, 0], [0.0, 1, 0, 0, 0, 0, 0]]
    append!(directions, [[0.0, 0, 1, 0, 0, 0, 0], [0.0, 0, 0, 1, 0, 0, 0]])
    if q == 1
        # Exchange symmetry at zero spin: perturb both signed spins together.
        directions = [directions[1], [0.0, 0, 1, 1, 0, 0, 0]]
    end
    errors = Vector{Float64}[]
    for v in directions
        ad = ForwardDiff.derivative(t -> waveform(x + t * v), 0.0)
        err = [relative(ad, fivepoint(waveform, x, v, eps)) for eps in STEPS]
        push!(errors, err)
        @test all(isfinite, ad)
        @test adjacent_pass(err, 1.0e-3)
    end
    record["intrinsic_derivative_errors"] = errors
    for j in 5:7
        v = zeros(7)
        v[j] = j == 5 ? 100.0 : 1.0
        # Match public extrinsic-only calls: intrinsic storage remains primal.
        extrinsic = function (t)
            p = (;
                mass_1 = x[1], mass_2 = x[2], chi_1 = x[3], chi_2 = x[4],
                luminosity_distance = x[5] + t * v[5],
                inclination = x[6] + t * v[6], phase = x[7] + t * v[7],
            )
            extrinsic_pol = S.evaluated_spectrum(F, p, G)
            return vcat(extrinsic_pol.plus, extrinsic_pol.cross)
        end
        ad = ForwardDiff.derivative(extrinsic, 0.0)
        eps = 1.0e-3
        fd = (-extrinsic(2eps) + 8extrinsic(eps) - 8extrinsic(-eps) + extrinsic(-2eps)) / (12eps)
        record["extrinsic_$j"] = relative(ad, fd)
        @test relative(ad, fd) < 1.0e-3
        j == 5 && @test ad ≈ -h / 5 rtol = 1.0e-10
    end
    tighter = SEOBNRv5HM(reltol = 1.0e-12, abstol = 1.0e-13)
    record["waveform_convergence"] = relative(waveform(x, tighter), h)
    @test record["waveform_convergence"] < 1.0e-4
    v = directions[1]
    ad = ForwardDiff.derivative(t -> waveform(x + t * v), 0.0)
    ad_tight = ForwardDiff.derivative(t -> waveform(x + t * v, tighter), 0.0)
    record["derivative_convergence"] = relative(ad_tight, ad)
    @test record["derivative_convergence"] < 1.0e-3
    gradient_eval = t -> waveform(x + t * v)
    ForwardDiff.derivative(gradient_eval, 0.0)
    timed_ad = @timed ForwardDiff.derivative(gradient_eval, 0.0)
    return record["gradient_seconds"], record["gradient_bytes"] = timed_ad.time, timed_ad.bytes
end

function composition_acceptance()
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    d1 = FrequencyDomain(; frequencies = F, epoch = 1126259460.0, channels = (:H1,))
    d2 = FrequencyDomain(; frequencies = F, epoch = 1126259460.0, channels = (:L1,))
    domain = FrequencyDomain(; frequencies = F, epoch = 1126259460.0, channels = (:H1, :L1))
    metric = FrequencyDomainMetric(domain, hcat(PSD, PSD) .* 1.0e-46; minimum_frequency = 20.0)
    project(x) = begin
        h = polarizations(x)
        hcat(project_fd(d1, h, sky, detectors.H1), project_fd(d2, h, sky, detectors.L1))
    end
    injected = [40.0, 20.0, 0.4, -0.3, 500.0, 0.7, 0.3]
    data = project(injected)
    x = injected + [0.1, -0.05, 0.001, -0.002, 1.0, 0.003, -0.002]
    loss(x) = log_likelihood_ratio(domain, project(x), data, metric)
    grad = ForwardDiff.gradient(loss, x)
    fd = [fivepoint(loss, x, [i == j ? 1.0 : 0.0 for i in eachindex(x)], 1.0e-3) for j in eachindex(x)]
    DIAGNOSTICS["composition_relative"] = relative(grad, fd)
    DIAGNOSTICS["composition_gradient"] = grad
    @test norm(grad) > 0
    return @test relative(grad, fd) < 1.0e-3
end

try
    @testset "SEOBNR complete acceptance (gates are not relaxed)" begin
        for M in (40.0, 60.0, 100.0), (q, s1, s2) in ((1.0, 0.0, 0.0), (2.0, 0.4, -0.3), (4.0, -0.6, 0.6))
            @testset "M=$M q=$q" begin
                try
                    case_acceptance(M, q, s1, s2)
                finally
                    mkpath(joinpath(ROOT, "reference/seobnr"))
                    open(joinpath(ROOT, "reference/seobnr/diagnostics.toml"), "w") do io
                        TOML.print(io, DIAGNOSTICS)
                    end
                end
            end
        end
        @testset "offset H1/L1 fixed-PSD likelihood" composition_acceptance()
    end
finally
    dest = joinpath(ROOT, "reference/seobnr/diagnostics.toml")
    mkpath(dirname(dest))
    open(dest, "w") do io
        TOML.print(io, DIAGNOSTICS)
    end
    @info "Acceptance diagnostics saved" dest
end
