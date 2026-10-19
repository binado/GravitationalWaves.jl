# Exploratory reverse-mode checks; failures are recorded rather than treated as support.
using GWSEOBNR, ForwardDiff, Enzyme, LinearAlgebra, TOML
using GWDetectors, GWMatchedFiltering, GWTimeFrequencyDomains
const S = GWSEOBNR
const GRID = 0:0.125:1024
const GENERATOR = SEOBNRv5HM()
const DRAW = [40.0, 20.0, 0.4, -0.3, 500.0, 0.7, 0.3]

function pol(x)
    p = (;
        mass_1 = x[1], mass_2 = x[2], chi_1 = x[3], chi_2 = x[4],
        luminosity_distance = x[5], inclination = x[6], phase = x[7],
    )
    return S.evaluated_spectrum(GRID, p, GENERATOR)
end

calibration_loss(x) = S.hamiltonian_coefficients(x...).dSO
function hamiltonian_loss(x)
    c = S.hamiltonian_coefficients(x...)
    return S.hamiltonian((12.0, 0.3), (-0.005, 4.0), c)[1]
end
function flux_loss(x)
    c = S.hamiltonian_coefficients(x...)
    d = S.hamiltonian_dynamics((12.0, 0.3), (-0.005, 4.0), c)
    wc = S.orbital_frequency((12.0, 0.3), (0.0, 4.0), c)
    return 1.0e5 * S.radiation_flux(4.0, d[4], wc, c.nu * d[5], c, S.flux_coefficients(c))
end
function initial_loss(x)
    u = S.initial_conditions(x..., π * (x[1] + x[2]) * S.MTSUN_SI * 16)
    return u[1] + 3u[3] + 2u[4]
end
function evolution_loss(x)
    e = S.evolve(x..., GENERATOR)
    return e.t_attach / 1000 + e.phi_ref
end
timing_loss(x) = S.continuous_root(t -> t^2 - x[1], 0.1, 3.0, t -> 2t)
fft_loss(x) = sum(abs2, S.real_fft(x))
tail_loss(x) = real(S.complex_tail(2, 2, x[1], x[2])) + 0.3imag(S.complex_tail(2, 2, x[1], x[2]))
function waveform_loss(x)
    h = pol(x)
    return 1.0e22 * sum(real(h.plus[k]) + 0.3imag(h.plus[k]) + 0.7real(h.cross[k]) for k in 161:37:4001)
end

const SKY = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
const H1_DOMAIN = FrequencyDomain(; frequencies = GRID, epoch = 1126259460.0, channels = (:H1,))
const L1_DOMAIN = FrequencyDomain(; frequencies = GRID, epoch = 1126259460.0, channels = (:L1,))
const DOMAIN = FrequencyDomain(; frequencies = GRID, epoch = 1126259460.0, channels = (:H1, :L1))
const PSD = [f >= 20 ? 1 + (40 / f)^4 + (f / 200)^2 : 1.0 for f in GRID]
const METRIC = FrequencyDomainMetric(DOMAIN, hcat(PSD, PSD) .* 1.0e-46; minimum_frequency = 20.0)
function signal(x)
    h = pol(x)
    return hcat(project_fd(H1_DOMAIN, h, SKY, detectors.H1), project_fd(L1_DOMAIN, h, SKY, detectors.L1))
end

function finite_difference(f, x, eps)
    return map(eachindex(x)) do j
        v = zeros(length(x)); v[j] = 1
        (-f(x + 2eps * v) + 8f(x + eps * v) - 8f(x - eps * v) + f(x - 2eps * v)) / (12eps)
    end
end
relative(a, b) = norm(a - b) / max(norm(b), eps(Float64))

function run_probe(name, mode_name, destination)
    objective, x = if name == "calibration"
        calibration_loss, DRAW[1:4]
    elseif name == "hamiltonian"
        hamiltonian_loss, DRAW[1:4]
    elseif name == "flux"
        flux_loss, DRAW[1:4]
    elseif name == "timing"
        timing_loss, [2.0]
    elseif name == "initial"
        initial_loss, DRAW[1:4]
    elseif name == "evolution"
        evolution_loss, DRAW[1:4]
    elseif name == "fft"
        fft_loss, [sin(0.7j) for j in 1:16]
    elseif name == "tail"
        tail_loss, [0.03, 0.98]
    elseif name == "waveform"
        waveform_loss, copy(DRAW)
    elseif name == "likelihood"
        injected = signal(DRAW)
        loss = y -> log_likelihood_ratio(DOMAIN, signal(y), injected, METRIC)
        loss, DRAW + [0.1, -0.05, 0.001, -0.002, 1.0, 0.003, -0.002]
    else
        error("unknown probe $name")
    end
    result = Dict{String, Any}(
        "case" => name, "mode" => mode_name,
        "julia" => string(VERSION), "enzyme" => string(pkgversion(Enzyme)), "parameters" => x
    )
    try
        println(stderr, "Reference gradients: $name"); flush(stderr)
        result["primal"] = objective(x)
        expected = ForwardDiff.gradient(objective, x)
        result["forwarddiff"] = expected
        steps = collect(10.0 .^ range(-2, -5; length = 5))
        fd = [finite_difference(objective, x, eps) for eps in steps]
        result["steps"] = steps
        result["forwarddiff_fd_errors"] = [relative(expected, g) for g in fd]
        println(stderr, "Enzyme reverse: $name ($mode_name)"); flush(stderr)
        mode_name == "relaxed_aliasing" && Enzyme.API.strictAliasing!(false)
        mode = mode_name == "runtime_activity" ? Enzyme.set_runtime_activity(Enzyme.Reverse) : Enzyme.Reverse
        actual = zeros(length(x))
        elapsed = @elapsed Enzyme.autodiff(mode, Enzyme.Const(objective), Enzyme.Active, Enzyme.Duplicated(copy(x), actual))
        result["gradient"] = actual
        result["seconds_including_compilation"] = elapsed
        result["enzyme_forwarddiff_error"] = relative(actual, expected)
        result["enzyme_fd_errors"] = [relative(actual, g) for g in fd]
        result["status"] = all(isfinite, actual) && relative(actual, expected) < 1.0e-3 ? "pass" : "incorrect_gradient"
        println(stderr, result["status"], " relative error = ", result["enzyme_forwarddiff_error"])
    catch err
        result["status"] = "error"
        result["error_type"] = string(typeof(err))
        result["error"] = sprint(showerror, err; context = :limit => true)
        showerror(stderr, err, catch_backtrace()); println(stderr)
    finally
        mkpath(dirname(destination))
        open(destination, "w") do io
            TOML.print(io, result)
        end
    end
    return result["status"] == "pass" ? 0 : 1
end

exit(run_probe(ARGS...))
