# Run each mode in a fresh process; full mode patches only that process's container.
using GWSEOBNR, SciMLSensitivity, Enzyme, ForwardDiff, LinearAlgebra, TOML, Test
import SciMLStructures as SS
include("scimlstructures_parameters.jl")
using .StructuredParametersPrototype: EOBParameters
const S = GWSEOBNR
const X = [40.0, 20.0, 0.4, -0.3]
const U = [12.0, 0.3, -0.005, 4.0]
parameters(x) = let c = S.hamiltonian_coefficients(x...)
    EOBParameters(c, S.flux_coefficients(c))
end
function rhs_loss(x)
    d = S.eob_rhs(U, parameters(x), 0.0)
    return d[1] + 2d[2] + 3d[3] + 4d[4]
end
function short_loss(x, structured = true)
    p = parameters(x)
    pars = structured ? p : (; p.c, p.coeffs)
    omega = π * (x[1] + x[2]) * S.MTSUN_SI * 16
    u0 = S.initial_conditions(x..., omega)
    sol = S.solve(
        S.ODEProblem(S.eob_rhs, u0, (0.0, 20.0), pars), S.Vern9();
        reltol = 1.0e-11, abstol = 1.0e-12, dtmax = 10.0, saveat = [20.0]
    )
    u = sol.u[end]
    return u[1] + 2u[2] + 3u[3] + 4u[4]
end
function fd(f, x, h)
    return map(eachindex(x)) do j
        v = zeros(length(x)); v[j] = h
        (-f(x + 2v) + 8f(x + v) - 8f(x - v) + f(x - 2v)) / (12h)
    end
end
relative(a, b) = norm(a - b) / max(norm(b), eps())
function record_gradient(f, x)
    result = Dict{String, Any}("primal" => f(x), "forwarddiff" => ForwardDiff.gradient(f, x))
    steps = 10.0 .^ range(-2, -5; length = 5)
    finite = [fd(f, x, h) for h in steps]
    result["steps"] = collect(steps)
    result["forwarddiff_fd_errors"] = [relative(result["forwarddiff"], g) for g in finite]
    try
        g = zeros(length(x))
        elapsed = @elapsed Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(f), Enzyme.Active, Enzyme.Duplicated(copy(x), g))
        result["gradient"] = g
        result["seconds_including_compilation"] = elapsed
        result["enzyme_forwarddiff_error"] = relative(g, result["forwarddiff"])
        result["enzyme_fd_errors"] = [relative(g, d) for d in finite]
        result["status"] = all(isfinite, g) && result["enzyme_forwarddiff_error"] < 1.0e-3 ? "pass" : "incorrect_gradient"
    catch err
        result["status"] = "error"
        result["error_type"] = string(typeof(err))
        result["error"] = sprint(showerror, err; context = :limit => true)
        showerror(stderr, err, catch_backtrace()); println(stderr)
    end
    return result
end

mode, destination = ARGS
if mode == "waveform" || mode == "evolution"
    StructuredParametersPrototype.install!()
    empty!(ARGS); append!(ARGS, [mode, "default", destination])
    include("probe_enzyme.jl")
else
    result = Dict{String, Any}(
        "julia" => string(VERSION), "enzyme" => string(pkgversion(Enzyme)),
        "scimlsensitivity" => string(pkgversion(SciMLSensitivity)), "scimlstructures" => string(pkgversion(SS)),
        "tunable_count" => StructuredParametersPrototype.COUNT, "container_assertions" => 102
    )
    @testset "SciMLStructures container" begin
        for M in (40.0, 60.0, 100.0), (q, a, b) in ((1.0, 0.0, 0.0), (2.0, 0.4, -0.3), (4.0, -0.6, 0.6))
            x = [M * q / (1 + q), M / (1 + q), a, b]
            p = parameters(x)
            v, repack, aliases = SS.canonicalize(SS.Tunable(), p)
            @test length(v) == result["tunable_count"]
            @test !aliases
            @test SS.isscimlstructure(p) && !SS.ismutablescimlstructure(p)
            @test !SS.hasportion(SS.Constants(), p)
            @test !SS.hasportion(SS.Caches(), p)
            @test SS.canonicalize(SS.Constants(), p) === (nothing, nothing, nothing)
            restored = repack(v)
            @test restored.c == p.c && restored.coeffs == p.coeffs
            @test SS.canonicalize(SS.Tunable(), SS.replace(SS.Tunable(), p, v .+ 1))[1] == v .+ 1
            @test SS.canonicalize(SS.Tunable(), p)[1] == v
            @test S.eob_rhs(U, p, 0.0) == S.eob_rhs(U, (; p.c, p.coeffs), 0.0)
            @test ForwardDiff.gradient(rhs_loss, x) ≈ ForwardDiff.gradient(
                x -> begin
                    p = parameters(x)
                    d = S.eob_rhs(U, (; p.c, p.coeffs), 0.0)
                    d[1] + 2d[2] + 3d[3] + 4d[4]
                end, x
            ) rtol = 1.0e-13
        end
        @test short_loss(X) == short_loss(X, false)
        @test ForwardDiff.gradient(short_loss, X) ≈ ForwardDiff.gradient(x -> short_loss(x, false), X) rtol = 1.0e-12
        @test_throws DimensionMismatch SS.replace(SS.Tunable(), parameters(X), zeros(2))
    end
    println(stderr, "Enzyme RHS"); flush(stderr)
    result["rhs"] = record_gradient(rhs_loss, X)
    println(stderr, "Enzyme short fixed-time ODE"); flush(stderr)
    result["short_ode"] = record_gradient(short_loss, X)
    mkpath(dirname(destination))
    open(destination, "w") do io
        TOML.print(io, result)
    end
    println("Results saved: ", destination)
    exit(all(result[k]["status"] == "pass" for k in ("rhs", "short_ode")) ? 0 : 1)
end
