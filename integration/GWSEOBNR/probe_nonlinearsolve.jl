# Run as a separate process for each sensitivity/activity configuration.
include("nonlinearsolve_initial.jl")
using .NonlinearInitialPrototype
using NonlinearSolve, SciMLSensitivity, Enzyme, ForwardDiff, LinearAlgebra, TOML
length(ARGS) == 3 || error("usage: probe_nonlinearsolve.jl default|dense default|runtime OUTPUT.toml")
ARGS[1] in ("default", "dense") || error("unknown sensitivity configuration")
ARGS[2] in ("default", "runtime") || error("unknown activity configuration")
const P = NonlinearInitialPrototype
const SENSEALG = ARGS[1] == "dense" ? SteadyStateAdjoint(; autojacvec = false) : nothing
coordinate_loss(x, j) = P.initial_conditions(x; sensealg = SENSEALG)[j]

function run_probe(activity, destination)
    result = Dict{String, Any}(
        "julia" => string(VERSION),
        "nonlinearsolve" => string(pkgversion(NonlinearSolve)),
        "scimlsensitivity" => string(pkgversion(SciMLSensitivity)),
        "enzyme" => string(pkgversion(Enzyme)), "sensitivity" => ARGS[1],
        "activity" => activity, "cases" => Dict{String, Any}[]
    )
    save() = open(io -> TOML.print(io, result), destination, "w")
    mkpath(dirname(destination))
    mode = activity == "runtime" ? Enzyme.set_runtime_activity(Enzyme.Reverse) : Enzyme.Reverse
    for M in (40.0, 60.0, 100.0), (q, chi1, chi2) in ((1.0, 0.0, 0.0), (2.0, 0.4, -0.3), (4.0, -0.6, 0.6))
        x = [M * q / (1 + q), M / (1 + q), chi1, chi2]
        case = Dict{String, Any}("parameters" => x, "total_mass" => M, "mass_ratio" => q, "status" => "running", "stage" => "primal")
        push!(result["cases"], case); save()
        println("M=$M q=$q: primal/ForwardDiff/FD"); flush(stdout)
        try
            sol = P.initial_solution(x; sensealg = SENSEALG)
            P.check_initial_solution(sol)
            u = P.initial_conditions(x; sensealg = SENSEALG)
            expected = ForwardDiff.jacobian(P.reference_initial, x)
            jac = ForwardDiff.jacobian(y -> P.initial_conditions(y; sensealg = SENSEALG), x)
            case["state"] = u
            reference_state = P.reference_initial(x)
            case["state_relative_error"] = P.relative(u, reference_state)
            case["state_component_errors"] = [P.relative(u[j], reference_state[j]) for j in (1, 3, 4)]
            case["state_components_match"] = all(isapprox.(u, reference_state; rtol = 1.0e-9, atol = 1.0e-13))
            case["residual_norm"] = norm(sol.resid, Inf)
            case["retcode"] = string(sol.retcode)
            case["forwarddiff_reference_error"] = P.relative(jac, expected)
            steps = collect(10.0 .^ range(-2, -5; length = 5))
            fd = [reduce(hcat, P.fivepoint(P.reference_initial, x, h)) for h in steps]
            case["fd_steps"] = steps
            case["forwarddiff_fd_errors"] = [P.relative(jac, g) for g in fd]
            case["stage"] = "enzyme"; save()
            println("M=$M q=$q: native Enzyme reverse"); flush(stdout)
            actual = zeros(4, 4)
            elapsed = @elapsed for j in (1, 3, 4)
                row = zeros(4)
                Enzyme.autodiff(
                    mode, Enzyme.Const(coordinate_loss), Enzyme.Active,
                    Enzyme.Duplicated(copy(x), row), Enzyme.Const(j)
                )
                actual[j, :] .= row
            end
            case["forwarddiff_jacobian_rows"] = collect.(eachrow(expected))
            case["enzyme_jacobian_rows"] = collect.(eachrow(actual))
            case["enzyme_seconds_including_compilation"] = elapsed
            case["enzyme_forwarddiff_error"] = P.relative(actual, expected)
            case["enzyme_row_errors"] = [P.relative(actual[j, :], expected[j, :]) for j in (1, 3, 4)]
            case["enzyme_fd_errors"] = [P.relative(actual, g) for g in fd]
            case["enzyme_fd_row_errors"] = [[P.relative(actual[j, :], g[j, :]) for g in fd] for j in (1, 3, 4)]
            adjacent = all(any(errors[k] < 1.0e-3 && errors[k + 1] < 1.0e-3 for k in 1:4) for errors in case["enzyme_fd_row_errors"])
            pass = case["state_components_match"] && case["state_relative_error"] < 1.0e-9 && case["residual_norm"] < 1.0e-12 &&
                case["forwarddiff_reference_error"] < 1.0e-3 && maximum(case["enzyme_row_errors"]) < 1.0e-3 && adjacent
            case["status"] = pass ? "pass" : "incorrect_gradient"
            case["stage"] = "complete"
        catch err
            case["status"] = "error"
            case["error_type"] = string(typeof(err))
            case["error"] = sprint(showerror, err; context = :limit => true)
            showerror(stdout, err, catch_backtrace());println();flush(stdout)
        end
        save()
        # An exception fails this configuration. Preserve the partial record;
        # remaining draws are untested, not counted as passing.
        case["status"] == "error" && return 1
    end
    x = [40.0, 20.0, 0.4, -0.3]
    actual = zeros(4)
    Enzyme.autodiff(mode, Enzyme.Const(P.initial_loss), Enzyme.Active, Enzyme.Duplicated(copy(x), actual))
    result["module_default_gradient_error"] = P.relative(actual, ForwardDiff.gradient(P.initial_loss, x))
    failed = P.initial_solution(x; maxiters = 0)
    result["maxiters_retcode"] = string(failed.retcode)
    result["maxiters_rejected"] = false
    try
        P.check_initial_solution(failed)
    catch err
        result["maxiters_rejected"] = err isa P.S.SolveError && !P.S.SciMLBase.successful_retcode(failed)
    end
    save()
    return all(c["status"] == "pass" for c in result["cases"]) &&
        result["module_default_gradient_error"] < 1.0e-3 && result["maxiters_rejected"] ? 0 : 1
end
exit(run_probe(ARGS[2], ARGS[3]))
