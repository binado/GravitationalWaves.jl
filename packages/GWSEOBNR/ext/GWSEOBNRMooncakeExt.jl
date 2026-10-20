module GWSEOBNRMooncakeExt

using GWSEOBNR: SEOBNRv5HM, real_fft, spectrum_vector, stopgrad, check_evolution_retcode, check_initial_solution, check_nqc_solution
using FFTW: plan_rfft, ESTIMATE
using ForwardDiff: ForwardDiff
using Mooncake: Mooncake
using ChainRulesCore: ChainRulesCore, NoTangent

# `stopgrad` is the identity on non-Dual reals, so only an explicit rule makes the
# implicit-function correction in `continuous_root` cancel to a pure first-order term.
Mooncake.@zero_adjoint Mooncake.DefaultCtx Tuple{typeof(stopgrad), Real}

# Discrete status validation has no differentiable output. Marking it zero-adjoint keeps
# the primal checks (and their errors) while skipping `cond`'s internal try/catch.
Mooncake.@zero_adjoint Mooncake.DefaultCtx Tuple{typeof(check_evolution_retcode), Any}
Mooncake.@zero_adjoint Mooncake.DefaultCtx Tuple{typeof(check_initial_solution), Any}
Mooncake.@zero_adjoint Mooncake.DefaultCtx Tuple{typeof(check_nqc_solution), Any, Any}

# Planning is fixed configuration; only the linear transform is differentiated.
# The pullback of a real-input FFT is the real part of the adjoint plan applied to the
# half-spectrum cotangent (AbstractFFTs handles normalization and the endpoints).
function ChainRulesCore.rrule(::typeof(real_fft), x::Vector{Float64})
    plan = plan_rfft(zeros(Float64, length(x)); flags = ESTIMATE)
    y = plan * x
    pullback(ȳ) = (NoTangent(), real(adjoint(plan) * collect(ChainRulesCore.unthunk(ȳ))))
    return y, pullback
end

Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(real_fft), Vector{Float64}}


# The ODE solve, its callbacks and the implicit initial conditions are not differentiated
# in reverse. The map from the 7 real parameters to the spectrum has a small Jacobian, so
# one forward-mode pass gives it exactly and the pullback is a pair of matrix products.
function real_spectrum(f, x, g)
    y = spectrum_vector(f, x, g)
    return vcat(real.(y), imag.(y))
end

function ChainRulesCore.rrule(::typeof(spectrum_vector), f::AbstractRange, x::Vector{Float64}, g::SEOBNRv5HM)
    n = 2 * length(f)
    result = ForwardDiff.DiffResults.JacobianResult(zeros(2n), x)
    ForwardDiff.jacobian!(result, z -> real_spectrum(f, z, g), x)
    value = ForwardDiff.DiffResults.value(result)
    J = ForwardDiff.DiffResults.jacobian(result)
    y = complex.(value[1:n], value[(n + 1):(2n)])
    function pullback(ȳ)
        ȳ = collect(ChainRulesCore.unthunk(ȳ))
        x̄ = J[1:n, :]' * real.(ȳ) + J[(n + 1):(2n), :]' * imag.(ȳ)
        return NoTangent(), NoTangent(), x̄, NoTangent()
    end
    return y, pullback
end

Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(spectrum_vector), AbstractRange, Vector{Float64}, SEOBNRv5HM}

end
