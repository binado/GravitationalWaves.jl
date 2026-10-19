module GWSEOBNREnzymeCoreExt

using GWSEOBNR: real_fft, check_evolution_retcode, check_initial_solution
using FFTW: plan_rfft, ESTIMATE
using EnzymeCore: Const, Duplicated
using EnzymeCore: EnzymeRules

# Discrete status validation has no numerical output or differentiable side
# effects. Inactivity preserves execution (including errors) in the primal.
# Do not mark the solve, trajectory, or continuous termination times inactive.
EnzymeRules.inactive(::typeof(check_evolution_retcode), retcode) = true
EnzymeRules.inactive(::typeof(check_initial_solution), sol) = true

# The constant-plan/adjoint construction follows the discussion in
# https://github.com/SciML/NonlinearSolve.jl/issues/597 . This rule is restricted
# to our out-of-place, one-dimensional real FFT, with no global FFTW annotations.
const FFTScalar = Union{Float32, Float64}

function EnzymeRules.augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, ::Const{typeof(real_fft)}, ::Type,
        x::Duplicated{Vector{T}},
    ) where {T <: FFTScalar}
    # Planning is fixed configuration. Use scratch storage even though ESTIMATE
    # preserves its input; never expose active waveform samples to the planner.
    plan = plan_rfft(zeros(T, length(x.val)); flags = ESTIMATE)
    y = plan * x.val
    shadow = EnzymeRules.needs_shadow(config) ? zero(y) : nothing
    primal = EnzymeRules.needs_primal(config) ? y : nothing
    return EnzymeRules.AugmentedReturn(primal, shadow, (plan, shadow))
end

function EnzymeRules.reverse(
        ::EnzymeRules.RevConfigWidth{1}, ::Const{typeof(real_fft)}, ::Type,
        tape, x::Duplicated{Vector{T}},
    ) where {T <: FFTScalar}
    plan, shadow = tape
    if shadow !== nothing
        # AbstractFFTs handles normalization and the real half-spectrum endpoints.
        x.dval .+= adjoint(plan) * shadow
        fill!(shadow, zero(eltype(shadow)))
    end
    return (nothing,)
end

end
