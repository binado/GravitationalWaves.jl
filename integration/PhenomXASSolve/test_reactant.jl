include("checks.jl")
using Reactant, Enzyme
Reactant.set_default_backend("cpu")

# Probe dense storage, with optional static paths; values and AD are separate.
# Compilation failures are reported and recorded, not silently counted as passes.
function reactant_probe(case, solver::F; values_only = false) where {F}
    z0 = perturbation(case, 0.0)
    tracked = Reactant.to_rarray(z0)
    fun(z) = Reactant.@allowscalar PhenomXASSolve.normalized_coefficients(case.region, z, case.x, case.scale, solver)
    loss(z) = Reactant.@allowscalar PhenomXASSolve.objective(case.region, z, case.x, case.scale, solver)
    println("Compiling values: ", case.name, " / ", solver)
    flush(stdout)
    value_exec = Reactant.@compile fun(tracked)
    draws = filter(c -> typeof(c.region) == typeof(case.region), fixtures())
    for draw in draws, amount in PERTURBATIONS
        # Change frequencies/RHS across all physical fixture draws using the
        # same executable. Intermediate systems include different LU pivots.
        z = (draw.x .* (1 .+ perturbation(draw, amount))) ./ case.x .- 1
        actual = Array(value_exec(Reactant.to_rarray(z)))
        checkclose(actual, fun(z); rtol = 1.0e-7)
    end
    values_only && return nothing
    println("Compiling gradient: ", case.name, " / ", solver)
    flush(stdout)
    gradient_exec = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, loss, tracked)
    for draw in draws, amount in PERTURBATIONS
        z = (draw.x .* (1 .+ perturbation(draw, amount))) ./ case.x .- 1
        gradient = Array(only(gradient_exec(Enzyme.Reverse, loss, Reactant.to_rarray(z))))
        checkclose(gradient, ForwardDiff.gradient(loss, z); rtol = 1.0e-6)
    end
    return nothing
end

failures = String[]
values_only = "--values-only" in ARGS
solvers = "--probe-static" in ARGS ?
    (PhenomXASSolve.native_dense_solve, PhenomXASSolve.static_solve, PhenomXASSolve.dense_solve) :
    (PhenomXASSolve.native_dense_solve,)
@testset "XAS phase systems: Reactant CPU probes" begin
    for solver in solvers
        # One compiled executable per region, then multiple parameter draws.
        for case in fixtures()[1:3]
            try
                reactant_probe(case, solver; values_only)
            catch err
                label = "$(case.name) / $solver"
                push!(failures, label)
                println(stderr, "UNSUPPORTED: ", label)
                showerror(stderr, err)
                println(stderr)
            end
        end
    end
end
println("Reactant unsupported paths: ", failures)
values_only && println("Gradient compilation was not requested (--values-only).")
isempty(failures) || exit(1)
