include("checks.jl")
using Reactant, Enzyme
Reactant.set_default_backend("cpu")

function assembly_loss(region, z, x0, Abar, bbar)
    A, b = PhenomXASSolve.dense_system(region, x0 .* (1 .+ z))
    return sum(A .* Abar) + sum(b .* bbar)
end

# Explicit first-order chain rule through the solve, followed by AD of the
# assembly only. This is a prototype of the rule body, not a registered rule
# for arbitrary calls to `\` and not a higher-order derivative implementation.
function explicit_gradient(region, z, x0, cscale)
    return Reactant.@allowscalar begin
        A, b = PhenomXASSolve.dense_system(region, x0 .* (1 .+ z))
        cbar = collect(1:length(cscale)) ./ (length(cscale) .* cscale)
        _, Abar, bbar = PhenomXASSolve.solve_pullback(A, b, cbar)
        builder(znew) = assembly_loss(region, znew, x0, Abar, bbar)
        only(Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(builder), z))
    end
end

function host_pullback_gradient(region, z, x0, cscale)
    A, b = PhenomXASSolve.dense_system(region, x0 .* (1 .+ z))
    cbar = collect(1:length(cscale)) ./ (length(cscale) .* cscale)
    _, Abar, bbar = PhenomXASSolve.solve_pullback(A, b, cbar)
    return ForwardDiff.gradient(znew -> assembly_loss(region, znew, x0, Abar, bbar), z)
end

@testset "Reactant gradient using solve pullback" begin
    for case in fixtures()[1:3]
        @testset "$(case.name)" begin
            z0 = perturbation(case, 0.0)
            grad(z) = explicit_gradient(case.region, z, case.x, case.scale)
            println("Compiling explicit solve pullback: ", case.name)
            flush(stdout)
            executable = Reactant.@compile grad(Reactant.to_rarray(z0))
            draws = filter(c -> typeof(c.region) == typeof(case.region), fixtures())
            for draw in draws, amount in PERTURBATIONS
                z = (draw.x .* (1 .+ perturbation(draw, amount))) ./ case.x .- 1
                loss(z) = PhenomXASSolve.objective(case.region, z, case.x, case.scale)
                expected = ForwardDiff.gradient(loss, z)
                # Validate the explicit chain rule on the host, independently
                # of Reactant lowering, then reuse the compiled gradient.
                checkclose(host_pullback_gradient(case.region, z, case.x, case.scale), expected; rtol = 1.0e-6)
                checkclose(Array(executable(Reactant.to_rarray(z))), expected; rtol = 1.0e-6)
            end
        end
    end
end
