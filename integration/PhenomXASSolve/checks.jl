using Test, LinearAlgebra, StaticArrays, ForwardDiff, FiniteDifferences, TOML
include("collocation.jl")
using .PhenomXASSolve

function fixtures()
    data = TOML.parsefile(joinpath(@__DIR__, "fixtures.toml"))
    return map(data["case"]) do c
        x = vcat(c["frequencies"], c["rhs"])
        scale = max.(abs.(c["solution"]), 1.0e-12)
        return (;
            name = c["name"], region = Val(Symbol(c["region"])), x, scale,
            matrix = permutedims(hcat(c["matrix"]...)), solution = c["solution"],
        )
    end
end

function implicit_jacobian(case, z)
    x = case.x .* (1 .+ z)
    A, b = PhenomXASSolve.system(case.region, x)
    n = length(b)
    # Independent high-precision solve, with derivatives of the matrix builder
    # only. This checks that AD through LU includes the dA contribution.
    c = Matrix{BigFloat}(A) \ Vector{BigFloat}(b)
    dA = ForwardDiff.jacobian(z -> vec(first(PhenomXASSolve.system(case.region, case.x .* (1 .+ z)))), z)
    db = ForwardDiff.jacobian(z -> last(PhenomXASSolve.system(case.region, case.x .* (1 .+ z))), z)
    J = hcat(
        [
            Matrix{BigFloat}(A) \ (BigFloat.(db[:, j]) - reshape(BigFloat.(dA[:, j]), n, n) * c)
                for j in eachindex(z)
        ]...
    )
    return Float64.(J) ./ case.scale
end

function checkclose(actual, expected; rtol = 1.0e-7, atol = 1.0e-9)
    @test all(isfinite, actual)
    @test norm(actual - expected, Inf) <= atol + rtol * norm(expected, Inf)
    return nothing
end

function check_primal(case)
    A, b = PhenomXASSolve.system(case.region, case.x)
    checkclose(Matrix(A), case.matrix; rtol = 1.0e-13, atol = 1.0e-13)
    c = PhenomXASSolve.static_solve(A, b)
    dense = PhenomXASSolve.coefficients(case.region, case.x, PhenomXASSolve.native_dense_solve)
    checkclose(dense ./ case.scale, c ./ case.scale; rtol = 1.0e-8)
    bigc = Matrix{BigFloat}(A) \ Vector{BigFloat}(b)
    # Check coefficients componentwise after scaling, not just largest sample.
    checkclose(c ./ case.scale, case.solution ./ case.scale; rtol = 1.0e-8)
    checkclose(c ./ case.scale, Float64.(bigc) ./ case.scale; rtol = 1.0e-8)
    backward_error = norm(A * c - b, Inf) / (norm(A, Inf) * norm(c, Inf) + norm(b, Inf))
    @test backward_error < 1.0e-14
    return (;
        condition = cond(Matrix(A)), backward_error,
        coefficient_error = maximum(abs, (c .- Float64.(bigc)) ./ case.scale),
    )
end

const PERTURBATIONS = (0.0, 0.01)
perturbation(case, amount) = amount .* sin.(collect(eachindex(case.x)))

function forward_checks()
    rows = []
    @testset "XAS phase systems: StaticArrays + ForwardDiff" begin
        for case in fixtures()
            @testset "$(case.name)" begin
                stats = check_primal(case)
                maxerror = 0.0
                for amount in PERTURBATIONS
                    z = perturbation(case, amount)
                    fun(z) = PhenomXASSolve.normalized_coefficients(case.region, z, case.x, case.scale)
                    J = ForwardDiff.jacobian(fun, z)
                    expected = implicit_jacobian(case, z)
                    checkclose(J, expected)
                    densefun(z) = PhenomXASSolve.normalized_coefficients(case.region, z, case.x, case.scale, PhenomXASSolve.native_dense_solve)
                    checkclose(ForwardDiff.jacobian(densefun, z), expected)
                    @test norm(expected[:, 1:4]) > 0 # matrix dependence is exercised
                    Jfd = only(FiniteDifferences.jacobian(central_fdm(5, 1), fun, z))
                    checkclose(J, Jfd; rtol = 3.0e-5, atol = 1.0e-6)
                    maxerror = max(maxerror, norm(J - expected, Inf) / max(1.0, norm(expected, Inf)))
                end
                push!(rows, (; name = case.name, stats..., jacobian_error = maxerror))
            end
        end
    end
    println("Maximum matrix condition number: ", maximum(r -> r.condition, rows))
    println("Maximum normalized coefficient error vs BigFloat: ", maximum(r -> r.coefficient_error, rows))
    println("Maximum relative Jacobian error vs implicit equation: ", maximum(r -> r.jacobian_error, rows))
    println("Maximum normwise backward error: ", maximum(r -> r.backward_error, rows))
    return rows
end
