include("checks.jl")
using Enzyme, Mooncake

@testset "XAS phase systems: Enzyme + Mooncake" begin
    for case in fixtures()
        @testset "$(case.name)" begin
            loss(z) = PhenomXASSolve.objective(case.region, z, case.x, case.scale)
            z0 = perturbation(case, 0.0)
            cache = Mooncake.prepare_gradient_cache(loss, z0)
            for amount in PERTURBATIONS
                z = perturbation(case, amount)
                expected = ForwardDiff.gradient(loss, z)
                enzyme = only(Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(loss), z))
                _, tangents = Mooncake.value_and_gradient!!(cache, loss, z)
                mooncake = tangents[2]
                checkclose(enzyme, expected)
                checkclose(mooncake, expected)
            end
        end
    end
end
