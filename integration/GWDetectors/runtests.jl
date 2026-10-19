using Test
using GWTimeFrequencyDomains
using GWDetectors
using ForwardDiff
using DifferentiationInterface
using Enzyme: Enzyme
using Mooncake: Mooncake
using Reactant

Reactant.set_default_backend("cpu")
include("../../packages/GWDetectors/test/ad_helpers.jl")

@testset "GWDetectors AD integration" begin
    samples = (AD_INPUT, AD_INPUT .+ [0.1, 0.02, -0.03, 0.125, 0.25, 0.01, -0.02, 0.04, -0.01, 0.001, 0.002, 10.0])
    @testset "Enzyme and Mooncake" begin
        for loss in (geometry_loss, antenna_loss, timing_loss, projection_loss)
            for backend in (AutoEnzyme(; function_annotation = Enzyme.Const), AutoMooncake(; config = nothing))
                for sample in samples
                    g = DifferentiationInterface.gradient(loss, backend, sample)
                    @test g ≈ ForwardDiff.gradient(loss, sample) rtol = 1.0e-9 atol = 1.0e-10
                end
            end
        end
    end
    @testset "Reactant CPU values and compiled Enzyme gradients" begin
        for loss in (geometry_loss, antenna_loss, timing_loss, projection_loss)
            # Scalar indexing is explicitly permitted; no materialization to host.
            objective(x) = Reactant.@allowscalar loss(x)
            tracked = Reactant.to_rarray(first(samples))
            executable = Reactant.@compile objective(tracked)
            gradient = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
            for sample in samples
                args = Reactant.to_rarray(sample)
                @test Reactant.to_number(executable(args)) ≈ loss(sample) rtol = 2.0e-10 atol = 2.0e-10
                g = Array(only(gradient(Enzyme.Reverse, objective, args)))
                @test g ≈ ForwardDiff.gradient(loss, sample) rtol = 2.0e-9 atol = 2.0e-9
            end
        end
    end
    include("test_reactant.jl")
end
