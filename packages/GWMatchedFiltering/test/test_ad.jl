using TestItems

@testsnippet Setup_ad begin
    include("ad_helpers.jl")
end

@testitem "Fixed-whitening strain derivatives" tags = [:ad] setup = [Common, Setup_ad] begin
    x = [0.8, 0.2]
    for loss in (diagonal_loss, correlated_loss, statistics_loss)
        actual = ForwardDiff.gradient(loss, x)
        delta = 1.0e-5
        expected = [(loss(x + delta * e) - loss(x - delta * e)) / (2delta) for e in eachcol(Matrix{Float64}(I, 2, 2))]
        @test actual ≈ expected rtol = 1.0e-8 atol = 1.0e-8
    end
end
