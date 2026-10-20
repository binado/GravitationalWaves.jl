@testset "Fixed Planck window ForwardDiff and Mooncake" begin
    td = TimeDomain(; epoch = 1126259460.125, nsamples = 17, sampling_frequency = 8, channels = (:H1,))
    w = planck_window(td; start_width = 0.5, end_width = 0.75)
    loss(h) = sum(abs2, w .* h)
    for h in (collect(range(-2.0, 3.0; length = 17)), reshape(collect(1.0:34), 17, 2))
        expected_value = sum(w[i]^2 * h[i, j]^2 for i in axes(h, 1), j in axes(h, 2))
        expected_gradient = 2 .* w .^ 2 .* h
        for backend in (AutoForwardDiff(), AutoMooncake(; config = nothing))
            value, gradient = DifferentiationInterface.value_and_gradient(loss, backend, h)
            @test value ≈ expected_value
            @test gradient ≈ expected_gradient
            @test size(gradient) == size(h)
            @test all(iszero, gradient[1, :]) && all(iszero, gradient[end, :])
        end
    end
end
