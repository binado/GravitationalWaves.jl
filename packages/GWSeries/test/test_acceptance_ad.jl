using TestItems

@testsnippet Setup_acceptance_ad begin
    using Test, GWSeries, FFTW, ForwardDiff, Mooncake, ChainRulesCore, DifferentiationInterface
end

@testitem "Physical Fourier AD" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_ad] begin
    for n in (7, 8), nc in (1, 2), dt in (0.125, 1 / 2048)
        a = reshape(sin.(1.0:(n * nc)), n, nc)
        saved = copy(a)
        meta = (; epoch = -0.2, sampling_interval = dt, channels = ntuple(i -> Symbol(:c, i), nc))
        fourier(x) = parent(to_frequencydomain(TimeSeries(x; meta...)))
        dft = [dt * cis(-2π * k * j / n) for k in 0:(n ÷ 2), j in 0:(n - 1)]
        @test fourier(a) ≈ dft * a
        for row in (1, 2, n ÷ 2 + 1)
            weights = zeros(ComplexF64, n ÷ 2 + 1, nc)
            weights[row, :] .= 0.7 + 0.3im
            loss(x) = sum(real, conj.(weights) .* fourier(x))
            expected = real.(adjoint(dft) * weights)
            for backend in (AutoForwardDiff(), AutoMooncake(; config = nothing))
                prep = prepare_gradient(loss, backend, a)
                for x in (a, a .* 1.01)
                    value, g = value_and_gradient(loss, prep, backend, x)
                    @test value ≈ loss(x) rtol = 1.0e-12 atol = 1.0e-14
                    @test g ≈ expected rtol = 1.0e-11 atol = 1.0e-14
                end
            end
        end
        @test a == saved
    end
end

@testitem "Mooncake through unrecorded spectrum metadata" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_ad] begin
    # Reading the grid records of a spectrum without a recorded length once crashed Mooncake's
    # rule compilation (a `Nothing`-typed field nested in the dims); keep this path compiled.
    spectrum(x) = FrequencySeries(reshape([x[1], x[2], x[1] * x[2]], :, 1); epoch = 0.0, frequency_spacing = 1.0, channels = (:H1,))
    loss(x) = (s = spectrum(x); sum(abs2, parent(s)) + something(time_samples(s), 0) + epoch(s))
    @test time_samples(spectrum([1.0, 2.0])) === nothing
    @test gradient(loss, AutoMooncake(; config = nothing), [1.5, 2.0]) ≈ gradient(loss, AutoForwardDiff(), [1.5, 2.0])
end
