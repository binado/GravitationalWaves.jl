using TestItems

@testsnippet Setup_acceptance_assembly begin
    using GWSeries
    using Test
    import AbstractWaveforms
    using AbstractWaveforms: mode_sum, sylm_m2, HigherModesWaveformGenerator, waveform_modes_fd, waveform_modes_td, hphc_fd, hphc_td
    using DifferentiationInterface
    using FiniteDifferences: central_fdm, grad
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake


    const BACKENDS = Pair{String, Any}[
        "ForwardDiff" => AutoForwardDiff(),
        "Mooncake" => AutoMooncake(; config = nothing),
    ]

    struct ADModes <: HigherModesWaveformGenerator end
    function AbstractWaveforms.waveform_modes_fd(f, p, ::ADModes)
        modes = complex.(p.x[1] .* collect(f), p.x[2] .* collect(f) .+ collect(f) .^ 2)
        return (; mode_indices = ((2, 0), (3, -1)), positive = hcat(modes, modes .* (1 + im)), negative = hcat(conj.(modes), modes .* 0.3))
    end
    function AbstractWaveforms.waveform_modes_td(::ADModes, p; kwargs...)
        return (; t = [0.0, 1.0], mode_indices = ((2, 0), (3, -1)), modes = complex.(reshape(p.x[1:4], 2, 2), reshape(p.x[2:5], 2, 2)))
    end

    # Compare the gradient of `f` at `x` from every backend with 5th-order finite differences.
    function check_gradients(f, x; backends = BACKENDS, rtol = 1.0e-6, atol = 1.0e-8)
        reference = grad(central_fdm(5, 1), f, x)[1]
        for (name, backend) in backends
            @testset "$name" begin
                g = DifferentiationInterface.gradient(f, backend, x)
                @test all(isfinite, g)
                @test isapprox(g, reference; rtol, atol)
            end
        end
        return nothing
    end
end

@testitem "shared assembly AD" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_assembly] begin
    labels = Tuple((l, m) for l in 2:5 for m in (-l):l)
    modes = [complex(sin(k + j), cos(2k - j)) for k in 1:3, j in eachindex(labels)]
    for theta in (0.0, 1.0e-10, 0.7, π - 1.0e-10, Float64(π))
        check_gradients(
            x -> begin
                h = mode_sum(modes .* complex(x[3], x[4]), labels, x[1], x[2])
                sum(real, h) + 0.3 * sum(imag, h)
            end,
            [theta, 0.4, 0.8, 0.2],
        )
    end
    for theta in (0.0, 0.7, Float64(π))
        check_gradients(
            x -> begin
                h = hphc_fd(0.0:1.0:4.0, (; x, theta_jn = x[3]), ADModes())
                sum(real, h[:, 1]) + 0.2 * sum(imag, h[:, 2])
            end,
            [0.7, 0.2, theta],
        )
        check_gradients(
            x -> begin
                h = hphc_td(ADModes(), (; x, theta_jn = x[5]); times = 0.0:1.0:1.0)
                sum(h[Chan = At(:plus)]) + 0.2 * sum(h[Chan = At(:cross)])
            end,
            [0.7, 0.2, 0.4, 0.3, theta],
        )
    end
end
