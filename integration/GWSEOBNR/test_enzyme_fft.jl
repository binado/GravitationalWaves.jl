using Test, GWSEOBNR, Enzyme, ForwardDiff, LinearAlgebra, TOML
const S = GWSEOBNR

function dft_matrix(n)
    return [cis(-2π * k * j / n) for k in 0:div(n, 2), j in 0:(n - 1)]
end

function enzyme_gradient(f, x; initial = zero(x))
    dx = copy(initial)
    Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(f), Enzyme.Active, Enzyme.Duplicated(copy(x), dx))
    return dx
end

function fd_gradient(f, x; step = 1.0e-2)
    return map(eachindex(x)) do j
        e = zeros(length(x)); e[j] = step
        (-f(x + 2e) + 8f(x + e) - 8f(x - e) + f(x - 2e)) / (12step)
    end
end

function mixed_loss(x, w)
    y = S.real_fft(x)
    value = zero(eltype(x))
    for k in eachindex(y)
        value += real(conj(w[k]) * y[k]) + oftype(value, 0.13) * abs2(y[k])
    end
    return value
end

records = Dict{String, Any}[]
@testset "local Enzyme real FFT reverse rule" begin
    @test Base.get_extension(S, :GWSEOBNREnzymeCoreExt) !== nothing
    for T in (Float64, Float32), n in (1, 2, 7, 8, 15, 16)
        x = T[sin(0.7j) + 0.2cos(1.3j) for j in 1:n]
        w = Complex{T}[complex(cos(0.3k), sin(0.9k)) for k in 1:(div(n, 2) + 1)]
        A = dft_matrix(n)
        tolerance = T === Float64 ? 2.0e-12 : 3.0e-6
        @test S.real_fft(x) ≈ A * x rtol = tolerance atol = tolerance
        before = copy(x)
        loss = y -> mixed_loss(y, w)
        actual = enzyme_gradient(loss, x)
        expected = real.(A' * (w + 0.26 .* (A * x)))
        @test x == before
        @test actual ≈ expected rtol = tolerance atol = tolerance
        @test actual ≈ ForwardDiff.gradient(loss, x) rtol = tolerance atol = tolerance
        # Five-point differences of an independent direct DFT, not FFTW.
        direct = y -> sum(real(conj(w[k]) * (A * y)[k]) + 0.13abs2((A * y)[k]) for k in eachindex(w))
        fd = fd_gradient(direct, Float64.(x))
        @test actual ≈ fd rtol = tolerance atol = tolerance
        # Existing input cotangents must accumulate, never be overwritten.
        initial = fill(T(1.7), n)
        @test enzyme_gradient(loss, x; initial) ≈ initial + expected rtol = tolerance atol = tolerance
        # Shared inputs, two FFT calls, and a direct dependence exercise accumulation.
        composed = y -> mixed_loss(y, w) + mixed_loss(2 .* y, w) + 0.21sum(abs2, y)
        @test enzyme_gradient(composed, x) ≈ ForwardDiff.gradient(composed, x) rtol = tolerance atol = tolerance
        # Individual endpoint and last-bin cotangents exercise odd/even conventions.
        for k in unique([1, div(n, 2) + 1]), component in (real, imag)
            bin = y -> component(S.real_fft(y)[k])
            target = component.(A[k, :])
            @test enzyme_gradient(bin, x) ≈ target rtol = tolerance atol = tolerance
        end
        push!(
            records, Dict(
                "type" => string(T), "length" => n,
                "direct_dft_gradient_error" => norm(actual - expected) / norm(expected),
                "finite_difference_error" => norm(actual - fd) / norm(fd)
            )
        )
    end
    # Physical normalization and time-origin phase remain active outside the rule.
    z = vcat([sin(0.7j) for j in 1:16], 0.02, -0.37)
    physical = y -> sum(real, S.physical_spectrum(y[1:16], y[17], y[18]))
    @test enzyme_gradient(physical, z) ≈ ForwardDiff.gradient(physical, z) rtol = 1.0e-11 atol = 1.0e-11
    @test enzyme_gradient(physical, z) ≈ fd_gradient(physical, z; step = 1.0e-6) rtol = 1.0e-7 atol = 1.0e-7
end
if !isempty(ARGS)
    open(ARGS[1], "w") do io
        TOML.print(io, Dict("julia" => string(VERSION), "enzyme" => string(pkgversion(Enzyme)), "cases" => records))
    end
end
