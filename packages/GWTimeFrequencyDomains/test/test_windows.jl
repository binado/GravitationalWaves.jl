@testset "Planck windows" begin
    domain(n; fs = 8, epoch = 0, channels = (:H1,)) =
        TimeDomain(; epoch, nsamples = n, sampling_frequency = fs, channels)
    weights(td; kwargs...) = PlanckWindow(td; kwargs...).weights
    td = domain(17)
    window = PlanckWindow(td; start_width = 0.5, end_width = 0.5)
    @test window isa PlanckWindow
    @test nsamples(window) == nsamples(td)
    @test step(window) == step(td)
    @test check_compatible(td, window) === nothing
    @test check_compatible(domain(17; epoch = 1126259460.125, channels = (:L1, :H1)), window) === nothing
    @test check_compatible(domain(17; fs = nextfloat(8.0)), window) === nothing
    @test_throws DimensionMismatch check_compatible(domain(16), window)
    @test_throws ArgumentError check_compatible(domain(17; fs = 16), window)
    w = weights(td; start_width = 0.5, end_width = 0.5)
    @test w isa Vector{Float64}
    @test length(w) == nsamples(td)
    @test w[1] === w[end] === 0.0
    @test w[3] === w[end - 2] === 0.5
    @test w[5:13] == ones(9)
    @test issorted(w[1:5]) && issorted(w[13:end]; rev = true)
    @test w == reverse(w)
    @test weights(td; start_width = 0.5, end_width = 0.5) !== w
    @test weights(domain(17; epoch = 1126259460.125); start_width = 0.5, end_width = 0.5) == w
    @test weights(domain(17; channels = (:L1, :H1)); start_width = 0.5, end_width = 0.5) == w

    asymmetric = weights(td; start_width = 0.5, end_width = 1)
    @test asymmetric[3] == asymmetric[13] == 0.5
    @test asymmetric[5:9] == ones(5)
    @test issorted(asymmetric[1:5]) && issorted(asymmetric[9:end]; rev = true)
    @test weights(td; start_width = 0, end_width = 0) == ones(17)
    @test weights(td; start_width = 0, end_width = 0.5)[1:13] == ones(13)
    @test weights(td; start_width = 0.5, end_width = 0)[5:end] == ones(13)
    @test_throws UndefKeywordError weights(td; start_width = 0)
    @test_throws UndefKeywordError weights(td; end_width = 0)

    # Continuous widths: the first interior sample is 2/3 through the transition.
    fractional = weights(td; start_width = 3 / 16, end_width = 0)
    @test fractional[2] ≈ inv(1 + exp(-1.5))
    @test fractional[3:end] == ones(15)
    @test weights(td; start_width = 1 / 32, end_width = 1 / 32) == [0; ones(15); 0]
    @test weights(domain(2); start_width = 1 / 32, end_width = 1 / 32) == zeros(2)
    touching = weights(domain(8); start_width = 7 / 16, end_width = 7 / 16)
    @test maximum(touching) < 1
    @test touching == reverse(touching)
    @test weights(domain(9); start_width = 0.5, end_width = 0.5)[5] == 1

    for n in (0, 1)
        @test weights(domain(n); start_width = 0, end_width = 0) == ones(n)
        @test_throws ArgumentError weights(domain(n); start_width = eps(), end_width = 0)
        @test_throws ArgumentError weights(domain(n); start_width = 0, end_width = eps())
    end
    for invalid in (-1, -Inf, Inf, NaN)
        @test_throws ArgumentError weights(td; start_width = invalid, end_width = 0)
        @test_throws ArgumentError weights(td; start_width = 0, end_width = invalid)
    end
    @test_throws ArgumentError weights(td; start_width = 1, end_width = nextfloat(1.0))
    @test_throws ArgumentError weights(td; start_width = duration(td), end_width = 0)
    for width in (nextfloat(0.0), 1 / 32, 3 / 16, 0.5, 1.0)
        samples = weights(td; start_width = width, end_width = width)
        @test all(x -> isfinite(x) && 0 <= x <= 1, samples)
    end
    # Interior exponents beyond Float64's exp range must remain finite.
    extreme = weights(domain(4097; fs = 4096); start_width = 0.5, end_width = 0.5)
    @test all(x -> isfinite(x) && 0 <= x <= 1, extreme)
    @test extreme[2] == extreme[end - 1] == 0
    @test weights(domain(17; fs = Float32(8)); start_width = Float32(0.5), end_width = 0) isa Vector{Float64}
    for big_td in (td, domain(17; fs = big"8"))
        big_w = weights(big_td; start_width = big"0.5", end_width = 0.5)
        @test big_w isa Vector{BigFloat}
        @test big_w[3] == big"0.5"
        @test Float64.(big_w) ≈ w
    end
    @test weights(domain(17; fs = big"8"); start_width = 0.5, end_width = 0) isa Vector{BigFloat}

    @testset "Window application" begin
        real_data = collect(1.0:17)
        complex_data = complex.(real_data, reverse(real_data))
        matrix_data = hcat(real_data, -real_data)
        tensor_data = reshape(collect(1.0:102), 17, 3, 2)
        for h in (
                real_data, complex_data, matrix_data, tensor_data,
                view(repeat(tensor_data; inner = (2, 1, 1)), 1:2:34, :, :),
                view(real_data, :), view(complex_data, :), view(matrix_data, :, 1:2),
                view(repeat(real_data; inner = 2), 1:2:34),
                view(repeat(matrix_data; inner = (2, 1)), 1:2:34, :),
            )
            original = copy(h)
            result = apply_window(window, h)
            @test size(result) == size(h)
            @test h == original
            out = similar(h)
            @test apply_window!(out, window, h) === out
            @test out == result
            @test h == original
            @test apply_window!(h, window) === h
            @test h == result
            h .= original
            @test apply_window!(h, window, h) === h
            @test h == result
            h .= original
            @test window.weights == w
            @test all(iszero, selectdim(result, 1, 1)) && all(iszero, selectdim(result, 1, length(w)))
            @test result == w .* original
        end
        @test eltype(apply_window(window, Float32.(real_data))) === Float64
        @test eltype(apply_window(window, ComplexF32.(complex_data))) === ComplexF64
    end
end
