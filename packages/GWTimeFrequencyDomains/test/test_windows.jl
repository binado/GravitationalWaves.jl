@testset "Planck windows" begin
    domain(n; fs = 8, epoch = 0, channels = (:H1,)) =
        TimeDomain(; epoch, nsamples = n, sampling_frequency = fs, channels)
    td = domain(17)
    w = planck_window(td; start_width = 0.5, end_width = 0.5)
    @test w isa Vector{Float64}
    @test length(w) == nsamples(td)
    @test w[1] === w[end] === 0.0
    @test w[3] === w[end - 2] === 0.5
    @test w[5:13] == ones(9)
    @test issorted(w[1:5]) && issorted(w[13:end]; rev = true)
    @test w == reverse(w)
    @test planck_window(td; start_width = 0.5, end_width = 0.5) !== w
    @test planck_window(domain(17; epoch = 1126259460.125); start_width = 0.5, end_width = 0.5) == w
    @test planck_window(domain(17; channels = (:L1, :H1)); start_width = 0.5, end_width = 0.5) == w

    asymmetric = planck_window(td; start_width = 0.5, end_width = 1)
    @test asymmetric[3] == asymmetric[13] == 0.5
    @test asymmetric[5:9] == ones(5)
    @test issorted(asymmetric[1:5]) && issorted(asymmetric[9:end]; rev = true)
    @test planck_window(td; start_width = 0, end_width = 0) == ones(17)
    @test planck_window(td; start_width = 0, end_width = 0.5)[1:13] == ones(13)
    @test planck_window(td; start_width = 0.5, end_width = 0)[5:end] == ones(13)
    @test_throws UndefKeywordError planck_window(td; start_width = 0)
    @test_throws UndefKeywordError planck_window(td; end_width = 0)

    # Continuous widths: the first interior sample is 2/3 through the transition.
    fractional = planck_window(td; start_width = 3 / 16, end_width = 0)
    @test fractional[2] ≈ inv(1 + exp(-1.5))
    @test fractional[3:end] == ones(15)
    @test planck_window(td; start_width = 1 / 32, end_width = 1 / 32) == [0; ones(15); 0]
    @test planck_window(domain(2); start_width = 1 / 32, end_width = 1 / 32) == zeros(2)
    touching = planck_window(domain(8); start_width = 7 / 16, end_width = 7 / 16)
    @test maximum(touching) < 1
    @test touching == reverse(touching)
    @test planck_window(domain(9); start_width = 0.5, end_width = 0.5)[5] == 1

    for n in (0, 1)
        @test planck_window(domain(n); start_width = 0, end_width = 0) == ones(n)
        @test_throws ArgumentError planck_window(domain(n); start_width = eps(), end_width = 0)
        @test_throws ArgumentError planck_window(domain(n); start_width = 0, end_width = eps())
    end
    for invalid in (-1, -Inf, Inf, NaN)
        @test_throws ArgumentError planck_window(td; start_width = invalid, end_width = 0)
        @test_throws ArgumentError planck_window(td; start_width = 0, end_width = invalid)
    end
    @test_throws ArgumentError planck_window(td; start_width = 1, end_width = nextfloat(1.0))
    @test_throws ArgumentError planck_window(td; start_width = duration(td), end_width = 0)
    for width in (nextfloat(0.0), 1 / 32, 3 / 16, 0.5, 1.0)
        weights = planck_window(td; start_width = width, end_width = width)
        @test all(x -> isfinite(x) && 0 <= x <= 1, weights)
    end
    # Interior exponents beyond Float64's exp range must remain finite.
    extreme = planck_window(domain(4097; fs = 4096); start_width = 0.5, end_width = 0.5)
    @test all(x -> isfinite(x) && 0 <= x <= 1, extreme)
    @test extreme[2] == extreme[end - 1] == 0
    @test planck_window(domain(17; fs = Float32(8)); start_width = Float32(0.5), end_width = 0) isa Vector{Float64}
    for big_td in (td, domain(17; fs = big"8"))
        big_w = planck_window(big_td; start_width = big"0.5", end_width = 0.5)
        @test big_w isa Vector{BigFloat}
        @test big_w[3] == big"0.5"
        @test Float64.(big_w) ≈ w
    end
    @test planck_window(domain(17; fs = big"8"); start_width = 0.5, end_width = 0) isa Vector{BigFloat}

    @testset "Broadcast application" begin
        real_data = collect(1.0:17)
        complex_data = complex.(real_data, reverse(real_data))
        matrix_data = hcat(real_data, -real_data)
        for h in (
                real_data, complex_data, matrix_data,
                view(real_data, :), view(complex_data, :), view(matrix_data, :, 1:2),
            )
            original = copy(h)
            result = w .* h
            @test size(result) == size(h)
            @test h == original
            @test all(iszero, result[1, :]) && all(iszero, result[end, :])
            for j in axes(h, 2), i in axes(h, 1)
                @test result[i, j] == w[i] * h[i, j]
            end
        end
        @test eltype(w .* Float32.(real_data)) === Float64
        @test eltype(w .* ComplexF32.(complex_data)) === ComplexF64
    end
end
