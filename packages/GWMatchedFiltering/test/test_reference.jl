using JSON

@testset "bilby and NumPy references" begin
    path = joinpath(@__DIR__, "..", "..", "..", "reference", "inference.json")
    if isfile(path)
        reference = JSON.parsefile(path)
        ascomplex(pair) = complex(pair...)
        df = 0.5
        h = ComplexF64[0 0; 1 + 2im 2 - im; 2im -1; -1 im; 0.5 2]
        d = ComplexF64[0 0; 2 - im -1 + 3im; im 2; 1 + 2im -im; 1 -0.5im]
        psds = Float64[0 0; 2 3; 3 4; 4 5; 5 6]
        W = whitening_diagonal(psds)
        @test complex_inner_product(h, d, W, df) ≈ ascomplex(reference["diagonal_overlap"]) rtol = 1.0e-14
        @test optimal_snr_squared(h, W, df) ≈ reference["diagonal_norm"] rtol = 1.0e-14
        @test matched_filter_snr(h[:, 1], d[:, 1], W[:, 1], df) ≈ ascomplex(reference["single_matched"]) rtol = 1.0e-14
        S = zeros(ComplexF64, 2, 2, 5)
        S[1, 1, :] = psds[:, 1]
        S[2, 2, :] = psds[:, 2]
        S[1, 2, 2:end] .= 0.2 + 0.3im
        S[2, 1, 2:end] .= 0.2 - 0.3im
        @test complex_inner_product(h, d, whitening_blocks(S), df) ≈ ascomplex(reference["correlated_overlap"]) rtol = 1.0e-14
    end
end
