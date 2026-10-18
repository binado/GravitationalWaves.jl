using JSON

@testset "bilby and NumPy references" begin
    path = joinpath(@__DIR__, "..", "..", "..", "reference", "inference.json")
    if !isfile(path)
        if get(ENV, "CI", "false") == "true"
            error("Missing inference fixtures; run python/generate_inference.py")
        end
        @warn "Skipping inference references; run python/generate_inference.py"
    else
        reference = JSON.parsefile(path)
        ascomplex(pair) = complex(pair...)
        f = 0.0:0.5:2.0
        h = ComplexF64[0 0; 1 + 2im 2 - im; 2im -1; -1 im; 0.5 2]
        d = ComplexF64[0 0; 2 - im -1 + 3im; im 2; 1 + 2im -im; 1 -0.5im]
        psds = Float64[0 0; 2 3; 3 4; 4 5; 5 6]
        metric = FrequencyDomainMetric(f, (:A, :B), psds)
        @test complex_inner_product(h, d, metric) ≈ ascomplex(reference["diagonal_overlap"]) rtol = 1.0e-14
        @test optimal_snr_squared(h, metric) ≈ reference["diagonal_norm"] rtol = 1.0e-14
        single = FrequencyDomainMetric(f, (:A,), psds[:, 1])
        @test matched_filter_snr(h[:, 1], d[:, 1], single) ≈ ascomplex(reference["single_matched"]) rtol = 1.0e-14
        S = zeros(ComplexF64, 5, 2, 2)
        S[:, 1, 1] = psds[:, 1]
        S[:, 2, 2] = psds[:, 2]
        S[2:end, 1, 2] .= 0.2 + 0.3im
        S[2:end, 2, 1] .= 0.2 - 0.3im
        correlated = FrequencyDomainMetric(f, (:A, :B), S)
        @test complex_inner_product(h, d, correlated) ≈ ascomplex(reference["correlated_overlap"]) rtol = 1.0e-14
    end
end
