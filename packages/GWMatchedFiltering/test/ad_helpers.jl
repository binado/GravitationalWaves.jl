const AD_F = 0.0:0.5:4.0
const AD_DF = step(AD_F)
const AD_PSD = fill(2.0, length(AD_F), 3)
const AD_DIAGONAL = whitening_diagonal(AD_PSD)
const AD_S0 = ComplexF64[2 0.2 + 0.1im 0.1; 0.2 - 0.1im 2 0.3im; 0.1 -0.3im 2]
const AD_BLOCKS = whitening_blocks(repeat(AD_S0, 1, 1, length(AD_F)))
const AD_DATA = reshape(complex.(sin.(1:27), cos.(1:27)), 9, 3)

ad_signal(x) = x[1] .* AD_DATA .* exp.(1.0im .* x[2] .* reshape(collect(AD_F), :, 1))
diagonal_loss(x) = log_likelihood_ratio(ad_signal(x), AD_DATA, AD_DIAGONAL, AD_DF)
correlated_loss(x) = log_likelihood_ratio(ad_signal(x), AD_DATA, AD_BLOCKS, AD_DF)

function statistics_loss(x)
    h = ad_signal(x)
    rho = matched_filter_snr(h, AD_DATA, AD_BLOCKS, AD_DF)
    overlap = complex_inner_product(h, AD_DATA, AD_BLOCKS, AD_DF)
    return real(rho) + imag(rho) + real(overlap) + imag(overlap) +
        optimal_snr(h, AD_BLOCKS, AD_DF) + residual_squared(AD_DATA, h, AD_BLOCKS, AD_DF)
end
