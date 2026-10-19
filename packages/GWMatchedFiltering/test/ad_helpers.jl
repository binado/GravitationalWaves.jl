const AD_F = 0.0:0.5:4.0
const AD_CHANNELS = (:ET1, :ET2, :ET3)
const AD_PSD = fill(2.0, length(AD_F), 3)
const AD_DIAGONAL_domain = FrequencyDomain(; frequencies = AD_F, epoch = 0.0, channels = AD_CHANNELS)
const AD_DIAGONAL = FrequencyDomainMetric(AD_DIAGONAL_domain, AD_PSD)
const AD_S0 = ComplexF64[2 0.2 + 0.1im 0.1; 0.2 - 0.1im 2 0.3im; 0.1 -0.3im 2]
const AD_BLOCKS_domain = FrequencyDomain(; frequencies = AD_F, epoch = 0.0, channels = AD_CHANNELS)
const AD_BLOCKS = FrequencyDomainMetric(AD_BLOCKS_domain, repeat(reshape(AD_S0, 1, 3, 3), length(AD_F), 1, 1))
const AD_DATA = reshape(complex.(sin.(1:27), cos.(1:27)), 9, 3)

ad_signal(x) = x[1] .* AD_DATA .* exp.(1.0im .* x[2] .* reshape(collect(AD_F), :, 1))
diagonal_loss(x) = log_likelihood_ratio(AD_DIAGONAL_domain, ad_signal(x), AD_DATA, AD_DIAGONAL)
correlated_loss(x) = log_likelihood_ratio(AD_BLOCKS_domain, ad_signal(x), AD_DATA, AD_BLOCKS)

function statistics_loss(x)
    h = ad_signal(x)
    rho = matched_filter_snr(AD_BLOCKS_domain, h, AD_DATA, AD_BLOCKS)
    overlap = complex_inner_product(AD_BLOCKS_domain, h, AD_DATA, AD_BLOCKS)
    return real(rho) + imag(rho) + real(overlap) + imag(overlap) +
        optimal_snr(AD_BLOCKS_domain, h, AD_BLOCKS) + residual_squared(AD_BLOCKS_domain, AD_DATA, h, AD_BLOCKS)
end
