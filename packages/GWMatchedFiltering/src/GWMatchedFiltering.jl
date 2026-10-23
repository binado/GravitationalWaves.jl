"""Whitened frequency-domain inner products, SNRs and Gaussian likelihood statistics on plain arrays."""
module GWMatchedFiltering

export complex_inner_product, noise_weighted_inner_product
export optimal_snr_squared, optimal_snr, matched_filter_snr, residual_squared
export likelihood_statistics, log_likelihood_ratio, check_template
export whiten, whiten!

include("whitened.jl")

end
