"""Joint frequency-domain noise metrics, inner products, and SNR statistics."""
module GWMatchedFiltering

import GWSeries
using GWSeries: FrequencySeries, check_same_grid, check_same_frequency_grid, frequencies, channels
import GWSeries: check_compatible, frequency_spacing

using LinearAlgebra: Hermitian, I, cholesky, ishermitian
using StaticArrays: SMatrix, SVector

export AbstractMetric, DiagonalMetric, BlockDiagonalMetric
export FrequencyDomainMetric
export complex_inner_product, noise_weighted_inner_product
export optimal_snr_squared, optimal_snr, matched_filter_snr, residual_squared
export likelihood_statistics, log_likelihood_ratio, check_template

include("types.jl")
include("metric.jl")
include("statistics.jl")

end
