"""Detached sampling domains, explicit validation, and physical Fourier conversion."""
module GWTimeFrequencyDomains

export TimeDomain, FrequencyDomain, epoch, channels, nsamples, times, frequencies
export sample_spacing, sampling_frequency, duration, frequency_spacing
export check_data, check_same_grid, check_same_frequency_grid, check_compatible, check_fourier_domain
export crop, pad, frequency_domain, time_domain, to_frequencydomain, to_timedomain
export plan_to_frequencydomain, plan_to_timedomain, output_domain, transform, transform!

include("domain_common.jl")
include("time_domain.jl")
include("frequency_domain.jl")
include("domain_operations.jl")
include("transform_interface.jl")

end
