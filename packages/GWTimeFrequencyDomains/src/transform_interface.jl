"""Physical dt*rfft of validated real CPU samples along dimension one. Load FFTW. No domain checks."""
function to_frequencydomain end
"""
    to_timedomain(domain, h)

Physical inverse FFT of validated CPU bins along dimension one. Load FFTW.
A source `FrequencyDomain` infers the even output length `2*(nsamples(domain)-1)`.
A target `TimeDomain` supplies the output length and sampling frequency explicitly,
including odd lengths. Bins must cover its complete zero-origin one-sided FFT grid.
Imaginary DC and even-length Nyquist parts are discarded. No domain checks.
"""
function to_timedomain end
"""Prepare a forward FFT on the host, validating the domain and prototype."""
function plan_to_frequencydomain end
"""Prepare an inverse FFT on the host, validating the domain and prototype."""
function plan_to_timedomain end
"""The prepared transform's immutable output domain."""
function output_domain end
"""Allocate and execute a prepared physical transform without domain validation."""
function transform end
"""Execute into separate destination storage; inputs must match the prepared prototype."""
function transform! end
