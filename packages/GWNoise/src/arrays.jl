# Plain-array noise covariance: preparation, whitening factors and sampling.
#
# Conventions (shared with GWMatchedFiltering):
# - Data are frequency-domain `(nf, nc)` complex arrays, one-sided from DC, in
#   `rfft(x) * dt` units. For even original length the last bin is Nyquist.
# - `S` is the one-sided PSD in strain²/Hz: `(nf, nc)` real for independent channels, or
#   `(nc, nc, nf)` Hermitian PSD per bin. `S == 0` marks an inactive bin; every array has
#   the full length `nf` (bands are zeros, never slices); DC and Nyquist must be zero.
# - Coefficients satisfy `E[n n†] = S / (2Δf)`.
# - `W` is the whitening factor `L⁻¹` with `L L† = S` per bin, zero wherever `S` is zero:
#   `(nf, nc)` real `1/√S`, or `(nc, nc, nf)` lower triangular.
# The inner product `⟨a, b⟩ = 4Δf Σₖ (W aₖ)† (W bₖ)` lives in GWMatchedFiltering.

"""
    band_psd(curve, f::AbstractRange; fmin = first(curve.frequencies), fmax = last(curve.frequencies))

One-sided PSD of `curve` on the grid `f`, as a `Vector{Float64}` of the same length.
The result is zero outside `[fmin, fmax]`, at DC and at the last bin, which is
taken to be the Nyquist bin of an even-length transform. The band must lie within the
coverage of `curve`.
"""
function band_psd(
        curve::TabulatedPSD, f::AbstractRange;
        fmin::Real = first(curve.frequencies), fmax::Real = last(curve.frequencies),
    )
    first(curve.frequencies) <= fmin <= fmax <= last(curve.frequencies) ||
        throw(ArgumentError("band [$fmin, $fmax] must lie within PSD coverage ($(first(curve.frequencies)), $(last(curve.frequencies)))"))
    length(f) >= 2 && first(f) >= 0 || throw(ArgumentError("need at least two nonnegative frequencies"))
    S = [fmin <= x <= fmax ? Float64(psd(x, curve)) : 0.0 for x in f]
    S[1] = 0.0
    S[end] = 0.0
    return S
end

check_diagonal_spectrum(S::AbstractMatrix{<:Real}) =
    all(x -> isfinite(x) && x >= 0, S) || throw(ArgumentError("PSDs must be finite and nonnegative"))

# Lower Cholesky factors `L` of the `(nc, nc, nf)` blocks; all-zero blocks give zero factors.
function cholesky_factors(S::AbstractArray{<:Number, 3})
    nc, nc2, nf = size(S)
    nc == nc2 || throw(DimensionMismatch("spectral matrices must have shape (channel, channel, frequency)"))
    T = typeof(float(zero(real(eltype(S)))))
    L = zeros(Complex{T}, nc, nc, nf)
    for k in 1:nf
        block = Matrix{Complex{T}}(S[:, :, k])
        all(isfinite, block) && ishermitian(block) ||
            throw(ArgumentError("spectral matrix at bin $k must be finite and Hermitian"))
        all(iszero, block) && continue
        factor = try
            cholesky(Hermitian(block))
        catch err
            err isa PosDefException || rethrow()
            throw(ArgumentError("spectral matrix at bin $k is not positive definite"))
        end
        L[:, :, k] = factor.L
    end
    return eltype(S) <: Real ? real.(L) : L
end

"""
    whitening_factor(S)

Whitening factor `W = L⁻¹` with `L L† = S` per bin, zero wherever `S` is zero. For
independent channels `S` is `(nf, nc)` real and `W = 1/√S`; for correlated channels `S` is
`(nc, nc, nf)` Hermitian and `W[:, :, k]` is the inverse of its lower Cholesky factor.
Inputs are validated once here (finite, non-negative, Hermitian, positive definite); a
rank-deficient nonzero block throws `ArgumentError`.
"""
function whitening_factor(S::AbstractMatrix{<:Real})
    check_diagonal_spectrum(S)
    return map(x -> iszero(x) ? zero(float(x)) : inv(sqrt(float(x))), S)
end
function whitening_factor(S::AbstractArray{<:Number, 3})
    L = cholesky_factors(S)
    W = similar(L)
    for k in axes(L, 3)
        nc = size(L, 1)
        W[:, :, k] = iszero(L[:, :, k]) ? zero(L[:, :, k]) : L[:, :, k] \ Matrix{eltype(L)}(I, nc, nc)
    end
    return W
end

"""
    sample_noise([rng], S, df) -> Matrix{ComplexF64}
    sample_noise!([rng], out, S, df)

Draw frequency-domain Gaussian noise with coefficient covariance `E[n n†] = S / (2df)` in
`rfft(strain) * dt` units: `n = L z / (2√df)` with `z = randn + i·randn` per bin, where
`L L† = S`. `S` is `(nf, nc)` real or `(nc, nc, nf)` Hermitian, and `out` is `(nf, nc)`.
Random numbers are drawn for every bin, active or not, so a stream stays aligned across
bands. The block Cholesky factorization is recomputed per call.
"""
sample_noise(S::AbstractArray, df::Real) = sample_noise(default_rng(), S, df)
function sample_noise(rng::AbstractRNG, S::AbstractArray, df::Real)
    out = Matrix{ComplexF64}(undef, size(S, ndims(S) == 3 ? 3 : 1), ndims(S) == 3 ? size(S, 1) : size(S, 2))
    return sample_noise!(rng, out, S, df)
end
sample_noise!(out::AbstractMatrix{<:Complex}, S::AbstractArray, df::Real) = sample_noise!(default_rng(), out, S, df)
function sample_noise!(rng::AbstractRNG, out::AbstractMatrix{<:Complex}, S::AbstractMatrix{<:Real}, df::Real)
    size(out) == size(S) || throw(DimensionMismatch("noise output has size $(size(out)), expected $(size(S))"))
    df > 0 || throw(ArgumentError("df must be positive"))
    check_diagonal_spectrum(S)
    for channel in axes(S, 2), k in axes(S, 1)
        scale = sqrt(S[k, channel] / df) / 2
        out[k, channel] = scale * complex(randn(rng), randn(rng))
    end
    return out
end
function sample_noise!(rng::AbstractRNG, out::AbstractMatrix{<:Complex}, S::AbstractArray{<:Number, 3}, df::Real)
    nc, _, nf = size(S)
    size(out) == (nf, nc) || throw(DimensionMismatch("noise output has size $(size(out)), expected ($nf, $nc)"))
    df > 0 || throw(ArgumentError("df must be positive"))
    L = cholesky_factors(S)
    scale = inv(2 * sqrt(df))
    z = Vector{ComplexF64}(undef, nc)
    for k in 1:nf
        for c in 1:nc
            z[c] = complex(randn(rng), randn(rng))
        end
        for i in 1:nc
            acc = zero(ComplexF64)
            for j in 1:i
                acc += L[i, j, k] * z[j]
            end
            out[k, i] = scale * acc
        end
    end
    return out
end
