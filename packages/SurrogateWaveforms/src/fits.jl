"""
    FitSettings(q_offset, q_slope, q_max_order, chi_max_order)

Model constants of the parametric fits: the mass-ratio parameter is rescaled to `[-1, 1]` with
`q_offset + q_slope * log(q)` and the basis-function orders are bounded by `q_max_order` and
`chi_max_order`.
"""
struct FitSettings
    q_offset::Float64
    q_slope::Float64
    q_max_order::Int
    chi_max_order::Int
end

"Fit settings shared by NRSur7dq4 and NRSur7dq4v2."
const NRSUR7DQ4_FIT_SETTINGS = FitSettings(-0.9857019407834238, 1.4298059216576398, 3, 2)

"Number of entries of the power table used by [`fit_powers!`](@ref)."
n_fit_powers(s::FitSettings) = s.q_max_order + 1 + 6 * (s.chi_max_order + 1)

# Table slot of the power `order` of fit parameter `j` (1 for q, 2:7 for the spin parameters).
fit_slot(s::FitSettings, j::Integer, order::Integer) =
    j == 1 ? order + 1 : s.q_max_order + 1 + (s.chi_max_order + 1) * (j - 2) + order + 1

"""
    FitSet(settings, orders, coefs)

A collection of parametric fits stored in flat, integer-indexed vectors (CSR style), so model
data can be read with plain loops: fit `k` owns the terms `offsets[k]+1 : offsets[k+1]`. Each
fit is a sum of `coef * prod_j x_j^n_j` over the seven fit parameters; `slots[:, t]` holds, for
term `t`, the positions of the seven powers in the table filled by [`fit_powers!`](@ref).

`orders` is a vector of `(7, n_terms)` integer matrices and `coefs` a vector of coefficient
vectors, one entry per fit.
"""
struct FitSet
    slots::Matrix{Int}
    coefs::Vector{Float64}
    offsets::Vector{Int}
end

function FitSet(settings::FitSettings, orders::AbstractVector, coefs::AbstractVector)
    length(orders) == length(coefs) || throw(DimensionMismatch("one orders matrix per coefficient vector"))
    n_terms = sum(length, coefs; init = 0)
    slots = Matrix{Int}(undef, 7, n_terms)
    flat = Vector{Float64}(undef, n_terms)
    offsets = zeros(Int, length(coefs) + 1)
    t = 0
    for (k, (o, c)) in enumerate(zip(orders, coefs))
        size(o) == (7, length(c)) || throw(DimensionMismatch("orders must have size (7, length(coefs))"))
        for i in eachindex(c)
            t += 1
            flat[t] = c[i]
            for j in 1:7
                n = o[j, i]
                max_order = j == 1 ? settings.q_max_order : settings.chi_max_order
                0 <= n <= max_order || throw(ArgumentError("basis function order $n out of range for parameter $j"))
                slots[j, t] = fit_slot(settings, j, n)
            end
        end
        offsets[k + 1] = t
    end
    return FitSet(slots, flat, offsets)
end

Base.length(fits::FitSet) = length(fits.offsets) - 1

"""
    fit_powers!(powers, settings, x)

Fill the buffer `powers` (length [`n_fit_powers`](@ref)) with the powers of the seven fit
parameters `x` that the fit terms index: powers `0:q_max_order` of the rescaled `x[1]`, then
powers `0:chi_max_order` of each of `x[2:7]`. Powers are built by successive multiplication, as
in `gwsurrogate`.
"""
function fit_powers!(powers::AbstractVector, s::FitSettings, x)
    qbase = s.q_offset + s.q_slope * x[1]
    powers[1] = one(qbase)
    for i in 1:s.q_max_order
        powers[i + 1] = powers[i] * qbase
    end
    for j in 2:7
        base = s.q_max_order + 1 + (s.chi_max_order + 1) * (j - 2)
        powers[base + 1] = one(qbase)
        for i in 1:s.chi_max_order
            powers[base + i + 1] = powers[base + i] * x[j]
        end
    end
    return powers
end

"""
    eval_fit(fits, k, powers)

Evaluate fit `k` of `fits` with the power table `powers` filled by [`fit_powers!`](@ref).
"""
function eval_fit(fits::FitSet, k::Integer, powers::AbstractVector)
    acc = zero(eltype(powers))
    for t in (fits.offsets[k] + 1):fits.offsets[k + 1]
        prod = powers[fits.slots[1, t]]
        for j in 2:7
            prod *= powers[fits.slots[j, t]]
        end
        acc += fits.coefs[t] * prod
    end
    return acc
end

"""
    NRSur7dq4Params
    NRSur7dq4v2Params

Fit-parameter conventions. [`fit_params`](@ref) maps the mass ratio and the two spin vectors to
the seven parameters the fits of that model were trained on.
"""
struct NRSur7dq4Params end
struct NRSur7dq4v2Params end

@inline function chi_hat(q, chi1z, chi2z)
    eta = q / (1 + q)^2
    chi_wt_avg = (q * chi1z + chi2z) / (1 + q)
    return (chi_wt_avg - 38 * eta / 113 * (chi1z + chi2z)) / (1 - 76 * eta / 113)
end

"""
    fit_params(params, q, chiA, chiB) -> NTuple{7}

Fit parameters for mass ratio `q` and spin 3-vectors `chiA`, `chiB`. NRSur7dq4 uses
`(log q, χAx, χAy, χ̂, χBx, χBy, (χAz - χBz)/2)`; NRSur7dq4v2 uses
`(log q, (χAx + χBx)/2, (χAx - χBx)/2, χ̂, (χAy + χBy)/2, (χAy - χBy)/2, χAz)`, where `χ̂` is the
effective spin of Eq. (3) of arXiv:1508.07253.
"""
function fit_params(::NRSur7dq4Params, q, chiA, chiB)
    return (log(q), chiA[1], chiA[2], chi_hat(q, chiA[3], chiB[3]), chiB[1], chiB[2], (chiA[3] - chiB[3]) / 2)
end

function fit_params(::NRSur7dq4v2Params, q, chiA, chiB)
    return (
        log(q), (chiA[1] + chiB[1]) / 2, (chiA[1] - chiB[1]) / 2, chi_hat(q, chiA[3], chiB[3]),
        (chiA[2] + chiB[2]) / 2, (chiA[2] - chiB[2]) / 2, chiA[3],
    )
end
