"""Complex loggamma with a local analytic first-order ForwardDiff adapter."""
tail_loggamma(z::Complex{<:AbstractFloat}) = loggamma(z)

function tail_loggamma(z::Complex{ForwardDiff.Dual{Tag, V, P}}) where {Tag, V <: AbstractFloat, P}
    z0 = complex(ForwardDiff.value(real(z)), ForwardDiff.value(imag(z)))
    value = loggamma(z0)
    slope = digamma(z0)
    dz = ntuple(P) do j
        slope * complex(ForwardDiff.partials(real(z))[j], ForwardDiff.partials(imag(z))[j])
    end
    return complex(
        ForwardDiff.Dual{Tag}(real(value), ntuple(j -> real(dz[j]), P)),
        ForwardDiff.Dual{Tag}(imag(value), ntuple(j -> imag(dz[j]), P)),
    )
end

function complex_tail(l, m, omega, H)
    k = m * omega
    khat = H * k
    return exp(
        tail_loggamma(complex(l + 1 + zero(khat), -2khat)) + π * khat +
            2im * khat * log(4k / sqrt(exp(1.0)))
    ) / factorial(l)
end
