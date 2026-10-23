# Test-only construction of `W` from `S`; GWNoise.whitening_factor is the production path.
whitening_diagonal(S::AbstractMatrix{<:Real}) = ifelse.(S .> 0, inv.(sqrt.(S)), zero(eltype(S)))
function whitening_blocks(S::AbstractArray{<:Number, 3})
    nc, _, nf = size(S)
    T = typeof(float(zero(eltype(S))))
    W = zeros(T, nc, nc, nf)
    for k in 1:nf
        block = Matrix{T}(S[:, :, k])
        all(iszero, block) && continue
        W[:, :, k] = cholesky(Hermitian(block)).L \ Matrix{T}(I, nc, nc)
    end
    return W
end
