# Prototype of XAS phase collocation solves, not a waveform implementation.
# Formulas follow Ripple revision 58c53024718697e0430b666ef5dd40ff3c2e029d.
# GPL-2.0-or-later, as for GWWaveforms and the underlying LAL model.
module PhenomXASSolve

using LinearAlgebra
using StaticArrays

const G4 = SVector(0.0, 0.25, 0.75, 1.0)
const G5 = SVector(0.0, 0.5 - 1 / (2 * sqrt(2.0)), 0.5, 0.5 + 1 / (2 * sqrt(2.0)), 1.0)
const DPHASE0 = 5.0 / (128.0 * Float64(π)^(5 / 3))

# x = [fMECO, fISCO, fRD, fDAMP, rhs...]. The RHS is already the
# assembled collocation RHS, including any Lorentzian subtraction.
function points(::Val{:inspiral}, x)
    xmin = 0.0026
    xmax = 1.02 * x[1]
    return SVector(ntuple(i -> G4[i] * (xmax - xmin) + xmin, Val(4)))
end

function points(::Val{:intermediate}, x)
    fmatch = 0.6 * (0.5 * x[3] + x[2])
    delta = (fmatch - x[1]) * 0.03
    xmin = x[1] - delta
    xmax = fmatch + 0.5 * delta
    return SVector(ntuple(i -> G5[i] * (xmax - xmin) + xmin, Val(5)))
end

function points(::Val{:ringdown}, x)
    xmin = 0.6 * (0.5 * x[3] + x[2])
    xmax = x[3] + 1.25 * x[4]
    p = SVector(ntuple(i -> G5[i] * (xmax - xmin) + xmin, Val(5)))
    return SVector(p[1], p[2], p[3], x[3], p[5])
end

basis(::Val{:inspiral}, f, j, x) = j == 1 ? one(f) : f^((j - 1) / 3)
basis(::Val{:intermediate}, f, j, x) = (x[3] / f)^(j - 1)
function basis(::Val{:ringdown}, f, j, x)
    j == 1 && return one(f)
    j == 2 && return f^(-1 / 3)
    j == 3 && return f^-2
    j == 4 && return f^-4
    return -DPHASE0 / (x[4]^2 + (f - x[3])^2)
end

function system(region, x)
    p = points(region, x)
    return assemble(region, x, p)
end

function assemble(region, x, p::StaticVector{N}) where {N}
    A = SMatrix{N, N}(ntuple(k -> basis(region, p[mod1(k, N)], cld(k, N), x), Val(N * N)))
    b = SVector{N}(ntuple(i -> x[4 + i], Val(N)))
    return A, b
end

static_solve(A, b) = A \ b
dense_solve(A, b) = Matrix(A) \ Vector(b)
native_dense_solve(A, b) = lu(A) \ b

# Real-valued solve pullback. A second factorization avoids depending on a
# backend's transpose-factorization wrapper; a production rule can reuse LU.
# These are ordinary solves: AD need not differentiate their factorizations.
function solve_pullback(A, b, cbar)
    c = native_dense_solve(A, b)
    # Match the RHS storage to b: a host vector of cotangents otherwise becomes
    # a vector of traced scalars instead of a traced array under Reactant.
    lambda = native_dense_solve(transpose(A), zero.(b) .+ cbar)
    return c, -(lambda * transpose(c)), lambda
end

# Assemble array columns directly instead of converting an SMatrix of traced
# scalars. This lets Reactant use its native dense linear-algebra lowering.
function dense_points(::Val{:inspiral}, x)
    return collect(G4) .* (1.02 * x[1] - 0.0026) .+ 0.0026
end

function dense_points(::Val{:intermediate}, x)
    fmatch = 0.6 * (0.5 * x[3] + x[2])
    delta = (fmatch - x[1]) * 0.03
    xmin = x[1] - delta
    return collect(G5) .* (fmatch + 0.5 * delta - xmin) .+ xmin
end

function dense_points(::Val{:ringdown}, x)
    xmin = 0.6 * (0.5 * x[3] + x[2])
    p = collect(G5) .* (x[3] + 1.25 * x[4] - xmin) .+ xmin
    p[4] = x[3]
    return p
end

dense_basis(::Val{:inspiral}, p, j, x) = p .^ ((j - 1) / 3)
dense_basis(::Val{:intermediate}, p, j, x) = (x[3] ./ p) .^ (j - 1)
function dense_basis(::Val{:ringdown}, p, j, x)
    j == 1 && return one.(p)
    j == 2 && return p .^ (-1 / 3)
    j == 3 && return p .^ -2
    j == 4 && return p .^ -4
    return -DPHASE0 ./ (x[4]^2 .+ (p .- x[3]) .^ 2)
end

function dense_system(region, x)
    p = dense_points(region, x)
    A = hcat(ntuple(j -> dense_basis(region, p, j, x), length(p))...)
    return A, x[5:(4 + length(p))]
end

function coefficients(region, x, ::typeof(native_dense_solve))
    A, b = dense_system(region, x)
    return native_dense_solve(A, b)
end

function coefficients(region, x, solver::F = static_solve) where {F}
    A, b = system(region, x)
    return solver(A, b)
end

# Dimensionless perturbations keep finite-difference steps and gradient scales
# comparable despite frequencies ~0.01 and some RHS entries ~1000.
function normalized_coefficients(region, z, x0, cscale, solver::F = static_solve) where {F}
    x = x0 .* (1 .+ z)
    return coefficients(region, x, solver) ./ cscale
end

# Linear loss gives independent sensitivity to every coefficient. We also test
# full Jacobians, so cancellation in this scalar cannot hide a bad derivative.
function objective(region, z, x0, cscale, solver::F = static_solve) where {F}
    c = normalized_coefficients(region, z, x0, cscale, solver)
    return sum(i -> (i / length(c)) * c[i], eachindex(c))
end

end
