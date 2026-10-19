"""
    CoorbitalComponent(node_indices, fits, basis)

One datapiece of the coorbital waveform surrogate: scalar fits at the empirical-interpolation
nodes `node_indices` (1-based positions on the coorbital time grid), contracted with the
empirical-interpolation `basis` of size `(n_times, n_nodes)`.
"""
struct CoorbitalComponent
    node_indices::Vector{Int}
    fits::FitSet
    basis::Matrix{Float64}
end

"""
    CoorbitalSurrogate(params, settings, t, masks, ell_max, components, starts)

Coorbital-frame waveform surrogate. `components` lists the datapieces in the order set by
`starts[ell - 1, m + 1]`, the index of the first datapiece of mode `(ell, m)` with `m >= 0`.
Each mode contributes 2 parts (`m = 0`: real, imaginary) or 4 parts (`m > 0`: `Re+`, `Re-`,
`Im+`, `Im-`), each split in `size(masks, 2)` subdomains that are blended with `masks`
(`n_times × n_subdomains`; a single all-ones column for NRSur7dq4).

The component layout and its integer indexing follow the AD design rules: no `Dict`, plain
loops over the flat model data.
"""
struct CoorbitalSurrogate{P}
    params::P
    settings::FitSettings
    t::Vector{Float64}
    masks::Matrix{Float64}
    ell_max::Int
    components::Vector{CoorbitalComponent}
    starts::Matrix{Int}
end

n_subdomains(s::CoorbitalSurrogate) = size(s.masks, 2)

"Index of the `part` (1-based) of mode `(ell, m)`, `m >= 0`, in subdomain `sub` in `s.components`."
component_index(s::CoorbitalSurrogate, ell, m, part, sub) = s.starts[ell - 1, m + 1] + (part - 1) * n_subdomains(s) + (sub - 1)

"Column of the mode `(ell, m)` in the output of [`coorbital_modes`](@ref)."
mode_column(ell, m) = ell * (ell + 1) - 3 + m

"Number of modes `(ell, m)` with `2 <= ell <= ell_max`."
n_modes(ell_max) = ell_max * ell_max + 2 * ell_max - 3

function eval_component!(
        out::AbstractVector, s::CoorbitalSurrogate, c::CoorbitalComponent, q, chiA, chiB, powers, nodes,
    )
    for (n, ni) in enumerate(c.node_indices)
        x = fit_params(s.params, q, (chiA[ni, 1], chiA[ni, 2], chiA[ni, 3]), (chiB[ni, 1], chiB[ni, 2], chiB[ni, 3]))
        fit_powers!(powers, s.settings, x)
        nodes[n] = eval_fit(c.fits, n, powers)
    end
    fill!(out, zero(eltype(out)))
    for n in eachindex(c.node_indices)
        w = nodes[n]
        for i in eachindex(out)
            out[i] += w * c.basis[i, n]
        end
    end
    return out
end

"""
    coorbital_modes(surrogate, q, chiA, chiB; ell_max = surrogate.ell_max)

Coorbital-frame modes `h_{ell m}(t)` on `surrogate.t` for mass ratio `q` and the time-dependent
coorbital-frame spins `chiA`, `chiB` (matrices of size `(n_times, 3)`). Returns a complex
`(n_times, n_modes)` matrix with columns ordered `(2, -2), …, (2, 2), (3, -3), …`, the same
order as `gwsurrogate`; modes absent from the model are zero.
"""
function coorbital_modes(s::CoorbitalSurrogate, q, chiA::AbstractMatrix, chiB::AbstractMatrix; ell_max::Integer = s.ell_max)
    2 <= ell_max <= s.ell_max || throw(ArgumentError("ell_max must be in 2:$(s.ell_max)"))
    size(chiA) == (length(s.t), 3) == size(chiB) || throw(DimensionMismatch("spins must have size (length(t), 3)"))
    T = promote_type(typeof(q), eltype(chiA), eltype(chiB), Float64)
    nt = length(s.t)
    nsub = n_subdomains(s)
    modes = zeros(Complex{T}, nt, n_modes(ell_max))
    powers = Vector{T}(undef, n_fit_powers(s.settings))
    nodes = Vector{T}(undef, maximum(c -> length(c.node_indices), s.components))
    parts = [Vector{T}(undef, nt) for _ in 1:4]
    for ell in 2:ell_max
        for m in 0:ell
            nparts = m == 0 ? 2 : 4
            for sub in 1:nsub
                for part in 1:nparts
                    c = s.components[component_index(s, ell, m, part, sub)]
                    eval_component!(parts[part], s, c, q, chiA, chiB, powers, nodes)
                end
                mask = view(s.masks, :, sub)
                if m == 0
                    col = mode_column(ell, 0)
                    for i in 1:nt
                        modes[i, col] += complex(parts[1][i], parts[2][i]) * mask[i]
                    end
                else
                    cp = mode_column(ell, m)
                    cn = mode_column(ell, -m)
                    for i in 1:nt
                        rep, rem, imp, imm = parts[1][i], parts[2][i], parts[3][i], parts[4][i]
                        modes[i, cp] += complex(rep - rem, imm - imp) * mask[i]
                        modes[i, cn] += complex(rep + rem, imp + imm) * mask[i]
                    end
                end
            end
        end
    end
    return modes
end
