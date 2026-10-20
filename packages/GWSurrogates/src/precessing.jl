"""
    PrecessingSurrogate(dynamics, coorbital)

A complete precessing surrogate: the [`DynamicsSurrogate`](@ref) that evolves the coprecessing
frame and the spins, and the [`CoorbitalSurrogate`](@ref) that gives the waveform modes in the
coorbital frame. [`inertial_modes`](@ref) combines them into the inertial-frame modes.
"""
struct PrecessingSurrogate{P}
    dynamics::DynamicsSurrogate{P}
    coorbital::CoorbitalSurrogate{P}
end

"Largest `ell` of the model."
ell_max(s::PrecessingSurrogate) = s.coorbital.ell_max

"""
    load_surrogate(model::Symbol; download = true) -> PrecessingSurrogate

Load `:NRSur7dq4` or `:NRSur7dq4v2` (see [`data_path`](@ref) for where the data comes from).
"""
function load_surrogate(model::Symbol; download::Bool = true)
    return load_surrogate(data_path(model; download), MODEL_FILES[model].params)
end

"""
    load_surrogate(path, params) -> PrecessingSurrogate

Load the surrogate from the HDF5 file `path` with the fit-parameter convention `params`
([`NRSur7dq4Params`](@ref) or [`NRSur7dq4v2Params`](@ref)).
"""
function load_surrogate(path::AbstractString, params)
    return PrecessingSurrogate(load_dynamics(path, params), load_coorbital(path, params))
end

# Hamilton product of quaternions (w, x, y, z), `multiplyQuats` in gwsurrogate.
@inline function quat_multiply(a, b)
    return (
        a[1] * b[1] - a[2] * b[2] - a[3] * b[3] - a[4] * b[4],
        a[3] * b[4] - b[3] * a[4] + a[1] * b[2] + b[1] * a[2],
        a[4] * b[2] - b[4] * a[2] + a[1] * b[3] + b[1] * a[3],
        a[2] * b[3] - b[2] * a[3] + a[1] * b[4] + b[1] * a[4],
    )
end

# Natural cubic spline of every column of `Y` (knots `t_in`), evaluated at `t_out`.
function spline_columns(t_out, t_in, Y::AbstractMatrix)
    issorted(t_in; lt = <=) || throw(ArgumentError("knots must be strictly increasing"))
    T = promote_type(eltype(Y), typeof(float(first(t_out))))
    out = Matrix{T}(undef, length(t_out), size(Y, 2))
    for j in axes(Y, 2)
        u = Y[:, j]
        m = spline_second_derivatives(t_in, u)
        for (i, x) in enumerate(t_out)
            out[i, j] = spline_eval(t_in, u, m, x)
        end
    end
    return out
end

# Rotate the coorbital-frame modes `h` (n_times × n_modes) into the inertial frame:
# h_inertial[ell, m′] = Σ_m D^ell_{m′ m}(R) h[ell, m] with the rotor R of each time.
function rotate_modes(h::AbstractMatrix, quat, ell_max::Integer)
    T = promote_type(eltype(h), eltype(quat))
    out = zeros(T, size(h))
    powers = SpinorPowers(complex(eltype(quat)), 2ell_max)
    for k in axes(h, 1)
        ra = complex(quat[k, 1], quat[k, 4])
        rb = complex(quat[k, 3], quat[k, 2])
        fill_powers!(powers, ra, rb)
        for l in 2:ell_max, mp in (-l):l
            acc = zero(T)
            for m in (-l):l
                acc += wigner_D_element(l, mp, m, powers) * h[k, mode_column(l, m)]
            end
            out[k, mode_column(l, mp)] = acc
        end
    end
    return out
end

"""
    inertial_modes(s, q, chiA0, chiB0; kwargs...) -> (t, modes)

Inertial-frame waveform modes `h_{ℓm}(t)` of the precessing surrogate in dimensionless units
(`G = c = M = 1`, `M` the total mass). `chiA0`, `chiB0` are the spin 3-vectors of the heavier and
lighter black hole at the reference epoch in the LAL source frame (`z` along the orbital angular
momentum, `x` from the lighter to the heavier hole). `modes` is complex `(n_times, n_modes)` with
columns `(2, -2), …, (2, 2), (3, -3), …` up to `ell_max`.

The reference epoch is set either by `t_ref` (preferred when differentiating) or by the
reference orbital frequency `f_ref` (in cycles/M, the orbital frequency `ω = π f_ref`), see
[`t_ref_from_omega`](@ref); without either the first node is used. The waveform starts at the time
where `ω = π f_low` when `f_low > 0`.

The output grid is the sparse model grid by default, a uniform grid of step `dt` starting at the
start time, or the explicit `times` (which must not start before the `f_low` start time). Specify at most one of `dt`, `times`; `init_orbphase` and
`init_quat` give the orbital phase and the coprecessing-to-inertial rotation at the reference epoch.
"""
function inertial_modes(
        s::PrecessingSurrogate, q, chiA0, chiB0;
        ell_max::Integer = GWSurrogates.ell_max(s), f_ref = nothing, f_low = nothing, t_ref = nothing,
        dt = nothing, times = nothing, init_orbphase = 0.0, init_quat = (1.0, 0.0, 0.0, 0.0),
    )
    2 <= ell_max <= GWSurrogates.ell_max(s) || throw(ArgumentError("ell_max must be in 2:$(GWSurrogates.ell_max(s))"))
    (dt === nothing || times === nothing) || throw(ArgumentError("specify at most one of dt and times"))
    use_f_ref = f_ref !== nothing && f_ref != 0
    (t_ref === nothing || !use_f_ref) || throw(ArgumentError("specify at most one of t_ref and f_ref"))
    dyn = s.dynamics
    coorb = s.coorbital
    t_coorb = coorb.t

    chiA_norm = sqrt(chiA0[1]^2 + chiA0[2]^2 + chiA0[3]^2)
    chiB_norm = sqrt(chiB0[1]^2 + chiB0[2]^2 + chiB0[3]^2)

    omega_ref = use_f_ref ? π * f_ref : nothing
    if omega_ref !== nothing
        t_ref = t_ref_from_omega(dyn, omega_ref, q, chiA0, chiB0; init_orbphase, init_quat)
    end
    t_low = nothing
    if f_low !== nothing && f_low != 0
        omega_low = π * f_low
        t_low = if omega_ref !== nothing && abs(omega_low - omega_ref) < 1.0e-10
            t_ref
        else
            t_ref_from_omega(dyn, omega_low, q, chiA0, chiB0; init_orbphase, init_quat)
        end
    end

    d = dynamics(dyn, q, chiA0, chiB0; init_quat, init_orbphase, t_ref)

    # Interpolate to the coorbital time grid, where the coorbital spins oscillate faster than the
    # coprecessing ones, and rebuild the frame quantities.
    chiA_copr = spline_columns(t_coorb, dyn.tds, d.chiA)
    chiB_copr = spline_columns(t_coorb, dyn.tds, d.chiB)
    n = length(t_coorb)
    for k in 1:n
        a = rescale_to_norm((chiA_copr[k, 1], chiA_copr[k, 2], chiA_copr[k, 3]), chiA_norm)
        b = rescale_to_norm((chiB_copr[k, 1], chiB_copr[k, 2], chiB_copr[k, 3]), chiB_norm)
        chiA_copr[k, :] .= a
        chiB_copr[k, :] .= b
    end
    orbphase = vec(spline_columns(t_coorb, dyn.tds, reshape(d.orbphase, :, 1)))
    quat = spline_columns(t_coorb, dyn.tds, d.quat)
    for k in 1:n
        nq = sqrt(quat[k, 1]^2 + quat[k, 2]^2 + quat[k, 3]^2 + quat[k, 4]^2)
        quat[k, :] ./= nq
    end

    # Spins in the coorbital frame: rotate about z by the orbital phase.
    T = promote_type(eltype(chiA_copr), eltype(chiB_copr))
    chiA_coorb = Matrix{T}(undef, n, 3)
    chiB_coorb = Matrix{T}(undef, n, 3)
    qfull = Matrix{T}(undef, n, 4)
    for k in 1:n
        sp = sin(orbphase[k])
        cp = cos(orbphase[k])
        chiA_coorb[k, 1] = chiA_copr[k, 1] * cp + chiA_copr[k, 2] * sp
        chiA_coorb[k, 2] = chiA_copr[k, 2] * cp - chiA_copr[k, 1] * sp
        chiA_coorb[k, 3] = chiA_copr[k, 3]
        chiB_coorb[k, 1] = chiB_copr[k, 1] * cp + chiB_copr[k, 2] * sp
        chiB_coorb[k, 2] = chiB_copr[k, 2] * cp - chiB_copr[k, 1] * sp
        chiB_coorb[k, 3] = chiB_copr[k, 3]
        # Frame of the orbit: coprecessing frame rotated by the orbital phase about z.
        qr = (cos(orbphase[k] / 2), zero(T), zero(T), sin(orbphase[k] / 2))
        qq = quat_multiply((quat[k, 1], quat[k, 2], quat[k, 3], quat[k, 4]), qr)
        qfull[k, 1], qfull[k, 2], qfull[k, 3], qfull[k, 4] = qq
    end

    h_coorb = coorbital_modes(coorb, q, chiA_coorb, chiB_coorb; ell_max)
    h = rotate_modes(h_coorb, qfull, ell_max)

    if dt === nothing && times === nothing
        keep = t_low === nothing ? (1:n) : findall(>=(t_low), t_coorb)
        return (t = t_coorb[keep], modes = h[keep, :])
    end
    t_out = if times !== nothing
        collect(times)
    else
        t0 = t_low === nothing ? t_coorb[1] : t_low
        # A plain loop: Enzyme cannot differentiate the float range `t0 .+ dt .* (0:N-1)`.
        nt = ceil(Int, (t_coorb[end] - t0) / dt)
        grid = Vector{typeof(t0 + dt)}(undef, nt)
        for k in 1:nt
            grid[k] = t0 + dt * (k - 1)
        end
        grid
    end
    t_out[end] <= t_coorb[end] + 0.01 || throw(ArgumentError("times extend beyond the end of the model"))
    t_out[1] >= t_coorb[1] || throw(ArgumentError("times start before the start of the model"))
    t_low === nothing || t_out[1] >= t_low || throw(ArgumentError("times start before the f_low start time"))
    return (t = t_out, modes = spline_columns(t_out, t_coorb, h))
end

# ₋₂Y_ℓm(θ, φ) = sqrt((2ℓ+1)/4π) d^ℓ_{m,2}(θ) e^{imφ}, matching `gwtools.harmonics.sYlm`.
@inline function sylm_m2(l::Integer, m::Integer, theta, phi)
    d = real(wigner_D_element(l, m, 2, complex(cos(theta / 2)), complex(sin(theta / 2))))
    return sqrt((2l + 1) / (4 * oftype(float(theta), π))) * d * cis(m * phi)
end

"""
    mode_sum(modes, ell_max, inclination, phi_ref) -> Vector

Complex strain `h = h₊ - i h×` from the mode matrix of [`inertial_modes`](@ref), evaluated on the
sky of the reference frame at polar angle `inclination` and azimuth `π/2 - phi_ref` (the LAL
convention).
"""
function mode_sum(modes::AbstractMatrix, ell_max::Integer, inclination, phi_ref)
    theta = inclination
    phi = oftype(float(phi_ref), π) / 2 - phi_ref
    T = promote_type(real(eltype(modes)), typeof(float(theta)), typeof(float(phi)))
    h = zeros(Complex{T}, size(modes, 1))
    for l in 2:ell_max, m in (-l):l
        y = sylm_m2(l, m, theta, phi)
        c = mode_column(l, m)
        for k in eachindex(h)
            h[k] += y * modes[k, c]
        end
    end
    return h
end

# Physical constants of gwtools (SI): solar mass in seconds, speed of light and the parsec.
const MSUN_IN_SECONDS = 4.925491025543576e-6
const SPEED_OF_LIGHT = 299792458.0
const PARSEC = 3.085677581491367e+16

"""
    PrecessingGenerator(surrogate; ell_max, f_ref, f_low)

Time-domain waveform generator around a [`PrecessingSurrogate`](@ref): the settings that are
fixed across draws. `f_ref` and `f_low` are the reference and starting frequencies of the
(2, 2) mode in Hz (`f_low = 0` returns the whole surrogate); the reference epoch is where the
coprecessing-frame orbital frequency equals `f_ref / 2`.
"""
struct PrecessingGenerator{S <: PrecessingSurrogate}
    surrogate::S
    ell_max::Int
    f_ref::Float64
    f_low::Float64
end

function PrecessingGenerator(
        surrogate::PrecessingSurrogate; ell_max::Integer = GWSurrogates.ell_max(surrogate), f_ref::Real, f_low::Real = 0.0,
    )
    2 <= ell_max <= GWSurrogates.ell_max(surrogate) || throw(ArgumentError("ell_max out of range"))
    f_ref >= f_low >= 0 || throw(ArgumentError("need f_ref >= f_low >= 0"))
    return PrecessingGenerator(surrogate, Int(ell_max), Float64(f_ref), Float64(f_low))
end

"""
    hphc_td(generator, params; dt = nothing, times = nothing) -> (t, hplus, hcross)

Plus and cross polarizations in SI units (seconds, strain) of the binary `params`, a named tuple
with `total_mass` (solar masses), `luminosity_distance` (Mpc), `mass_ratio` (`q >= 1`),
`spin_1`, `spin_2` (3-vectors in the LAL source frame at the reference epoch), `inclination`
(angle between the orbital angular momentum at the reference epoch and the line of sight) and
`phi_ref`; optional `init_orbphase` and `init_quat`. Time is zero at the peak of the waveform.
Specify the grid with the step `dt` (s) or the array `times` (s); the default is the model grid.
"""
function hphc_td(g::PrecessingGenerator, params::NamedTuple; dt = nothing, times = nothing)
    M = params.total_mass
    t_scale = MSUN_IN_SECONDS * M
    amp_scale = M * MSUN_IN_SECONDS * SPEED_OF_LIGHT / (1.0e6 * params.luminosity_distance * PARSEC)
    r = inertial_modes(
        g.surrogate, params.mass_ratio, params.spin_1, params.spin_2;
        ell_max = g.ell_max, f_ref = g.f_ref * t_scale, f_low = g.f_low * t_scale,
        dt = dt === nothing ? nothing : dt / t_scale, times = times === nothing ? nothing : times ./ t_scale,
        init_orbphase = get(params, :init_orbphase, 0.0), init_quat = get(params, :init_quat, (1.0, 0.0, 0.0, 0.0)),
    )
    h = mode_sum(r.modes, g.ell_max, params.inclination, params.phi_ref)
    return (t = r.t .* t_scale, hplus = real.(h) .* amp_scale, hcross = -imag.(h) .* amp_scale)
end
