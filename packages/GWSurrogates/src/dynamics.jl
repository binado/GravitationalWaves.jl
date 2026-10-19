"""
    DynamicsSurrogate(params, settings, t, fits, omega_ref_max)

Dynamics surrogate of the precessing models: it evolves the 11-variable state
`y = (q0, qx, qy, qz, φ, χA…, χB…)` (coprecessing-frame quaternion, orbital phase and the two
spin vectors in the coprecessing frame) on the fixed node grid `t` with RK4 start-up steps and a
variable-step Adams–Bashforth 4 scheme, exactly as `gwsurrogate` does. `fits` holds nine
parametric fits per node, in the order `Ω⁰`, `Ω¹` (the coorbital-frame orbital angular velocity
components), `ω`, `χ̇A` (3), `χ̇B` (3).

The first six steps of `t` are half steps (three RK4 steps of doubled width); the integration
writes `length(t) - 3` rows, sampled at `output_times(s)`.
"""
struct DynamicsSurrogate{P}
    params::P
    settings::FitSettings
    t::Vector{Float64}
    fits::FitSet
    omega_ref_max::Float64
    tds::Vector{Float64}
    diff_t::Vector{Float64}
    dt_ab4::Vector{Float64}
end

function DynamicsSurrogate(params, settings::FitSettings, t::Vector{Float64}, fits::FitSet, omega_ref_max::Real)
    length(fits) == 9 * length(t) || throw(DimensionMismatch("expected nine fits per node"))
    d = diff(t)
    for i in 1:3
        d[2i - 1] == d[2i] || throw(ArgumentError("AB4 needs three half-step pairs at the start of t"))
    end
    tds = vcat(t[1:2:5], t[7:end])
    dt_ab4 = vcat(2 .* d[1:2:5], d[7:end])
    return DynamicsSurrogate(params, settings, t, fits, Float64(omega_ref_max), tds, d, dt_ab4)
end

"Times at which [`dynamics`](@ref) samples the state: the full nodes of `t`, skipping the half-step nodes."
output_times(s::DynamicsSurrogate) = s.tds

const State = SVector{11}

# Fit inputs: q and the spins rotated into the coorbital frame by the orbital phase.
@inline function ds_fit_x(q, y)
    sp = sin(y[5])
    cp = cos(y[5])
    return (
        q, y[6] * cp + y[7] * sp, -1 * y[6] * sp + y[7] * cp, y[8],
        y[9] * cp + y[10] * sp, -1 * y[9] * sp + y[10] * cp, y[11],
    )
end

# Fit parameters and filled power table at the state `y`.
@inline function ds_fit_powers!(powers, s::DynamicsSurrogate, q, y)
    x = ds_fit_x(q, y)
    fp = fit_params(s.params, x[1], (x[2], x[3], x[4]), (x[5], x[6], x[7]))
    return fit_powers!(powers, s.settings, fp)
end

"""
    node_derivative(s, node, q, y, powers)

Time derivative of the state `y` from the fits of the 0-based grid node `node` (the ODE
right-hand side). `powers` is a buffer of length `n_fit_powers(s.settings)`.
"""
function node_derivative(s::DynamicsSurrogate, node::Integer, q, y, powers)
    ds_fit_powers!(powers, s, q, y)
    f = 9 * node
    r1 = eval_fit(s.fits, f + 1, powers)
    r2 = eval_fit(s.fits, f + 2, powers)
    r3 = eval_fit(s.fits, f + 3, powers)
    r4 = eval_fit(s.fits, f + 4, powers)
    r5 = eval_fit(s.fits, f + 5, powers)
    r6 = eval_fit(s.fits, f + 6, powers)
    r7 = eval_fit(s.fits, f + 7, powers)
    r8 = eval_fit(s.fits, f + 8, powers)
    r9 = eval_fit(s.fits, f + 9, powers)
    cp = cos(y[5])
    sp = sin(y[5])
    # Rotate the coorbital-frame quantities to the coprecessing frame.
    ox = r1 * cp - r2 * sp
    oy = r1 * sp + r2 * cp
    return State(
        (-0.5) * y[2] * ox - 0.5 * y[3] * oy,
        (-0.5) * y[4] * oy + 0.5 * y[1] * ox,
        0.5 * y[4] * ox + 0.5 * y[1] * oy,
        0.5 * y[2] * oy - 0.5 * y[3] * ox,
        r3,
        r4 * cp - r5 * sp, r4 * sp + r5 * cp, r6,
        r7 * cp - r8 * sp, r7 * sp + r8 * cp, r9,
    )
end

"Orbital frequency `ω` of the fits of 0-based node `node` at the state `y`."
function node_omega(s::DynamicsSurrogate, node::Integer, q, y, powers)
    ds_fit_powers!(powers, s, q, y)
    return eval_fit(s.fits, 9 * node + 3, powers)
end

# Natural cubic spline through four knots, evaluated at `x`; the data are vectors.
@inline function interp4(t1, t2, t3, t4, u1, u2, u3, u4, x)
    h1 = t2 - t1
    h2 = t3 - t2
    h3 = t4 - t3
    r2 = 6 * ((u3 - u2) / h2 - (u2 - u1) / h1)
    r3 = 6 * ((u4 - u3) / h3 - (u3 - u2) / h2)
    a22 = 2 * (h1 + h2)
    a33 = 2 * (h2 + h3)
    det = a22 * a33 - h2 * h2
    m2 = (r2 * a33 - h2 * r3) / det
    m3 = (a22 * r3 - h2 * r2) / det
    z = zero(m2)
    seg = 1 + (x >= t2) + (x >= t3)
    ts = (t1, t2, t3, t4)
    us = (u1, u2, u3, u4)
    ms = (z, m2, m3, z)
    h = ts[seg + 1] - ts[seg]
    a = (ts[seg + 1] - x) / h
    b = (x - ts[seg]) / h
    return a * us[seg] + b * us[seg + 1] + ((a^3 - a) * ms[seg] + (b^3 - b) * ms[seg + 1]) * (h^2 / 6)
end

"""
    time_derivative(s, t, q, y, powers)

State derivative at an arbitrary time `t` in the node range: the derivatives at the four
nearest nodes (all evaluated at `y`) interpolated with a natural cubic spline, as in
`gwsurrogate`.
"""
function time_derivative(s::DynamicsSurrogate, t, q, y, powers)
    tn = s.t
    n = length(tn)
    i0 = 1
    best = abs(tn[1] - t)
    for i in 2:n
        d = abs(tn[i] - t)
        if d < best
            best = d
            i0 = i
        end
    end
    imin = ifelse(t > tn[i0], i0 - 1, i0 - 2)
    imin = min(max(1, imin), n - 3)
    d1 = node_derivative(s, imin - 1, q, y, powers)
    d2 = node_derivative(s, imin, q, y, powers)
    d3 = node_derivative(s, imin + 1, q, y, powers)
    d4 = node_derivative(s, imin + 2, q, y, powers)
    return interp4(tn[imin], tn[imin + 1], tn[imin + 2], tn[imin + 3], d1, d2, d3, d4, t)
end

# Adams–Bashforth 4 increment with variable steps (`precessing_utils.c`, `ab4_dy`).
@inline function ab4_dy(k1, k2, k3, k4, dt1, dt2, dt3, dt4)
    dt12 = dt1 + dt2
    dt123 = dt12 + dt3
    dt23 = dt2 + dt3
    D1 = dt1 * dt12 * dt123
    D2 = dt1 * dt2 * dt23
    D3 = dt2 * dt12 * dt3
    B41 = dt3 * dt23 / D1
    B42 = -1 * dt3 * dt123 / D2
    B43 = dt23 * dt123 / D3
    B4 = B41 + B42 + B43
    C41 = (dt23 + dt3) / D1
    C42 = -1 * (dt123 + dt3) / D2
    C43 = (dt123 + dt23) / D3
    C4 = C41 + C42 + C43
    A = k4
    B = k4 * B4 - k1 * B41 - k2 * B42 - k3 * B43
    C = k4 * C4 - k1 * C41 - k2 * C42 - k3 * C43
    D = (k4 - k1) / D1 - (k4 - k2) / D2 + (k4 - k3) / D3
    return dt4 * (A + dt4 * (0.5 * B + dt4 * (C / 3.0 + dt4 * 0.25 * D)))
end

"""
    normalize_state(y, normA, normB)

Renormalize the quaternion to unit length and the spin vectors to the magnitudes `normA` and
`normB` of the initial spins. Zero spins are left unchanged ([`rescale_to_norm`](@ref)).
"""
@inline function normalize_state(y, normA, normB)
    quat_norm = sqrt(y[1] * y[1] + y[2] * y[2] + y[3] * y[3] + y[4] * y[4])
    cA = rescale_to_norm((y[6], y[7], y[8]), normA)
    cB = rescale_to_norm((y[9], y[10], y[11]), normB)
    return State(
        y[1] / quat_norm, y[2] / quat_norm, y[3] / quat_norm, y[4] / quat_norm, y[5],
        cA[1], cA[2], cA[3], cB[1], cB[2], cB[3],
    )
end

# Row r (0-based) of the output is node t[r + 1] for r < 3 doubled steps and t[r + 4] after.
@inline node_time_index(r) = ifelse(r < 3, 2r, r + 3)

# One RK4 step from row `r`, forward or backward (`direction` = ±1), with spline-interpolated
# derivatives; returns the new state and the first stage derivative.
function rk4_step(s::DynamicsSurrogate, q, y, r, direction, normA, normB, powers)
    it = node_time_index(r)
    t1 = s.t[it + 1]
    if direction > 0
        t2 = s.t[ifelse(r < 3, it + 2, it + 1) + 1]
    else
        t2 = s.t[ifelse(r <= 3, it - 2, it - 1) + 1]
    end
    half_dt = 0.5 * (t2 - t1)
    y0 = y[r + 1]
    k1 = time_derivative(s, t1, q, y0, powers)
    k2 = time_derivative(s, t1 + half_dt, q, y0 + half_dt * k1, powers)
    k3 = time_derivative(s, t1 + half_dt, q, y0 + half_dt * k2, powers)
    k4 = time_derivative(s, t2, q, y0 + 2 * half_dt * k3, powers)
    ynext = y0 + (half_dt / 3.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    y[r + direction + 1] = normalize_state(ynext, normA, normB)
    return k1
end

function integrate_forward!(s::DynamicsSurrogate, q, y, normA, normB, i0, k_ab4, dt_ab4, powers)
    i0 >= 3 || throw(ArgumentError("i0 must be at least 3"))
    k1, k2, k3 = k_ab4
    dt1, dt2, dt3 = dt_ab4
    d = s.diff_t
    for i in 0:(length(d) - (i0 + 3) - 1)
        dt4 = d[i0 + 3 + i + 1]
        r = i0 + i
        k4 = node_derivative(s, r + 3, q, y[r + 1], powers)
        ynext = y[r + 1] + ab4_dy(k1, k2, k3, k4, dt1, dt2, dt3, dt4)
        y[r + 2] = normalize_state(ynext, normA, normB)
        k1, k2, k3 = k2, k3, k4
        dt1, dt2, dt3 = dt2, dt3, dt4
    end
    return y
end

function integrate_backward!(s::DynamicsSurrogate, q, y, normA, normB, i0, k_ab4, dt_ab4, powers)
    i0 <= length(s.t) - 7 || throw(ArgumentError("i0 must be <= length(t) - 7"))
    k1, k2, k3 = k_ab4
    dt1, dt2, dt3 = dt_ab4
    for r in (i0 - 1):-1:0
        node = ifelse(r < 2, 2 + 2r, r + 4)
        dt4 = s.dt_ab4[r + 1]
        k4 = node_derivative(s, node, q, y[r + 2], powers)
        ynext = y[r + 2] - ab4_dy(k1, k2, k3, k4, dt1, dt2, dt3, dt4)
        y[r + 1] = normalize_state(ynext, normA, normB)
        k1, k2, k3 = k2, k3, k4
        dt1, dt2, dt3 = dt2, dt3, dt4
    end
    return y
end

# Rotate the (x, y) components of a spin about z by `phase` (`rotate_spin`).
@inline function rotate_spin(chi, phase)
    sp = sin(phase)
    cp = cos(phase)
    return (chi[1] * cp + chi[2] * sp, chi[2] * cp - chi[1] * sp, chi[3])
end

"""
    t_ref_from_omega(s, omega_ref, q, chiA0, chiB0; init_orbphase = 0, init_quat = (1, 0, 0, 0))

Reference time at which the orbital frequency reaches `omega_ref`, by linear interpolation of
`ω` between nodes (as `gwsurrogate`). The result is only C⁰ in `omega_ref` and the spins: its
derivative jumps at every node crossing, so gradient-based samplers should fix `t_ref` instead.
"""
function t_ref_from_omega(
        s::DynamicsSurrogate, omega_ref, q, chiA0, chiB0;
        init_orbphase = 0.0, init_quat = (1.0, 0.0, 0.0, 0.0),
    )
    omega_ref <= s.omega_ref_max || throw(ArgumentError("omega_ref = $omega_ref exceeds the model limit $(s.omega_ref_max)"))
    a = rotate_spin(chiA0, -1 * init_orbphase)
    b = rotate_spin(chiB0, -1 * init_orbphase)
    T = promote_type(typeof(q), eltype(a), eltype(b), typeof(init_orbphase), Float64)
    y0 = State{T}(init_quat[1], init_quat[2], init_quat[3], init_quat[4], init_orbphase, a..., b...)
    powers = Vector{T}(undef, n_fit_powers(s.settings))
    omega0 = node_omega(s, 0, q, y0, powers)
    omega_ref >= omega0 || throw(ArgumentError("omega_ref = $omega_ref is below omega_0 = $omega0"))
    # The three half-step nodes (1, 3, 5) are skipped, as in LAL.
    nodes = [i for i in 0:(length(s.t) - 1) if i != 1 && i != 3 && i != 5]
    imax = 2
    omega_min = omega0
    omega_max = node_omega(s, nodes[imax], q, y0, powers)
    while omega_max <= omega_ref
        imax += 1
        omega_min = omega_max
        omega_max = node_omega(s, nodes[imax], q, y0, powers)
    end
    t_min = s.t[nodes[imax - 1] + 1]
    t_max = s.t[nodes[imax] + 1]
    return (t_min * (omega_max - omega_ref) + t_max * (omega_ref - omega_min)) / (omega_max - omega_min)
end

"""
    dynamics(s, q, chiA0, chiB0; init_quat = (1, 0, 0, 0), init_orbphase = 0, t_ref = nothing)

Evolve the dynamics surrogate from the initial conditions at the reference time `t_ref`
(default: the first node). `chiA0`, `chiB0` are the spin 3-vectors in the LAL source frame at
the reference epoch. Returns a named tuple `(t, quat, orbphase, chiA, chiB)` sampled at
[`output_times`](@ref): `quat` is `(n, 4)`, the spins are `(n, 3)` in the coprecessing frame.

If `t_ref` is not a node, one forward-Euler step to the nearest node is taken first. Magnitudes
of the spins are preserved by the per-step renormalization.
"""
function dynamics(
        s::DynamicsSurrogate, q, chiA0, chiB0;
        init_quat = (1.0, 0.0, 0.0, 0.0), init_orbphase = 0.0, t_ref = nothing,
    )
    a = rotate_spin(chiA0, -1 * init_orbphase)
    b = rotate_spin(chiB0, -1 * init_orbphase)
    normA = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    normB = sqrt(b[1]^2 + b[2]^2 + b[3]^2)
    max(normA, normB) <= 1.001 || throw(ArgumentError("spin magnitude $(max(normA, normB)) exceeds 1"))
    T = promote_type(typeof(q), eltype(a), eltype(b), typeof(init_orbphase), eltype(init_quat), Float64)
    if t_ref !== nothing
        s.t[1] <= t_ref <= s.t[end] || throw(ArgumentError("t_ref outside the time domain"))
    end
    powers = Vector{T}(undef, n_fit_powers(s.settings))
    y0 = State{T}(init_quat[1], init_quat[2], init_quat[3], init_quat[4], init_orbphase, a..., b...)
    y = Vector{State{T}}(undef, length(s.t) - 3)

    i0 = 0
    if t_ref === nothing
        y[1] = y0
    else
        i0 = 0
        best = abs(s.tds[1] - t_ref)
        for r in 1:length(s.tds)
            d = abs(s.tds[r] - t_ref)
            if d < best
                best = d
                i0 = r - 1
            end
        end
        dydt0 = time_derivative(s, t_ref, q, y0, powers)
        y[i0 + 1] = normalize_state(y0 + (s.tds[i0 + 1] - t_ref) * dydt0, normA, normB)
    end

    if i0 == 0
        d = s.diff_t
        k_ab4 = ntuple(_ -> zero(State{T}), 3)
        for i in 0:2
            dt = d[2i + 1]
            yi = y[i + 1]
            k1 = node_derivative(s, 2i, q, yi, powers)
            k2 = node_derivative(s, 2i + 1, q, yi + dt * k1, powers)
            k3 = node_derivative(s, 2i + 1, q, yi + dt * k2, powers)
            k4 = node_derivative(s, 2i + 2, q, yi + 2 * dt * k3, powers)
            ynext = yi + (dt / 3.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            y[i + 2] = normalize_state(ynext, normA, normB)
            k_ab4 = Base.setindex(k_ab4, k1, i + 1)
        end
        dt_ab4 = (s.dt_ab4[1], s.dt_ab4[2], s.dt_ab4[3])
        integrate_forward!(s, q, y, normA, normB, 3, k_ab4, dt_ab4, powers)
    elseif i0 > 2
        k_ab4 = ntuple(_ -> zero(State{T}), 3)
        for i in 0:2
            k_ab4 = Base.setindex(k_ab4, rk4_step(s, q, y, i0 - i, -1, normA, normB, powers), i + 1)
        end
        dt_ab4 = (s.dt_ab4[i0], s.dt_ab4[i0 - 1], s.dt_ab4[i0 - 2])
        integrate_backward!(s, q, y, normA, normB, i0 - 3, k_ab4, dt_ab4, powers)
        tmp_k = node_derivative(s, i0, q, y[i0 - 3 + 1], powers)
        integrate_forward!(
            s, q, y, normA, normB, i0, (tmp_k, k_ab4[3], k_ab4[2]), (dt_ab4[3], dt_ab4[2], dt_ab4[1]), powers,
        )
    else
        k_ab4 = ntuple(_ -> zero(State{T}), 3)
        for i in 0:2
            k_ab4 = Base.setindex(k_ab4, rk4_step(s, q, y, i0 + i, 1, normA, normB, powers), i + 1)
        end
        dt_ab4 = (s.dt_ab4[i0 + 1], s.dt_ab4[i0 + 2], s.dt_ab4[i0 + 3])
        integrate_forward!(s, q, y, normA, normB, i0 + 3, k_ab4, dt_ab4, powers)
        tmp_k = node_derivative(s, i0 + 3, q, y[i0 + 3 + 1], powers)
        integrate_backward!(
            s, q, y, normA, normB, i0, (tmp_k, k_ab4[3], k_ab4[2]), (dt_ab4[3], dt_ab4[2], dt_ab4[1]), powers,
        )
    end

    n = length(y)
    return (
        t = s.tds,
        quat = T[y[r][c] for r in 1:n, c in 1:4],
        orbphase = T[y[r][5] for r in 1:n],
        chiA = T[y[r][5 + c] for r in 1:n, c in 1:3],
        chiB = T[y[r][8 + c] for r in 1:n, c in 1:3],
    )
end
