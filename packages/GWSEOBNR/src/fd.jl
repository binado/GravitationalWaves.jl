"""
    hphc_fd(f, params, generator::SEOBNRv5HM)

Experimental native (2,±2) inspiral–merger–ringdown on a fixed physical FFT grid.
Inputs are checked. Insufficient onset/ringdown coverage raises `SolveError`.
First-order ForwardDiff support is experimental until complete acceptance passes.
"""
function hphc_fd(f::AbstractRange, params::NamedTuple, g::SEOBNRv5HM)
    validate_parameters(f, params, g)
    return evaluated_spectrum(f, aligned_parameters(params), g)
end

# Internal signed-spin continuation used by acceptance finite differences at the
# envelope boundaries. Public evaluations always pass the host validation above.
function evaluated_spectrum(f, p, g)
    x = collect(promote(p.mass_1, p.mass_2, p.chi_1, p.chi_2, p.luminosity_distance, p.inclination, p.phase))
    y = spectrum_vector(f, x, g)
    n = length(f)
    return (; plus = y[1:n], cross = y[(n + 1):(2n)])
end

# Parameter order of `spectrum_vector`. Gives reverse-mode backends one differentiable
# boundary: 7 real inputs, `[plus; cross]` out. Backend extensions replace it with a
# Jacobian-based rule instead of differentiating the ODE solve.
function spectrum_vector(f, x::AbstractVector, g)
    p = (;
        mass_1 = x[1], mass_2 = x[2], chi_1 = x[3], chi_2 = x[4],
        luminosity_distance = x[5], inclination = x[6], phase = x[7],
    )
    T = 1 / step(f)
    N = 2(length(f) - 1)
    dt = T / N
    g.post_merger_duration >= 0.25 || throw(SolveError("post_merger_duration must include the final 0.25-second taper"))
    e = evolve(p.mass_1, p.mass_2, p.chi_1, p.chi_2, g)
    nqc = nqc_coefficients(e)
    mr = merger_coefficients(e, nqc)
    tpeak = peak_time(e, nqc, mr)
    Msec = (p.mass_1 + p.mass_2) * MTSUN_SI
    tstart = -T + g.post_merger_duration
    tstart <= -primal(tpeak * Msec) || throw(SolveError("FFT segment does not contain the full start transition"))
    ringdown_start = (e.t_attach - tpeak) * Msec
    tau = -Msec / imag(mr.omega)
    g.post_merger_duration - 0.25 >= primal(ringdown_start + 30tau) ||
        throw(SolveError("FFT segment needs 30 damping times before the final taper"))
    real(mr.omega) / Msec <= 2π * last(f) || throw(SolveError("grid does not resolve the ringdown frequency"))
    scale = -(p.mass_1 + p.mass_2) * MRSUN_SI / (p.luminosity_distance * MPC_SI)
    cosi = cos(p.inclination)
    Yplus = sqrt(5 / (64π)) * (1 + cosi)^2
    Yminus = sqrt(5 / (64π)) * (1 - cosi)^2
    rotation = cis(2 * (π / 2 - p.phase))
    function sample(j)
        trel = tstart + j * dt
        t = trel / Msec + tpeak
        if primal(t) < e.sol.t[1]
            return zero(complex(scale + p.inclination + p.phase + p.chi_1 + p.chi_2))
        end
        h22 = complete_mode(e, nqc, mr, t) * rotation
        onset = if primal(t) <= primal(e.t_attach)
            h0 = inspiral_mode(e, t)
            freq = -time_jet(s -> angle(inspiral_mode(e, s) / h0), t)[2] / (2π * Msec)
            quintic_step((freq - 0.8g.minimum_frequency) / (0.2g.minimum_frequency))
        else
            one(scale)
        end
        ending = quintic_step((g.post_merger_duration - trel) / 0.25)
        return scale * onset * ending * (Yplus * h22 + Yminus * conj(h22))
    end
    strain = [sample(j) for j in 0:(N - 1)]
    plus = physical_spectrum(real.(strain), dt, tstart)
    cross = physical_spectrum(-imag.(strain), dt, tstart)
    fmin, fmax = frequency_bounds(g, f)
    for i in eachindex(f)
        if !(fmin <= f[i] <= fmax)
            plus[i] = zero(eltype(plus))
            cross[i] = zero(eltype(cross))
        end
    end
    return vcat(plus, cross)
end
