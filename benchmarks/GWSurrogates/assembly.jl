# Run with --project=packages/GWSurrogates. Model data is needed only for generation timings.
using GWSurrogates
using AbstractGWWaveforms: mode_sum, hphc_td, waveform_modes_td

function previous_sum(modes, ell_max, theta, phi_ref)
    phi = π / 2 - phi_ref
    h = zeros(ComplexF64, size(modes, 1))
    for l in 2:ell_max, m in (-l):l
        d = real(GWSurrogates.wigner_D_element(l, m, 2, complex(cos(theta / 2)), complex(sin(theta / 2))))
        y = sqrt((2l + 1) / (4π)) * d * cis(m * phi)
        for k in eachindex(h)
            h[k] += y * modes[k, GWSurrogates.mode_column(l, m)]
        end
    end
    return h
end

function measure(f; repeats = 20)
    f()
    times = [@elapsed f() for _ in 1:repeats]
    return (; seconds = sort(times)[cld(repeats, 2)], bytes = @allocated(f()))
end

function previous_strain(g, params; dt)
    t_scale = GWSurrogates.MSUN_IN_SECONDS * params.total_mass
    amp_scale = params.total_mass * GWSurrogates.MSUN_IN_SECONDS * GWSurrogates.SPEED_OF_LIGHT /
        (1.0e6 * params.luminosity_distance * GWSurrogates.PARSEC)
    r = GWSurrogates.inertial_modes(
        g.surrogate, params.mass_ratio, params.spin_1, params.spin_2;
        ell_max = g.ell_max, f_ref = g.f_ref * t_scale, f_low = g.f_low * t_scale, dt = dt / t_scale
    )
    h = previous_sum(r.modes, g.ell_max, params.inclination, params.phi_ref)
    return (; t = r.t .* t_scale, hplus = real.(h) .* amp_scale, hcross = -imag.(h) .* amp_scale)
end

println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME)
for ell_max in (4, 5), n in (1, 2000, 20000)
    labels = Tuple((l, m) for l in 2:ell_max for m in (-l):l)
    modes = [complex(sin(k + j), cos(k - j)) for k in 1:n, j in eachindex(labels)]
    old = () -> previous_sum(modes, ell_max, 0.9, 0.4)
    shared = () -> mode_sum(modes, labels, 0.9, π / 2 - 0.4)
    @assert old() ≈ shared()
    println((; ell_max, n, previous = measure(old), shared = measure(shared)))
end
for model in (:NRSur7dq4, :NRSur7dq4v2)
    surrogate = GWSurrogates.load_surrogate(model; download = false)
    g = GWSurrogates.PrecessingGenerator(surrogate; f_ref = 25.0)
    params = (;
        total_mass = 60.0, luminosity_distance = 400.0, mass_ratio = 2.2,
        spin_1 = (0.3, 0.2, 0.4), spin_2 = (-0.1, 0.5, -0.3), inclination = 0.7, phi_ref = 0.3,
    )
    old = () -> previous_strain(g, params; dt = 1 / 4096)
    shared = () -> hphc_td(g, params; dt = 1 / 4096)
    @assert old().hplus ≈ shared().hplus
    @assert old().hcross ≈ shared().hcross
    println(
        (;
            model, previous = measure(old), modes = measure(() -> waveform_modes_td(g, params; dt = 1 / 4096)),
            strain = measure(() -> hphc_td(g, params; dt = 1 / 4096)),
        )
    )
end
