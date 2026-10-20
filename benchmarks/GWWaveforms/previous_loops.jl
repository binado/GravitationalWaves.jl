# Direct-frequency loops before scalar-kernel extraction, retained as the benchmark baseline.
module PreviousLoops

using GWWaveforms: IMRPhenomXAS, phenom_frequency_slice, phenom_reference_frequency,
    F_CUT, MTSUN_SI, phenomd_model, phenom_amplitude_scale, useful_powers,
    _phenomd_amplitude, _phenomd_phase, canonical_amplitude_phase, phenomxas_model,
    phenomxas_amplitude, phenomxas_phase

function baseline_phenomd(
        ::Type{T}, f::AbstractRange, lo::Integer, hi::Integer, p, g, f_start, f_max
    ) where {T}
    amplitude = zeros(T, length(f))
    phase = zeros(T, length(f))
    lo > hi && return (; amplitude, phase)
    m1 = p.mass_1
    m2 = p.mass_2
    fRef = phenom_reference_frequency(g, f_start)
    M = m1 + m2
    fCut = F_CUT / (M * MTSUN_SI)
    f_max_eff = ifelse(fCut < f_max, fCut, f_max)

    model = phenomd_model(m1, m2, p.chi_1, p.chi_2, fRef, p.phase)
    amp0 = phenom_amplitude_scale(m1, m2, p.luminosity_distance)

    deltaF = Float64(step(f))
    stop = floor(f_max_eff / deltaF) * deltaF
    for i in lo:hi
        frequency = Float64(f[i])
        if frequency < stop
            Mf = model.M_sec * frequency
            pm = useful_powers(Mf)
            amp = amp0 * _phenomd_amplitude(model.amp, Mf, pm)
            phi = _phenomd_phase(model, Mf, pm) - (model.t0 * (Mf - model.MfRef) + model.phi_precalc)
            amplitude[i], phase[i] = canonical_amplitude_phase(amp, -phi)
        end
    end
    return (; amplitude, phase)
end

function baseline_phenomxas(g::IMRPhenomXAS, ::Type{T}, f, p, f_start, f_min, f_max) where {T}
    lo, hi = phenom_frequency_slice(f, f_start, f_min, f_max)
    amplitude = zeros(T, length(f))
    phase = zeros(T, length(f))
    lo > hi && return (; amplitude, phase)
    fref = phenom_reference_frequency(g, f_start)
    m = phenomxas_model(p.mass_1, p.mass_2, p.chi_1, p.chi_2, fref, p.phase)
    amp0 = phenom_amplitude_scale(p.mass_1, p.mass_2, p.luminosity_distance)
    norm = amp0 * sqrt(2m.r.eta / 3) * Float64(π)^(-1 / 6)
    cutoff = 0.3 / m.M_sec
    for i in lo:hi
        frequency = Float64(f[i])
        if frequency <= cutoff
            Mf = m.M_sec * frequency
            amp = norm * phenomxas_amplitude(m.amplitude, Mf)
            phi = phenomxas_phase(m.phase, Mf) + m.linb * Mf + m.phiref
            amplitude[i], phase[i] = canonical_amplitude_phase(amp, phi)
        end
    end
    return (; amplitude, phase)
end

end
