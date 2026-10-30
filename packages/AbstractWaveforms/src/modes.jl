"""
    HigherModesWaveformGenerator <: AbstractWaveformGenerator

Generator exposing fixed-frame spin-weight -2 modes through [`waveform_modes_td`](@ref)
or [`waveform_modes_fd`](@ref), inheriting polarization assembly in the supported domain.
No transform between domains is implied. Backends document their frame and reference phase.
"""
abstract type HigherModesWaveformGenerator <: AbstractWaveformGenerator end

supports_mode_array(::HigherModesWaveformGenerator) = true

"""
    waveform_modes_fd(f, params, generator) -> (; mode_indices, positive, negative)

Distance-normalized, inclination-independent complex modes on the nonnegative uniform grid
`f`. `mode_indices` is an ordered tuple of unique signed `(l,m)` labels with `2 ≤ l ≤ 5`.
Matrices have shape `(length(f), length(mode_indices))`; `positive[i,j]` is h̃ⱼ(f[i])
and `negative[i,j]` is h̃ⱼ(-f[i]), in the same row order. DC entries agree. Both branches
vanish outside configured frequency bounds and intrinsic model support. Reference phase
belongs in the modes by default; a backend may override [`mode_angles`](@ref).
"""
function waveform_modes_fd end

"""
    waveform_modes_td(generator, params; kwargs...) -> (; t, mode_indices, modes)

Fixed inertial/source-frame spherical modes in physical strain, excluding observer inclination.
`t` is a vector of seconds; complex `modes` has shape `(length(t), length(mode_indices))`.
Labels follow [`waveform_modes_fd`](@ref). The backend documents its reference-phase convention.
"""
function waveform_modes_td end

"""
    hphc_td(generator, params; kwargs...) -> TimeSeries

Time-domain polarizations in seconds and physical strain. Higher-mode generators evaluate
[`waveform_modes_td`](@ref) once and assemble H = Σ ₋₂Yₗₘ hₗₘ, hplus = real(H),
hcross = -imag(H), with angles from [`mode_angles`](@ref).
Times, mode labels and matrix dimensions are assumed to satisfy the [`waveform_modes_td`](@ref)
contract and are not revalidated during assembly.
"""
function hphc_td end

"""
    mode_angles(params, generator) -> (theta, phi)

Observer angles in the mode frame. Defaults to `(params.theta_jn, 0)`; the reference
phase is included in the modes. Backends may override this convention.
"""
mode_angles(params, ::HigherModesWaveformGenerator) = (params.theta_jn, zero(params.theta_jn))

"""
    mode_sum(modes, mode_indices, theta, phi)

Assemble complex strain H = Σ ₋₂Yₗₘ(theta,phi) hₗₘ from `(sample, mode)` complex
modes. Compute angles once and accumulate without weighted-array intermediates.
Mode labels and matrix dimensions are assumed to satisfy the waveform mode contract.
The azimuth `phi` is the actual mode-frame viewing angle.
"""
function mode_sum(modes::AbstractMatrix, labels, theta, phi)
    theta, phi = promote(float(theta), float(phi))
    coefficients = mode_coefficients(labels, theta, phi)
    return sum_modes(modes, coefficients)
end

function sum_modes(modes, coefficients)
    T = promote_type(eltype(modes), eltype(coefficients))
    h = zeros(T, size(modes, 1))
    for j in eachindex(coefficients)
        y = coefficients[j]
        for k in eachindex(h)
            h[k] += y * modes[k, j]
        end
    end
    return h
end

"""Evaluate plus/cross columns at backend-supported, peak-relative time coordinates."""
function hphc_td(t::AbstractVector{<:Real}, params::NamedTuple, g::HigherModesWaveformGenerator; kwargs...)
    r = waveform_modes_td(g, params; times = t, kwargs...)
    return assemble_td(r, params, g)
end

function assemble_td(r, params, g)
    theta, phi = mode_angles(params, g)
    h = sum_modes(r.modes, mode_coefficients(r.mode_indices, theta, phi))
    return hcat(real.(h), -imag.(h))
end

"""Sample polarizations at a time prototype's peak-relative coordinates, without detector placement."""
function hphc_td(prototype::TimeSeries, params::NamedTuple, g::AbstractWaveformGenerator; kwargs...)
    data = hphc_td(times(prototype), params, g; kwargs...)
    return TimeSeries(data; epoch = epoch(prototype), sampling_interval = sampling_interval(prototype), channels = (:plus, :cross))
end

function hphc_td(g::HigherModesWaveformGenerator, params::NamedTuple; dt = nothing, times = nothing, kwargs...)
    (dt === nothing) != (times === nothing) || throw(ArgumentError("assembled TD waveforms require exactly one of dt and a uniform times range"))
    if times !== nothing
        times isa AbstractRange{<:Real} && step(times) > 0 && !isempty(times) || throw(ArgumentError("times must be a nonempty increasing uniform range"))
        data = hphc_td(times, params, g; kwargs...)
        return TimeSeries(data; epoch = first(times), sampling_interval = step(times), channels = (:plus, :cross))
    end
    isfinite(dt) && dt > 0 || throw(ArgumentError("dt must be finite and positive"))
    r = waveform_modes_td(g, params; dt, kwargs...)
    return TimeSeries(assemble_td(r, params, g); epoch = first(r.t), sampling_interval = dt, channels = (:plus, :cross))
end

function hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, g::HigherModesWaveformGenerator)
    r = waveform_modes_fd(f, params, g)
    theta, phi = mode_angles(params, g)
    coefficients = mode_coefficients(r.mode_indices, theta, phi)
    hp = sum_modes(r.positive, coefficients)
    hn = sum_modes(r.negative, coefficients)
    return hcat((hp .+ conj.(hn)) ./ 2, (im / 2) .* (hp .- conj.(hn)))
end
