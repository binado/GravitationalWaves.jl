"""
    SEOBNRv5HM(; minimum_frequency=20, reference_frequency=20,
                maximum_frequency=nothing, reltol=1e-11, abstol=1e-12,
                post_merger_duration=1)

Experimental aligned-spin SEOBNRv5HM prototype restricted to (2,±2) output.
Configuration and sampling must remain fixed during differentiation.
`is_differentiable` remains false until complete first-order ForwardDiff CPU
acceptance. Higher modes contribute to the flux, but are not output modes.
"""
struct SEOBNRv5HM <: AbstractWaveformGenerator
    minimum_frequency::Float64
    reference_frequency::Float64
    maximum_frequency::Union{Nothing, Float64}
    reltol::Float64
    abstol::Float64
    post_merger_duration::Float64
    function SEOBNRv5HM(;
            minimum_frequency::Real = 20,
            reference_frequency::Real = 20,
            maximum_frequency::Union{Nothing, Real} = nothing,
            reltol::Real = 1.0e-11, abstol::Real = 1.0e-12,
            post_merger_duration::Real = 1,
        )
        isfinite(minimum_frequency) && minimum_frequency > 0 ||
            throw(ArgumentError("minimum_frequency must be finite and positive"))
        reference_frequency == minimum_frequency ||
            throw(ArgumentError("reference_frequency must equal minimum_frequency"))
        isnothing(maximum_frequency) ||
            (isfinite(maximum_frequency) && maximum_frequency >= minimum_frequency) ||
            throw(ArgumentError("maximum_frequency must be finite and at least minimum_frequency"))
        all(x -> isfinite(x) && x > 0, (reltol, abstol, post_merger_duration)) ||
            throw(ArgumentError("tolerances and post_merger_duration must be finite and positive"))
        return new(
            Float64(minimum_frequency), Float64(reference_frequency),
            isnothing(maximum_frequency) ? nothing : Float64(maximum_frequency),
            Float64(reltol), Float64(abstol), Float64(post_merger_duration),
        )
    end
end

approximant(::SEOBNRv5HM) = "SEOBNRv5HM"
minimum_frequency(g::SEOBNRv5HM) = g.minimum_frequency
reference_frequency(g::SEOBNRv5HM) = g.reference_frequency
maximum_frequency(g::SEOBNRv5HM) = g.maximum_frequency
spin_support(::SEOBNRv5HM) = :aligned
supports_mode_array(::SEOBNRv5HM) = false
is_differentiable(::SEOBNRv5HM) = false
