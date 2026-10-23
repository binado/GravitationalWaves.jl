# Peak-normalized error avoids division by tiny individual samples. Zero outputs
# are checked exactly instead of normalizing by zero.
function check_pointwise(actual, expected; rtol, exact_support = false)
    @test size(actual) == size(expected)
    @test all(isfinite, actual)
    @test all(isfinite, expected)
    peak = maximum(abs, expected)
    @test iszero(peak) ? all(iszero, actual) : maximum(abs, actual .- expected) <= rtol * peak
    exact_support && @test iszero.(actual) == iszero.(expected)
    return nothing
end

# FD generators are evaluated on `grid` (frequencies); TD generators on a time grid, then
# `prepare(::TimeSeries)` (window, padding, ...) and the physical rfft bring them to FD.
function evaluate_polarizations(g, params, grid, domain::Symbol, prepare)
    if domain === :fd
        return hphc_fd(grid, params, g), grid
    elseif domain === :td
        prototype = TimeSeries(
            zeros(length(grid), 2); epoch = first(grid), sampling_interval = step(grid), channels = (:plus, :cross),
        )
        spectrum = to_frequencydomain(prepare(hphc_td(prototype, params, g)))
        return parent(spectrum), frequencies(spectrum)
    end
    throw(ArgumentError("domain must be :fd or :td"))
end

# Intersection of the two generators' supported bands, clipped to the grid.
function default_band(f, gs...)
    low = max(first(f), maximum(minimum_frequency, gs))
    highs = filter(!isnothing, map(maximum_frequency, gs))
    high = min(last(f), highs...)
    return low, high
end

resolve_metric(metric::FrequencyDomainMetric, f, gs) = metric
resolve_metric(metric, f, gs) = metric(f)
function resolve_metric(::Nothing, f, gs)
    low, high = default_band(f, gs...)
    return flat_metric(f; minimum_frequency = low, maximum_frequency = high)
end

"""
    test_reference(target, reference, cases = applicable_cases(target, reference); grid, max_mismatch,
        parameters = fd_parameters, target_domain = :fd, reference_domain = target_domain,
        prepare = identity, metric = nothing, maximize = :time_phase, max_shift = 0.05,
        rtol = nothing, exact_support = false)

Compare a `target` generator with a trusted `reference` generator (for example `LALWaveforms` or the
PyCall `gwsurrogate` adapter) on every case, requiring `mismatch <= max_mismatch`.

`cases` are `(id, case)` pairs and `parameters` maps a case to the generators' parameter `NamedTuple`.
`grid` is the frequency grid for FD generators, or the time grid (a uniform range) for TD generators;
a `:td` side goes through `prepare` and the physical rfft, and both sides must end on the same
frequency grid. `metric` is `nothing` (flat PSD over the generators' common band), a
`FrequencyDomainMetric`, or a function of the frequency grid returning one.

When `rtol` is set, each polarization must additionally agree pointwise to `rtol` relative to the
reference peak; `exact_support` also requires identical zero patterns. This keeps strong tests for
ports meant to reproduce the reference bit-for-bit-ish, since mismatch is quadratic in the error.
"""
function test_reference(
        target, reference, cases = applicable_cases(target, reference);
        grid, max_mismatch, parameters = fd_parameters, target_domain::Symbol = :fd,
        reference_domain::Symbol = target_domain, prepare = identity, metric = nothing,
        maximize::Symbol = :time_phase, max_shift = 0.05, rtol = nothing, exact_support::Bool = false,
    )
    isempty(cases) && throw(ArgumentError("no cases are supported by both generators"))
    return @testset "reference: $(approximant(target)) vs $(approximant(reference))" begin
        for (id, case) in cases
            @testset "$id" begin
                params = parameters(case)
                h_target, f_target = evaluate_polarizations(target, params, grid, target_domain, prepare)
                h_reference, f_reference = evaluate_polarizations(reference, params, grid, reference_domain, prepare)
                @test size(h_target) == size(h_reference)
                @test all(isfinite, h_target) && all(isfinite, h_reference)
                @test collect(f_target) ≈ collect(f_reference)
                m = resolve_metric(metric, f_target, (target, reference))
                @test mismatch(h_target, h_reference, m; maximize, max_shift) <= max_mismatch
                if rtol !== nothing
                    for channel in axes(h_target, 2)
                        check_pointwise(h_target[:, channel], h_reference[:, channel]; rtol, exact_support)
                    end
                end
            end
        end
    end
end
