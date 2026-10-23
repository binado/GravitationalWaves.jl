"""
    FixedObjective(f)

Wrap a callable whose captures are fixed configuration (grids, generators, metrics), never aliases of
the differentiated input. The Mooncake extension marks it non-differentiable.
"""
struct FixedObjective{F}
    f::F
end
(f::FixedObjective)(x) = f.f(x)

finite_difference_backend() = AutoFiniteDifferences(; fdm = central_fdm(5, 1))

relative(g, scale) = g ./ scale

"""
    test_gradients(objective, x; backends, rtol, finite_difference = true, fd_rtol = 1e-4,
        reuse_index = nothing, atol = 1e-8)

Self-calibrating gradient checks of a scalar `objective` at `x`, always through DifferentiationInterface:
each AD backend must reproduce the value, agree with the others, and (if `finite_difference`) with
`FiniteDifferences` at `fd_rtol`. Gradients are compared componentwise after dividing by the objective
value. If `reuse_index` is set, a gradient prepared at `x` is reused after scaling that coordinate by
1.001 and must match an unprepared evaluation. Returns the first backend's gradient.
"""
function test_gradients(
        objective, x::AbstractVector; backends, rtol, finite_difference::Bool = true, fd_rtol = 1.0e-4,
        reuse_index = nothing, atol = 1.0e-8,
    )
    loss = FixedObjective(objective)
    reference = loss(x)
    @test isfinite(reference) && !iszero(reference)
    gradients = map(backends) do backend
        prep = prepare_gradient(loss, backend, x)
        value, g = value_and_gradient(loss, prep, backend, x)
        @test value ≈ reference rtol = 1.0e-10
        @test all(isfinite, g)
        if reuse_index !== nothing
            nearby = copy(x)
            nearby[reuse_index] *= 1.001
            reused_value, reused = value_and_gradient(loss, prep, backend, nearby)
            expected_value, expected = value_and_gradient(loss, backend, nearby)
            @test reused_value ≈ expected_value rtol = 1.0e-10
            @test all(isfinite, reused)
            @test all(isapprox.(relative(reused, reference), relative(expected, reference); rtol, atol))
        end
        g
    end
    for g in gradients[2:end]
        @test all(isapprox.(relative(gradients[1], reference), relative(g, reference); rtol, atol))
    end
    if finite_difference
        fd = gradient(loss, finite_difference_backend(), x)
        @test all(isfinite, fd)
        @test all(isapprox.(relative(fd, reference), relative(gradients[1], reference); rtol = fd_rtol, atol = 1.0e-7))
    end
    return gradients[1]
end

"""
    test_snr_gradients(strain, metric, x; kwargs...)

[`test_gradients`](@ref) for the optimal SNR of `strain(x)` (a `(frequency, channel)` array) in
`metric`. The last three coordinates must be distance, inclination and phase with distance scaled to one,
so the distance derivative of `rho` is exactly `-rho` and prepared gradients are reused at a nearby distance.
"""
function test_snr_gradients(strain, metric, x::AbstractVector; kwargs...)
    loss = FixedObjective(z -> sqrt(optimal_snr_squared(strain(z), metric)))
    g = test_gradients(loss, x; reuse_index = length(x) - 2, kwargs...)
    @test g[end - 2] / loss(x) ≈ -1 rtol = get(kwargs, :rtol, 1.0e-6)
    return g
end

"""
    test_phase_gradients(strain, x; backends, rtol)

A phase-sensitive overlap with a phase-shifted copy of `strain(x)` checks the reference-phase
derivative (the last coordinate) that a norm such as the SNR cannot see.
"""
function test_phase_gradients(strain, x::AbstractVector; backends, rtol, delta = 0.2)
    shifted = copy(x)
    shifted[end] += delta
    template = conj.(strain(shifted))
    scale = sum(abs2, template)
    loss = FixedObjective(z -> sum(real, strain(z) .* template) / scale)
    gradients = map(backend -> gradient(loss, backend, x), backends)
    foreach(g -> @test(all(isfinite, g)), gradients)
    for g in gradients[2:end]
        @test all(isapprox.(gradients[1], g; rtol, atol = 1.0e-8))
    end
    fd = gradient(loss, finite_difference_backend(), x)
    @test gradients[1][end] ≈ fd[end] rtol = 1.0e-5 atol = 1.0e-7
    @test abs(gradients[1][end]) > 1.0e-6
    return nothing
end

"""
    test_waveform_jacobian(generator, x, adapter; grid, backends, rtol, fd_rtol = 1e-5)

Jacobians of `[real; imag]` of the flattened `hphc_fd(grid, adapter(x), generator)` with respect to `x`:
the AD backends must agree with each other and with `FiniteDifferences`, componentwise relative to
the largest entry. Needs no PSD and no reference generator.
"""
function test_waveform_jacobian(generator, x::AbstractVector, adapter; grid, backends, rtol, fd_rtol = 1.0e-5)
    function flat(z)
        h = hphc_fd(grid, adapter(z), generator)
        return vcat(real.(vec(h)), imag.(vec(h)))
    end
    fixed = FixedObjective(flat)
    jacobians = map(backend -> jacobian(fixed, backend, x), backends)
    scale = maximum(abs, jacobians[1])
    @test all(isfinite, jacobians[1]) && scale > 0
    for J in jacobians[2:end]
        @test maximum(abs, J .- jacobians[1]) <= rtol * scale
    end
    fd = jacobian(fixed, finite_difference_backend(), x)
    @test maximum(abs, fd .- jacobians[1]) <= fd_rtol * scale
    return jacobians[1]
end
