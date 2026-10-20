using Test
using AbstractGWWaveforms: AlignedSpinQuadrupolarWaveformGenerator, amplitude_phase_fd

# Shared checks for every `AbstractWaveformGenerator`. Backends `include` this file from their
# tests, with the interface (`hphc_fd`, the accessors, ...) in scope.
#
# `params` must be valid for `g`, `f` a valid grid, and the support of `g` must lie inside `f`
# with some of it strictly inside, so the masking can be checked on both sides.
function test_generator_contract(g; params, f)
    return @testset "contract: $(sprint(show, g))" begin
        @testset "accessors" begin
            @test approximant(g) isa String
            @test reference_frequency(g) isa Real
            @test minimum_frequency(g) isa Real
            @test maximum_frequency(g) isa Union{Nothing, Real}
            @test spin_support(g) in (:none, :aligned, :precessing)
            @test supports_mode_array(g) isa Bool
            @test is_differentiable(g) isa Bool
            @test occursin(approximant(g), sprint(show, g))
        end

        @testset "parameter validation" begin
            @test validate_parameters(f, params, g) === nothing
            @test_throws ArgumentError validate_parameters(f, merge(params, (; not_a_parameter = 1.0)), g)
        end

        if g isa AlignedSpinQuadrupolarWaveformGenerator
            test_amplitude_phase_fd_contract(g; params, f)
        end

        @testset "output" begin
            h = hphc_fd(f, params, g)
            @test h.plus isa AbstractVector{<:Complex}
            @test h.cross isa AbstractVector{<:Complex}
            @test length(h.plus) == length(f)
            @test length(h.cross) == length(f)

            f_min = minimum_frequency(g)
            f_max = something(maximum_frequency(g), Float64(last(f)))
            outside = [!(f_min <= fi <= f_max) for fi in f]
            inside = .!outside
            @test any(outside) && any(inside)
            for x in (h.plus, h.cross)
                @test all(iszero, x[outside])
                @test any(!iszero, x[inside])
                @test all(isfinite, x)
            end
        end
    end
end

# Additional contract for the optional aligned-spin quadrupole interface.
function test_amplitude_phase_fd_contract(g; params, f)
    return @testset "amplitude/phase contract" begin
        a = amplitude_phase_fd(f, params, g)
        @test keys(a) == (:amplitude, :phase)
        @test a.amplitude isa AbstractVector{<:Real}
        @test a.phase isa AbstractVector{<:Real}
        @test length(a.amplitude) == length(a.phase) == length(f)
        @test all(isfinite, a.amplitude) && all(isfinite, a.phase)
        @test all(>=(0), a.amplitude)
        @test iszero(a.amplitude[1]) && iszero(a.phase[1])
        @test all(iszero, a.phase[iszero.(a.amplitude)])
        f_min = minimum_frequency(g)
        f_max = something(maximum_frequency(g), Float64(last(f)))
        outside = [!(f_min <= fi <= f_max) for fi in f]
        @test all(iszero, a.amplitude[outside]) && all(iszero, a.phase[outside])
        face = amplitude_phase_fd(f, merge(params, (; theta_jn = 0.0)), g)
        @test a.amplitude == face.amplitude
        @test a.phase == face.phase
        h = hphc_fd(f, params, g)
        strain = a.amplitude .* cis.(a.phase)
        c = cos(params.theta_jn)
        @test h.plus ≈ ((1 + c * c) / 2) .* strain
        @test h.cross ≈ (-im * c) .* strain
    end
end

# Domain-specific contracts: modes in both domains have explicit, ordered signed labels.
# Backends choose alternate viewing parameters through `viewing_params` when necessary.
function test_modes_td_contract(g; params, viewing_params = merge(params, (; theta_jn = 0.2)), kwargs...)
    return @testset "TD mode contract" begin
        r = AbstractGWWaveforms.waveform_modes_td(g, params; kwargs...)
        @test r.mode_indices isa Tuple
        @test !isempty(r.mode_indices)
        @test all(label -> label isa Tuple{Integer, Integer}, r.mode_indices)
        @test all(((l, m),) -> 2 <= l <= 5 && abs(m) <= l, r.mode_indices)
        @test allunique(r.mode_indices)
        @test r.t isa AbstractVector{<:Real}
        @test r.modes isa AbstractMatrix{<:Complex}
        @test size(r.modes) == (length(r.t), length(r.mode_indices))
        @test all(isfinite, r.t) && all(isfinite, r.modes)
        other = AbstractGWWaveforms.waveform_modes_td(g, viewing_params; kwargs...)
        @test r == other
        theta, phi = AbstractGWWaveforms.mode_angles(params, g)
        h = AbstractGWWaveforms.mode_sum(r.modes, r.mode_indices, theta, phi)
        out = AbstractGWWaveforms.hphc_td(g, params; kwargs...)
        @test out.t == r.t
        @test out.hplus ≈ real.(h)
        @test out.hcross ≈ -imag.(h)
    end
end

function test_modes_fd_contract(g; params, f, viewing_params = merge(params, (; theta_jn = 0.2)))
    return @testset "FD mode contract" begin
        check_frequency_grid(f)
        r = AbstractGWWaveforms.waveform_modes_fd(f, params, g)
        @test r.mode_indices isa Tuple
        @test !isempty(r.mode_indices)
        @test all(label -> label isa Tuple{Integer, Integer}, r.mode_indices)
        @test all(((l, m),) -> 2 <= l <= 5 && abs(m) <= l, r.mode_indices)
        @test allunique(r.mode_indices)
        @test r.positive isa AbstractMatrix{<:Complex}
        @test r.negative isa AbstractMatrix{<:Complex}
        @test size(r.positive) == size(r.negative) == (length(f), length(r.mode_indices))
        @test r.positive[1, :] == r.negative[1, :]
        lo, hi = AbstractGWWaveforms.frequency_bounds(g, f)
        outside = [!(lo <= fi <= hi) for fi in f]
        for branch in (r.positive, r.negative)
            @test all(isfinite, branch)
            @test all(iszero, branch[outside, :])
        end
        @test r == AbstractGWWaveforms.waveform_modes_fd(f, viewing_params, g)
    end
end
