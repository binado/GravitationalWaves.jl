using Reactant
using Enzyme: Enzyme

Reactant.set_default_backend("cpu")

_host(x::AbstractArray) = Array(x)
_host(x::Union{Tuple, NamedTuple}) = map(_host, x)

# Relative to the largest sample: the waveform spans many orders of magnitude.
function _close(actual, expected; rtol = 1.0e-11)
    return @test maximum(abs, actual .- expected) <= rtol * maximum(abs, expected)
end

# One executable must serve several parameter draws, with only the scalars tracked.
const RX_GRID = range(0, 256; step = 0.25)
const RX_BASE = (
    mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 450.0, theta_jn = 0.7, phase = 1.1,
    a_1 = 0.3, a_2 = 0.2, tilt_1 = 0.0, tilt_2 = Float64(π),
)
const RX_DRAWS = (
    RX_BASE,
    merge(RX_BASE, (; mass_1 = 12.0, mass_2 = 31.0, a_1 = 0.9, tilt_1 = Float64(π), tilt_2 = 0.0)),
    merge(RX_BASE, (; mass_1 = 60.0, mass_2 = 55.0, theta_jn = 2.4, phase = 0.2, luminosity_distance = 1200.0)),
)

@testset "Reactant CPU" begin
    @testset "TaylorF2 hphc_fd" begin
        for kw in (
                (;),
                (; minimum_frequency = 25.0, maximum_frequency = 180.0),
                (; minimum_frequency = 20.0, pn_amplitude_order = 3),
                (; minimum_frequency = 30.0, pn_phase_order = 4, pn_spin_order = 3, reference_frequency = 0.0),
            )
            @testset "$kw" begin
                g = TaylorF2(; kw...)
                tracked = Reactant.to_rarray(first(RX_DRAWS); track_numbers = true)
                executable = Reactant.@compile hphc_fd(RX_GRID, tracked, g)
                for params in RX_DRAWS
                    validate_parameters(RX_GRID, params, g)
                    actual = _host(executable(RX_GRID, Reactant.to_rarray(params; track_numbers = true), g))
                    expected = hphc_fd(RX_GRID, params, g)
                    _close(actual.plus, expected.plus)
                    _close(actual.cross, expected.cross)
                    @test all(iszero, actual.plus[RX_GRID .< minimum_frequency(g)])
                end
            end
        end

        @testset "other mass parameters" begin
            g = TaylorF2(; minimum_frequency = 25.0)
            m1, m2 = RX_BASE.mass_1, RX_BASE.mass_2
            alt = (;
                chirp_mass = (m1 * m2)^(3 / 5) / (m1 + m2)^(1 / 5), mass_ratio = m2 / m1,
                luminosity_distance = 450.0, theta_jn = 0.7, phase = 1.1,
            )
            tracked = Reactant.to_rarray(alt; track_numbers = true)
            executable = Reactant.@compile hphc_fd(RX_GRID, tracked, g)
            _close(_host(executable(RX_GRID, tracked, g)).plus, hphc_fd(RX_GRID, alt, g).plus)
        end

        @testset "the ISCO cutoff is not supported" begin
            g = TaylorF2(; maximum_frequency = 0)
            tracked = Reactant.to_rarray(RX_BASE; track_numbers = true)
            @test_throws ArgumentError Reactant.@compile hphc_fd(RX_GRID, tracked, g)
        end
    end

    @testset "IMRPhenomD hphc_fd" begin
        # `fCut = 0.2 / M` falls inside `RX_GRID` (256 Hz) for M ≳ 160 M☉, and the same
        # executable must serve draws on both sides of it.
        draws = (
            RX_DRAWS...,
            merge(RX_BASE, (; mass_1 = 130.0, mass_2 = 100.0, a_1 = 0.7, a_2 = 0.5)),
            merge(RX_BASE, (; mass_1 = 90.0, mass_2 = 95.0, a_1 = 0.9, tilt_1 = Float64(π), tilt_2 = Float64(π))),
            merge(RX_BASE, (; mass_1 = 200.0, mass_2 = 30.0, a_1 = 0.0, a_2 = 0.9, tilt_2 = 0.0)),
        )
        for kw in (
                (;),
                (; minimum_frequency = 25.0, maximum_frequency = 180.0),
                (; minimum_frequency = 20.1, maximum_frequency = 200.3, reference_frequency = 0.0),
                (; minimum_frequency = 30.0, maximum_frequency = 140.0),
            )
            @testset "$kw" begin
                g = IMRPhenomD(; kw...)
                tracked = Reactant.to_rarray(first(draws); track_numbers = true)
                executable = Reactant.@compile hphc_fd(RX_GRID, tracked, g)
                for params in draws
                    validate_parameters(RX_GRID, params, g)
                    actual = _host(executable(RX_GRID, Reactant.to_rarray(params; track_numbers = true), g))
                    expected = hphc_fd(RX_GRID, params, g)
                    _close(actual.plus, expected.plus)
                    _close(actual.cross, expected.cross)
                    # The mask is the same on both sides, including the zero at the cutoff.
                    @test iszero.(actual.plus) == iszero.(expected.plus)
                    @test all(iszero, actual.plus[RX_GRID .< minimum_frequency(g)])
                end
            end
        end

        @testset "the cutoff depending on the masses is not supported" begin
            g = IMRPhenomD(; maximum_frequency = 0)
            tracked = Reactant.to_rarray(RX_BASE; track_numbers = true)
            @test_throws ArgumentError Reactant.@compile hphc_fd(RX_GRID, tracked, g)
        end
    end

    @testset "TaylorF2 Enzyme gradient" begin
        generator = TaylorF2(; minimum_frequency = 25.0)
        x_template = AD_X0 .* [1.01, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        template = 1.0e22 .* hphc_fd(AD_GRID, ad_params(x_template), generator).cross
        loss = ad_loss(generator, template)
        objective(x) = Reactant.@allowscalar loss(x)

        executable = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, Reactant.to_rarray(AD_X0))
        for x in (AD_X0, AD_X0 .* [1.05, 0.9, 1.0, 1.0, 1.1, 0.8, 0.5])
            actual = Array(only(executable(Enzyme.Reverse, objective, Reactant.to_rarray(x))))
            expected = ForwardDiff.gradient(loss, x)
            @test all(isfinite, actual)
            @test isapprox(actual, expected; rtol = 1.0e-8)
        end
    end

    @testset "IMRPhenomD Enzyme gradient" begin
        generator = IMRPhenomD(; minimum_frequency = 25.0)
        # The second start has `fCut` (about 195 Hz) inside `AD_GRID`: the mask moves with the masses.
        heavy = [90.0, 0.8, 0.3, 0.2, 450.0, 0.7, 1.1]
        for x0 in (AD_X0, heavy)
            x_template = x0 .* [1.01, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
            template = 1.0e22 .* hphc_fd(AD_GRID, ad_params(x_template), generator).cross
            loss = ad_loss(generator, template)
            objective(x) = Reactant.@allowscalar loss(x)

            executable = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, Reactant.to_rarray(x0))
            for x in (x0, x0 .* [1.05, 0.9, 1.0, 1.0, 1.1, 0.8, 0.5])
                actual = Array(only(executable(Enzyme.Reverse, objective, Reactant.to_rarray(x))))
                expected = ForwardDiff.gradient(loss, x)
                @test all(isfinite, actual)
                @test isapprox(actual, expected; rtol = 1.0e-6)
            end
        end
    end
end
