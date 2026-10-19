using DifferentiationInterface
using FiniteDifferences: central_fdm, grad
using ForwardDiff: ForwardDiff

# Shared gradient checks for the default ForwardDiff tests and the opt-in integration suite.

const AD_GRID = range(0, 256; step = 0.5)
const AD_X0 = [28.0, 0.8, 0.3, 0.2, 450.0, 0.7, 1.1]

# (chirp_mass, mass_ratio, a_1, a_2, luminosity_distance, theta_jn, phase), spin 1 aligned
# and spin 2 anti-aligned with the orbital angular momentum.
function ad_params(x)
    return (
        chirp_mass = x[1], mass_ratio = x[2], a_1 = x[3], a_2 = x[4], luminosity_distance = x[5],
        theta_jn = x[6], phase = x[7], tilt_1 = 0.0, tilt_2 = Float64(π),
    )
end

# A real scalar that keeps the phase in the gradient: the norm of the plus polarisation plus
# the overlap of the cross polarisation with a fixed template. Rescaled to order one.
function ad_loss(generator, template)
    return function (x)
        h = hphc_fd(AD_GRID, ad_params(x), generator)
        scale = 1.0e22
        return scale^2 * (sum(abs2, h.plus) + sum(real, h.cross .* conj(template)))
    end
end

function test_waveform_gradients(backends; generators = (TaylorF2, IMRPhenomD, IMRPhenomXAS))
    return @testset "automatic differentiation" begin
        x_template = AD_X0 .* [1.01, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        for G in generators
            @testset "$G" begin
                generator = G(; minimum_frequency = 25.0)
                template = 1.0e22 .* hphc_fd(AD_GRID, ad_params(x_template), generator).cross
                loss(x) = ad_loss(generator, template)(x)
                fd = central_fdm(5, 1)
                g_fd = grad(fd, loss, AD_X0)[1]
                @test all(isfinite, g_fd)

                gradients = Dict{String, Vector{Float64}}()
                for (name, backend) in backends
                    @testset "$name" begin
                        g = DifferentiationInterface.gradient(loss, backend, AD_X0)
                        gradients[name] = g
                        @test all(isfinite, g)
                        @test isapprox(g, g_fd; rtol = 1.0e-5)
                    end
                end
                @testset "backends agree" begin
                    reference = gradients["ForwardDiff"]
                    for (name, g) in gradients
                        @test isapprox(g, reference; rtol = 1.0e-10)
                    end
                end
                if G === IMRPhenomD
                    @testset "equal-mass gradients" begin
                        function component_loss(x)
                            p = (;
                                mass_1 = x[1], mass_2 = x[2], a_1 = x[3], a_2 = x[4],
                                luminosity_distance = 500.0, theta_jn = 0.4, phase = 0.2,
                                tilt_1 = 0.0, tilt_2 = Float64(π),
                            )
                            h = hphc_fd(AD_GRID, p, generator)
                            return 1.0e22 * (sum(real, h.plus) + sum(imag, h.cross))
                        end
                        for mass in (30.0, 60.0), spins in ((0.0, 0.0), (0.3, 0.3), (0.6, 0.2))
                            x = [mass, mass, spins...]
                            expected = grad(central_fdm(5, 1), component_loss, x)[1]
                            @test all(isfinite, expected)
                            for (name, backend) in backends
                                @testset "$name, mass=$mass, spins=$spins" begin
                                    g = DifferentiationInterface.gradient(component_loss, backend, x)
                                    @test all(isfinite, g)
                                    @test isapprox(g, expected; rtol = 1.0e-5)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
