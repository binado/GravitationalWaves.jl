using Test
using DifferentiationInterface
using Enzyme: Enzyme
using FiniteDifferences: central_fdm, grad
using ForwardDiff: ForwardDiff
using GWSurrogates
using Mooncake: Mooncake

const G = GWSurrogates

const BACKENDS = Pair{String, Any}[
    "ForwardDiff" => AutoForwardDiff(),
    "Mooncake" => AutoMooncake(; config = nothing),
    "Enzyme" => AutoEnzyme(; function_annotation = Enzyme.Const),
]

# Enzyme with runtime activity: the documented fallback when static analysis rejects a change.
const ENZYME_RUNTIME_ACTIVITY = AutoEnzyme(;
    mode = Enzyme.set_runtime_activity(Enzyme.Reverse), function_annotation = Enzyme.Const,
)

# Compare the gradient of `f` at `x` from every backend with 5th-order finite differences.
function check_gradients(f, x; backends = BACKENDS, rtol = 1.0e-6, atol = 1.0e-8)
    reference = grad(central_fdm(5, 1), f, x)[1]
    for (name, backend) in backends
        @testset "$name" begin
            g = DifferentiationInterface.gradient(f, backend, x)
            @test all(isfinite, g)
            @test isapprox(g, reference; rtol, atol)
        end
    end
    return nothing
end

@testset "GWSurrogates AD" begin
    @testset "kernels" begin
        t = [0.0, 0.3, 0.9, 1.4, 2.2, 3.0, 3.1, 4.5]
        u0 = @. sin(2t) + 0.1t^2
        # The functional kernels: the `NaturalCubicSpline` struct returns the constant knots together
        # with the differentiated data, which Enzyme's static activity analysis rejects.
        @testset "spline, data" begin
            check_gradients(u -> G.spline_eval(t, u, G.spline_second_derivatives(t, u), 1.234), u0)
        end
        @testset "spline, evaluation point" begin
            m0 = G.spline_second_derivatives(t, u0)
            check_gradients(x -> G.spline_eval(t, u0, m0, x[1]), [1.234])
        end
        @testset "spline, knots" begin
            check_gradients(tt -> G.spline_eval(tt, u0, G.spline_second_derivatives(tt, u0), 1.234), t)
        end
        @testset "Wigner-D" begin
            check_gradients(
                q -> begin
                    n = sqrt(sum(abs2, q))
                    D = G.wigner_D(3, q[1] / n, q[2] / n, q[3] / n, q[4] / n)
                    sum(real, D) + sum(imag, D .* (1:7))
                end,
                [0.3, -0.2, 0.7, 0.4],
            )
        end
        @testset "guarded norms" begin
            check_gradients(c -> sum(abs2, G.rescale_to_norm((c[1], c[2], c[3]), 0.7) .* (1.0, 2.0, 3.0)), [0.3, -0.4, 0.2])
            # At zero spin the norm has a finite (zero) gradient on every backend.
            for (name, backend) in BACKENDS
                @testset "$name zero spin" begin
                    g = DifferentiationInterface.gradient(c -> G.guarded_sqrt(sum(abs2, c)), backend, zeros(3))
                    @test g == zeros(3)
                end
            end
        end
    end

    for model in (:NRSur7dq4, :NRSur7dq4v2)
        @testset "$model" begin
            sur = G.load_surrogate(model)
            a2 = (0.3, 0.2, 0.4)
            b2 = (-0.1, 0.5, -0.3)

            @testset "coorbital modes" begin
                nt = length(sur.coorbital.t)
                s = range(0, 1; length = nt)
                base = [0.3 .* cos.(3 .* s) 0.2 .* sin.(2 .* s) 0.4 .* ones(nt)]
                check_gradients(
                    p -> begin
                        h = G.coorbital_modes(sur.coorbital, p[1], p[2] .* base, p[3] .* base; ell_max = 3)
                        sum(abs2, h) + sum(real, h[:, 5])
                    end,
                    [2.2, 1.0, 0.7],
                )
            end

            @testset "dynamics, t_ref = $t_ref" for t_ref in (nothing, -100.0)
                check_gradients(
                    p -> begin
                        out = G.dynamics(sur.dynamics, p[1], (p[2], p[3], p[4]), (p[5], 0.5, -0.3); t_ref)
                        sum(out.quat[end, :] .* (1, 2, 3, 4)) + out.orbphase[end] + sum(out.chiA[end, :] .* (1, -1, 2)) +
                            sum(out.chiB[end, :])
                    end,
                    [2.0, 0.3, 0.2, 0.4, -0.1],
                )
            end

            @testset "inertial modes and strain, fixed t_ref" begin
                objective(p) = begin
                    out = G.inertial_modes(
                        sur, p[1], (p[2], p[3], p[4]), (p[5], 0.5, -0.3);
                        ell_max = 3, t_ref = -100.0, times = range(-1000.0, 90.0; length = 40),
                    )
                    h = G.mode_sum(out.modes, 3, p[6], p[7])
                    sum(real, h) + 2 * sum(imag, h)
                end
                x = [2.2, 0.3, 0.2, 0.4, -0.1, 0.7, 0.3]
                check_gradients(objective, x)
                @testset "Enzyme runtime activity" begin
                    reference = DifferentiationInterface.gradient(objective, AutoForwardDiff(), x)
                    @test DifferentiationInterface.gradient(objective, ENZYME_RUNTIME_ACTIVITY, x) ≈ reference rtol = 1.0e-8
                end
            end

            @testset "SI strain with f_ref and dt" begin
                gen = G.PrecessingGenerator(sur; ell_max = 3, f_ref = 25.0)
                check_gradients(
                    p -> begin
                        params = (
                            total_mass = p[1], luminosity_distance = 400.0, mass_ratio = p[2],
                            spin_1 = (p[3], p[4], p[5]), spin_2 = (p[6], 0.5, -0.3), inclination = p[7], phi_ref = p[8],
                        )
                        out = G.hphc_td(gen, params; dt = 1 / 1024)
                        n = length(out.hplus)
                        acc = zero(eltype(out.hplus))
                        for k in 1:n
                            acc += out.hplus[k] * (1 + (k - 1) / (n - 1)) + 2 * out.hcross[k]
                        end
                        1.0e21 * acc
                    end,
                    [60.0, 2.2, 0.3, 0.2, 0.4, -0.1, 0.7, 0.3];
                    rtol = 1.0e-5, atol = 1.0e-5,
                )
            end
        end
    end
end
