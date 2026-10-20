using ForwardDiff
using LinearAlgebra

@testset "IMRPhenomXAS internals and boundaries" begin
    @testset "small pivoted solves" begin
        for n in (3, 4, 5), shift in 0:(n - 1)
            A = [sin(i * j) + (i == j ? 3.0 : 0.0) for i in 1:n, j in 1:n]
            A = A[circshift(collect(1:n), shift), :]
            b = cos.(collect(1:n))
            rows = ntuple(i -> Tuple(A[i, :]), n)
            c = collect(GWWaveforms.xas_solve(rows, Tuple(b)))
            @test c ≈ A \ b rtol = 1.0e-12
            @test norm(A * c - b) < 1.0e-12
        end
    end

    @testset "ISCO at zero final spin" begin
        isco = GWWaveforms.xas_isco
        @test isco(0.0) ≈ 1 / (6sqrt(6.0) * Float64(π))
        @test ForwardDiff.derivative(isco, 0.0) ≈ 11 / (216Float64(π)) rtol = 1.0e-12
        for a in (-1.0e-7, 0.0, 1.0e-7, -0.7, 0.7)
            @test isfinite(ForwardDiff.derivative(isco, a))
        end
    end

    @testset "analytic frequency derivatives and C1 matching" begin
        for (m1, m2, c1, c2) in ((36.0, 29.0, 0.3, -0.2), (30.0, 30.0, 0.6, -0.3), (100.0, 1.0, 0.2, -0.1))
            m = GWWaveforms.phenomxas_model(m1, m2, c1, c2, 50.0, 0.2)
            p = m.phase
            for (fun, dfun, coeff, f) in (
                    (GWWaveforms.xas_inspiral_phase, GWWaveforms.xas_dinspiral_phase, p.ins, p.f1),
                    (GWWaveforms.xas_int_phase, GWWaveforms.xas_dint_phase, p.int, p.f2),
                    (GWWaveforms.xas_rd_phase, GWWaveforms.xas_drd_phase, p.rd, m.r.fRD),
                    (GWWaveforms.xas_inspiral_amp, GWWaveforms.xas_dinspiral_amp, m.amplitude.ins, m.amplitude.f1),
                    (GWWaveforms.xas_rd_amp, GWWaveforms.xas_drd_amp, m.amplitude.rd, m.amplitude.f4),
                )
                @test dfun(coeff, f) ≈ ForwardDiff.derivative(x -> fun(coeff, x), f) rtol = 1.0e-10 atol = 1.0e-12
            end
            @test GWWaveforms.xas_inspiral_phase(p.ins, p.f1) ≈ GWWaveforms.xas_int_phase(p.int, p.f1) + p.alpha0 + p.alpha1 * p.f1 rtol = 1.0e-11
            @test GWWaveforms.xas_dinspiral_phase(p.ins, p.f1) ≈ GWWaveforms.xas_dint_phase(p.int, p.f1) + p.alpha1 rtol = 1.0e-11
            @test GWWaveforms.xas_int_phase(p.int, p.f2) + p.alpha0 + p.alpha1 * p.f2 ≈ GWWaveforms.xas_rd_phase(p.rd, p.f2) + p.beta0 + p.beta1 * p.f2 rtol = 1.0e-11
            @test GWWaveforms.xas_dint_phase(p.int, p.f2) + p.alpha1 ≈ GWWaveforms.xas_drd_phase(p.rd, p.f2) + p.beta1 rtol = 1.0e-11
        end
    end

    @testset "mass ordering and equal-mass derivatives" begin
        f = range(0, 512; step = 0.5)
        p = (
            mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.7, phase = 1.1,
            a_1 = 0.3, a_2 = 0.2, tilt_1 = 0.0, tilt_2 = Float64(π),
        )
        @test_throws ArgumentError validate_parameters(
            f, merge(p, (mass_1 = p.mass_2, mass_2 = p.mass_1, a_1 = p.a_2, a_2 = p.a_1, tilt_1 = p.tilt_2, tilt_2 = p.tilt_1)),
            IMRPhenomXAS(),
        )
        for equal in (false, true)
            x = [36.0, equal ? 36.0 : 29.0, 0.3, 0.2]
            loss(x) = sum(real, hphc_fd(f, merge(p, (mass_1 = x[1], mass_2 = x[2], a_1 = x[3], a_2 = x[4])), IMRPhenomXAS()).plus) * 1.0e22
            g = ForwardDiff.gradient(loss, x)
            fd = [(loss(x .+ (1:4 .== i) .* 1.0e-5) - loss(x .- (1:4 .== i) .* 1.0e-5)) / 2.0e-5 for i in 1:4]
            @test all(isfinite, g)
            @test g ≈ fd rtol = 2.0e-5
        end
    end

    @testset "inclusive cutoff and exact frequency limits" begin
        p = (mass_1 = 180.0, mass_2 = 120.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 0.2)
        cutoff = 0.3 / ((p.mass_1 + p.mass_2) * GWWaveforms.MTSUN_SI)
        f = range(0; step = cutoff / 64, length = 130)
        h = hphc_fd(f, p, IMRPhenomXAS())
        @test !iszero(h.plus[65])
        @test all(iszero, h.plus[66:end])
        @test_throws ArgumentError IMRPhenomXAS(pn_amplitude_order = 1)
        @test_throws MethodError IMRPhenomXAS(mode_array = [(2, 2)])
        @test validate_parameters(f, p, IMRPhenomXAS()) === nothing
        @test_throws ArgumentError validate_parameters(f, merge(p, (mass_1 = 1001.0, mass_2 = 1.0)), IMRPhenomXAS())
    end
end
