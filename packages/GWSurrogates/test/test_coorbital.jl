using FiniteDifferences
using ForwardDiff
using NPZ

const FIXTURE_DIR = joinpath(@__DIR__, "..", "..", "..", "reference", "surrogates")
const MODEL_PARAMS = (NRSur7dq4 = GWSurrogates.NRSur7dq4Params(), NRSur7dq4v2 = GWSurrogates.NRSur7dq4v2Params())

@testset "fits" begin
    s = GWSurrogates.FitSettings(-0.5, 2.0, 3, 2)
    @test GWSurrogates.n_fit_powers(s) == 4 + 6 * 3

    @testset "power table and evaluation" begin
        x = (0.3, 0.1, -0.2, 0.4, 0.5, -0.6, 0.7)
        powers = GWSurrogates.fit_powers!(zeros(GWSurrogates.n_fit_powers(s)), s, x)
        qb = -0.5 + 2.0 * 0.3
        @test powers[1:4] ≈ [1.0, qb, qb^2, qb^3]
        @test powers[5:7] ≈ [1.0, 0.1, 0.01]
        @test powers[(end - 2):end] ≈ [1.0, 0.7, 0.49]
        # Two fits: 1.5 * q^2 * x2 (plus a zero term) and 2.0 - 0.5 * q * x3^2 * x6.
        o1 = [2 0; 1 0; 0 0; 0 0; 0 0; 0 0; 0 0]
        o2 = [0 0; 0 0; 0 2; 0 0; 0 0; 0 0; 0 0]
        o2[1, 2] = 1
        o2[2, 2] = 0
        o2[6, 2] = 1
        fits = GWSurrogates.FitSet(s, [o1, o2], [[1.5, 0.0], [2.0, -0.5]])
        @test length(fits) == 2
        @test GWSurrogates.eval_fit(fits, 1, powers) ≈ 1.5 * qb^2 * 0.1
        @test GWSurrogates.eval_fit(fits, 2, powers) ≈ 2.0 - 0.5 * qb * (-0.2)^2 * (-0.6)
    end

    @testset "argument checks" begin
        @test_throws ArgumentError GWSurrogates.FitSet(s, [fill(4, 7, 1)], [[1.0]])
        @test_throws DimensionMismatch GWSurrogates.FitSet(s, [zeros(Int, 7, 2)], [[1.0]])
    end

    @testset "fit parameters" begin
        q = 2.0
        a = (0.1, 0.2, 0.3)
        b = (-0.4, 0.5, -0.6)
        hat = GWSurrogates.chi_hat(q, 0.3, -0.6)
        @test GWSurrogates.fit_params(GWSurrogates.NRSur7dq4Params(), q, a, b) ==
            (log(q), 0.1, 0.2, hat, -0.4, 0.5, (0.3 + 0.6) / 2)
        @test GWSurrogates.fit_params(GWSurrogates.NRSur7dq4v2Params(), q, a, b) ==
            (log(q), (0.1 - 0.4) / 2, (0.1 + 0.4) / 2, hat, (0.2 + 0.5) / 2, (0.2 - 0.5) / 2, 0.3)
        # Equal masses and aligned spins: chiHat reduces to the mean spin.
        @test GWSurrogates.chi_hat(1.0, 0.5, 0.5) ≈ (0.5 - 38 * 0.25 / 113 * 1.0) / (1 - 76 * 0.25 / 113)
    end
end

@testset "data_path" begin
    @test_throws ArgumentError GWSurrogates.data_path(:NotAModel)
    mktempdir() do dir
        withenv("GWSURROGATES_DATA_DIR" => dir) do
            @test_throws ErrorException GWSurrogates.data_path(:NRSur7dq4)
            touch(joinpath(dir, "NRSur7dq4.h5"))
            @test GWSurrogates.data_path(:NRSur7dq4) == joinpath(dir, "NRSur7dq4.h5")
        end
    end
end

@testset "coorbital surrogate" begin
    ci = get(ENV, "CI", "false") == "true"
    for (model, params) in pairs(MODEL_PARAMS)
        fixture = joinpath(FIXTURE_DIR, "$(model)_coorbital.npz")
        if !isfile(fixture)
            msg = "coorbital reference fixture $fixture not found; generate it with python/generate_surrogates.py"
            ci ? error(msg) : @warn msg
            continue
        end
        @testset "$model" begin
            ref = npzread(fixture)
            sur = GWSurrogates.load_coorbital(model)
            @test sur.params == params
            @test sur.t == ref["t_coorb"]
            @test sur.ell_max == ref["ell_max"]
            for i in 0:(ref["n_cases"] - 1)
                q = ref["case$(i)_q"]
                chiA = ref["case$(i)_chiA"]
                chiB = ref["case$(i)_chiB"]
                expected = permutedims(ref["case$(i)_modes"])
                modes = GWSurrogates.coorbital_modes(sur, q, chiA, chiB)
                @test size(modes) == size(expected)
                scale = maximum(abs, expected)
                @test maximum(abs, modes - expected) < 1.0e-12 * scale
            end
            # Lower ell_max returns the leading modes.
            q, chiA, chiB = ref["case1_q"], ref["case1_chiA"], ref["case1_chiB"]
            low = GWSurrogates.coorbital_modes(sur, q, chiA, chiB; ell_max = 3)
            @test low == GWSurrogates.coorbital_modes(sur, q, chiA, chiB)[:, 1:size(low, 2)]
            @test_throws ArgumentError GWSurrogates.coorbital_modes(sur, q, chiA, chiB; ell_max = 1)
            @test_throws DimensionMismatch GWSurrogates.coorbital_modes(sur, q, chiA[1:10, :], chiB)

            # Gradients through the fits and the empirical-interpolation contraction.
            fdm = central_fdm(5, 1)
            objective(q) = sum(abs2, GWSurrogates.coorbital_modes(sur, q, chiA, chiB; ell_max = 2))
            @test ForwardDiff.derivative(objective, q) ≈ fdm(objective, q) rtol = 1.0e-7
            scaled(a) = sum(abs2, GWSurrogates.coorbital_modes(sur, q, a .* chiA, a .* chiB; ell_max = 3))
            @test ForwardDiff.derivative(scaled, 1.0) ≈ fdm(scaled, 1.0) rtol = 1.0e-7
            @test (@inferred GWSurrogates.coorbital_modes(sur, q, chiA, chiB; ell_max = 2)) isa Matrix{ComplexF64}
        end
    end
end
