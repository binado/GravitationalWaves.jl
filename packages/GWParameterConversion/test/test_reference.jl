using JSON

@testset "reference (bilby)" begin
    file = joinpath(@__DIR__, "..", "..", "..", "reference", "conversion.json")
    if !isfile(file)
        msg = "bilby conversion fixtures missing ($file); " *
            "run `uv run --project python python python/generate.py`"
        get(ENV, "CI", "") == "true" ? error(msg) : @warn msg
    else
        fixtures = JSON.parsefile(file)
        cases = fixtures["masses"]
        # Julia and Python can bind different libm `pow`s, so allow an ulp-level difference.
        rtol = 1.0e-14
        for c in cases
            @testset "m1=$(c["mass_1"]) m2=$(c["mass_2"])" begin
                m1, m2 = c["mass_1"], c["mass_2"]
                mc, mt, q, eta = c["chirp_mass"], c["total_mass"], c["mass_ratio"], c["symmetric_mass_ratio"]
                @test isapprox(chirp_mass(m1, m2), mc; rtol)
                @test total_mass(m1, m2) == mt
                @test mass_ratio(m1, m2) == q
                @test isapprox(symmetric_mass_ratio(m1, m2), eta; rtol)
                @test all(isapprox.(total_mass_and_mass_ratio_to_component_masses(mt, q), c["from_total_mass_mass_ratio"]; rtol))
                @test all(isapprox.(chirp_mass_and_mass_ratio_to_component_masses(mc, q), c["from_chirp_mass_mass_ratio"]; rtol))
                @test isapprox(chirp_mass_and_mass_ratio_to_total_mass(mc, q), c["total_mass_from_chirp_mass_mass_ratio"]; rtol)
                @test isapprox(chirp_mass_and_total_mass_to_symmetric_mass_ratio(mc, mt), c["eta_from_chirp_mass_total_mass"]; rtol)
                @test isapprox(symmetric_mass_ratio_to_mass_ratio(eta), c["mass_ratio_from_eta"]; rtol)
                @test isapprox(mass_1_and_chirp_mass_to_mass_ratio(m1, mc), c["mass_ratio_from_mass_1_chirp_mass"]; rtol)
                expected = c["mass_ratio_from_mass_2_chirp_mass"]
                if expected === nothing
                    @test_throws DomainError mass_2_and_chirp_mass_to_mass_ratio(m2, mc)
                else
                    @test isapprox(mass_2_and_chirp_mass_to_mass_ratio(m2, mc), expected; rtol)
                end
            end
        end
        @testset "spin frame" begin
            for c in fixtures["spins"]
                p = NamedTuple(Symbol(k) => v for (k, v) in c if k in ("mass_1", "mass_2", "theta_jn", "phase", "a_1", "a_2", "tilt_1", "tilt_2", "phi_12", "phi_jl"))
                s = spins_to_lal_frame(p, c["f_ref"])
                # bilby runs LAL's own transform; the port should agree to rounding.
                @test collect(values(s)) ≈ c["spins"] atol = 1.0e-13
            end
        end
    end
end
