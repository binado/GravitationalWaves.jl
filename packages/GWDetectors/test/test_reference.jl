using JSON

@testset "LAL and Bilby fixtures" begin
    path = joinpath(@__DIR__, "../../../reference/detectors.json")
    if !isfile(path)
        if get(ENV, "CI", "false") == "true"
            error("Missing detector fixtures: run python/generate_detectors.py")
        end
        @warn "Skipping detector references; run uv run --project python python python/generate_detectors.py"
    else
        fixtures = JSON.parsefile(path)
        rows(x) = Tuple(Tuple(row) for row in x)
        for g in fixtures["detectors"]
            d = Detector(g["name"])
            @test d.vertex == Tuple(g["vertex"])
            @test d.xarm == Tuple(g["xarm"])
            @test d.yarm == Tuple(g["yarm"])
            @test d.arm_lengths == Tuple(g["arm_lengths"])
            @test all(isapprox.(collect(Iterators.flatten(detector_tensor(d))), collect(Iterators.flatten(g["native_tensor"])); atol = 2.0e-16))
        end
        for g in fixtures["geodetic"]
            names = (:latitude, :longitude, :elevation, :xarm_azimuth, :yarm_azimuth, :xarm_tilt, :yarm_tilt, :arm_lengths)
            kwargs = NamedTuple{names}(map(k -> g[String(k)], names))
            d = Detector("geodetic"; kwargs...)
            @test collect(d.vertex) ≈ g["vertex"] rtol = 2.0e-15 atol = 2.0e-9
            @test collect(d.xarm) ≈ g["xarm"] atol = 4.0e-16
            @test collect(d.yarm) ≈ g["yarm"] atol = 4.0e-16
        end
        for c in fixtures["responses"]
            args = (c["ra"], c["dec"], c["psi"], c["gmst"])
            lal = GWDetectors._antenna(args..., rows(c["rounded_tensor"]))
            native = antenna_pattern(args..., Detector(c["name"]))
            @test collect(values(lal)) ≈ c["lal"] atol = 5.0e-16
            @test collect(values(native)) ≈ c["native"] atol = 8.0e-16
            g = only(filter(g -> g["name"] == c["name"], fixtures["detectors"]))
            cached = GWDetectors._antenna(args..., rows(g["cached_tensor"]))
            @test collect(values(cached)) ≈ c["cached"] atol = 5.0e-16
        end
        for c in fixtures["gmst"]
            @test greenwich_mean_sidereal_time(c["gps"]) ≈ c["gmst"] rtol = 0 atol = 3.0e-11
        end
        complex_vector(v) = [complex(x...) for x in v]
        for c in fixtures["projections"]
            s = c["sky"]
            sky = (; ra = s["ra"], dec = s["dec"], psi = s["psi"], geocent_time = s["geocent_time"])
            d = Detector(c["name"])
            p = (; plus = complex_vector(c["plus"]), cross = complex_vector(c["cross"]))
            actual = project_fd(c["f"], p, sky, d; start_time = c["start_time"], reference_time = c["reference_time"])
            expected = complex_vector(c["strain"])
            # Bilby reduces sidereal angles internally; timing differences of
            # ~1e-13 s accumulate phase at the largest frequency.
            @test maximum(abs, actual - expected) / maximum(abs, expected) < 1.0e-9
            @test time_delay_from_geocenter(sky, d) ≈ c["delay"] atol = 2.0e-14
            a = antenna_pattern(merge(sky, (; geocent_time = c["reference_time"])), d)
            @test collect(values(a)) ≈ c["antenna"] atol = 2.0e-11
        end
    end
end
