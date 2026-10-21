using JSON
using NPZ

@testset "reference (bilby)" begin
    LW = LALWaveforms
    root = normpath(@__DIR__, "..", "..", "..")
    spec = JSON.parsefile(joinpath(root, "python", "cases.json"))
    refdir = joinpath(root, "reference")
    df = spec["frequency"]["delta_f"]
    fmax = spec["frequency"]["f_max"]
    f = range(0, fmax; step = df)
    # Julia loads libopenlibm globally, so on Linux LAL's libm calls bind to it rather than
    # glibc. Results are bit-identical to Python on macOS but differ at the ulp level on
    # Linux (amplified by sensitive models), hence the optional per-case `rtol`.

    missing_names = [
        c["name"] for c in spec["cases"] if !isfile(joinpath(refdir, c["name"] * ".npz"))
    ]
    if !isempty(missing_names)
        msg = "bilby reference fixtures missing in $refdir ($(join(missing_names, ", "))); " *
            "run `uv run --project python python python/generate.py`"
        if get(ENV, "CI", "") == "true"
            error(msg)
        else
            @warn msg
        end
    end

    for c in spec["cases"]
        file = joinpath(refdir, c["name"] * ".npz")
        isfile(file) || continue
        @testset "$(c["name"])" begin
            ref = npzread(file)
            # `julia_params` lets a case feed Julia a different (equivalent) parameterisation than bilby.
            params = NamedTuple(
                Symbol(k) => Float64(v) for (k, v) in get(c, "julia_params", c["params"])
            )
            kw = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in c["settings"])
            if haskey(kw, :mode_array)
                kw[:mode_array] = [(Int(m[1]), Int(m[2])) for m in kw[:mode_array]]
            end
            settings = LALWaveformGenerator(c["approximant"]; kw...)

            h = hphc_fd(f, params, settings)
            rtol = Float64(get(c, "rtol", spec["default_rtol"]))
            @test isapprox(h[:, :plus], ref["plus"]; rtol)
            @test isapprox(h[:, :cross], ref["cross"]; rtol)

            fref = get(kw, :reference_frequency, 50.0)
            s = LW.bilby_to_lalsimulation_spins(LW.bilby_parameters(params), fref)
            @test collect(values(s)) ≈ ref["spins"] atol = 1.0e-14
        end
    end
end
