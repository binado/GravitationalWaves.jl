using LALWaveforms, NPZ, JSON
spec = JSON.parsefile("test/python/cases.json")
f = range(0, 1024; step = 0.25)
function rep(name, a, b)
    d = abs.(a .- b)
    peak = maximum(abs, b)
    i = argmax(d)
    println("  $name: normrel=", sqrt(sum(abs2, a .- b) / sum(abs2, b)), " maxabs/peak=", maximum(d) / peak,
        " at f=", f[i], " nbad=", count(d .> 1e-12 * peak), " ident=", a == b)
end
for c in spec["cases"]
    println(c["name"])
    ref = npzread("test/reference/$(c["name"]).npz")
    ref1 = npzread("/tmp/ref1/$(c["name"]).npz")
    rep("py-vs-py(2 procs)", ref["plus"], ref1["plus"])
    params = NamedTuple(Symbol(k) => Float64(v) for (k, v) in c["params"])
    kw = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in c["settings"])
    haskey(kw, :mode_array) && (kw[:mode_array] = [(Int(m[1]), Int(m[2])) for m in kw[:mode_array]])
    s = WaveformSettings(c["approximant"]; kw...)
    h1 = hphc_fd(f, params, s); h2 = hphc_fd(f, params, s)
    rep("jl-vs-jl", h1.plus, h2.plus)
    rep("jl-vs-py", h1.plus, ref["plus"])
    rep("jl-vs-py1", h1.plus, ref1["plus"])
end
