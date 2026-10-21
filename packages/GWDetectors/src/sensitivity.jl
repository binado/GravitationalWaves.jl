"""
    sensitivities

Named tuple mapping detector names to their default bundled GWNoise curve names.
Load, for example, `GWNoise.read_psd(GWNoise.psd_path(sensitivities.H1))`.
These sensitivity choices are independent of the detector geometry.
"""
const sensitivities = let
    path = joinpath(@__DIR__, "..", "data", "sensitivity.toml")
    include_dependency(path)
    catalog = TOML.parsefile(path)["detectors"]
    names = keys(detectors)
    curves = map(names) do name
        filename = catalog[String(name)]["psd_reference"]
        Symbol(replace(chopsuffix(filename, "_psd.txt"), '-' => '_'))
    end
    NamedTuple{names}(curves)
end
