const _PSD_DIRECTORY = joinpath(@__DIR__, "..", "data", "noise_curves")
const _PSD_FILES = (
    AplusDesign = "AplusDesign_psd.txt",
    ET_000A_18_ETD_hlf = "ET-000A-18_ETD_hlf_psd.txt",
    ET_COBA_10km = "ET_COBA_10km_psd.txt",
    ET_COBA_15km = "ET_COBA_15km_psd.txt",
    ET_COBA_20km = "ET_COBA_20km_psd.txt",
    avirgo_O5high_NEW = "avirgo_O5high_NEW_psd.txt",
    avirgo_O5low_NEW = "avirgo_O5low_NEW_psd.txt",
    cosmic_explorer_20km_pm_strain = "cosmic_explorer_20km_pm_strain_psd.txt",
    cosmic_explorer_20km_strain = "cosmic_explorer_20km_strain_psd.txt",
    cosmic_explorer_40km_lf_strain = "cosmic_explorer_40km_lf_strain_psd.txt",
    cosmic_explorer_strain = "cosmic_explorer_strain_psd.txt",
    kagra_128Mpc = "kagra_128Mpc_psd.txt",
)
"""
    available_psds()

Return the names of all bundled one-sided PSD curves, including alternative designs.
"""
available_psds() = keys(_PSD_FILES)

"""
    psd_path(name::Symbol)

Return the path to a bundled two-column frequency (Hz), strain PSD (strain²/Hz)
text file. `name` is a curve name from `available_psds()`.
"""
function psd_path(name::Symbol)
    haskey(_PSD_FILES, name) || throw(ArgumentError("unknown PSD curve: $name"))
    return joinpath(_PSD_DIRECTORY, _PSD_FILES[name])
end
psd_path(name::AbstractString) = psd_path(Symbol(name))
