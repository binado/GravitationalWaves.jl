using Downloads: Downloads
using HDF5: HDF5
using Scratch: @get_scratch!
using SHA: sha256

const ZENODO_RECORD = "https://zenodo.org/records/22257361/files"

# File name, SHA-256 and size of the model data (CC-BY-4.0, doi:10.5281/zenodo.22257361).
const MODEL_FILES = (
    NRSur7dq4 = (
        file = "NRSur7dq4.h5",
        sha256 = "3fc2f58f6ed7a4fff58856b6c01a7dc9155954d2b9fa5d58126a7807826d45d5",
        params = NRSur7dq4Params(),
    ),
    NRSur7dq4v2 = (
        file = "NRSur7dq4v2.h5",
        sha256 = "62f3afb2baeb837ea19692ff6ca8737808e336cfbd9881c474147ae133efc9de",
        params = NRSur7dq4v2Params(),
    ),
)

"""
    data_path(model::Symbol; download = true) -> String

Path of the HDF5 file of `model` (`:NRSur7dq4` or `:NRSur7dq4v2`). Looks in the directory named
by the environment variable `GWSURROGATES_DATA_DIR` first; otherwise the file is cached in a
scratch space and, when missing and `download` is true, fetched from Zenodo and checked against
its SHA-256. Throws when the file is not available.
"""
function data_path(model::Symbol; download::Bool = true)
    haskey(MODEL_FILES, model) || throw(ArgumentError("unknown model $model; choose from $(keys(MODEL_FILES))"))
    info = MODEL_FILES[model]
    override = get(ENV, "GWSURROGATES_DATA_DIR", "")
    if !isempty(override)
        path = joinpath(override, info.file)
        isfile(path) || error("$(info.file) not found in GWSURROGATES_DATA_DIR=$override")
        return path
    end
    path = joinpath(@get_scratch!("model-data"), info.file)
    if !isfile(path)
        download || error("$(info.file) is not cached; call data_path($(repr(model))) with download = true")
        tmp = path * ".part"
        Downloads.download("$ZENODO_RECORD/$(info.file)", tmp)
        checksum = bytes2hex(open(sha256, tmp))
        if checksum != info.sha256
            rm(tmp)
            error("checksum mismatch for $(info.file): expected $(info.sha256), got $checksum")
        end
        mv(tmp, path)
    end
    return path
end

function read_datapiece(group, settings::FitSettings)
    node_indices = Int.(HDF5.read_dataset(group, "nodeIndices")) .+ 1
    basis = Matrix{Float64}(HDF5.read_dataset(group, "EIBasis"))
    nodes = HDF5.open_group(group, "nodeModelers")
    n = length(node_indices)
    orders = [Matrix{Int}(HDF5.read_dataset(nodes, "bfOrders_$i")) for i in 0:(n - 1)]
    coefs = [Vector{Float64}(HDF5.read_dataset(nodes, "coefs_$i")) for i in 0:(n - 1)]
    close(nodes)
    return CoorbitalComponent(node_indices, FitSet(settings, orders, coefs), basis)
end

"""
    load_coorbital(path, params; settings = NRSUR7DQ4_FIT_SETTINGS) -> CoorbitalSurrogate

Read the coorbital waveform surrogate of an NRSur7dq4 (`params = NRSur7dq4Params()`) or
NRSur7dq4v2 (`NRSur7dq4v2Params()`) HDF5 file. The file layout (one or two time subdomains) is
detected from the group names.
"""
function load_coorbital(path::AbstractString, params; settings::FitSettings = NRSUR7DQ4_FIT_SETTINGS)
    return HDF5.h5open(path, "r") do f
        multidomain = haskey(f, "hCoorb_2_2_Re+_subdomain_0")
        nsub = multidomain ? 2 : 1
        suffix(sub) = multidomain ? "_subdomain_$(sub - 1)" : ""
        ell_max = 2
        while haskey(f, "hCoorb_$(ell_max + 1)_$(ell_max + 1)_Re+" * suffix(1))
            ell_max += 1
        end
        t = Vector{Float64}(HDF5.read_dataset(f, "t_coorb"))
        masks = multidomain ? Matrix{Float64}(HDF5.read_dataset(f, "masks")) : ones(length(t), 1)
        components = CoorbitalComponent[]
        starts = zeros(Int, ell_max - 1, ell_max + 1)
        for ell in 2:ell_max, m in 0:ell
            starts[ell - 1, m + 1] = length(components) + 1
            names = m == 0 ? ("real", "imag") : ("Re+", "Re-", "Im+", "Im-")
            for name in names, sub in 1:nsub
                group = HDF5.open_group(f, "hCoorb_$(ell)_$(m)_$name" * suffix(sub))
                push!(components, read_datapiece(group, settings))
                close(group)
            end
        end
        return CoorbitalSurrogate(params, settings, t, masks, ell_max, components, starts)
    end
end

"""
    load_coorbital(model::Symbol; download = true) -> CoorbitalSurrogate

Load the coorbital surrogate of `:NRSur7dq4` or `:NRSur7dq4v2` from [`data_path`](@ref).
"""
function load_coorbital(model::Symbol; download::Bool = true)
    return load_coorbital(data_path(model; download), MODEL_FILES[model].params)
end
