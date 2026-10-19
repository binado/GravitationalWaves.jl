"""
    GWSurrogates

Differentiable Julia ports of the NRSur7dq4 family of numerical-relativity surrogates from
`gwsurrogate` (MIT). The model data is CC-BY-4.0 (Zenodo 10.5281/zenodo.22257361).

The coorbital waveform surrogate loads its model data (HDF5) lazily from a cache, downloading it
from Zenodo on first use unless `GWSURROGATES_DATA_DIR` points at a local copy.

The kernels (natural cubic splines, closed-form Wigner-D matrices, guarded norms) are written
for ForwardDiff, Enzyme and Mooncake: no `try`/`catch`, no `Dict`, and branches only through
`ifelse` on values.
"""
module GWSurrogates

using StaticArrays: SVector

public NaturalCubicSpline, spline_eval, spline_second_derivatives, wigner_D, wigner_D_element, guarded_sqrt, rescale_to_norm
public FitSettings, FitSet, fit_powers!, eval_fit, fit_params, NRSur7dq4Params, NRSur7dq4v2Params
public CoorbitalSurrogate, coorbital_modes, load_coorbital, data_path
public DynamicsSurrogate, dynamics, output_times, t_ref_from_omega, load_dynamics
public PrecessingSurrogate, load_surrogate, inertial_modes, mode_sum, PrecessingGenerator, hphc_td

include("spline.jl")
include("wignerd.jl")
include("norms.jl")
include("fits.jl")
include("coorbital.jl")
include("dynamics.jl")
include("precessing.jl")
include("data.jl")

end
