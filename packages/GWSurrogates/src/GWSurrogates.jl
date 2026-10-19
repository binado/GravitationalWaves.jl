"""
    GWSurrogates

Differentiable Julia ports of the NRSur7dq4 family of numerical-relativity surrogates from
`gwsurrogate` (MIT). The model data is CC-BY-4.0 (Zenodo 10.5281/zenodo.22257361).

The kernels (natural cubic splines, closed-form Wigner-D matrices, guarded norms) are written
for ForwardDiff, Enzyme and Mooncake: no `try`/`catch`, no `Dict`, and branches only through
`ifelse` on values.
"""
module GWSurrogates

public NaturalCubicSpline, spline_eval, wigner_D, wigner_D_element, guarded_sqrt, rescale_to_norm

include("spline.jl")
include("wignerd.jl")
include("norms.jl")

end
