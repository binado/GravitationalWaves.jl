module GWSEOBNRSciMLSensitivityExt
using SciMLSensitivity: SteadyStateAdjoint
import GWSEOBNR: initial_condition_sensealg

# Out-of-place automatic selection uses Zygote, which cannot differentiate our
# mutable coefficient construction. Use dense ForwardDiff residual Jacobians
# inside the native implicit adjoint, as validated by the nine-case prototype.
initial_condition_sensealg(::Type{<:Real}) = SteadyStateAdjoint(; autojacvec = false)
end
