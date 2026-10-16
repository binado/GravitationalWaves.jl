module GWParameterConversionForwardDiffExt

using ForwardDiff: Dual, value
import GWParameterConversion: _coordinate_value

# Strip nested duals only for domain checks. Arithmetic still uses the original
# inputs, including derivatives, at the equal-mass boundary.
_coordinate_value(x::Dual) = _coordinate_value(value(x))

end
