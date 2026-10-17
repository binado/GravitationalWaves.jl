module GWParameterConversionChangesOfVariablesExt

using GWParameterConversion: chirp_mass_and_mass_ratio_to_component_masses,
    total_mass_and_mass_ratio_to_component_masses,
    component_masses_to_chirp_mass_and_mass_ratio,
    component_masses_to_total_mass_and_mass_ratio
import ChangesOfVariables: with_logabsdet_jacobian

# |∂(m₁,m₂)/∂(Mc,q)| = Mc (1+q)^(2/5) q^(-6/5).
_chirp_logjac(mc, q) = log(mc) + (2 * log1p(q) - 6 * log(q)) / 5
# |∂(m₁,m₂)/∂(M,q)| = M / (1+q)^2.
_total_logjac(mt, q) = log(mt) - 2 * log1p(q)

function with_logabsdet_jacobian(
        f::typeof(chirp_mass_and_mass_ratio_to_component_masses), p::NamedTuple
    )
    y = f(p)
    return y, _chirp_logjac(p.chirp_mass, p.mass_ratio)
end

function with_logabsdet_jacobian(
        f::typeof(total_mass_and_mass_ratio_to_component_masses), p::NamedTuple
    )
    y = f(p)
    return y, _total_logjac(p.total_mass, p.mass_ratio)
end

function with_logabsdet_jacobian(
        f::typeof(component_masses_to_chirp_mass_and_mass_ratio), p::NamedTuple
    )
    y = f(p)
    return y, -_chirp_logjac(y.chirp_mass, y.mass_ratio)
end

function with_logabsdet_jacobian(
        f::typeof(component_masses_to_total_mass_and_mass_ratio), p::NamedTuple
    )
    y = f(p)
    return y, -_total_logjac(y.total_mass, y.mass_ratio)
end

end
