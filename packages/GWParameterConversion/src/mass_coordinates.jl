# NamedTuple coordinate changes have explicit physical domains. The positional
# conversions retain bilby's arithmetic and their existing permissive behavior.
_coordinate_value(x::Real) = x

function _check_coordinate_keys(p::NamedTuple, keys::Tuple)
    length(p) == length(keys) && all(k -> haskey(p, k), keys) ||
        throw(ArgumentError("expected exactly the coordinate keys $keys"))
    all(x -> x isa Real, values(p)) || throw(ArgumentError("coordinates must be Real"))
    return nothing
end

function _check_positive_mass(m::Real)
    v = _coordinate_value(m)
    isfinite(v) && v > 0 || throw(DomainError(m, "mass must be finite and positive"))
    return nothing
end

function _check_mass_ratio(q::Real)
    v = _coordinate_value(q)
    isfinite(v) && 0 < v <= 1 || throw(DomainError(q, "mass_ratio must satisfy 0 < q ≤ 1"))
    return nothing
end

function _check_component_coordinates(p::NamedTuple)
    _check_coordinate_keys(p, (:mass_1, :mass_2))
    _check_positive_mass(p.mass_1)
    _check_positive_mass(p.mass_2)
    _coordinate_value(p.mass_1) >= _coordinate_value(p.mass_2) ||
        throw(DomainError(p, "mass_1 must be ≥ mass_2"))
    return nothing
end

"""
    chirp_mass_and_mass_ratio_to_component_masses(p::NamedTuple)

Change `(; chirp_mass, mass_ratio)` to `(; mass_1, mass_2)`. Accepts exactly these
input keys in either order, with finite positive chirp mass and `0 < mass_ratio ≤ 1`.
The two-positional-argument method instead returns an ordinary tuple.

Loading ChangesOfVariables enables `with_logabsdet_jacobian(f, p)` for this method;
loading InverseFunctions registers [`component_masses_to_chirp_mass_and_mass_ratio`](@ref)
as its inverse on this domain.
"""
function chirp_mass_and_mass_ratio_to_component_masses(p::NamedTuple)
    _check_coordinate_keys(p, (:chirp_mass, :mass_ratio))
    _check_positive_mass(p.chirp_mass)
    _check_mass_ratio(p.mass_ratio)
    mass_1, mass_2 = chirp_mass_and_mass_ratio_to_component_masses(p.chirp_mass, p.mass_ratio)
    return (; mass_1, mass_2)
end

"""
    total_mass_and_mass_ratio_to_component_masses(p::NamedTuple)

Change `(; total_mass, mass_ratio)` to `(; mass_1, mass_2)`. Accepts exactly these
input keys in either order, with finite positive total mass and `0 < mass_ratio ≤ 1`.
The two-positional-argument method instead returns an ordinary tuple.

Loading ChangesOfVariables enables `with_logabsdet_jacobian(f, p)` for this method;
loading InverseFunctions registers [`component_masses_to_total_mass_and_mass_ratio`](@ref)
as its inverse on this domain.
"""
function total_mass_and_mass_ratio_to_component_masses(p::NamedTuple)
    _check_coordinate_keys(p, (:total_mass, :mass_ratio))
    _check_positive_mass(p.total_mass)
    _check_mass_ratio(p.mass_ratio)
    mass_1, mass_2 = total_mass_and_mass_ratio_to_component_masses(p.total_mass, p.mass_ratio)
    return (; mass_1, mass_2)
end

"""
    component_masses_to_chirp_mass_and_mass_ratio(mass_1, mass_2) -> (chirp_mass, mass_ratio)
    component_masses_to_chirp_mass_and_mass_ratio(p::NamedTuple)

Change component masses to chirp mass and mass ratio using the existing scalar
conversions. The NamedTuple method accepts exactly `mass_1` and `mass_2` in either
order, requires finite `mass_1 ≥ mass_2 > 0`, and returns `(; chirp_mass, mass_ratio)`.
The positional method returns an ordinary tuple without domain validation.

The NamedTuple method supports ChangesOfVariables Jacobians and is the registered
InverseFunctions inverse of [`chirp_mass_and_mass_ratio_to_component_masses`](@ref).
"""
function component_masses_to_chirp_mass_and_mass_ratio(mass_1::Real, mass_2::Real)
    return chirp_mass(mass_1, mass_2), mass_ratio(mass_1, mass_2)
end

function component_masses_to_chirp_mass_and_mass_ratio(p::NamedTuple)
    _check_component_coordinates(p)
    mc, q = component_masses_to_chirp_mass_and_mass_ratio(p.mass_1, p.mass_2)
    return (; chirp_mass = mc, mass_ratio = q)
end

"""
    component_masses_to_total_mass_and_mass_ratio(mass_1, mass_2) -> (total_mass, mass_ratio)
    component_masses_to_total_mass_and_mass_ratio(p::NamedTuple)

Change component masses to total mass and mass ratio using the existing scalar
conversions. The NamedTuple method accepts exactly `mass_1` and `mass_2` in either
order, requires finite `mass_1 ≥ mass_2 > 0`, and returns `(; total_mass, mass_ratio)`.
The positional method returns an ordinary tuple without domain validation.

The NamedTuple method supports ChangesOfVariables Jacobians and is the registered
InverseFunctions inverse of [`total_mass_and_mass_ratio_to_component_masses`](@ref).
"""
function component_masses_to_total_mass_and_mass_ratio(mass_1::Real, mass_2::Real)
    return total_mass(mass_1, mass_2), mass_ratio(mass_1, mass_2)
end

function component_masses_to_total_mass_and_mass_ratio(p::NamedTuple)
    _check_component_coordinates(p)
    mt, q = component_masses_to_total_mass_and_mass_ratio(p.mass_1, p.mass_2)
    return (; total_mass = mt, mass_ratio = q)
end
