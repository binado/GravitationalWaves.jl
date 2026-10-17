# Mass conversions ported from bilby/gw/conversion.py. Operation order and exponent
# literals mirror the Python source so results agree bit for bit where libm does.
# Python's `**` always calls `pow`, whereas Julia lowers `x^3` to `x*x*x`, so powers above
# 2 use float literals (`3.0`, `5.0`) where bilby writes `** 3`, `** 5`.

"""
    chirp_mass(mass_1, mass_2)

Chirp mass `(m₁m₂)^(3/5) / (m₁+m₂)^(1/5)`.
"""
chirp_mass(mass_1::Number, mass_2::Number) = (mass_1 * mass_2)^0.6 / (mass_1 + mass_2)^0.2

"""
    total_mass(mass_1, mass_2)
"""
total_mass(mass_1::Number, mass_2::Number) = mass_1 + mass_2

"""
    mass_ratio(mass_1, mass_2)

`mass_2 / mass_1`. Not reordered: it is ≤ 1 only when `mass_1 ≥ mass_2`.
"""
mass_ratio(mass_1::Number, mass_2::Number) = mass_2 / mass_1

"""
    symmetric_mass_ratio(mass_1, mass_2)

`m₁m₂ / (m₁+m₂)²`, clamped to its maximum 1/4 as bilby does to absorb rounding.
"""
symmetric_mass_ratio(mass_1::Number, mass_2::Number) = min((mass_1 * mass_2) / (mass_1 + mass_2)^2, 1 / 4)

"""
    reduced_mass(mass_1, mass_2)
"""
reduced_mass(mass_1::Number, mass_2::Number) = mass_1 * mass_2 / (mass_1 + mass_2)

"""
    total_mass_and_mass_ratio_to_component_masses(total_mass, mass_ratio) -> (mass_1, mass_2)

Note the argument order is the reverse of bilby's signature, which is keyword-called.
"""
function total_mass_and_mass_ratio_to_component_masses(total_mass::Number, mass_ratio::Number)
    mass_1 = total_mass / (1 + mass_ratio)
    mass_2 = mass_1 * mass_ratio
    return mass_1, mass_2
end

"""
    chirp_mass_and_mass_ratio_to_total_mass(chirp_mass, mass_ratio)
"""
function chirp_mass_and_mass_ratio_to_total_mass(chirp_mass::Number, mass_ratio::Number)
    return chirp_mass * (1 + mass_ratio)^1.2 / mass_ratio^0.6
end

"""
    chirp_mass_and_mass_ratio_to_component_masses(chirp_mass, mass_ratio) -> (mass_1, mass_2)
"""
function chirp_mass_and_mass_ratio_to_component_masses(chirp_mass::Number, mass_ratio::Number)
    total = chirp_mass_and_mass_ratio_to_total_mass(chirp_mass, mass_ratio)
    return total_mass_and_mass_ratio_to_component_masses(total, mass_ratio)
end

"""
    chirp_mass_and_total_mass_to_symmetric_mass_ratio(chirp_mass, total_mass)
"""
function chirp_mass_and_total_mass_to_symmetric_mass_ratio(chirp_mass::Number, total_mass::Number)
    return (chirp_mass / total_mass)^(5 / 3)
end

"""
    symmetric_mass_ratio_to_mass_ratio(symmetric_mass_ratio)

Root with `mass_ratio ≤ 1`. NaN for `symmetric_mass_ratio > 1/4`, as in bilby.
"""
function symmetric_mass_ratio_to_mass_ratio(symmetric_mass_ratio::Number)
    temp = 1 / symmetric_mass_ratio / 2 - 1
    return temp - (temp^2 - 1)^0.5
end

"""
    mass_1_and_chirp_mass_to_mass_ratio(mass_1, chirp_mass)

Solve `chirp_mass = mass_1 q^(3/5) / (1+q)^(1/5)` for `q` (closed-form cubic root).
"""
function mass_1_and_chirp_mass_to_mass_ratio(mass_1::Number, chirp_mass::Number)
    temp = (chirp_mass / mass_1)^5.0
    root = 3^0.5 * (27 * temp^2 - 4 * temp^3.0)^0.5 + 9 * temp
    return (2 / 3 / root)^(1 / 3) * temp + (root / (2 * 3^2))^(1 / 3)
end

"""
    mass_2_and_chirp_mass_to_mass_ratio(mass_2, chirp_mass)

Same solve with `mass_2`, which yields `1/q` because the chirp mass is symmetric.
"""
function mass_2_and_chirp_mass_to_mass_ratio(mass_2::Number, chirp_mass::Number)
    return 1 / mass_1_and_chirp_mass_to_mass_ratio(mass_2, chirp_mass)
end

"""
    source_to_detector_mass(mass_source, redshift)

`mass_source * (1 + redshift)`.
"""
source_to_detector_mass(mass_source::Number, redshift::Number) = mass_source * (1 + redshift)

"""
    detector_to_source_mass(mass_detector, redshift)

`mass_detector / (1 + redshift)`.
"""
detector_to_source_mass(mass_detector::Number, redshift::Number) = mass_detector / (1 + redshift)
