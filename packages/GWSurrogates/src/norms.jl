"""
    guarded_sqrt(s)

`sqrt(s)` for `s > 0` and `zero(s)` otherwise, with a finite derivative everywhere. The argument
of the discarded `sqrt` is replaced by one, so reverse-mode backends never evaluate the infinite
slope at zero.
"""
@inline function guarded_sqrt(s)
    positive = s > zero(s)
    return ifelse(positive, sqrt(ifelse(positive, s, one(s))), zero(s))
end

"""
    rescale_to_norm(v, target)

Rescale the tuple or static vector `v` to Euclidean norm `target`. A zero vector is returned
unchanged instead of producing `NaN`; the true gradient has a cone-shaped kink there, so
gradient-based samplers should not treat exactly-zero spin as a regular point.
"""
@inline function rescale_to_norm(v, target)
    s = sum(abs2, v)
    positive = s > zero(s)
    scale = ifelse(positive, target / sqrt(ifelse(positive, s, one(s))), one(s))
    return v .* scale
end
