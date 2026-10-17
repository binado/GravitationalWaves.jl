# XLALLeapSeconds.h: GPS seconds at 23:59:60 UTC (1981-07-01 through 2017-01-01).
const LEAP_SECONDS = (
    46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
    346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
    599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017,
)

"""
    greenwich_mean_sidereal_time(gps_time)

LAL-compatible, unwrapped Greenwich mean sidereal angle in radians relative to
J2000. GPS input is in seconds, including fractions. Supports GPS epochs from 0
through 2147483647 inclusive. The pinned leap table ends at 2017-01-01; subsequent
epochs assume GPS−UTC remains 18 seconds. No UT1 or equation-of-equinoxes correction.

Whole and fractional seconds are kept separate. Derivatives are meaningful away
from integer-second rounding boundaries and leap transitions. Input is assumed
finite and real, in the supported epoch range; argument values are not validated.
"""
function greenwich_mean_sidereal_time(gps_time)
    # LAL uses double precision even when the supplied epoch has a narrower
    # storage type. Arithmetic promotion preserves dual and traced values.
    gps = gps_time + 0.0
    whole = floor(gps)
    fraction = gps - whole
    # At 23:59:60 LAL represents the civil second as 60, so it maps to the
    # same Julian day as the next second. Strict > reproduces that convention.
    offset = zero(whole)
    for leap in LEAP_SECONDS
        offset += ifelse(whole > leap, one(whole), zero(whole))
    end
    utc = whole - offset
    day = floor(utc / 86400)
    seconds = utc - day * 86400
    # Match XLALConvertCivilTimeToJD's integer-day plus fractional-day ordering.
    jd = (2444245 + day) + (seconds / 86400 - 0.5)
    t_hi = (jd - 2451545) / 36525
    t_lo = fraction / (36525 * 86400.0)
    t = t_hi + t_lo
    sidereal = (-6.2e-6 * t + 0.093104) * t * t + 67310.54841
    sidereal += 8640184.812866 * t_lo
    sidereal += 3155760000.0 * t_lo
    sidereal += 8640184.812866 * t_hi
    sidereal += 3155760000.0 * t_hi
    return sidereal * Float64(π) / 43200
end
