# Pinned LALDetectors.h constants (LAL 7.7.1 / lalsuite 7.26.16).
# Keep published precision: do not renormalize the arm directions.
const _lal_detectors = (
    H1 = Detector((-2.16141492636e+6, -3.83469517889e+6, 4.60035022664e+6), (-0.22389266154, 0.79983062746, 0.55690487831), (-0.91397818574, 0.02609403989, -0.40492342125), (3995.084, 3995.044)),
    L1 = Detector((-7.42760447238e+4, -5.49628371971e+6, 3.22425701744e+6), (-0.95457412153, -0.1415807734, -0.26218911324), (0.29774156894, -0.48791033647, -0.82054461286), (3995.15, 3995.15)),
    V1 = Detector((4.546374099e+6, 8.42989697626e+5, 4.37857696241e+6), (-0.70045821479, 0.20848948619, 0.68256166277), (-0.05379255368, -0.96908180549, 0.24080451708), (3000.0, 3000.0)),
    K1 = Detector((-3777336.024, 3484898.411, 3765313.697), (-0.375904, -0.8361583, 0.3994189), (0.7164378, 0.01114076, 0.697562), (3026.507, 3023.222)),
)

# Additional sites use the supplied geodetic catalog. Existing LAL sites retain
# their published Cartesian precision. Lengths are nominal design metadata.
include_dependency(joinpath(@__DIR__, "..", "data", "geometry.toml"))

function _catalog_detector(name, length_m)
    g = TOML.parsefile(joinpath(@__DIR__, "..", "data", "geometry.toml"))["detectors"][name]["geometry"]
    return Detector(;
        latitude = g["latitude_rad"], longitude = g["longitude_rad"],
        elevation = g["elevation_m"], xarm_azimuth = g["xarm_azimuth_rad"],
        yarm_azimuth = g["yarm_azimuth_rad"],
        xarm_tilt = get(g, "xarm_tilt_rad", 0.0),
        yarm_tilt = get(g, "yarm_tilt_rad", 0.0),
        arm_lengths = (length_m, length_m),
    )
end

"""
    detectors

Named tuple of predefined `Detector` instances: `H1`, `L1`, `V1`, `K1`,
`E1`–`E3`, `S1`, `S2`, `R1`, `R2`, `C1`, and `C2`.
Access a detector with `detectors.H1` or `detectors[:H1]`. The original four sites use pinned LAL Cartesian geometry; additional sites
use the bundled geodetic catalog and nominal design arm lengths.
"""
const detectors = merge(
    _lal_detectors, (
        E1 = _catalog_detector("E1", 10000.0),
        E2 = _catalog_detector("E2", 10000.0),
        E3 = _catalog_detector("E3", 10000.0),
        S1 = _catalog_detector("S1", 15000.0),
        S2 = _catalog_detector("S2", 15000.0),
        R1 = _catalog_detector("R1", 15000.0),
        R2 = _catalog_detector("R2", 15000.0),
        C1 = _catalog_detector("C1", 40000.0),
        C2 = _catalog_detector("C2", 40000.0),
    )
)
