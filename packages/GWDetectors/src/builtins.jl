# Pinned LALDetectors.h constants (LAL 7.7.1 / lalsuite 7.26.16).
# Keep published precision: do not renormalize the arm directions.
# Additional geodetic sites are transcribed from data/geometry.toml;
# arm lengths are nominal design metadata.
"""
    detectors

Named tuple of predefined `Detector` instances: `H1`, `L1`, `V1`, `K1`,
`E1`–`E3`, `S1`, `S2`, `R1`, `R2`, `C1`, and `C2`.
Access a detector with `detectors.H1` or `detectors[:H1]`. The original four sites use pinned LAL Cartesian geometry; additional sites
use the bundled geodetic catalog and nominal design arm lengths.
"""
const detectors = (
    H1 = Detector((-2.16141492636e+6, -3.83469517889e+6, 4.60035022664e+6), (-0.22389266154, 0.79983062746, 0.55690487831), (-0.91397818574, 0.02609403989, -0.40492342125), (3995.084, 3995.044)),
    L1 = Detector((-7.42760447238e+4, -5.49628371971e+6, 3.22425701744e+6), (-0.95457412153, -0.1415807734, -0.26218911324), (0.29774156894, -0.48791033647, -0.82054461286), (3995.15, 3995.15)),
    V1 = Detector((4.546374099e+6, 8.42989697626e+5, 4.37857696241e+6), (-0.70045821479, 0.20848948619, 0.68256166277), (-0.05379255368, -0.96908180549, 0.24080451708), (3000.0, 3000.0)),
    K1 = Detector((-3777336.024, 3484898.411, 3765313.697), (-0.375904, -0.8361583, 0.3994189), (0.7164378, 0.01114076, 0.697562), (3026.507, 3023.222)),
    E1 = Detector(;
        latitude = 0.761487152645126,
        longitude = 0.1832595714594046,
        elevation = 52,
        xarm_azimuth = 0.3391174736624982,
        yarm_azimuth = 5.575105229645487,
        arm_lengths = (10000.0, 10000.0),
    ),
    E2 = Detector(;
        latitude = 0.7630579489719208,
        longitude = 0.18413223608540177,
        elevation = 60,
        xarm_azimuth = 4.527907678448889,
        yarm_azimuth = 3.4807101272522916,
        arm_lengths = (10000.0, 10000.0),
    ),
    E3 = Detector(;
        latitude = 0.7627088831215221,
        longitude = 0.18186330805780915,
        elevation = 60,
        xarm_azimuth = 2.4335125760556933,
        yarm_azimuth = 1.3863150248590959,
        arm_lengths = (10000.0, 10000.0),
    ),
    S1 = Detector(;
        latitude = 0.7071492352663692,
        longitude = 0.1643518378961327,
        elevation = 300,
        xarm_azimuth = 1.5707963267948966,
        yarm_azimuth = 0.0,
        arm_lengths = (15000.0, 15000.0),
    ),
    S2 = Detector(;
        latitude = 0.7071492352663692,
        longitude = 0.1643518378961327,
        elevation = 300,
        xarm_azimuth = 1.5707963267948966,
        yarm_azimuth = 0.0,
        arm_lengths = (15000.0, 15000.0),
    ),
    R1 = Detector(;
        latitude = 0.8852843261164459,
        longitude = 0.1033331879916865,
        elevation = 300,
        xarm_azimuth = 1.5707963267948966,
        yarm_azimuth = 0.0,
        arm_lengths = (15000.0, 15000.0),
    ),
    R2 = Detector(;
        latitude = 0.8852843261164459,
        longitude = 0.1033331879916865,
        elevation = 0,
        xarm_azimuth = 0.7853981633974483,
        yarm_azimuth = 5.497787143782138,
        arm_lengths = (15000.0, 15000.0),
    ),
    C1 = Detector(;
        latitude = 0.8107952638302022,
        longitude = -2.084056769166123,
        elevation = 142.554,
        xarm_azimuth = 5.65487724843714,
        yarm_azimuth = 4.084080921642243,
        xarm_tilt = -1.0812314716104872e-5,
        yarm_tilt = 2.181661564992912e-7,
        arm_lengths = (40000.0, 40000.0),
    ),
    C2 = Detector(;
        latitude = 0.533423135055784,
        longitude = -1.584309370782929,
        elevation = -6.574,
        xarm_azimuth = 4.403177723455114,
        yarm_azimuth = 2.8323813966602174,
        xarm_tilt = -5.447172595474303e-6,
        yarm_tilt = -1.065872574192937e-5,
        arm_lengths = (40000.0, 40000.0),
    ),
)
