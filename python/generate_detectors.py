"""Generate independent LAL/Bilby detector oracles in reference/detectors.json.

No Julia results are used. Native tensors and REAL4 LAL tensors are recorded
separately so cached-tensor rounding cannot hide a response-kernel error.
"""
from importlib.metadata import version
import json
from pathlib import Path

import bilby
import lal
import numpy as np
from bilby.gw.detector import get_empty_interferometer
from bilby_cython.geometry import get_polarization_tensor, calculate_arm

OUT = Path(__file__).resolve().parent.parent / "reference" / "detectors.json"
PREFIXES = dict(H1="LHO_4K", L1="LLO_4K", V1="VIRGO", K1="KAGRA")
LEAPS = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
         346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
         599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017]


def pair(h):
    return [[float(z.real), float(z.imag)] for z in h]


def main():
    epochs = [0.0, 630763213.0, 1126259462.123456, 1400000000.25, 2147483647.0]
    epochs += [t + delta for t in LEAPS for delta in (-1.25, -0.25, 0.0, 0.25, 1.0, 1.25)]
    data = dict(versions={p: version(p) for p in ("bilby", "lalsuite", "numpy")},
                lal_version=lal.__version__,
                gmst=[dict(gps=t, gmst=lal.GreenwichMeanSiderealTime(lal.LIGOTimeGPS(t))) for t in epochs],
                detectors=[], geodetic=[], responses=[], projections=[])
    for lat, lon, height, az, opening, tx, ty in (
        (0., 0., 0., 0., np.pi/3, 0., 0.),
        (0.6, -1.1, 50., 0.2, 1.1, 0.03, -0.02),
        (-1.2, 2.7, -30., 4.5, np.pi/2, -0.001, 0.002),
    ):
        # WGS-84 radius evaluated independently; Bilby's vertex helper uses a
        # rounded polar radius, so it is deliberately not the vertex oracle.
        b_over_a = 1 - 1 / 298.257223563
        radius = 6378137 / np.sqrt(np.cos(lat)**2 + b_over_a**2 * np.sin(lat)**2)
        vertex = [(radius+height)*np.cos(lat)*np.cos(lon),
                  (radius+height)*np.cos(lat)*np.sin(lon),
                  (b_over_a**2*radius+height)*np.sin(lat)]
        # Bilby measures azimuth counterclockwise from east; our geodetic
        # convention (like LAL) is clockwise from north.
        u = calculate_arm(tx, np.pi/2-az, lon, lat)
        v = calculate_arm(ty, np.pi/2-(az+opening), lon, lat)
        data['geodetic'].append(dict(latitude=lat, longitude=lon, elevation=height,
            xarm_azimuth=az, yarm_azimuth=az+opening, xarm_tilt=tx, yarm_tilt=ty,
            arm_lengths=[3000.,4000.], vertex=vertex, xarm=list(u), yarm=list(v)))
    rng = np.random.default_rng(749)
    for name, prefix in PREFIXES.items():
        def constants(stem):
            return np.array([getattr(lal, f"{prefix}_{stem}_{axis}") for axis in "XYZ"])
        vertex = np.array([getattr(lal, f"{prefix}_VERTEX_LOCATION_{axis}_SI") for axis in "XYZ"])
        u, v = constants("ARM_X_DIRECTION"), constants("ARM_Y_DIRECTION")
        lengths = [2 * getattr(lal, f"{prefix}_DETECTOR_ARM_{axis}_MIDPOINT_SI") for axis in "XY"]
        native = (np.outer(u, u) - np.outer(v, v)) / 2
        cached = next(d for d in lal.CachedDetectors if d.frDetector.prefix == name)
        data["detectors"].append(dict(name=name, vertex=vertex.tolist(), xarm=u.tolist(), yarm=v.tolist(),
                                      arm_lengths=lengths, native_tensor=native.tolist(), cached_tensor=cached.response.tolist()))
        for _ in range(16):
            ra, dec, psi, gmst = rng.uniform(0, 2*np.pi), rng.uniform(-np.pi/2, np.pi/2), rng.uniform(-4, 4), rng.uniform(-10, 10)
            rounded = native.astype(np.float32)
            plus, cross = lal.ComputeDetAMResponse(rounded, ra, dec, psi, gmst)
            # Contract double geometry with a separately implemented wave basis.
            # Bilby's public basis accepts GPS; fixed-GMST basis below follows
            # its theta/phi definition without involving GPS conversion.
            theta, phi = np.pi/2-dec, ra-gmst
            m = np.array([-np.sin(phi), np.cos(phi), 0.])
            n = np.array([-np.cos(theta)*np.cos(phi), -np.cos(theta)*np.sin(phi), np.sin(theta)])
            X, Y = -np.cos(psi)*m + np.sin(psi)*n, np.sin(psi)*m + np.cos(psi)*n
            native_response = [np.sum(native*(np.outer(X,X)-np.outer(Y,Y))), np.sum(native*(np.outer(X,Y)+np.outer(Y,X)))]
            data["responses"].append(dict(name=name, ra=ra, dec=dec, psi=psi, gmst=gmst,
                rounded_tensor=rounded.tolist(), lal=[plus,cross], native=native_response,
                cached=list(lal.ComputeDetAMResponse(cached.response,ra,dec,psi,gmst))))
        ifo = get_empty_interferometer(name)
        # Override all geometry caches: compare identical Cartesian geometry.
        ifo.geometry._x, ifo.geometry._y = u, v
        ifo.geometry._x_updated = ifo.geometry._y_updated = True
        ifo.geometry._vertex = vertex
        ifo.geometry._detector_tensor = native
        ifo.geometry._vertex_updated = True
        ifo.geometry._detector_tensor_updated = True
        np.testing.assert_array_equal(ifo.geometry.detector_tensor, native)
        np.testing.assert_array_equal(ifo.geometry.vertex, vertex)
        for grid in (list(12.5 + np.arange(5)*23.), list(np.arange(9)*0.25)):
            f = np.array(grid)
            hp = rng.normal(size=len(f)) + 1j*rng.normal(size=len(f))
            hc = rng.normal(size=len(f)) + 1j*rng.normal(size=len(f))
            sky = dict(ra=1.23, dec=-0.42, psi=0.71, geocent_time=1126259462.123456)
            for reference_time in (sky['geocent_time'], sky['geocent_time'] + 1234.5):
                start = sky['geocent_time']-2.25
                ifo.strain_data.start_time = start
                ifo.reference_time = reference_time
                h = ifo.get_detector_response(dict(plus=hp, cross=hc), sky, frequencies=f)
                delay = lal.TimeDelayFromEarthCenter(vertex, sky['ra'], sky['dec'], lal.LIGOTimeGPS(sky['geocent_time']))
                data['projections'].append(dict(name=name, f=grid, plus=pair(hp), cross=pair(hc), sky=sky,
                    reference_time=reference_time, start_time=start, strain=pair(h), delay=delay,
                    antenna=[float(np.sum(native*get_polarization_tensor(sky['ra'],sky['dec'],reference_time,sky['psi'],mode))) for mode in ('plus','cross')]))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(data, indent=1) + '\n')
    print(f"wrote {OUT}: {len(data['responses'])} responses, {len(data['projections'])} projections")


if __name__ == '__main__':
    main()
