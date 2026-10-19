"""Apply the prototype's physical sampling, continuous peak and conditioning to oracle modes."""
from pathlib import Path
import lal
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from generate import ROOT, CASES, verify_versions


def smoothstep(x):
    x = np.clip(x, 0, 1)
    return x**3 * (10 + x * (-15 + 6*x))


def condition(path, mass, df=0.0625, fmax=2048):
    data = np.load(path)
    time, mode = data["t"], data["h22"]
    amp = CubicSpline(time, np.abs(mode))
    i = np.argmax(np.abs(mode))
    # The upstream (2,2) matching target has zero amplitude slope at the
    # continuous attachment time. Recover its coordinate after the sampled
    # upstream output-time shift, instead of defining timing by that sample.
    peak = float(data["t_attach"]) + float(time[0])
    # Interpolate amplitude and unwrapped phase to avoid high-frequency complex interpolation error.
    phase = CubicSpline(time, np.unwrap(np.angle(mode)))
    Msec = mass * lal.MTSUN_SI
    N = round(2*fmax/df)
    dt, t_start = 1/(2*fmax), 1 - 1/df
    t = (t_start + np.arange(N)*dt) / Msec + peak
    inside = (t >= time[0]) & (t <= time[-1])
    h = np.zeros(t.shape, dtype=complex)
    h[inside] = amp(t[inside]) * np.exp(1j * phase(t[inside]))
    freq = np.zeros(t.shape)
    freq[inside] = -phase.derivative()(t[inside]) / (2*np.pi*Msec)
    onset = np.ones(t.shape)
    onset[(t < data['t_attach']) & inside] = smoothstep((freq[(t < data['t_attach']) & inside] - 16)/4)
    trel = t_start + np.arange(N)*dt
    ending = smoothstep((1 - trel)/.25)
    h *= onset*ending
    cosi = np.cos(.7)
    yp, ym = np.sqrt(5/(64*np.pi))*(1+cosi)**2, np.sqrt(5/(64*np.pi))*(1-cosi)**2
    rot = np.exp(2j*(np.pi/2-.3))
    scale = -mass*lal.MRSUN_SI/(500*lal.PC_SI*1e6)
    strain = scale * (yp*h*rot + ym*np.conj(h*rot))
    f = np.arange(N//2+1)*df
    correction = np.exp(-2j*np.pi*f*t_start)*dt
    hp, hc = np.fft.rfft(strain.real)*correction, np.fft.rfft(-strain.imag)*correction
    hp[f<20] = 0
    hc[f<20] = 0
    suffix = '_fd.npz' if df == .0625 and fmax == 2048 else '_short_fd.npz'
    dest = path.with_name(path.stem + suffix)
    np.savez(dest, plus=hp, cross=hc, peak=peak, t_attach=data['t_attach'],
             t_isco=data['t_isco'], nqc=data['nqc'])
    print(dest, flush=True)
    return f, hp, hc


if __name__ == '__main__':
    verify_versions()
    for mass, q, _, _ in CASES:
        condition(ROOT / f'reference/seobnr/M{mass}_q{q}.npz', mass)
