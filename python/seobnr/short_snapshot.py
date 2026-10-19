"""Compact spectral anchors for a short, non-Python package regression test."""
import numpy as np
import tomlkit
from condition import condition
from generate import ROOT
f, hp, hc = condition(ROOT/'reference/seobnr/M40_q4.npz', 40, df=.125, fmax=1024)
i = np.unique(np.round(np.geomspace(160, 8192, 80)).astype(int))
snapshot = dict(mass=40, q=4, chi1=-.6, chi2=.6, df=.125, fmax=1024,
                indices=(i+1).tolist(), plus_real=hp[i].real.tolist(), plus_imag=hp[i].imag.tolist(),
                cross_real=hc[i].real.tolist(), cross_imag=hc[i].imag.tolist(),
                plus_norm=float(np.linalg.norm(hp)), cross_norm=float(np.linalg.norm(hc)))
(ROOT/'packages/GWSEOBNR/test/short_waveform.toml').write_text(tomlkit.dumps(snapshot))
