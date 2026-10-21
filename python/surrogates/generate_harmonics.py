"""Regenerate the small, committed harmonic oracle with gwtools 1.2.0."""
from pathlib import Path
import math
from gwtools.harmonics import sYlm

out = Path(__file__).resolve().parents[2] / "packages/AbstractWaveforms/test/harmonics.csv"
with out.open("w") as f:
    f.write("# gwtools 1.2.0 sYlm(-2,l,m,theta,phi): l,m,theta,phi,real,imag\n")
    for theta in (0., 1e-4, .37, 1.4, math.pi - 1e-4, math.pi):
        for l in range(2, 6):
            for m in range(-l, l + 1):
                y = sYlm(-2, l, m, theta, .63)
                f.write(f"{l},{m},{theta:.17g},0.63,{y.real:.17g},{y.imag:.17g}\n")
