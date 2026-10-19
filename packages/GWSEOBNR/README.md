# GWSEOBNR (experimental)

Native Julia aligned-spin SEOBNRv5HM prototype, tied to pySEOBNR **0.3.7**.
**Only (2,±2) waveform output is implemented.** The radiation-reaction flux
retains all upstream multipoles through ℓ=8.

This is an unvalidated research prototype. `is_differentiable` returns **false**
until every primal, derivative, convergence, and composition acceptance gate
passes. ForwardDiff CPU first derivatives are the intended contract; no Hessian,
reverse-mode, Reactant, GPU, precession, eccentricity, or performance guarantee.

Inputs use the shared bilby parameter names and GWParameterConversion mass
conversions. Detector-frame total masses are restricted to 40–100 solar masses,
ordered mass ratios mass_1/mass_2 to 1–4, spin magnitudes to 0–0.6 with aligned
or anti-aligned tilts, positive distances, and zero eccentricity. Fixed nonzero
spin magnitudes can be differentiated with the tilt held at 0 or π.

```
using GWSEOBNR
f = 0:0.0625:2048
p = (mass_1=40., mass_2=20., luminosity_distance=500.,
     theta_jn=0.7, phase=0.3, a_1=0.4, a_2=0.3, tilt_2=π)
g = SEOBNRv5HM()
validate_parameters(f, p, g)
h = hphc_fd(f, p, g)
```

Defaults: minimum/reference frequency 20 Hz (initially required equal),
maximum frequency `nothing` (last grid bin), relative/absolute solve tolerances
1e-11/1e-12, and post-merger duration 1 second. Keep configuration, sampling,
and tolerances fixed during differentiation.

The waveform starts evolving at 0.8 f_min. Orbital phase is referenced to f_ref;
time zero is the continuous amplitude maximum. For K frequency samples,
N=2(K−1), T=1/Δf, and Δt=T/N. Samples cover [−T+t_post,t_post).
A quintic onset spans instantaneous frequency 0.8 f_min to f_min; a final
0.25-second taper starts only after at least 30 QNM damping times. Inadequate
segments, failed roots/solves, missing ISCOs, and attachment clamping requirements
raise explicit errors. The physical spectrum is Δt*rfft with the start-time
phase correction. Shared FFT APIs are unchanged.

The Hamiltonian uses the upstream analytic gradient and Hessian; root partials
use implicit differentiation. Evolution uses Vern9 dense output. The complex
tail and FFT have package-local ForwardDiff adapters. Timing and matching use
continuous interpolants, with fixed five-point time derivatives for matching.
Physical piecewise formulas are retained.

Loading EnzymeCore enables an optional reverse rule for the local real FFT on
dense Float32/Float64 vectors (one reverse lane). Planning uses scratch storage;
the pullback applies the constant plan's adjoint. This component rule does not
establish reverse-mode waveform support: root and ODE compilation blockers remain.
Run `just test-integration-seobnr-fft` for direct-DFT and gradient checks.

Small committed component/QNM snapshots run without Python. The separate locked
oracle environment is documented in [python/seobnr](../../python/seobnr/README.md).
Full waveform/gradient acceptance belongs to `integration/GWSEOBNR`; large
fixtures are gitignored. See `ACCEPTANCE.md` for measured status and limitations.

GPL-3.0-or-later; see [NOTICE.md](NOTICE.md), [LICENSE](LICENSE), and
[provenance.toml](provenance.toml). Existing workspace package licenses are unchanged.
