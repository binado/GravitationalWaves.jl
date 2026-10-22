# Waveform numerical acceptance

`parameters.csv` is shared by public-output oracle comparisons and the heavy
optimal-SNR AD suites. It is data only: `cases.jl` validates the fixed header,
unique IDs, finite values, ordered positive detector-frame masses, spin magnitudes,
and positive distances. Masses are solar masses, distances Mpc, angles radians,
and spins dimensionless Cartesian components in the LAL source frame at the
reference epoch. Inclination is measured from the orbital angular momentum.
The reference phase is passed to each backend's documented phase convention;
these suites compare each implementation to its own oracle, not different models.

FD adapters require zero transverse spins and use magnitude/tilt for the LAL
interface. Surrogates use total mass, q=m1/m2 (the inverse of bilby's mass ratio),
Cartesian spins, and observer azimuth pi/2-phi_ref. Generator settings and explicit
applicability lists live in Julia. The common equal/unequal/zero/aligned/antialigned
BBHs enter both suites; the precessing row is surrogate-only. FD extremes retain
separate model-appropriate oracle coverage. No precessing row is projected onto
an aligned model.

## Entry points

- `Pkg.test()` for FrequencyDomainWaveforms compares both unwindowed polarizations
  with LAL, including existing PN options, cutoffs and grid variants.
- `just test-oracle-surrogates` compares both surrogate models with the pinned
  gwsurrogate 1.2.0 PyCall oracle in its own environment, without Mooncake. CI runs
  this as a required job whenever SurrogateWaveforms is selected.
- `just test-integration-waveforms` and `just test-integration-surrogates` run DI
  ForwardDiff/Mooncake optimal-SNR tests locally.
- `just test-integration-series` validates production physical Fourier AD against
  an explicit DFT, including DC, interior bins, odd/even final bins, multiple
  channels, both real and imaginary cotangents, and prepared-gradient reuse.

All oracle comparisons precede windowing. FD peak-relative tolerance is 1e-10;
surrogate tolerance is 1e-11. Identically zero oracle outputs require exact zero.
The sole surrogate exclusion is the known upstream last-interval spline defect:
only times after the penultimate coorbital knot are omitted from numerical
agreement; shapes and finiteness still cover every sample. Correct interpolation
is checked independently against DataInterpolations and the four-knot regression.

## SNR convention

The FD grid is 0:0.5:1024 Hz; generators use 20–1024 Hz. The TD grid has 432 samples
starting at -0.2 seconds with dt=1/2048 s. The reference is the first dynamics node
(`f_ref=0`), and spins are specified in that reference frame. The window has start
and end Planck widths of 0.02 and 0.005 seconds; 80 zero samples are appended after
windowing, for N=512, df=4 Hz. There is no taper-power renormalization. The Fourier
convention is dt*rfft along rows, with no epoch phase factor. An overall epoch
phase would leave this objective unchanged. Native FD outputs are not tapered.

Fixed response coefficients (Fplus,Fcross)=(0.7,0.3) combine both polarizations.
Each model/case is checked with both bundled ET-D (`ET_000A_18_ETD_hlf`) and CE
reference (`cosmic_explorer_strain`) curves. These are one-sided PSDs in strain²/Hz,
read through GWNoise with log-log interpolation. The analysis band is 20–512 Hz
intersected with table/grid coverage. Only active bins are interpolated; no PSD
extrapolation is permitted. GWMatchedFiltering supplies 4df normalization and
excludes DC; the TD metric explicitly receives N and excludes even-length Nyquist.
Native FD grids have no inferred time length.

Configuration, windows, padding buffers, model data, grids, and PSD metrics are
prepared outside differentiation. Only waveform parameters are differentiated.
Mass/distance coordinates are scaled by their fixture values, spins/angles by one.
FD tilt signs are fixed at preparation, so the local signed-spin coordinate is
smooth at zero. Gradients are compared componentwise after division by base SNR,
with rtol 1e-8 (FD) / 1e-6 (TD) and atol 1e-8. Values use rtol 1e-10. Finite-difference
spot checks use smooth aligned FD and generic-precession TD cases; boundary and
zero-spin cases compare backend conventions without asserting differentiability
where the physical parameterization has a kink. Every case tests cache reuse at
a nearby distance and the analytic inverse-distance derivative. A phase-sensitive
overlap checks the phase derivative that an SNR norm can miss.

## PSD provenance follow-up

Reuse the maintainer-supplied GWNoise tables verbatim; do not duplicate them here.
Their exact external source/version and licensing attribution were not supplied.
Resolving that attribution remains an explicit follow-up, not a completed claim
of this test migration. GWNoise documents units and validates frequency coverage.
