# Heavy waveform acceptance suites

These suites are too heavy for CI and run locally from the repository root. They use the
generic checks of [AbstractWaveformsTests](../packages/AbstractWaveformsTests): reference-generator
comparisons by mismatch (`test_reference`) and self-calibrating gradient checks
(`test_gradients`, `test_snr_gradients`, `test_phase_gradients`, `test_waveform_jacobian`), with
`FiniteDifferences` as the oracle and every call going through DifferentiationInterface.

The shared physical parameters are `packages/AbstractWaveformsTests/data/parameters.csv`
(solar masses, Mpc, radians, dimensionless Cartesian spins in the LAL source frame at the
reference epoch; inclination is measured from the orbital angular momentum). FD generators use
magnitude/tilt spins (zero transverse spin required); surrogates use total mass, q = m1/m2 (the
inverse of bilby's mass ratio), Cartesian spins and observer azimuth pi/2 - phi_ref. Each suite
lists the case IDs a model covers; no precessing row is projected onto an aligned model.

## Entry points

- `just test-reference-waveforms`: `FrequencyDomainWaveforms` against `LALWaveforms`
  (`FrequencyDomainWaveforms/test_reference.jl`; needs `LALWAVEFORMS_LIBDIR`). Unwindowed
  polarizations must agree to mismatch 1e-14 and pointwise to 1e-10 relative to the peak, with
  identical zero support.
- `just test-reference-surrogates`: both surrogate models against upstream gwsurrogate through
  `SurrogatePyCallWaveforms` (`SurrogateWaveforms/gwsurrogate_reference`; needs the pinned Python
  environment). Modes agree to 1e-11 relative to the peak; the sole exclusion is the known
  upstream last-interval spline defect (times after the penultimate coorbital knot).
- `just test-integration-waveforms` and `just test-integration-surrogates`: optimal-SNR gradient
  suites with ForwardDiff and Mooncake, plus forward-mode waveform Jacobians for the FD models.
- `just test-integration-series` validates production physical Fourier AD against
  an explicit DFT, including DC, interior bins, odd/even final bins, multiple
  channels, both real and imaginary cotangents, and prepared-gradient reuse.

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
smooth at zero. The shared PSD setup is `common/snr.jl`. Gradients are compared componentwise after division by base SNR,
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
