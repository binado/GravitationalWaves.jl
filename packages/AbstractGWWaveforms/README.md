# AbstractGWWaveforms

The interface shared by the gravitational-wave waveform backends in this repository,
[`GWWaveforms`](../GWWaveforms) (pure Julia, differentiable) and [`LALWaveforms`](../LALWaveforms)
(LALSimulation through `ccall`). It has no dependencies.

It defines `AbstractWaveformGenerator` and the generic function

```julia
hphc_fd(f, params, generator) -> (; plus, cross)
```

where `f` is a uniform frequency grid starting at 0, `params` a `NamedTuple` with bilby's
parameter names, and `generator` a concrete generator such as `TaylorF2()`, `IMRPhenomD()` or
`LALWaveformGenerator("IMRPhenomXPHM")`. Code written against this package works with any backend.

and `validate_parameters(f, params, generator)`, which throws `ArgumentError` on invalid input and
returns `nothing` otherwise. Backends that compilers can trace (`GWWaveforms` with Reactant) do not
check parameter values in `hphc_fd`; callers run `validate_parameters` once on the host instead.
The fallback accepts everything, so generic code can always call it.

Every generator provides `approximant`, `reference_frequency`, `minimum_frequency` and
`maximum_frequency`, and may answer the capability queries `spin_support` (`:none`, `:aligned` or
`:precessing`), `supports_mode_array` and `is_differentiable`. `check_frequency_grid` and
`frequency_bounds` hold the grid checks and masking bounds backends share.

`test/contract.jl` defines `test_generator_contract(g; params, f)`, which backends `include` from
their tests to check a generator against the contract.

MIT; see [LICENSE](LICENSE).

## Aligned-spin quadrupole interface

Generators for aligned-spin, quadrupole-only radiation may subtype
`AlignedSpinQuadrupolarWaveformGenerator <: AbstractWaveformGenerator` and implement

```julia
amplitude_phase(f, params, generator) -> (; amplitude, phase)
```

The returned real vectors describe the inclination-independent strain
`amplitude .* exp.(im .* phase)`. Amplitude is nonnegative and includes distance normalization;
phase is analytic and unwrapped, including reference-phase alignment and model time shifts.
Both vectors are zero at DC and outside configured or intrinsic model support. Phase is zero
wherever amplitude vanishes. Frequencies are in Hz and phase is in radians.

The default `hphc_fd` evaluates `amplitude_phase` once and applies the inclination factors
`(1 + cos(theta_jn)^2)/2` for plus and `-im * cos(theta_jn)` for cross.
Generators continue to store fixed configuration.


## Higher modes in either domain

Subtype `HigherModesWaveformGenerator` and implement either or both:

```julia
waveform_modes_td(generator, params; kwargs...) -> (; t, mode_indices, modes)
waveform_modes_fd(f, params, generator) -> (; mode_indices, positive, negative)
```

The default `hphc_td` and `hphc_fd` evaluate the corresponding method once. No automatic
FFT conversion is provided. Modes are standard spin-weight -2 spherical modes in a fixed
inertial/source frame, normalized by luminosity distance and excluding observer inclination.
Each backend documents its frame and reference-phase convention. TD times are seconds and
modes are physical strain. `mode_indices` is a nonempty ordered tuple of unique signed `(l,m)`
labels, with `2 ≤ l ≤ 5`, including `m=0`; matrices have shape `(sample, mode)`.

FD `positive[i,j]` is h̃ⱼ(f[i]) and `negative[i,j]` is h̃ⱼ(-f[i]), in the same row order.
DC values agree, and both branches vanish outside configured or intrinsic support. The existing
nonnegative, uniform grid contract applies. The assembler independently sums H(f) and H(-f),
then returns `(H(f) + conj(H(-f)))/2` and `im*(H(f) - conj(H(-f)))/2`. It does not assume
any symmetry between signed modes. TD assembly returns `real(H)` and `-imag(H)`.

Backend helpers `mode_angles`, `mode_sum`, and `sylm_m2` are public.
`mode_angles(params, generator)` defaults to `(params.theta_jn, 0)` with reference phase already
in the modes. Override it for another viewing convention. The half-angle harmonic polynomial
supports finite derivatives at the poles and preserves promoted numeric types. Angular
coefficients are computed once per call, independently of sample count.

`test/contract.jl` also provides `test_modes_td_contract` and `test_modes_fd_contract`;
`viewing_params` supplies a second observer configuration to check mode independence.
The small committed `test/harmonics.csv` oracle is regenerated with
`uv run --project python/surrogates python python/surrogates/generate_harmonics.py`.
