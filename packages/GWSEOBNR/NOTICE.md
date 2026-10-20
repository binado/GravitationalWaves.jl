# Attribution and modifications

GWSEOBNR is GPL-3.0-or-later. Its package-local LICENSE contains the full
license. Other packages in this monorepo retain their existing licenses.

The aligned Hamiltonian, analytic gradients and Hessian, waveform/flux
coefficients, calibrated input-value and merger fits, NQC equations, and
radiation reaction are derived from **pySEOBNR 0.3.7**, by Serguei Ossokine,
Lorenzo Pompili, Deyan Mihaylov, Antoni Ramos Buades, Michael Puerrer, Hector
Estelles, Raffi Enficiaud, and contributors. Upstream:
https://git.ligo.org/waveforms/software/pyseobnr . Source archive and individual
SHA-256 checksums are recorded in `provenance.toml`.

Modifications: translated Python/Cython equations to generic Julia arithmetic;
normalized component masses; precomputed draw-dependent flux coefficients;
replaced GSL evolution by OrdinaryDiffEqVerner; implicit first-order derivatives
of initial-condition roots and continuous timing roots; continuous trajectory
matching and amplitude-peak timing; fixed physical FFT sampling and conditioning.
Initial-condition equations are solved as a coupled three-variable Newton system
with NonlinearSolve rather than sequential custom Newton/bisection solves. Their
physical residual equations are unchanged.
These changes are experimental and require the acceptance checks in integration.
The scalar remnant fits originate in LAL's `lalinference.imrtgr.nrutils`
(LAL 7.24, February 2025), as identified by pySEOBNR's `nr_utils.py`.

`qnm_data.jl` contains piecewise cubic interpolation coefficients exported from
**qnm 0.4.4**, by Leo C. Stein, licensed under MIT; see `qnm-LICENSE`.
Positive remnant spins use (s,l,m,n)=(-2,2,2,0); negative spins use
(-2,2,-2,0) evaluated at the absolute spin, as in the pinned oracle.
The exporter is `python/seobnr/export_qnm.py`; Julia never downloads QNM data.
Please cite Stein, JOSS 4(42), 1683 (2019), doi:10.21105/joss.01683,
and the SEOBNRv5HM model papers listed by pySEOBNR.
