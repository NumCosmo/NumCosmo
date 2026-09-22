# APES centre shrinkage: scripts and how to rerun

Analysis behind the `apes-center-shrink` branch. `review.md` has the full write-up with
every table; this file is the map and the restart instructions.

The MCMC catalogs produced while developing this are not here (a few hundred MB, all
reproducible). Every script writes its own catalog and is safe to rerun.

## What the branch changes

1. Centre shrinkage in `NcmStatsDist` (`center-shrink` property, default off), exposed on
   `NcmFitESMCMCWalkerAPES`, in `create_esmcmc`, `create_stats_dist` and
   `numcosmo catalog calibrate --center-shrink`.
2. Kernel variance factor `ncm_stats_dist_kernel_get_var_factor()`, needed because the
   shrinkage factor is `a = 1 / sqrt(1 + kappa h^2 s^2)`. Kernels with no finite
   covariance (Student-t with nu <= 2, so Cauchy) are refused.
3. `ncm_stats_dist_get_href()` now returns the bandwidth actually applied, `a * h`.
4. `numcosmo run mcmc apes --center-shrink`, and `numcosmo generate sampler-test`, which
   writes the synthetic benchmark targets as experiment files.

Bug fixes that are independent of the feature and worth keeping regardless:

* `set_method` / `set_k_type` silently dropped `local_frac` and the covariance type.
* `METHOD_KDE` built a VKDE.
* `NcmNNLS` aborted on an empty passive set and on rank-deficient QR. The rank-deficient
  path had a leftover `g_assert (ret <= 0)` in front of its own DGELSD fallback, so the
  fallback could never run; `numcosmo run mcmc apes --center-shrink` with a small ensemble
  dumped core. Now only a negative LAPACK info, which is a programming error, asserts.
* `prepare_interp` asserted when the fit returned no positive weight.
* The Python `APES.run_mcmc` ignored its `initial_sample`.
* `nc_data_snia_cov_set_abs_mag_set` double counted, which aborted threaded SNIa runs.

## Scripts

| script | what it does |
|--------|--------------|
| `apes_validate.py` | correctness: rebuilds the mixture density in numpy from the public API and compares with `eval`, KS-tests projected draws against the exact 1d mixture law, checks the covariance identity |
| `apes_interp_bench.py` | one-step proposal quality: expected independence-MH acceptance from an i.i.d. training sample |
| `apes_e2e.py` | end-to-end chain on a correlated Gaussian, stationary start |
| `apes_hard.py` | end-to-end on funnel, bimodal, banana, Student-t, Rosenbrock, with target-specific diagnostics |
| `apes_cosmo.py` | end-to-end on real cheap cosmology: JLA SNIa + BAO + H(z), XCDM or spline w(z) |
| `apes_paper50.py` | the paper's truncated Gaussian at d = 50, starting from the prior |
| `apes_gap.py` | scores the walker's own prepared proposal, to tell a bad proposal from a bad setup |
| `apes_ensemble_quality.py` | proposal built from the chain ensemble vs from an i.i.d. sample |
| `apes_burnin.py` | ensemble expansion from a collapsed start (never run) |
| `planck_cost.py` | times one Planck 2018 TT likelihood call |
| `planck_phaseA.py` | Planck baseline chain, adapts an ensemble and writes a catalog |
| `planck_phaseB.py` | seeds from that catalog and measures one candidate configuration |
| `summarize.py` | runs `numcosmo catalog analyze` twice per catalog, taking the burn-in cut from the Constant Break statistic, and tabulates tau / R-1 / ESS / HW |
| `lambda_probe3.py` | NNLS-weight-weighted covariance against the true one; bias by kernel and shrink setting |
| `lambda_sweep.py` | `r = tr(C_lambda)/tr(C)` along the bandwidth |
| `lambda_accept4.py` | out-of-sample one-step acceptance along h, equal weights or NNLS, shrink on or off |
| `lambda_banana.py` | the same on the Haario banana, both weightings side by side |
| `param_ceiling.py` | acceptance of a single fitted Gaussian, the h -> infinity limit |
| `lambda_bimodal.py` | bimodal target: equal weights vs NNLS, plus the single-Gaussian ceiling |
| `run_scalar.sh` | run a probe against a worktree build of the pre-change (scalar shrinkage) library |
| `lambda_rosenbrock.py` | NumCosmo's Rosenbrock and k independent pairs (d = 2k), exact sampler; equal vs NNLS, shrink on/off (`SHRINK`, `HGRID`, `N` env) |
| `lambda_highd.py` | correlated Gaussian at d = 50 and 100, equal vs NNLS, plus the single-Gaussian ceiling |
| `var_trace.py` | Var(-2lnL) per ensemble iteration over 2n for a list of catalogs, and the iteration at which it settles |

Always report from `summarize.py`, not from a bare autocorrelation call: tau is only
meaningful after the Constant Break cut, and needs ~130 iterations beyond it.

## Rerunning

```sh
# real cosmology, 15 and 23 parameters (cheap, ~56 likelihood evals/s per core)
python apes_cosmo.py wspline8  1500 300 cauchy 1.0 0 0.3 1 base
python apes_cosmo.py wspline8  1500 300 gauss  3.0 1 0.3 1 shrink
python apes_cosmo.py wspline16 3000 250 cauchy 1.0 0 0.3 1 base23
python apes_cosmo.py wspline16 3000 250 gauss  3.0 1 0.3 1 shrink23
python summarize.py "cosmo_*.fits"

# synthetic, for the dimension scaling
python apes_hard.py gauss 50 12000 200 cauchy 1.0 0 0.2 1
python apes_hard.py gauss 50 12000 200 gauss  3.0 1 0.2 1

# Planck 2018 TT + LCDM, 21 parameters, ~9.3 s per call on a laptop core
OMP_NUM_THREADS=<ncores> python planck_phaseA.py 300 4 0.35
OMP_NUM_THREADS=<ncores> python planck_phaseB.py <catalog.fits> 300 1 0.35 gauss 3.0 1 shrink
```

Argument order for `apes_cosmo.py` and `apes_hard.py`:
`<target> <d|dataset> <nwalkers> <niter> <kernel> <over_smooth> <shrink 0|1> <local_frac> <interp 0|1> [tag]`.

## Command line

The whole benchmark now runs through the CLI, which exercises the same code path as a
production run instead of the Python `APES` wrapper:

```sh
numcosmo generate sampler-test gauss25.yaml --target gauss-constraint --dim 25
numcosmo run mcmc apes gauss25.yaml --nwalkers 10000 --nsamples 200 --parallel threads \
    --interpolation-kernel gauss --over-smooth 4.0 --local-fraction 0.4 --center-shrink \
    -o gauss25_shrink.fits
numcosmo catalog analyze gauss25_shrink.fits
```

`generate sampler-test` writes any of the four synthetic targets used to benchmark the
sampler (`gauss-constraint`, `funnel`, `rosenbrock`, `gaussmix2d`); `gauss-constraint`
reproduces the published target exactly, same seed and same covariance. `--center-shrink`
is exposed on `numcosmo run mcmc apes` and on every experiment runner in
`numcosmo_py.experiments`.

## Left unfinished

* **d = 50 at the published walker count.** The published runs used 30000 walkers at
  d = 50; the synthetic runs here used 2400, which handicaps the baseline. The published
  rows above d = 40 have acceptance below 0.5 % and tau that scatters by a factor of 10
  between neighbouring walker counts, so they are not converged either and cannot serve as
  a baseline as they stand.
* 3000-walker runs at 15 and 23 parameters, killed mid-run.
* The cosmology baselines (`apes_cosmo.py`) were all taken at Cauchy over-smooth 1. The
  d = 25 sweep showed h = 1 is the optimum for every kernel, so the bandwidth part of the
  concern is settled, but ST3 and a larger local_frac were worth 1.6x there and have not
  been tried on `wspline8`.
* Planck: set up and costed, never started.
* `apes_burnin.py` never run. Shrinkage removes the over-dispersion that helps a collapsed
  ensemble expand, so it could slow burn-in; the Constant Break cuts measured at d = 15
  went the other way (71 with shrinkage against 106 without), but that is one target and
  not a test.

## Plan

`plan.md` lists the queued studies and their evidence.

## Next feature

`auto_tuning.md` has the design. Short version: the block that is not being updated is a
free validation set with its `-2 ln L` already computed, so the bandwidth can be chosen by
maximising a out-of-sample estimate of the acceptance itself,

    A(h) = (1/m) sum_j sum_k w_k min(1, e^(r_j - r_k)),   w_k proportional to e^(r_k)

with `r_j = ln p~(x_j) + 0.5 m2lnL_j`. Measured against the one-step acceptance over 18
bandwidth scans it loses a median of 2.3 % of the attainable acceptance against 7.1 % for
the variance objective the review originally proposed, and its worst case is 30 % against
92 %. Caching the three h-independent Mahalanobis pieces makes a 20-point scan cost a few
per cent of one preparation, so it can run at every half step. The same scan selects the
kernel, which is the other choice that has to follow the target.

The anisotropic shrinkage matrix `A = U_C^T U_M^{-T}` is implemented (2026-09-15); measured
against the scalar on five targets it changes acceptance by less than one standard error.
See `equal_weights_vs_nnls.md`, section 7.
`bench/`: the sampler benchmark protocol (`PROTOCOL.md`), its decisions log (`DECISIONS.md`),
`bench_a.py` (acceptance benchmark, out-of-sample acceptance) and `bench_b.py` (convergence benchmark, short chains through the
CLI: fixtures, three chain states, tau and settle iteration). Tables live in
`~/data/apes/bench/decisions/`.
