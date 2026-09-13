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

Bug fixes that are independent of the feature and worth keeping regardless:

* `set_method` / `set_k_type` silently dropped `local_frac` and the covariance type.
* `METHOD_KDE` built a VKDE.
* `NcmNNLS` aborted on an empty passive set and on rank-deficient QR.
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

## Left unfinished

Do the first item before quoting any number from `review.md` anywhere else.

* **The fairness control, and it is the one that can move the conclusions.** Every
  comparison here is against the Cauchy kernel at over-smooth 1. The notebooks use 1.1
  and 2.0, so the baseline may be mis-tuned and the reported ratios too generous. Sweep
  the baseline's over-smooth at d = 25 (and on `wspline8`) and requote against the best
  Cauchy, not against h = 1. If a tuned baseline reaches, say, tau 10 instead of 34 at
  d = 25, the honest factor is 3x rather than 11x.
* 3000-walker runs at 15 and 23 parameters, killed mid-run.
* d = 50 synthetic at a walker count matching the paper's scaling (it used 30000 at
  d = 50; the runs here used 2400, which handicaps the baseline).
* Planck: set up and costed, never started.
* `apes_burnin.py` never run. Shrinkage removes the over-dispersion that helps a
  collapsed ensemble expand, so it could slow burn-in; the Constant Break cuts measured
  at d = 15 went the other way (71 with shrinkage against 106 without), but that is one
  target and not a test.

## Next feature

The walker hardcodes `NCM_STATS_DIST_CV_NONE` at all four construction sites, while
`numcosmo catalog calibrate` already auto-tunes the over-smooth through
`CrossValidationMethod.SPLIT_NOFIT`. Exposing the cross-validation type on the walker is
most of the automatic-bandwidth work, and without it the gains measured here need
hand-tuning to h ~ 3. After that, the anisotropic shrinkage matrix
`A = Sigma^(1/2) (Sigma + kappa h^2 Cbar)^(-1/2)` is the natural refinement: the scalar
factor matches the covariance only in trace, which is visible on the bimodal target.
