# APES review (2026-09-12)

Code reviewed: `ncm_fit_esmcmc_walker_apes.c`, `ncm_stats_dist.c`, `ncm_stats_dist_kde.c`,
`ncm_stats_dist_vkde.c`, `ncm_stats_dist_kernel_st.c`, `ncm_nnls.c`, and the acceptance
step in `ncm_fit_esmcmc.c`. Paper: arXiv:2303.13667v2.

Benchmark script: `apes_interp_bench.py` (this directory). It draws an i.i.d. sample of
size n from the target (the converged half ensemble), builds the NcmStatsDist proposal
exactly as the walker does, and measures on independent points:

* `std_r`: scatter of ln(pi_tilde/pi) at points x ~ pi,
* `acc`: expected acceptance of the independence step, E[min(1, pi(y) pt(x) / (pi(x) pt(y)))],
  x ~ pi, y ~ pi_tilde. This is the quantity APES optimises implicitly.

All CSVs are in this directory; columns are self-describing.

## 1. Correctness of the sampler itself

* MH ratio, forward/backward proposal densities, truncated-normal random-walk mixture
  normalisation, block structure: all correct. No detailed-balance issue found.
* `random-walk-scale` property is stored but never used; `_prepare_random_walk` hard-codes
  0.25 (`ncm_fit_esmcmc_walker_apes.c:685-688`). Bug.
* `_setup` fills `m2lnL_s0`/`m2lnL_s1` twice (dead first loop, T = 1).
* `exploration` mode drops the proposal correction (greedy). Fine for burn-in only; should
  be documented as "not a valid MCMC step".
* VKDE size guard requires only k >= 2 neighbours. With k <= d the local covariance is
  singular and nearPD fallbacks produce garbage: d = 10, local_frac 0.02, n = 500 gives
  k = 10 and acceptance 0.000 (E_localfrac10.csv). Guard should be k >= ~2d, or use a
  shrinkage covariance estimator (see 4).

## 2. What the benchmark says about the current configuration

Gaussian target, d = 10, n = 500 (1000 walkers), local_frac 0.05, interpolation on
(A_gauss10_interp.csv). Acceptance as a function of h and kernel:

| h    | Cauchy | ST3  | ST10 | Gauss |
|------|--------|------|------|-------|
| 0.2  | 0.05   | 0.01 | 0.00 | 0.00  |
| 0.5  | 0.13   | 0.10 | 0.02 | 0.01  |
| 1.0  | 0.19   | 0.28 | 0.30 | 0.25  |
| 1.5  | 0.18   | 0.30 | 0.37 | 0.39  |
| 2.0  | 0.15   | 0.20 | 0.18 | 0.17  |

Rosenbrock, n = 160 (D_rosen_interp.csv): best Cauchy 0.49 (h 0.7), ST3 0.58, Gauss 0.60.
Funnel 10D, n = 1500 (D_funnel10_*.csv): best with interpolation 0.11 / 0.18 / 0.20;
without interpolation 0.16 / 0.20 / 0.20 (Cauchy / ST3 / Gauss).

Observations:

1. Cauchy is never the best kernel once h is tuned; it is the most *forgiving* kernel.
   Its acceptance varies by 4x over h in [0.1, 2], Gauss varies by 100x. With a fixed,
   untuned h (0.2 in the Python wrapper) Cauchy wins. With a tuned h, ST3/Gauss win by
   20-60 %.
2. Cauchy wastes proposals. `std_rY` (scatter of ln(pt/pi) at proposed points) is 1e4-1e6
   for Cauchy versus 1-10 for Gauss: the single radial scale sqrt(nu/chi2_nu) with nu = 1
   sends ~40 % of proposals to >= 2 sigma in *every* direction simultaneously, which is a
   certain rejection for d >= 5. This cost is independent of interpolation quality.
3. The NNLS interpolation adds little on Gaussians (0.19 -> 0.28 for ST3 at h = 1) and
   hurts on the funnel: the fit keeps only 3 % of the weights non-zero (nz = 0.031) for
   h <= 0.75. The collocation problem is ill-posed where kernels overlap strongly and
   the non-negativity constraint then zeroes most of the solution.
4. Above d ~ 25 nothing matters: std_r is ~2.5 (d = 25, n = 1000) and ~3.3 (d = 50,
   n = 2000) for every kernel and every h, acceptance <= 0.02. This is the paper's
   "acceptance below 5 % at d = 50". The error is set by n, not by h: the family
   "sum of n kernels centred on the sample" cannot represent a 25D Gaussian.
5. local_frac: larger is monotonically better on Gaussians (d = 10, ST3, h = 1: 0.27 at
   0.05, 0.37 at 0.1, 0.41 at 0.2, 0.45 at 0.4). The default 0.05 is too small for d >= 10.
6. `shrink` (mixing with uniform weights) has no measurable effect (F_shrink10.csv).
7. Automatic h works: CV_SPLIT_NOFIT with split_frac 0.5 picks h ~ 1.1-1.2 and gets
   0.24 (ST3) / 0.31 (Gauss) with half the kernels; CV_SPLIT (levmar on the interpolation
   residual) gets 0.34 for Gauss. Both are within 10-20 % of the best hand-tuned value.

## 3. The big win: shrink the kernel centres (West 1993 / Liu-West)

A KDE with bandwidth h has covariance (1 + h^2 s^2) Sigma, where s^2 Sigma is the
average kernel covariance. In high d that inflation alone destroys the acceptance.
The classical fix is to move each centre toward the mean,

    x_i' = mu + a (x_i - mu),   a = 1 / sqrt(1 + h^2 s^2),

so that the mixture covariance equals the ensemble covariance for every h. As h -> large
the proposal becomes a single Gaussian/Student-t with the ensemble covariance; as h -> 0
it is the ordinary KDE. Nothing else in the sampler changes (the mixture density is still
exact, MH stays exact). Implemented in the benchmark as `--cshrink` (two passes: first
build measures s^2 from the local covariances, then rebuild on the shrunk centres).
Results, uniform weights (no interpolation):

| target, n            | kernel | best acc, current | best acc, centre-shrink |
|----------------------|--------|-------------------|-------------------------|
| Gauss d=10, 400      | Gauss  | 0.30 (h 1)        | 0.75 (h 3)              |
| Gauss d=10, 400      | ST3    | 0.21 (h 1)        | 0.54 (h 2)              |
| Gauss d=25, 1000     | Gauss  | 0.02 (h 1)        | 0.55 (h 3)              |
| Gauss d=25, 1000     | ST3    | 0.02 (h 1)        | 0.28 (h 3)              |
| Gauss d=50, 2000     | Gauss  | 0.00              | 0.36 (h 3, still rising) |
| Gauss d=50, 2000     | ST3    | 0.00              | 0.13 (h 3)              |
| Rosenbrock, 160      | ST3    | 0.51 (h 0.7)      | 0.50 (h 0.7)            |
| Funnel d=10, 1500    | Gauss  | 0.20 (h 1)        | 0.26 (h 1.5)            |
| Funnel d=10, 1500    | ST3    | 0.20 (h 1)        | 0.25 (h 1.5)            |

Files: G_vkde_d*.csv, G_vkde_d*--cshrink.csv, G_rosen_cshrink_*.csv, G_funnel_*.csv.
The non-Gaussian targets do not get worse (the optimum h moves), the Gaussian-like ones
improve by 2-25x. Cosmological posteriors are much closer to the Gaussian column.

Caveat: combining centre shrinkage with the NNLS interpolation at large h is bad
(Rosenbrock, h >= 1: 1 % non-zero weights, acceptance 0.01). Use uniform weights with
the shrinkage, or fit the weights at the *original* points (needs a C change, see 5).

## 4. Recommendations, in order of expected payoff / effort

1. **Centre shrinkage in NcmStatsDistKDE/VKDE** (new property, default on). ~50 lines:
   compute mu, s^2 in `prepare_kernel`, store shrunk centres for evaluation and sampling.
   Expected: high-d acceptance from ~0 to 0.3-0.5 for Gaussian-like posteriors.
2. **Automatic bandwidth from the other half of the ensemble.** The walker already holds
   the block not being updated *with its -2 ln L values*: an independent validation set
   for free. Minimise Var[ln pt(x_j) - ln pi(x_j)] over the validation block (a constant
   offset does not matter for MH). For VKDE the chi^2 matrix scales as 1/h^2, so store
   the unscaled chi^2 once and re-apply only the kernel function per trial h; the scan is
   then O(n^2) per trial instead of O(n^2 d^2). Replaces the user-tuned `over_smooth`
   (0.2, 1.1, 2.0 in the three notebooks) and makes the light-tailed kernels safe.
3. **Kernel choice.** Make the default ST with nu depending on d (nu ~ 3 for d <= 5,
   nu ~ d for large d), or add a defensive-mixture kernel: K = (1 - eps) K_Gauss +
   eps K_Cauchy with eps ~ 0.05-0.1, same covariance. Keeps the bounded pi/pt ratio that
   makes Cauchy robust, without sending 40 % of the proposals into the tails. Trivial as
   a NcmStatsDistKernel subclass (sample by picking a component).
4. **Local covariance.** Raise the guard to k >= 2d, scale the default local_frac so that
   k >= max(0.05 n, 2d), and shrink each local covariance toward the global one
   (Ledoit-Wolf style, C_i = (1 - lambda) S_i + lambda s_i^2 Sigma_global). Removes the
   nearPD fallback path and the singular-covariance failure mode.
5. **Interpolation fit.** Three cheap changes to `_ncm_stats_dist_prepare_interp`:
   a. Use n_kernels < n_obs by default (split_frac ~ 0.5): the overdetermined NNLS is
      better conditioned, does not collapse to 3 % of the weights, halves the cost of
      sampling/evaluation, and lost only ~10 % acceptance in the tests.
   b. Switch NNLS from `NCM_NNLS_UMETHOD_NORMAL` (Cholesky of A^T A, squares the
      condition number) to QR or DGELSD, both already implemented in NcmNNLS.
   c. Warm-start the NNLS active set from the previous iteration (the ensemble changes
      slowly); currently `Pset` is reset to "all columns" each solve, making the first
      unconstrained solve an O(n^3) dense factorisation every half step.
   Also consider a regularised objective ||K w / pi - 1||^2 + lambda ||w - 1/n||^2 in
   place of the post-hoc `shrink` mixing, and add the fit residual (already computed) to
   the run log as an approximation-quality diagnostic.
6. **Fix `random-walk-scale`** (unused property) and remove the duplicated loop in
   `_setup`.
7. **Cost.** For d >= 20 the Cauchy/ST kernel K(r)/K(0) = (1 + r^2/nu)^-(nu+d)/2 is
   effectively compact (1e-25 at r = 3, d = 50). The IM and the evaluation sums can be
   truncated with the kd-tree that VKDE already builds, turning O(n^2 d^2) into
   O(n k d^2). Only worth doing after 1-3.

## 5. Ideas beyond the current framework

* **Use the history.** Each iteration discards L/2 evaluated (x, pi(x)) pairs. Pooling
  the last m iterations of the opposite block gives m x more points for the fit. Validity
  needs care: proposals then depend on past states, so either freeze the pooled proposal
  after burn-in (standard adaptive-MCMC argument, diminishing adaptation) or use the pool
  only for the scalar bandwidth/shrinkage choices, which is unproblematic.
* **Semi-parametric proposal.** With centre shrinkage, h -> large gives a Student-t fitted
  to the ensemble. Going one step further, fit a Gaussian mixture with few components
  (k-means or EM on the ensemble, e.g. sqrt(n) components with full covariances) and use
  the kernel mixture only as the small-h limit. In d = 50 with n = 2000, 40 components x
  50 points each is a well-posed covariance estimate; 2000 kernels x 100 neighbours is not.
* **Normalising-flow proposal** trained on the ensemble (as in flowMC / pocoMC) is the
  modern alternative for d > 30. It is a large implementation effort and loses the
  "no training loop" property of APES; the items above should be exhausted first.

## 6. Implementation and end-to-end validation (branch apes-center-shrink)

Implemented as NcmStatsDist:center-shrink (contract centres by a about the mean and
apply bandwidth a h), exposed on the APES walker. Native results reproduce the emulation
exactly (N_vkde_d10_native.csv vs N_vkde_d10_python.csv). Fitting the NNLS weights at
the original points (what the C code does) removes the collapse seen when the emulation
fitted at the shrunk centres.

End-to-end chains (apes_e2e.py, correlated Gaussian, stationary start, tau from the
catalog), 600 walkers at d = 10, 1600 walkers at d = 25, local_frac 0.2 at d = 25:

| d  | kernel | h | shrink | acceptance | tau_mean | tau_max |
|----|--------|---|--------|-----------|----------|---------|
| 10 | Cauchy | 1 | off (today's default) | 0.116 | 15.2 | 19.6 |
| 10 | Cauchy | 1 | on  | 0.078 | 27.0 | 41.9 |
| 10 | Cauchy | 2 | on  | 0.196 | 11.9 | 19.7 |
| 10 | Gauss  | 1 | off | 0.105 | 23.0 | 48.7 |
| 10 | Gauss  | 2 | on  | 0.316 |  6.1 | 11.2 |
| 10 | ST3    | 2 | on  | 0.298 |  6.2 |  8.8 |
| 25 | Cauchy | 1 | off | 0.026 | 26.9 | 56.7 |
| 25 | Cauchy | 1 | on  | 0.038 | 19.8 | 49.4 |
| 25 | Gauss  | 2 | off | 0.008 | 40.6 | 91.8 |
| 25 | Gauss  | 2 | on  | 0.056 | 25.8 | 50.2 |
| 25 | ST3    | 2 | on  | 0.057 | 21.2 | 40.8 |

(d = 25 chains have only 120 iterations, so tau there is indicative, not converged.)

Lessons that the one-step benchmark could not show:

* Shrinkage moves the optimum bandwidth up (h ~ 2 in units of the local covariance).
  Switched on with the old defaults (Cauchy, h = 1) it makes the chain worse, so the
  property is off by default; pairing it with h = 2 and a Gauss/ST3 kernel gives
  2.5x lower tau at d = 10.
* The pre-existing runs started far from stationarity because the Python APES wrapper
  ignored its initial_sample argument (inserted after start_run had already generated
  points). Fixed on the branch. Burn-in ensembles trigger the m2lnp cut and, with
  narrow kernels, produce all-zero interpolation columns that the NNLS did not tolerate
  (two crash paths, both fixed on the branch).
* Stuck walkers, not bulk acceptance, dominate tau_max: a walker sitting where
  pi_tilde << pi cannot leave. Heavy tails limit this, which is the real reason the
  Cauchy kernel was the safe choice. A delayed-rejection local move would remove the
  failure mode independently of the kernel.

## 7. Extended validation, focused on high dimension (branch apes-center-shrink)

Script `apes_validate.py`. Three independent checks per configuration, over a grid of
target (Gaussian, funnel), d (5, 10, 25, 50), method (KDE, VKDE), kernel (Gauss, ST3,
ST10), bandwidth, and interpolation on/off:

1. **exact** : the mixture density is rebuilt in numpy from the public API
   (peek_center_array, peek_cov_decomp, peek_weights, get_lnnorm, get_href,
   get_center_shrink_factor) and compared with eval() and eval_m2lnp().
2. **ks** : draws are projected on random directions and compared with the exact
   one-dimensional law of the mixture along that direction (marginals of Gaussian and
   Student-t kernels keep their family). This tests that sample() draws from the density
   eval() returns; an inconsistency there breaks detailed balance. A global
   importance-sampling normalisation check was tried first and is useless above
   d ~ 10, its variance swamps the signal.
3. **cov** : covariance of the draws against the sample covariance.

Results over 88 configurations:

* exact: **0 failures**, max relative error 1e-14 at d = 25 and 6e-14 at d = 50.
* ks: 4 marginal exceedances of the 1 % threshold out of 264 direction tests, 3 of them
  on the shrink-off (unmodified) path. Consistent with the nominal false-positive rate.
* cov: exposed a real bug, see below.

### Bug found: the kernel variance factor was missing

The mixture covariance is $a^2 (1 + \kappa h^2 s^2) \Sigma$ where $\kappa$ is the ratio
between the kernel covariance and its scale matrix. The first implementation assumed
$\kappa = 1$, true only for the Gaussian kernel. Measured ratio of mixture covariance to
sample covariance with shrinkage on, before the fix:

| kernel | kappa | d = 5 | d = 25 |
|--------|-------|-------|--------|
| Gauss  | 1     | 1.00  | 1.00   |
| ST3    | 3     | 3.69  | 2.01   |
| ST10   | 1.25  | -     | 1.11   |

After adding `ncm_stats_dist_kernel_get_var_factor()` all kernels land at 1.00 +- 0.05.
The Cauchy kernel has no covariance at all ($\kappa = \infty$), so centre shrinkage is
undefined for it and is now refused with an explicit error in both NcmStatsDist and the
APES walker. That also explains the earlier end-to-end observation that shrinkage with
the default Cauchy kernel made the chain worse.

Related fix: `ncm_stats_dist_get_href()` returned the nominal bandwidth while
`get_lnnorm()` and `get_Ki()` used the applied one; with shrinkage the two differ by the
factor a. The getter now returns the applied bandwidth.

### Proposal quality in high dimension (one-step acceptance, i.i.d. training sample)

Gaussian target, n = 40 d training points, local_frac 0.2, uniform weights:

| d  | kernel | best without shrinkage | best with shrinkage |
|----|--------|------------------------|---------------------|
| 25 | Cauchy | 0.011 (h 1)            | not applicable      |
| 25 | ST3    | 0.018 (h 1)            | 0.173 (h 1.5)       |
| 25 | Gauss  | 0.019 (h 1)            | **0.572 (h 4)**     |
| 50 | any    | 0.000                  | -                   |
| 50 | ST3    | 0.000                  | 0.085 (h 2)         |
| 50 | Gauss  | 0.000                  | **0.387 (h 4)**     |

Interpolation changes these by less than 10 % either way at high d; the NNLS is not the
limiting factor there.

Funnel target, d = 25, n = 1500:

| kernel | without shrinkage | with shrinkage |
|--------|-------------------|----------------|
| Cauchy | 0.029 (h 1.5)     | not applicable |
| ST3    | 0.032 (h 2)       | 0.068 (h 3)    |
| Gauss  | 0.024 (h 2)       | 0.034 (h 3)    |

Two lessons beyond the Gaussian case:

* On the funnel the shrinkage factor stays near 1 (0.68 to 0.98) because the local
  covariances are far smaller than the global one, so shrinkage is nearly inert there.
  It does not distort a strongly non-Gaussian target, and it does not rescue one either.
* The Gaussian kernel is disastrous on the funnel (scatter of ln(pi_tilde/pi) of 130 to
  3400 versus 3 to 6 for ST3). Kernel choice has to follow the target: light tails plus
  shrinkage for Gaussian-like posteriors, heavy tails for funnels.
* Student-t kernels saturate in high dimension even with the covariance matched (0.17 at
  d = 25, 0.085 at d = 50). The single radial scale sqrt(nu / chi2_nu) multiplies every
  direction at once, and the density ratio picks up that factor to the power d.

### Pre-existing pathology confirmed, quantified

With interpolation on, the funnel mixture covariance is 6 % to 29 % of the sample
covariance at d = 10 and 25 (cov_frob ~ 0.9, i.e. essentially the wrong matrix). The NNLS
concentrates the weight in the neck of the funnel, and at d = 10, h = 2 the proposal
density underflows to zero at typical target draws. This is on the shrink-off path and is
independent of this branch's changes.

## 8. Two pre-existing walker bugs found while chasing a benchmark discrepancy

The one-step benchmark predicted 0.57 acceptance for Gauss + shrinkage at d = 25 while the
chain measured 0.19. Scoring the walker's OWN prepared proposal object (peek_sds) gave
0.364 against a chain acceptance of 0.367 at d = 10, while an independently rebuilt
NcmStatsDist with nominally identical settings gave 0.61. The walker was therefore not
using the settings it had been given.

1. `ncm_fit_esmcmc_walker_apes_set_method()` and `set_k_type()` call `_set_sys()`, which
   destroys and recreates both NcmStatsDist objects. `over_smooth`, `shrink`, `use_threads`
   and `center_shrink` are reapplied afterwards; **`local_frac` and the covariance type are
   not**, so they silently revert to 0.05 and COV_TYPE_SAMPLE. `create_esmcmc()` and the
   Python APES class both call `set_local_frac` and `set_cov_robust` BEFORE `set_k_type`,
   so every run with a non-default kernel or method used local_frac 0.05 and a
   non-robust covariance whatever the user asked for. At d = 25 with 1600 walkers that is
   40 neighbours for a 25-dimensional local covariance, i.e. the singular regime flagged in
   section 1. Runs with the default Cauchy kernel are unaffected, because then `set_k_type`
   does not change the method/kernel pair and no rebuild happens.
2. `NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE` created a `ncm_stats_dist_vkde_new()`, so the
   fixed-bandwidth method ran the variable-bandwidth estimator while `desc()` reported
   "KDE". Present since commit 372c802e (#92). Any KDE vs VKDE comparison made through the
   walker was VKDE against itself. The notebooks that build NcmStatsDistKDE directly are
   not affected.

Both fixed on the branch: the walker now stores local_frac, cov type and fixed covariance
and reapplies them after a rebuild, and METHOD_KDE builds a NcmStatsDistKDE.

All end-to-end numbers in section 6 for the ST3 and Gaussian kernels were taken with
local_frac 0.05 instead of the requested 0.2, so they understate the benefit of centre
shrinkage; the Cauchy baselines are unaffected. Section 9 repeats them on the fixed code.

## 9. End-to-end on the fixed code (d = 25 Gaussian, 1600 walkers, 200 iterations)

Stationary start, tau from the catalog, local_frac 0.2 now actually applied.

| kernel | h | shrink | interp | acceptance | tau_mean | tau_max |
|--------|---|--------|--------|-----------|----------|---------|
| Cauchy | 1 | off    | on     | 0.026 | 33.9 | 102.5 |
| ST3    | 2 | on     | on     | 0.153 | 12.3 |  26.1 |
| Gauss  | 3 | on     | on     | 0.485 |  4.0 |  10.9 |
| Gauss  | 3 | on     | off    | 0.499 |  3.6 |   5.5 |
| Gauss  | 4 | on     | on     | 0.526 |  3.0 |   4.3 |

Against the current default (first row): **11x lower mean autocorrelation and 24x lower
worst-walker autocorrelation**. The chain acceptance (0.50) now agrees with the one-step
benchmark (0.57); the earlier 3x gap was entirely the local_frac bug of section 8.

## 10. Harder targets, end to end (d = 10)

| target   | config                        | acceptance | tau_mean | note |
|----------|-------------------------------|-----------|----------|------|
| bimodal  | Cauchy h=1, no shrink (default) | 0.121 | 14.4 | modes resolved, 48 % of walkers change mode |
| bimodal  | Gauss h=3, no shrink            | 0.013 | 49.8 | proposal 9x too wide, mode changes drop to 12 % |
| bimodal  | Gauss h=3, shrink               | 0.286 |  5.8 | modes resolved, 51 % change mode |
| bimodal  | Gauss h=3, shrink, uniform weights | 0.310 |  7.0 | same as with interpolation |
| funnel   | Cauchy h=1, no shrink (default) | 0.044 | 90.3 | tau_v 35 |
| funnel   | ST3 h=1.5, no shrink            | 0.121 | 41.5 | tau_v 36 |
| funnel   | ST3 h=1.5, shrink               | 0.149 | 41.9 | tau_v 21, a = 0.78 |
| banana   | Cauchy h=1, no shrink (default) | 0.145 | 13.6 | |
| banana   | Gauss h=3, shrink               | 0.239 | 10.5 | |
| banana   | ST3 h=2, shrink                 | 0.131 | 17.9 | WORSE than the default |

Bimodal at d = 25 (1200 walkers): default 0.017 acceptance, tau 51.3, only 26 % of walkers
change mode; Gauss h=3 with shrinkage 0.193, tau 11.2, 48 % change mode. Shrinkage improves
both the efficiency and the mode mixing, and both keep the modes correctly placed.

Two modes separated by 8 sigma, so the global covariance is dominated by the separation
(Sigma_00 = 17 against a within-mode variance of 1). This is the case where scalar centre
shrinkage is least justified, and the measured factor is indeed aggressive, a = 0.386,
contracting the mode centres from +-4 to +-1.5. The chain nevertheless samples the target
correctly (mean |x0| 3.98 against a true 4.02, gap occupancy 1.3 % against a true 2.2 %)
and mixes between modes better than the default. Two reasons: MH corrects the proposal
exactly, and the interpolation weights partly compensate for the contraction by
up-weighting the outermost kernels.

The scalar factor matches the covariance only in TRACE. When the local covariances have a
very different shape from the global one (bimodal, funnel) the match is not per-direction:
for the bimodal target the shrunk mixture is about twice too narrow along the separation
and correspondingly too wide in the other nine directions. The natural fix is an
anisotropic shrinkage matrix A solving A (Sigma + kappa h^2 Cbar) A^T = Sigma, that is
A = Sigma^{1/2} (Sigma + kappa h^2 Cbar)^{-1/2}, at the cost of one extra d x d
factorisation per preparation. Left as future work.

### What the hard targets say overall

Reduction in mean autocorrelation against the current default (Cauchy, h = 1, interpolation):

| target                | factor | best configuration |
|-----------------------|--------|--------------------|
| Gaussian d = 25       | 11x    | Gauss h = 4, shrink |
| bimodal d = 25        | 4.6x   | Gauss h = 3, shrink |
| bimodal d = 10        | 2.5x   | Gauss h = 3, shrink |
| funnel d = 10         | 2.2x   | ST3 h = 1.5 (shrinkage adds little beyond the kernel change) |
| banana d = 10         | 1.3x   | Gauss h = 3, shrink |

The gain is largest exactly where cosmological posteriors live, close to Gaussian and in
many dimensions, and it degrades gracefully on curved and multimodal targets. It is not
universal: ST3 with shrinkage on the banana is worse than the default (tau 17.9 against
13.6), so the kernel must be chosen for the target rather than switched on blindly. On the
funnel most of the gain comes from replacing Cauchy with ST3, not from shrinkage, which is
nearly inert there (a = 0.78).
