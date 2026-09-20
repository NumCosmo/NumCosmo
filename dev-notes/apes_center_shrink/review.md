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

## 11. The fairness control (2026-09-13)

Every number in sections 9 and 10 was quoted against the Cauchy kernel at over-smooth 1,
which is the library default but not necessarily a well-tuned baseline. Sweeping the
baseline (d = 25 Gaussian, 1600 walkers, 200 iterations, interpolation on, no shrinkage,
same binary and same session as the shrinkage runs):

| kernel | h    | local_frac | acceptance | tau_mean | tau_max |
|--------|------|-----------|-----------|----------|---------|
| Cauchy | 0.5  | 0.2 | 0.019 | 40.7 |  92.1 |
| Cauchy | 1.0  | 0.2 | 0.025 | 31.7 |  77.3 |
| Cauchy | 1.5  | 0.2 | 0.024 | 33.8 |  88.4 |
| Cauchy | 2.0  | 0.2 | 0.019 | 43.4 | 115.1 |
| Cauchy | 3.0  | 0.2 | 0.008 | 65.3 | 123.5 |
| Cauchy | 1.0  | 0.4 | 0.033 | 41.3 | 101.5 |
| ST3    | 0.5  | 0.2 | 0.016 | 37.2 |  68.2 |
| ST3    | 0.75 | 0.2 | 0.031 | 26.1 |  48.3 |
| ST3    | 1.0  | 0.2 | 0.044 | 22.5 |  39.6 |
| ST3    | 1.5  | 0.2 | 0.035 | 30.3 |  63.0 |
| ST3    | 2.0  | 0.2 | 0.019 | 39.1 |  87.2 |
| ST3    | 3.0  | 0.2 | 0.005 | 55.2 | 141.8 |
| ST3    | 1.0  | 0.4 | 0.055 | 21.8 |  34.1 |
| ST3    | 1.0  | 0.6 | 0.060 | **20.0** |  42.4 |
| Gauss  | 0.5  | 0.2 | 0.004 | 80.2 | 216.2 |
| Gauss  | 0.75 | 0.2 | 0.008 | 34.0 |  77.6 |
| Gauss  | 1.0  | 0.2 | 0.054 | 30.0 |  72.6 |
| Gauss  | 2.0  | 0.2 | 0.004 | 76.7 | 184.9 |
| Gauss  | 3.0  | 0.2 | 0.004 | 67.3 | 156.9 |
| Gauss  | 1.0  | 0.4 | 0.095 | 23.2 |  46.3 |

and the shrinkage configurations re-measured in the same batch:

| kernel | h | local_frac | acceptance | tau_mean | tau_max |
|--------|---|-----------|-----------|----------|---------|
| Gauss  | 3 | 0.2 | 0.486 | 3.5 | 6.0 |
| Gauss  | 4 | 0.2 | 0.525 | 3.4 | 5.1 |
| Gauss  | 4 | 0.4 | 0.580 | **3.0** | **4.5** |

Conclusions:

* **h = 1 is the optimum for every kernel without shrinkage.** The baseline was not
  mis-tuned in h; the 4x span in tau across h in [0.5, 3] is real but h = 1 sits at the
  minimum. The notebooks' 1.1 is the same point.
* The baseline *is* improved by the other two knobs: ST3 instead of Cauchy (31.7 -> 22.5)
  and a larger local_frac (22.5 -> 20.0). Best tuned baseline: ST3, h = 1, local_frac 0.6,
  tau_mean 20.0; best tau_max 34.1 at local_frac 0.4.
* Against that best tuned baseline the centre-shrinkage configuration gives **6.7x lower
  mean autocorrelation and 7.6x lower worst-walker autocorrelation**, not the 11x and 24x
  quoted against the untuned default. The conclusion survives; the size of the worst-walker
  gain was mostly the untuned baseline.
* Run-to-run scatter on tau_mean at these settings is about 7 % (Cauchy h = 1 gave 33.9 in
  section 9 and 31.7 here), so differences below ~15 % in these tables are not significant.

### The published baseline, from the paper's own runs

The APES paper's `gauss_constraint` catalogs (dimensions 2 to 50, walkers-per-dimension
100 to 600) give the baseline at its own walker count, which is the comparison that matters
for any claim made outside this directory. Reading tau from the stored `.ess` diagnostics:

| d  | walkers | w/d | acceptance | tau_mean | tau_max |
|----|---------|-----|-----------|----------|---------|
| 25 |  2500   | 100 | 0.026 | 106.1 | 232.6 |
| 25 |  5000   | 200 | 0.089 |  21.2 |  50.1 |
| 25 | 10000   | 400 | 0.164 |  17.2 |  22.4 |
| 25 | 15000   | 600 | 0.210 |  10.8 |  15.0 |
| 50 | 15000   | 300 | 0.003 |  99.9 | 276.4 |
| 50 | 20000   | 400 | 0.004 | 111.8 | 478.0 |
| 50 | 30000   | 600 | 0.005 |  36.2 |  72.0 |

Two things to note. First, the acceptance of the published runs at d = 25 with 2500 walkers
(0.026) matches the baseline measured here at 1600 walkers (0.025), so the two setups agree
where they overlap despite the different Gaussian target. Second, at d >= 40 the published
acceptance is below 0.5 % and tau scatters by a factor of 10 between neighbouring walker
counts; those rows are not converged and should not be quoted as a baseline.

The honest cross-comparison at d = 25: the shrinkage configuration reaches tau_mean 3.0
with 1600 walkers, against tau_mean 10.8 for the published configuration with 15000. That
is 3.6x in tau at 9.4x fewer likelihood evaluations per iteration. It is not a like-for-like
run (different Gaussian target, different code revision), which is what the
`numcosmo generate sampler-test` command added on this branch is for.

## 12. Automatic tuning, and what the gain actually depends on (2026-09-13)

### The tuner already existed; the walker could not reach it

`NCM_STATS_DIST_CV_SPLIT_NOFIT` already splits each block into kernel centres and a disjoint
out-of-sample set and minimises the out-of-sample `-2 ln p~` over `ln over_smooth`, with the NNLS
weights fitted afterwards at the chosen bandwidth. The APES walker hardcoded
`NCM_STATS_DIST_CV_NONE` at all four construction sites, so none of it ever ran inside a
chain. That is the whole reason the centre-shrinkage results needed hand-tuning.

Section 4 of this review proposed minimising `Var[ln(p~/pi)]` instead. Measured over 18
bandwidth scans, that is the worst of the three candidates and should not be built:

| objective | picks the grid optimum | median acceptance lost | worst |
|-----------|------------------------|------------------------|-------|
| `mean[r]`, i.e. what `CV_SPLIT_NOFIT` already does | 10 / 18 | 0.0 % | 38.2 % |
| `Var[r]`, proposed in section 4 | 7 / 18 | 7.1 % | 92.2 % |
| out-of-sample importance-sampling estimate of the acceptance | 8 / 18 | 2.3 % | 29.7 % |

### The kernel is a continuous parameter, not a choice

The Student-t kernel `(1 + chi2/nu)^(-(nu+d)/2)` is the Cauchy kernel at `nu = 1` and tends
to the Gaussian one as `nu` grows, so kernel and bandwidth are a single two-parameter fit
over `(ln over_smooth, ln nu)`. Bounds: `nu > 2` with centre shrinkage, since the kernel
covariance is `nu / (nu - 2)`; and `nu <= 1e4` at the other end, where the kernel is already
within 5e-5 of the Gaussian and beyond which the `(1 + chi2/nu)` form loses precision rather
than gaining accuracy. The Gaussian kernel itself is installed when the fit reaches the
ceiling. See `auto_tuning.md`.

### End to end

d = 25, 1600 walkers, 200 iterations:

| configuration | over-smooth | acceptance | tau_mean | time |
|---------------|-------------|-----------|----------|------|
| Cauchy h = 1, the default | 1.0 fixed | 0.024 | 32.2 | 170 s |
| ST3 h = 1, local_frac 0.4, best tuned baseline | 1.0 fixed | 0.056 | 19.9 | 192 s |
| Gauss + shrink, **auto**, local_frac 0.4, split 0.9 | -> 3.04 | 0.564 | **2.7** | 231 s |
| Gauss h = 4 + shrink, local_frac 0.4, hand-tuned | 4.0 fixed | 0.583 | 2.6 | 189 s |

d = 50, 2400 walkers, 150 iterations:

| configuration | over-smooth | acceptance | tau_mean | tau_max | time |
|---------------|-------------|-----------|----------|---------|------|
| Cauchy h = 1, the default | 1.0 fixed | 0.002 | 40.7 | 109.9 | 1674 s |
| ST3 h = 1, local_frac 0.4 | 1.0 fixed | 0.002 | 51.0 | 180.2 | 1870 s |
| Gauss + shrink, **auto**, local_frac 0.4, split 0.8 | -> 4.34 | 0.328 | **5.8** | 11.4 | **1588 s** |
| Gauss h = 4 + shrink, local_frac 0.4, hand-tuned | 4.0 fixed | 0.333 | 4.7 | 10.5 | 1757 s |

The tuner reaches hand-tuned quality without the knob: 2.7 against 2.6 at d = 25, 5.8
against 4.7 at d = 50. It is start-independent, landing on the same bandwidth from h = 1 and
from h = 4. At d = 50 it is also the fastest configuration in the table, because
`split_frac 0.8` builds the approximation from 20 % fewer kernels and that more than pays
for the fit. Against the default that is 7.0x lower tau at 0.95x the wall-clock.

ST3 is worse than Cauchy at d = 50 (51.0 against 40.7), which is the saturation of section 7
showing up end to end: a fixed heavy-tailed kernel is the wrong default in high dimension,
and that is the argument for fitting `nu` rather than choosing it.

### The gain depends on walkers per dimension, and that changes the claim

Repeating the comparison on the published `gauss_constraint` target at its own operating
point, d = 25 with 10000 walkers (400 per dimension), 200 iterations:

| configuration | acceptance | tau |
|---------------|-----------|-----|
| the published configuration, Gauss h = 1.1 | 0.219 | 40.7 |
| the same with centre shrinkage at h = 1.1 | 0.116 | 62.0 |
| Gauss h = 4 + shrink, local_frac 0.4 | **0.327** | **26.2** |
| ST3 h = 1, local_frac 0.4 | 0.204 | 41.5 |

That is 1.55x, not the 6.7x measured at 1600 walkers. These four chains are short (the
Constant Break cut leaves 100 to 161 iterations, R-1 is 0.086 and the ESS is 6), so the
absolute tau values are not converged and only the relative comparison is usable. The
direction is nevertheless clear and it is what the mechanism predicts: the rule-of-thumb
bandwidth falls as `n^(-1/(d+4))`, so the covariance inflation `1 + kappa h^2 s^2` shrinks as
the ensemble grows and there is less for the shrinkage to remove.

**The honest claim is therefore not lower tau at a fixed walker count, it is the same
quality with far fewer walkers.** At d = 25 the shrinkage configuration reaches tau 3.0 with
1600 walkers where the published configuration needs 15000 for tau 10.8. The second row is
also a reminder that shrinkage at the old bandwidth is worse than no shrinkage at all.

### d = 100, beyond where the paper stops

4800 walkers (48 per dimension), 120 iterations, correlated Gaussian:

| configuration | fitted h | acceptance | tau_mean | tau_max | time |
|---------------|----------|-----------|----------|---------|------|
| Cauchy h = 1, the default | 1.0 | **0.000** | 3.3 | 135.5 | 12684 s |
| Gauss + shrink, auto | 5.52 | 0.176 | 11.1 | 49.3 | 20929 s |
| Gauss h = 4 + shrink | 4.0 | 0.137 | 11.6 | 39.6 | 13746 s |
| Gauss h = 6 + shrink | 6.0 | 0.187 | **9.9** | 30.2 | 13731 s |

The default does not sample: acceptance is zero to three decimals. Its `tau_mean` of 3.3 is
an artefact and a warning about reading tau alone, a frozen chain has no autocorrelation to
measure; `tau_max = 135.5` against a 120-iteration chain is the real signal, and
`Var(-2lnL)` would have said "collapsed" immediately.

The fitted bandwidth follows roughly `h ~ sqrt(d)`: 2.0, 3.0, 3.9, 5.5 at `d = 10, 25, 50,
100`, so no fixed default can be right across dimensions. The tuner lands within 12 % of the
best hand-tuned tau without being told anything about the dimension.

Cost of the tuner at `d = 100` is 52 % of wall-clock (20929 s against 13731 s), because each
objective evaluation does a triangular solve per out-of-sample point and kernel. That matters on
a synthetic target and not at all on a real posterior: one Planck likelihood is 2.2 s, which
puts the same tuner at about 1 %.

### Robustness: where the default collapses, the tuned configuration does not

Same target at `d = 10`, varying the ensemble size only:

| walkers | default: Var/2n | default tau | auto: Var/2n | auto tau |
|---------|-----------------|-------------|--------------|----------|
| 100  | 1.53 | 61.8 | 1.05 | 13.7 |
| 200  | **0.00** | 129 | 0.99 | 4.8 |
| 400  | **0.00** | 192 | 1.01 | 2.5 |
| 800  | 1.02 | 10.7 | 1.00 | 1.9 |
| 1600 | 1.00 | 7.8  | 1.00 | 1.4 |

The shipped default collapses the ensemble at 20 to 40 walkers per dimension, which
`Var(-2lnL) -> 0` detects at once while tau is non-monotonic and easy to misread. The tuned
configuration holds `Var/2n = 1.00` down to 10 walkers per dimension and is 5.6x better in
tau at the top of the range. The gain is not only speed; it is that the sampler still works
in an ensemble regime where the default fails outright.
