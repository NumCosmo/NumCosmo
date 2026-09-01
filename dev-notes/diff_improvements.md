# NcmDiff improvements (2026-09, branch diff_improvements)

Three work items, one commit each, on top of a committed baseline benchmark.
Benchmark harness: `dev-notes/diff_benchmark.py`, 15 scalar cases (smooth,
oscillatory, near-singularity, sharp features centered on and off the eval
point, a 1e10 constant offset for cancellation), d1 and d2, recording actual
error against analytic references, the estimated error, honesty
(estimate >= actual), and function evaluations. JSON snapshots per phase:
`diff_benchmark_{baseline,dual,spectral}.json`.

## 1. Reorganization (behavior-preserving)

The per-component convergence state machine, duplicated between the vector
driver and the Hessian driver, is now one `NcmDiffConv` + `_ncm_diff_conv_update`.
The first-order and general Richardson accumulation branches were unified
(identical FP operation order), step evaluation and extrapolation are shared
helpers, and the internal table constructors are static. Verified bit-identical
to master on all rf/rc/Hessian paths (probe script comparing every output value).

## 2. Dual-series scheme (property `dual-series`, default off)

Two Richardson ladders from initial steps h0 and h0/sqrt(rs) share the
extrapolation tables; their same-order disagreement |E_A - E_B| replaces the
consecutive-order difference as the truncation estimate. It measures the error
at the current order instead of lagging one order behind, so the pads drop
from 3e4 to 10.

Details that mattered:
- Keep refining while the round-off is far below the best error (the legacy
  headroom rule): without it, functions whose structure appears only at small
  steps (atan(1000(x-1)) at x=1) stall pre-asymptotically and return garbage
  flagged only by the estimate.
- Return the smaller-step ladder when truncation dominates, the ladder with
  less round-off when noise dominates (recovers the legacy value on the
  offset-1e10 d2 case exactly).

Measured vs single series: rc_d1 geometric-mean error 5x smaller
(1.1e-14 vs 5.8e-14), estimates ~100-1000x tighter (est/act ~10-1e3 vs
~1e4-1e5), all 45 cells honest, +28% evaluations.

## 3. Spectral derivatives (`ncm_diff_sc_*`, property `spectral-window`)

Chebyshev fit on a window [x-R, x+R], derivative taken in coefficient space at
the center. R is found by hill descent on a dyadic ladder, scored by the
estimated derivative error of a fixed-N=17 probe fit; refinement then runs on
nested 17/33/65 Lobatto grids (direct DCT, self-contained; no FFTW/NcmSpectral
dependency).

Design decisions forced by failures observed along the way:
- Convergence is judged on the derivative between refinement levels, never on
  the coefficient norm: a large constant offset dominates the norm while the
  derivative-carrying coefficients are still noise (Sandro's warning; confirmed).
- Truncation is estimated from the tail of the ORIGINAL coefficients propagated
  through the differentiation. The differentiated series of a non-converging
  (step-like) fit is bottom-heavy, so its own tail looks converged while the
  value is garbage: measured on the atan1000 case, which was 100% wrong with a
  confident estimate before this change.
- Resolvedness (q_fit = coefficient mass fraction in the top quarter) gates the
  walk: while no resolved window has been seen, shrinking never gives up, and
  expansion is disallowed. For saturated step data the whole t-space picture is
  window-independent, so the 1/R scaling of the error estimate otherwise makes
  ever-wider windows score better (observed: walked out to R=256 on atan1000).
- Round-off floor eps * sum|c_k| per coefficient, propagated with the same
  recurrence on magnitudes (amplification ~k^2 per derivative order).

Measured (15 cases, d1/d2, all honest): geometric mean comparable to Richardson
central (1.5e-13 vs 5.8e-14 d1); worst case 18x better on d1 (1.0e-5 vs 1.9e-4,
the cancellation case; 20-300x better on offset cases specifically); cost ~5x
evaluations (median ~84 per derivative). The window search is the point: it
finds features at scales far from |x| (atan1000: 3e-13 relative where a fixed
window would return noise).

## numdifftools 0.11.0 comparison (2026-09-01, same battery, defaults both sides)

nd d1: geo-mean 6.3e-11, worst 1.0 (total failure), 31 evals/call, 5/15 cells
with estimate < actual error. nd d2: geo-mean 8.7e-09, worst 1.0, 5/15
dishonest. Failures: inv@1e-3 and sqrt@1e-4 (its default step ladder is not
scaled to |x|, steps across the singularity: NaN or 100% error with confident
estimates), sin100 d2, atan1000 d1, gauss_n d2 (fixed ladder misses the
feature scale). On easy smooth cases it matches NcmDiff (~1e-14) with tight
estimates (est/act ~2-10). All NcmDiff variants: 0 dishonest cells, no value
failures, on this battery. Run via a scratch venv; numdifftools is not
installed system-wide.

## Higher-order derivatives from the spectral fit (prototype, numpy)

One Chebyshev fit yields every derivative order by repeating the coefficient
recurrence; the marginal cost per extra order is zero function evaluations.
Measured (same window-scan logic, N = 65): d3/d4 actual relative errors
exp@1 2e-12/6e-11, sin@1.3 9e-14/1e-12, sin100 2e-13/8e-13, log@1e-2 1e-8/6e-8;
noise-limited offset1e10 degrades to ~1e-2 as expected (amplification ~N^2/R
per order). NcmDiff currently has no d3+ at all, and sc_d1/sc_d2 each redo the
scan+fit: a shared-fit multi-order entry point would halve the d1+d2 cost and
add d^n. The window should be scored on the highest requested order.

## Augur + derivkit assessment (2026-09-01)

Augur (origin/master f8ab7cf; local checkout left untouched on its stale
branch) offers three derivative backends for its Fisher matrices: a five-point
stencil with a FIXED absolute step (default 0.01, not scaled per parameter),
numdifftools Jacobian (same fixed step), and derivkit (LSSTDESC, PyPI 1.2.1)
with defaults n_points=27, spacing='1%', base_abs=1e-3, ridge=1e-8.

derivkit's "adaptive" method is one ridge-regularized least-squares polynomial
of degree order+2 (cap 8) on a Chebyshev grid of fixed half-width (1% of |x0|,
floor base_abs); the returned "error" is the fit's relative RMS residual, in
function units. On the 15-case battery (d1/d2):

- dk defaults: d1 geo-mean 8.2e-6, worst 2.5e+1; d2 worst 5.8e+1; the residual
  proxy is below the actual error in 8/15 (d1) and 13/15 (d2) cells;
  sqrt@1e-4 is NaN (base_abs floor steps across the singularity).
- dk Augur config: the ridge puts a hard floor of ~7e-9 relative on EVERY
  smooth d1 (even an exact cubic); offset1e10 d2 comes out 7.7e+5 relative;
  no scale adaptation (inv@1e-3 is 100% wrong, atan1000 67%).
- Its one strength is realized in the noisy-pipeline regime it was built for:
  with 1e-6 relative white noise it gets ~2e-5..8e-5 honestly at 28 evals.
  But NcmDiff spectral is 3-400x more accurate on the same noisy battery
  (wide window averages the noise), dual stays honest with tight estimates,
  and both add scale discovery and derivative-unit error estimates that
  derivkit lacks. NcmDiff variants: 0 dishonest cells under noise too
  (single's estimate is honest but uselessly loose there).

Summary: the "smoothing polynomial = stability" idea is our spectral method
minus the window search, minus the refinement, minus honest error estimation,
plus a ridge bias floor. Its advantages are fixed low cost (28 evals) and
embarrassingly parallel evaluation.

## Follow-ups (not done)

- Spectral Hessian off-diagonals (needs 2D fits).
- The spectral round-off model over-bounds by ~1e3 (sum of amplified
  magnitudes; the true DCT noise partially cancels). A quadrature-sum model
  would tighten est/act but needs a honesty re-audit.
- numdifftools comparison not re-run (not installed in this env).
- Dual-series Hessian and spectral methods could share the window between d1
  and d2 requests for the same point.
