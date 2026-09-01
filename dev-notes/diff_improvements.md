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

## Follow-ups (not done)

- Spectral Hessian off-diagonals (needs 2D fits).
- The spectral round-off model over-bounds by ~1e3 (sum of amplified
  magnitudes; the true DCT noise partially cancels). A quadrature-sum model
  would tighten est/act but needs a honesty re-audit.
- numdifftools comparison not re-run (not installed in this env).
- Dual-series Hessian and spectral methods could share the window between d1
  and d2 requests for the same point.
