# Exact outer quadrature for NcXcor

Status: `NC_XCOR_METHOD_KERNEL_EXACT` computes the outer `k` integral exactly,
by GL(5) over the **common refinement of each pair's own knot sets**. An earlier
implementation instead sampled all kernels jointly onto one shared abscissa;
that was slower and less accurate and has been removed. Sections 1-7 record why,
so the joint design is not re-derived. Sections 8-10 are unchanged and still
apply.

## 1. The outer k integral is exactly integrable

Each component of a spline-backed `NcXcorKernelIntegrand` is a **cubic
not-a-knot spline in `k`** — `ncm_spline_vec_eval (sid->spline_vec, k, ...)`,
no log transform. The outer integral is

    C_l ~ int d(ln k) k^3 W_i W_j  =  int dk k^2 W_i W_j

so on any interval where both splines are a single cubic piece, the integrand is
a polynomial of degree `2 + 3 + 3 = 8`. Gauss-Legendre with `n` nodes is exact
through degree `2n - 1`, hence:

> **GL(5) integrates the outer k integral exactly, on any panel where both
> members of the pair are a single cubic.**

It must be done **in `k`, not `ln k`**: under `u = ln k` the integrand is a cubic
in `e^u`, not in `u`, and GL is no longer exact. This also explains why a linear
`k` grid beat a log `k` grid at l=20 by 10-25x in an earlier trapezoid study.

## 2. Those panels are the common refinement — merge, do not co-sample

Two kernels sampled independently do not share an abscissa, so the panels with
the property §1 needs are the panels of the **union of the two knot vectors**.
Both are sorted, so this is a linear merge, done per pair inside
`_nc_xcor_kernel_integrate_block_exact()`.

**An earlier version of this note argued the union was useless**, on the
measurement that `union / sum-of-separate = 1.000` — independently adapted
splines resolve the same oscillation at nearby but distinct abscissas, so the
union is as large as the sum. That measurement is correct but it is the wrong
comparison. The union is not competing with "one kernel's knots"; it is
competing with a **shared** set, and a shared set is a max-density set: it
forces every pair to carry all `N` kernels' knots, not just its own two.
Measured over 7 bins at `l = 0`: 28 per-pair unions total 69356 panels against
93072 for 28 sweeps of the shared set.

## 3. The sampling dominates, not the quadrature

This is the fact that decides the design. Measured, 7 bins, 28 pairs,
`scaled_abstol = 1e-6`:

| step | cost at l=0 |
|------|-------------|
| sampling, 7 separate builds | 0.085 s |
| sampling, 1 joint build     | 0.247 s (**2.91x**) |
| outer quadrature, either    | ~0.02 s |

Joint sampling costs ~1.85-2.9x more Levin evaluations to build the shared
abscissa, and buys a cheaper version of a step that is under 10% of the runtime.
Optimizing the quadrature at the sampling's expense is the wrong trade, and that
is what the joint design did.

## 4. Result

`KERNEL_EXACT` reuses the per-kernel closures `KERNEL_CUBATURE` already builds —
and which `NcXcorSolver` already caches once per kernel per ell block — and
differs from it *only* in the outer-integration call. Measured over 28 pairs of
7 bins:

| | l = 0 | l = 0..26 | agreement vs CUBATURE |
|---|-------|-----------|------------------------|
| `KERNEL_EXACT`, joint abscissa (removed) | 0.271 s | 4.53 s | 7.5e-5 |
| `KERNEL_EXACT`, per-pair union (current) | **0.099 s** | **1.96 s** | **2.4e-8** |
| `KERNEL_CUBATURE` | 0.104 s | 2.32 s | — |

So the exact method is now slightly *faster* than the adaptive one, and agrees
with it three and a half orders of magnitude better.

## 5. Why the accuracy improved too

The joint sampler's shared domain is the union of the kernels' domains, so a
kernel was splined past its own adaptive cutoff, where the spline carries
extrapolation rather than signal. For a mid-z bin the true `W` at `k~ = 2703` is
`2.05e-5` while the spline gives `14.9`. Weighted by `k^2` that junk dominated
the auto spectrum (`S(1,1) = 4.46e-4` against the correct `5.36e-5`).

That forced a per-component support mask (`|y| >= epsilon * max|y|`), which then
*clipped* the opposite case: for well-separated bins at high multipole the cross
spectrum lives on the overlap edge — tail times tail — and a per-component
criterion cannot know that a *pair* needs the tail. Measured: bins z=0.1-0.2 and
z=0.6-0.7 at l=9 gave `-3.8e-11` against cubature's `-7.6e-11`.

Both problems are artifacts of the shared domain. With per-kernel knots each
spline is only ever evaluated inside its own fitted domain, the pair range is
the intersection of the two `get_range()`s — exactly what `KERNEL_CUBATURE`
uses — and neither the extrapolation blow-up nor the tail clipping arises. The
support machinery was deleted with the joint sampler.

## 6. Share across kernels? No. Across ell blocks? Also no.

Knot count grows steeply with `l` (3023 at l=0 to 17755 at l=100, tracking
`k_max`), so a set shared across ell blocks makes every block pay the worst
block's count: ~160k evaluations for 9 blocks versus 66k for separate sets,
~2.4x worse. Sharing across ell blocks was never done. Sharing across kernels
within a block was tried, and §3 is why it was undone.

## 7. What KERNEL_EXACT is for

- **no tolerance parameter**, and no adaptive outer step that can fail;
- **deterministic cost**, set by the knot count rather than by how hard the
  integrand turns out to be;
- it collapses the `reltol` / `scaled_abstol` interaction to a single knob,
  since the quadrature contributes no error at all;
- and it is not slower, so there is no accuracy-for-speed trade to weigh.

Note that the `pcubature` abort that originally motivated an exact method is
*not* an intrinsic property of top-hat kernels — it is induced by
`scaled_abstol`, measured at `l = 32-39` over 7 bins:

| `scaled_abstol` | knots/kernel | `KERNEL_CUBATURE` |
|-----------------|--------------|-------------------|
| `1e-4` (default used by `test_kernel.py`) | 493 | OK, 0.57 s |
| `1e-6` | 2367 | OK, 2.63 s |
| `1e-7` | 4573 | OK, 5.01 s |
| `1e-8` | ~8759 | **abort** |

Cubature is fine on top-hats until the spline is refined enough to resolve the
ringing tail of §9, which makes the outer integrand genuinely oscillatory. The
exact rule removes that failure mode by construction, which is what makes
`scaled_abstol = 1e-8` usable.

## 8. Where the sampled domain is actually decided

`ncm_function_sample_set_expand_domain()` owns it, not the xcor kernel code:

    norm_y = ncm_vector_dnrm2 (y);
    F_ref  = ncm_function_sample_set_get_absmaxF_l2_norm (fss);
    if (norm_y < epsilon * F_ref)   /* pointwise amplitude vs running max */

The per-component `left_boundary_found` / `right_boundary_found` tracking in
`_component_states_compute_non_limber()` looks like it decides the range but
does not: it only controls when a component switches to exponential tail
*extrapolation*. An attempt to make the truncation cumulative by editing that
tracking left `k_max` identical to four digits, while shifting when
extrapolation engages -- which *increased* the knot count (3023 -> 3613) and the
runtime (0.36 s -> 0.63 s). Reverted.

So the tail problem of §9 is fixable only in `expand_domain`'s convergence test,
which is generic `NcmFunctionSampleSet` code shared with other callers. That
wants either an opt-in mode or an audit of every caller, not an in-place change.

## 9. The top-hat tail: why these integrands are hard

The other xcor kernels (`gal`, `wl`, `cmb_lens`, `tsz`) are smooth in `z`, so
their `U(k)` decays fast: amplitude and integrand fall off together, the
pointwise criterion of §8 cuts in the right place, and `KERNEL_CUBATURE` is
happy at `l = 800` after years of use. A top-hat bin is the pathological case.
Its radial integral is over a compact interval, so `W(k)` *rings* -- decaying
as a slow power law -- while the integrand `k^3 W^2` decays fast.

Measured on bin `z = 0.1-0.2` at `l = 0` (2327 knots over `k~ = [0.0044, 3248]`):

| fraction of the integral | reached by | knots used |
|--------------------------|------------|------------|
| 90%                      | `k~ = 54`  | 5.3%       |
| 99%                      | `k~ = 191` | 12.5%      |
| 99.9%                    | `k~ = 489` | 28%        |
| 99.99%                   | `k~ = 1037`| 49%        |

**71% of the knots lie above `k~ = 500`, where under 0.1% of the integral is.**
Local slopes: `|W| ~ k^-2.0` to `k^-3.1`, `k^3 W^2 ~ k^-1.1` to `k^-3.2`.

Candidate truncation criteria, and what each yields:

| criterion                          | `k_max` |
|------------------------------------|---------|
| `\|W\|` (what §8 uses)             | 1668    |
| `k^{3/2} \|W\|`                    | 3248 — *wider* |
| `k^3 W^2` (the integrand itself)   | 3236 — *wider* |
| remaining tail integral < `1e-3`   | **489** |
| remaining tail integral < `1e-4`   | **1037**|

The lesson: *any* pointwise "amplitude below `epsilon` times peak" test
truncates late on a power-law tail, at roughly `(1/epsilon)^(1/p)` times the
peak location, whatever weight it carries. Only a cumulative test -- comparing
the tail still to come against what has accumulated -- tracks the integral.
Worth 2-3.5x in knots, hence in runtime.

It would also remove the `pcubature` abort of §7 at its root, since that abort
is what happens when this tail is resolved well enough to be oscillatory: the
effort spent resolving it is what breaks the outer integration. Nothing here is
urgent, though -- at the tolerances actually in use (`scaled_abstol >= 1e-6`)
the current criterion works and top-hats validate at `l = 800`.

## 10. Related open items

- `ncm_integral_nd_eval` hard-`g_error`s when `pcubature` runs out of
  Clenshaw-Curtis levels. This blocks `scaled_abstol = 1e-8` on masked runs and
  would kill an MCMC chain mid-flight. The exact GL(5) path removes it by
  construction.
- `_nc_xcor_kernel_gsl` targets `reltol * 1e-2` and aborts on GSL roundoff for
  strongly cancelling cross integrands; it converges fine at `reltol = 1e-4`.
