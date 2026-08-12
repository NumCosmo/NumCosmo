# Shared knot sets and exact outer quadrature for NcXcor

Status: knot accessor, joint sampler and per-component support masking all
landed and tested. Written up so the measurements below are not repeated and the
retracted claim in §5 is not re-derived. Remaining work is in §7.

## 1. The outer k integral is exactly integrable

Each component of a spline-backed `NcXcorKernelIntegrand` is a **cubic
not-a-knot spline in `k`** — `ncm_spline_vec_eval (sid->spline_vec, k, ...)`,
no log transform. The outer integral is

    C_l ~ int d(ln k) k^3 W_i W_j  =  int dk k^2 W_i W_j

so on every knot panel the integrand is a polynomial of degree `2 + 3 + 3 = 8`.
Gauss-Legendre with `n` nodes is exact through degree `2n - 1`, hence:

> **GL(5) per spline panel integrates the outer k integral exactly.**

Measured, refining every panel 8x:

| l  | pair  | panels | evals  | rel. err |
|----|-------|--------|--------|----------|
| 0  | (0,0) |   3022 |  15110 | 2.4e-15  |
| 0  | (0,1) |   6027 |  30135 | 1.8e-15  |
| 0  | (0,6) |   9136 |  45680 | 4.3e-14  |
| 20 | (0,6) |  18177 |  90885 | 1.6e-10  |

It must be done **in `k`, not `ln k`**: under `u = ln k` the integrand is a cubic
in `e^u`, not in `u`, and GL is no longer exact. This also explains why a linear
`k` grid beat a log `k` grid at l=20 by 10-25x in an earlier trapezoid study.

Consequence: the outer quadrature never needs to be more accurate than
`scaled_abstol * peak`, so its node count follows from the inner tolerance with
nothing left to tune — and no adaptive outer rule that can fail to converge
(see the `pcubature` abort, §6).

## 2. Knot sets must be built jointly, never unioned

Independently adapted splines resolve the same oscillation at *nearby but
distinct* abscissas, so unioning them buys nothing:

    union over 9 ells = 65920 knots, sum of separate = 65945  ->  union/sum = 1.000
    union across 7 bins at fixed ell -> 5.2x the per-bin max (5.6x at l=100)

Built **jointly** (`nc_xcor_kernel_get_eval_vectorized_joint()`), the shared set
is a max-density set, not a sum:

| l  | per-bin max | sum   | JOINT | joint/max | joint/sum |
|----|-------------|-------|-------|-----------|-----------|
| 0  |        6897 | 35851 |  9455 |     1.37x |     0.26x |
| 8  |       10137 | 53103 | 12301 |     1.21x |     0.23x |
| 32 |       16409 | 86973 | 18629 |     1.14x |     0.21x |

Cost: ~1.85x more Levin evaluations at l=0 (1.50x at l=32) than separate
builds. What it buys is exact quadrature plus every pair from one `U^T W U`
(47270 nodes, 0.05 s for all 28 pairs of 7 bins, versus ~2 s for 28 separate
exact quadratures on per-pair unions).

## 3. Share across kernels, NOT across ell blocks

The two axes behave oppositely. Knot count grows steeply with `l` (3023 at l=0
to 17755 at l=100, tracking `k_max`), so a set shared across ell blocks makes
every block pay the worst block's count: ~160k evaluations for 9 blocks versus
66k for separate sets, **~2.4x worse with no compensating benefit** (the matrix
product is per-ell regardless). Sharing across kernels within a block is the
win; sharing across blocks is not.

## 4. Per-component support masking (landed)

`NcmSpline` deliberately does **not** range-check on evaluation — it is a sharp
edge kept for hot-loop performance. The NumCosmo convention is therefore
*evaluate the domain first, then evaluate only inside it*.

The joint sampler currently violates this: its shared domain is the union of the
kernels' own domains, so a kernel is sampled and splined past its own adaptive
cutoff, where the spline carries extrapolation junk. For a mid-z bin the true
`W` at `k~ = 2703` is `2.05e-5` while the spline gives `14.9` — off by `7e5`.
Weighted by `k^2` that junk dominated the auto spectrum (`S(1,1) = 4.46e-4`
against the correct `5.36e-5`); crosses survived because the other kernel damps
the region.

Fix, implemented: keep the knot set shared, keep the **support per component**.
`nc_xcor_kernel_integrand_get_component_range()` reports each component's own
`[k_min, k_max]`, recorded after sampling as the outermost knots where
`|y| >= epsilon * max|y|` — the same criterion the per-kernel boundary search
converges to, applied to the final data so it does not depend on the order the
sampler explored in. A pair `(i,j)` is integrated over
`[max(kmin_i, kmin_j), min(kmax_i, kmax_j)]`; those bounds are knots of the
shared set, so the panels stay spline panels and the quadrature stays exact.
Separate (single-kernel) integrands carry supports too, so callers need no
special case.

With this, GL(5) on the shared panels reproduces `KERNEL_CUBATURE` for every
pair, and the whole xcor suite passes (924 tests).

## 5. Retracted: "the per-kernel cutoff measures the wrong quantity"

An earlier reading of the above argued that the adaptive cutoff is decided from
`|U_i|` against its own peak while the integrated quantity is `k^2 U_i U_j`, and
that the cutoff therefore had to become pair-driven. **That was wrong** — it was
built on spline values read outside the fitted domain. The cutoff is correctly
placed:

- `W` matches an independent recomputation to `<= 1.6e-3` across the whole
  adaptive range, and the disagreement appears only outside it.
- The integrand `k^3 W^2` peaks at `k~ = 46` and decays to `0.005` by
  `k~ = 1.3e4`; **100.00%** of the mass lies below the cutoff at `k~ = 2459`.
- GL(5) on each kernel's own panels reproduces `KERNEL_CUBATURE` and
  `KERNEL_GSL` to 7 digits on every pair tested.

Domain expansion does not need to become pair-driven. Autoknots do need joint
orchestration (done); supports do not.

## 6. Idea kept for later: a global orchestrator

Worth revisiting once the above is stable. Because the sampler walks `k` by `k`
and kernel by kernel and stitches, it *could* in principle choose the sampled
range from what the requested **pair set** actually needs, rather than from each
kernel independently. That would help when only a particular product is of
interest, and would hurt in tail-times-tail cases, where the pair criterion
demands more range than either kernel alone.

**There is now a concrete case that motivates it.** `KERNEL_CUBATURE` intersects
the two integrands' full *ranges*; `KERNEL_FIXED` intersects their per-component
*supports*, which are narrower. For well-separated bins at high multipole the
cross spectrum lives exactly on that overlap edge — a tail-times-tail case — and
the two methods diverge by a factor of ~2 on entries three orders below the
block maximum (bins z=0.1-0.2 and z=0.6-0.7 at l=9: `-7.6e-11` vs `-3.8e-11`).
The per-component support criterion (`|y| >= epsilon * max|y|`, §4) cannot know
that a *pair* needs the tail, so it clips it.

This does not affect anything measured so far — the entries involved are far
below the absolute floor that `scaled_abstol` sets anyway, and the SSC use case
never looks at them. But it is the first evidence that a pair-aware criterion
buys accuracy, not just speed, and it identifies where to look for it: high `l`,
widely separated kernels.

## 7. KERNEL_FIXED: what it is and is not worth

`NC_XCOR_METHOD_KERNEL_FIXED` samples the kernels jointly and integrates the
outer `k` integral by exact GL(5) over the knot panels, masking each component
to its own support. `NcXcorSolver` drives it with **one joint build per ell
block** covering every registered kernel, reading each requested pair out of it.

It agrees with `KERNEL_CUBATURE` to `1.4e-5` of the block maximum over all 28
pairs of 7 J-PAS bins, for `l = 0` and for `l = 0..26`.

**It is slower, not faster**: 0.49x at `l = 0`, 0.72x over `l = 0..26`. An
earlier estimate of a 4-9x *win* was wrong, and the reason is worth recording so
it is not re-derived: the `N^2 -> N` argument does not apply to the sampling,
because `KERNEL_CUBATURE` already caches one integrand per kernel per block in a
hash table, so it too samples `N` times, not `N(N+1)/2`. The only `N^2` was in
the outer integration, which is cheap (spline evaluations). Joint sampling then
costs ~1.85x more Levin evaluations (§2) with no sampling redundancy left to
remove.

What it does buy:

- **no tolerance parameter**, and no adaptive outer step that can fail;
- **deterministic cost**, set by the knot count rather than by how hard the
  integrand turns out to be;
- it collapses the `reltol` / `scaled_abstol` interaction to a single knob,
  since the quadrature contributes no error at all.

**Its niche is narrower than first claimed.** The `pcubature` abort that
motivated it is *not* an intrinsic property of top-hat kernels -- it is induced
by `scaled_abstol`, measured at `l = 32-39` over 7 bins:

| `scaled_abstol` | knots/kernel | `KERNEL_CUBATURE` |
|-----------------|--------------|-------------------|
| `1e-4` (default used by `test_kernel.py`) | 493 | OK, 0.57 s |
| `1e-6` | 2367 | OK, 2.63 s |
| `1e-7` | 4573 | OK, 5.01 s |
| `1e-8` | ~8759 | **abort** |

`kernel_cluster_tophat` is exercised at `l = 200/500/800` in `test_kernel.py`
and passes, at `scaled_abstol = 1e-4`. Cubature is fine on top-hats; it breaks
only once the spline is refined enough to actually resolve the ringing tail of
§9, which makes the outer integrand genuinely oscillatory. In other words the
abort is manufactured by resolving, at great cost, a tail worth under 0.1% of
the answer -- which is exactly the §9 problem, and further evidence that the
loud abort is a useful diagnostic rather than a fragility.

So `KERNEL_FIXED` is worth keeping as a reference implementation, and for the
narrow case of wanting the far off-diagonal entries resolved without breaking
the outer integration. It is not needed at sane tolerances.

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
