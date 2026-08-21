# The ISW Limber step: why GH #297 aborts, and the shape of the fix

Branch `xcor-isw-aborts`, off `d2fd03b3`. Everything below is measured, not
estimated. Visual write-up of the same material:
<https://claude.ai/code/artifact/6d88e2fa-5b01-411f-a0de-bb5cf1f0d2af>

## 1. The math

The outer integral, as the code computes it (`nc_xcor.c`, variable `ln k`):

```
C_l = 2 / (pi RH^3) * INT d ln k  k^3 W_l^(1)(k) W_l^(2)(k)
```

Exactly, the kernel is

```
W_l(k) = INT_{xi_min}^{xi_max} d xi  K(xi, k) j_l(k xi)
```

which is **smooth in k for every k**: k enters through a smooth `K` and an
entire `j_l`, over a *fixed* interval. The true `C_l` has no feature at any
special k. Every discontinuity discussed here is created by the approximation.

Limber collapses that integral to a single point at `xi = nu / k`, `nu = l+1/2`:

```
W_l^Limber(k) = sqrt (pi / (2 nu)) * (1/k) * K (nu/k, k)
```

The kernel has compact support, so `K = 0` for `xi > xi_max`, and under Limber
that becomes a condition on k:

```
nu/k > xi_max   <=>   k < k_min_l = nu / xi_max      (likewise k_max_l = nu / xi_min)
```

So `W_l` is identically zero below `k_min_l` and jumps there by an amount set by
one number: **the kernel's value at the edge of its own support**. And
`k_min_l` is *linear in l*, so a multipole batch produces a comb of steps, ~5%
apart at l ~ 20.

## 2. Why ISW and not weak lensing

`K(xi, k=10)` approaching `xi_max`:

| 1 - xi/xi_max | weak lensing | ISW |
|---|---|---|
| 1e-02 | 2.87e-17 | 5.12e+00 |
| 1e-04 | 1.55e-21 | 1.24e+01 |
| 1e-06 | 1.54e-25 | 1.25e+01 |
| 0 | **0** | **1.25e+01** |

The lensing efficiency vanishes at its support edge, so its Limber integrand is
continuous and batching never creates a step. The ISW kernel is still at its
plateau when the integration is truncated at recombination
(`xi_max ~ 3.18 RH ~ 13600 Mpc`) -- a modelling truncation that becomes a cliff.

`W_l(k)` for the batch l = 20..23, showing the comb:

| k/k_min - 1 | l=20 | l=21 | l=22 | l=23 |
|---|---|---|---|---|
| 0 | 2.64e+00 | 0 | 0 | 0 |
| 3e-02 | 3.18e-01 | 7.1e-223 | 0 | 0 |
| 1e-01 | 3.49e-02 | **1.25e-01** | **1.57e+00** | 1.3e-244 |
| 2e-01 | 6.93e-03 | 1.31e-02 | 3.04e-02 | **1.03e-01** |

## 3. Cross-correlation is immune, and why

Cross-spectra integrate the *intersection* of the two k-ranges
(`k_min = MAX(k1,k2)`, `k_max = MIN(k1,k2)`). A low-z tracer has a far larger
`k_min`, which cuts the ISW comb off entirely:

| spectrum | k range | interior breakpoints | result |
|---|---|---|---|
| ISW x ISW | [6.45, 402] | 3 | **aborts** |
| ISW x WL | [24.18, 402] | none (all below 24.18) | works |

So ISW is only in trouble in **auto-spectra**, or in any pairing whose
intersection dips below ~`nu/xi_max`.

## 4. Why the code extends beyond the natural interval at all

Structural, not physical. `_nc_xcor_kernel_build_spline_integrand` fits **one**
`NcmSplineVec` for the whole multipole batch over **one** k-grid. A vector
spline needs a value for every component at every knot, but each l has a
different natural interval, so the shared grid spans the union and every
not-yet-active component needs something filled in below its own start. That
"something" is the `DECAY_RATE` extrapolation in
`_component_states_compute_limber`.

`DECAY_RATE = 1e10` makes `exp (-(1e10 dk/k)^2)` fall from 1 to underflow across
`dk/k ~ 1e-9`. It is not a taper, it is a step at machine resolution -- and on
the way down it passes through the subnormal range, emitting sign-flipping
garbage that the spline then interpolates.

Verified with a marker test (temporarily returning 0.5 from the decay): the raw
Limber value is exactly 0 where `_spline_integrand_eval` returns
+/-1e-226 ... 1e-206. The noise is real and is what the quadrature is asked to
fit at reltol 1e-6 with `abstol = 0`.

## 5. What was tried, and why each failed

| attempt | result |
|---|---|
| `CUBATURE_H_V` instead of `_P_V` | fixes #297; WL agrees to ~1e-7; but h-adaptive must *discover* a 220-order step by bisection, and its error estimate straddling a step is not trustworthy. ISW cost 64x WL. |
| Split at analytic breakpoints (B) | works as designed: WL **bit-identical** and zero added cost (the collector rejects edges where the kernel vanishes); ISW splits at exactly the right k; two methods agree to ~1e-8. **But** it manufactures pieces where some multipoles contribute nothing, and `abstol = 0` makes those unconvergeable. |
| abstol = 1e-16 * `k_max^3 max|W|^2 * width` | far too generous (|W| peaks at small k, k^3 at large k). Fixed #297, broke 2 tests. |
| abstol = 1e-16 * sampled max integrand * width | correctly tight, the 2 tests recover, but no longer rescues #297. **Not diagnosed -- open thread.** |

**Correction to an earlier claim in this effort:** the abstol is *not* orthogonal
to the split. Splitting manufactures the zero-component problem; B depends on
the tolerance rather than composing with it.

## 6. The shape the fix should take (not built)

Zeroing components exactly outside their own interval is necessary but not
sufficient: it turns the comb into genuine hard steps, and splitting at them
puts the jump on a *shared endpoint*, where the two pieces need different
one-sided values. Clenshaw-Curtis weights endpoints, so whichever piece gets the
wrong one picks up a spurious contribution that moves under refinement.

The clean resolution is to **stop sharing an interval across multipoles**:
integrate each l over its own `[nu/xi_max, nu/xi_min]`, intersected with the
fitted spline domain (and the second kernel's interval for a cross-spectrum).
The discontinuity is then exactly *at the integration limit*, where it is not a
discontinuity of the integrand at all -- the function is smooth over the whole
interval and the endpoint is the correct one-sided limit.

Consequences: `abstol = 0` stays, `pcubature` stays, no breakpoint machinery, no
masking, no floor to size. `DECAY_RATE` and `last_values_left/right` become dead
code.

The batching that this gives up is batching of the part that costs nothing --
§9.6 of the batching plan measured the outer k-integral at **~0%, noise-level**;
all cost is the closure build, which stays once per block.

Three things to check before building it:

1. `k_max_l = nu/xi_min` with `xi_min ~ 1e-6` gives `k_max_l ~ 2e7`, far outside
   the fitted spline domain, so it clips to the spline's `k_max`. Confirm the
   clip does not reintroduce an edge.
2. Cost: §9.6 measured 23 l's in *one* call; 23 calls carry per-call setup.
   Measure it -- the same mistake made twice already in this effort was
   asserting a cost without measuring.
3. WL must stay bit-identical (its per-l intervals are wider than the fitted
   domain at the bottom, so it should clip to today's range -- a prediction to
   verify, not a claim).

## 7. State

- `0be5fdea` on `xcor-isw-aborts`: ISW type registration + the two umbrella
  headers. Green, 906 xcor tests.
- WIP for the split/decay/abstol attempts: branch `xcor-isw-split-wip`. Does not
  work; kept only as a record of the measurements above.
- GH #297, #298 remain open. **#298 has not been touched** -- it reports
  `roundoff error` from `KERNEL_GSL`, a different path, and may or may not share
  this root cause. Measure before assuming.

## 8. Resolution (branch `xcor-isw-integration-fixes`)

Built as §6 prescribes, and its three checks came out as follows.

**#297.** `NcXcorKernelIntegrand` now reports a per-component range
(`nc_xcor_kernel_integrand_get_range_comp()`): for a Limber closure, the
multipole's own band taken over the kernel's components and clipped to the
fitted domain; for a non-Limber one, the fitted domain itself. The cubature
driver integrates *runs* of consecutive multipoles that share a range, so the
non-Limber case is still one call over the whole block and only Limber blocks
split. A run evaluates only its own components -- one index lookup on the
shared abscissa, then `ncm_spline_eval_idx` per component -- which is what
keeps the split from costing anything.

Measured on the issue's own reproducer (ISW auto, tier 2, l = 20..27,
`ell_batch_size = 4`), for l = 22:

| lower limit | value | reported error | true error | evaluations |
| --- | --- | --- | --- | --- |
| shared | 7.94006 | 2.1e-07 | 1.4e-02 | 1048952 |
| its own | 8.05126 | 2.4e-07 | 3.9e-08 | 257 |

The old failure was never only the abort: after PR #307 gave
`ncm_integral_nd_eval()` an h-adaptive retry, `hcubature` converged *by its own
estimator* while missing part of the spike -- a region whose nodes all land
where the component vanishes reports no error and is never subdivided.

Check 1 (the `k_max_l = nu/xi_min` clip): clips to the spline's `k_max`, which
is where every non-Limber component already sits, so it introduces no edge.
Check 2 (cost), l = 2..1025: ISW tier 2 at block 16 went 8.056 s -> 3.135 s,
i.e. from 2.6x `KERNEL_FIXED` to level with it; CMB lensing at block 64 went
0.227 s -> 0.235 s. Check 3 (WL bit-identical) was wrong as stated: a block's
fitted domain starts at its *lowest* multipole's band edge, so every higher
multipole in the block has a band edge of its own and splits too. CMB lensing
moves in its last two digits (7.7917110410e-08 -> 7.7917110417e-08), both
converged.

**#298 does not share the cause.** `KERNEL_GSL` builds a closure per multipole,
so it never sees the comb. Its integrand is a cubic spline over ~380 knots:
C^2, with a jump in the third derivative at every knot. A Gauss-Kronrod rule
spanning several knots applies a smooth-function error model to a function that
is not, and its estimate stops falling under bisection -- which is exactly what
QUADPACK reports as roundoff. Measured for ISW tier 3, l = 20, at the failing
`reltol * 1e-2` request:

- 159 subintervals, each claiming 4-5e-9 of error whether it contributes 7.5e-2
  or 1.6e-5 to the integral -- an estimate uncorrelated with the contribution.
- Total claimed 1.8e-7 against a true error of 3.2e-9.
- One 4-knot subinterval claims 75-88x more error whole than split on its knots.
- The whole range, panel by panel: claimed 2.3e-14, exact to every digit.

So `_nc_xcor_kernel_gsl()` now integrates with `gsl_integration_qagp()` on the
merged knots of the pair when both integrands report knots, leaving every
subinterval a single cubic piece: machine precision (3e-15 against the exact
quadrature) instead of a stall near 1e-6, for about twice the evaluations of an
outer integral that costs nothing next to the per-multipole closure build. The
`reltol * 1e-2` margin is gone from that path -- measured, it only moved the
stall -- and a non-success status is fatal only when the achieved error misses
NcXcor:reltol.

**A third defect, newly visible.** With both methods working, the tier-2
CMB lensing x ISW cross sits 1.05e-3 below tier 1, flat under refinement of the
spline tolerance and collapsing to 1.5e-6 when `NcXcorKernel:adaptive-epsilon`
is tightened to 1e-8. The closure is fitted only out to where *its own*
amplitude dies; an auto spectrum never pays for that (its integrand carries W
squared) but a cross pays once. No per-kernel criterion can fix it, since a
closure is built without knowing what it will be paired with; a criterion
relative to the kernel's own peak would at least bound the loss by
eps * max|W_1| * max|W_2|, i.e. relative to sqrt(C_1 C_2) -- the amplitude that
sets the cross spectrum's sample variance -- rather than relative to C_x. Here
the two are 7% correlated, so the 1.05e-3 is 7e-5 of that scale. Recorded in
`test_xcor_kernel_methods`, not fixed.
