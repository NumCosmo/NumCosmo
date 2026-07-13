# UltraLevin ↔ Xcor: mechanics, correctness gap, and a batched interface plan

Status: draft plan (2026-07-13). Written before any code changes, per decision
below (§5, item 1). Nothing in this document has been implemented yet.

## 1. What "UltraLevin" is

Not a class name in the codebase — the term refers to the combination of two
objects introduced together in "New xcor" (#240):

- [`NcmSBesselOdeSolver`](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.c) —
  solves the modified spherical Bessel ODE
  `x²y'' + 2xy' + (x² - ℓ(ℓ+1))y = f(x)` on an interval `[a,b]` by ultraspherical
  (Chebyshev→Gegenbauer, `NcmSpectral`) spectral collocation, producing a banded
  operator that can be reused for multiple right-hand sides.
- [`NcmSBesselIntegratorLevin`](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c) —
  a Levin-type collocation integrator for `∫ K(x,k) jℓ(kx) dx`, built entirely on
  top of that ODE solver (`sbilv->ode_solver` / `ode_operator`,
  [ncm_sbessel_integrator_levin.c:75-76](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L75-L76)).
  It implements the generic `NcmSBesselIntegrator` interface, alongside two
  older backends (`GL`, `FFTL`) used for benchmarking/cross-checks.

### 1.1 Why it is fast: three independent axes of reuse

The object's internal state is deliberately structured so that the expensive
part (panel operator factorization + spherical-Bessel-at-knots cache) is
computed once and reused along three separate axes:

1. **Ell-block axis.** `ncm_sbessel_integrator_set_ell_range(sbi, ell_min,
   ell_max)` is the only thing that triggers a rebuild — and only if the range
   actually changed:
   [ncm_sbessel_integrator_levin.c:850](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L850)
   (`if ((ell_min != sbilv->alloc_ell_min) || (ell_max != sbilv->alloc_ell_max))`).
   A rebuild recomputes the `jl_knots` cache
   ([`_ncm_sbessel_integrator_levin_prepare_ell_cache`](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L116))
   and the per-panel operators
   ([`_ncm_sbessel_integrator_levin_prepare_knots_operators`](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L118))
   for *every* ℓ in `[ell_min, ell_max]` at once — "batched mode is
   automatically selected based on `n_ℓ`" (header docstring,
   [ncm_sbessel_ode_solver.h:69-70](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.h#L69-L70)).
   Calling `set_ell_range` with the *same* range repeatedly is free.
2. **Knot/panel axis.** Operators are built once per log-spaced `y = kx` knot
   panel ([`knots`](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L95),
   bounds `y_knots_min`/`y_knots_max`,
   [ncm_sbessel_integrator_levin.h:63-64](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.h#L63-L64))
   and stored in `operators`
   ([ncm_sbessel_integrator_levin.c:97](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L97)).
   Because the panels live in the *dimensionless* `y = kx` variable, a single
   panel set serves **any** physical `k` and **any** physical `[a,b]` interval
   that map into it — this is what lets one integrator instance answer many
   `k` queries during an adaptive k-sampling without re-factorizing anything.
3. **RHS axis.** For a fixed panel + ell-block, only the right-hand side
   (`_compute_rhs`, driven by the caller's `K(x,k)` at the new `k`) changes
   between calls; the factorized operator is solved against it via
   `ncm_sbessel_ode_operator_solve*`
   ([ncm_sbessel_ode_solver.h:108-109](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.h#L108-L109)).

### 1.2 The one rule that falls out of this

**Configure the ell range once per block, then hammer the same integrator
instance with every k-query and every kernel that needs that block, before
touching a different ell range.** Creating a fresh integrator per kernel, or
calling `set_ell_range` once per single ℓ, or once per kernel pair, throws away
all three reuse axes and pays the full panel-factorization + jℓ-cache cost
every time.

## 2. How Xcor currently uses it — and where reuse is being thrown away

`NcXcorKernel` owns one `NcmSBesselIntegrator*`
([nc_xcor_kernel.c:73](numcosmo/nc/xcor/nc_xcor_kernel.c#L73), settable
per-kernel via the `integrator` construct property). Non-Limber evaluation
goes through `nc_xcor_kernel_get_eval_vectorized(xclk, cosmo, lmin, lmax)`
([nc_xcor_kernel.c:1090-1123](numcosmo/nc/xcor/nc_xcor_kernel.c#L1090-L1123)),
which:

- calls `ncm_sbessel_integrator_set_ell_range(self->sbi, lmin, lmax)` — fine,
  cheap if unchanged;
- then unconditionally rebuilds the **entire** k-space spline for that ℓ-block
  via `_nc_xcor_kernel_build_spline_integrand`
  ([nc_xcor_kernel.c:995-1065](numcosmo/nc/xcor/nc_xcor_kernel.c#L995-L1065)) —
  an adaptive `NcmFunctionSampleSet` construction that queries the (expensive)
  Levin integrator at many `k` points. **There is no cache keyed on
  `(kernel, ℓ-block)`.**

`NcXcor`'s kernel-space methods (`KERNEL_GSL`/`KERNEL_CUBATURE`,
[nc_xcor.c:838-1059](numcosmo/nc/xcor/nc_xcor.c#L838-L1059)) call
`nc_xcor_kernel_get_eval_vectorized` once per kernel *per pair* — see the loop
building `xcor_arg->xclki1`/`xclki2` inside `_nc_xcor_kernel_cubature`
([nc_xcor.c:916, 943, 998-999, 1030-1031](numcosmo/nc/xcor/nc_xcor.c#L916)).
`NcDataXcor` then drives this pair-by-pair over the full tomographic set:

```c
/* nc_data_xcor.c:530-556 */
for (a) {
  nc_xcor_compute (dxc->xc, xcl1, NULL, cosmo, 0, ..., cl_th_0_aa);   /* auto */
  for (b) {
    nc_xcor_compute (dxc->xc, xcl1, xcl2, cosmo, 0, ..., cl_th_0_ab); /* cross */
  }
}
```

**Consequence:** for `N` kernels and the `N(N+1)/2` auto+cross pairs a full
tomographic analysis needs, each kernel's expensive k-space spline for a given
ℓ-block gets rebuilt once *per pair it participates in* (up to `N-1` times
more than necessary), instead of once per `(kernel, ℓ-block)` and shared across
every pair that touches it in that block. This is the gap you flagged: the
current interface has no place to say "here are all the kernels and all the
Cℓ's I want — figure out the cheapest order," because `nc_xcor_compute()` is a
single-pair, single-block call with no memory of sibling calls.

## 3. Known correctness bug (independent of the batching issue, must land first)

Numeric smoke test (CMB-lensing auto-Cℓ, ℓ=2–50, comparing
`LIMBER_Z_CUBATURE` vs `KERNEL_GSL` vs `KERNEL_CUBATURE`): the two kernel-space
methods are wrong by **15-27 orders of magnitude** and disagree with *each
other* by ℓ-dependent factors (ruling out a missing global constant).

Root cause (traced, not fixed): `nc_xcor_kernel_get_eval_vectorized` builds a
`NcmSplineCubicNotaknot` over a domain set by adaptive epsilon-expansion
(`sid->k_min/k_max = ncm_function_sample_set_get_x_min/max`,
[nc_xcor_kernel.c:1048-1049](numcosmo/nc/xcor/nc_xcor_kernel.c#L1048-L1049)).
But `_nc_xcor_kernel_gsl`/`_nc_xcor_kernel_cubature` pick their *outer*
k-integration bounds from `nc_xcor_kernel_get_k_range()`
([nc_xcor_kernel.c:1388-1429](numcosmo/nc/xcor/nc_xcor_kernel.c#L1388-L1429)),
a completely independent calculation (`ν/ξ_max ≤ k ≤ ν/ξ_min` Limber band ∩
hard component limits) with no guarantee of matching the fitted spline domain.
`_spline_integrand_eval` never clamps
([nc_xcor_kernel.c:441-453](numcosmo/nc/xcor/nc_xcor_kernel.c#L441-L453)), so
GSL/cubature almost certainly evaluate the not-a-knot cubic spline far outside
its fitted range, where cubic extrapolation diverges — consistent with the
observed magnitudes and the GSL-vs-cubature disagreement (different adaptive
sampling ⇒ different extrapolation garbage).

Test coverage never caught this: `test_xcor_compute_methods`
([test_integration.py:93-133](tests/python/nc/xcor/test_integration.py#L93-L133))
claims (in its docstring) to check `KERNEL_GSL` but never calls it;
`KERNEL_CUBATURE` appears in no test file.

**This must be fixed before the batching redesign is validated** — there is no
point optimizing the reuse pattern of a computation that returns garbage.

## 4. Proposed two-phase Xcor interface

Goal: let the caller declare *everything it wants* up front, so Xcor can pick
an ℓ-block tiling and, for each block, build every distinct kernel's spline
exactly once and reuse it across every pair that needs it, before moving to
the next block.

### 4.1 Phase 1 — registration (configuration only, no computation)

```c
guint nc_xcor_register_kernel   (NcXcor *xc, NcXcorKernel *xclk);
void  nc_xcor_request_cl        (NcXcor *xc, guint kernel_id_1, guint kernel_id_2,
                                  guint lmin, guint lmax);
```

- `register_kernel` is idempotent by pointer identity (registering the same
  kernel twice returns the same id).
- `request_cl` declares one desired output block (a pair + ℓ-range); call it
  once per entry in the data vector / covariance block the caller ultimately
  needs (this mirrors what `NcDataXcor` already knows per `xcab[a][b]`,
  including `ell_th_cut_off`).

### 4.2 Phase 2 — solve (all computation happens here)

```c
void nc_xcor_solve (NcXcor *xc, NcHICosmo *cosmo);
NcmVector *nc_xcor_get_cl (NcXcor *xc, guint kernel_id_1, guint kernel_id_2);
```

`nc_xcor_solve` internally:

1. Computes the union of ℓ-blocks needed across all registered requests
   (tiling `[ℓmin,ℓmax]` at a configurable batch size, consistent with the
   existing `ell_batch_size`/`MAX_ELL_BLOCK`/`ell_cache_max` constants).
2. For each ℓ-block, in order:
   a. For each *distinct* kernel referenced by a request touching this block,
      build its k-space spline once
      (`nc_xcor_kernel_get_eval_vectorized`), and cache it for the block's
      lifetime only.
   b. For each request whose ℓ-range intersects this block, compute the
      k-integral of the (cached) kernel-product and write the corresponding
      slice of that pair's output vector.
   c. Drop the block's cached splines before moving to the next block.
3. Results become retrievable via `nc_xcor_get_cl`.

This turns the cost from `O(N_pairs × N_blocks)` kernel-spline builds into
`O(N_kernels × N_blocks)` — the actual win you're after.

`nc_xcor_compute()` stays as today's convenience wrapper (register both
kernels + one request + solve, for the single-pair case), documented as the
non-optimal path once more than a couple of kernels are involved.

## 5. Decisions needed before implementation starts

| # | Question | Notes |
|---|----------|-------|
| 1 | Fix the k-range bug (§3) first, in isolation, or alongside the batching redesign? | Recommend: fix first, with its own regression tests, so the batched interface is validated against known-correct numbers. |
| 2 | Where does the new stateful API live — on `NcXcor` itself, or a new companion object (e.g. `NcXcorSolver`)? | `NcXcor` today is close to stateless per-call; bolting registration/request state onto it changes its contract. A companion object keeps `nc_xcor_compute()`'s existing callers untouched. |
| 3 | How to tile ℓ-blocks when different pairs need different ℓmin/ℓmax (as `nc_data_xcor.c`'s per-pair `ell_th_cut_off` already implies)? | Union-then-mask (compute the full block, use only the requested sub-range per pair) vs. per-pair block grids. Union-then-mask is simpler but may compute unused ℓ's for some pairs. |
| 4 | Per-block spline cache lifetime beyond "duration of one `nc_xcor_solve` call"? | Splines depend on the current cosmology, which changes every likelihood/MCMC step, so cross-step caching is not applicable — only cross-pair-within-one-`solve()` caching matters. Confirm no case needs longer-lived caching. |
| 5 | Threading model for block processing? | Building the `N_kernels` splines within a block are independent of each other and could run in parallel; need to confirm `NcmSBesselIntegratorLevin`/`NcmSBesselOdeSolver` instances are one-per-thread (they hold mutable working arrays, so likely each kernel needs its own integrator instance already — check whether that's already true or whether kernels currently share one integrator). |
| 6 | Migration of `NcDataXcor` | The call site at [nc_data_xcor.c:530-556](numcosmo/nc/data/nc_data_xcor.c#L530-L556) is the one production consumer that benefits; plan its switch to register/request/solve explicitly rather than leaving it on the single-pair path. |

## 6. Suggested order of work (once the above is confirmed)

1. Fix the k-range/domain mismatch (§3); add `KERNEL_GSL`/`KERNEL_CUBATURE`
   regression tests (agreement with each other, convergence to Limber at high
   ℓ, re-run against the CCL comparison currently Limber-only).
2. Land the companion object / new API (§4) purely as an additive layer,
   `nc_xcor_compute()` unchanged and re-expressed in terms of it.
3. Migrate `NcDataXcor` to register/request/solve.
4. Expose the non-Limber knobs (`integrator`, `l-limber`, `reltol`,
   `adaptive-*`) in `numcosmo_py/app/xcor` production config, not just the
   `view.py` debug tool.
5. Performance pass + docs (theory page for the ODE-solver/Levin combination,
   `docs/theory/`, once behavior is verified — no page exists today).
