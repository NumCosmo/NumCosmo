# NcXcor non-Limber infrastructure plan (UltraLevin)

Status: plan, nothing implemented yet. Documentation-first by design — this
is the spec the implementation work will follow, kept up to date as decisions
are made and as pieces land.

## 0. Goal

`NcXcor` computes angular auto/cross power spectra `C_ℓ^{AB}` from pairs of
projection kernels. Three ways of computing it coexist in the code today, at
different stages of completion (§2). This plan brings the two unfinished ones
to production quality:

1. **Correctness** — the non-Limber kernel-space methods
   (`KERNEL_GSL`/`KERNEL_CUBATURE`) currently return numbers wrong by 15-27
   orders of magnitude. Fix the actual bug (§3).
2. **Performance** — once correct, make multi-kernel tomographic runs cheap
   by batching kernel-spline construction per ℓ-block instead of per pair
   (§4-§5), which is the entire point of UltraLevin's reuse design (§1).
3. **Production readiness** — get both of the above into `NcDataXcor` and
   `numcosmo_py/app/xcor` config, with test coverage that would have caught
   the bug in (1) (§6-§7).

The classical Limber-in-z path (`LIMBER_Z_GSL`/`LIMBER_Z_CUBATURE`) is
explicitly **out of scope for restructuring** — it stays exactly as it is,
a separate, cheaper, well-tested formula. See §2 for why, and the non-goals
in §8.

## 1. What "UltraLevin" is, and why it exists

Not a class name in the codebase — the term refers to two objects introduced
together in "New xcor" (#240), used by the non-Limber tiers:

- [`NcmSBesselOdeSolver`](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.c) —
  solves the modified spherical Bessel ODE
  `x²y'' + 2xy' + (x² - ℓ(ℓ+1))y = f(x)` on an interval `[a,b]` by
  ultraspherical (Chebyshev→Gegenbauer, `NcmSpectral`) spectral collocation,
  producing a banded operator reusable for multiple right-hand sides.
- [`NcmSBesselIntegratorLevin`](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c) —
  a Levin-type collocation integrator for `∫ K(x,k) jℓ(kx) dx`, built on top
  of that ODE solver (`sbilv->ode_solver`/`ode_operator`,
  [ncm_sbessel_integrator_levin.c:75-76](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L75-L76)).
  Implements the generic `NcmSBesselIntegrator` interface alongside two older
  benchmarking backends (`GL`, `FFTL`).

Reviewed the object itself before trusting anything built on top of it: 447
passing unit tests, including `scipy`+truth-table accuracy checks up to
ℓ=500. The panel-selection code
([ncm_sbessel_integrator_levin.c:934-1031](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L934-L1031))
correctly falls back to a direct ODE solve when `y=kx` falls outside the
precomputed knot grid rather than silently extrapolating. **The numerics
here are solid — every bug and gap in this plan is in the xcor glue layer
around it, not in `NcmSBesselIntegratorLevin` itself.**

### 1.1 Why it is fast: three independent axes of reuse

The expensive part — panel operator factorization + spherical-Bessel-at-knots
cache — is computed once and reused along three axes:

1. **Ell-block axis.** `ncm_sbessel_integrator_set_ell_range(sbi, ell_min,
   ell_max)` triggers a rebuild only if the range actually changed
   ([ncm_sbessel_integrator_levin.c:850](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L850)).
   A rebuild recomputes the `jl_knots` cache and per-panel operators for
   *every* ℓ in `[ell_min, ell_max]` at once — batched mode is automatically
   selected based on `n_ℓ`
   ([ncm_sbessel_ode_solver.h:69-70](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.h#L69-L70)).
   Calling `set_ell_range` with the *same* range repeatedly is free.
2. **Knot/panel axis.** Operators are built once per log-spaced `y = kx` knot
   panel ([ncm_sbessel_integrator_levin.c:95](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L95),
   bounds `y_knots_min`/`y_knots_max`,
   [ncm_sbessel_integrator_levin.h:63-64](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.h#L63-L64)).
   Because panels live in the *dimensionless* `y = kx` variable, one panel
   set serves any physical `k` and any physical `[a,b]` — this is what lets
   one integrator instance answer many `k` queries during adaptive
   k-sampling without re-factorizing anything.
3. **RHS axis.** For a fixed panel + ell-block, only the right-hand side
   (`_compute_rhs`, driven by the caller's `K(x,k)` at the new `k`) changes
   between calls; the factorized operator is solved against it via
   `ncm_sbessel_ode_operator_solve*`.

**The rule that falls out of this:** configure the ell range once per block,
then hammer the same integrator instance with every k-query and every kernel
that needs that block, before touching a different ell range. Creating a
fresh integrator per kernel, or calling `set_ell_range` once per single ℓ or
once per kernel pair, throws away all three reuse axes and pays the full
panel-factorization + jℓ-cache cost every time. §4-§5 design the batched
interface around this rule.

### 1.2 Hard limits and untested scale

- **`ell_cache_max` is a hard, non-catchable `g_error`**, default 1200
  ([ncm_sbessel_integrator_levin.h:66](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.h#L66)),
  enforced the moment `set_ell_range` exceeds it. This collides with
  production defaults already in the repo — `KernelCMBLensingConfig`/
  `KernelCMBISWConfig` both default `lmax=3000`
  ([kernels.py:89,122](numcosmo_py/app/xcor/kernels.py#L89)). The block
  planner (§5) **must** tile any request into `≤ get_ell_cache_max()` chunks
  as a hard constraint, not something discovered when a real config crashes
  the process.
- **Current xcor batch granularity is `MAX_ELL_BLOCK = 64`**
  ([nc_xcor_kernel.c:491](numcosmo/nc/xcor/nc_xcor_kernel.c#L491)), a fixed
  stack-array cap independent of and much smaller than the integrator's own
  `ell_cache_max` ceiling of 1200. The original version of this plan assumed
  the new interface's default block size should therefore target close to
  `ell_cache_max`, on the theory that a bigger block amortizes the panel
  setup cost over more ℓ's. **That assumption turned out to be wrong — see
  §1.3, which empirically found the opposite: 8, not 1200, is the sweet
  spot, for a specific, now-understood reason.**
- **Setup cost is kernel-independent but not shareable by construction.**
  Each `NcmSBesselIntegratorLevin` privately owns its own
  `NcmSBesselOdeSolver` ([ncm_sbessel_integrator_levin.c:134](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L134)),
  with no injection hook and no clone/dup constructor, and
  `NcXcorKernel` gives each kernel its own `integrator` property — N kernels
  pay panel/jℓ-cache setup N times today. Addressed by injection, not by
  making the object reentrant (§6).
- **Untested at production scale before this plan — now checked, see §1.3.**
  Everything exercising the Levin integrator previously used toy ℓ ranges
  (`Ncm.SBesselIntegratorLevin.new(0, 8)` in `view.py`; one unit test up to
  1000). §1.3's benchmark is the first time it's been run end-to-end near
  the 1200-3000 range production kernel configs already default to.

### 1.3 Empirical block-size sweet spot: 8, not `ell_cache_max`

Benchmarked directly (`dev-notes/xcor_ultralevin_batching_bench.py`, not
committed to the build — a standalone script): for a fixed total ℓ range
`[0, lmax]`, tile it into contiguous blocks of size `B` and call
`get_eval_vectorized` once per kernel per block (the "prepare all kernels for
this block" half of §5.2 step 2a — no outer k-integral involved, that cost is
identical regardless of blocking and would only dilute the signal). 10
LSST-Y1 kernels (5 `XcorKernelGal` lens bins + 5 `XcorKernelWeakLensing`
source bins), tier 3 (`l_limber=-1`, true non-Limber).

| block `B` | lmax=47, N=5 | lmax=400, N=5 | lmax=400, N=10 |
|---|---|---|---|
| 2  | 16.6s | 70.9s  | 122.8s |
| 4  | 14.4s | 58.9s  | 103.8s |
| **8**  | **14.1s** | **55.0s**  | **94.6s**  |
| 16 | 20.0s | 63.5s  | 107.0s |
| 32 | 23.5s | 77.2s  | 119.2s |
| 64 | 30.5s | 103.2s | 145.2s |

`B=8` is the best or tied-for-best at every `(lmax, N)` combination tested —
a clean U-shape, not a monotonic trend in either direction. Below it, the
fixed per-block overhead (paid once per block regardless of size) dominates,
since smaller `B` means more blocks for the same `lmax`. Above it, a
mechanism found via `perf record` (software events, no root needed;
`_ncm_sbessel_ode_operator_diagonalize_batched` and
`_solve_endpoints_batched_internal` account for 75-91% of samples) explains
the cost growth precisely, not just correlates with it:
`_ncm_sbessel_check_convergence_batched`
([ncm_sbessel_ode_solver.c:1262-1296](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.c#L1262-L1296))
gates the whole batch's convergence on a **max-reduction across all ℓ's in
it** — the number of spectral columns needed (`solution_order`) is set by
the single hardest-to-resolve ℓ in the block, not an average. Combined with
`_ncm_sbessel_apply_rotations_batched` doing one Givens rotation per ℓ per
column (linear in `n_ell`,
[ncm_sbessel_ode_solver.c:1815](numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.c#L1815)),
total cost per call is `O(n_ell × solution_order(n_ell))`, and
`solution_order` itself grows with block width because a wider ℓ-span is
statistically more likely to contain a harder ℓ. That's why
`solve_endpoints`'s profile share overtakes `diagonalize`'s as `B` grows
(45%/30% at `B=8` → 35%/56% at `B=64`): more of the cost sits in re-running
the shared convergence loop out to the block's worst member.

**Threading check.** The above ran with `OMP_NUM_THREADS=3`/
`OPENBLAS_NUM_THREADS=3` (NumCosmo links `libopenblas.so.0`, which is
multi-threaded; the profile's `dptts2_` symbol confirms some LAPACK work is
in the mix). Re-ran the two extremes (`B=8`, `B=64`) at `lmax=400` fully
serial (`OMP_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1`) to check this wasn't a
threading artifact:

| | threaded (3) | serial (1) |
|---|---|---|
| B=8,  N=5  | 55.0s  | 65.8s  |
| B=64, N=5  | 103.2s | 108.9s |
| ratio (64/8) | 1.88× | 1.66× |
| B=8,  N=10 | 94.6s  | 111.6s |
| B=64, N=10 | 145.2s | 153.3s |
| ratio (64/8) | 1.54× | 1.37× |

The `B=8` win holds serially too — the ratio shrinks a little (large `B`'s
BLAS-touching share parallelizes slightly better than small `B`'s), meaning
multi-threaded production deployment makes the `B=8` recommendation, if
anything, understated rather than an artifact of the 3-thread environment
this was benchmarked in. The dominant hot-path functions
(`diagonalize_batched`/`solve_endpoints_batched_internal`) use `#pragma omp
simd` only, never `omp parallel` — architecturally single-threaded
regardless of `OMP_NUM_THREADS` — consistent with the mechanism above being
the real driver, not a scheduling effect.

**Conclusion for §5's block planner: default block size is 8, hard-coded
from this data, not derived from `ell_cache_max`.** `ell_cache_max` (1200)
remains relevant only as the hard ceiling no block may exceed, alongside
`MAX_ELL_BLOCK` (64) — 8 comfortably satisfies both, with wide margin.

**Out of scope here, logged for later:** the max-reduction convergence gate
is a property of the shared `ncm_sbessel_ode_solver.c` primitive, not xcor
glue. An early-exit per-ℓ (stop rotating/counting columns for an ℓ once
*it* converges, instead of running every ℓ in the batch out to the slowest
one) could in principle recover some of the loss at large `B` — but that's
an UltraLevin-level optimization, a separate effort from this plan, and
`B=8` already performs well without it.

**Corrected — the mechanism above is wrong, and the per-ℓ early exit is not
worth doing.** Both claims were measured directly by shadow-tracking, inside
`_ncm_sbessel_check_convergence_batched`, the column at which each ℓ would have
converged on its own, against the column the batch actually ran to (~47k
diagonalizations, N=5 kernels, ℓ=[2,60]):

| B | wasted column-work, mean / p90 | per-ℓ column spread (max/min) |
|---|---|---|
| 8 | 0.2% / 0.0% | 1.01× |
| 16 | 0.4% / 0.3% | 1.02× |
| 32 | 0.7% / 1.3% | 1.04× |
| 64 | 1.2% / 4.1% | 1.08× |

Every ℓ in a block converges at essentially the same column, so a perfect
per-ℓ early exit recovers ~1% even at `B=64` — not worth restructuring the
rotation storage, the `last_n_cols` reuse path and the back-substitution.

`solution_order` also does not grow with block width; it slightly shrinks:

| B | n_ell | mean `solution_order` | calls | total `n_ell × cols` |
|---|---|---|---|---|
| 8 | 8 | 68.7 | 47388 | 24.4M |
| 16 | 16 | 67.1 | 31368 | 31.6M |
| 32 | 32 | 63.9 | 20297 | 38.9M |
| 64 | 59 | 62.4 | 16118 | **59.3M** |

The real reason wide blocks cost more is in the last two columns: per-call cost
is linear in `n_ell`, but widening the block 8× reduces the call count only
**2.9×**, so total work grows 2.4×. Calls are driven by how many **k-samples**
the adaptive closure needs, not by ℓ, and a wider ℓ-block needs more of them —
it spans a wider k-range and a harder function to represent. `B=8` therefore
wins because ℓ-batching stops paying once the extra k-sampling outruns the
amortisation, not because it dodges a convergence-gate penalty.

Consequence: the cost of a wide ℓ-range lives in
`_nc_xcor_kernel_build_spline_integrand`'s adaptive k-sampling, i.e. ordinary
xcor glue, not in the delicate `ncm_sbessel_ode_solver.c`. That is where to
look next. Also unexplained and possibly a better target: the long tail of hard
diagonalizations — mean `solution_order` is 69 but p90 is 133 and the maximum
414.

## 2. Current architecture: three evaluation tiers

`NcXcor` does not have a binary Limber/non-Limber split — it has three, and
they differ in *which* nested integral gets Limber-approximated.

| Tier | Selector | Radial (component) integral | Outer k-integral | Code |
|---|---|---|---|---|
| 1 — Limber-z | `xc->meth == LIMBER_Z_GSL/CUBATURE` (`NcXcor` property, default) | N/A — kernel evaluated as `W(z)` directly, no k-space spline | Replaced by the classical Limber substitution `k = (ℓ+½)/χ(z)`; only the z-integral is numerical | `nc_xcor_kernel_eval_limber_z`, wholly separate from tiers 2/3 |
| 2 — kernel-Limber | `xc->meth == KERNEL_GSL/CUBATURE`, kernel's `l_limber == 0` (`NcXcorKernel` property, **default**, [nc_xcor_kernel.c:120](numcosmo/nc/xcor/nc_xcor_kernel.c#L120)) | Limber-approximated: `_component_states_compute_limber` substitutes `x = ν/k`, giving `W_ℓ(k) ≈ √(π/2ν)·K(ν/k,k)/k` | Exact, numerical (GSL/cubature over k) | `_nc_xcor_kernel_build_limber_integrand` |
| 3 — non-Limber | `xc->meth == KERNEL_GSL/CUBATURE`, kernel's `l_limber < 0` | Exact: `_component_states_compute_non_limber` calls `ncm_sbessel_integrator_integrate` — the real UltraLevin radial Bessel transform | Exact, numerical (GSL/cubature over k) | `_nc_xcor_kernel_build_non_limber_integrand` |

**Tier 1 stays a permanently separate path** — a genuinely different, cheaper
formula (pure z-quadrature, no k-space spline or Bessel machinery at all),
it's what production (`numcosmo_py/app/xcor/`) actually uses today, and it's
well-tested. Everything in §4-§7 unifies tiers 2/3 only; tier 1's code and
behavior are untouched by this plan.

Tiers 2 and 3 share their outer-integral code (`_nc_xcor_kernel_gsl`/
`_nc_xcor_kernel_cubature` in `nc_xcor.c`) — only the per-k callback used to
build the k-space spline differs. §3's correctness bug lives in that shared
code, and hits tier 3 harder than tier 2: the reason to pay for tier 3 is to
capture the diffraction tails *outside* the naive Limber k-band (what the
exponential-decay extrapolation in `_component_states_compute_non_limber`
models); constraining tier 3's outer integral to the Limber band throws away
the thing tier 3 exists to compute, on top of the extrapolation-blowup
failure mode in §3.

**Production is not exposed to any of this today.** `numcosmo_py/app/xcor/kernels.py`
never touches `l_limber`, `meth`, or `integrator` — production always runs
tier 1. Tiers 2/3 are exercised only by the debug tool `view.py`
([view.py:668-693](numcosmo_py/app/xcor/view.py#L668-L693)), the only place
`l_limber` is ever set away from its default. No fire drill, but tiers 2/3
have effectively zero real-world mileage, reinforcing §1.2's scale-testing
gap and §7's production-exposure milestone.

**Known doc/code drift, same "unfinished transition" family:**
`nc_xcor_kernel_get_eval_vectorized`'s docstring
([nc_xcor_kernel.c:1462-1465](numcosmo/nc/xcor/nc_xcor_kernel.c#L1462-L1465))
claims a third branch — "otherwise falls back to single-l get_eval for
lmin" — for a block whose `[lmin,lmax]` straddles the kernel's `l_limber`
threshold. That branch does not exist; only two are implemented
([nc_xcor_kernel.c:1474-1477](numcosmo/nc/xcor/nc_xcor_kernel.c#L1474-L1477)),
so a straddling block today silently goes entirely tier-3 for every ℓ in it,
including the ones ≥ `l_limber`. §5's block planner resolves this by
construction rather than implementing the missing branch: treat each
kernel's `l_limber` as a hard block-tiling boundary, alongside
`ell_cache_max` (§1.2) and the fixed-size `MAX_ELL_BLOCK`/`MAX_COMP_BLOCK`
caps ([nc_xcor_kernel.c:491-492](numcosmo/nc/xcor/nc_xcor_kernel.c#L491-L492)).
Once no block can straddle `l_limber`, the stale docstring branch is
provably unreachable and should be deleted.

## 3. Correctness fix (tiers 2/3)

Numeric smoke test (CMB-lensing auto-Cℓ, ℓ=2–50, comparing tier 1 vs
`KERNEL_GSL`/`KERNEL_CUBATURE`): the kernel-space methods are wrong by
**15-27 orders of magnitude** and disagree with *each other* by
ℓ-dependent factors (ruling out a missing global constant).

**Root cause.** `_nc_xcor_kernel_build_spline_integrand` builds a
`NcmSplineCubicNotaknot` over a domain set by adaptive epsilon-expansion
(`sid->k_min/k_max = ncm_function_sample_set_get_x_min/max`,
[nc_xcor_kernel.c:1048-1049](numcosmo/nc/xcor/nc_xcor_kernel.c#L1048-L1049)).
But `_nc_xcor_kernel_gsl`/`_nc_xcor_kernel_cubature`
([nc_xcor.c:811-812](numcosmo/nc/xcor/nc_xcor.c#L811-L812),
[nc_xcor.c:850-874](numcosmo/nc/xcor/nc_xcor.c#L850-L874)) pick their *outer*
k-integration bounds from `nc_xcor_kernel_get_k_range()`
([nc_xcor_kernel.c:1388-1432](numcosmo/nc/xcor/nc_xcor_kernel.c#L1388-L1432)),
a completely independent Limber-band calculation (`ν/ξ_max ≤ k ≤ ν/ξ_min` ∩
hard component limits) with no guarantee of matching the fitted spline
domain — and, per §2, actively wrong in spirit for tier 3.
`_spline_integrand_eval` never clamps
([nc_xcor_kernel.c:441-453](numcosmo/nc/xcor/nc_xcor_kernel.c#L441-L453)), so
GSL/cubature almost certainly evaluate the not-a-knot cubic spline far
outside its fitted range, where cubic extrapolation diverges — consistent
with the observed magnitudes and the GSL-vs-cubature disagreement (different
adaptive sampling ⇒ different extrapolation garbage).

Test coverage never caught this: `test_xcor_compute_methods`
([test_integration.py:93-133](tests/python/nc/xcor/test_integration.py#L93-L133))
claims (in its docstring) to check `KERNEL_GSL` but never calls it;
`KERNEL_CUBATURE` appears in no test file.

**Fix: make the integrand the source of truth, not `get_k_range()`.** The
machinery for this already exists and is unwired — someone started the fix
and didn't finish. `NcXcorKernelIntegrand` carries a `get_range_func` slot
([nc_xcor_kernel.h:69-98](numcosmo/nc/xcor/nc_xcor_kernel.h#L69-L98)), and the
spline-backed implementation (`_spline_integrand_get_range`,
[nc_xcor_kernel.c:455-462](numcosmo/nc/xcor/nc_xcor_kernel.c#L455-L462)) already
reports the correct thing — the actual fitted domain, `sid->k_min/k_max`. But
`nc_xcor_kernel_integrand_get_range()` ([nc_xcor_kernel.h:218](numcosmo/nc/xcor/nc_xcor_kernel.h#L218))
has **zero callers** outside its own definition.

**Status: DONE.**

1. In `_nc_xcor_kernel_gsl`/`_nc_xcor_kernel_cubature`, build both kernels'
   integrands *before* computing the outer bound (inverts the old
   bound-then-build order): `xclki1 = nc_xcor_kernel_get_eval(...)`,
   `xclki2 = nc_xcor_kernel_get_eval(...)`, then
   `nc_xcor_kernel_integrand_get_range(xclki1, &k1_min, &k1_max)` / same for
   `xclki2`, then intersect exactly as before:
   `k_min = max(k1_min,k2_min)`, `k_max = min(k1_max,k2_max)`. Landed, plus a
   pre-existing integrand refcount leak fixed along the way (`get_eval`/
   `get_eval_vectorized` are `transfer full`; the old code never unreffed).
2. `nc_xcor_kernel_get_k_range()` stays as public API — still useful as a
   cheap seed/estimate (e.g. §1's k-seed computation, or a caller that wants
   a bound without paying for a full spline build) — but stopped being used
   as the outer-integral bound for tiers 2/3.
3. Regression tests landed (`test_xcor_kernel_methods`,
   `tests/python/nc/xcor/test_integration.py`): `KERNEL_GSL`/`KERNEL_CUBATURE`,
   tier 2 and tier 3 separately, agreement with each other and convergence to
   tier 1 at moderate ℓ. Two separate, pre-existing bugs found along the way
   (not fixed, out of scope for this pass — see §9.0): `KERNEL_CUBATURE`
   fatal-asserts for `kernel_cmb_isw` at every tier; `KERNEL_GSL`
   fatal-errors on GSL roundoff for it at tier 3. Both skipped in the test
   with a documented reason, not worked around silently.
4. External validation against CCL done — see §9 for the full comparison
   (correctness, "who's righter", speed, and a deep precision/speed-lever
   investigation this raised). Short version: NumCosmo's tier-3 answer
   agrees with CCL's independent FKEM non-Limber computation to ~0.15% at
   default settings, and convergence testing suggests NumCosmo is the one
   closer to the true converged value.

This landing, and being validated, was a hard prerequisite for trusting any
part of the batching redesign below — there was no point optimizing the
reuse pattern of a computation that returned garbage.

## 4. Where reuse is currently thrown away

`NcXcorKernel` owns one `NcmSBesselIntegrator*`
([nc_xcor_kernel.c:73](numcosmo/nc/xcor/nc_xcor_kernel.c#L73), settable
per-kernel via the `integrator` construct property). Tier 2/3 evaluation goes
through `nc_xcor_kernel_get_eval_vectorized(xclk, cosmo, lmin, lmax)`, which
calls `ncm_sbessel_integrator_set_ell_range` (cheap if unchanged) and then
unconditionally rebuilds the **entire** k-space spline for that ℓ-block via
`_nc_xcor_kernel_build_spline_integrand` — an adaptive `NcmFunctionSampleSet`
construction that queries the (expensive) Levin integrator at many `k`
points. **There is no cache keyed on `(kernel, ℓ-block)`.**

`_nc_xcor_kernel_gsl`/`_nc_xcor_kernel_cubature` call
`nc_xcor_kernel_get_eval_vectorized` once per kernel *per pair*.
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
tomographic analysis needs, each kernel's expensive k-space spline for a
given ℓ-block gets rebuilt once *per pair it participates in* (up to `N-1`
times more than necessary), instead of once per `(kernel, ℓ-block)` and
shared across every pair that touches it in that block. The current
interface has no place to say "here are all the kernels and all the Cℓ's I
want — figure out the cheapest order," because `nc_xcor_compute()` is a
single-pair, single-block call with no memory of sibling calls. §5 fixes
this.

## 5. Batched interface for tiers 2/3

Goal: let the caller declare *everything it wants* up front, so Xcor can pick
an ℓ-block tiling and, for each block, build every distinct kernel's spline
exactly once and reuse it across every pair that needs it, before moving to
the next block.

### 5.1 Phase 1 — registration (configuration only, no computation)

```c
guint nc_xcor_register_kernel   (NcXcor *xc, NcXcorKernel *xclk);
void  nc_xcor_request_cl        (NcXcor *xc, guint kernel_id_1, guint kernel_id_2,
                                  guint lmin, guint lmax);
```

- `register_kernel` is idempotent by pointer identity (registering the same
  kernel twice returns the same id).
- `request_cl` declares one desired output block (a pair + ℓ-range); call it
  once per entry in the data vector / covariance block the caller ultimately
  needs (mirrors what `NcDataXcor` already knows per `xcab[a][b]`, including
  `ell_th_cut_off`).

### 5.2 Phase 2 — solve (all computation happens here)

```c
void nc_xcor_solve (NcXcor *xc, NcHICosmo *cosmo);
NcmVector *nc_xcor_get_cl (NcXcor *xc, guint kernel_id_1, guint kernel_id_2);
```

`nc_xcor_solve` internally:

1. Computes the union of ℓ-blocks needed across all registered requests,
   tiling `[ℓmin,ℓmax]` at a **default block size of 8** (§1.3, empirically
   determined — not derived from `ell_cache_max`), always
   `≤ ncm_sbessel_integrator_levin_get_ell_cache_max()` and
   `≤ MAX_ELL_BLOCK` as hard ceilings, and never straddling any referenced
   kernel's `l_limber` threshold (§2). Union-then-mask (compute the full
   block, use only the requested sub-range per pair) is the simplest tiling
   strategy when different pairs need different `[ℓmin,ℓmax]`; it may
   compute some unused ℓ's for a given pair, traded for not needing per-pair
   block grids.
2. For each ℓ-block, in order:
   a. For each *distinct* kernel referenced by a request touching this
      block, build its k-space spline once
      (`nc_xcor_kernel_get_eval_vectorized`), cached for the block's
      lifetime only. Splines depend on the current cosmology (changes every
      likelihood/MCMC step), so caching never needs to outlive one
      `nc_xcor_solve` call.
   b. For each request whose ℓ-range intersects this block, compute the
      k-integral of the (cached) kernel-product and write the corresponding
      slice of that pair's output vector.
   c. Drop the block's cached splines before moving to the next block.
3. Results become retrievable via `nc_xcor_get_cl`.

This turns the cost from `O(N_pairs × N_blocks)` kernel-spline builds into
`O(N_kernels × N_blocks)` — the actual optimization target.

`nc_xcor_compute()` stays as today's convenience wrapper (register both
kernels + one request + solve, for the single-pair case), documented as the
non-optimal path once more than a couple of kernels are involved.

**Open question:** does this stateful register/request/solve API live on
`NcXcor` itself, or a new companion object (e.g. `NcXcorSolver`)? `NcXcor`
today is close to stateless per-call; bolting registration/request state
onto it changes its contract. A companion object keeps `nc_xcor_compute()`'s
existing callers untouched — leaning toward this, to be confirmed before
§7 milestone 3 starts.

## 6. Integrator injection and concurrency

Step 5.2.2a is where per-kernel spline construction happens. The first draft
of this section recommended giving every concurrently-evaluated kernel its
own fresh `NcmSBesselIntegratorLevin` (one slot per kernel in an
`NcmMemoryPool`) and parallelizing across all of them. **§6.1 below found
that recommendation wrong for at least one important kernel type** — it
trades away a real, large reuse benefit that the naive "just parallelize
everything" test happened to obscure. Read §6.1 before implementing
anything here; the design that follows it (§6.2) is the corrected one.

### 6.1 Cross-kernel reuse is real, large, and *order-sensitive*

`_ncm_sbessel_ode_operator_diagonalize_batched` has an early-out
(`op->last_n_cols > 0` → `_ncm_sbessel_apply_all_stored_rotations_batched`)
that reuses a previously-computed Givens-rotation decomposition *from the
same integrator instance* — cheap (docstring: "4 FLOPs") when the stored
decomposition already covers the new kernel's RHS, falling through to
expensive extension or full diagonalization when it doesn't. This is a
property of the **operator** (ℓ-block, panel), not of any specific kernel's
`K(x,k)` — so in principle it should transfer across *different* kernels
sharing one integrator instance, not just across repeated k-queries within
one kernel's own adaptive sampling.

Tested directly (10 distinct kernels, block `[0,7]`, tier 3; first with a
mixed lens+source group, which — see below — was itself a misleading test):

| group | shared integrator (serial) | fresh integrator (serial) | ratio |
|---|---|---|---|
| 5 lens bins | 5321ms | 3724ms | fresh **0.70×** (faster) |
| 5 source bins | 443ms | 1241ms | fresh **2.80×** (slower) |

Opposite conclusions depending on kernel type — worth distrusting on sight,
so re-checked with per-kernel timings, fresh `dndz` splines per kernel (no
shared-object confound), and interleaved run order (rules out simple
warm-up bias):

```
source, shared: [272, 50, 47, 41, 32] ms   <- kernel 1 full cost, 2-5 ~8x cheaper
source, fresh:   [272, 244, 243, 255, 234] ms  <- every kernel full cost

lens, shared (ascending bin order): [273, 346, 217, 1809, 2535] ms
lens, fresh:                        [277, 508, 395, 999, 1862] ms
```

Source confirms the reuse mechanism working exactly as designed: kernel 1
pays full diagonalization, kernels 2-5 land on the cheap stored-rotation
path. Lens looked like reuse *hurts* — until the per-kernel numbers showed
why: lens bins get harder (need more spectral columns) as `bin_idx`
increases, so processing them in ascending order means every "harder" bin
has to *extend* a too-small stored decomposition, which costs more than
just diagonalizing that bin fresh. **Reuse only helps when a later kernel's
requirement is a subset of what's already stored** — extending upward is
expensive; reusing downward (or sideways, same difficulty) is nearly free.

Confirmed by fixing the order: processing lens bins **hardest-first**
(descending, `[4,3,2,1,0]`) instead of ascending drops the shared-integrator
total from 5102ms to 3517ms — beating even the fresh/no-reuse total
(~4041ms). So this isn't "lens kernels don't benefit from reuse," it's "lens
kernels benefit from reuse *only with the right processing order*," same as
source.

**Revised conclusion: don't give every kernel its own integrator.** Group
kernels that are likely to share a decomposition well, process each group
**serially against one shared integrator, hardest-member-first**, and
parallelize **across groups** (§6.2), not across every individual kernel.
Blanket per-kernel parallelism (the original §6 design, briefly validated in
an earlier draft of this section with a since-corrected 10-kernel mixed test
showing an apparent 1.76-2.69× speedup) is a straightforward net loss for
kernel types like weak lensing, where within-group reuse is large: parallel
fresh-integrator evaluation of the 5 source bins alone (956-1045ms) was
*slower* than serial shared-integrator evaluation of the same 5 (464ms) —
nearly 2× worse, not better.

**Open, unresolved by this plan:** how does the real implementation know
which kernel is "hardest" (needs the most spectral columns) *before*
evaluating any of them, in order to schedule it first? Candidate approaches,
none tested here: (a) a cheap proxy heuristic (e.g. dn/dz width/sharpness,
redshift range) computed without a full closure build; (b) remember which
kernel was hardest on the *previous* `nc_xcor_solve()` call (e.g. the
previous MCMC step) and schedule it first on the next one, since cosmology
changes only incrementally between calls; (c) accept the ascending-order
penalty as a one-time cost per group and not optimize scheduling order at
all, if profiling (§7 milestone 6) shows it's not the dominant cost anyway.
This needs real design work, not just an `NcmMemoryPool` of interchangeable
integrators.

### 6.2 Revised design: grouped sharing, parallel across groups

```c
void nc_xcor_set_integrator        (NcXcor *xc, NcmSBesselIntegrator *sbi);
void nc_xcor_set_integrator_pool   (NcXcor *xc, NcmMemoryPool *pool);
```

- `set_integrator`: a single shared instance, used serially across *all*
  kernels in a block, in whatever order the caller registered them — today's
  implicit behavior, and the right default for a one-off script. Per §6.1,
  this is not just "the simple case" — with good ordering it can beat naive
  per-kernel parallelism outright.
- `set_integrator_pool`: an `NcmMemoryPool` (same idiom already used for GSL
  workspaces in [ncm_integrate.c:75-88](numcosmo/ncm/integration/ncm_integrate.c#L75-L88)),
  but pool slots now correspond to **kernel groups**, not individual
  kernels. Each slot's integrator is shared serially (hardest-first, §6.1)
  by every kernel routed to that group; different groups run concurrently
  against different pool slots. Grouping policy (by `GType`, by some
  cheap difficulty proxy, or something else) is exactly the open question
  from §6.1 — not resolved here. A single-group, single-slot pool
  degenerates to exactly the `set_integrator` case for free.

Don't make `NcmSBesselIntegratorLevin` itself reentrant — its scratch state
is deeply threaded through the hot path, not worth it. Xcor stays agnostic
to *how much* concurrency exists; it only needs an integrator when step
5.2.2a needs one, and asks the injected source for it.

Two mechanical findings from testing this, still valid regardless of the
grouping policy above: PyGObject releases the GIL during these C calls (so
real threading works, verified via genuine wall-time reduction, not just
GIL-serialized busywork), and the correct pattern for combining with
NumCosmo's linked `libopenblas.so.0` (itself multi-threaded) is **outer
parallelism, inner serial BLAS** — pin `OPENBLAS_NUM_THREADS=1` per worker
to avoid oversubscription, don't leave it implicit.

### 6.3 Block-level parallelism: the primary axis, not kernel-level

Different ℓ-blocks share **no** reusable operator state with each other —
`set_ell_range` to a different block is a full rebuild regardless (the ODE
coefficient depends on ℓ) — so unlike kernel-level parallelism (§6.1),
parallelizing across blocks has no reuse-vs-parallelism tradeoff to reason
about. Each concurrently-run block needs its own dedicated integrator (and
kernel instances bound to it — can't have two blocks racing to mutate one
integrator's `set_ell_range` state), but that's the only constraint.

Tested directly: 8 blocks of size 8, each doing the full 10-kernel (hardest
lens bins first, per §6.1) closure build, correctness checked (per-block
k-ranges must match a serial run over the same 8 blocks):

| workers | wall time | speedup | correct? |
|---|---|---|---|
| serial | 19427ms | 1.0× | — |
| 2 | 12409ms | 1.57× | ✓ |
| 4 | 8331ms | 2.33× | ✓ |
| 8 | 7517ms | 2.58× | ✓ |

Clean win, no caveats beyond the usual (sandbox-capped core count, so the
real ceiling is untested here). More importantly, this is the axis with far
more headroom: a production run has `ell_cache_max`-bounded but numerous
blocks (51 at `lmax=400`/`block=8`, §1.3) vs. only ~10-20 kernels — block
parallelism alone offers more exploitable concurrency than kernel-level ever
could, without any of §6.1's hardest-first ordering complexity.

**Revised recommendation:** parallelize primarily **across blocks** (maps
naturally onto an OpenMP loop over step 5.1's block list); keep processing
**within** a block simple — serial, one shared integrator, hardest-kernel-
first (§6.1) — rather than *also* parallelizing kernels within a block. Only
reconsider kernel-level parallelism (§6.1/6.2) if profiling (§7 milestone 6)
shows too few blocks to saturate available cores (e.g. a narrow `lmax`
range) — the common case, more blocks than cores, doesn't need it.

Whether tier-2 kernel-spline-building or tier-3 pair-integration dominates
overall wall time at realistic N is still a profiling question for after a
correct, batched baseline exists (§7 milestone 6).

## 7. Milestones

1. **Correctness fix (§3) — DONE.** Reordered both `_nc_xcor_kernel_gsl` and
   `_nc_xcor_kernel_cubature` to build the integrand(s) first and read the
   outer bound off `nc_xcor_kernel_integrand_get_range()`. `KERNEL_GSL`/
   `KERNEL_CUBATURE` regression tests landed for tier 2 and tier 3 separately
   (agreement with each other, convergence to tier 1). External CCL
   validation done, see §9. **Still open:** delete the stale third-branch
   docstring on `get_eval_vectorized` once milestone 3's block planner makes
   it provably unreachable; two separate pre-existing bugs found for
   `kernel_cmb_isw` specifically (§3 point 3, §9.0) remain unfixed.
2. **Scale validation (§1.3) — mostly done.** Timing+stability confirmed up
   to `lmax=400` (§1.3's table) and the block-size question is answered
   (default 8). Along the way, found and fixed a real, separate bug this
   validation would otherwise have hit immediately: `ncm_function_sample_set`
   (a shared NCM primitive, not xcor-specific) could never converge when any
   component of a multi-ℓ block is identically zero — true for weak lensing
   at ℓ=0,1 (spin-2 field, no monopole/dipole) — causing exponential sample-
   count blowup and OOM. Fixed in `ncm_function_sample_set_refine`/
   `get_absmaxF_min` with 5 new regression tests
   (`tests/python/ncm/stats/test_function_sample_set.py`); verified against
   the exact originally-crashing calls. *Accuracy* (not just timing/
   stability) checked against CCL, see §9 — agreement to ~0.15%, and
   convergence testing suggests NumCosmo is the more accurate of the two.
   **Still open:** the full `lmax=1200-3000` production range — `lmax=400`
   was chosen for benchmark turnaround time, not because 1200-3000 is
   expected to behave differently in kind (block=8's tiling makes per-call
   cost independent of the overall `lmax` ceiling anyway).
3. **Batched interface (§5-§6) — DONE.** Landed as `NcXcorSolver`, a
   companion object (§5's open placement question resolved this way):
   `register_kernel`/`request_cl`/`plan_blocks`/`solve`/`get_result`;
   `nc_xcor_compute()` itself is unchanged, still the non-batched
   convenience path. `solve()`'s `KERNEL_CUBATURE` requests share each
   kernel's k-space closure once per ℓ-block across every request needing
   it in that block (the actual optimization target); every other method
   delegates to `nc_xcor_compute()` per request, correct but unoptimized
   (tier 1 untouched by design, `KERNEL_GSL` already has no block-shared
   closure to reuse). Blocks are parallelized primarily across an OpenMP
   team (§6.3), each thread working from its own `NcmSerialize`-duplicated
   kernel clones so no thread touches another's integrator state; §6.1's
   kernel-grouping/hardest-first-ordering question stays unresolved/
   deferred (only matters if too few blocks exist to saturate cores).
   Regression-tested against `nc_xcor_compute()`, including multi-kernel
   closure reuse and multi-block stitching, stress-tested for the parallel
   path at 1/2/4/8 OMP threads.
4. **Migrate `NcDataXcor` — SKIPPED.** Investigated and found this class to
   be untested, dead code: no tests exist for it anywhere in the repo, and
   nothing in `numcosmo_py` (including `app/xcor`) constructs or uses it —
   the plan's assumption that it was "the one production consumer that
   benefits" doesn't hold in this codebase. Migrating an untested class
   with no way to verify the result end-to-end isn't worth the risk for a
   performance win nothing currently exercises; revisit if/when something
   comes to depend on `NcDataXcor` being fast.
5. **Production config exposure — done via milestone 1's `view.py` work**
   (§7 milestone 1 exposed `l-limber`/`integrator`
   reltol/cheb-reltol/max-order as CLI options). No separate non-`view.py`
   "production config" path exists to extend, since `NcDataXcor` (the
   would-be production consumer) isn't wired into `numcosmo_py` at all —
   see milestone 4.
6. **Performance pass — DONE, see §9.8** for the CCL scaling comparison
   against the real `NcXcorSolver` (not a Python-level simulation):
   14.83×→7.58× (N=5→10, default precision, trustworthy accuracy 0.3-0.6%),
   confirming the `O(N)` vs `O(N²)` scaling shape. Along the way, found and
   fixed a real heap-corruption bug in `_ensure_operator_size`
   (`ncm_sbessel_ode_solver.c`) exposed by tier-3 kernel duplication across
   a shrinking ℓ-block, and corrected §9.4's "tuned config" claim (doesn't
   generalize, was wrong by up to 23-29% outside the one case it was
   measured on). **Still open:** N beyond 10; block-level OpenMP effect on
   these specific ratios (measured single-threaded here); docs (a theory
   page for the ODE-solver/Levin combination, `docs/theory/`) not started.

## 8. Non-goals

- **No changes to tier 1** (`LIMBER_Z_GSL`/`LIMBER_Z_CUBATURE`) — code,
  behavior, or performance. It stays a separate path (§2).
- **No reentrant `NcmSBesselIntegratorLevin`.** Concurrency is handled by
  injection (§6), not by making the object thread-safe internally.
- **No cross-`nc_xcor_solve()`-call spline caching.** Splines are keyed to a
  cosmology that changes every likelihood/MCMC step; caching beyond one
  `solve()` call has no use case identified so far.

## 9. External validation (CCL) and precision/speed tuning

Uses `numcosmo_py.ccl.nc_ccl.create_nc_obj(ccl_cosmo, ...)` — existing
paired-object infrastructure — so NumCosmo and CCL cosmological parameters
match exactly. Test case throughout: a synthetic Gaussian-dndz weak-lensing
kernel, tier 3 (`l_limber=-1`), auto-correlation, `KERNEL_GSL`/
`KERNEL_CUBATURE`.

### 9.0 Two more pre-existing bugs found, not fixed here

Same family as §3's — found while building this section's tests, unrelated
to the §3 fix (which supplies a valid, sensible outer bound in every case
here too):

- **`KERNEL_CUBATURE` aborts the process** for `kernel_cmb_isw` specifically,
  at every ℓ tier: fatal `g_assert` in `ncm_integral_nd_eval`
  (`"assertion failed (ret == 0)"` — the underlying cubature call itself
  returns a nonzero/error status, which `ncm_integral_nd_eval` treats as
  fatal). `wl`/`gal`/`cmb_lens`/`tsz` all work fine with `KERNEL_CUBATURE`.
- **`KERNEL_GSL` fatal-errors on `kernel_cmb_isw` at tier 3** specifically
  (works fine at tier 2): `_nc_xcor_kernel_gsl: roundoff error` — GSL's
  `qag` hits `GSL_EROUND`, and `_nc_xcor_kernel_gsl` treats *any*
  non-`GSL_SUCCESS` status as fatal via `g_error`. `GSL_EROUND` often just
  means the result can't be formally certified to the requested tolerance,
  not that it's wrong — this fatal-on-any-error policy is itself fragile,
  independent of whether ISW's kernel shape is the real cause.

Both are process aborts, not catchable exceptions. Both newly visible only
because neither method had ever been exercised against ISW in any test
before this session. Regression tests (`test_xcor_kernel_methods`) skip
these combinations explicitly, with the reason documented in the test file
rather than silently worked around. TODO, out of scope for this plan:
root-cause both for ISW's kernel shape; reconsider whether
`_nc_xcor_kernel_gsl` should tolerate `GSL_EROUND` instead of aborting.

A third bug was found *and* fixed: kernel fixtures in
`tests/python/fixtures_xcor.py` share ONE session-wide, `functools.cache`d
`NcmSBesselIntegratorLevin` across every test file. Exercising tier 3 (real
Levin/ODE-solve) against it for the first time left it in a state that
crashed an unrelated, later test in a different file (malloc-corruption
abort inside `ncm_function_sample_set_expand_domain`, only reproducible
running the full `nc/xcor` test directory, or `test_integration.py` followed
by `test_kernel.py`). Fixed by swapping in dedicated, test-local integrators
for the duration of the new tests instead of touching the shared one.

### 9.1 Correctness: NumCosmo vs. CCL/FKEM non-Limber

`ells=[2,4,8,16,24]`, default precision both sides:

| ℓ | CCL non-Limber | NumCosmo non-Limber | ratio | CCL Limber | NumCosmo Limber |
|---|---|---|---|---|---|
| 2 | 1.772849e-08 | 1.775498e-08 | 1.00149 | 1.523274e-08 | 1.523273e-08 |
| 4 | 2.425802e-08 | 2.426623e-08 | 1.00034 | 2.312258e-08 | 2.312257e-08 |
| 8 | 2.318400e-08 | 2.318996e-08 | 1.00026 | 2.286964e-08 | 2.286961e-08 |
| 16 | 1.646916e-08 | 1.647263e-08 | 1.00021 | 1.640978e-08 | 1.640954e-08 |
| 24 | 1.174147e-08 | 1.174400e-08 | 1.00022 | 1.172276e-08 | 1.172272e-08 |

Tier 1 (Limber) agrees with CCL's Limber to 1.5e-5 — expected, both compute
the identical classical approximation, and this is the same comparison
`test_ccl_comparison.py` already covers. Tier 3 (non-Limber) agrees with
CCL's independent FKEM non-Limber computation to **~0.15%** across ℓ=2-24 —
this is the real, new validation: not internal self-consistency (GSL vs.
cubature agreeing with each other could both be wrong the same way), but
agreement with a *different codebase* using a *different non-Limber
algorithm* (FKEM vs. UltraLevin/Chebyshev-Gegenbauer spectral collocation).

### 9.2 Who's righter: convergence study

Tightened each method's own precision knob independently and watched which
answer moved:

| | default | tightened | shift |
|---|---|---|---|
| NumCosmo (kernel `reltol` 1e-4→1e-8) | 1.775498e-08 | 1.775498e-08 | **1.2e-7** (none) |
| CCL/FKEM (`fkem_Nchi` default→2000→8000) | 1.772849e-08 | 1.772622e-08 | **1.3e-4**, still moving |

NumCosmo's answer is completely stable across 4 orders of magnitude of
tolerance tightening — default settings were already fully converged.
CCL/FKEM's answer is still drifting even after quadrupling its grid
(`Nchi` 2000→8000), and drifting *toward* NumCosmo's value. Cross-checking
both at their tightest tested settings still leaves ~0.16% disagreement,
which doesn't shrink at anywhere near the rate you'd expect if NumCosmo were
the one that was wrong. **Read: NumCosmo is righter.** FKEM is a fixed,
finite log-spaced-grid method (its precision-control story is literally a
`CCLWarning: Nchi must be a positive integer. Setting to match tracer with
large chi samples x2` heuristic, not a continuous adaptive tolerance)
carrying real, still-shrinking discretization error at the settings tested;
NumCosmo's adaptive reltol-based refinement has actually converged.

### 9.3 Speed: CCL is 70-250x faster at default settings, and it isn't parallelism

**This section is single-kernel, single-pair — §9.7 shows the gap narrows
sharply (to ~5-10x and closing) once realistic kernel/pair counts are
in play. Read this as "per-pair cost," not "the real-world gap."**

| | time (3 ells) |
|---|---|
| CCL/FKEM, default `Nchi` | 1.6-1.9 ms |
| CCL/FKEM, `Nchi=8000` | 6.1-7.2 ms |
| NumCosmo `KERNEL_GSL`, default precision, unbatched | ~415-431 ms |

CCL's C backend (`ccl_cls.c`) does have real `#pragma omp parallel`/
`omp for schedule(dynamic)` (most likely over the classic-Limber ℓ loop),
linked against `libgomp`. Forcing `OMP_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1`
for the whole process left CCL's FKEM timing essentially unchanged — its
actual non-Limber code path (`pyccl/nonlimber/_nonlimber_FKEM.py`) appears
to be a largely separate Python/NumPy routine, not the OpenMP-parallelized
Limber C code. The speed gap is real, not a hidden-parallelism artifact.

### 9.4 NumCosmo's precision-control knobs: full inventory, and which ones matter

The precision/speed tradeoff is spread across **four separate objects**,
found by enumerating every `g_param_spec` in the relevant files rather than
assuming:

- `NcXcorKernel`: `reltol`, `scaled-abstol` (k-space closure adaptive-
  midpoint refinement), `adaptive-epsilon`, `adaptive-boundary-tries`,
  `max-border-expansions`, `expansion-factor` (domain-expansion/boundary-
  finding phase), `max-iter`.
- `NcXcor`: `reltol` (outer k-integral, default `NC_XCOR_PRECISION = 1e-6`,
  further scaled by `1e-2` inside `_nc_xcor_kernel_gsl`).
- **`NcmSBesselIntegratorLevin`**: `reltol` (default **1e-13**),
  `cheb-reltol` (default **1e-8**), `max-order` (default **16384**),
  `cheb-min-order`, `n-knots`, `y-knots-min`/`y-knots-max`. Near machine
  precision by design — this is the object the plan's own §1.1 identifies
  as the expensive part (panel factorization).
- `NcmSBesselOdeSolver`: `tolerance` — not directly reachable from Python
  (no injection hook, per §1.2/§9.0's "no reentrant" non-goal), but
  confirmed in `ncm_sbessel_integrator_levin_set_reltol()`
  ([ncm_sbessel_integrator_levin.c:1214](numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.c#L1214))
  that `NcmSBesselIntegratorLevin.reltol` propagates to it internally, so
  it's controlled indirectly.

Measured (single kernel, 23 ℓ's `[2,24]`, unbatched `KERNEL_GSL` unless
noted; deviation relative to fully-converged reference):

| knob | result |
|---|---|
| `NcXcorKernel.reltol` 1e-4→1e-1 (1000x looser) | 415→385 ms (**~7%**), dev ≤7.8e-6 — **not a real lever**, cost is elsewhere |
| `NcmSBesselIntegratorLevin.reltol`/`cheb-reltol`/`max-order`, default→(1e-1,1e-1,32) | 421→164.5 ms (**2.6x**), dev 1.2e-4 — **the real lever** |
| `NcXcorKernel.adaptive-epsilon` 1e-3→1e-2→1e-1 | 333→302→271 ms, dev 4.9e-4→6.3e-3→**5.2e-2** — works, but accuracy degrades faster than it's worth |
| `NcXcorKernel.max-border-expansions` 500(default)→200→50 | no effect at all (dev=0) — default is already non-binding | 
| `NcXcorKernel.max-border-expansions`→10 | 271 ms but dev=**0.20** (20% error) — **not a real lever**, just truncates the domain and gives a wrong answer, not a genuine accuracy/speed tradeoff |

**Conclusion: the integrator's own `reltol`/`cheb-reltol`/`max-order` is the
only knob that trades meaningful speed for negligible accuracy loss.**
Everything at the `NcXcorKernel` level either does nothing (kernel `reltol`,
`max-border-expansions` below its natural convergence point) or trades
accuracy faster than it saves time (`adaptive-epsilon`).

**Correction (§9.8): the "negligible accuracy loss" part of this
conclusion does not generalize.** The specific settings tried here
((1e-1,1e-1,32), dev 1.2e-4) were measured on one easy single-kernel case
only; §9.8 found the similar-but-looser (1e-2,1e-2,64) config used
downstream in §9.5/9.7 gives up to 23-29% error on other kernels/ℓ. Treat
any `reltol`/`cheb-reltol`/`max-order` looser than default as unvalidated
until checked per-case.

### 9.5 ℓ-batching helps a single kernel's own closure too, not just cross-kernel reuse

§1.3/§6.3 established batching for *reuse across different kernels* sharing
a block. Separately: `KERNEL_GSL` builds each ℓ's closure independently
(`nc_xcor_kernel_get_eval`, i.e. `get_eval_vectorized(l,l)` — effectively
`block=1`, unbatched) where `KERNEL_CUBATURE` batches via
`ell_batch_size`/`get_eval_vectorized(lmin,lmax)`. For **one kernel's own**
23-ℓ non-Limber closure, batching alone (no precision change) already helps:

| `ell_batch_size` | time | dev from unbatched |
|---|---|---|
| 1 (`KERNEL_GSL`, effectively) | 421 ms | — |
| 8 | 304-307 ms | 7.6e-6 |
| 16 | 270 ms | 9.8e-6 |
| 23 (all ℓ's at once) | 270-271 ms | 9.9e-6 |

Diminishing returns past 16 for this single-kernel case (unlike §1.3's
cross-kernel finding where 8 was the clear minimum against a U-shaped
curve) — makes sense, there's no "worst-ℓ-gates-convergence" cost inflation
here since it's the same kernel throughout, just less fixed per-block
overhead to amortize as blocks get bigger.

**Combining both real levers** (`ell_batch_size=23` + integrator
`reltol=1e-2`/`cheb-reltol=1e-2`/`max-order=64`): **151.6 ms**, dev only
**1.1e-5** — the best configuration found, ~2.8x faster than the naive
default with excellent accuracy retained (much better than any single lever
alone, and far better accuracy-per-ms than loosening `adaptive-epsilon`).

### 9.6 The outer k-integral is negligible — "fixed knots" there would not help

Isolated cleanly (fresh kernel + fresh, never-used integrator for each
measurement, to rule out the cross-call decomposition-reuse confound from
§6.1):

| | time |
|---|---|
| closure build alone (23 ℓ's) | 153.76 ms |
| full compute (closure build + outer k-integral) | 153.27 ms |
| implied outer-k-integral share | **~0%**, noise-level |

All measured cost is in the closure build (the Levin/ODE-solve machinery);
the outer cubature/GSL integral over the already-built spline costs
essentially nothing in comparison — consistent with §1.3's `perf` profiling
(`diagonalize_batched`/`solve_endpoints_batched_internal` at 75-91% of
samples). **A fixed-knots (non-adaptive) scheme for the outer k-integral
would not help** — there is nothing to gain there. `NcmIntegralND`'s four
methods (`CUBATURE_H`/`_P`/`_H_V`/`_P_V`) are all adaptive regardless (the
"_V" suffix means vectorized-integrand, already in use via `KERNEL_CUBATURE`
— not a fixed-node mode); implementing genuine fixed knots would mean
bypassing `NcmIntegralND` for this integral entirely, and §9.6's isolation
says that effort has no payoff here.

### 9.7 §9.1-9.6 were all single-kernel, single-pair — the gap narrows sharply at realistic N

Everything through §9.6 measured one kernel correlated with itself: one
closure build, one pair, one outer integral. That's not what a tomographic
analysis needs — `N` kernels, `N(N+1)/2` auto+cross pairs. §9.6 already
implies the shape of the answer: NumCosmo's real cost is **per kernel**
(closure build), and its outer integral is **per pair but ~free** — so total
NumCosmo cost should scale ~`O(N)`. Whether CCL's cost is `O(N)` or
`O(N(N+1)/2)` depends on whether `_nonlimber_FKEM` caches anything across
calls sharing a tracer — checked: no `@lru_cache`/memoization decorator
anywhere in `pyccl/nonlimber/_nonlimber_FKEM.py`, so it recomputes from
scratch every `pyccl.angular_cl()` call — total CCL cost should scale
~`O(N(N+1)/2)`.

Tested directly, `N=5` and `N=10` distinct weak-lensing kernels (shifted
Gaussian dndz's as tomographic-bin stand-ins), all auto+cross pairs, best
NumCosmo config (§9.4/9.5: integrator `reltol=1e-2`/`cheb-reltol=1e-2`/
`max-order=64`), pair integral via a **fixed, log-spaced 48-node** Simpson
sum evaluating all 23 ℓ's per k-node in one `eval_array` call (per §9.6,
native C already does this essentially for free, but no C API exists today
to reuse a *pre-built* closure across multiple pairs — §4/§5's whole point —
so this is a Python-level prototype of that capability, not the real thing):

| N | pairs | NumCosmo closures | NumCosmo pair integrals | NumCosmo total | CCL total | **ratio** |
|---|---|---|---|---|---|---|
| 5 | 15 | 1171.3 ms | 4.3 ms (0.28 ms/pair) | 1175.6 ms | 120.2 ms | **9.78×** |
| 10 | 55 | 2176.3 ms | 14.0 ms (0.25 ms/pair) | 2190.3 ms | 424.1 ms | **5.16×** |

The gap **nearly halves as N doubles** — exactly the `O(N)` vs. `O(N²)`
scaling predicted, and a fundamentally different picture from §9.3's
single-pair 70-250×. Extrapolating: N=20 → ~2.6×, and it's plausible the
ratio crosses below 1 (NumCosmo faster *in total*) somewhere in the
tens-of-kernels range real tomographic analyses actually use. Pair-integral
cost is genuinely negligible here too (0.25-0.28ms/pair) — confirms §9.6
was not an artifact of the single-pair case.

**Methodology lesson, worth keeping:** the first attempt at this measurement
used `scipy.integrate.quad` called separately per ℓ (23 adaptive integrals
per pair, each doing many Python/GI-marshaled `eval_array` calls) and got
793.7ms for the same 15 pairs — 180x worse than the fixed-knot version,
and it *looked like* the outer integral was expensive, flatly contradicting
§9.6. It wasn't a real cost, it was per-sample-point Python/FFI overhead
compounding across 23 separate adaptive-refinement loops. Switching to fixed
log-spaced knots evaluated once per k-node for *all* ℓ's together (matching
what the native C path already does) made it match §9.6 immediately. Do not
trust a Python-level cost measurement that contradicts an already-verified
native-C one without checking for exactly this kind of per-call overhead
compounding first.

**Not yet resolved:** accuracy at N=5/N=10 came out at 1.46e-2 (worse than
the single-pair 1.5e-3, and identical between the two N's, suggesting one
specific — likely widest-separated — pair dominates the error) — the fixed
48-knot count may be under-resolving that pair; needs tuning/verification
before trusting the accuracy number quantitatively, though it doesn't change
the timing/scaling conclusion. A real C-level "cached closures, cheap C
pair integral" implementation (i.e. actually building §4/§5) is the only way
to get a trustworthy absolute number here rather than a Python prototype.

**Untested:** `NcmSBesselIntegratorLevin`'s `n-knots` (default 21, panel
count) and `y-knots-min`/`y-knots-max` (default 1e-4 to 1e6 — 10 decades,
the dimensionless-`y=kx` panel range). Given §9.4's finding that the
integrator-level knobs are where the real cost lives, these are the next
candidates to probe if further speedup is wanted.

### 9.8 §9.7 redone with the real `NcXcorSolver`, and a correction to §9.4's "tuned config"

Milestone 3 (§7) is now done — a real C-level `NcXcorSolver.solve()` exists
(register/request/plan_blocks/solve, per-block shared closures, OpenMP
across blocks), so §9.7's "not yet resolved" call for "a real C-level
cached-closures implementation... the only way to get a trustworthy
absolute number" can finally be answered directly instead of prototyped.

**§9.4's "tuned config" doesn't generalize — corrected here.** Rerunning
§9.7's exact N=5 case through the real `solve()` with §9.4's tuned
integrator settings (`reltol=1e-2`/`cheb-reltol=1e-2`/`max-order=64`)
gives **29% max deviation** from `nc_xcor_compute()` on the *same*
kernels/request — not the ~1.5% §9.7 flagged as "not yet resolved," an
order of magnitude worse. Isolated to a single kernel (N=1, no cross-kernel
sharing involved at all, so not a `solve()`-specific bug): at §9.4's tuned
settings, `solve()` and `nc_xcor_compute()` individually diverge from each
other by up to 23% at the high-ℓ end of a 24-ℓ block; at **default**
precision, the same comparison agrees to **1.2e-12** (machine precision).
§9.4's own dev number (1.2e-4) was measured on one easy single-kernel case
and never re-checked elsewhere — it does not bound the error generally.
**Correction: don't use `reltol=1e-2`/`cheb-reltol=1e-2`/`max-order=64`
for anything where the answer matters — only default precision has been
shown trustworthy so far.** §9.5's "best configuration" (§9.5's closing
paragraph) inherits the same caveat.

**Real `solve()` vs CCL, default precision** (same setup as §9.7: N
shifted-Gaussian-dndz weak-lensing kernels, all auto+cross pairs, ℓ=[2,25]
— 24 ℓ's, an exact multiple of the block size to sidestep GH #299 — single
thread, `OMP_NUM_THREADS=1`):

| N | pairs | NumCosmo `solve()` total | CCL total | ratio | max dev vs CCL |
|---|---|---|---|---|---|
| 5 | 15 | 1763.7 ms (352.7 ms/kernel) | 119.0 ms (7.93 ms/pair) | **14.83×** | 3.0e-3 |
| 10 | 55 | 3272.8 ms (327.3 ms/kernel) | 431.6 ms (7.85 ms/pair) | **7.58×** | 5.7e-3 |

**Superseded by §9.9** — both the ratios (NumCosmo has since got ~3.5× faster
on this exact benchmark) and the accuracy column (a max over ℓ=[2,25], which
ℓ=2 sets on its own).

The gap **nearly halves as N doubles** (14.83 → 7.58, ratio 1.96) — the
same `O(N)` vs `O(N²)` shape §9.7's prototype predicted, now confirmed with
the actual C-level implementation rather than a hand-rolled pair-integral
standing in for it, and at accuracy (0.3-0.6%) an order of magnitude better
than §9.7's flagged-untrustworthy 1.46e-2. Absolute ratios are worse than
§9.7's prototype numbers (14.83×/7.58× vs. 9.78×/5.16×) — not a regression,
the comparison basis changed: §9.7 used the unsafe tuned config (fast but,
per the correction above, wrong by up to 23% for some pairs) baked into its
NumCosmo side, where this table uses default precision throughout. These
default-precision numbers are the trustworthy ones; §9.7's are superseded
for absolute-ratio purposes (its scaling *shape* finding still stands,
independently reconfirmed here).

**Still open:** whether a properly *validated* looser-than-default
config (unlike §9.4's) exists that's both fast and safe across kernel
shapes/ℓ ranges — worth exploring later with per-pair accuracy checks
before trusting any number from it, not just a single easy case. Also
still open: block-level OpenMP parallelism's effect on these ratios (all
of the above is single-threaded; §7 measured up to 2.58× on 8 workers for
a different, non-CCL benchmark) and N beyond 10.


### 9.9 §9.8 redone: the ratios were stale and the accuracy column was ℓ=2

§9.8's table was re-measured with the same script, same setup and the same
`OMP_NUM_THREADS=1`. CCL's side is unchanged (117 ms vs 119.0 ms at N=5,
434 ms vs 431.6 ms at N=10), so the comparison basis held and the whole
difference is on NumCosmo's side:

| N | §9.8 NumCosmo | now | §9.8 ratio | now (cold) |
|---|---|---|---|---|
| 5 | 1763.7 ms | **497.5 ms** | 14.83× | **4.28×** |
| 10 | 3272.8 ms | **863.8 ms** | 7.58× | **1.99×** |

**~3.5× at N=5 and ~3.8× at N=10**, from work landed after §9.8 was written:
FFTW wisdom I/O caching, the spectral/abstol convergence work, and the
per-block integrators (one `NcmSBesselIntegratorLevin` per ℓ-block, pinned to
its range and reused across `solve()` calls, with the registered kernels
evaluated directly instead of duplicated per thread).

§9.8's `O(N)` vs `O(N²)` shape still holds: the ratio roughly halves as N
doubles (4.28 → 1.99).

**Threads, which §9.8 left open.** All of §9.8 was single-threaded:

| | N=5 cold | N=5 warm | N=10 cold | N=10 warm |
|---|---|---|---|---|
| 1 thread | 4.26× | 5.77× | 1.97× | 2.91× |
| 12 threads | 1.66× | 2.24× | **0.82×** | 1.16× |

At N=10 on 12 threads NumCosmo computes a cold batch faster than CCL. The warm
columns are *less* favourable than the cold ones: CCL gains more from repetition
than we do (118 → 72 ms at N=5), so a cold-only comparison, which is all §9.8
measured, flatters us. Warm is the number that describes an MCMC.

**This benchmark cannot show threading, and that is a property of its ℓ range.**
ℓ=[2,25] at block size 8 is **3 blocks**, and blocks are the unit of parallelism
(§6.3), so nothing here can use more than 3 threads no matter what is set. Going
from 1 to 12 threads moves N=5 warm only 414 → 185 ms, and 3 → 12 threads not at
all. Read the ratios above as a 3-thread result.

Measured separately on a workload with enough blocks to spread (N=5 galaxy
kernels, ℓ=[2,200] — 25 blocks, warm `solve()`):

| threads | 1 | 2 | 3 | 4 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|---|
| warm solve | 1.74 s | 0.92 s | 0.65 s | 0.52 s | 0.37 s | 0.31 s | **0.25 s** |
| speedup | 1.00× | 1.89× | 2.68× | 3.35× | 4.70× | 5.61× | **6.96×** |

Clean scaling to 12 threads, no plateau — 58% efficiency at 12, consistent with
25 blocks quantising into 3 waves. §6.3's block-parallel axis holds.

**Watch `OMP_THREAD_LIMIT` when measuring.** This machine's default environment
(`set_threads.sh 3`) exports `OMP_THREAD_LIMIT=3` alongside `OMP_NUM_THREADS=3`,
and `OMP_THREAD_LIMIT` caps the whole program regardless of what
`OMP_NUM_THREADS` is set to per-run. An earlier draft of this section reported a
"plateau past 4 threads" that was entirely this: every run labelled 4, 8 or 12
threads was executing on 3. Override both, or re-source the script with the
count you want.

**A production-width ℓ range, for the threading the ℓ=[2,25] case cannot show**
(N=5, same weak-lensing kernels, `--lmax=200` — 25 blocks):

| | NumCosmo cold | NumCosmo warm | CCL cold | CCL warm | ratio cold | ratio warm |
|---|---|---|---|---|---|---|
| 1 thread | 10370.8 ms | 8850.9 ms | 964.3 ms | 587.7 ms | 10.75× | 15.06× |
| 12 threads | 2616.3 ms | 2156.7 ms | 1031.4 ms | 628.9 ms | 2.54× | 3.43× |

Threading now does what §6.3 predicted — 4.10× on 12 threads (8850.9 → 2156.7
ms). But the ratio against CCL is **worse** at ℓ=[2,200] than at ℓ=[2,25]
(15.06× vs 5.77× single-threaded): our cost per ℓ grows faster than CCL's.
Going 25 → 200 multipoles is 8× the ℓ's for 11.6× our time and 8.7× CCL's. That
is the scaling axis to attack next, and it is not visible at all in the ℓ=[2,25]
benchmark §9.8 used.

**There is no ℓ ceiling — corrected.** An earlier version of this section claimed
forced non-Limber broke down above ℓ ≈ 200, because `--lmax=210` aborts in
`_ncm_spectral_compute_chebyshev_coeffs_adaptive_internal`. That diagnosis was
wrong. Measured directly:

- `numcosmo xcor kernel view --kernel "weak-lensing survey=LSST-Y1 bin_idx=0"
  --ell=990` works at default precision.
- With this benchmark's *own* narrow Gaussian kernel, a single-kernel closure
  build succeeds at ℓ = [2,9], [194,201], [202,209], [210,217], [500,507] and
  [990,997].

The failure is not a function of ℓ alone. Isolating through `NcXcorSolver` at
`lmax=210` first looked like a **kernel-count threshold** (1 OK, 2 OK, 3+ fail),
but that was an artifact of how the reproducer built its list: the kernels are
μ = 0.2, 0.3, 0.4, 0.5, 0.6 taken in order, so "3 kernels" is the first list
that contains μ = 0.4. Registering each kernel *alone* through the solver:

| μ | 0.2 | 0.3 | **0.4** | 0.5 | 0.6 |
|---|---|---|---|---|---|
| result | OK | OK | **fails** | OK | OK |

One kernel is enough. A finer μ scan at `lmax=210` is knife-edge and
non-monotonic — 0.34 OK, 0.36 fails, 0.38 OK, 0.39 fails, 0.40 fails, 0.41 OK —
and the same μ = 0.4 kernel is fine at `lmax` ≤ 200 and fails at 210, 250, 300,
400. Root cause and fix in §9.10.

Also ruled out by measurement, each independently: the dn/dz tail truncation
(±8σ, ±5σ, ±4σ, ±3σ all fail), the kernel width (σ = 0.08, 0.15, 0.25 all
fail), the block geometry (including a final block of `n_ell = 1`, which works),
and ℓ in isolation. Checking out the pre-per-block-integrator versions of the
four xcor files and rerunning reproduces the failure, so it predates that work
and is not a regression from it.

**Measurement note:** `numcosmo xcor kernel view` ends in `plt.show()`
(view.py:901), which blocks on an interactive backend. Set `MPLBACKEND=Agg`
when running it non-interactively, or it will appear to hang with no output
(stdout is block-buffered, so nothing reaches the log either).

**The accuracy column was one multipole.** §9.8's "0.3-0.6%" is
`max(|nc/ccl - 1|)` over ℓ=[2,25] and over all pairs, which ℓ=2 sets by
itself. Resolved per ℓ (N=5, default precision):

| ℓ | 2 | 3 | 5 | 8 | 15 | 25 |
|---|---|---|---|---|---|---|
| max dev | 3.0e-3 | 7.1e-4 | 5.0e-4 | 4.2e-4 | 3.3e-4 | 3.0e-4 |

| | max | median |
|---|---|---|
| ℓ < 8 | 3.0e-3 | 3.2e-4 |
| ℓ ≥ 8 | **4.2e-4** | **1.6e-4** |

So agreement away from the lowest multipoles is an order of magnitude better
than §9.8 reported. The 0.3% figure is a real bound, but it describes ℓ=2, not
the method. The bench script now prints the per-ℓ breakdown and the cold/warm
split so this cannot recur.

**§9.4/§9.8's tuned config is now refused in code.** §9.8 concluded
`reltol=1e-2`/`cheb-reltol=1e-2` "should not be used for anything where the
answer matters". `_nc_xcor_check_kernel_tolerance()` now fails loudly when
`NcXcor:reltol` asks for more precision than the kernel's closure carries,
which is exactly that configuration, so `--tuned` no longer runs.

### 9.10 §9.9's abort: root cause, fix, and a second (unfixed) defect

Found by running the failing case under gdb rather than by reasoning about it.
The abort is one specific panel, and the backtrace names every actor:

```
_ncm_spectral_compute_chebyshev_coeffs_adaptive_internal (a=99.999999999999957,
    b=101.13342580488548, tol=1e-08, abstol=5.4259111550197567e-31)
_ncm_sbessel_integrator_levin_integrate_panel (a_p_idx=12, b_p_idx=-1, k=127.40274063204693,
    ell_min=210, ell_max=210)
_component_states_compute_non_limber (k=127.40274063204693)
ncm_function_sample_set_expand_domain (x_min_hard=0.00428, x_max_hard=4282749.4)
```

So: the **last** panel, `[last_knot, y_max] = [100, 101.133]` in y-space, for the
single-multipole block ℓ = [210,210], at k = 127.4 — a k the *domain-expansion*
phase is probing, not one the answer depends on. Two things are wrong there.

**(a) The abstol bound throws away 45 orders of magnitude.** The panel's
contribution is `b_p·j_ℓ(b_p)·u'(b_p) − a_p·j_ℓ(a_p)·u'(a_p)`, and
`_ncm_sbessel_integrator_levin_panel_abstol()` converted a bound on the result
into a bound on the Chebyshev coefficients by dividing by `b_p` — i.e. using
`|j_ℓ| ≤ 1`. But y = 100 is far *below* the ℓ = 210 turning point:
`j_210(101.13) = 1.4e-47`. The panel provably cannot move the accumulated
result (5.5e-13) no matter what the RHS looks like, yet it was being asked for a
1e-8 *relative* Chebyshev fit of the integrand, and refined to N = 65537 trying
to get one.

Fixed: the scale is now `max over the batch of max(b_p·|j_ℓ(b_p)|,
a_p·|j_ℓ(a_p)|)` — the same inequality, not bounded away. Since `|j_ℓ| ≤ 1` that
scale is never larger than `b_p`, so the floor is **never tighter** than before:
the change can only skip work that provably cannot matter, never demand less
accuracy where it can. The computation also moved into
`_ncm_sbessel_integrator_levin_solve_and_accumulate()`, the one place that has
the `j_ℓ` arrays; that removes three duplicated call sites and incidentally
covers the no-knots single-panel branch, which never set a floor at all.

Results:

- every previously-aborting case passes — μ ∈ {0.36, 0.39, 0.40} at lmax 210 and
  400, a 15-point μ sweep from 0.15 to 0.90 at lmax = 260, and the original
  reproducer at all of N = 1…5 kernels;
- answers unchanged: over 15 pairs × 199 multipoles the **median** deviation
  from the pre-fix code is exactly 0 and the max is 3.2e-12, four orders below
  the integrator's own 1e-8 reltol;
- agreement with CCL is unchanged (ℓ ≥ 8: max 4.2e-4, median 1.624e-4 — the same
  figures as §9.9);
- ~6% faster, consistently: 8.29 → 7.79 s at lmax = 200 and 17.28 → 16.24 s at
  lmax = 300 (best of 3, `OMP_NUM_THREADS=1`, N = 5, 15 pairs).

**(b) `nc_xcor_lensing_efficiency_eval()` extrapolated below its first knot.**
This is *why* the panel was not smooth. Sampling the weak-lensing component's
own integrand toward `xi_max`:

| 1 − ξ/ξ_max | 1e-6 | 1e-7 | 1e-8 | 1e-10 | 1e-12 | 0 |
|---|---|---|---|---|---|---|
| before | 9.8e-28 | 9.8e-30 | **6.1e-24** | **1.3e-20** | **1.3e-20** | **1.3e-20** |
| after | 9.8e-28 | 9.8e-30 | 9.8e-32 | 9.8e-36 | 9.8e-40 | **0** |

F vanishes quadratically, as it should — and then turned around and climbed
back to a floor of 1.3e-20, a 40× jump at the endpoint. `ncm_spline_eval()` on
a cubic below the first knot takes `ncm_spline_get_index() == 0` and evaluates
the first cubic with a *negative* `delx`: silent extrapolation, no clamp and no
error. `nc_xcor_lensing_efficiency_prepare()` starts CVODE at
`mz_ini = -zmax + dz`, `dz = zmax·sqrt(DBL_EPSILON) = 1.55e-8`, so the last `dz`
of the range had no spline over it. The noise region matched `dz` exactly, and
a Chebyshev series cannot converge relatively across it.

**The `dz` offset is necessary — a first attempt to remove it was wrong.** The
ODE *is* regular at `mz = -zmax`: `ydot0 = -y1/E = 0` and
`ydot1 = -W_src/d_T ≠ 0`, so `y = (0,0)` is a legitimate initial condition and
the offset IC is just its second-order Taylor expansion. Starting there anyway
fails, for a different reason — CVODE's error weight is
`1/(reltol·|y| + abstol)`, and with `y = (0,0)` and `abstol = 1e-50` that is
`1e50`, an effectively infinite accuracy demand:

```
[CVode] Internal t = -1.04 and h = 5.38761266744418e-26 are such that
        t + h = t on the next step. The solver will continue anyway.
```

The offset exists to give the state a nonzero scale, not to step away from a
singularity.

So the repair goes in the evaluator instead. `prepare()` now also stores
`head_zmax`, `head_dz` and `head_c2 = f/dz²`, and
`nc_xcor_lensing_efficiency_eval()` uses the second-order expansion the initial
condition was already built from on `[zmax − dz, zmax]`, returning exactly 0 at
and above `zmax`. It agrees with the spline at the first knot by construction —
same `f` — and the `f` expression is left bit-identical so the ODE start is
untouched. Inside the spline domain the sampled values are unchanged to every
digit; the Cl impact is pure roundoff (max 8.5e-16, median 0), full suite 86/86
and the pytest lane 893 passed.

**The two fixes are independently sufficient, not jointly necessary.** Found
while trying to build a regression test, by rebuilding all four combinations
and rerunning the μ = 0.4 / lmax = 210 reproducer:

| | lens_eff old | lens_eff new |
|---|---|---|
| **levin old** | **aborts** | OK |
| **levin new** | OK | OK |

So (b) is the root cause and (a) is independent robustness — either one alone
removes this failure. The practical consequence is that *this* failure cannot
serve as a regression test for (a): with (b) in place the trigger is gone, and
such a test would pass with (a) reverted.

Each fix is therefore tested through its own lever, both in ~0.5 s:

- `TestPanelAbstolEvanescent`
  (`tests/python/ncm/specfunc/test_sbessel_integrator_levin.py`): a deliberately
  unresolvable integrand over a range entirely below the ℓ-th turning point,
  with a caller-supplied `abstol` — the exact quantity (a) rescales. Core dumps
  without the fix, returns in 11 ms with it. Paired with the opposite case,
  ℓ = 10 with the turning point *inside* the range, checked against
  `scipy.quad` to rtol 1e-8, so the floor is shown to be a relaxation and not a
  licence to under-resolve.
- `test_weak_lensing_kernel_decays_to_zero_at_xi_max`
  (`tests/python/nc/xcor/test_kernel_component.py`): samples the WL component
  across the first spline knot on all five LSST Y1 source bins. Without (b):
  `kernel grows toward xi_max: 1.356403e-23 at 1-xi/xi_max=1.0e-10 exceeds
  2.906403e-34 at 1.0e-08`.

Three earlier candidate tests were discarded because they **passed on the
pre-fix build**: a `max()` floor (a kink converges algebraically, so the
level-doubling test clears it), the reproducer shrunk to `xc.compute` on the
single failing block (wrong code path), and the reproducer with a narrowed ℓ
range (which is what exposed the 2×2 above). Always rebuild the pre-fix code
and confirm the candidate fails before keeping it.
