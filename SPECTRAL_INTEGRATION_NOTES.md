# Spectral (Chebyshev) integration for NcDataClusterWL — assessment & plan

Notes by Claude, 2026-06-17, at Sandro's request. Context: we just finished hardening
the cluster-WL redshift integration (FIXED_NODES / LNINT / CUBATURE), found and fixed
the kink-at-z_cl + effective-support bug, and built a deterministic truth-table test.
Sandro prepared `NcmSpectral` (Chebyshev toolkit) and `tests/python/math/test_spectral_weighted.py`
and asks whether to add a *fourth* integration method based on it.

This file is an untracked working note — move/keep/delete as you like.

## 1. What NcmSpectral gives us (from the test + header)

- `compute_chebyshev_coeffs_adaptive(F,a,b,k_min,tol)` → standard Chebyshev coeffs of F,
  refining k until the tail is below `tol` (N = 2^k + 1).
- `compute_chebyshev_coeffs_adaptive_weighted(...)` → coeffs of F(x(t))·√(1−t²); the key
  property tested is **a₀·π = ∫_a^b F dx** (a Clenshaw–Curtis-style quadrature read off the
  first weighted coefficient).
- **Product integral** (TestSpectralProductIntegral): with `a_k` = *weighted* coeffs of F and
  `b_k` = *standard* coeffs of G,
  ∫_a^b F·G dx = π·a₀·b₀ + (π/2)·Σ_{k>0} a_k·b_k.
- Plus `chebyshev_eval/deriv`, and projection / multiplication-by-x / **differentiation**
  matrices (`get_d_matrix`, …) — i.e. a real spectral-operator toolkit, not just a quadrature.

Exponential convergence for smooth integrands; algebraic for ones with limited smoothness
(e.g. a C² cubic spline) or a kink.

## 2. How it maps onto the WL per-galaxy integral

Per galaxy the likelihood is
  P_gal = ∫ dz  P(z | z_obs) · S(z; cosmo, M, profile)
where **P(z) is data** (gauss kernel or pz cubic spline — fixed once the catalog is loaded)
and **S(z) is the model-dependent** reduced-shear shape factor (changes every model eval /
fit step). This is *exactly* the F·G product structure: F = P (data), G = S (model).

So a spectral method would:
- precompute `a_k` = weighted Chebyshev coeffs of P(z) **once per galaxy at prepare**;
- per model eval, compute `b_k` = standard coeffs of S(z) (evaluate S at Chebyshev nodes);
- combine with the dot-product formula above.

## 3. The honest catch: FIXED already exploits the same separation

`_nc_data_cluster_wl_eval_m2lnP_fixed` + `make_fixed_nodes` already:
- bake the per-galaxy P(z) into the fixed-node **weights at prepare** (data part, amortized),
- recompute only S(z) via `nc_galaxy_sd_shape_eval_at_nodes` per eval (model part),
- and (since this PR) split at z_cl so no panel straddles the reduced-shear kink.

And we measured FIXED at 20 nodes to be accurate to ~1.4e-8 (gauss) / ~2e-6 (pz) vs the
320-node truth, at trivial cost. So the spectral product method's *structural* advantage
(amortize data, cheap model combine) is **already captured** by FIXED. For the smooth shear
factor S(z), 20 Gauss–Legendre nodes are effectively exact, so Chebyshev's exponential
convergence would not cut the node count materially. The pz spline P(z) is the only
algebraically-converging piece, and it is precomputed in *both* schemes.

Two caveats that do **not** go away with spectral:
- **The kink at z_cl** breaks smoothness → spectral would still need the same panel split at
  z_cl (otherwise adaptive refinement explodes), exactly like FIXED/LNINT/CUBATURE.
- **pz P(z)** is only C² → algebraic convergence of `a_k`; fine (precomputed once) but it is
  not where spectral shines.

## 4. Where spectral *could* genuinely win (worth a benchmark, not assumed)

1. **Analytic fit gradients.** The MC/fit path uses `NCM_FIT_GRAD_NUMDIFF_FORWARD` (finite
   differences → N_param extra full likelihood evals per gradient). The spectral
   differentiation matrices let you represent P_gal(θ) and differentiate in spectral space,
   potentially giving **analytic ∂P_gal/∂θ** cheaply. That would speed up *fits* (the
   expensive monte_carlo test, and real MCMC) far more than shaving integration nodes — this
   is the most promising angle by far, and it is orthogonal to "which quadrature."
2. **A regime where S(z) is expensive or wiggly** (e.g. richer shear models, magnification,
   non-trivial reduced-shear calibration) where exponential convergence in few S-evaluations
   beats a fixed node count. Not the case for the current smooth NFW shear, but could matter
   for future models.
3. **Representing P_gal as a smooth function of a parameter** (e.g. mass) to reuse across a
   grid — niche.

## 5. Recommendation

- **Do not add a 4th production `NcDataClusterWLIntegMethod` just for the integral itself.**
  FIXED is now validated, exact-enough, and already does the data/model split; a spectral
  *quadrature* would add a method to maintain and validate (incl. a 4th truth-table column)
  for little measurable gain on the current integrand.
- **Do pursue the spectral toolkit for the gradient angle** as a separate investigation:
  prototype `∂(−2lnP)/∂θ` via NcmSpectral and benchmark a fit against the current
  numdiff Nelder-Mead / numdiff-gradient path. If it gives a clear fit speedup (it plausibly
  will), that is the high-value use of the work — and it is independent of the integration
  method enum.
- Keep `NcmSpectral` + `test_spectral_weighted.py` as a first-class math utility regardless;
  it is well-tested and useful on its own (and for #4.1/#4.2).

## 6. If we *do* add it as an integration method later — integration checklist

- Add `NC_DATA_CLUSTER_WL_INTEG_METHOD_SPECTRAL` to the enum; dispatch in
  `_nc_data_cluster_wl_eval_m2lnP` / `_eval_m2lnP_gal`.
- Precompute weighted `a_k` of P(z) per galaxy in `_nc_data_cluster_wl_prepare`, **split at
  z_cl** into two panels exactly as the fixed-node grid is (reuse `fixed_nodes_zcl` logic).
- Per eval, compute standard `b_k` of S(z) (needs `eval_at_nodes` at Chebyshev abscissas —
  a new shape entry point, or reuse with Chebyshev z-nodes), combine via the product formula,
  recombine the two panels in the linear domain (the `_nc_data_cluster_wl_integ_combine`
  helper already does logsumexp panel combination).
- Validate: it drops straight into `test_nc_data_cluster_wl_truth.c` as a 4th method vs the
  FIXED-320 golden, at a tolerance set from its measured floor (same pattern as the others).
- Beware: the product formula returns the *linear* integral; for the lnint-style summation
  you'd still take −2 ln of it. Make sure the weighted-coeff normalization (the `·π` factor)
  and the [a,b]→[−1,1] scaling are handled per panel.

## 7. One-line verdict

The spectral product structure is elegant and a genuinely good fit for the data×model form of
the integral — but FIXED already captures that structure and is exact and cheap, so the
payoff is not in replacing the quadrature. The real opportunity is **spectral derivatives for
analytic fit gradients**; chase that, benchmark it, and only then decide on a 4th method.
