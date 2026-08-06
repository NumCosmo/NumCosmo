# WL Shape-Factor / Ellipticity: Design History

This is the consolidated development history for the
`nc_galaxy_shape_factor_*` classes, `nc_wl_ellipticity`, and
`nc_galaxy_shape_intrinsic_mode` — bugs found and fixed, approaches tried
and rejected, and the verification evidence behind non-obvious design
choices. The technical docs (gtk-doc comments, the `wl_shape_marginalization_*.qmd`
pages) state current facts and point here for the "why" behind them; this
file is not rendered into the Quarto site and is not meant to be read
top-to-bottom — it's a reference for "why does the code do X" questions.

## `NcGalaxyShapeFactorSeries` (noise-side) and `NcGalaxyShapeFactorKnots`: removed

Both classes were removed from the codebase (kept: `Quad`, `FixedQuad`,
`SeriesLensed`, `VarAdd`, `Laplace`).

- `NcGalaxyShapeFactorSeries` expanded the *noise kernel* in powers of $g$.
  Its Taylor coefficients scale like $1/\sigma_\mathrm{noise}^{2n}$ and
  combinatorially blow up whenever a galaxy's per-object `std_noise`
  (observational data, unbounded, no prior) is small relative to $g$ — not
  roundoff (reproduced with 50-digit `mpmath`, still wrong at high order),
  a genuine slow-convergence catastrophe, the same cancellation mechanism as
  reconstructing $e^{-x}$ from its own Taylor series.
  `NcGalaxyShapeFactorSeriesLensed` (kept) sidesteps this by expanding the
  *population* term instead, whose coefficients scale with
  $1/\sigma_\mathrm{pop}^2$ — a population/prior parameter this project
  already hard-constrains to $(0.2,0.4)$, never per-galaxy data.
  `NcGalaxyShapeFactorFixedQuad` (kept) avoids the whole failure family
  structurally, since it has no series in $g$ at all.
- `NcGalaxyShapeFactorKnots` did fixed-node quadrature in $\chi_I$-space via
  a $(u,v)$-plane compactification, requiring the box to be re-centered on
  every evaluation because the integrand's peak moves with $g$. It was
  explicitly documented as "not a recommended likelihood backend... meant to
  be poked at", with no accuracy validation. `NcGalaxyShapeFactorFixedQuad`
  is the validated version of the same underlying idea (fixed-node
  quadrature), with a $g$-independent domain that needs no re-centering.

Removal (this pass): deleted `.c`/`.h`, build registrations, Python test
files; `numcosmo/ncm/algebra/ncm_laurent_series.c` (used by `SeriesLensed`)
and `numcosmo/nc/lss/galaxy/nc_galaxy_shape_intrinsic_mode.h` (used by
`Laplace`) only referenced the removed classes in doc comments, no
functional dependency.

## Noise-side `Series`: the arbitrary-order engine

Both ellipticity conventions' own map series
$f_g(\chi_I) = \chi_I + \sum_n F_n(\chi_I) g^n$ satisfy a short linear
recursion (verified against `mpmath` 50-digit Taylor coefficients of the
exact maps through order 8):

$$
\text{TRACE\_DET:}\quad F_n = -\chi_I F_{n-1}\ (n\ge2), \quad F_1 = 1-\chi_I^2.
$$
$$
\text{TRACE:}\quad F_n = -2sF_{n-1} - F_{n-2}\ (n\ge3,\ s=\mathrm{Re}\,\chi_I),
\quad F_0=\chi_I,\ F_1=2-2s\chi_I,\ F_2=\bar\chi_I-2sF_1-F_0.
$$

This recursion is what made an arbitrary truncation order practical in the
shipped `NcGalaxyShapeFactorSeries`: $\Delta$ and $\exp\Delta$ were built
generically from Laurent-in-$\theta$ coefficient arrays (same
add/scale/convolve/conjugate operations for both conventions, only the
$F_n$ recursion differing), with $\exp\Delta$'s coefficients from the
standard formal-power-series-exponentiation recursion,
$c_m = \tfrac1m\sum_{k=1}^m k\,\tilde\delta_k\,c_{m-k}$
($\tilde\delta_k=\delta_k/2\sigma_\mathrm{noise}^2$, $c_0=1$) — see
`verify_Fn_recursions`, `F_list_eps`/`F_list_chi`, `marginal_generic` in
`dev-notes/wl_shape_series_marginalization_derivation.py`.

Higher truncation order helped substantially, for both conventions
(disc-integrated, $\sigma_\mathrm{noise}=0.30$, $\sigma_\mathrm{pop}=0.25$):

| $g$ | $\epsilon$, $N{=}4$ | $\epsilon$, $N{=}8$ | $\epsilon$, $N{=}20$ | $\chi$, $N{=}4$ | $\chi$, $N{=}8$ | $\chi$, $N{=}20$ |
|---:|---:|---:|---:|---:|---:|---:|
| 0.3 | $6.0\times10^{-3}$ | $3.7\times10^{-6}$ | $3.0\times10^{-15}$ | $3.1\times10^{-1}$ | $2.5\times10^{-2}$ | $1.0\times10^{-8}$ |
| 0.5 | $9.2\times10^{-2}$ | $4.5\times10^{-4}$ | $1.9\times10^{-10}$ | $6.3\times10^{0}$ | $5.0\times10^{0}$ | $6.3\times10^{-3}$ |
| 0.8 | $2.5\times10^{0}$ | $9.0\times10^{-2}$ | $1.3\times10^{-5}$ | $1.3\times10^{2}$ | $9.1\times10^{2}$ | (not converged by $N{=}20$) |

$\epsilon$ converged cleanly and fast — usably accurate even at $g=0.8$ by
$N\sim20$. $\chi$ converged too (poles at exactly $|g|=1$, the same radius
as $\epsilon$'s), but far more slowly, needing roughly $2$–$3\times$ the
order for comparable accuracy at $g\lesssim0.5$; at $g=0.8$ it had not
converged by $N=20$ in double precision, with non-monotonic intermediate
orders — a genuinely convergent but very slowly-convergent series being
evaluated where floating-point roundoff competes with the terms' own size,
not a settled result. Not chased further since it was far outside this
project's actual $g$ range.

The shipped class dispatched both conventions once at construction, had a
`trunc-order` construct property (default 4) and an `n-nodes` property for
the fixed Gauss-Legendre $\rho$-quadrature (default 60), windowed to
$\rho\in[\max(0,R-8\sigma_\mathrm{noise}),\min(1,R+8\sigma_\mathrm{noise})]$
around the noise envelope's peak at $\rho=R$. Per-order radial integrals
were cached per galaxy (depend only on $(R,\sigma_\mathrm{noise})$ and the
population, never on $g$), turning a whole `NcDataClusterWLFactor`
$z$-integral's cost from *(per-node full radial integral)* into *(one
radial integral per galaxy) + (per-node cheap polynomial evaluation)*.
Verified against the Python reference to $\sim10^{-12}$ for both
conventions and several orders before removal.

## `NcGalaxyShapePopBeta`: alpha's floor loosened below 1

$\alpha$'s floor was originally $\ge 1$, matching $\beta$'s (both avoid a
divergence, at $x=0$ resp. $x=1$). Loosened to $\ge 0.5001$ after a real
fit's posterior concentrated its mass right at $\alpha=1$ under
`FixedQuad`: a hard floor sitting exactly where the data want the posterior
to be would truncate/bias the inferred $\alpha$ one-sidedly, worse than
allowing the (`FixedQuad`-safe) divergent regime the data are actually
pointing to. `SeriesLensed` still needs $\alpha\ge1$ in practice — its
Taylor expansion's radius of convergence shrinks below that, since the
population stops being analytic at $x=0$ — but that is a `SeriesLensed`
limitation, not a reason to keep the class-wide floor at 1.

## `NcGalaxyShapePopBeta`: reparametrized to $r=|\chi_I|$ (2026)

Switched the model's own variable from $x=|\chi_I|^2\sim\mathrm{Beta}(\alpha,\beta)$
to $r=|\chi_I|\sim\mathrm{Beta}(\alpha,\beta)$ directly. Motivation: nearly
every weak-lensing paper reporting/fitting an ellipticity distribution does
so against $|\chi_I|$, not its square, so this way $\alpha$/$\beta$ are
directly comparable to those fits and $\mathrm{mode}(r)=(\alpha-1)/(\alpha+\beta-2)$
is literally the "peak ellipticity" number such papers quote — and a real
fit's posterior wanting $\alpha_\mathrm{old}<1$ (see the section above) was
itself a symptom of fitting the wrong (squared) coordinate: under the new
parametrization that same population lands at
$\alpha_\mathrm{new}\approx2\alpha_\mathrm{old}\ge1$ (this class's own
default moved $0.7\to1.4$ accordingly), inside $r$-space's ordinary,
non-divergent region.

`eval_p(x)` still reports the density w.r.t. $x$ (the base class's vfunc
contract), composing $P(r)$ with the $r=\sqrt{x}$ Jacobian:
$P(x)=x^{\alpha/2-1}(1-\sqrt{x})^{\beta-1}/[2B(\alpha,\beta)]$. Two
consequences worth recording, both found only by testing against a
brute-force/numerical cross-check rather than trusting the derivation on
its face (see `verify-closed-form-derivations-numerically` house rule):

- **The $x$-divergence threshold moved from $\alpha<1$ to $\alpha<2$.**
  $P(r)$ itself is non-divergent for $\alpha\ge1$ (why the floor could go
  back to the clean, symmetric $\alpha,\beta\ge1$, no more loosening below
  1), but the extra $x^{-1/2}$ Jacobian factor means $P(x)$ — what
  `eval_p()`/every numerical consumer actually integrates — still has a
  genuine pole at $x=0$ for any $\alpha<2$. This is not a wider bug surface
  than before (the *same* previously-divergent populations, e.g. old
  $\alpha=0.7$, are still exactly the ones divergent now, at
  $\alpha_\mathrm{new}=1.4$) — it only looks wider if the raw threshold
  numbers are compared without rescaling.
- **`SeriesLensed`'s usable range shrank from $\alpha\ge1$ (practical
  guidance) to $\alpha\ge2$ (much stronger).** `eval_p_rho2_g_series()`
  needs $\sqrt{x(g)}$ unconditionally now, even at integer $\alpha$ where
  the old $x^{\alpha-1}$ term was an entire function with no singularity at
  all (e.g. $\alpha=4$ used to compose to arbitrary order with no issue).
  $\sqrt{\cdot}$'s branch point at $x=0$ shrinks the $g$-Taylor series'
  radius of convergence to wherever $x(g)$ first reaches 0 in the complex
  $g$-plane, which empirically collapses to unusably small for any
  $\alpha<2$ — increasing `trunc_order` makes a previously-passing
  regression test (`alpha=beta=3`, chosen specifically as a "safe, no
  branch point" case pre-reparametrization) diverge by 16 orders of
  magnitude rather than converge, the classic non-roundoff
  outside-the-radius-of-convergence signature. This class's own default
  ($\alpha=1.4$) is now *inside* the broken region for `SeriesLensed` —
  `FixedQuad`/`Quad` have no such restriction (no series in $g$, or an
  entire-function exponent respectively) and remain the ones to trust for
  $\alpha<2$ Beta populations with this shape factor.

## `NcGalaxyShapeFactorQuad`: known accuracy bug for alpha<2 Beta populations

For a Beta population with $\alpha<2$ ($P(x)$ diverges at $x=0$ — see the
reparametrization entry above for why the threshold is 2, not 1), `Quad`'s
adaptive Divonne cubature loses accuracy against an independent scipy
oracle in a $g\sim[0.14,0.19]$ window, while `FixedQuad` stays accurate
throughout at the same configuration. Suspected cause: the singularity
isn't resolved by Divonne's adaptive subdivision of the generic box `Quad`
integrates over, unlike `FixedQuad`'s fixed lens-domain nodes. Not yet
fixed; `test_marginal_alpha_below_one_known_accuracy_bug`
(`test_galaxy_shape_factor_quad.py`) and
`test_marginal_matches_scipy_truth_table_beta_alpha_below_one_g_scan`
(`test_galaxy_shape_factor_fixed_quad.py`) pin the current behavior against
regression. `FixedQuad` is the one to trust for this population today.

## `NcGalaxyShapeFactorQuad`: rejected implementations

Before settling on Divonne cubature over the lensed frame $\chi_L$ with
explicit peak hints, two earlier approaches were tried and rejected:

1. Evaluating $P_\mathrm{pop}$ directly at $\chi_I$ (not $\chi_L$) via
   h-adaptive `NcmIntegralND` cubature.
2. Cuba's Cuhre (fixed-degree base rule, no explicit peak hints) with a box
   scaled to $\sigma_\mathrm{noise}$.

Both silently returned confidently wrong results — off by orders of
magnitude, at tight `reltol` — whenever the population was narrow relative
to the integration box and off-center: neither has a mechanism to notice an
isolated feature it never samples near. Divonne's explicit peak hints solve
this directly by telling the search where to look.

A related mistake: the second peak hint ($f_g(0)$, the population's peak
location under the disc) was briefly hardcoded as $g$. This is only correct
in the TRACE_DET convention ($f_g(0)=g$); the TRACE convention's map gives
$2g/(1+|g|^2)$, and the hardcoded version silently missed the peak for
narrow, off-center TRACE-convention populations. Fixed by computing $f_g(0)$
via the actual forward map.

## `NcGalaxyShapeFactorFixedQuad`: the "no overlap" boundary bug

The originally-shipped domain construction had a fourth branch, "no
overlap" ($d\ge R_1+R_2$, i.e. the fixed $n_\sigma$ noise-disk window misses
the unit disc entirely), returning a fixed `1e-300` floor regardless of how
close $d$ was to that threshold. This produced an arbitrary, discontinuous
jump in the marginal: at `std_noise=0.03`, `R=1.239999 → R=1.240000` (an
infinitesimal change) jumped `-2lnL` by **over 1000** for that one galaxy —
a real risk for any MC-style analysis where galaxies in different
realizations land on either side of this arbitrary, unphysical threshold
(the true noise kernel has infinite support).

**First fix attempt, wrong**: fall back to full-disc quadrature
(`_regen_unit_contained`) near that boundary. Looked fine at
`std_noise=0.03` (checked, a few percent error) but failed badly at smaller
`std_noise`: off by **~2350×** at `std_noise=0.01`, **~10^14** at
`std_noise=0.005`, at the same node count that looked fine at 0.03. The
surviving probability there is a deep `exp(-x^2)` tail concentrated in a
narrow sliver near the disc edge closest to the observation — narrower the
smaller `std_noise` is — and a full-disc grid with a modest, uniformly-spaced
node count has no reason to sample near an exponentially-localized peak it
was never aimed at.

**Lesson**: validating a fix at only one parameter value is not validation,
especially near a domain boundary or an asymptotic/tail regime. Sweep the
actual parameter range that matters before declaring a fix good.

**Correct fix**: don't switch quadrature schemes — grow the lens branch's
own effective noise-disk radius (`R2_eff`, used only for domain
construction, not the fixed window used for branch-selection checks) so it
always reaches at least `NSIGMA_TAIL` (=8) sigma of noise-kernel tail depth
into the unit disc: `R2_eff = max(R2, (d-R1) + NSIGMA_TAIL*std_noise)`. Since
`R1+R2_eff > d` whenever `std_noise>0`, the old threshold is never reached,
so the separate "no overlap" branch was removed entirely (three branches,
not four). Guard: if this growth would itself make `R2_eff` fully contain
the unit disc, fall through to the full-disc branch instead. Verified
smooth and accurate (sub-percent to a few percent, down to `~1e-130`-scale
deep-tail values) against scipy across the actual `std_noise` range this
project uses (0.03–0.3); residual jump at the old boundary reduced to a few
tenths of a unit in `-2lnL`.

The genuine-lens branch's own parametrization degenerates at exactly
$d=R_1+R_2$ (both two-arc half-angles go to zero, collapsing the $(u,v)$
grid onto a single point regardless of node count — a parametrization
degeneracy, not a resolution one); the `R2_eff` fix keeps this threshold
permanently out of reach rather than trying to resolve through it.

## `NcGalaxyShapeFactorFixedQuad`: narrow-population limitation

A fixed grid cannot resolve a population much narrower than its node
spacing ($\sigma_\mathrm{pop}\lesssim0.05$, or a sharply concentrated Beta
population). Tried `ncm_spectral_compute_chebyshev_coeffs_adaptive_weighted`
(adaptive weighted Chebyshev quadrature, used elsewhere in the project) as a
possible fix: it fails identically to the plain fixed grid (converges
confidently to an answer wrong by dozens of orders of magnitude), because
it's *passive* refinement (doubles resolution once its own coefficients
haven't decayed enough), not *active* peak-seeking like Divonne's explicit
hints — if the initial grid never samples near an isolated narrow peak, it
sees "nothing here" and stops. Not pursued further; use `Quad` for
narrower/more exotic populations. Production only ever uses Gaussian
populations with $\sigma_\mathrm{pop}\in(0.2,0.4)$, well inside `FixedQuad`'s
validated regime.

## `NcGalaxyShapeFactorFixedQuad`/`SeriesLensed`: branch-cost cliff

At `std_noise` roughly in $[0.09, 0.2]$ (where the noise-disk radius
$R_2=8\sigma_\mathrm{noise}$ is close to the unit-disc radius $R_1=1$),
nearly every galaxy in a realistic population lands in `FixedQuad`'s
expensive genuine-lens branch (1681 nodes at the default `n-lens=41`)
instead of the cheap contained branches (~225–961 nodes) — measured
directly: 198/200 galaxies at `std_noise=0.12`. This project's actual
production regime, `std_noise~0.3`, is outside this range and is in fact
`FixedQuad`'s cheapest, cleanest regime.

An adaptive-`n_lens` calibration (minimal node count reaching a target
`reltol`, ported from the `pz_auto_nodes` branch's `ncm_integral_fixed_calibrate`
pattern) was prototyped in Python across 21 configs in the expensive middle
range: calibrated `n_lens` ranged 23–29 (529–841 nodes) vs. the fixed
default 41 (1681 nodes) — mean 2.70× fewer nodes. Prototype scripts:
`dev-notes/wl_fixed_quad_lens_nlens_calibration.py`,
`dev-notes/wl_fixed_quad_lens_domain_prototype.py`.

Shipped as the opt-in `auto-lens-nodes`/`lens-node-reltol` properties on
`NcGalaxyShapeFactorFixedQuad` (default off, so the existing `rtol=2e-4`
`n_lens=41` behavior is unchanged unless a caller opts in). A full
mass-recovery sweep with `auto-lens-nodes` enabled shows 0 crashed cells
across the previously-affected `std_noise` range, at 1.3-1.75× lower cost.

## `SeriesLensed`/`Series`: large-$g$ crash (fixed, commit `f3f6db7b`)

Found via a broad mass-recovery Monte-Carlo sweep (true cluster mass up to
$3\times10^{15}\,M_\odot$, `std_noise`$\in[0.03,0.3]$,
$N_\mathrm{gal}\in[3000,15000]$): both `SeriesLensed` and the (now-removed)
noise-side `Series` could fatally abort the whole fitting process — not
just return an inaccurate number — for a small tail of galaxies at high
enough cluster mass.

Root cause: $J(g)$ is evaluated by Horner's rule on cached, $g$-independent
coefficients — an exact degree-$N$ polynomial in $g$. Any nonconstant
polynomial satisfies $J(g)\to\pm\infty$ as $g\to\infty$; if $J(g)$ starts
positive at $g=0$ (it must) and trends downward as $g$ grows (the common
case), it is *guaranteed* to cross zero at some finite $g$ and diverge to
$-\infty$ beyond that. At the crash-triggering configuration
($R=0.3825$, $\sigma_\mathrm{noise}=0.03$, $\sigma_\mathrm{pop}=0.28$, $N=4$),
$J(g)$ crosses zero at $g\approx0.29$–$0.30$. Once $J(g)<0$,
`eval_ln_marginal`'s $\log J(g)$ is `NaN`, which reaches
`NcDataClusterWLFactor`'s adaptive redshift integral and aborts the whole
process (`ncm_integral1d_eval_lnint`'s `g_error()` path, uncatchable from
Python). Checked directly: the crossing is not an artifact of the
$\rho$-integration — a single isolated $(\rho,\theta)=(R,\phi)$ term of the
sum crosses zero at essentially the same $g$, since the noise envelope's
weight is concentrated right around $\rho\approx R$.

**A candidate fix that didn't pan out**: reparametrizing $g=\tanh(\zeta/2)$
(a "rapidity", linearizing shear composition the way it linearizes
relativistic velocity addition). Tested at a single isolated near-disc-edge
point ($\rho=0.99$): a real, order-of-magnitude improvement. But pushed
through the actual radial quadrature at the real crash configuration, it is
*slightly worse*, not better — the rapidity substitution's benefit lives at
the disc edge, but the integral's weight lives near $R$, which for this
(and most) configurations isn't there.

**The fix actually shipped**: guard $J(g)$'s sign (and NaN) directly —
`J(g) > 0 ? J(g) : 1e-300`, a comfortably positive floor ($\log\approx-691$,
never $-\infty$). A first attempt used a smooth window in $g$ instead
(worried a bare clamp's slope discontinuity would defeat the caller's
adaptive integrator), but that broke the test asserting a *higher*
truncation order stays accurate further out than a lower one — a fixed
$g$-window isn't order-aware, since the crossing point depends on both $N$
and $(R,\phi,\sigma_\mathrm{noise})$. The bare sign/NaN guard is
order-agnostic and was applied identically to both `SeriesLensed` and the
(now-removed) `Series`.

## `SeriesLensed` derivation: rejected approaches and a bug

- A first naive attempt at the population-side Taylor coefficients used a
  $g$-Taylor extraction via fixed-step central differences; turned out to
  be numerically unstable at high order for an unrelated reason ($h^n$
  cancellation, not a convergence failure of the true series — confirmed
  with `mpmath`'s adaptive-precision differentiation). Fixed by using the
  closed-form $\chi_I(\chi_L,g)$/Jacobian series instead (no finite
  differences anywhere).
- The guard on $J(g)$'s sign is checked against a fixed threshold in
  principle simpler than a $g$-magnitude cutoff — see the crash section
  above for why a fixed $g$-window doesn't work (order-dependent crossing
  point).

## Noise-side `Series` derivation: bugs found during derivation

- **Rotation-covariance bug**: the first version of the noise-kernel series
  (Step 2 in `docs/theory/wl_shape_marginalization_series.qmd`) assumed
  $f_g$ commutes with rotating $\chi_I$ (equivalently, rotating
  $\epsilon_\mathrm{obs}$ to lie on the real axis before expanding) while
  keeping $g$ real and fixed. It does not — full rotation covariance
  requires $g$ to rotate too. The mistake produced $O(g)$ truncation
  "errors" against the brute-force integral that did *not* shrink as
  $g\to0$ — the tell that something was structurally wrong, not just
  imprecise. Fixed by keeping $\epsilon_\mathrm{obs}$ as a general complex
  parameter throughout instead of rotating it away. Lesson: this is exactly
  the class of mistake that demands checking every closed-form derivation
  step against a brute-force numerical reference before trusting it.
- **Sympy extraction bug** (TRACE/$\chi$ convention derivation): a first
  attempt used `sp.fraction()` + `Poly(den,w).degree()` to auto-detect how
  much to clear denominators by, which silently dropped the
  $\sigma_\mathrm{noise}^2$ part of a mixed $(w,\sigma^2)$ denominator for
  $m\ge2$ — the $w$-structure of the printed coefficients looked plausible,
  but the $1/\sigma^{2m}$ factor was missing. Caught by an explicit
  cross-check ($c_1$ recomputed two ways), not by inspection. Fixed by
  clearing with an explicit $w^{2m}$ instead of auto-detecting.

## `nc_galaxy_shape_intrinsic_mode.c` (TRACE closed-form mode-finder): off-by-factor-of-2

The TRACE-convention closed-form theta-profile solves a quartic in
$t=\tan(\theta_B/4)$ (a Weierstrass-substitution variable). The actual angle
is $\theta_B/2$, not $\theta_B$ itself — i.e. the correct recovery is
$2\arctan(t)+\text{branch}\cdot\pi$, not $4\arctan(t)+\text{branch}\cdot2\pi$.
The wrong (factor-of-2) version caused a ~45% NaN-mismatch rate on the
first C implementation; caught by re-deriving from the verified Python
reference.

## `NcGalaxyShapeFactorFixedQuad`/shared header: performance work

Two rounds of profiling (`perf record`) on `FixedQuad`'s lens-branch hot
loop found and fixed three real costs, combined ~3.2× fewer perf samples
(26384 → 8164), bit-identical output throughout:

1. `exp(lndet_jac(...))`: the shared `nc_wl_ellipticity_lndet_jac_*_c()`
   computes the Jacobian via `log1p` then `FixedQuad` immediately
   exponentiated back to linear scale (it sums linearly, not in log-space)
   — roughly a quarter of runtime. Fixed by adding direct (non-log) variants
   `nc_wl_ellipticity_det_jac_{trace,trace_det}_c()` to the shared header
   (verified to match `exp(lndet_jac(...))` to `1e-9` relative precision).
2. `hypot`/`atan2` were recomputed on every `eval_marginal` call just to
   check cache validity, even though `eps_obs` is fixed per galaxy across a
   fit. Fixed by comparing the raw `epsilon_obs_1`/`epsilon_obs_2` doubles
   directly, computing `R`/`phi` only inside the actual domain-rebuild path.
3. `apply_shear_inv`'s `cabs(g)<=1.0` branch check (a `hypot`-style sqrt)
   runs every call since `g` changes every call and can't be cached like the
   Jacobian round-trip. Fixed by comparing `|g|^2<=1.0` instead (identical
   since `sqrt` is monotonic) — in the shared header, so this also benefited
   `Quad`, `VarAdd`, `Laplace`, and the since-removed `Knots`.

Separately, `nc_wl_ellipticity.c`'s naive complex division
(`_nc_wl_ellipticity_cdiv`, skipping the overflow/underflow guards libc's
`__divdc3` provides) was added because those guards dominated runtime in
`nc_galaxy_shape_intrinsic_mode.c`'s per-galaxy mode search, which calls
`apply_shear`/`apply_shear_inv` repeatedly inside a finite-difference Newton
search — unneeded here since $g$/$\chi$ values are always within or near the
unit disk.

## `NcGalaxyShapeFactor`: rotate the shear, not the observed ellipticity

The per-galaxy cache used to store `epsilon_obs_t`/`epsilon_obs_x` (the
observed ellipticity pre-rotated into the tangential/cross frame) and
`c1_rot`/`c2_rot` (the calibration bias, same rotation). Both were rebuilt
on every `marginal()` call because the rotation angle $\phi$ depends on the
cluster's `ra`/`dec` — a free parameter in a real fit — so `epsilon_obs`,
which should be a fixed per-galaxy datum, silently changed with every model
state. This forced `FixedQuad`'s marginal-spline cache (keyed on
`epsilon_obs`) to rebuild from scratch on every evaluation instead of once
per thread.

Fixed by rotating the reduced shear into the data frame instead: the
marginal is exactly rotation-covariant
($P(\epsilon_\mathrm{obs}e^{ia}\mid ge^{ia})=P(\epsilon_\mathrm{obs}\mid g)$),
so rotating `gt+bias` by $e^{2i\phi}$ gives the identical answer as the old
convention, while `data->epsilon_obs_1/2` never changes. Verified against
the frozen legacy-oracle tests (rtol=1e-12), which caught a bug in an
earlier draft: only `gt` (tangential-native) should rotate — the
calibration bias is already expressed in `data->coord` and must be added
unrotated.

## Mass-recovery sweep: representative numbers

An initial full sweep (`N_\mathrm{gal}\times\text{std\_noise}\times\text{true\_mass}\times\text{method}`),
run before `NcDataClusterWLFactor` had `FIXED_NODES`/`auto-nodes` support
and before the `auto-lens-nodes` branch-cost-cliff fix above, found 109/120
cells converged cleanly and 11 timed out (all `n_gal>=12000`,
`std_noise=0.09–0.12`, the branch-cost-cliff zone above), with `LNINT` as
the only available redshift-integration method. At `std_noise=0.3`
(production regime), averaged across `n_gal`: VarAdd 3.83s/-27.0% bias,
`Series` (removed) 6.08s/-9.8%, `SeriesLensed` 8.68s/-9.8%, `Laplace`
23.37s/-9.2%, `FixedQuad` 59.32s/-7.6% (best accuracy, but slowest of the
well-behaved methods at this specific regime — a nuance against
overselling `FixedQuad` purely on the `Quad`-comparison speedup headline,
since `Quad` itself is infeasible at this scale: a feasibility demo found
no detectable mass signal out of `Quad` at `std_noise=0.3` even at 3000
galaxies, ~57 minutes of compute).

`NcDataClusterWLFactor` now also implements `FIXED_NODES` (with the
`auto-nodes` calibration) and `CUBATURE`, matching `NcDataClusterWL`'s own
integ-method options. Switching the sweep to `FIXED_NODES` (verified to
match `LNINT` to machine precision) and applying the `auto-lens-nodes` fix
above brought the full 480-cell sweep to 0 crashed cells across all four
methods, with `SeriesLensed` the standout cost/accuracy pick (~1.9x
VarAdd's cost, matching `Laplace`/`FixedQuad`'s accuracy).

## `NcGalaxyShapeFactorFixedQuad`: replaced with a two-branch psi/native-chi_I design (2026)

The domain-intersection design above (three geometry branches: noise-
contained, unit-contained, genuine-lens, plus the effective-radius boundary
fix) is superseded. `Quad` was independently redesigned first (same
session) to integrate $\psi=f_h^{-1}(\chi_L)$ in exact finite polar
coordinates with a puncture correction, fixing a real `Quad` accuracy bug
for $\alpha<2$ Beta populations near $g\sim0.18$ (see above). Porting the
same $\psi$-reparametrization to `FixedQuad` as a single fixed polar grid
(replacing all three geometry branches) fixed that class' own long-standing
narrow-population/branch-cost-cliff issues above as a side effect — but
introduced a NEW failure at $|\epsilon_\mathrm{obs}|$ close to 1 combined
with broad $\sigma_\mathrm{noise}$ (worst case ~16% at $n=15$), traced to
the $\psi$-reparametrization's own Jacobian developing a huge
($\sim10^6\times$) dynamic range near the disc boundary there — a
structural property of the reparametrization itself, not fixable by more
radial resolution alone.

A radial two-panel split (exact split radius via a Möbius circumcircle
construction — three boundary points of the physical noise disc mapped
through $f_h^{-1}$ and circumscribed, since Möbius maps send circles to
circles exactly) fully fixed the "far from pop pull" / small-noise regime,
but not this new one: the difficulty there is the Jacobian itself, not
where the noise kernel's mass sits within the disc. The fix was a second,
complementary scheme instead: integrate directly in $\chi_I$'s own polar
coordinates when the noise kernel is broad enough to cover the whole disc
($1+|\epsilon_\mathrm{obs}|\le8\sigma_\mathrm{noise}$) — no
reparametrization, no Jacobian for the population term (the quadrature
measure's own $r\,\mathrm{d}r$ cancels $P_\mathrm{pop}(r)/(2\pi r)$'s $1/r$
analytically), and, found while implementing it, no puncture correction
either: unlike $\psi$-space, this grid is centered on $\chi_I=0$ itself, so
the noise kernel can be evaluated plainly and safely with no unrelated
singular point to protect against. Confirmed empirically: the plain
(uncorrected) form matches the corrected form's accuracy on the very case
that originally motivated the correction (Beta $\alpha=1.55$), converging
cleanly with node count — the correction was only ever needed because the
OLD geometry-branch design gridded in $\chi_L$-space and inverted to
$\chi_I$, unevenly sparsifying nodes near the singularity; gridding
natively in $\chi_I$ from the start never has that problem.

An 840-case sweep (7 Beta $\alpha$ values, both ellipticity conventions, 5
$|\epsilon_\mathrm{obs}|$, 4 $\sigma_\mathrm{noise}$, 3 shears) confirmed
the two-branch switch closes essentially the whole gap: single-panel and
two-panel alone both peak at 15–17% error, always in the same
$|\epsilon_\mathrm{obs}|\approx0.95,\sigma_\mathrm{noise}=0.35$-type
corner; the switch brings the worst case down to 2.0% (a Beta
$\alpha=3$/tiny-$\sigma_\mathrm{noise}$ resolution edge, not the puncture
problem), median error $\sim10^{-6}$. A third candidate (stereographic
$(u,v)$-plane substitution applied to $\psi$, hypothesized to have its own
Jacobian cancel the reparametrization's) was tried and rejected: it neither
fixes the boundary regime cleanly nor preserves the small-noise regime's
accuracy (up to 600% error on unrelated cases).

See `docs/theory/wl_shape_marginalization_fixed_quad.qmd` for the shipped
design.
