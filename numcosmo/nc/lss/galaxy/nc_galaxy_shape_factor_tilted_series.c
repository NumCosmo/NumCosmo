/***************************************************************************
 *            nc_galaxy_shape_factor_tilted_series.c
 *
 *  Thu Sep 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_tilted_series.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyShapeFactorTiltedSeries:
 *
 * Exponential-tilt evaluation of the intrinsic-ellipticity marginal.
 *
 * A closed-form, table-free alternative to #NcGalaxyShapeFactorMomentSeries'
 * Gaussian clamp. `MomentSeries` carries the *exact*-in-$g$ mean and
 * covariance of the marginal but represents the density itself by a
 * Gaussian, which has no fourth cumulant while the disc-truncated population
 * does: its posterior is unbiased for the shear point estimate, but its
 * credible intervals come out too narrow. This class matches the *same*
 * exact moments and replaces the base measure instead: rather than the
 * maximum-entropy density relative to Lebesgue measure (a Gaussian), it uses
 * the exponential tilt of the exact zero-shear marginal
 * $P_0 = P_\mathrm{pop} * \mathcal{N}(\sigma_\nu)$,
 * $$
 * \ln P(\chi_\mathrm{obs}\mid g) = \ln P_0(|\chi_\mathrm{obs}|)
 *   + \lambda(g)\cdot T(\chi_\mathrm{obs}) - W(g),
 *   \qquad T = (x,\ x^2,\ y^2),
 * $$
 * with $\lambda$ fixed by requiring $\mathbb{E}_\lambda[T]$ to reproduce that
 * exact mean and covariance. Because the shear map fixes the unit circle,
 * $P_0$ already carries the disc's hard edge and every non-Gaussian cumulant
 * the truncation generates; only the smooth shear-induced deformation is
 * expanded, and the result is positive by construction (an exponential
 * family), unlike a density expanded directly in $g$, which goes negative on
 * the counter-shear side.
 *
 * ## The ordering variable is the distortion, not the reduced shear
 *
 * The exact marginal is invariant under $g\to1/g$ -- the classical local
 * degeneracy of weak lensing, which survives the observed-plane noise
 * convolution because the noise kernel is isotropic, so it is a property of
 * the per-galaxy likelihood and not only of the noise-free push-forward. A
 * function of $g$ invariant under $g\to1/g$ is a function of $g+1/g$, hence
 * of the distortion
 * $$ \delta = \frac{2g}{1+g^2} = \frac{2}{g+1/g}, $$
 * and every coefficient of $\lambda(g)$ and $W(g)$ is therefore a function of
 * $\delta$ alone. A truncated polynomial in $g$ cannot represent such a
 * function: the only $g$-polynomials invariant under $g\to1/g$ are the
 * constants, and the exact answer is bounded on $g\in[0,\infty)$ while a
 * non-constant polynomial is not.
 *
 * That is not academic, because $|g|>1$ is reachable on real data:
 * $g=\gamma/(1-\kappa)$ has a pole at $\kappa=1$, and it is crossed both at
 * the innermost fit radii once the sampler visits a high enough mass and
 * through the source-redshift quadrature, which integrates over source
 * planes far enough behind an ordinary cluster to reach it. Ordered in $g$,
 * the truncated $\lambda$ grows without bound there and the model returns
 * spuriously large probabilities -- a few such galaxies are enough to put a
 * cluster likelihood's global maximum at the mass prior's upper edge.
 *
 * So this class orders every series in $\delta$: the coefficients are the
 * $\delta$-coefficients of the exact moments, and the argument the series is
 * evaluated at is $\delta$ itself, which is bounded by 1 and invariant by
 * construction. The model then obeys the $g\to1/g$ degeneracy identically
 * rather than approximately, and is bounded a priori by
 * $\max_{[0,1]}|\lambda|$. $\delta$ is also smooth in $g$ everywhere and
 * stationary at $g=1$, so the model inherits the exact stationarity of the
 * truth at the self-dual point -- the fold $\hat g=\min(g,1/g)$ does not,
 * having a corner there.
 *
 * That stationarity has a price, and it bounds what this class can do. It
 * is $\mathrm{d}\delta/\mathrm{d}\hat g = 0$ at $\hat g=1$, so the
 * inverse $\hat g(\delta) = (1-\sqrt{1-\delta^2})/\delta$ carries a
 * square-root branch point at $\delta=1$. $\lambda$ is analytic in
 * $\hat g$ there and is therefore *not* analytic in $\delta$: its
 * $\delta$-coefficients decay only algebraically near the critical curve,
 * however many of them are kept. The endpoint correction below pins the
 * value at $\delta=1$, which is what restores admissibility, but it cannot
 * repair the local behaviour approaching it -- so accuracy near
 * $|g|\simeq1$ stays limited by the ordering itself, and raising
 * `trunc-order` buys little there. An exact, non-perturbative solve
 * interpolated in $\hat g$ is the construction that escapes this; it is a
 * different class, not a deeper truncation of this one.
 *
 * ## The endpoint correction
 *
 * Ordering in $\delta$ is necessary but not sufficient. The exact tangential
 * variance has a $\delta$-expansion whose coefficients are positive beyond
 * the second, so every truncation discards a positive tail and undershoots
 * $C_t$ at every $\delta$; and the exact $C_t$ *vanishes* at $\delta=1$,
 * where the disc collapses onto the single point $\chi=1$ and the marginal
 * is pure noise. The undershoot therefore drives $C_t$ negative near the
 * critical curve. That is not an inaccuracy but an impossibility: $C_t>0$
 * and $C_\times>0$ hold if and only if the target is the moment vector of
 * some distribution, so a negative $C_t$ leaves the tilt equations with no
 * solution at all, at any order. Raising the order does not help -- the
 * undershoot is structural.
 *
 * The repair uses the one piece of exact information a truncation cannot
 * have on its own. The individual omitted coefficients are unknown, but
 * their *sum* is known, because the value at $\delta=1$ is known in closed
 * form. Adding
 * $$ \bigl(\text{endpoint} - \textstyle\sum_{p\le N}\text{retained}\bigr)\,
 *    \delta^{q} $$
 * at the first order $q>N$ of the right parity leaves every retained
 * coefficient untouched -- so the small-$\delta$ accuracy of order $N$ is
 * preserved exactly -- while making the truncation exact at $\delta=1$. See
 * `_tilted_series_target_series()` for the endpoint values and the identity
 * behind the construction.
 *
 * `trunc-order` is the number of true moment orders $N$ retained. The two
 * orders above it carry the endpoint terms and no new moment information, so
 * the internal solve runs at $N+2$; the cost of the solve is set by that
 * larger order.
 *
 * ## Solving for lambda
 *
 * $\lambda(\delta)$ is solved as a power series in $\delta$, one fixed
 * $3\times3$ linear solve per order: no Newton iteration, no
 * population-indexed interpolation table. This is possible because
 * $\mathbb{E}_{P_0}[x^2] = M_2/2+\sigma_\nu^2$ identically equals the exact
 * zero-shear covariance, so $\lambda(0)=\bm0$ exactly and the moment
 * conditions can be expanded order by order. Each order's linear solve uses
 * the closed-form $V=\mathrm{Cov}_{P_0}(T)$; the right-hand side at every
 * order comes from a noise-integrated generating function
 * $Z(\lambda)=\mathbb{E}_{P_0}[e^{\lambda\cdot T}]$, evaluated as a formal
 * power series through the population's own radial moments
 * ($M_{2k}=$ nc_galaxy_shape_pop_moment_2k()) -- an exact, arbitrary-order,
 * population-generic route that never hardcodes $\lambda^{(p)}$ as a
 * symbolic rational function of a fixed set of moments.
 *
 * At order $p$ in the order-by-order solve, only the $\delta^p$ coefficient
 * of each generating-function series is ever consumed, and every series
 * operation used here (product, reciprocal, log) is lower-triangular --
 * coefficient $p$ of the output depends only on input coefficients $\le p$.
 * `_tilted_series_solve()` exploits that by truncating every series computed
 * inside the loop to the *current* order $p$ rather than the full solve
 * order, which is an exact optimisation, not an approximation.
 *
 * Three further exact optimisations sit on top of that, all of them
 * reorganisations that reproduce the unoptimised coefficients bit for bit:
 *
 * - **Valuation-aware products.** $\lambda_1$ is odd in $\delta$ and
 *   $\lambda_2,\lambda_3$ even, so in `_tilted_series_compute_D()`
 *   $\alpha=\lambda_1/a_2$ has valuation 1 and $\beta,\gamma$ valuation 2:
 *   $\alpha^{p}\beta^{q}\gamma^{s}$ is structurally zero below
 *   $\delta^{p+2q+2s}$. `series_mul_v()`/`series_axpy_v()` skip those zeros
 *   instead of multiplying them, which at the top of the triple sum turns a
 *   full Cauchy product into a single term.
 * - **Precomputed combinatorial weights.** The $S_j/(p!q!s!)$ factors of the
 *   triple sum depend only on the population moments $M_{2k}$, never on
 *   $\lambda$, so they are built once per solve
 *   (`_tilted_series_Stab_new()`) instead of being recomputed inside every
 *   `_tilted_series_compute_D()` pass.
 * - **No post-loop pass at odd order.** $W$ is even, so `_eval()` never
 *   reads past $2\lfloor N_\mathrm{solve}/2\rfloor$; at odd solve order the
 *   loop's own last pass already left $a_2,a_3,D_0$ correct to that order.
 *
 * All the series scratch is carved from one allocation per solve rather than
 * the ~100 small ones the straightforward form needs.
 *
 * The target moments are the *same* ones `MomentSeries` builds from the
 * population's radial moments, so this class shares that build
 * (nc_galaxy_shape_factor_moment_series_private.h) rather than duplicating
 * it. That build produces $g$-coefficients; converting them to
 * $\delta$-coefficients is linear and population-independent, so it is
 * folded into the shared tables once at construction and never repeated per
 * galaxy (`_tilted_series_delta_tables()`). Both ellipticity conventions
 * come for free, since boundary invariance holds for `TRACE_DET` too and
 * both maps are the identity at $g=0$.
 *
 * The per-galaxy constant $\ln P_0(|\chi_\mathrm{obs}|)$ is $g$-independent,
 * so any quadrature error in it is a per-galaxy *additive* constant in
 * $\ln P$: it cancels identically from $\partial_g\ln P$ and never biases
 * the shear point estimate or the sandwich variance. This is what licenses
 * evaluating it by one population-generic, localized, log-space
 * Gauss-Legendre quadrature (never a table, and never the Marcum
 * $Q$-function closed form of a disc-truncated Gaussian population --
 * unneeded because a per-galaxy population such as
 * #NcGalaxyShapePopGaussLocal has a per-galaxy, not fitted, width, so this
 * quantity is computed at most once per galaxy; see
 * #NcGalaxyShapeFactorMomentSeries' own note on
 * nc_galaxy_shape_factor_update_data_pop() for when that stops being true).
 *
 * Unlike `MomentSeries`, both `eval_marginal()` and `eval_ln_marginal()` are
 * equally expensive here: the fixed-nodes evaluation path
 * (nc_galaxy_shape_factor_eval_at_nodes()) calls `eval_marginal`, so an
 * `exp()` is unavoidable on that path regardless of which hook the caller
 * nominally wants -- this class does not claim to be cheaper than a Gaussian
 * evaluation.
 *
 * ## Sanity guards
 *
 * Two tests run at evaluation time, not silent clamps, exactly like
 * `MomentSeries`' own covariance-positivity guard: $Z(\lambda)$'s natural
 * domain $\max(\lambda_2,\lambda_3)<1/2\sigma_\nu^2$, and an exact ceiling
 * on the marginal itself. The marginal is a convolution of a probability
 * density with the noise kernel, so $P\le1/2\pi\sigma_\nu^2$ for every $g$
 * and every $\chi_\mathrm{obs}$; the second test trips well above that. The
 * ceiling is the one that catches a $\lambda$ far too *negative*, which
 * makes $P$ a spurious spike -- no test on $\max(\lambda_2,\lambda_3)$ can
 * see that, since the exact $\lambda_2,\lambda_3$ are themselves negative.
 *
 * Whether raising `trunc-order` helps when a guard fires depends on which
 * of two things went wrong, and they pull in opposite directions.
 *
 * - On a narrow, low-noise population the model's peak excess over the
 *   ceiling *grows* with order, steadily and without the order-by-order
 *   solve ever failing. Raising the order there makes the very quantity the
 *   ceiling watches worse, so it is not a remedy.
 * - At low order the target itself can still be far from its endpoint, which
 *   puts a large coefficient in the single correction term and distorts the
 *   target just below $\delta=1$. A $\lambda_2$ over the domain bound near
 *   $|g|=1$ is the signature, and there raising the order does clear it,
 *   because the retained sum reaches the endpoint on its own.
 *
 * Neither can touch the boundedness, which is structural. Read the reported
 * $\lambda_2$ against the bound before deciding: over the bound at a low
 * order points at the second, an excess over the ceiling alone at the first.
 *
 * An MCMC walker reaching that region is a routine event -- walkers are
 * initialised over the whole prior box -- so the default is to report it,
 * not to abort the process. The guard returns $0$ from `eval_marginal()`
 * (and $-\infty$ from `eval_ln_marginal()`, keeping the two consistent),
 * which routes into #NcDataClusterWLFactor's existing `NC_GALAXY_LOW_PROB`
 * path and pushes the sampler back out of the region. Occurrences are
 * counted and readable through
 * nc_galaxy_shape_factor_tilted_series_get_domain_error_count(); a non-zero
 * count means the model returned something no probability density can
 * return, and should be reported rather than worked around by widening
 * bounds. Set #NcGalaxyShapeFactorTiltedSeries:strict-domain to restore the
 * fatal behaviour.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_tilted_series.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moment_series_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * ---- Formal power-series arithmetic in delta, truncated at order N ----
 *
 * Every series below is a plain gdouble[N+1] array, coefficient of delta^p
 * at index p. delta = 2g/(1+g^2) is this class' ordering variable throughout
 * (see the class doc comment); the only place a g-series ever appears is the
 * one-off change of basis in _tilted_series_delta_tables(), which converts
 * the shared moment tables from g to delta once, at construction.
 *
 * These are the generic building blocks for the noise-integrated
 * generating-function solve; none of them is population- or order-specific.
 */

static inline void
series_zero (gdouble *out, guint order)
{
  memset (out, 0, (order + 1) * sizeof (gdouble));
}

static inline void
series_copy (gdouble *out, const gdouble *in, guint order)
{
  memcpy (out, in, (order + 1) * sizeof (gdouble));
}

/* out += scale * in */
static inline void
series_axpy (gdouble *out, const gdouble *in, gdouble scale, guint order)
{
  guint p;

  for (p = 0; p <= order; p++)
    out[p] += scale * in[p];
}

/* out = a * b (Cauchy product, truncated); out must not alias a or b. */
static void
series_mul (gdouble *out, const gdouble *a, const gdouble *b, guint order)
{
  guint p;

  for (p = 0; p <= order; p++)
  {
    gdouble s = 0.0;
    guint k;

    for (k = 0; k <= p; k++)
      s += a[k] * b[p - k];

    out[p] = s;
  }
}

/* out = a * b where a is known to have valuation @va and b valuation @vb
 * (every coefficient below the valuation is a structural zero). The result
 * is identical to series_mul() on the same inputs -- the Cauchy sums simply
 * do not visit the zeros. out must not alias a or b. */
static void
series_mul_v (gdouble *out, const gdouble *a, guint va, const gdouble *b, guint vb, guint order)
{
  const guint v = va + vb;
  guint p;

  if (v > order)
  {
    memset (out, 0, (order + 1) * sizeof (gdouble));

    return;
  }

  memset (out, 0, v * sizeof (gdouble));

  for (p = v; p <= order; p++)
  {
    const guint khi = p - vb;
    gdouble s       = 0.0;
    guint k;

    for (k = va; k <= khi; k++)
      s += a[k] * b[p - k];

    out[p] = s;
  }
}

/* out += scale * in, where in has valuation @v. */
static inline void
series_axpy_v (gdouble *out, const gdouble *in, gdouble scale, guint v, guint order)
{
  guint p;

  for (p = v; p <= order; p++)
    out[p] += scale * in[p];
}

/* out = 1/a, requires a[0] != 0; out must not alias a. */
static void
series_recip (gdouble *out, const gdouble *a, guint order)
{
  guint p;

  out[0] = 1.0 / a[0];

  for (p = 1; p <= order; p++)
  {
    gdouble s = 0.0;
    guint k;

    for (k = 1; k <= p; k++)
      s += a[k] * out[p - k];

    out[p] = -s / a[0];
  }
}

/* out = ln(f), requires f[0] > 0; out must not alias f. Via
 * (ln f)' = f'/f, integrated term by term: reuses series_recip/series_mul,
 * no new machinery beyond what the D_u solve already needs. */
static void
series_log (gdouble *out, const gdouble *f, guint order)
{
  gdouble *fp   = g_new0 (gdouble, order + 1);
  gdouble *finv = g_new0 (gdouble, order + 1);
  gdouble *q    = g_new0 (gdouble, order + 1);
  guint p;

  for (p = 0; p < order; p++)
    fp[p] = (p + 1) * f[p + 1];

  series_recip (finv, f, order);
  series_mul (q, fp, finv, order);

  out[0] = log (f[0]);

  for (p = 1; p <= order; p++)
    out[p] = q[p - 1] / p;

  g_free (fp);
  g_free (finv);
  g_free (q);
}

/* table[n] = base^n, n=0..maxn, written into caller-provided storage @flat
 * as consecutive length-(order+1) rows (row n at @flat + n*(order+1)), so
 * the solve's hot path allocates nothing. @base has valuation @base_val,
 * hence base^n has valuation n*@base_val. */
static void
series_pow_table (gdouble *flat, const gdouble *base, guint base_val, guint maxn, guint order)
{
  const guint row = order + 1;
  guint n;

  memset (flat, 0, (maxn + 1) * row * sizeof (gdouble));
  flat[0] = 1.0;

  for (n = 1; n <= maxn; n++)
    series_mul_v (&flat[n * row], &flat[(n - 1) * row], (n - 1) * base_val, base, base_val, order);
}

static inline gdouble
factorial_d (guint n)
{
  gdouble r = 1.0;
  guint k;

  for (k = 2; k <= n; k++)
    r *= k;

  return r;
}

/* Double factorial, n!! with the convention (-1)!!=0!!=1, valid for both
 * parities of n -- the isotropic-moment formula below uses both. */
static inline gdouble
dfact (gint n)
{
  gdouble r = 1.0;

  while (n > 1)
  {
    r *= n;
    n -= 2;
  }

  return r;
}

/*
 * <(Re chi_I)^i_pow (Im chi_I)^j_pow> for an isotropic intrinsic population,
 * zero unless both powers are even. Writing chi_I = |chi_I| e^{i phi} with
 * phi uniform, the angular average of cos^{2ii} sin^{2jj} factorises out as
 * (2ii-1)!!(2jj-1)!!/(2ii+2jj)!!, leaving the radial moment
 * M_{2(ii+jj)} = <|chi_I|^{2(ii+jj)}>:
 *   = M_{2(ii+jj)} * (2ii-1)!!(2jj-1)!!/(2ii+2jj)!!, ii=i_pow/2, jj=j_pow/2.
 * @M holds M_2k at M[k]; @n_M must cover k=ii+jj (asserted).
 */
static inline gdouble
iso_moment (const gdouble *M, guint n_M, guint i_pow, guint j_pow)
{
  if (((i_pow % 2) != 0) || ((j_pow % 2) != 0))
    return 0.0;

  {
    const gint ii = (gint) (i_pow / 2);
    const gint jj = (gint) (j_pow / 2);
    const guint k = (guint) (ii + jj);

    g_assert_cmpuint (k, <, n_M);

    return M[k] * dfact (2 * ii - 1) * dfact (2 * jj - 1) / dfact (2 * ii + 2 * jj);
  }
}

/*
 * ---- The noise-integrated generating-function solve ----
 *
 * Computes a2(d)=1-2*lambda_2(d)*sn2, a3(d)=1-2*lambda_3(d)*sn2 and their
 * reciprocals, then D_u(d)=E_chis[chis_x^u chis_y^v e^Q] for the four (u,v)
 * in {(0,0),(1,0),(2,0),(0,2)} needed by the t1/t2/t3 model moments, with
 * Q = (lambda_1 sx + lambda_2 sx^2)/a2 + lambda_3 sy^2/a3.
 *
 * The construction: E_lambda[T] is a ratio of Gaussian integrals against
 * P_0 = P_pop * N(sn), and the noise integral can be done in closed form
 * for fixed intrinsic chi_I, which is what leaves a2, a3 and Q behind.
 * Expanding e^Q in alpha=lambda_1/a2, beta=lambda_2/a2, gamma=lambda_3/a3
 * turns what remains into a finite triple sum over (pp,qq,ss) with
 * pp+2qq+2ss<=order, each term an isotropic moment of chi_I
 * (population-generic, from nc_galaxy_shape_pop_moment_2k()) times a fixed
 * combinatorial factor -- no quadrature, no table.
 *
 * @order truncates every series computed here (may be less than the lengths
 * of the @lam1/@lam2/@lam3/output arrays, which are always caller-sized for
 * the full solve order): coefficients above @order are left untouched. This
 * is exact, not approximate -- every series operation below is
 * lower-triangular, so a caller that only reads coefficients <=@order back
 * out gets the same values truncation at the full order would give. The
 * order-by-order solve in _tilted_series_solve() exploits this by calling
 * with @order=p inside its loop.
 *
 * @Stab holds the triple sum's combinatorial weights, built once per solve
 * by _tilted_series_Stab_new() with stride @half_N (which is half the
 * caller's FULL order, NOT @order/2 -- the table is shared across every call
 * whatever @order each one uses). @scratch is caller-provided working
 * storage of at least (5 + (order+1) + 2*(order/2+1)) * (order+1) doubles.
 */
static void
_tilted_series_compute_D (const gdouble *lam1, const gdouble *lam2, const gdouble *lam3,
                          gdouble sn2, guint order,
                          const gdouble *Stab, guint half_N, gdouble *scratch,
                          gdouble *a2, gdouble *a3, gdouble *a2_inv, gdouble *a3_inv,
                          gdouble *D0, gdouble *D1v, gdouble *D2v, gdouble *D3v)
{
  const guint row        = order + 1;
  const guint half_order = order / 2;
  const guint stride     = half_N + 1;
  gdouble *alpha         = scratch;
  gdouble *beta          = scratch + row;
  gdouble *gamma         = scratch + 2 * row;
  gdouble *tmp1          = scratch + 3 * row;
  gdouble *tmp2          = scratch + 4 * row;
  gdouble *alpha_pows    = scratch + 5 * row;
  gdouble *beta_pows     = alpha_pows + (order + 1) * row;
  gdouble *gamma_pows    = beta_pows + (half_order + 1) * row;
  guint pp, qq, ss, k;

  series_zero (a2, order);
  series_zero (a3, order);
  a2[0] = 1.0;
  a3[0] = 1.0;

  for (k = 1; k <= order; k++)
  {
    a2[k] = -2.0 * sn2 * lam2[k];
    a3[k] = -2.0 * sn2 * lam3[k];
  }

  series_recip (a2_inv, a2, order);
  series_recip (a3_inv, a3, order);

  series_mul (alpha, lam1, a2_inv, order);
  series_mul (beta, lam2, a2_inv, order);
  series_mul (gamma, lam3, a3_inv, order);

  /* lambda_1 is odd in delta and lambda_2/lambda_3 even, so alpha has
   * valuation 1 and beta/gamma valuation 2 (see the class doc comment). */
  series_pow_table (alpha_pows, alpha, 1, order, order);
  series_pow_table (beta_pows, beta, 2, half_order, order);
  series_pow_table (gamma_pows, gamma, 2, half_order, order);

  series_zero (D0, order);
  series_zero (D1v, order);
  series_zero (D2v, order);
  series_zero (D3v, order);

  for (pp = 0; pp <= order; pp++)
  {
    for (qq = 0; pp + 2 * qq <= order; qq++)
    {
      const guint v12 = pp + 2 * qq;
      const gdouble *t1;

      /* beta^0 = 1 and gamma^0 = 1: skip the identity products rather
       * than running a full Cauchy product against a constant series. */
      if (qq == 0)
      {
        t1 = &alpha_pows[pp * row];
      }
      else
      {
        series_mul_v (tmp1, &alpha_pows[pp * row], pp, &beta_pows[qq * row], 2 * qq, order);
        t1 = tmp1;
      }

      for (ss = 0; v12 + 2 * ss <= order; ss++)
      {
        const guint v    = v12 + 2 * ss;
        const gdouble *S = &Stab[4 * (((pp * stride) + qq) * stride + ss)];
        const gdouble *t2;

        if ((S[0] == 0.0) && (S[1] == 0.0) && (S[2] == 0.0) && (S[3] == 0.0))
          continue;

        if (ss == 0)
        {
          t2 = t1;
        }
        else
        {
          series_mul_v (tmp2, t1, v12, &gamma_pows[ss * row], 2 * ss, order);
          t2 = tmp2;
        }

        if (S[0] != 0.0)
          series_axpy_v (D0, t2, S[0], v, order);

        if (S[1] != 0.0)
          series_axpy_v (D1v, t2, S[1], v, order);

        if (S[2] != 0.0)
          series_axpy_v (D2v, t2, S[2], v, order);

        if (S[3] != 0.0)
          series_axpy_v (D3v, t2, S[3], v, order);
      }
    }
  }
}

/*
 * The (pp,qq,ss) combinatorial weights S_j/(pp! qq! ss!) of
 * _tilted_series_compute_D()'s triple sum. They are functions of the
 * population moments @M alone -- lambda enters the sum only through
 * alpha/beta/gamma -- so they are built once per solve here instead of being
 * recomputed inside each compute_D() pass. Layout: four consecutive doubles
 * (S0,S1,S2,S3)/fac at index ((pp*stride + qq)*stride + ss),
 * stride = N/2 + 1. Entries outside the triple sum's own
 * (pp + 2qq + 2ss <= N) region are left zero and never read.
 *
 * @n_M must cover the largest radial moment any entry needs. The tightest
 * index comes from S[2] (and S[3]) at even pp, where iso_moment() is asked
 * for k = pp/2 + qq + ss + 1 <= floor(N/2) + 1; see the n_M derivation in
 * _nc_galaxy_shape_factor_tilted_series_constructed().
 */
static gdouble *
_tilted_series_Stab_new (const gdouble *M, guint n_M, guint N)
{
  const guint half_N = N / 2;
  const guint stride = half_N + 1;
  gdouble *Stab      = g_new0 (gdouble, 4 * (N + 1) * stride * stride);
  guint pp, qq, ss;

  for (pp = 0; pp <= N; pp++)
  {
    for (qq = 0; (qq <= half_N) && (pp + 2 * qq <= N); qq++)
    {
      for (ss = 0; (ss <= half_N) && (pp + 2 * qq + 2 * ss <= N); ss++)
      {
        const guint sx_pow = pp + 2 * qq;
        const guint sy_pow = 2 * ss;
        const gdouble inv  = 1.0 / (factorial_d (pp) * factorial_d (qq) * factorial_d (ss));
        gdouble *S         = &Stab[4 * (((pp * stride) + qq) * stride + ss)];

        S[0] = iso_moment (M, n_M, sx_pow, sy_pow) * inv;
        S[1] = iso_moment (M, n_M, sx_pow + 1, sy_pow) * inv;
        S[2] = iso_moment (M, n_M, sx_pow + 2, sy_pow) * inv;
        S[3] = iso_moment (M, n_M, sx_pow, sy_pow + 2) * inv;
      }
    }
  }

  return Stab;
}

/*
 * Order-by-order solve of lambda(delta) and the closed-form log-normaliser
 * W(delta) via the recursion V*lambda^(p) = Delta^(p) - R^(p), where R^(p)
 * is the delta^p coefficient of grad W evaluated at the trial series with
 * lambda^(p) itself omitted. @Delta1/@Delta2/@Delta3 are this galaxy's
 * target moment series in delta (_tilted_series_target_series());
 * @lam1/@lam2/@lam3/@Wser (each length N+1) receive the result.
 *
 * @N is the SOLVE order, i.e. trunc-order plus the two orders carrying the
 * endpoint corrections -- not trunc-order itself.
 *
 * V = [[A,0,0],[0,B,kappa],[0,kappa,B]] is the closed-form covariance of
 * T=(x,x^2,y^2) under P_0; it is positive-definite unconditionally (as
 * sn->0, 2*kappa+A^2 -> Var(|chi_I|^2)/4 > 0 and kappa+A^2 -> M_4/8 > 0, and
 * sn^2 only raises A), so the guard below is an internal invariant, not a
 * reachable user-facing failure mode (contrast MomentSeries' own covariance
 * guard).
 *
 * After the order-by-order loop, one more pass recomputes a2, a3, D0 with
 * the now-complete lambda, and W is read off in closed form,
 * W = -0.5 ln a2 - 0.5 ln a3 + lambda_1^2 sn2/(2 a2) + ln D0 -- "free" given
 * a2/a3/D0, rather than by integrating W' = t.lambda' term by term (that
 * identity is a check on this routine, not its code path).
 */
static void
_tilted_series_solve (guint N, const gdouble *M, guint n_M, gdouble sn2,
                      const gdouble *Delta1, const gdouble *Delta2, const gdouble *Delta3,
                      gdouble *lam1, gdouble *lam2, gdouble *lam3, gdouble *Wser)
{
  const gdouble A     = 0.5 * M[1] + sn2;
  const gdouble kappa = 0.125 * (M[2] - 2.0 * M[1] * M[1]);
  const gdouble B     = 3.0 * kappa + 2.0 * A * A;
  const gdouble det2  = B * B - kappa * kappa;

  const guint row = N + 1;

  /* One arena for every temporary: the 17 series below, then
   * _tilted_series_compute_D()'s own scratch (5 series plus the three flat
   * power tables alpha^0..alpha^N, beta^0..beta^{N/2}, gamma^0..gamma^{N/2},
   * sized for the largest @order any pass uses, which is N). */
  const guint n_scr  = (5 + (N + 1) + 2 * (N / 2 + 1)) * row;
  gdouble *arena     = g_new (gdouble, 17 * row + n_scr);
  gdouble *a2        = arena;
  gdouble *a3        = arena + 1 * row;
  gdouble *a2_inv    = arena + 2 * row;
  gdouble *a3_inv    = arena + 3 * row;
  gdouble *D0        = arena + 4 * row;
  gdouble *D1v       = arena + 5 * row;
  gdouble *D2v       = arena + 6 * row;
  gdouble *D3v       = arena + 7 * row;
  gdouble *D0_inv    = arena + 8 * row;
  gdouble *r1s       = arena + 9 * row;
  gdouble *r2s       = arena + 10 * row;
  gdouble *r3s       = arena + 11 * row;
  gdouble *tmp_a     = arena + 12 * row;
  gdouble *tmp_b     = arena + 13 * row;
  gdouble *t1_model  = arena + 14 * row;
  gdouble *t2_model  = arena + 15 * row;
  gdouble *t3_model  = arena + 16 * row;
  gdouble *scratch   = arena + 17 * row;
  gdouble *Stab      = _tilted_series_Stab_new (M, n_M, N);
  const guint half_N = N / 2;
  /* W is even in delta, so this is the highest coefficient _eval() reads. */
  const guint W_ord = 2 * (N / 2);
  guint p;

  /* Internal invariant, not a reachable guard -- see this function's own
   * doc comment. */
  g_assert (det2 > 0.0);

  series_zero (lam1, N);
  series_zero (lam2, N);
  series_zero (lam3, N);

  for (p = 1; p <= N; p++)
  {
    gdouble rhs1, rhs2, rhs3;

    /* lam1/lam2/lam3 here carry only orders < p (order p and above are
     * still zero), which is exactly R^(p).
     *
     * Only the delta^p coefficient of rhs1/rhs2/rhs3 is read below, and
     * every series op here is lower-triangular, so truncating the whole loop
     * body at @p instead of the full N is exact -- see the class doc comment
     * and _tilted_series_compute_D()'s. Cost drops from N*C(N) to
     * sum_{p=1}^N C(p). */
    _tilted_series_compute_D (lam1, lam2, lam3, sn2, p, Stab, half_N, scratch,
                              a2, a3, a2_inv, a3_inv, D0, D1v, D2v, D3v);

    series_recip (D0_inv, D0, p);
    series_mul (r1s, D1v, D0_inv, p);
    series_mul (r2s, D2v, D0_inv, p);
    series_mul (r3s, D3v, D0_inv, p);

    /* t1_model = (r1s + sn2*lam1) * a2_inv */
    series_copy (tmp_a, r1s, p);
    series_axpy (tmp_a, lam1, sn2, p);
    series_mul (t1_model, tmp_a, a2_inv, p);

    /* t2_model = (r2s + 2*sn2*lam1*r1s + sn2^2*lam1^2) * a2_inv^2 + sn2*a2_inv */
    series_mul (tmp_a, lam1, r1s, p);
    series_copy (tmp_b, r2s, p);
    series_axpy (tmp_b, tmp_a, 2.0 * sn2, p);
    series_mul (tmp_a, lam1, lam1, p);
    series_axpy (tmp_b, tmp_a, sn2 * sn2, p);
    series_mul (tmp_a, a2_inv, a2_inv, p);
    series_mul (t2_model, tmp_b, tmp_a, p);
    series_axpy (t2_model, a2_inv, sn2, p);

    /* t3_model = r3s * a3_inv^2 + sn2*a3_inv */
    series_mul (tmp_a, a3_inv, a3_inv, p);
    series_mul (t3_model, r3s, tmp_a, p);
    series_axpy (t3_model, a3_inv, sn2, p);

    /* The model moments t*_model carry the noise variance in their p=0
     * coefficient; the targets are increments from delta=0, where that same
     * constant sits. At p>=1 the two therefore compare directly. */
    rhs1 = Delta1[p] - t1_model[p];
    rhs2 = Delta2[p] - t2_model[p];
    rhs3 = Delta3[p] - t3_model[p];

    lam1[p] = rhs1 / A;
    lam2[p] = (B * rhs2 - kappa * rhs3) / det2;
    lam3[p] = (B * rhs3 - kappa * rhs2) / det2;
  }

  /* W is even in delta: lambda_2/lambda_3 are even, so a2 and a3 are;
   * lambda_1^2 is; and D0 keeps only the even powers of alpha (its S0 weight
   * vanishes unless pp + 2qq is even, i.e. unless pp is even). So every odd
   * coefficient of Wser is an exact zero and _eval() never reads past
   * W_ord = 2*floor(N/2).
   *
   * At odd N, W_ord = N-1 and the p = N pass of the loop above already left
   * a2, a3, a2_inv and D0 correct to that order: only lambda^(N) was missing
   * from its inputs, and every series operation here is lower-triangular, so
   * it can perturb coefficient N alone. The extra full-order pass is
   * therefore needed at even N only. */
  if (W_ord == N)
    _tilted_series_compute_D (lam1, lam2, lam3, sn2, N, Stab, half_N, scratch,
                              a2, a3, a2_inv, a3_inv, D0, D1v, D2v, D3v);

  series_zero (Wser, N);

  series_log (tmp_a, a2, W_ord);
  series_axpy (Wser, tmp_a, -0.5, W_ord);

  series_log (tmp_a, a3, W_ord);
  series_axpy (Wser, tmp_a, -0.5, W_ord);

  series_mul (tmp_a, lam1, lam1, W_ord);
  series_mul (tmp_b, tmp_a, a2_inv, W_ord);
  series_axpy (Wser, tmp_b, 0.5 * sn2, W_ord);

  series_log (tmp_a, D0, W_ord);
  series_axpy (Wser, tmp_a, 1.0, W_ord);

  g_free (Stab);
  g_free (arena);
}

/*
 * ---- Change of basis: the moment tables from g to delta ----
 *
 * Taylor coefficients of the inverse map g(delta), odd in delta:
 *
 *   g(delta) = (1 - sqrt(1-delta^2))/delta = delta/2 + delta^3/8 + ...,
 *   c_1 = 1/2,  c_{k+1} = c_k (2k-1)/(2(k+1)).
 *
 * This inverts delta = 2g/(1+g^2) on the branch vanishing at the origin. The
 * quadratic's two roots are g and 1/g, so choosing a branch is exactly
 * choosing a representative of the g -> 1/g degeneracy; the series returns
 * min(g,1/g), which is why it is bounded by 1 however large |g| gets.
 */
static void
_tilted_series_g_of_delta (gdouble *gd, guint order)
{
  gdouble ck = 0.5;
  guint k;

  series_zero (gd, order);

  for (k = 1; 2 * k - 1 <= order; k++)
  {
    gd[2 * k - 1] = ck;
    ck           *= (2.0 * k - 1.0) / (2.0 * (k + 1.0));
  }
}

/*
 * Rewrites the shared moment tables from g-coefficients to
 * delta-coefficients, in place, once at construction.
 *
 * _nc_galaxy_shape_factor_moment_series_build_tables() produces tables in
 * g: tab_m row l holds the g^{2l+1} coefficient of the mean and tab_v/tab_w
 * row l the g^{2l} coefficients of the second moments, each as n_moments
 * weights to be contracted against the population's radial moments M_2k.
 * This class needs the same quantities as coefficients of delta.
 *
 * Both the change of basis and the per-galaxy contraction against M are
 * linear, so they commute: converting the tables once here is identical to
 * converting each galaxy's contracted series, and costs nothing per galaxy.
 *
 * The row layout survives untouched, which is what makes this a drop-in.
 * g(delta) is odd with no constant term, so it maps odd exponents to odd and
 * even to even, and the row index is the same function of the exponent in
 * both variables. Two consequences worth stating, because the rest of the
 * class relies on them:
 *
 *  - delta-coefficient p of f(g(delta)) depends only on the g-coefficients
 *    <= p, since g(delta) has valuation 1. Rows 0..N therefore come out as
 *    the EXACT delta-coefficients of the exact moments -- the change of
 *    basis approximates nothing and needs no row the build did not produce.
 *  - Row 0 of tab_v/tab_w is zeroed rather than converted: t(0) is handled
 *    separately by the lambda(0)=0 theorem, and g(delta)^0 = 1 contributes
 *    to delta-coefficient 0 alone, so dropping it cannot disturb any row
 *    above it.
 *
 * The fold also mixes rows only, never columns, so a moment column that was
 * zero throughout stays zero and @n_moments remains a valid bound.
 */
static void
_tilted_series_delta_tables (guint N, guint n_m, guint n_v, guint n_moments,
                             gdouble *tab_m, gdouble *tab_v, gdouble *tab_w)
{
  const guint row  = N + 1;
  const gsize sz_m = (gsize) n_m * n_moments;
  const gsize sz_v = (gsize) n_v * n_moments;
  gdouble *gd      = g_new0 (gdouble, row);
  /* gpow[n] = g(delta)^n as a delta-series, n = 0..N. */
  gdouble *gpow  = g_new0 (gdouble, (gsize) (N + 1) * row);
  gdouble *src_m = g_new (gdouble, sz_m);
  gdouble *src_v = g_new (gdouble, sz_v);
  gdouble *src_w = g_new (gdouble, sz_v);
  guint l, a, k;

  memcpy (src_m, tab_m, sz_m * sizeof (gdouble));
  memcpy (src_v, tab_v, sz_v * sizeof (gdouble));
  memcpy (src_w, tab_w, sz_v * sizeof (gdouble));

  _tilted_series_g_of_delta (gd, N);
  series_pow_table (gpow, gd, 1, N, N);

  /* Mean: target exponent p = 2l+1, fed by source exponents 2k+1 <= p. */
  for (l = 0; l < n_m; l++)
  {
    const guint p = 2 * l + 1;

    for (a = 0; a < n_moments; a++)
    {
      gdouble s = 0.0;

      for (k = 0; k <= l; k++)
        s += gpow[(2 * k + 1) * row + p] * src_m[k * n_moments + a];

      tab_m[l * n_moments + a] = s;
    }
  }

  /* Second moments: target exponent p = 2l, fed by source exponents
   * 2k <= p with k >= 1. */
  for (l = 0; l < n_v; l++)
  {
    const guint p = 2 * l;

    for (a = 0; a < n_moments; a++)
    {
      gdouble sv = 0.0;
      gdouble sw = 0.0;

      for (k = 1; k <= l; k++)
      {
        const gdouble b = gpow[(2 * k) * row + p];

        sv += b * src_v[k * n_moments + a];
        sw += b * src_w[k * n_moments + a];
      }

      tab_v[l * n_moments + a] = sv;
      tab_w[l * n_moments + a] = sw;
    }
  }

  g_free (gd);
  g_free (gpow);
  g_free (src_m);
  g_free (src_v);
  g_free (src_w);
}

/*
 * ---- Target moment series in delta, with the endpoint correction ----
 *
 * Delta^(p) is the delta^p coefficient of t = (mu, C_t + mu^2, C_x), the
 * same exact moments MomentSeries carries, as an increment from delta = 0
 * (the p = 0 term is handled separately by the lambda(0) = 0 theorem).
 * Reading them off needs no polynomial algebra beyond the tab_m/tab_v/tab_w
 * contraction against the population's radial moments: mu has only odd
 * powers and C_t + mu^2 / C_x only even ones. The mu^2 that MomentSeries
 * itself subtracts to turn <(Re S)^2> into Var(Re S) cancels exactly against
 * the mu^2 added back in here -- both are the self-convolution of the same
 * mu series -- so this class never performs that subtraction at all. The
 * tables arrive already expressed in delta; see
 * _tilted_series_delta_tables().
 *
 * THE ENDPOINT CORRECTION
 *
 * Truncating the delta-series is inadmissible near delta = 1, and the repair
 * has to happen here because this is where the target is built. The
 * discarded tail of C_t is positive term by term, so the truncation
 * undershoots C_t at every delta; and the exact C_t vanishes at delta = 1,
 * so the undershoot drives it negative there. A target with C_t < 0 is not
 * the moment vector of any distribution, so the solve then has no solution
 * at all rather than an inaccurate one -- and raising the order does not
 * help, the undershoot being structural.
 *
 * The individual omitted coefficients are unknown, but their SUM is not. At
 * delta = 1 the shear map sends the whole unit disc to the single point
 * chi = 1 -- the critical curve -- so the marginal there is just the noise
 * kernel centred on (1,0), and the increments are
 *
 *   Delta1 -> 1,   Delta2 -> 1 - M_2/2,   Delta3 -> -M_2/2,
 *
 * population-independent but for M_2. (At delta = 0, isotropy gives
 * <(Re S)^2> = <(Im S)^2> = M_2/2, which is what is subtracted off; the
 * noise variance is delta-independent and cancels from both ends.) These
 * values hold for BOTH ellipticity conventions, so they need no branch on
 * ellip-conv: at g = 1 the TRACE map
 * (chi_s + 2g + g^2 conj(chi_s))/(1 + |g|^2 + 2 Re(g conj(chi_s))) and the
 * TRACE_DET map (eps_s + g)/(1 + conj(g) eps_s) both collapse to 1
 * identically, whatever the source ellipticity.
 *
 * Adding (endpoint - sum of retained coefficients) at the first UNUSED order
 * of each parity leaves every retained coefficient untouched -- so the
 * accuracy of order N at small delta is preserved exactly, and nothing is
 * recomputed -- while making the truncation exact at delta = 1. Parity fixes
 * where it lands: Delta1 is odd and Delta2/Delta3 even, so at odd N the
 * corrections sit at N+2 and N+1 respectively, and the other way round at
 * even N. @N_ext = N + 2 is the order the caller allocates and solves to.
 * The tables stop at N, so both correction slots land in an all-zero region:
 * the writes below are assignments into empty space, not overwrites of a
 * partially computed coefficient.
 *
 * Equivalently -- and this is why the correction is a truncation of
 * something exact rather than an ad hoc patch -- it is the truncation of the
 * factored form
 *
 *   Delta(delta) = E + (1 - delta^2) H(delta),
 *
 * with E the endpoint value and H truncated. Telescoping
 * (1-d^2) sum_{m<=M} H_m d^2m with H_m = sum_{k<=m} (Delta - E)_2k
 * reproduces the retained coefficients plus a single leftover term
 * -H_M d^(2M+2), and -H_M = E - sum_{k<=M} Delta_2k is exactly the
 * correction added below. Pulling out an exact factor of (1 - delta^2) from
 * the delta-dependent part and truncating only the cofactor is the same
 * operation as pinning the endpoint.
 *
 * Two things this does NOT do. It corrects Delta2 = <x^2>, the coordinate
 * the solve consumes, not C_t directly; C_t = <x^2> - mu^2 then differs from
 * a directly corrected C_t beyond order N, but conservatively -- it comes
 * out larger, hence more comfortably feasible -- and both are exactly 0 at
 * delta = 1. And it restores admissibility, not accuracy: it leaves the
 * retained orders alone by construction, so it cannot improve them, and at
 * moderate delta it can be marginally worse than the uncorrected truncation.
 * What it guarantees is that the target stays a possible moment vector all
 * the way to the critical curve.
 */
static void
_tilted_series_target_series (const gdouble *tab_m, guint n_m,
                              const gdouble *tab_v, const gdouble *tab_w, guint n_v,
                              const gdouble *M, guint n_moments, guint N, guint N_ext,
                              gdouble *Delta1, gdouble *Delta2, gdouble *Delta3)
{
  const guint p_odd  = ((N % 2) == 1) ? N + 2 : N + 1;
  const guint p_even = ((N % 2) == 1) ? N + 1 : N + 2;
  gdouble e1         = 1.0;
  gdouble e2         = 1.0 - 0.5 * M[1];
  gdouble e3         = -0.5 * M[1];
  guint l, p;

  g_assert_cmpuint (p_odd, <=, N_ext);
  g_assert_cmpuint (p_even, <=, N_ext);

  series_zero (Delta1, N_ext);
  series_zero (Delta2, N_ext);
  series_zero (Delta3, N_ext);

  for (l = 0; l < n_m; l++)
  {
    const guint pl = 2 * l + 1;
    gdouble s      = 0.0;
    guint a;

    if (pl > N)
      break;

    for (a = 0; a < n_moments; a++)
      s += tab_m[l * n_moments + a] * M[a];

    Delta1[pl] = s;
  }

  /* Row 0 is the delta = 0 term, which the lambda(0) = 0 theorem handles
   * and _tilted_series_delta_tables() already zeroed; start past it. */
  for (l = 1; l < n_v; l++)
  {
    const guint pl = 2 * l;
    gdouble sv     = 0.0;
    gdouble sw     = 0.0;
    guint a;

    if (pl > N)
      break;

    for (a = 0; a < n_moments; a++)
    {
      sv += tab_v[l * n_moments + a] * M[a];
      sw += tab_w[l * n_moments + a] * M[a];
    }

    Delta2[pl] = sv;
    Delta3[pl] = sw;
  }

  for (p = 0; p <= N; p++)
  {
    e1 -= Delta1[p];
    e2 -= Delta2[p];
    e3 -= Delta3[p];
  }

  Delta1[p_odd]  = e1;
  Delta2[p_even] = e2;
  Delta3[p_even] = e3;
}

/*
 * ---- ln P_0 by a localized, log-space Gauss-Legendre quadrature ----
 *
 * P_0(R) = (1/2*pi*sn^2) * int_0^1 P_pop^NC(r) * exp(-(R-r)^2/2sn^2) *
 * I0_scaled(R*r/sn^2) dr, where P_pop^NC is NumCosmo's own radial-marginal
 * convention (nc_galaxy_shape_pop_eval_p(): the disc measure 2*pi*r already
 * folded in), NOT a 2D area density. Converting between the two is exactly
 * what turns the unscaled exp(-(R^2+r^2)/2sn^2)*I0(Rr/sn^2) of the 2D form
 * into exp(-(R-r)^2/2sn^2)*I0_scaled(Rr/sn^2) here (the exponents combine
 * because I0(z) = exp(z)*I0_scaled(z)): getting this conversion wrong
 * silently returns the density for the wrong radial convention, not an
 * error.
 *
 * Any quadrature error here is a per-galaxy additive constant in ln P and
 * cancels from every derivative in g (see the class doc comment), so the
 * requirement is a smooth, finite, good-enough-for-ln-L quadrature -- not a
 * shear-unbiased one. A single fixed-node Gauss-Legendre panel over the
 * *whole* [0,1] would need enormously many nodes to resolve a Gaussian of
 * width sn much smaller than 1, so the integration window is localized
 * around where the integrand's mass actually is:
 *   R<=1: r in [R-sn*sqrt(2*DELTA), R+sn*sqrt(2*DELTA)] cap [0,1]
 *         (a symmetric window around the interior peak at r=R)
 *   R>1:  r in [R-sqrt((R-1)^2+2*sn^2*DELTA), 1]
 *         (the integrand instead peaks at the r=1 endpoint, with decay
 *         length sn^2/(R-1) -- a symmetric window sized for the interior
 *         case would badly under-resolve this regime)
 * DELTA=40 puts e^{-DELTA} ~ 4e-18 of a unit Gaussian's mass outside the
 * window, negligible next to the quadrature's own node-count error.
 *
 * The dynamic range that forces a log-space combination lives entirely in
 * the exponential factors: at small sn and R far from the disc, individual
 * terms underflow to a genuine (not clamped) double-precision zero well
 * before the sum does, so accumulating the density itself, not its log,
 * returns exactly zero. Only the *exponent* needs that treatment, though:
 * the quadrature weight and the population density are positive O(1)
 * prefactors, so they are carried in linear space and only
 * exp(a_i - max a) is formed. That removes two of the three logarithms per
 * node with no loss of range -- the terms actually summed are still bounded
 * by 1 -- and the result agrees with the all-logarithms form to 1.8e-15
 * relative over the whole (R, sn) box.
 *
 * The Gauss-Legendre nodes are computed once for [-1,1] and affinely mapped
 * per call: the table depends only on the (compile-time) node count, so
 * allocating one per galaxy would be pure waste.
 */
#define NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES 64
#define NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_DELTA 40.0

/* Gauss-Legendre nodes and weights on [-1,1], built once. */
static gdouble _tilted_series_lnP0_x[NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES];
static gdouble _tilted_series_lnP0_w[NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES];

static void
_tilted_series_lnP0_nodes_init (void)
{
  static gsize init = 0;

  if (g_once_init_enter (&init))
  {
    const guint n_nodes                  = NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES;
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (n_nodes);
    guint i;

    for (i = 0; i < n_nodes; i++)
      gsl_integration_glfixed_point (-1.0, 1.0, i,
                                     &_tilted_series_lnP0_x[i], &_tilted_series_lnP0_w[i], table);

    gsl_integration_glfixed_table_free (table);
    g_once_init_leave (&init, 1);
  }
}

static gdouble
_tilted_series_ln_P0 (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, gdouble R, gdouble sn)
{
  const guint n_nodes = NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES;
  const gdouble sn2   = sn * sn;
  gdouble r_lo, r_hi, half, mid;
  GArray *r_arr;
  GArray *p_arr = NULL;
  gdouble *r_data;
  gdouble expo[NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES];
  gdouble max_expo = -G_MAXDOUBLE;
  gdouble sum      = 0.0;
  guint i;

  _tilted_series_lnP0_nodes_init ();

  if (R <= 1.0)
  {
    const gdouble half_width = sn * sqrt (2.0 * NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_DELTA);

    r_lo = MAX (0.0, R - half_width);
    r_hi = MIN (1.0, R + half_width);
  }
  else
  {
    const gdouble d          = R - 1.0;
    const gdouble half_width = sqrt (d * d + 2.0 * sn2 * NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_DELTA);

    r_lo = MAX (0.0, R - half_width);
    r_hi = 1.0;
  }

  half  = 0.5 * (r_hi - r_lo);
  mid   = 0.5 * (r_hi + r_lo);
  r_arr = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_nodes);
  g_array_set_size (r_arr, n_nodes);
  r_data = (gdouble *) r_arr->data;

  for (i = 0; i < n_nodes; i++)
    r_data[i] = mid + half * _tilted_series_lnP0_x[i];

  nc_galaxy_shape_pop_eval_p_array (pop, pop_data, r_arr, &p_arr);

  for (i = 0; i < n_nodes; i++)
  {
    const gdouble r  = r_data[i];
    const gdouble dr = R - r;
    const gdouble z  = R * r / sn2;

    expo[i]  = -0.5 * dr * dr / sn2 + log (gsl_sf_bessel_I0_scaled (z));
    max_expo = MAX (max_expo, expo[i]);
  }

  {
    const gdouble *p_data = (const gdouble *) p_arr->data;

    for (i = 0; i < n_nodes; i++)
      sum += _tilted_series_lnP0_w[i] * p_data[i] * exp (expo[i] - max_expo);
  }

  g_array_unref (r_arr);
  g_array_unref (p_arr);

  return max_expo + log (half * sum) - log (2.0 * M_PI * sn2);
}

struct _NcGalaxyShapeFactorTiltedSeries
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorTiltedSeriesPrivate
{
  /* N: the number of true moment orders retained -- the `trunc-order`
   * property. */
  guint trunc_order;

  /* N + 2: the order lambda(delta) is actually solved to. The two extra
   * orders carry only the endpoint corrections of
   * _tilted_series_target_series(), no new moment information, but the
   * solve runs at this order and its cost is set by it. Derived in
   * constructed(), alongside n_M. */
  guint solve_order;

  /* Target-series tables, shared with MomentSeries: same layout as
   * NcGalaxyShapeFactorMomentSeriesPrivate's own fields, but rewritten from
   * g-coefficients to delta-coefficients at construction by
   * _tilted_series_delta_tables(). */
  guint n_m;
  guint n_v;
  guint n_moments;
  gdouble *tab_m;
  gdouble *tab_v;
  gdouble *tab_w;

  /* Radial moments M_2k fetched per galaxy, M[k] = M_2k. The
   * generating-function solve and the shared tables want different numbers
   * of them, and which dominates depends on the order; n_M covers both. See
   * constructed() for the derivation. */
  guint n_M;

  guint64 pop_hash;

  /* Domain-guard bookkeeping; see the class doc comment. @strict_domain
   * restores fatal behaviour. @domain_error_count is the running number of
   * out-of-domain evaluations and @domain_warned the one-shot flag for the
   * warning accompanying the first of them; both are gint and touched only
   * through g_atomic_int_*, since a shared factor can be evaluated from
   * several walker threads at once and an undercount here would understate
   * exactly the diagnostic the caller is reading. Neither is reset
   * internally -- the caller owns that, through
   * nc_galaxy_shape_factor_tilted_series_reset_domain_error_count(). */
  gboolean strict_domain;
  gint domain_error_count;
  gint domain_warned;
} NcGalaxyShapeFactorTiltedSeriesPrivate;

/*
 * Per-galaxy scratch, split into two independently-validated groups since
 * they depend on different values -- lambda/W solve only from
 * (pop_hash, sn), never from R (_tilted_series_solve()'s Delta/M inputs
 * carry no R dependence at all), while ln_P0 additionally depends on R.
 * Rebuilding lambda/W on every R change alone (as a single combined cache
 * used to) is a wasted 3x3-solve-per-order re-run on any path where R varies
 * at fixed sn -- e.g. a multi-R scan at fixed catalogue noise, or
 * nc_galaxy_shape_factor_gen() sweeping epsilon_obs at fixed std_noise.
 *
 * Both groups are still refreshed when the population generation moved or a
 * new catalog row was read (mirrors MomentSeries' own two invalidation axes;
 * see ldata_read_row() below, which invalidates both -- a new row may move
 * sn, R, or both) or -- unlike MomentSeries -- when std_noise or the
 * observed radius changed, since (unlike MomentSeries' m/v/w) this cache
 * depends on both: nc_galaxy_shape_factor_gen() and
 * nc_galaxy_shape_factor_data_set() both write epsilon_obs_1/2 and std_noise
 * without invoking ldata_read_row, so pop_hash/row invalidation alone is not
 * enough here.
 *
 * The four series are sized to the SOLVE order, not to trunc-order: the
 * endpoint corrections live above trunc-order and the solve carries them.
 */
typedef struct _NcGalaxyShapeFactorTiltedSeriesLData
{
  /* Keyed on (pop_hash, sn) only. */
  gdouble *lam1; /* solve_order+1 */
  gdouble *lam2; /* solve_order+1 */
  gdouble *lam3; /* solve_order+1 */
  gdouble *Wser; /* solve_order+1 */
  gdouble lam_bound;
  gdouble sn_seen;
  guint64 pop_hash_seen_lw;
  gboolean lw_valid;

  /* Keyed on (pop_hash, sn, R2) in addition to the above. */
  gdouble ln_P0;
  gdouble R2_seen;
  guint64 pop_hash_seen_p0;
  gboolean p0_valid;
} NcGalaxyShapeFactorTiltedSeriesLData;

enum
{
  PROP_0,
  PROP_TRUNC_ORDER,
  PROP_STRICT_DOMAIN,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorTiltedSeries, nc_galaxy_shape_factor_tilted_series, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
nc_galaxy_shape_factor_tilted_series_init (NcGalaxyShapeFactorTiltedSeries *gsfts)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (gsfts);

  /* Mirrors the `trunc-order` pspec default. The property is
   * construct-only, so GObject applies that default here anyway; keeping
   * the two in step means _init() alone leaves a coherent object. */
  self->trunc_order = 5;
  self->solve_order = 7;
  self->n_m         = 0;
  self->n_v         = 0;
  self->n_moments   = 0;
  self->n_M         = 0;
  self->tab_m       = NULL;
  self->tab_v       = NULL;
  self->tab_w       = NULL;
  self->pop_hash    = 0;

  self->strict_domain      = FALSE;
  self->domain_error_count = 0;
  self->domain_warned      = 0;
}

static void
_nc_galaxy_shape_factor_tilted_series_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorTiltedSeries *gsfts              = NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (object);
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (gsfts);

  switch (prop_id)
  {
    case PROP_TRUNC_ORDER:
      self->trunc_order = g_value_get_uint (value);
      break;
    case PROP_STRICT_DOMAIN:
      self->strict_domain = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_tilted_series_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorTiltedSeries *gsfts              = NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (object);
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (gsfts);

  switch (prop_id)
  {
    case PROP_TRUNC_ORDER:
      g_value_set_uint (value, self->trunc_order);
      break;
    case PROP_STRICT_DOMAIN:
      g_value_set_boolean (value, self->strict_domain);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_tilted_series_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_tilted_series_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactor *gsf                            = NC_GALAXY_SHAPE_FACTOR (object);
    NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (object));
    const NcGalaxyWLObsEllipConv ellip_conv             = nc_galaxy_shape_factor_get_ellip_conv (gsf);

    /* Everything derived from trunc-order is derived here, in one place:
     * the property is construct-only, so this runs exactly once and after
     * the property has been set. */
    self->solve_order = self->trunc_order + 2;

    /* Same population-independent, build-once-at-construction table this
     * class shares with MomentSeries -- see the class doc comment and
     * nc_galaxy_shape_factor_moment_series_private.h. It is built at
     * trunc-order, not at solve_order: trunc-order is by definition the
     * number of true moment orders, and the two orders above it carry only
     * the endpoint corrections, which are not table-derived. */
    _nc_galaxy_shape_factor_moment_series_build_tables (self->trunc_order, ellip_conv,
                                                        &self->n_m, &self->n_v, &self->n_moments,
                                                        &self->tab_m, &self->tab_v, &self->tab_w);

    /* The shared build works in g; this class orders everything in delta.
     * Rewrite the tables once, here, so no galaxy ever pays for it. */
    _tilted_series_delta_tables (self->trunc_order, self->n_m, self->n_v, self->n_moments,
                                 self->tab_m, self->tab_v, self->tab_w);

    /* Radial moments M_2k needed, as the max of three requirements:
     *
     *  - the shared tables contract against M[0..n_moments-1];
     *  - _tilted_series_Stab_new() asks iso_moment() for indices up to
     *    floor(solve_order/2) + 1, so it needs solve_order/2 + 2 entries;
     *  - _tilted_series_solve() unconditionally reads M[1] and M[2] to
     *    build the closed-form V = Cov_{P_0}(T) regardless of order, so
     *    n_M must be >= 3 even at trunc-order = 1, where neither of the
     *    first two reaches 3. */
    self->n_M = MAX (MAX (self->n_moments, self->solve_order / 2 + 2), 3);
  }
}

static void
_nc_galaxy_shape_factor_tilted_series_finalize (GObject *object)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (object));

  g_clear_pointer (&self->tab_m, g_free);
  g_clear_pointer (&self->tab_v, g_free);
  g_clear_pointer (&self->tab_w, g_free);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_tilted_series_parent_class)->finalize (object);
}

static void
_nc_galaxy_shape_factor_tilted_series_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorTiltedSeriesLData *ldata = (NcGalaxyShapeFactorTiltedSeriesLData *) p;

  g_free (ldata->lam1);
  g_free (ldata->lam2);
  g_free (ldata->lam3);
  g_free (ldata->Wser);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_tilted_series_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

/* A per-galaxy population reads its moments from the catalog row, so a new
 * row invalidates the cache without any model pkey moving -- mirrors
 * MomentSeries' own ldata_read_row(). This alone is not sufficient here (see
 * NcGalaxyShapeFactorTiltedSeriesLData's own comment on the sn/R2
 * value-keying); it is kept anyway since it is the axis that also
 * invalidates pop_data->e_rms itself. */
static void
_nc_galaxy_shape_factor_tilted_series_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorTiltedSeriesLData *ldata = (NcGalaxyShapeFactorTiltedSeriesLData *) data->ldata;

  ldata->lw_valid = FALSE;
  ldata->p0_valid = FALSE;
}

static void
_nc_galaxy_shape_factor_tilted_series_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_tilted_series_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (gsf));
  NcGalaxyShapeFactorTiltedSeriesLData *ldata         = g_new0 (NcGalaxyShapeFactorTiltedSeriesLData, 1);
  const guint N_ext                                   = self->solve_order;

  /* g_new0 leaves @lw_valid/@p0_valid FALSE, so the first evaluation
   * populates both. */
  ldata->lam1 = g_new0 (gdouble, N_ext + 1);
  ldata->lam2 = g_new0 (gdouble, N_ext + 1);
  ldata->lam3 = g_new0 (gdouble, N_ext + 1);
  ldata->Wser = g_new0 (gdouble, N_ext + 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_tilted_series_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_tilted_series_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_tilted_series_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_tilted_series_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_tilted_series_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (gsf));

  /* No capability gate, as MomentSeries: this class needs only
   * nc_galaxy_shape_pop_moment_2k() and nc_galaxy_shape_pop_eval_p_array(),
   * which every NcGalaxyShapePop provides, so it works with every
   * population, NcGalaxyShapePopBeta included. */
  self->pop_hash = nc_galaxy_shape_factor_get_pop_hash (gsf);
}

/*
 * Refreshes this galaxy's cached {lambda, W, lambda bound} and {ln P_0}
 * independently (see NcGalaxyShapeFactorTiltedSeriesLData's own comment):
 * the former on population generation, catalog row, or std_noise changes,
 * the latter additionally on the observed radius.
 */
static inline void
_nc_galaxy_shape_factor_tilted_series_peek_coeffs (NcGalaxyShapeFactorTiltedSeriesPrivate * const self,
                                                   NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                   const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                                   const gdouble **lam1_out, const gdouble **lam2_out,
                                                   const gdouble **lam3_out, const gdouble **Wser_out,
                                                   gdouble *ln_P0_out, gdouble *lam_bound_out)
{
  NcGalaxyShapeFactorTiltedSeriesLData *ldata = (NcGalaxyShapeFactorTiltedSeriesLData *) data->ldata;
  const gdouble sn                            = data->std_noise;

  /* R = |epsilon_obs| is rotation-invariant, so it is read from the
   * epsilon_obs_1/2 ARGUMENTS eval_marginal()/eval_ln_marginal() receive
   * (the contract every NcGalaxyShapeFactor subclass evaluates against --
   * see e.g. MomentSeries' own _eval()), not from data->epsilon_obs_1/2:
   * the two coincide on the fixed-nodes pipeline path (which hoists et/ex
   * from @data before calling), but nothing guarantees it in general
   * (nc_galaxy_shape_factor_integ()'s own integrand always passes
   * data->epsilon_obs_1/2 through unchanged, but a caller may not). */
  const gdouble R2 = gsl_pow_2 (epsilon_obs_1) + gsl_pow_2 (epsilon_obs_2);

  /* lambda/W depend only on (pop_hash, sn) -- ln_P0 additionally depends on
   * R -- so the two groups are validated and rebuilt independently (this
   * struct's own doc comment): an R-only change (fixed sn) never re-runs the
   * 3x3-solve-per-order loop. */
  if (G_UNLIKELY (!ldata->lw_valid || (ldata->pop_hash_seen_lw != self->pop_hash) ||
                  (ldata->sn_seen != sn)))
  {
    const guint N     = self->trunc_order;
    const guint N_ext = self->solve_order;
    const gdouble sn2 = sn * sn;
    gdouble *M        = g_new (gdouble, self->n_M);
    gdouble *Delta1   = g_new (gdouble, N_ext + 1);
    gdouble *Delta2   = g_new (gdouble, N_ext + 1);
    gdouble *Delta3   = g_new (gdouble, N_ext + 1);
    guint a;

    for (a = 0; a < self->n_M; a++)
      M[a] = nc_galaxy_shape_pop_moment_2k (pop, data->pop_data, a);

    _tilted_series_target_series (self->tab_m, self->n_m, self->tab_v, self->tab_w, self->n_v,
                                  M, self->n_moments, N, N_ext, Delta1, Delta2, Delta3);

    _tilted_series_solve (N_ext, M, self->n_M, sn2, Delta1, Delta2, Delta3,
                          ldata->lam1, ldata->lam2, ldata->lam3, ldata->Wser);

    ldata->lam_bound = 1.0 / (2.0 * sn2);

    g_free (M);
    g_free (Delta1);
    g_free (Delta2);
    g_free (Delta3);

    ldata->pop_hash_seen_lw = self->pop_hash;
    ldata->sn_seen          = sn;
    ldata->lw_valid         = TRUE;

    /* sn moved, so any previously cached ln_P0 (computed against the old
     * sn) is stale too, even if R2 did not change. */
    ldata->p0_valid = FALSE;
  }

  if (G_UNLIKELY (!ldata->p0_valid || (ldata->pop_hash_seen_p0 != self->pop_hash) ||
                  (ldata->R2_seen != R2)))
  {
    ldata->ln_P0 = _tilted_series_ln_P0 (pop, data->pop_data, sqrt (R2), sn);

    ldata->pop_hash_seen_p0 = self->pop_hash;
    ldata->R2_seen          = R2;
    ldata->p0_valid         = TRUE;
  }

  *lam1_out      = ldata->lam1;
  *lam2_out      = ldata->lam2;
  *lam3_out      = ldata->lam3;
  *Wser_out      = ldata->Wser;
  *ln_P0_out     = ldata->ln_P0;
  *lam_bound_out = ldata->lam_bound;
}

/*
 * Gauge-fixes (g, eps_obs) together by -arg(g) (exact, the same rotation
 * MomentSeries' own _eval uses), forms the ordering variable
 * delta = 2|g|/(1+|g|^2), Horner-evaluates lambda_1 (odd powers only) and
 * lambda_2/lambda_3/W (even powers only) in u = delta^2, applies the two
 * sanity guards, and returns
 * ln P_0 + lambda_1 x + lambda_2 x^2 + lambda_3 y^2 - W.
 *
 * Both eval_marginal and eval_ln_marginal route through here (@want_log
 * picks whether the final exp() runs), since
 * nc_galaxy_shape_factor_eval_at_nodes() calls eval_marginal on its hot
 * path -- see the class doc comment on why this class is not cheaper than a
 * Gaussian evaluation.
 *
 * The Horner bounds run to the SOLVE order, so the endpoint-correction
 * orders are evaluated along with the true moment orders: they are what
 * makes the model finite and admissible as delta approaches 1.
 */
static gdouble
_nc_galaxy_shape_factor_tilted_series_eval (NcGalaxyShapeFactorTiltedSeriesPrivate * const self,
                                            const gdouble *lam1, const gdouble *lam2, const gdouble *lam3,
                                            const gdouble *Wser, const gdouble ln_P0, const gdouble lam_bound,
                                            const gdouble g_1, const gdouble g_2,
                                            const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                            const gboolean want_log)
{
  const gdouble g_mag  = hypot (g_1, g_2);
  const gdouble cos_pg = (g_mag > 0.0) ? g_1 / g_mag : 1.0;
  const gdouble sin_pg = (g_mag > 0.0) ? g_2 / g_mag : 0.0;
  const gdouble x      = epsilon_obs_1 * cos_pg + epsilon_obs_2 * sin_pg;
  const gdouble y      = -epsilon_obs_1 * sin_pg + epsilon_obs_2 * cos_pg;
  const guint N        = self->solve_order;

  /* The ordering variable. Bounded by 1, invariant under g -> 1/g, and odd
   * in g, so the odd/even split below is the parity of the coefficients in
   * delta just as it was in g. */
  const gdouble d      = 2.0 * g_mag / (1.0 + g_mag * g_mag);
  const gdouble u      = d * d;
  const gint kmax_odd  = (gint) ((N - 1) / 2);
  const gint kmax_even = (gint) (N / 2);
  gdouble l1, l2, l3, Wg, lnP;
  gint k;

  l1 = lam1[2 * kmax_odd + 1];

  for (k = kmax_odd - 1; k >= 0; k--)
    l1 = l1 * u + lam1[2 * k + 1];

  l1 *= d;

  l2 = lam2[2 * kmax_even];
  l3 = lam3[2 * kmax_even];
  Wg = Wser[2 * kmax_even];

  for (k = kmax_even - 1; k >= 0; k--)
  {
    l2 = l2 * u + lam2[2 * k];
    l3 = l3 * u + lam3[2 * k];
    Wg = Wg * u + Wser[2 * k];
  }

  lnP = ln_P0 + l1 * x + l2 * x * x + l3 * y * y - Wg;

  /* Two independent sanity tests, in one branch.
   *
   * The first is the natural domain of the generating function Z(lambda),
   * outside which the noise integral defining it diverges. The second is an
   * exact ceiling: the marginal is a convolution of a probability density
   * with the noise kernel, so P <= max(kernel) = 1/(2 pi sn^2) for every g
   * and every eps_obs, i.e. lnP <= log(lam_bound / pi). The model is an
   * approximation and can sit a little above that, so the test trips only
   * far above it: this is a "returned nonsense" tripwire, not an accuracy
   * test.
   *
   * The ceiling is the one that catches the failure mode a domain test
   * cannot see. A truncated lambda far too NEGATIVE makes P a spurious
   * spike, and no test on max(lambda_2, lambda_3) notices, because the
   * exact lambda_2 and lambda_3 are themselves negative. Testing |lambda|
   * instead is not an option: the exact |lambda_3| comes close to lam_bound
   * at small sn, so a two-sided bound would fire on legitimate
   * configurations.
   */
  if (G_UNLIKELY ((MAX (l2, l3) >= lam_bound) ||
                  (lnP > log (lam_bound / M_PI) + 20.0)))
  {
    if (G_UNLIKELY (self->strict_domain))
      g_error ("NcGalaxyShapeFactorTiltedSeries: the tilt left its sanity "
               "bounds (lambda_2=%g, lambda_3=%g, bound=%g, lnP=%g, "
               "ceiling=%g) at trunc-order=%u, |g|=%g, delta=%g. The series "
               "is ordered in delta and its target is pinned at delta=1, so "
               "this is a truncation artefact rather than a reachable "
               "regime: either the target moment series handed to the solve "
               "was not the moment vector of any distribution -- check the "
               "population parameters -- or trunc-order is wrong for this "
               "population in one of two opposite ways: too low, leaving the "
               "target far from its endpoint and lambda_2 over the bound near "
               "|g|=1, or too high, the excess over the ceiling growing with "
               "order on narrow, low-noise populations.",
               l2, l3, lam_bound, lnP, log (lam_bound / M_PI),
               self->trunc_order, g_mag, d);

    /* One warning per instance: the sampler can revisit this region for
     * many galaxies over many walkers, and the count below is the number
     * that matters, not one line per occurrence. */
    if (G_UNLIKELY (g_atomic_int_compare_and_exchange (&self->domain_warned, 0, 1)))
      g_warning ("NcGalaxyShapeFactorTiltedSeries: the tilt left its sanity "
                 "bounds (lambda_2=%g, lambda_3=%g, bound=%g, lnP=%g, "
                 "ceiling=%g) at trunc-order=%u, |g|=%g, delta=%g; returning "
                 "zero probability, which routes into the caller's "
                 "NC_GALAXY_LOW_PROB path. The series is ordered in delta and "
                 "its target is pinned at delta=1, so this is a truncation "
                 "artefact rather than a reachable regime: check the "
                 "population parameters, and see this class' description "
                 "before changing trunc-order -- raising it clears a lambda_2 "
                 "over the bound, but worsens an excess over the ceiling. "
                 "Warned once per instance -- read the running total with "
                 "nc_galaxy_shape_factor_tilted_series_get_domain_error_count().",
                 l2, l3, lam_bound, lnP, log (lam_bound / M_PI),
                 self->trunc_order, g_mag, d);

    g_atomic_int_inc (&self->domain_error_count);

    /* exp(GSL_NEGINF) == 0.0 exactly, so the two hooks stay consistent. */
    return want_log ? GSL_NEGINF : 0.0;
  }

  return want_log ? lnP : exp (lnP);
}

static gdouble
_nc_galaxy_shape_factor_tilted_series_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (gsf));
  const gdouble *lam1, *lam2, *lam3, *Wser;
  gdouble ln_P0, lam_bound;

  _nc_galaxy_shape_factor_tilted_series_peek_coeffs (self, pop, data, epsilon_obs_1, epsilon_obs_2, &lam1, &lam2, &lam3, &Wser, &ln_P0, &lam_bound);

  return _nc_galaxy_shape_factor_tilted_series_eval (self, lam1, lam2, lam3, Wser, ln_P0, lam_bound,
                                                     g_1, g_2, epsilon_obs_1, epsilon_obs_2, FALSE);
}

static gdouble
_nc_galaxy_shape_factor_tilted_series_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES (gsf));
  const gdouble *lam1, *lam2, *lam3, *Wser;
  gdouble ln_P0, lam_bound;

  _nc_galaxy_shape_factor_tilted_series_peek_coeffs (self, pop, data, epsilon_obs_1, epsilon_obs_2, &lam1, &lam2, &lam3, &Wser, &ln_P0, &lam_bound);

  return _nc_galaxy_shape_factor_tilted_series_eval (self, lam1, lam2, lam3, Wser, ln_P0, lam_bound,
                                                     g_1, g_2, epsilon_obs_1, epsilon_obs_2, TRUE);
}

static void
nc_galaxy_shape_factor_tilted_series_class_init (NcGalaxyShapeFactorTiltedSeriesClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_tilted_series_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_tilted_series_get_property;
  object_class->constructed  = &_nc_galaxy_shape_factor_tilted_series_constructed;
  object_class->finalize     = &_nc_galaxy_shape_factor_tilted_series_finalize;

  /**
   * NcGalaxyShapeFactorTiltedSeries:trunc-order:
   *
   * Number of true moment orders $N$ retained in the $\delta$-series for
   * $\lambda$.
   *
   * The internal solve runs two orders higher, at $N+2$: those two orders
   * carry the endpoint corrections that keep the target admissible at the
   * critical curve (see this class' description) and no new moment
   * information. The cost of the solve is set by the larger order.
   *
   * Default 5, matching #NcGalaxyShapeFactorMomentSeries' own default: a
   * bias comparison across $N=5/7/9$ found the tilt series' remaining bias
   * numerically negligible at every order once solved correctly, so $N$
   * only buys back a fraction of a percent of calibration in the hardest
   * (small-$\sigma_\nu$) corner of the catalogue box
   * (`docs/theory/wl_shape_factor_history.md`). Whether raising it helps
   * when a sanity guard fires depends on which guard and why -- see this
   * class' description. Independent
   * of `MomentSeries`' own `trunc-order`, even though this class reuses
   * `MomentSeries`' target-series build at the same $N$.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRUNC_ORDER,
                                   g_param_spec_uint ("trunc-order",
                                                      "Truncation order",
                                                      "Number of true moment orders N retained in the delta-series for lambda",
                                                      1, G_MAXUINT, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorTiltedSeries:strict-domain:
   *
   * Whether leaving the evaluation-time sanity bounds -- $Z(\lambda)$'s
   * natural domain $\max(\lambda_2,\lambda_3) < 1/2\sigma_\nu^2$, or the
   * exact ceiling on the marginal -- is fatal.
   *
   * %FALSE (the default) returns zero probability and counts the
   * occurrence, so a sampler that wanders past $|g|=1$ is pushed back out
   * by #NcDataClusterWLFactor's `NC_GALAXY_LOW_PROB` penalty instead of
   * killing the process mid-chain. %TRUE aborts, which is useful in any
   * batch job that would rather fail than quietly penalise.
   *
   * Deliberately NOT %G_PARAM_CONSTRUCT: that flag would have GObject write
   * the pspec default over whatever _init() set, at construction time, for
   * every instance.
   */
  g_object_class_install_property (object_class,
                                   PROP_STRICT_DOMAIN,
                                   g_param_spec_boolean ("strict-domain",
                                                         "Strict domain",
                                                         "Abort instead of returning zero probability when the tilt leaves its sanity bounds",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_tilted_series_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_tilted_series_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_tilted_series_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_tilted_series_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_tilted_series_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 * @trunc_order: number of true moment orders $N$ retained in the
 * $\delta$-series for $\lambda$, $N\ge1$
 *
 * Creates a new #NcGalaxyShapeFactorTiltedSeries.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorTiltedSeries.
 */
NcGalaxyShapeFactorTiltedSeries *
nc_galaxy_shape_factor_tilted_series_new (NcGalaxyWLObsEllipConv ellip_conv, guint trunc_order)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_TILTED_SERIES,
                       "ellip-conv", ellip_conv,
                       "trunc-order", trunc_order,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_tilted_series_ref:
 * @gsfts: a #NcGalaxyShapeFactorTiltedSeries
 *
 * Increases the reference count of @gsfts by one.
 *
 * Returns: (transfer full): @gsfts.
 */
NcGalaxyShapeFactorTiltedSeries *
nc_galaxy_shape_factor_tilted_series_ref (NcGalaxyShapeFactorTiltedSeries *gsfts)
{
  return g_object_ref (gsfts);
}

/**
 * nc_galaxy_shape_factor_tilted_series_free:
 * @gsfts: a #NcGalaxyShapeFactorTiltedSeries
 *
 * Decreases the reference count of @gsfts by one.
 *
 */
void
nc_galaxy_shape_factor_tilted_series_free (NcGalaxyShapeFactorTiltedSeries *gsfts)
{
  g_object_unref (gsfts);
}

/**
 * nc_galaxy_shape_factor_tilted_series_clear:
 * @gsfts: a #NcGalaxyShapeFactorTiltedSeries
 *
 * Decreases the reference count of *@gsfts by one, and sets the pointer
 * *@gsfts to NULL.
 *
 */
void
nc_galaxy_shape_factor_tilted_series_clear (NcGalaxyShapeFactorTiltedSeries **gsfts)
{
  g_clear_object (gsfts);
}

/**
 * nc_galaxy_shape_factor_tilted_series_get_domain_error_count:
 * @gsfts: a #NcGalaxyShapeFactorTiltedSeries
 *
 * Number of evaluations that left the tilt's sanity bounds since
 * construction, or since the last
 * nc_galaxy_shape_factor_tilted_series_reset_domain_error_count().
 *
 * Each one returned zero probability rather than a truncated-series value
 * (see this class' description), so a non-zero count means part of the
 * likelihood was replaced by #NcDataClusterWLFactor's `NC_GALAXY_LOW_PROB`
 * penalty. On a converged chain that is first a signal to restrict the
 * shear range -- in a cluster fit, to tighten the mass prior. See this
 * class' description before reaching for
 * #NcGalaxyShapeFactorTiltedSeries:trunc-order instead: it helps or hurts
 * depending on which bound was crossed.
 *
 * Returns: the number of out-of-bounds evaluations.
 */
guint
nc_galaxy_shape_factor_tilted_series_get_domain_error_count (NcGalaxyShapeFactorTiltedSeries *gsfts)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate *self;

  g_return_val_if_fail (NC_IS_GALAXY_SHAPE_FACTOR_TILTED_SERIES (gsfts), 0);

  self = nc_galaxy_shape_factor_tilted_series_get_instance_private (gsfts);

  return (guint) g_atomic_int_get (&self->domain_error_count);
}

/**
 * nc_galaxy_shape_factor_tilted_series_reset_domain_error_count:
 * @gsfts: a #NcGalaxyShapeFactorTiltedSeries
 *
 * Zeroes the counter read by
 * nc_galaxy_shape_factor_tilted_series_get_domain_error_count() and re-arms
 * the once-per-instance warning, so a caller can attribute out-of-bounds
 * evaluations to a single likelihood evaluation, chain segment, or fit.
 *
 */
void
nc_galaxy_shape_factor_tilted_series_reset_domain_error_count (NcGalaxyShapeFactorTiltedSeries *gsfts)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate *self;

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR_TILTED_SERIES (gsfts));

  self = nc_galaxy_shape_factor_tilted_series_get_instance_private (gsfts);

  g_atomic_int_set (&self->domain_error_count, 0);
  g_atomic_int_set (&self->domain_warned, 0);
}

