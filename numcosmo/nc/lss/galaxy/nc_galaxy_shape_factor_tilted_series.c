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
 * A closed-form, table-free replacement for #NcGalaxyShapeFactorMomentSeries'
 * Gaussian clamp. `MomentSeries` carries the *exact*-in-$g$ mean and
 * covariance of the marginal but represents the density itself by a
 * Gaussian, which has no fourth cumulant while the disc-truncated
 * population does: its posterior is provably unbiased for the shear point
 * estimate but its credible intervals are 15%-18% too narrow (see
 * `MOMENT_SERIES.md` sec. 11 and `TILT_SERIES.md`/`tilt_series.tex` at the
 * repository root, this class' own design notes). This class matches the
 * *same* exact moments but replaces the base measure: instead of the
 * maximum-entropy density relative to Lebesgue measure (a Gaussian), it
 * uses the exponential tilt of the exact zero-shear marginal
 * $P_0 = P_\mathrm{pop} * \mathcal{N}(\sigma_\nu)$,
 * $$
 * \ln P(\chi_\mathrm{obs}\mid g) = \ln P_0(|\chi_\mathrm{obs}|)
 *   + \lambda(g)\cdot T(\chi_\mathrm{obs}) - W(g),
 *   \qquad T = (x,\ x^2,\ y^2),
 * $$
 * with $\lambda(g)$ fixed by requiring $\mathbb{E}_\lambda[T]$ to match the
 * *same* exact-in-$g$ mean/covariance `MomentSeries` already computes.
 * Because the shear map fixes the unit circle, $P_0$ carries the disc's
 * hard edge and every non-Gaussian cumulant the truncation generates; only
 * the smooth shear-induced deformation is expanded, and the resulting
 * density is positive by construction (an exponential family), unlike the
 * $g$-ordered density expansions that go negative on the counter-shear
 * side.
 *
 * $\lambda(g)$ is solved as a convergent power series in $g$
 * (`trunc-order`, default 9 -- higher than `MomentSeries`' default of 5,
 * because the tilt series converges an order of magnitude more slowly than
 * the moment series it consumes, and the two truncation orders are set
 * independently), one fixed $3\times3$ linear solve per order: no Newton
 * iteration, no population-indexed interpolation table. This is possible
 * because $\mathbb{E}_{P_0}[x^2] = M_2/2+\sigma_\nu^2$ identically equals
 * the exact zero-shear covariance, so $\lambda(0)=\bm 0$ exactly and the
 * moment conditions can be expanded order by order (`TILT_SERIES.md`
 * sec. 3). Each order's linear solve uses the closed-form
 * $V=\mathrm{Cov}_{P_0}(T)$; the right-hand side at every order is obtained
 * from a noise-integrated generating function
 * $Z(\lambda)=\mathbb{E}_{P_0}[e^{\lambda\cdot T}]$, evaluated as a formal
 * power series in $g$ via the population's own radial moments
 * ($M_{2k}=$ nc_galaxy_shape_pop_moment_2k(), to $k\lesssim(N{+}2)/2$) --
 * an exact, arbitrary-order, population-generic route that never
 * hardcodes $\lambda^{(p)}$ as a symbolic rational function of a fixed set
 * of moments the way `MomentSeries`' own $m_j/v_j/w_j$ tables do (a
 * departure from `TILT_SERIES.md` sec. 8's stated production plan --
 * hardcoding would freeze both the truncation order and the population's
 * functional form; see `docs/theory/wl_shape_factor_history.md`).
 *
 * The target moment series ($\Delta^{(p)}$, tilt_series.tex sec. 5) is the
 * *same* $g$-series `MomentSeries` builds from the population's radial
 * moments, so this class shares that build
 * (nc_galaxy_shape_factor_moment_series_private.h) rather than duplicating
 * it: both ellipticity conventions come for free, since boundary
 * invariance holds for `TRACE_DET` too and both maps are the identity at
 * $g=0$.
 *
 * The per-galaxy constant $\ln P_0(|\chi_\mathrm{obs}|)$ is $g$-independent,
 * so any quadrature error in it is a per-galaxy *additive* constant in
 * $\ln P$: it cancels identically from $\partial_g\ln P$ and never biases
 * the shear point estimate or the sandwich variance. This is what licenses
 * evaluating it by one population-generic, localized, log-space
 * Gauss-Legendre quadrature (never a table, and never the Marcum
 * $Q$-function closed form of a disc-truncated Gaussian population --
 * unneeded because the pipeline's `NcGalaxyShapePopGaussLocal` has a
 * per-galaxy, not fitted, width, so this quantity is computed at most once
 * per galaxy; see #NcGalaxyShapeFactorMomentSeries' own note on
 * `nc_galaxy_shape_factor_update_data_pop()` for when that stops being
 * true).
 *
 * Unlike `MomentSeries`, both `eval_marginal()` and `eval_ln_marginal()`
 * are equally expensive here: the fixed-nodes evaluation path
 * (nc_galaxy_shape_factor_eval_at_nodes()) calls `eval_marginal`, so an
 * `exp()` is unavoidable on that path regardless of which hook the caller
 * nominally wants -- this class does not claim to be cheaper than a
 * Gaussian evaluation.
 *
 * The tilt parameter must stay inside $Z(\lambda)$'s natural domain,
 * $\max(\lambda_2,\lambda_3) < 1/2\sigma_\nu^2$ (`TILT_SERIES.md' sec. 7):
 * asserted at evaluation time, not silently clamped, exactly like
 * `MomentSeries`' own covariance-positivity guard -- but unlike that guard,
 * this one is not reachable at any physical shear/noise combination the
 * catalogue exercises (the margin is $>300\times$ everywhere; the closed
 * form for the second-order coefficients shows $V\succ0$ unconditionally).
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
 * ---- Formal power-series arithmetic in g, truncated at order N ----
 *
 * Every series below is a plain gdouble[N+1] array, coefficient of g^p at
 * index p. These are the generic building blocks for D1's
 * noise-integrated generating-function solve (see this file's class doc
 * comment and TILT_SERIES.md sec. 3); none of them is population- or
 * order-specific.
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

/* table[n] = base^n, n=0..maxn, each a length-(order+1) series. */
static gdouble **
series_pow_table_new (const gdouble *base, guint maxn, guint order)
{
  gdouble **table = g_new (gdouble *, maxn + 1);
  guint n;

  for (n = 0; n <= maxn; n++)
    table[n] = g_new0 (gdouble, order + 1);

  table[0][0] = 1.0;

  for (n = 1; n <= maxn; n++)
    series_mul (table[n], table[n - 1], base, order);

  return table;
}

static void
series_pow_table_free (gdouble **table, guint maxn)
{
  guint n;

  for (n = 0; n <= maxn; n++)
    g_free (table[n]);

  g_free (table);
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
 * parities of n -- tilt_series.tex eq. (2.13)'s isotropic-moment formula
 * uses both. */
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
 * <(Re chi_I)^i_pow (Im chi_I)^j_pow> under isotropy (A4), zero unless
 * both powers are even (tilt_series.tex eq. 2.13):
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
 * ---- D1: the noise-integrated generating-function solve ----
 *
 * Computes a2(g)=1-2*lambda_2(g)*sn2, a3(g)=1-2*lambda_3(g)*sn2 and their
 * reciprocals, then D_u(g)=E_chis[chis_x^u chis_y^v e^Q] for the four
 * (u,v) in {(0,0),(1,0),(2,0),(0,2)} needed by the t1/t2/t3 model moments,
 * with Q = (lambda_1 sx + lambda_2 sx^2)/a2 + lambda_3 sy^2/a3
 * (TILT_SERIES.md sec. "Solving for lambda(g) in closed form", this
 * file's own generating-function derivation). Expanding e^Q in
 * alpha=lambda_1/a2, beta=lambda_2/a2, gamma=lambda_3/a3 gives a finite
 * triple sum over (pp,qq,ss) with pp+2qq+2ss<=N, each term an isotropic
 * moment of chi_I (population-generic, from nc_galaxy_shape_pop_moment_2k())
 * times a fixed combinatorial factor -- no quadrature, no table.
 */
static void
_tilted_series_compute_D (const gdouble *lam1, const gdouble *lam2, const gdouble *lam3,
                          const gdouble *M, guint n_M, gdouble sn2, guint N,
                          gdouble *a2, gdouble *a3, gdouble *a2_inv, gdouble *a3_inv,
                          gdouble *D0, gdouble *D1v, gdouble *D2v, gdouble *D3v)
{
  gdouble *alpha = g_new (gdouble, N + 1);
  gdouble *beta  = g_new (gdouble, N + 1);
  gdouble *gamma = g_new (gdouble, N + 1);
  gdouble *tmp1  = g_new (gdouble, N + 1);
  gdouble *tmp2  = g_new (gdouble, N + 1);
  const guint half_N = N / 2;
  gdouble **alpha_pows, **beta_pows, **gamma_pows;
  guint pp, qq, ss, k;

  series_zero (a2, N);
  series_zero (a3, N);
  a2[0] = 1.0;
  a3[0] = 1.0;

  for (k = 1; k <= N; k++)
  {
    a2[k] = -2.0 * sn2 * lam2[k];
    a3[k] = -2.0 * sn2 * lam3[k];
  }

  series_recip (a2_inv, a2, N);
  series_recip (a3_inv, a3, N);

  series_mul (alpha, lam1, a2_inv, N);
  series_mul (beta, lam2, a2_inv, N);
  series_mul (gamma, lam3, a3_inv, N);

  alpha_pows = series_pow_table_new (alpha, N, N);
  beta_pows  = series_pow_table_new (beta, half_N, N);
  gamma_pows = series_pow_table_new (gamma, half_N, N);

  series_zero (D0, N);
  series_zero (D1v, N);
  series_zero (D2v, N);
  series_zero (D3v, N);

  for (pp = 0; pp <= N; pp++)
  {
    for (qq = 0; pp + 2 * qq <= N; qq++)
    {
      series_mul (tmp1, alpha_pows[pp], beta_pows[qq], N);

      for (ss = 0; pp + 2 * qq + 2 * ss <= N; ss++)
      {
        const guint sx_pow = pp + 2 * qq;
        const guint sy_pow = 2 * ss;
        const gdouble fac  = factorial_d (pp) * factorial_d (qq) * factorial_d (ss);
        const gdouble S0   = iso_moment (M, n_M, sx_pow, sy_pow);
        const gdouble S1   = iso_moment (M, n_M, sx_pow + 1, sy_pow);
        const gdouble S2   = iso_moment (M, n_M, sx_pow + 2, sy_pow);
        const gdouble S3   = iso_moment (M, n_M, sx_pow, sy_pow + 2);

        if ((S0 == 0.0) && (S1 == 0.0) && (S2 == 0.0) && (S3 == 0.0))
          continue;

        series_mul (tmp2, tmp1, gamma_pows[ss], N);

        if (S0 != 0.0)
          series_axpy (D0, tmp2, S0 / fac, N);

        if (S1 != 0.0)
          series_axpy (D1v, tmp2, S1 / fac, N);

        if (S2 != 0.0)
          series_axpy (D2v, tmp2, S2 / fac, N);

        if (S3 != 0.0)
          series_axpy (D3v, tmp2, S3 / fac, N);
      }
    }
  }

  series_pow_table_free (alpha_pows, N);
  series_pow_table_free (beta_pows, half_N);
  series_pow_table_free (gamma_pows, half_N);
  g_free (alpha);
  g_free (beta);
  g_free (gamma);
  g_free (tmp1);
  g_free (tmp2);
}

/*
 * Order-by-order solve of lambda(g) and the closed-form log-normaliser
 * W(g) (TILT_SERIES.md's recursion V*lambda^(p) = Delta^(p) - R^(p)).
 * @Delta1/@Delta2/@Delta3 are this galaxy's target moment series (this
 * class' own _tilted_series_target_series(), reusing MomentSeries' table
 * build); @lam1/@lam2/@lam3/@Wser (each length N+1) receive the result.
 *
 * V = [[A,0,0],[0,B,kappa],[0,kappa,B]] is the closed-form covariance of
 * T=(x,x^2,y^2) under P_0 (tilt_series.tex eq. 5.9); it is
 * positive-definite unconditionally (as sn->0, 2*kappa+A^2 ->
 * Var(|chi_I|^2)/4 > 0 and kappa+A^2 -> M_4/8 > 0, and sn^2 only raises A),
 * so the guard below is an internal invariant, not a reachable
 * user-facing failure mode (contrast MomentSeries' own covariance guard).
 *
 * After the order-by-order loop, one more pass recomputes a2, a3, D0 with
 * the now-complete lambda(g), and W(g) is read off in closed form,
 * W = -0.5 ln a2 - 0.5 ln a3 + lambda_1^2 sn2/(2 a2) + ln D0 -- "free"
 * given a2/a3/D0, rather than by integrating W'=t.lambda' term by term
 * (that identity is a test of this routine, not its code path; see
 * tests/python/nc/lss/galaxy/test_galaxy_shape_factor_tilted_series.py).
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

  gdouble *a2       = g_new (gdouble, N + 1);
  gdouble *a3       = g_new (gdouble, N + 1);
  gdouble *a2_inv   = g_new (gdouble, N + 1);
  gdouble *a3_inv   = g_new (gdouble, N + 1);
  gdouble *D0       = g_new (gdouble, N + 1);
  gdouble *D1v      = g_new (gdouble, N + 1);
  gdouble *D2v      = g_new (gdouble, N + 1);
  gdouble *D3v      = g_new (gdouble, N + 1);
  gdouble *D0_inv   = g_new (gdouble, N + 1);
  gdouble *r1s      = g_new (gdouble, N + 1);
  gdouble *r2s      = g_new (gdouble, N + 1);
  gdouble *r3s      = g_new (gdouble, N + 1);
  gdouble *tmp_a    = g_new (gdouble, N + 1);
  gdouble *tmp_b    = g_new (gdouble, N + 1);
  gdouble *t1_model = g_new (gdouble, N + 1);
  gdouble *t2_model = g_new (gdouble, N + 1);
  gdouble *t3_model = g_new (gdouble, N + 1);
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
     * still zero), which is exactly R^(p): the g^p coefficient of
     * nabla W evaluated at the trial series with lambda^(p) omitted. */
    _tilted_series_compute_D (lam1, lam2, lam3, M, n_M, sn2, N,
                              a2, a3, a2_inv, a3_inv, D0, D1v, D2v, D3v);

    series_recip (D0_inv, D0, N);
    series_mul (r1s, D1v, D0_inv, N);
    series_mul (r2s, D2v, D0_inv, N);
    series_mul (r3s, D3v, D0_inv, N);

    /* t1_model = (r1s + sn2*lam1) * a2_inv */
    series_copy (tmp_a, r1s, N);
    series_axpy (tmp_a, lam1, sn2, N);
    series_mul (t1_model, tmp_a, a2_inv, N);

    /* t2_model = (r2s + 2*sn2*lam1*r1s + sn2^2*lam1^2) * a2_inv^2 + sn2*a2_inv */
    series_mul (tmp_a, lam1, r1s, N);
    series_copy (tmp_b, r2s, N);
    series_axpy (tmp_b, tmp_a, 2.0 * sn2, N);
    series_mul (tmp_a, lam1, lam1, N);
    series_axpy (tmp_b, tmp_a, sn2 * sn2, N);
    series_mul (tmp_a, a2_inv, a2_inv, N);
    series_mul (t2_model, tmp_b, tmp_a, N);
    series_axpy (t2_model, a2_inv, sn2, N);

    /* t3_model = r3s * a3_inv^2 + sn2*a3_inv */
    series_mul (tmp_a, a3_inv, a3_inv, N);
    series_mul (t3_model, r3s, tmp_a, N);
    series_axpy (t3_model, a3_inv, sn2, N);

    rhs1 = Delta1[p] - t1_model[p];
    rhs2 = Delta2[p] - t2_model[p];
    rhs3 = Delta3[p] - t3_model[p];

    lam1[p] = rhs1 / A;
    lam2[p] = (B * rhs2 - kappa * rhs3) / det2;
    lam3[p] = (B * rhs3 - kappa * rhs2) / det2;
  }

  _tilted_series_compute_D (lam1, lam2, lam3, M, n_M, sn2, N,
                            a2, a3, a2_inv, a3_inv, D0, D1v, D2v, D3v);

  series_log (tmp_a, a2, N);
  series_zero (Wser, N);
  series_axpy (Wser, tmp_a, -0.5, N);

  series_log (tmp_a, a3, N);
  series_axpy (Wser, tmp_a, -0.5, N);

  series_mul (tmp_a, lam1, lam1, N);
  series_mul (tmp_b, tmp_a, a2_inv, N);
  series_axpy (Wser, tmp_b, 0.5 * sn2, N);

  series_log (tmp_a, D0, N);
  series_axpy (Wser, tmp_a, 1.0, N);

  g_free (a2);
  g_free (a3);
  g_free (a2_inv);
  g_free (a3_inv);
  g_free (D0);
  g_free (D1v);
  g_free (D2v);
  g_free (D3v);
  g_free (D0_inv);
  g_free (r1s);
  g_free (r2s);
  g_free (r3s);
  g_free (tmp_a);
  g_free (tmp_b);
  g_free (t1_model);
  g_free (t2_model);
  g_free (t3_model);
}

/*
 * ---- D2: target moment series, reusing MomentSeries' tables ----
 *
 * Delta^(p) (p=1..N) = the g^p coefficient of t(g) = (mu, C_t+mu^2, C_x),
 * exactly the exact-in-g moments MomentSeries already carries. Reading
 * off Delta1/Delta2/Delta3 needs no polynomial algebra beyond
 * MomentSeries' own tab_m/tab_v/tab_w contraction: mu(g) has only odd
 * powers (tab_m row l gives the g^{2l+1} coefficient), and Ct(g)+mu(g)^2
 * / Cx(g) both reduce to the *raw*, pre-mu^2-subtraction tab_v/tab_w
 * contraction -- the mu^2 term that MomentSeries itself subtracts to turn
 * <(Re S)^2> into Var(Re S) cancels exactly against the mu(g)^2 added back
 * in here (both are the self-convolution of the same mu(g) series), so
 * this class never performs that subtraction at all.
 */
static void
_tilted_series_target_series (const gdouble *tab_m, guint n_m,
                              const gdouble *tab_v, const gdouble *tab_w, guint n_v,
                              const gdouble *M, guint n_moments, guint N,
                              gdouble *Delta1, gdouble *Delta2, gdouble *Delta3)
{
  guint l;

  series_zero (Delta1, N);
  series_zero (Delta2, N);
  series_zero (Delta3, N);

  for (l = 0; l < n_m; l++)
  {
    const guint p = 2 * l + 1;
    gdouble s = 0.0;
    guint a;

    if (p > N)
      break;

    for (a = 0; a < n_moments; a++)
      s += tab_m[l * n_moments + a] * M[a];

    Delta1[p] = s;
  }

  for (l = 0; l < n_v; l++)
  {
    const guint p = 2 * l;
    gdouble sv = 0.0, sw = 0.0;
    guint a;

    if (p > N)
      break;

    if (p == 0)
      continue; /* t(0) is handled separately (Theorem: lambda(0)=0). */

    for (a = 0; a < n_moments; a++)
    {
      sv += tab_v[l * n_moments + a] * M[a];
      sw += tab_w[l * n_moments + a] * M[a];
    }

    Delta2[p] = sv;
    Delta3[p] = sw;
  }
}

/*
 * ---- D3: ln P_0 by a localized, log-space Gauss-Legendre quadrature ----
 *
 * P_0(R) = (1/2*pi*sn^2) * int_0^1 P_pop^NC(r) * exp(-(R-r)^2/2sn^2) *
 * I0_scaled(R*r/sn^2) dr, where P_pop^NC is NumCosmo's own radial-marginal
 * convention (nc_galaxy_shape_pop_eval_p(): the disc-measure 2*pi*r
 * already folded in), NOT the tex's 2D area density -- converting between
 * the two is exactly what turns the tex's un-scaled exp(-(R^2+r^2)/2sn^2)*
 * I0(Rr/sn^2) into exp(-(R-r)^2/2sn^2)*I0_scaled(Rr/sn^2) here (the
 * exponents combine because I0(z) = exp(z)*I0_scaled(z)): getting this
 * conversion wrong silently returns the density for the wrong radial
 * convention, not an error.
 *
 * Any quadrature error here is a per-galaxy additive constant in ln P and
 * cancels from every derivative in g (see this file's own class doc
 * comment), so the requirement is a smooth, finite, good-enough-for-ln-L
 * quadrature -- not a shear-unbiased one. A single fixed-node
 * Gauss-Legendre panel over the *whole* [0,1] would need enormously many
 * nodes to resolve a Gaussian of width sn much smaller than 1, so the
 * integration window is localized around where the integrand's mass
 * actually is:
 *   R<=1: r in [R-sn*sqrt(2*DELTA), R+sn*sqrt(2*DELTA)] cap [0,1]
 *         (a symmetric window around the interior peak at r=R)
 *   R>1:  r in [R-sqrt((R-1)^2+2*sn^2*DELTA), 1]
 *         (the integrand instead peaks at the r=1 endpoint, with decay
 *         length sn^2/(R-1) -- a symmetric window sized for the interior
 *         case would badly under-resolve this regime)
 * DELTA=40 puts e^{-DELTA} ~ 4e-18 of a unit Gaussian's mass outside the
 * window, negligible next to the quadrature's own node-count error.
 *
 * Node values are combined in log space with the maximum exponent
 * factored out first: at small sn and R far from the disc, individual
 * terms underflow to a genuine (not clamped) double-precision zero well
 * before the sum does (e.g. R=1.99, sn=0.0033 gives ln P_0 ~ -45000), so
 * accumulating the density itself, not its log, returns exactly zero.
 */
#define NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES 64
#define NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_DELTA 40.0

static gdouble
_tilted_series_ln_P0 (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, gdouble R, gdouble sn)
{
  const guint n_nodes = NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES;
  const gdouble sn2   = sn * sn;
  gdouble r_lo, r_hi;
  gsl_integration_glfixed_table *table;
  GArray *r_arr;
  GArray *p_arr = NULL;
  gdouble *r_data;
  gdouble weight[NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES];
  gdouble log_terms[NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_LNP0_NNODES];
  gdouble max_log = -G_MAXDOUBLE;
  gdouble sum     = 0.0;
  guint i;

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

  table = gsl_integration_glfixed_table_alloc (n_nodes);
  r_arr = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_nodes);
  g_array_set_size (r_arr, n_nodes);
  r_data = (gdouble *) r_arr->data;

  for (i = 0; i < n_nodes; i++)
    gsl_integration_glfixed_point (r_lo, r_hi, i, &r_data[i], &weight[i], table);

  nc_galaxy_shape_pop_eval_p_array (pop, pop_data, r_arr, &p_arr);

  {
    const gdouble *p_data = (const gdouble *) p_arr->data;

    for (i = 0; i < n_nodes; i++)
    {
      const gdouble r  = r_data[i];
      const gdouble dr = R - r;
      const gdouble z  = R * r / sn2;

      log_terms[i] = log (weight[i]) + log (p_data[i]) - 0.5 * dr * dr / sn2 + log (gsl_sf_bessel_I0_scaled (z));

      max_log = MAX (max_log, log_terms[i]);
    }
  }

  for (i = 0; i < n_nodes; i++)
    sum += exp (log_terms[i] - max_log);

  g_array_unref (r_arr);
  g_array_unref (p_arr);
  gsl_integration_glfixed_table_free (table);

  return max_log + log (sum) - log (2.0 * M_PI * sn2);
}

struct _NcGalaxyShapeFactorTiltedSeries
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorTiltedSeriesPrivate
{
  guint trunc_order; /* N */

  /* Target-series tables, shared with MomentSeries (D2): same layout as
   * NcGalaxyShapeFactorMomentSeriesPrivate's own fields. */
  guint n_m;
  guint n_v;
  guint n_moments;
  gdouble *tab_m;
  gdouble *tab_v;
  gdouble *tab_w;

  /* Moments needed by D1's generating-function solve alone can exceed
   * n_moments at low trunc-order and fall short of it at high trunc-order
   * (TILT_SERIES.md's derived bound is k<=(N+2)/2, MomentSeries' own table
   * needs more at large N): n_M = max(n_moments, N/2+2) covers both, and
   * M[0..n_M-1] is fetched once per galaxy build. */
  guint n_M;

  guint64 pop_hash;
} NcGalaxyShapeFactorTiltedSeriesPrivate;

/*
 * Per-galaxy scratch: lambda(g)'s three components, W(g), ln P_0 and the
 * domain-guard bound, refreshed when the population generation moved, a
 * new catalog row was read (mirrors MomentSeries' own two invalidation
 * axes) or -- unlike MomentSeries -- when std_noise or the observed radius
 * changed, since (unlike MomentSeries' m/v/w) this cache depends on both:
 * nc_galaxy_shape_factor_gen() and nc_galaxy_shape_factor_data_set() both
 * write epsilon_obs_1/2 and std_noise without invoking ldata_read_row, so
 * pop_hash/row invalidation alone is not enough here.
 */
typedef struct _NcGalaxyShapeFactorTiltedSeriesLData
{
  gdouble *lam1; /* N+1 */
  gdouble *lam2; /* N+1 */
  gdouble *lam3; /* N+1 */
  gdouble *Wser; /* N+1 */
  gdouble ln_P0;
  gdouble lam_bound;
  gdouble sn_seen;
  gdouble R2_seen;
  guint64 pop_hash_seen;
  gboolean valid;
} NcGalaxyShapeFactorTiltedSeriesLData;

enum
{
  PROP_0,
  PROP_TRUNC_ORDER,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorTiltedSeries, nc_galaxy_shape_factor_tilted_series, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
nc_galaxy_shape_factor_tilted_series_init (NcGalaxyShapeFactorTiltedSeries *gsfts)
{
  NcGalaxyShapeFactorTiltedSeriesPrivate * const self = nc_galaxy_shape_factor_tilted_series_get_instance_private (gsfts);

  self->trunc_order = 9;
  self->n_m         = 0;
  self->n_v         = 0;
  self->n_moments   = 0;
  self->n_M         = 0;
  self->tab_m       = NULL;
  self->tab_v       = NULL;
  self->tab_w       = NULL;
  self->pop_hash    = 0;
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

    /* Same population-independent, build-once-at-construction table this
     * class shares with MomentSeries -- see this class' own doc comment
     * and nc_galaxy_shape_factor_moment_series_private.h. */
    _nc_galaxy_shape_factor_moment_series_build_tables (self->trunc_order, ellip_conv,
                                                        &self->n_m, &self->n_v, &self->n_moments,
                                                        &self->tab_m, &self->tab_v, &self->tab_w);

    self->n_M = MAX (self->n_moments, self->trunc_order / 2 + 2);
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

/* A per-galaxy population (NcGalaxyShapePopGaussLocal) reads its moments
 * from the catalog row, so a new row invalidates the cache without any
 * model pkey moving -- mirrors MomentSeries' own ldata_read_row(). This
 * alone is not sufficient here (see NcGalaxyShapeFactorTiltedSeriesLData's
 * own comment on the sn/R2 value-keying); it is kept anyway since it is
 * the axis that also invalidates pop_data->e_rms itself. */
static void
_nc_galaxy_shape_factor_tilted_series_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorTiltedSeriesLData *ldata = (NcGalaxyShapeFactorTiltedSeriesLData *) data->ldata;

  ldata->valid = FALSE;
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
  const guint N = self->trunc_order;

  /* g_new0 leaves @valid FALSE, so the first evaluation populates it. */
  ldata->lam1 = g_new0 (gdouble, N + 1);
  ldata->lam2 = g_new0 (gdouble, N + 1);
  ldata->lam3 = g_new0 (gdouble, N + 1);
  ldata->Wser = g_new0 (gdouble, N + 1);

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
 * Refreshes this galaxy's cached {lambda(g), W(g), ln P_0, lambda bound}
 * when the population model generation moved, a new catalog row was read,
 * or (see NcGalaxyShapeFactorTiltedSeriesLData's own comment) std_noise or
 * the observed radius changed.
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
  const gdouble sn = data->std_noise;

  /* R = |epsilon_obs| is rotation-invariant, so it is read from the
   * epsilon_obs_1/2 ARGUMENTS eval_marginal()/eval_ln_marginal() receive
   * (the contract every NcGalaxyShapeFactor subclass evaluates against --
   * see e.g. MomentSeries' own _eval()), not from data->epsilon_obs_1/2:
   * the two coincide on the fixed-nodes pipeline path (which hoists
   * et/ex from @data before calling), but nothing guarantees it in
   * general (nc_galaxy_shape_factor_integ()'s own integrand always passes
   * data->epsilon_obs_1/2 through unchanged, but a caller may not). */
  const gdouble R2 = gsl_pow_2 (epsilon_obs_1) + gsl_pow_2 (epsilon_obs_2);

  if (G_UNLIKELY (!ldata->valid || (ldata->pop_hash_seen != self->pop_hash) ||
                  (ldata->sn_seen != sn) || (ldata->R2_seen != R2)))
  {
    const guint N     = self->trunc_order;
    const gdouble sn2 = sn * sn;
    gdouble *M        = g_new (gdouble, self->n_M);
    gdouble *Delta1   = g_new (gdouble, N + 1);
    gdouble *Delta2   = g_new (gdouble, N + 1);
    gdouble *Delta3   = g_new (gdouble, N + 1);
    guint a;

    for (a = 0; a < self->n_M; a++)
      M[a] = nc_galaxy_shape_pop_moment_2k (pop, data->pop_data, a);

    _tilted_series_target_series (self->tab_m, self->n_m, self->tab_v, self->tab_w, self->n_v,
                                  M, self->n_moments, N, Delta1, Delta2, Delta3);

    _tilted_series_solve (N, M, self->n_M, sn2, Delta1, Delta2, Delta3,
                          ldata->lam1, ldata->lam2, ldata->lam3, ldata->Wser);

    ldata->ln_P0     = _tilted_series_ln_P0 (pop, data->pop_data, sqrt (R2), sn);
    ldata->lam_bound = 1.0 / (2.0 * sn2);

    g_free (M);
    g_free (Delta1);
    g_free (Delta2);
    g_free (Delta3);

    ldata->pop_hash_seen = self->pop_hash;
    ldata->sn_seen        = sn;
    ldata->R2_seen        = R2;
    ldata->valid          = TRUE;
  }

  *lam1_out      = ldata->lam1;
  *lam2_out      = ldata->lam2;
  *lam3_out      = ldata->lam3;
  *Wser_out      = ldata->Wser;
  *ln_P0_out     = ldata->ln_P0;
  *lam_bound_out = ldata->lam_bound;
}

/*
 * Gauge-fixes (g,eps_obs) together by -arg(g) (exact, same rotation as
 * MomentSeries' own _eval), Horner-evaluates lambda_1 (odd powers of g
 * only), lambda_2/lambda_3/W (even powers only) in u=g^2, asserts the
 * domain guard (this class' own doc comment / TILT_SERIES.md sec. 7), and
 * returns ln P_0 + lambda_1 x + lambda_2 x^2 + lambda_3 y^2 - W. Both
 * eval_marginal and eval_ln_marginal route through here (want_log picks
 * whether the final exp() runs) since nc_galaxy_shape_factor_eval_at_nodes()
 * calls eval_marginal on its hot path -- see this file's own class doc
 * comment on why this class is not cheaper than a Gaussian evaluation.
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
  const guint N        = self->trunc_order;
  const gdouble u      = g_mag * g_mag;
  const gint kmax_odd  = (gint) ((N - 1) / 2);
  const gint kmax_even = (gint) (N / 2);
  gdouble l1, l2, l3, Wg, lnP;
  gint k;

  l1 = lam1[2 * kmax_odd + 1];

  for (k = kmax_odd - 1; k >= 0; k--)
    l1 = l1 * u + lam1[2 * k + 1];

  l1 *= g_mag;

  l2 = lam2[2 * kmax_even];
  l3 = lam3[2 * kmax_even];
  Wg = Wser[2 * kmax_even];

  for (k = kmax_even - 1; k >= 0; k--)
  {
    l2 = l2 * u + lam2[2 * k];
    l3 = l3 * u + lam3[2 * k];
    Wg = Wg * u + Wser[2 * k];
  }

  if (G_UNLIKELY (MAX (l2, l3) >= lam_bound))
    g_error ("NcGalaxyShapeFactorTiltedSeries: the tilt parameter left the "
             "natural domain (lambda_2=%g, lambda_3=%g, bound=%g) at "
             "trunc-order=%u, |g|=%g. This marks a truncation breakdown "
             "outside this order's valid range -- raise trunc-order or "
             "restrict the shear range.",
             l2, l3, lam_bound, N, g_mag);

  lnP = ln_P0 + l1 * x + l2 * x * x + l3 * y * y - Wg;

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
   * Truncation order $N$ of the $g$-power series for $\lambda(g)$. Default
   * 9, higher than #NcGalaxyShapeFactorMomentSeries' default of 5: the
   * tilt series converges roughly a decade slower per two orders than the
   * moment series it consumes (`TILT_SERIES.md` sec. 3, "Convergence"), so
   * the two truncation orders are independent knobs even though this
   * class reuses `MomentSeries`' own target-series build at the same $N$.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRUNC_ORDER,
                                   g_param_spec_uint ("trunc-order",
                                                      "Truncation order",
                                                      "Truncation order N of the g-power series for lambda(g)",
                                                      1, G_MAXUINT, 9,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_tilted_series_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_tilted_series_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_tilted_series_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_tilted_series_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_tilted_series_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 * @trunc_order: truncation order $N$ of the $g$-power series for
 * $\lambda(g)$, $N\ge1$
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
