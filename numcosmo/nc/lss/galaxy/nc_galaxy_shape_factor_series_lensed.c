/***************************************************************************
 *            nc_galaxy_shape_factor_series_lensed.c
 *
 *  Wed Jul 8 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_series_lensed.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxyShapeFactorSeriesLensed:
 *
 * Lensed-frame series evaluation of the intrinsic-ellipticity marginal, both
 * ellipticity conventions, arbitrary truncation order, Gaussian population
 * only (v1).
 *
 * Changes to the LENSED frame (same substitution #NcGalaxyShapeFactorQuad
 * uses) and Taylor-expands the *population* term in $g$, keeping the noise
 * kernel exact:
 *
 *   marginal(eps_obs) = Int d^2chi_L [P_pop(f_{-g}(chi_L)) |Jac(f_{-g})(chi_L)|]
 *                        * N_2(eps_obs - chi_L; sigma_noise^2)
 *
 * The bracketed term $F_g(\chi_L)$ is Taylor-expanded in $g$ at fixed
 * $\chi_L$; its coefficients scale with $1/\sigma_\mathrm{pop}^2$, and
 * $\sigma_\mathrm{pop}$ is a population/prior parameter this project
 * hard-constrains away from being pathologically small (see
 * docs/theory/wl_shape_factor_history.md for why this matters relative to
 * expanding the noise kernel instead).
 *
 * Both conventions' $\chi_I(\chi_L,g)$ and Jacobian have closed forms (no
 * finite differences anywhere):
 *
 * - TRACE_DET (eps): $\chi_I=(\chi_L-g)/(1-g\chi_L)$, holomorphic, so
 *   $\mathrm{Jac}=|f'(\chi_L,-g)|^2$; both closed-form power series in $g$
 *   with no recursion needed at all.
 * - TRACE (chi): $\chi_I(\chi_L,g)=f_\chi(\chi_L,-g)$ is a ratio of two
 *   QUADRATICS in $g$ ($\chi_L,\bar\chi_L$ fixed), so its Taylor series
 *   follows a 2-term power-series-division recursion; the chi map is NOT
 *   holomorphic in $\chi_L$ (involves $\bar\chi_L$ too), so the Jacobian
 *   needs the full real $2\times2$ determinant, itself a closed-form
 *   reciprocal-series-cubed construction. Both closed forms match the
 *   already-shipped
 *   nc_wl_ellipticity_apply_shear_inv_trace_c/nc_wl_ellipticity_lndet_jac_trace_c
 *   to machine precision (see dev-notes'
 *   verify_inverse_map_and_jacobian_match_production).
 *
 * $P_\mathrm{pop}=\exp(-|\chi_I|^2/2\sigma_\mathrm{pop}^2)$'s composition
 * with $\chi_I(\chi_L,g)$'s own series uses the standard exp-of-power-series
 * recursion, with the $g^0$ term ($\chi_I(g=0)=\chi_L\neq0$) factored out as
 * a scalar prefactor first.
 *
 * $\theta$-integral done EXACTLY via Jacobi-Anger
 * (ncm_laurent_series_jacobi_anger_reduce()); $\rho$-integral a fixed,
 * windowed Gauss-Legendre quadrature (see
 * #NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA). Per-galaxy
 * $J[0..\mathrm{order}]$ cached, keyed to $(R,\phi,\sigma_\mathrm{noise})$
 * AND nc_galaxy_shape_factor_get_pop_hash(), since these coefficients depend
 * on the population's own parameters (e.g. $\sigma_\mathrm{pop}$), not just
 * per-galaxy data.
 *
 * Gaussian population only in this version: uses
 * nc_galaxy_shape_pop_get_sigma() (the same capability-based accessor
 * #NcGalaxyShapeFactorVarAdd already relies on), checked once at prepare()
 * time, erroring clearly for any population that doesn't support it.
 * General-population support (Taylor-composing an arbitrary population's own
 * local Taylor series in $x=|\chi_I|^2$ with the closed-form $\chi_I(g)$
 * series) is future work.
 *
 * See <a href="../../theory/wl_shape_marginalization_series.html">A
 * Small-Shear Series Marginalization for the Shape Likelihood</a> and
 * `dev-notes/wl_shape_series_marginalization_derivation.py` sections 9-11
 * for the full derivation and symbolic/numeric verification this class
 * mirrors.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_series_lensed.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/algebra/ncm_laurent_series.h"
#include "ncm/core/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* ===========================================================================
 * Generic bivariate (power-series-in-g of Laurent series) helpers. Plain
 * NcmLaurentSeries** arrays of length order+1, every #NcmLaurentSeries drawn
 * from @arena instead of freshly allocated -- see the "arena" struct and
 * _arena_next() below. The exact sequence and sizes requested here are
 * driven purely by @order (fixed for an instance's whole life) and which
 * chi_taylor/jac_taylor pair is selected (also fixed), never by rho/R/phi
 * or any other per-node/per-galaxy value -- so @arena grows only during the
 * very first node ever computed for a given workspace, and is pure reuse
 * (zero allocation) on every call after that. See _compute_J for where
 * @arena's backing pool lives (a per-instance NcmMemoryPool of
 * pre-populated workspaces) and where its bump index gets reset to 0 once
 * per node.
 * ===========================================================================
 */

typedef struct _NcGalaxyShapeFactorSeriesLensedArena
{
  GPtrArray *pool; /* NcmLaurentSeries*, persists for the owning workspace's whole life */
  guint idx;       /* bump index, reset to 0 at the start of each node's computation */
} NcGalaxyShapeFactorSeriesLensedArena;

static NcmLaurentSeries *
_arena_next (NcGalaxyShapeFactorSeriesLensedArena *arena)
{
  NcmLaurentSeries *a;

  if (arena->idx < arena->pool->len)
  {
    a = g_ptr_array_index (arena->pool, arena->idx);
  }
  else
  {
    a = ncm_laurent_series_new (0, 0);
    g_ptr_array_add (arena->pool, a);
  }

  arena->idx++;

  return a;
}

static void
_bivar_conv (NcGalaxyShapeFactorSeriesLensedArena *arena, NcmLaurentSeries * const *a, NcmLaurentSeries * const *b, guint order, NcmLaurentSeries **out)
{
  guint m;

  for (m = 0; m <= order; m++)
  {
    NcmLaurentSeries *acc = _arena_next (arena);
    guint k;

    ncm_laurent_series_reset (acc, 0, 0);

    for (k = 0; k <= m; k++)
    {
      NcmLaurentSeries *term = _arena_next (arena);
      NcmLaurentSeries *acc2 = _arena_next (arena);

      ncm_laurent_series_conv_into (term, a[k], b[m - k]);
      ncm_laurent_series_add_into (acc2, acc, term, 1.0);
      acc = acc2;
    }

    out[m] = acc;
  }
}

static void
_bivar_scale (NcGalaxyShapeFactorSeriesLensedArena *arena, NcmLaurentSeries * const *a, complex double s, guint order, NcmLaurentSeries **out)
{
  guint m;

  for (m = 0; m <= order; m++)
  {
    NcmLaurentSeries *o = _arena_next (arena);

    ncm_laurent_series_scale_c_into (o, a[m], s);
    out[m] = o;
  }
}

static void
_bivar_conj (NcGalaxyShapeFactorSeriesLensedArena *arena, NcmLaurentSeries * const *a, guint order, NcmLaurentSeries **out)
{
  guint m;

  for (m = 0; m <= order; m++)
  {
    NcmLaurentSeries *o = _arena_next (arena);

    ncm_laurent_series_conj_into (o, a[m]);
    out[m] = o;
  }
}

/* ===========================================================================
 * Convention-specific closed-form chi_I(chi_L,g) and Jac(chi_L,g)
 * Taylor-in-g coefficients -- see the class docs above and
 * dev-notes/wl_shape_series_marginalization_derivation.py sections 10-11
 * (verify_chi_closed_form_pieces, verify_eps_closed_form_pieces) for the
 * derivation and numeric verification.
 * ===========================================================================
 */

typedef void (*NcGalaxyShapeFactorSeriesLensedChiTaylorFunc) (NcGalaxyShapeFactorSeriesLensedArena *arena, gdouble rho, guint order, NcmLaurentSeries **c);
typedef void (*NcGalaxyShapeFactorSeriesLensedJacTaylorFunc) (NcGalaxyShapeFactorSeriesLensedArena *arena, gdouble rho, guint order, NcmLaurentSeries **jac);

/* eps: chi_I(chi_L,g)=(chi_L-g)/(1-g*chi_L), closed form, no recursion:
 * c_0=chi_L={1:rho}; c_k=chi_L^(k-1)*(chi_L^2-1) = rho^(k+1)*w^(k+1) -
 * rho^(k-1)*w^(k-1) for k>=1 -- a 2-term Laurent series at every order. */
static void
_chi_taylor_eps (NcGalaxyShapeFactorSeriesLensedArena *arena, gdouble rho, guint order, NcmLaurentSeries **c)
{
  guint k;
  NcmLaurentSeries *c0 = _arena_next (arena);

  ncm_laurent_series_set_single_into (c0, 1, rho);
  c[0] = c0;

  for (k = 1; k <= order; k++)
  {
    gdouble rho_km1      = pow (rho, (gdouble) (k - 1));
    gdouble rho_kp1      = rho_km1 * rho * rho;
    NcmLaurentSeries *ck = _arena_next (arena);

    ncm_laurent_series_reset (ck, (gint) k - 1, (gint) k + 1);
    ncm_laurent_series_set_c (ck, (gint) k - 1, -rho_km1);
    ncm_laurent_series_set_c (ck, (gint) k + 1, rho_kp1);
    c[k] = ck;
  }
}

/* eps: f_eps'(chi_L,-g)=(1-g^2)/(1-g*chi_L)^2 Taylor-in-g coefficients --
 * closed form: 1/(1-g*chi_L)^2 = sum_n(n+1)*chi_L^n*g^n (literal binomial
 * series, single-term at harmonic n each), times (1-g^2). Jac=|that|^2 is
 * this series convolved with its own bivariate conjugate. */
static void
_jac_taylor_eps (NcGalaxyShapeFactorSeriesLensedArena *arena, gdouble rho, guint order, NcmLaurentSeries **jac)
{
  NcmLaurentSeries **binom      = g_new0 (NcmLaurentSeries *, order + 3);
  NcmLaurentSeries **deriv      = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **deriv_conj = g_new0 (NcmLaurentSeries *, order + 1);
  guint n;

  for (n = 0; n <= order + 2; n++)
  {
    NcmLaurentSeries *b = _arena_next (arena);

    ncm_laurent_series_set_single_into (b, (gint) n, (n + 1.0) * pow (rho, (gdouble) n));
    binom[n] = b;
  }

  for (n = 0; n <= order; n++)
  {
    if (n >= 2)
    {
      NcmLaurentSeries *d = _arena_next (arena);

      ncm_laurent_series_add_into (d, binom[n], binom[n - 2], -1.0);
      deriv[n] = d;
    }
    else
    {
      deriv[n] = binom[n];
    }
  }

  _bivar_conj (arena, deriv, order, deriv_conj);
  _bivar_conv (arena, deriv, deriv_conj, order, jac);

  g_free (binom);
  g_free (deriv);
  g_free (deriv_conj);
}

/* chi: chi_I(chi_L,g)=f_chi(chi_L,-g) is Num(g)/Den(g), ratio of two
 * quadratics in g (chi_L,chibar_L fixed): Num=chi_L-2g+chibar_L*g^2,
 * Den=1-(chi_L+chibar_L)*g+g^2 -- 2-term power-series-division recursion,
 * c_0=chi_L=n_0; c_k=n_k-d1*c_{k-1}-d2*c_{k-2} (k>=1), n_1=-2, n_2=chibar_L,
 * n_k=0 otherwise, d1=-(chi_L+chibar_L), d2=1 (so the d2 term is just
 * c_{k-2} itself, no convolution needed). */
static void
_chi_taylor_chi (NcGalaxyShapeFactorSeriesLensedArena *arena, gdouble rho, guint order, NcmLaurentSeries **c)
{
  NcmLaurentSeries *chi_L    = _arena_next (arena);
  NcmLaurentSeries *chibar_L = _arena_next (arena);
  NcmLaurentSeries *two_s    = _arena_next (arena);
  NcmLaurentSeries *d1       = _arena_next (arena);
  guint k;

  ncm_laurent_series_set_single_into (chi_L, 1, rho);
  ncm_laurent_series_set_single_into (chibar_L, -1, rho);
  ncm_laurent_series_add_into (two_s, chi_L, chibar_L, 1.0);
  ncm_laurent_series_scale_c_into (d1, two_s, -1.0);

  c[0] = chi_L;

  for (k = 1; k <= order; k++)
  {
    NcmLaurentSeries *base;
    NcmLaurentSeries *conv_term;
    NcmLaurentSeries *acc;

    if (k == 1)
    {
      base = _arena_next (arena);
      ncm_laurent_series_set_single_into (base, 0, -2.0);
    }
    else if (k == 2)
    {
      base = chibar_L;
    }
    else
    {
      base = _arena_next (arena);
      ncm_laurent_series_reset (base, 0, 0);
    }

    conv_term = _arena_next (arena);
    ncm_laurent_series_conv_into (conv_term, d1, c[k - 1]);

    acc = _arena_next (arena);
    ncm_laurent_series_add_into (acc, base, conv_term, -1.0);

    if (k >= 2)
    {
      NcmLaurentSeries *acc2 = _arena_next (arena);

      ncm_laurent_series_add_into (acc2, acc, c[k - 2], -1.0);
      acc = acc2;
    }

    c[k] = acc;
  }
}

/* chi: Jac(chi_L,g) = (1-g^2)^3 / D(g)^3, D(g) the same denominator as
 * _chi_taylor_chi (the chi map is not holomorphic in chi_L, so the full real
 * 2x2 Jacobian is needed -- derived by hand, confirmed to machine precision
 * against the already-shipped nc_wl_ellipticity_lndet_jac_trace_c, see the
 * class docs). */
static void
_jac_taylor_chi (NcGalaxyShapeFactorSeriesLensedArena *arena, gdouble rho, guint order, NcmLaurentSeries **jac)
{
  NcmLaurentSeries *chi_L    = _arena_next (arena);
  NcmLaurentSeries *chibar_L = _arena_next (arena);
  NcmLaurentSeries *two_s    = _arena_next (arena);
  NcmLaurentSeries *d1       = _arena_next (arena);
  NcmLaurentSeries **r       = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **r2      = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **r3      = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **num     = g_new0 (NcmLaurentSeries *, order + 1);
  guint k;

  ncm_laurent_series_set_single_into (chi_L, 1, rho);
  ncm_laurent_series_set_single_into (chibar_L, -1, rho);
  ncm_laurent_series_add_into (two_s, chi_L, chibar_L, 1.0);
  ncm_laurent_series_scale_c_into (d1, two_s, -1.0);

  /* r = 1/D(g), D(g)=1+d1*g+g^2 (d0=1, d2=1) */
  r[0] = _arena_next (arena);
  ncm_laurent_series_set_single_into (r[0], 0, 1.0);

  for (k = 1; k <= order; k++)
  {
    NcmLaurentSeries *conv_term = _arena_next (arena);
    NcmLaurentSeries *acc       = _arena_next (arena);

    ncm_laurent_series_conv_into (conv_term, d1, r[k - 1]);
    ncm_laurent_series_scale_c_into (acc, conv_term, -1.0);

    if (k >= 2)
    {
      NcmLaurentSeries *acc2 = _arena_next (arena);

      ncm_laurent_series_add_into (acc2, acc, r[k - 2], -1.0);
      acc = acc2;
    }

    r[k] = acc;
  }

  _bivar_conv (arena, r, r, order, r2);
  _bivar_conv (arena, r2, r, order, r3);

  /* (1-g^2)^3 = 1 - 3g^2 + 3g^4 - g^6 */
  for (k = 0; k <= order; k++)
  {
    gdouble coeff;
    NcmLaurentSeries *n;

    if (k == 0)
      coeff = 1.0;
    else if (k == 2)
      coeff = -3.0;
    else if (k == 4)
      coeff = 3.0;
    else if (k == 6)
      coeff = -1.0;
    else
      coeff = 0.0;

    n = _arena_next (arena);
    ncm_laurent_series_set_single_into (n, 0, coeff);
    num[k] = n;
  }

  _bivar_conv (arena, num, r3, order, jac);

  g_free (r);
  g_free (r2);
  g_free (r3);
  g_free (num);
}

/* ===========================================================================
 * Population composition (Gaussian only, convention-independent) and the
 * final F_g(chi_L) = P_pop(chi_I(chi_L,g)) * Jac(chi_L,g) product.
 * ===========================================================================
 */

/* exp(A(g)) as a g-power series, A possibly having a NONZERO g^0 term
 * (unlike the noise-side scheme's exp(Delta), where Delta(theta,0)=0
 * always) -- pull out exp(a0) as a scalar prefactor, apply the standard
 * c_m recursion to the g>=1 remainder. */
static void
_exp_series (NcGalaxyShapeFactorSeriesLensedArena *arena, NcmLaurentSeries * const *A, guint order, NcmLaurentSeries **out)
{
  complex double a0 = ncm_laurent_series_get_c (A[0], 0);
  complex double prefactor;
  NcmLaurentSeries **c = g_new0 (NcmLaurentSeries *, order + 1);
  guint m;

  g_assert_cmpint (ncm_laurent_series_hmin (A[0]), ==, 0);
  g_assert_cmpint (ncm_laurent_series_hmax (A[0]), ==, 0);

  prefactor = cexp (a0);

  c[0] = _arena_next (arena);
  ncm_laurent_series_set_single_into (c[0], 0, 1.0);

  for (m = 1; m <= order; m++)
  {
    NcmLaurentSeries *acc = _arena_next (arena);
    guint k;

    ncm_laurent_series_reset (acc, 0, 0);

    for (k = 1; k <= m; k++)
    {
      NcmLaurentSeries *term = _arena_next (arena);
      NcmLaurentSeries *acc2 = _arena_next (arena);

      ncm_laurent_series_conv_into (term, A[k], c[m - k]);
      ncm_laurent_series_add_into (acc2, acc, term, (gdouble) k);
      acc = acc2;
    }

    c[m] = _arena_next (arena);
    ncm_laurent_series_scale_c_into (c[m], acc, 1.0 / m);
  }

  for (m = 0; m <= order; m++)
  {
    NcmLaurentSeries *o = _arena_next (arena);

    ncm_laurent_series_scale_c_into (o, c[m], prefactor);
    out[m] = o;
  }

  g_free (c);
}

static void
_F_g_taylor (NcGalaxyShapeFactorSeriesLensedArena *arena, NcGalaxyShapeFactorSeriesLensedChiTaylorFunc chi_taylor, NcGalaxyShapeFactorSeriesLensedJacTaylorFunc jac_taylor,
             gdouble rho, gdouble sigma_pop, guint order, NcmLaurentSeries **F)
{
  NcmLaurentSeries **c        = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **cbar     = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **abs_sq   = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **exponent = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **pop      = g_new0 (NcmLaurentSeries *, order + 1);
  NcmLaurentSeries **jac      = g_new0 (NcmLaurentSeries *, order + 1);

  chi_taylor (arena, rho, order, c);
  _bivar_conj (arena, c, order, cbar);
  _bivar_conv (arena, c, cbar, order, abs_sq);
  _bivar_scale (arena, abs_sq, -1.0 / (2.0 * sigma_pop * sigma_pop), order, exponent);
  _exp_series (arena, exponent, order, pop);

  jac_taylor (arena, rho, order, jac);

  _bivar_conv (arena, pop, jac, order, F);

  g_free (c);
  g_free (cbar);
  g_free (abs_sq);
  g_free (exponent);
  g_free (pop);
  g_free (jac);
}

/* ===========================================================================
 * GObject boilerplate
 * ===========================================================================
 */

struct _NcGalaxyShapeFactorSeriesLensed
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorSeriesLensedPrivate
{
  NcGalaxyShapeFactorSeriesLensedChiTaylorFunc chi_taylor;
  NcGalaxyShapeFactorSeriesLensedJacTaylorFunc jac_taylor;

  guint trunc_order;
  guint n_nodes;

  /* Pool of per-computation workspaces (NcGalaxyShapeFactorSeriesLensedWS,
   * defined further down, right before _compute_J), each bundling a
   * reusable NcmLaurentSeries arena plus the rho_nodes/weights/GL-table
   * triple _compute_J needs. NcmMemoryPool (ncm_memory_pool.h) is this
   * project's own thread-safe, auto-growing checkout pool -- same pattern
   * as ncm_integral_get_workspace() and NcmStatsDistKDE's mp_eval_vars:
   * _compute_J checks out one workspace per call and returns it when done,
   * so concurrent calls (were this ever parallelized) get independent
   * workspaces rather than racing on shared state. */
  NcmMemoryPool *mp_ws;
} NcGalaxyShapeFactorSeriesLensedPrivate;

enum
{
  PROP_0,
  PROP_TRUNC_ORDER,
  PROP_N_NODES,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorSeriesLensed, nc_galaxy_shape_factor_series_lensed, NC_TYPE_GALAXY_SHAPE_FACTOR)

/* Defined near _compute_J, forward-declared here for _init()'s use. */
static gpointer _nc_galaxy_shape_factor_series_lensed_ws_new (gpointer userdata);
static void _nc_galaxy_shape_factor_series_lensed_ws_free (gpointer p);

static void
nc_galaxy_shape_factor_series_lensed_init (NcGalaxyShapeFactorSeriesLensed *gsfsl)
{
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);

  self->chi_taylor  = NULL;
  self->jac_taylor  = NULL;
  self->trunc_order = 4;
  self->n_nodes     = 60;

  /* Lazily populated (ncm_memory_pool_new itself allocates no workspace),
   * so this runs safely before trunc-order/n-nodes's CONSTRUCT properties
   * are set -- _ws_new (defined near _compute_J) reads them from @gsfsl's
   * own private data only when a workspace is actually first requested,
   * by which point construction has completed. */
  self->mp_ws = ncm_memory_pool_new (&_nc_galaxy_shape_factor_series_lensed_ws_new, gsfsl, &_nc_galaxy_shape_factor_series_lensed_ws_free);
}

static void
_nc_galaxy_shape_factor_series_lensed_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorSeriesLensed *gsfsl              = NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (object);
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);

  switch (prop_id)
  {
    case PROP_TRUNC_ORDER:
      self->trunc_order = g_value_get_uint (value);
      break;
    case PROP_N_NODES:
      self->n_nodes = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_series_lensed_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorSeriesLensed *gsfsl              = NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (object);
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);

  switch (prop_id)
  {
    case PROP_TRUNC_ORDER:
      g_value_set_uint (value, self->trunc_order);
      break;
    case PROP_N_NODES:
      g_value_set_uint (value, self->n_nodes);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_series_lensed_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_series_lensed_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactor *gsf                            = NC_GALAXY_SHAPE_FACTOR (object);
    NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (object));
    const NcGalaxyWLObsEllipConv ellip_conv             = nc_galaxy_shape_factor_get_ellip_conv (gsf);

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->chi_taylor = &_chi_taylor_chi;
        self->jac_taylor = &_jac_taylor_chi;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->chi_taylor = &_chi_taylor_eps;
        self->jac_taylor = &_jac_taylor_eps;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

/* Gaussian-only guard, checked once here (a property of the pop MODEL,
 * identical for every galaxy sharing it), matching
 * NcGalaxyShapeFactorVarAdd's own prepare() override exactly. */
static void
_nc_galaxy_shape_factor_series_lensed_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapePop *pop              = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  NcGalaxyShapePopData *tmp_pop_data = nc_galaxy_shape_pop_data_new (pop);
  const gboolean has_sigma           = tmp_pop_data->ldata_get_sigma != NULL;

  nc_galaxy_shape_pop_data_unref (tmp_pop_data);

  if (!has_sigma)
    g_error ("NcGalaxyShapeFactorSeriesLensed: this scheme's v1 Gaussian-only "
             "population composition requires a population supporting "
             "nc_galaxy_shape_pop_get_sigma(), got %s.", G_OBJECT_TYPE_NAME (pop));
}

static void
_nc_galaxy_shape_factor_series_lensed_finalize (GObject *object)
{
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (object));

  if (self->mp_ws != NULL)
  {
    ncm_memory_pool_free (self->mp_ws, TRUE);
    self->mp_ws = NULL;
  }

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_series_lensed_parent_class)->finalize (object);
}

typedef struct _NcGalaxyShapeFactorSeriesLensedData
{
  gboolean cache_valid;
  gdouble cached_R;
  gdouble cached_phi;
  gdouble cached_std_noise;
  guint64 cached_pop_hash;
  gdouble *J; /* size trunc_order+1 */
} NcGalaxyShapeFactorSeriesLensedData;

static void
_nc_galaxy_shape_factor_series_lensed_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorSeriesLensedData *ldata = (NcGalaxyShapeFactorSeriesLensedData *) p;

  g_free (ldata->J);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_series_lensed_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_series_lensed_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_series_lensed_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (gsf));
  NcGalaxyShapeFactorSeriesLensedData *ldata          = g_new0 (NcGalaxyShapeFactorSeriesLensedData, 1);

  ldata->cache_valid = FALSE;
  ldata->J           = g_new0 (gdouble, self->trunc_order + 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_series_lensed_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_series_lensed_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_series_lensed_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_series_lensed_ldata_required_columns;
}

/* Gaussian-envelope half-width margin (in units of sigma_noise) used to
 * window the rho integral around its peak at rho=R -- see _compute_J.
 * Same value and rationale as the noise-side scheme's own window
 * (independently chosen here, not shared code): 8 sigma leaves a tail of
 * exp(-32) ~ 1e-14, already below the truncation error of everything else
 * in this pipeline. */
#define NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA (8.0)

/* Comfortably positive (avoids log(0)=-Infinity, which risks Inf-Inf=NaN
 * arithmetic inside the caller's adaptive z-integral) but small enough
 * (-2*log(...) ~ 1381) to represent "this scheme's truncated approximation
 * has broken down here" for any practical purpose -- see the guard in
 * _nc_galaxy_shape_factor_series_lensed_marginal. */
#define NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_MIN_MARGINAL (1.0e-300)

/* Per-computation workspace, pooled via #NcmMemoryPool (see mp_ws's own
 * comment on the private struct above): bundles the NcmLaurentSeries arena
 * every function in the _F_g_taylor cascade draws its scratch objects from,
 * plus the rho_nodes/weights/GL-table triple _compute_J needs. n_nodes is
 * CONSTRUCT_ONLY (fixed for the instance's whole life), so rho_nodes/
 * weights/table are sized correctly here, once, and never resized; the
 * arena starts empty and is grown lazily by _arena_next() the first time
 * this workspace is actually used (see the arena's own comment above --
 * that growth is a one-time, first-node-only event, stable forever after).
 */
typedef struct _NcGalaxyShapeFactorSeriesLensedWS
{
  GPtrArray *arena;
  gdouble *rho_nodes;
  gdouble *weights;
  gsl_integration_glfixed_table *table;
} NcGalaxyShapeFactorSeriesLensedWS;

static gpointer
_nc_galaxy_shape_factor_series_lensed_ws_new (gpointer userdata)
{
  NcGalaxyShapeFactorSeriesLensed *gsfsl              = NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (userdata);
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);
  NcGalaxyShapeFactorSeriesLensedWS *ws               = g_new0 (NcGalaxyShapeFactorSeriesLensedWS, 1);

  ws->arena     = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_laurent_series_free);
  ws->rho_nodes = g_new (gdouble, self->n_nodes);
  ws->weights   = g_new (gdouble, self->n_nodes);
  ws->table     = gsl_integration_glfixed_table_alloc (self->n_nodes);

  return ws;
}

static void
_nc_galaxy_shape_factor_series_lensed_ws_free (gpointer p)
{
  NcGalaxyShapeFactorSeriesLensedWS *ws = (NcGalaxyShapeFactorSeriesLensedWS *) p;

  g_ptr_array_unref (ws->arena);
  g_free (ws->rho_nodes);
  g_free (ws->weights);
  gsl_integration_glfixed_table_free (ws->table);
  g_free (ws);
}

/*
 * Per-galaxy radial integral: fills J[0..order] with the (g-independent)
 * Bessel-reduced radial integrals at fixed (R,phi,sig2,sigma_pop), one fixed
 * Gauss-Legendre pass over rho in a window around the envelope's peak at
 * rho=R. Checks out one pooled workspace for the whole call (all n_nodes
 * iterations), returning it at the end -- same get/return granularity as
 * ncm_integral_locked_a_b's own bracketing of ncm_integral_get_workspace().
 */
static void
_compute_J (NcGalaxyShapeFactorSeriesLensedPrivate * const self, gdouble sigma_pop, gdouble pop_norm,
            gdouble R, gdouble phi, gdouble sig2, gdouble *J)
{
  const guint order                          = self->trunc_order;
  const guint n_nodes                        = self->n_nodes;
  const gdouble sigma_n                      = sqrt (sig2);
  const gdouble rho_lo                       = fmax (0.0, R - NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA * sigma_n);
  const gdouble rho_hi                       = fmin (1.0, R + NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA * sigma_n);
  NcGalaxyShapeFactorSeriesLensedWS **ws_ptr = ncm_memory_pool_get (self->mp_ws);
  NcGalaxyShapeFactorSeriesLensedWS *ws      = *ws_ptr;
  NcGalaxyShapeFactorSeriesLensedArena arena = { ws->arena, 0 };
  guint i, m;

  for (i = 0; i < n_nodes; i++)
    gsl_integration_glfixed_point (rho_lo, rho_hi, i, &ws->rho_nodes[i], &ws->weights[i], ws->table);

  for (m = 0; m <= order; m++)
    J[m] = 0.0;

  for (i = 0; i < n_nodes; i++)
  {
    const gdouble rho       = ws->rho_nodes[i];
    const gdouble w         = ws->weights[i];
    const gdouble z         = R * rho / sig2;
    const gdouble envelope  = exp (-gsl_pow_2 (R - rho) / (2.0 * sig2)) / (2.0 * M_PI * sig2);
    const gdouble prefactor = w * rho * envelope * pop_norm;
    NcmLaurentSeries **F    = g_new0 (NcmLaurentSeries *, order + 1);
    gint maxk, k;
    gdouble *Ik;

    arena.idx = 0;
    _F_g_taylor (&arena, self->chi_taylor, self->jac_taylor, rho, sigma_pop, order, F);

    maxk = 0;

    for (m = 0; m <= order; m++)
      maxk = MAX (maxk, ncm_laurent_series_hmax (F[m]));

    Ik = g_new (gdouble, maxk + 1);

    for (k = 0; k <= maxk; k++)
      Ik[k] = gsl_sf_bessel_In_scaled (k, z);

    for (m = 0; m <= order; m++)
      J[m] += prefactor * ncm_laurent_series_jacobi_anger_reduce (F[m], phi, Ik, maxk + 1);

    g_free (Ik);
    g_free (F);
  }

  ncm_memory_pool_return (ws_ptr);
}

static gdouble
_nc_galaxy_shape_factor_series_lensed_marginal (NcGalaxyShapeFactorSeriesLensed *gsfsl, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                const gdouble g_1, const gdouble g_2,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);
  NcGalaxyShapeFactorSeriesLensedData *ldata          = (NcGalaxyShapeFactorSeriesLensedData *) data->ldata;
  const guint64 pop_hash                              = nc_galaxy_shape_factor_get_pop_hash (NC_GALAXY_SHAPE_FACTOR (gsfsl));

  /* Gauge-fix: rotate (g,eps_obs) together by -arg(g) so g becomes real and
   * non-negative -- the chi_taylor/jac_taylor closed forms above are
   * real-g-only. Harmless when g=0. */
  const gdouble g_mag = hypot (g_1, g_2);
  const gdouble phi_g = (g_mag > 0.0) ? atan2 (g_2, g_1) : 0.0;
  const gdouble cphi  = cos (phi_g);
  const gdouble sphi  = sin (phi_g);
  const gdouble e1r   = cphi * epsilon_obs_1 + sphi * epsilon_obs_2;
  const gdouble e2r   = -sphi * epsilon_obs_1 + cphi * epsilon_obs_2;
  const gdouble R     = hypot (e1r, e2r);
  const gdouble phi   = atan2 (e2r, e1r);
  const gdouble sig2  = gsl_pow_2 (data->std_noise);
  gdouble result;
  gint m;

  if (!ldata->cache_valid || (ldata->cached_R != R) || (ldata->cached_phi != phi) ||
      (ldata->cached_std_noise != data->std_noise) || (ldata->cached_pop_hash != pop_hash))
  {
    const gdouble sigma_pop = nc_galaxy_shape_pop_get_sigma (pop, data->pop_data);

    /* NcGalaxyShapePopGauss (and its per-galaxy-sigma siblings) normalize
     * their density over the TRUNCATED disc (see that class' own "norm"
     * factor, nc_galaxy_shape_pop_gauss.c), not the plain
     * exp(-rho^2/2sigma_pop^2) this scheme's own Taylor-series composition
     * assumes -- read the actual normalization back from the population
     * object itself (eval_p_rho2 at rho^2=0 gives exactly that constant,
     * for any of the Gaussian-family siblings, all sharing the same
     * functional form per nc_galaxy_shape_pop_gauss.c's own docs) rather
     * than hardcoding the truncated-Gaussian formula here. The disc-polar
     * 1/pi convention matches the noise-side scheme's own
     * p_pop=eval_p(...)/M_PI convention (nc_galaxy_shape_pop_eval_p's own
     * docs). */
    const gdouble pop_norm = nc_galaxy_shape_pop_eval_p_rho2 (pop, data->pop_data, 0.0) / M_PI;

    _compute_J (self, sigma_pop, pop_norm, R, phi, sig2, ldata->J);
    ldata->cached_R         = R;
    ldata->cached_phi       = phi;
    ldata->cached_std_noise = data->std_noise;
    ldata->cached_pop_hash  = pop_hash;
    ldata->cache_valid      = TRUE;
  }

  result = ldata->J[self->trunc_order];

  for (m = (gint) self->trunc_order - 1; m >= 0; m--)
    result = result * g_mag + ldata->J[m];

  /* result is a degree-trunc_order polynomial in g_mag -- only a LOCAL
   * approximation around g=0. For configurations where the true P_pop*Jac
   * envelope decreases with g (the common case), a low-order polynomial fit
   * is *guaranteed* to eventually cross zero and diverge to -Infinity as
   * g_mag grows (an elementary property of nonconstant polynomials),
   * regardless of the true underlying function staying positive for any
   * physical |g|<1. Where exactly this crossing lands depends on both
   * trunc_order and (R,phi,std_noise) -- e.g. order=8 stays accurate well
   * past where order=4 has already broken down -- so this guards on the
   * sign of the actual computed value rather than a fixed g_mag threshold
   * (see docs/theory/wl_shape_factor_history.md for why a fixed threshold
   * does not work here). Left unguarded, a negative result reaches
   * eval_ln_marginal's log() as log(negative)=NaN, which poisons the
   * caller's adaptive z-integral and aborts the whole process once GSL's
   * quadrature fails to converge (ncm_integral1d_eval_lnint's own
   * g_error() path) -- refuse instead, matching the "vanishingly small but
   * always finite, always positive" contract this function's callers
   * expect. Checked as !(result>0.0), not result<=0.0, so it also catches
   * NaN (any comparison with NaN is false). */
  if (!(result > 0.0))
    return NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_MIN_MARGINAL;

  return result;
}

static gdouble
_nc_galaxy_shape_factor_series_lensed_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_series_lensed_marginal (NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

static gdouble
_nc_galaxy_shape_factor_series_lensed_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return log (_nc_galaxy_shape_factor_series_lensed_marginal (NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

static void
nc_galaxy_shape_factor_series_lensed_class_init (NcGalaxyShapeFactorSeriesLensedClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_series_lensed_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_series_lensed_get_property;
  object_class->constructed  = &_nc_galaxy_shape_factor_series_lensed_constructed;
  object_class->finalize     = &_nc_galaxy_shape_factor_series_lensed_finalize;

  /**
   * NcGalaxyShapeFactorSeriesLensed:trunc-order:
   *
   * Truncation order $N$ of the $g$-power series (both conventions'
   * closed forms generate any order equally easily). Default 4.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRUNC_ORDER,
                                   g_param_spec_uint ("trunc-order",
                                                      "Truncation order",
                                                      "Truncation order N of the g-power series",
                                                      1, G_MAXUINT, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorSeriesLensed:n-nodes:
   *
   * Number of fixed Gauss-Legendre nodes used for the per-galaxy radial
   * ($\rho$) integral. Default 60.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_NODES,
                                   g_param_spec_uint ("n-nodes",
                                                      "Number of radial nodes",
                                                      "Number of fixed Gauss-Legendre nodes for the rho integral",
                                                      2, G_MAXUINT, 60,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_series_lensed_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_series_lensed_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_series_lensed_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_series_lensed_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_series_lensed_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 * @trunc_order: truncation order $N$ of the $g$-power series, $N\ge1$
 *
 * Creates a new #NcGalaxyShapeFactorSeriesLensed.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorSeriesLensed.
 */
NcGalaxyShapeFactorSeriesLensed *
nc_galaxy_shape_factor_series_lensed_new (NcGalaxyWLObsEllipConv ellip_conv, guint trunc_order)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_SERIES_LENSED,
                       "ellip-conv", ellip_conv,
                       "trunc-order", trunc_order,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_series_lensed_ref:
 * @gsfsl: a #NcGalaxyShapeFactorSeriesLensed
 *
 * Increases the reference count of @gsfsl by one.
 *
 * Returns: (transfer full): @gsfsl.
 */
NcGalaxyShapeFactorSeriesLensed *
nc_galaxy_shape_factor_series_lensed_ref (NcGalaxyShapeFactorSeriesLensed *gsfsl)
{
  return g_object_ref (gsfsl);
}

/**
 * nc_galaxy_shape_factor_series_lensed_free:
 * @gsfsl: a #NcGalaxyShapeFactorSeriesLensed
 *
 * Decreases the reference count of @gsfsl by one.
 *
 */
void
nc_galaxy_shape_factor_series_lensed_free (NcGalaxyShapeFactorSeriesLensed *gsfsl)
{
  g_object_unref (gsfsl);
}

/**
 * nc_galaxy_shape_factor_series_lensed_clear:
 * @gsfsl: a #NcGalaxyShapeFactorSeriesLensed
 *
 * Decreases the reference count of *@gsfsl by one, and sets the pointer
 * *@gsfsl to NULL.
 *
 */
void
nc_galaxy_shape_factor_series_lensed_clear (NcGalaxyShapeFactorSeriesLensed **gsfsl)
{
  g_clear_object (gsfsl);
}

