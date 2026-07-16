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
 * ellipticity conventions, arbitrary truncation order, any population model
 * implementing #NcGalaxyShapePop's eval_p_rho2_g_series vfunc (currently
 * #NcGalaxyShapePopGauss, #NcGalaxyShapePopGaussLocal and
 * #NcGalaxyShapePopBeta).
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
 * This class is a thin orchestrator over two collaborators, each owning its
 * own closed-form derivation and doing its own construction-time
 * allocation:
 *
 * - #NcWLEllipticitySeriesTrace / #NcWLEllipticitySeriesTraceDet
 *   (nc_wl_ellipticity_series.h) supply $\chi_I(\chi_L,g)$'s Jacobian
 *   series and $|\chi_I(\chi_L,g)|^2$'s series (population-independent,
 *   pure shear-map output), one object per ellipticity convention, chosen
 *   once per workspace based on #NcGalaxyWLObsEllipConv.
 * - #NcGalaxyShapePop's eval_p_rho2_g_series vfunc composes
 *   $P_\mathrm{pop}(\chi_I(\chi_L,g))$'s own $g$-Taylor series from that
 *   $|\chi_I|^2$ series -- e.g. #NcGalaxyShapePopGauss composes
 *   $\exp(-|\chi_I|^2/2\sigma_\mathrm{pop}^2)$ via the standard
 *   exp-of-power-series recursion, with the $g^0$ term
 *   ($\chi_I(g=0)=\chi_L\neq0$) factored out as a scalar prefactor first.
 *
 * $F_g(\chi_L)$ itself is just $P_\mathrm{pop}(\chi_I(\chi_L,g))$'s series
 * convolved with the Jacobian series (ncm_laurent_series_tps_conv()).
 *
 * $\theta$-integral done EXACTLY via Jacobi-Anger
 * (ncm_laurent_series_jacobi_anger_reduce()); $\rho$-integral a fixed,
 * windowed Gauss-Legendre quadrature (see
 * #NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA). Caching is
 * two-layered, per galaxy: the Jacobi-Anger reduction's own
 * $\phi$-independent harmonic sums $H$ (the expensive part: Bessel
 * functions and Laurent-series algebra over the whole $\rho$ quadrature)
 * are keyed to $(R,\sigma_\mathrm{noise},\text{pop\_hash})$ -- deliberately
 * NOT to $\phi$, which is the argument of
 * $g=(1+m)g_\mathrm{true}+c_1+ic_2$ and so drifts with $g$ whenever the
 * galaxy's additive bias is nonzero. $J[0..\mathrm{order}]$, the actual
 * Horner-method polynomial coefficients, is instead recomputed from $H$ on
 * every call via a cheap trig-polynomial evaluation (_finalize_J()),
 * tracking $\phi$'s drift at negligible cost.
 *
 * A population that does not implement eval_p_rho2_g_series errors clearly
 * the first time this class actually needs it (base class default), rather
 * than at prepare() time.
 *
 * The $g$-Taylor series' own radius of convergence is set by eval_p's
 * nearest non-analytic point in $x=|\chi_I|^2$: entire for the Gaussian
 * family (#NcGalaxyShapePopGauss, #NcGalaxyShapePopGaussLocal), so any $g$
 * this project's own small-shear regime allows works; but
 * #NcGalaxyShapePopBeta's $P(x)\propto x^{\alpha-1}(1-x)^{\beta-1}$ has a
 * branch point at $x=0$ (if $\alpha<1$) and/or $x=1$ (if $\beta<1$), so the
 * series can diverge at even modest $g$ once the $\rho$ quadrature comes
 * close enough to that point -- confirmed by raising trunc_order making
 * such cases *worse*, the standard signature of evaluating a Taylor series
 * outside its own disk of convergence, rather than better as truncation
 * error alone would predict. $\alpha,\beta>1$ (an interior-peaked,
 * everywhere-smooth density) has no such restriction.
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
#include "nc/lss/wl/nc_wl_ellipticity_series.h"
#include "ncm/algebra/ncm_laurent_series.h"
#include "ncm/core/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

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
  NcGalaxyWLObsEllipConv ellip_conv;

  guint trunc_order;
  guint n_nodes;

  /* Pool of per-computation workspaces (NcGalaxyShapeFactorSeriesLensedWS,
   * defined further down, right before _compute_H), each bundling the
   * ellip_conv-appropriate #NcWLEllipticitySeriesTrace/TraceDet, this
   * class's own pop_coeffs/F, and the rho_nodes/weights/GL-table triple
   * _compute_H needs. NcmMemoryPool (ncm_memory_pool.h) is this project's
   * own thread-safe, auto-growing checkout pool -- same pattern as
   * ncm_integral_get_workspace() and NcmStatsDistKDE's mp_eval_vars:
   * _compute_H checks out one workspace per call and returns it when done,
   * so concurrent calls (were this ever parallelized) get independent
   * workspaces -- and independent, non-shared mutable series state --
   * rather than racing on shared state. */
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

/* Defined near _compute_H, forward-declared here for _init()'s use. */
static gpointer _nc_galaxy_shape_factor_series_lensed_ws_new (gpointer userdata);
static void _nc_galaxy_shape_factor_series_lensed_ws_free (gpointer p);

static void
nc_galaxy_shape_factor_series_lensed_init (NcGalaxyShapeFactorSeriesLensed *gsfsl)
{
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);

  self->trunc_order = 4;
  self->n_nodes     = 60;

  /* Lazily populated (ncm_memory_pool_new itself allocates no workspace),
   * so this runs safely before trunc-order/n-nodes's CONSTRUCT properties
   * are set -- _ws_new (defined near _compute_H) reads them from @gsfsl's
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

    self->ellip_conv = nc_galaxy_shape_factor_get_ellip_conv (gsf);

    switch (self->ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
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

/* Two independent per-galaxy caches, layered:
 *
 * - eps cache: epsilon_obs' own polar form (R,phi_eps). Depends only on
 *   epsilon_obs, which is fixed for the whole redshift-node loop.
 * - harmonics cache: H[m][k], the phi-INDEPENDENT Jacobi-Anger harmonic
 *   content of J[0..order] (see _compute_H). Depends on (R,std_noise,
 *   pop_hash) -- NOT on phi, since phi is the argument of
 *   (1+m)*g_true+(c1+i*c2) and so drifts with g whenever the galaxy's
 *   additive bias (c1,c2) is nonzero, while R/std_noise/pop_hash don't.
 *
 * J itself is recomputed from H on every call (cheap, see _finalize_J), so
 * it tracks phi's drift without needing to be part of either cache. */
typedef struct _NcGalaxyShapeFactorSeriesLensedData
{
  /* eps cache */
  gboolean eps_valid;
  gdouble cached_eps1, cached_eps2;
  gdouble R;
  gdouble phi_eps;

  /* harmonics cache */
  gboolean harmonics_valid;
  gdouble cached_std_noise;
  guint64 cached_pop_hash;
  guint maxk;        /* harmonic degree populated in H; a property of trunc_order
                      * and the ellipticity convention alone, never of rho */
  guint h_alloc;     /* allocated per-m stride of H (>= maxk+1); grows, never shrinks */
  complex double *H; /* flattened [m][k], size (trunc_order+1)*h_alloc */

  /* recomputed every call from H by _finalize_J */
  gdouble *J; /* size trunc_order+1: final per-phi polynomial coefficients */
} NcGalaxyShapeFactorSeriesLensedData;

static void
_nc_galaxy_shape_factor_series_lensed_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorSeriesLensedData *ldata = (NcGalaxyShapeFactorSeriesLensedData *) p;

  g_free (ldata->J);
  g_free (ldata->H);
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

  /* g_new0 already zero-initializes both caches to invalid/empty. */
  ldata->J = g_new0 (gdouble, self->trunc_order + 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_series_lensed_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_series_lensed_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_series_lensed_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_series_lensed_ldata_required_columns;
}

/* Gaussian-envelope half-width margin (in units of sigma_noise) used to
 * window the rho integral around its peak at rho=R -- see _compute_H.
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
 * comment on the private struct above). @series_trace/@series_trace_det:
 * only one is ever non-NULL, chosen at _ws_new() time from @self's own
 * @ellip_conv (resolved once, at constructed() time) and never revisited.
 * @pop_coeffs/@F are this class's own (population composed with the
 * Jacobian series); n_nodes is CONSTRUCT_ONLY (fixed for the instance's
 * whole life), so rho_nodes/weights/table are sized correctly here, once,
 * and never resized. Every one of these is owned, allocated exactly once
 * at workspace-creation time -- _compute_H's per-rho-node loop never
 * allocates anything, it only ever refills already-existing storage. */
typedef struct _NcGalaxyShapeFactorSeriesLensedWS
{
  NcWLEllipticitySeriesTrace *series_trace;
  NcWLEllipticitySeriesTraceDet *series_trace_det;
  NcmLaurentSeriesTPS *pop_coeffs;
  NcmLaurentSeriesTPS *F;
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

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      ws->series_trace = nc_wl_ellipticity_series_trace_new (self->trunc_order);
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      ws->series_trace_det = nc_wl_ellipticity_series_trace_det_new (self->trunc_order);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }

  ws->pop_coeffs = ncm_laurent_series_tps_new (self->trunc_order);
  ws->F          = ncm_laurent_series_tps_new (self->trunc_order);
  ws->rho_nodes  = g_new (gdouble, self->n_nodes);
  ws->weights    = g_new (gdouble, self->n_nodes);
  ws->table      = gsl_integration_glfixed_table_alloc (self->n_nodes);

  return ws;
}

static void
_nc_galaxy_shape_factor_series_lensed_ws_free (gpointer p)
{
  NcGalaxyShapeFactorSeriesLensedWS *ws = (NcGalaxyShapeFactorSeriesLensedWS *) p;

  nc_wl_ellipticity_series_trace_clear (&ws->series_trace);
  nc_wl_ellipticity_series_trace_det_clear (&ws->series_trace_det);
  ncm_laurent_series_tps_clear (&ws->pop_coeffs);
  ncm_laurent_series_tps_clear (&ws->F);
  g_free (ws->rho_nodes);
  g_free (ws->weights);
  gsl_integration_glfixed_table_free (ws->table);
  g_free (ws);
}

/*
 * Fills ldata->H[m][k] with the phi-independent Jacobi-Anger harmonic
 * content of J[0..order] at fixed (ldata->R,sig2,sigma_pop): a fixed
 * Gauss-Legendre pass over rho in a window around the envelope's peak at
 * rho=R. This is ncm_laurent_series_jacobi_anger_reduce()'s own sum
 *   Re(c_0)*Ik[0] + sum_{k=1}^{maxk} 2*Ik[k]*Re(c_k*exp(i*k*phi))
 * with the exp(i*k*phi) factor -- the only place phi enters -- pulled out
 * of the rho-node sum: H[m][k] accumulates prefactor*Ik[k]*c_k (k>=1) or
 * prefactor*Re(c_0)*Ik[0] (k=0); phi is applied afterward by _finalize_J.
 *
 * maxk (the harmonic degree F[m] actually populates) is a structural
 * property of trunc_order and the ellipticity convention alone -- every
 * harmonic slot the WL series objects and ws->F populate is set via an
 * integer loop index, never derived from rho's numeric value -- so it is
 * computed once, from the first rho node, and asserted stable for the
 * rest.
 *
 * Checks out one pooled workspace for the whole call (all n_nodes
 * iterations), returning it at the end -- same get/return granularity as
 * ncm_integral_locked_a_b's own bracketing of ncm_integral_get_workspace().
 */
static void
_compute_H (NcGalaxyShapeFactorSeriesLensedPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, gdouble pop_norm,
            gdouble sig2, NcGalaxyShapeFactorSeriesLensedData *ldata)
{
  const guint order                          = self->trunc_order;
  const guint n_nodes                        = self->n_nodes;
  const gdouble R                            = ldata->R;
  const gdouble sigma_n                      = sqrt (sig2);
  const gdouble rho_lo                       = fmax (0.0, R - NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA * sigma_n);
  const gdouble rho_hi                       = fmin (1.0, R + NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA * sigma_n);
  NcGalaxyShapeFactorSeriesLensedWS **ws_ptr = ncm_memory_pool_get (self->mp_ws);
  NcGalaxyShapeFactorSeriesLensedWS *ws      = *ws_ptr;
  guint i, m;
  guint h_stride = ldata->h_alloc;

  for (i = 0; i < n_nodes; i++)
    gsl_integration_glfixed_point (rho_lo, rho_hi, i, &ws->rho_nodes[i], &ws->weights[i], ws->table);

  for (i = 0; i < n_nodes; i++)
  {
    const gdouble rho       = ws->rho_nodes[i];
    const gdouble w         = ws->weights[i];
    const gdouble z         = R * rho / sig2;
    const gdouble envelope  = exp (-gsl_pow_2 (R - rho) / (2.0 * sig2)) / (2.0 * M_PI * sig2);
    const gdouble prefactor = w * rho * envelope * pop_norm;
    NcmLaurentSeriesTPS *x_series;
    NcmLaurentSeriesTPS *jac_series;
    gint maxk, k;
    gdouble *Ik;

    if (ws->series_trace != NULL)
    {
      nc_wl_ellipticity_series_trace_eval (ws->series_trace, rho);
      x_series   = nc_wl_ellipticity_series_trace_get_abs_sq (ws->series_trace);
      jac_series = nc_wl_ellipticity_series_trace_get_jac (ws->series_trace);
    }
    else
    {
      nc_wl_ellipticity_series_trace_det_eval (ws->series_trace_det, rho);
      x_series   = nc_wl_ellipticity_series_trace_det_get_abs_sq (ws->series_trace_det);
      jac_series = nc_wl_ellipticity_series_trace_det_get_jac (ws->series_trace_det);
    }

    nc_galaxy_shape_pop_eval_p_rho2_g_series (pop, pop_data, x_series, ws->pop_coeffs);
    ncm_laurent_series_tps_conv (ws->F, ws->pop_coeffs, jac_series);

    maxk = 0;

    for (m = 0; m <= order; m++)
      maxk = MAX (maxk, ncm_laurent_series_hmax (ncm_laurent_series_tps_get (ws->F, m)));

    if (i == 0)
    {
      const guint needed_stride = (guint) maxk + 1;

      if (ldata->h_alloc < needed_stride)
      {
        g_free (ldata->H);
        ldata->H       = g_new0 (complex double, (order + 1) * needed_stride);
        ldata->h_alloc = needed_stride;
      }
      else
      {
        guint idx;

        for (idx = 0; idx < (order + 1) * ldata->h_alloc; idx++)
          ldata->H[idx] = 0.0;
      }

      ldata->maxk = (guint) maxk;
      h_stride    = ldata->h_alloc;
    }
    else
    {
      /* See this function's own doc comment: maxk must be identical at
       * every rho node for a given trunc_order/convention. */
      g_assert_cmpint (maxk, ==, (gint) ldata->maxk);
    }

    Ik = g_new (gdouble, maxk + 1);

    for (k = 0; k <= maxk; k++)
      Ik[k] = gsl_sf_bessel_In_scaled (k, z);

    for (m = 0; m <= order; m++)
    {
      complex double *Hm   = &ldata->H[m * h_stride];
      NcmLaurentSeries *Fm = ncm_laurent_series_tps_get (ws->F, m);

      Hm[0] += prefactor * creal (ncm_laurent_series_get_c (Fm, 0)) * Ik[0];

      for (k = 1; k <= maxk; k++)
      {
        complex double v = ncm_laurent_series_get_c (Fm, k);

        if (v != 0.0)
          Hm[k] += prefactor * Ik[k] * v;
      }
    }

    g_free (Ik);
  }

  ncm_memory_pool_return (ws_ptr);
}

/* Cheap per-call finalization: applies the phi-dependent trig-polynomial
 * evaluation ncm_laurent_series_jacobi_anger_reduce() would have done
 * per-rho-node, but against the already-rho-summed harmonics in ldata->H:
 * J[m] = 2*pi*(Re(H[m][0]) + sum_k 2*Re(H[m][k]*exp(i*k*phi))), an
 * O(order*maxk) trig-polynomial evaluation instead of the full
 * Bessel/Laurent-series recomputation. */
static void
_finalize_J (NcGalaxyShapeFactorSeriesLensedData *ldata, guint order, gdouble phi)
{
  const guint h_stride = ldata->h_alloc;
  guint m;

  for (m = 0; m <= order; m++)
  {
    const complex double *Hm = &ldata->H[m * h_stride];
    gdouble term             = creal (Hm[0]);
    guint k;

    for (k = 1; k <= ldata->maxk; k++)
      term += 2.0 * creal (Hm[k] * cexp (I * (gdouble) k * phi));

    ldata->J[m] = 2.0 * M_PI * term;
  }
}

/* Refreshes the eps cache (R,phi_eps) when epsilon_obs itself changed. */
static void
_ensure_eps_cache (NcGalaxyShapeFactorSeriesLensedData *ldata, gdouble epsilon_obs_1, gdouble epsilon_obs_2)
{
  if (ldata->eps_valid && (ldata->cached_eps1 == epsilon_obs_1) && (ldata->cached_eps2 == epsilon_obs_2))
    return;

  ldata->R           = hypot (epsilon_obs_1, epsilon_obs_2);
  ldata->phi_eps     = atan2 (epsilon_obs_2, epsilon_obs_1);
  ldata->cached_eps1 = epsilon_obs_1;
  ldata->cached_eps2 = epsilon_obs_2;
  ldata->eps_valid   = TRUE;
}

/* Refreshes the harmonics cache (ldata->H) when (R,std_noise,pop_hash)
 * changed. pop_norm is a fixed 1/pi disc-measure factor, population-
 * agnostic: eval_p_rho2_g_series() composes the population's own fully
 * normalized density P(x(g)) directly (see its doc comment in
 * nc_galaxy_shape_pop.h), so no population-specific normalization lookup
 * is needed here. */
static void
_ensure_harmonics_cache (NcGalaxyShapeFactorSeriesLensedPrivate * const self, NcGalaxyShapeFactorSeriesLensedData *ldata,
                         NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, guint64 pop_hash, gdouble sig2)
{
  if (ldata->harmonics_valid && (ldata->cached_std_noise == data->std_noise) && (ldata->cached_pop_hash == pop_hash))
    return;

  {
    const gdouble pop_norm = 1.0 / M_PI;

    _compute_H (self, pop, data->pop_data, pop_norm, sig2, ldata);
  }

  ldata->cached_std_noise = data->std_noise;
  ldata->cached_pop_hash  = pop_hash;
  ldata->harmonics_valid  = TRUE;
}

/* Horner-evaluates the degree-trunc_order polynomial J[0..order] at
 * g_mag, clamping non-positive results to a tiny positive floor. This is
 * only a LOCAL approximation around g=0: whenever the true envelope
 * decreases with g, a low-order polynomial fit is guaranteed to eventually
 * cross zero and diverge to -Infinity as g_mag grows, regardless of the
 * true function staying positive -- caught here on the sign of the actual
 * result (not a fixed g_mag threshold, since where the crossing lands
 * depends on trunc_order/R/phi/std_noise; !(result>0.0) rather than
 * result<=0.0 so it also catches NaN). Left unclamped, log(negative)=NaN
 * would poison the caller's adaptive z-integral. */
static gdouble
_horner_eval_clamped (const gdouble *J, guint order, gdouble g_mag)
{
  gdouble result = J[order];
  gint m;

  for (m = (gint) order - 1; m >= 0; m--)
    result = result * g_mag + J[m];

  return (result > 0.0) ? result : NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_MIN_MARGINAL;
}

static gdouble
_nc_galaxy_shape_factor_series_lensed_marginal (NcGalaxyShapeFactorSeriesLensed *gsfsl, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                const gdouble g_1, const gdouble g_2,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSeriesLensedPrivate * const self = nc_galaxy_shape_factor_series_lensed_get_instance_private (gsfsl);
  NcGalaxyShapeFactorSeriesLensedData *ldata          = (NcGalaxyShapeFactorSeriesLensedData *) data->ldata;
  const guint64 pop_hash                              = nc_galaxy_shape_factor_get_pop_hash (NC_GALAXY_SHAPE_FACTOR (gsfsl));
  const gdouble sig2                                  = gsl_pow_2 (data->std_noise);

  /* Gauge-fix: (g,eps_obs) rotated together by -arg(g) so g becomes real
   * and non-negative -- the WL series objects are real-g-only. Rotation
   * preserves norm and is additive in angle, so the rotated eps_obs' polar
   * form is obtained directly from eps_obs' own (R,phi_eps): R is
   * unchanged by the rotation, phi=phi_eps-phi_g. */
  const gdouble g_mag = hypot (g_1, g_2);
  const gdouble phi_g = (g_mag > 0.0) ? atan2 (g_2, g_1) : 0.0;

  _ensure_eps_cache (ldata, epsilon_obs_1, epsilon_obs_2);
  _ensure_harmonics_cache (self, ldata, pop, data, pop_hash, sig2);
  _finalize_J (ldata, self->trunc_order, ldata->phi_eps - phi_g);

  return _horner_eval_clamped (ldata->J, self->trunc_order, g_mag);
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

  /* prepare: base default (no-op) is right too -- population-composition
   * capability is now checked by NcGalaxyShapePop's own
   * eval_p_rho2_g_series vfunc dispatch (see _compute_H), the first time
   * it is actually invoked, instead of a class-specific Gaussian-only
   * guard here. */
  gsf_class->data_init        = &_nc_galaxy_shape_factor_series_lensed_data_init;
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

