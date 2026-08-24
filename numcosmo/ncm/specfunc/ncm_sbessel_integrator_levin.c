/***************************************************************************
 *            ncm_sbessel_integrator_levin.c
 *
 *  Sat January 25 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_levin.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmSBesselIntegratorLevin:
 *
 * Levin-Bessel method for spherical Bessel function integration.
 *
 * This class implements integration of functions multiplied by spherical Bessel
 * functions using a Levin-type method for low multipoles and vector cubature
 * integration for high multipoles.
 *
 * The integral is written in the dimensionless variable $y = k x$, so that
 * $\int K(x, k) j_\ell(k x) \mathrm{d}x = \int F(y) j_\ell(y) \mathrm{d}y$ with
 * $F(y) = K(y / k, k) / k$. This is what allows one panel set to serve every $k$.
 *
 * For low ell values, the contribution of a panel $[a, b]$ is obtained by solving
 * $y^2 w''(y) + 2 y w'(y) + (y^2 - \ell(\ell+1)) w(y) = F(y)$ with boundary conditions
 * $w(a) = w(b) = 0$. Since $j_\ell$ solves the homogeneous equation, the combination
 * $y^2 (w' j_\ell - j_\ell' w)$ has derivative $F j_\ell$, so the panel integral is
 * the boundary term $b^2 w'(b) j_\ell(b) - a^2 w'(a) j_\ell(a)$.
 *
 * In practice the equation is solved for $u = y w$, which removes the first-derivative
 * term and yields $y^2 u'' + (y^2 - \ell(\ell+1)) u = y F(y)$ --- the forcing handed to
 * #NcmSBesselOdeSolver is the weighted $y F(y)$. Because $w$ vanishes at the endpoints,
 * $u'(a) = a w'(a)$ and $u'(b) = b w'(b)$ there, and the panel contribution actually
 * evaluated is $b j_\ell(b) u'(b) - a j_\ell(a) u'(a)$.
 *
 * For high ell values, it uses vector cubature integration where the integrand evaluates
 * $f(x)$ and all $j_\ell(x)$ values simultaneously for efficiency.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "ncm/specfunc/ncm_sbessel_ode_solver.h"
#include "ncm/sphere/ncm_spectral.h"
#include "ncm/specfunc/ncm_sf_sbessel.h"
#include "ncm/algebra/ncm_lapack.h"
#include "ncm/core/ncm_c.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "ncm/core/ncm_dtuple.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_integration.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * Per-panel diagnostic record. Off by default and zero cost when off.
 *
 * The integral is assembled as a sum of per-panel boundary terms that cancel
 * heavily, so the attainable relative accuracy is bounded by the cancellation
 * ratio sum|contrib| / |total| times the machine epsilon. That ratio cannot be
 * inferred from the result; recording the individual contributions is the only
 * way to measure it.
 */
typedef struct _NcmSBesselIntegratorLevinPanelRec
{
  gdouble a;
  gdouble b;
  gint ell;
  gdouble contrib;
} NcmSBesselIntegratorLevinPanelRec;

struct _NcmSBesselIntegratorLevin
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  guint max_order;
  gdouble reltol;
  guint cheb_min_order;
  gdouble cheb_reltol;
  NcmSBesselOdeSolver *ode_solver;
  NcmSBesselOdeOperator *ode_operator;
  NcmSFSBesselArray *sba; /* Allocation tracking */
  guint alloc_max_order;
  guint alloc_ell_min;
  guint alloc_ell_max;
  gboolean constructed;
  /* Pre-allocated working arrays */
  GArray *cheb_coeffs;
  GArray *edge_cheb_coeffs;
  gdouble *edge_transform_work;
  gsize edge_transform_work_len;
  GArray *gegen_coeffs;
  GArray *rhs;
  GArray *values_result;
  gdouble *j_array_a;
  gdouble *j_array_b;
  GArray *endpoints_result;
  gdouble *jl_arr;
  gboolean record_panels; /* Diagnostic panel recording (off by default) */
  GArray *panel_records;  /* NcmSBesselIntegratorLevinPanelRec, valid when recording */
  /* Knots-based paneling */
  gdouble y_knots_min;
  gdouble y_knots_max;
  guint n_knots;
  guint ell_cache_max;                        /* Maximum ell for precomputed j_l at knots */
  GArray *knots;                              /* Log-spaced knots array */
  gdouble *jl_knots;                          /* Precomputed j_l at knots: [n_knots * (ell_cache_max + 1)] */
  GPtrArray *operators;                       /* Operators for each panel between consecutive knots */
  GHashTable *edge_operators;                 /* Dyadic fixed-cell operators used by moving edge panels */
  gdouble panel_abstol;                       /* Absolute coefficient floor for the current panel's RHS (0 = relative only) */
  NcmSBesselOdeOperator *ode_operator_temp_a; /* Temporary operator for [a, smallest_knot > a] */
  NcmSBesselOdeOperator *ode_operator_temp_b; /* Temporary operator for [largest_knot < b, b] */
  gboolean ode_operator_temp_a_valid;         /* True when temp_a matches the cached panel */
  gboolean ode_operator_temp_b_valid;         /* True when temp_b matches the cached panel */
  gdouble ode_operator_temp_a_a;
  gdouble ode_operator_temp_a_b;
  gdouble ode_operator_temp_b_a;
  gdouble ode_operator_temp_b_b;
  guint ode_operator_temp_a_ell_min;
  guint ode_operator_temp_a_ell_max;
  guint ode_operator_temp_b_ell_min;
  guint ode_operator_temp_b_ell_max;
};

enum
{
  PROP_0,
  PROP_MAX_ORDER,
  PROP_RELTOL,
  PROP_CHEB_MIN_ORDER,
  PROP_CHEB_RELTOL,
  PROP_Y_KNOTS_MIN,
  PROP_Y_KNOTS_MAX,
  PROP_N_KNOTS,
  PROP_ELL_CACHE_MAX,
};

static void _ncm_sbessel_integrator_levin_prepare_knots_array (NcmSBesselIntegratorLevin *sbilv);
static void _ncm_sbessel_integrator_levin_prepare_ell_cache (NcmSBesselIntegratorLevin *sbilv);
static void _ncm_sbessel_integrator_levin_ensure_prepared (NcmSBesselIntegratorLevin *sbilv, guint max_order, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_prepare_knots_operators (NcmSBesselIntegratorLevin *sbilv, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_compute_rhs (NcmSBesselIntegratorLevin *sbilv, NcmSpectral *spectral, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, gpointer user_data);
static void _ncm_sbessel_integrator_levin_solve_and_accumulate (NcmSBesselIntegratorLevin *sbilv, NcmSpectral *spectral, NcmSBesselOdeOperator *operator, NcmSBesselIntegratorF F, gdouble a_p, gdouble b_p, const gdouble *j_a_p, const gdouble *j_b_p, gdouble k, guint ell_min, guint ell_max, gdouble *result_data, gpointer user_data);
static void _ncm_sbessel_integrator_levin_integrate_panel (NcmSBesselIntegratorLevin *sbilv, gint a_p_idx, gint b_p_idx, gdouble a_p, gdouble b_p, NcmSpectral *spectral, NcmSBesselIntegratorF F, gdouble k, guint ell_min, guint ell_max, gdouble *result_data, gpointer user_data);
static void _ncm_sbessel_integrator_levin_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data);

G_DEFINE_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_levin_init (NcmSBesselIntegratorLevin *sbilv)
{
  sbilv->max_order               = 0;
  sbilv->reltol                  = 0.0;
  sbilv->cheb_min_order          = 0;
  sbilv->cheb_reltol             = 0.0;
  sbilv->ode_solver              = ncm_sbessel_ode_solver_new ();
  sbilv->ode_operator            = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, 2, 2);
  sbilv->sba                     = ncm_sf_sbessel_array_new ();
  sbilv->alloc_max_order         = 0;
  sbilv->panel_abstol            = 0.0;
  sbilv->alloc_ell_min           = -1;
  sbilv->alloc_ell_max           = -1;
  sbilv->cheb_coeffs             = NULL;
  sbilv->edge_cheb_coeffs        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->edge_transform_work     = NULL;
  sbilv->edge_transform_work_len = 0;
  sbilv->gegen_coeffs            = NULL;
  sbilv->rhs                     = NULL;
  sbilv->values_result           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->j_array_a               = NULL;
  sbilv->j_array_b               = NULL;
  sbilv->endpoints_result        = NULL;
  sbilv->jl_arr                  = NULL;
  sbilv->constructed             = FALSE;
  sbilv->record_panels           = FALSE;
  sbilv->panel_records           = g_array_new (FALSE, FALSE, sizeof (NcmSBesselIntegratorLevinPanelRec));
  /* Knots-based paneling */
  sbilv->y_knots_min    = 0.0;
  sbilv->y_knots_max    = 0.0;
  sbilv->n_knots        = 0;
  sbilv->ell_cache_max  = 0;
  sbilv->knots          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->jl_knots       = NULL;
  sbilv->operators      = NULL;
  sbilv->edge_operators = g_hash_table_new_full (g_int64_hash, g_int64_equal, g_free,
                                                 (GDestroyNotify) ncm_sbessel_ode_operator_unref);
  sbilv->ode_operator_temp_a       = NULL;
  sbilv->ode_operator_temp_b       = NULL;
  sbilv->ode_operator_temp_a_valid = FALSE;
  sbilv->ode_operator_temp_b_valid = FALSE;
}

static void
_ncm_sbessel_integrator_levin_dispose (GObject *object)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  ncm_sbessel_ode_solver_clear (&sbilv->ode_solver);
  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator);
  ncm_sf_sbessel_array_clear (&sbilv->sba);
  g_clear_pointer (&sbilv->panel_records, g_array_unref);
  g_clear_pointer (&sbilv->cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->edge_cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->edge_transform_work, g_free);
  g_clear_pointer (&sbilv->gegen_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->rhs, g_array_unref);
  g_clear_pointer (&sbilv->endpoints_result, g_array_unref);

  /* Clear knots-based paneling resources */
  g_clear_pointer (&sbilv->knots, g_array_unref);

  if (sbilv->jl_knots != NULL)
  {
    g_free (sbilv->jl_knots);
    sbilv->jl_knots = NULL;
  }

  if (sbilv->operators != NULL)
  {
    g_ptr_array_unref (sbilv->operators);
    sbilv->operators = NULL;
  }

  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator_temp_a);
  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator_temp_b);
  g_clear_pointer (&sbilv->edge_operators, g_hash_table_unref);

  g_clear_pointer (&sbilv->values_result, g_array_unref);

  if (sbilv->j_array_a != NULL)
  {
    g_free (sbilv->j_array_a);
    sbilv->j_array_a = NULL;
  }

  if (sbilv->j_array_b != NULL)
  {
    g_free (sbilv->j_array_b);
    sbilv->j_array_b = NULL;
  }

  if (sbilv->jl_arr != NULL)
  {
    g_free (sbilv->jl_arr);
    sbilv->jl_arr = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_levin_parent_class)->dispose (object);
}

static void
_ncm_sbessel_integrator_levin_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_levin_parent_class)->finalize (object);
}

static void
_ncm_sbessel_integrator_levin_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_sbessel_integrator_levin_parent_class)->constructed (object);

  {
    NcmSBesselIntegrator *sbi        = NCM_SBESSEL_INTEGRATOR (object);
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);
    guint ell_min, ell_max;

    /* Prepare knots array and precompute spherical Bessel functions at construction time
     * since knots parameters and ell_cache_max are CONSTRUCT_ONLY */
    _ncm_sbessel_integrator_levin_prepare_knots_array (sbilv);
    _ncm_sbessel_integrator_levin_prepare_ell_cache (sbilv);

    /* Mark as constructed and trigger set_ell_range to initialize operators */
    sbilv->constructed = TRUE;
    ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
    ncm_sbessel_integrator_set_ell_range (sbi, ell_min, ell_max);
  }
}

static void
_ncm_sbessel_integrator_levin_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_LEVIN (object));

  switch (prop_id)
  {
    case PROP_MAX_ORDER:
      ncm_sbessel_integrator_levin_set_max_order (sbilv, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      ncm_sbessel_integrator_levin_set_reltol (sbilv, g_value_get_double (value));
      break;
    case PROP_CHEB_MIN_ORDER:
      ncm_sbessel_integrator_levin_set_cheb_min_order (sbilv, g_value_get_uint (value));
      break;
    case PROP_CHEB_RELTOL:
      ncm_sbessel_integrator_levin_set_cheb_reltol (sbilv, g_value_get_double (value));
      break;
    case PROP_Y_KNOTS_MIN:
      sbilv->y_knots_min = g_value_get_double (value);
      break;
    case PROP_Y_KNOTS_MAX:
      sbilv->y_knots_max = g_value_get_double (value);
      break;
    case PROP_N_KNOTS:
      sbilv->n_knots = g_value_get_uint (value);
      break;
    case PROP_ELL_CACHE_MAX:
      sbilv->ell_cache_max = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_integrator_levin_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_LEVIN (object));

  switch (prop_id)
  {
    case PROP_MAX_ORDER:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_max_order (sbilv));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_reltol (sbilv));
      break;
    case PROP_CHEB_MIN_ORDER:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_cheb_min_order (sbilv));
      break;
    case PROP_CHEB_RELTOL:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_cheb_reltol (sbilv));
      break;
    case PROP_Y_KNOTS_MIN:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_y_knots_min (sbilv));
      break;
    case PROP_Y_KNOTS_MAX:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_y_knots_max (sbilv));
      break;
    case PROP_N_KNOTS:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_n_knots (sbilv));
      break;
    case PROP_ELL_CACHE_MAX:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_ell_cache_max (sbilv));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_sbessel_integrator_levin_class_init (NcmSBesselIntegratorLevinClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmSBesselIntegratorClass *parent_class = NCM_SBESSEL_INTEGRATOR_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_integrator_levin_set_property;
  object_class->get_property = &_ncm_sbessel_integrator_levin_get_property;
  object_class->constructed  = &_ncm_sbessel_integrator_levin_constructed;
  object_class->dispose      = &_ncm_sbessel_integrator_levin_dispose;
  object_class->finalize     = &_ncm_sbessel_integrator_levin_finalize;

  /**
   * NcmSBesselIntegratorLevin:max-order:
   *
   * Maximum order of Chebyshev decomposition for the Levin method. Higher order may
   * give better accuracy but is more expensive. Default is 2^14.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum Chebyshev order",
                                                      2, G_MAXUINT, 1 << 14,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:reltol:
   *
   * Relative tolerance for convergence.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:cheb-min-order:
   *
   * Minimum order of Chebyshev decomposition used when computing the RHS for the
   * Levin method.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_MIN_ORDER,
                                   g_param_spec_uint ("cheb-min-order",
                                                      NULL,
                                                      "Minimum Chebyshev order for RHS",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:cheb-reltol:
   *
   * Relative tolerance for Chebyshev decomposition of the integrand when computing
   * the RHS for the Levin method.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_RELTOL,
                                   g_param_spec_double ("cheb-reltol",
                                                        NULL,
                                                        "Chebyshev decomposition relative tolerance",
                                                        0.0, 1.0, 1.0e-8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:y-knots-min:
   *
   * Minimum value for knots in log-spaced grid. Set to 0 to disable knots-based
   * paneling. This property can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_Y_KNOTS_MIN,
                                   g_param_spec_double ("y-knots-min",
                                                        NULL,
                                                        "Minimum knot value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:y-knots-max:
   *
   * Maximum value for knots in log-spaced grid. Set to 0 to disable knots-based
   * paneling. This property can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_Y_KNOTS_MAX,
                                   g_param_spec_double ("y-knots-max",
                                                        NULL,
                                                        "Maximum knot value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:n-knots:
   *
   * Number of knots in the log-spaced grid. The knots will be equally spaced in log
   * space between y-knots-min and y-knots-max. Set to 0 to disable knots-based
   * paneling. This property can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_KNOTS,
                                   g_param_spec_uint ("n-knots",
                                                      NULL,
                                                      "Number of knots",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:ell-cache-max:
   *
   * Maximum ell value for precomputed spherical Bessel functions at knots. The
   * integrator will precompute j_ell(knot) for all knots and all ell from 0 to
   * ell-cache-max. This enables fast lookup during integration when the requested
   * ell values are within the cached range. For ell values beyond ell-cache-max,
   * the integrator will compute spherical Bessel functions on-the-fly. This property
   * can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_CACHE_MAX,
                                   g_param_spec_uint ("ell-cache-max",
                                                      NULL,
                                                      "Maximum ell for cache",
                                                      0, G_MAXUINT, 500,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->set_ell_range = &_ncm_sbessel_integrator_levin_set_ell_range;
  parent_class->integrate     = &_ncm_sbessel_integrator_levin_integrate;
}

static void
_ncm_sbessel_integrator_levin_ensure_prepared (NcmSBesselIntegratorLevin *sbilv, guint max_order, guint ell_min, guint ell_max)
{
  if (sbilv->alloc_max_order == max_order)
    return;

  /* Free existing allocations */
  g_clear_pointer (&sbilv->cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->gegen_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->rhs, g_array_unref);
  g_clear_pointer (&sbilv->endpoints_result, g_array_unref);

  /* Allocate arrays for spectral coefficients */
  sbilv->cheb_coeffs  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order);
  sbilv->gegen_coeffs = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order);
  sbilv->rhs          = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order + 2);

  /* Allocate result matrix for batched endpoint computation (max block size is 8) */
  sbilv->endpoints_result = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 8 * 3);

  sbilv->alloc_max_order = max_order;
}

/**
 * _ncm_sbessel_integrator_levin_prepare_knots_array:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Creates the log-spaced knots array from y_knots_min to y_knots_max.
 * This only needs to be done when knot parameters change, not when ell range changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_knots_array (NcmSBesselIntegratorLevin *sbilv)
{
  if (sbilv->n_knots < 2)
    return;

  g_assert_cmpfloat (sbilv->y_knots_min, >, 0.0);
  g_assert_cmpfloat (sbilv->y_knots_max, >, sbilv->y_knots_min);
  g_assert_cmpuint (sbilv->n_knots, <, G_MAXUINT / 2); /* Prevent overflow in log spacing calculation */
  g_array_set_size (sbilv->knots, sbilv->n_knots);

  {
    const gdouble ln_y_min = log (sbilv->y_knots_min);
    const gdouble ln_y_max = log (sbilv->y_knots_max);
    const gdouble L        = ln_y_max - ln_y_min;
    const gdouble dL       = L / (sbilv->n_knots - 1.0);
    const gdouble expm1_dL = expm1 (dL);
    gdouble y0             = exp (ln_y_min);
    guint i;

    g_array_index (sbilv->knots, gdouble, 0) = y0;

    for (i = 1; i < sbilv->n_knots; i++)
    {
      const gdouble dy = y0 * expm1_dL;
      const gdouble y  = y0 + dy;

      g_array_index (sbilv->knots, gdouble, i) = y;
      y0                                       = y;
    }
  }
}

/**
 * _ncm_sbessel_integrator_levin_prepare_ell_cache:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * Precomputes spherical Bessel functions j_ell at all knots for ell from 0 to ell_cache_max.
 * This is only needed when ell_cache_max or knot parameters change, not when ell range changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_ell_cache (NcmSBesselIntegratorLevin *sbilv)
{
  const guint n_ell = sbilv->ell_cache_max + 1;
  guint i;

  g_assert_cmpuint (sbilv->ell_cache_max, <, G_MAXUINT / 2);

  /* Allocate flat array: jl_knots[knot_idx * n_ell + ell] */
  sbilv->jl_knots = g_new (gdouble, sbilv->n_knots * n_ell);

  for (i = 0; i < sbilv->n_knots; i++)
  {
    const gdouble y     = g_array_index (sbilv->knots, gdouble, i);
    gdouble *jl_at_knot = &sbilv->jl_knots[i * n_ell];

    /* Compute all j_ell(y) for ell = 0, 1, ..., ell_cache_max */
    ncm_sf_sbessel_array_eval (sbilv->sba, sbilv->ell_cache_max, y, jl_at_knot);
  }

  /* Allocate arrays for spherical Bessel functions */
  sbilv->j_array_a = g_new0 (gdouble, n_ell);
  sbilv->j_array_b = g_new0 (gdouble, n_ell);
  sbilv->jl_arr    = g_new (gdouble, n_ell);
}

/**
 * _ncm_sbessel_integrator_levin_prepare_knots_operators:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 *
 * Prepares ODE operators for the knots-based paneling system:
 * - Pre-allocated ODE operators for each panel between consecutive knots
 * - Two temporary operators for edge panels [a, smallest_knot > a] and [largest_knot < b, b]
 *
 * Uses ncm_sbessel_ode_operator_reset() to efficiently update operators when ell range changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_knots_operators (NcmSBesselIntegratorLevin *sbilv, guint ell_min, guint ell_max)
{
  if (sbilv->n_knots < 2)
  {
    return;
  }
  else
  {
    const gboolean need_create = (sbilv->operators == NULL);
    const gboolean need_reset  = !need_create && ((sbilv->alloc_ell_min != ell_min) || (sbilv->alloc_ell_max != ell_max));
    guint i;

    if (!need_create && !need_reset)
      return;

    if (need_create)
    {
      sbilv->operators = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_sbessel_ode_operator_unref);

      for (i = 0; i < sbilv->n_knots - 1; i++)
      {
        const gdouble y_a = g_array_index (sbilv->knots, gdouble, i);
        const gdouble y_b = g_array_index (sbilv->knots, gdouble, i + 1);
        NcmSBesselOdeOperator *op;

        op = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, y_a, y_b, ell_min, ell_max);
        g_ptr_array_add (sbilv->operators, op);
      }

      sbilv->ode_operator_temp_a       = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_b       = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_a_valid = FALSE;
      sbilv->ode_operator_temp_b_valid = FALSE;
    }
    else if (need_reset)
    {
      g_hash_table_remove_all (sbilv->edge_operators);

      for (i = 0; i < sbilv->operators->len; i++)
      {
        NcmSBesselOdeOperator *op = g_ptr_array_index (sbilv->operators, i);
        const gdouble y_a         = g_array_index (sbilv->knots, gdouble, i);
        const gdouble y_b         = g_array_index (sbilv->knots, gdouble, i + 1);

        ncm_sbessel_ode_operator_reset (op, y_a, y_b, ell_min, ell_max);
      }

      /* Reset temporary operators */
      ncm_sbessel_ode_operator_reset (sbilv->ode_operator_temp_a, 0.0, 1.0, ell_min, ell_max);
      ncm_sbessel_ode_operator_reset (sbilv->ode_operator_temp_b, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_a_valid = FALSE;
      sbilv->ode_operator_temp_b_valid = FALSE;
    }

    sbilv->alloc_ell_min = ell_min;
    sbilv->alloc_ell_max = ell_max;
  }
}

/* Wrapper function that transforms K(x, k) -> f(y) = K(y/k, k)/k */
typedef struct _NcmSBesselIntegratorLevinWrapper
{
  NcmSBesselIntegratorF K;
  gdouble k;
  gpointer user_data;
} NcmSBesselIntegratorLevinWrapper;

static gdouble
_ncm_sbessel_integrator_levin_wrapper_func (gpointer data, gdouble y)
{
  NcmSBesselIntegratorLevinWrapper *wrapper = (NcmSBesselIntegratorLevinWrapper *) data;

  const gdouble x     = y / wrapper->k;
  const gdouble K_val = wrapper->K (wrapper->user_data, x, wrapper->k);

  return x * K_val;
}

static void
_ncm_sbessel_integrator_levin_build_rhs (NcmSBesselIntegratorLevin *sbilv)
{
  ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);
  g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
  {
    gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
    const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

    rhs_data[0] = 0.0;
    rhs_data[1] = 0.0;
    memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
  }
}

/**
 * _ncm_sbessel_integrator_levin_compute_rhs:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @spectral: spectral methods object
 * @F: integrand function K(x, k)
 * @a: lower integration bound in y-space
 * @b: upper integration bound in y-space
 * @k: wave number parameter
 * @user_data: user data for integrand
 *
 * Computes the RHS for the Levin ODE by:
 * 1. Computing Chebyshev coefficients for f(y) = K(y/k, k)/k
 * 2. Converting to Gegenbauer C^(2) basis
 * 3. Setting up RHS with homogeneous boundary conditions
 */
static void
_ncm_sbessel_integrator_levin_compute_rhs (NcmSBesselIntegratorLevin *sbilv,
                                           NcmSpectral *spectral,
                                           NcmSBesselIntegratorF F,
                                           gdouble a, gdouble b,
                                           gdouble k,
                                           gpointer user_data)
{
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};

  ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, &_ncm_sbessel_integrator_levin_wrapper_func,
                                                       a, b, sbilv->cheb_min_order, sbilv->cheb_reltol, sbilv->panel_abstol,
                                                       &sbilv->cheb_coeffs, &wrapper);

  _ncm_sbessel_integrator_levin_build_rhs (sbilv);
}

/**
 * _ncm_sbessel_integrator_levin_get_panel_resources:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @a_p_idx: knot index for left endpoint (-1 if not on a knot)
 * @b_p_idx: knot index for right endpoint (-1 if not on a knot)
 * @a_p: lower bound of panel
 * @b_p: upper bound of panel
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @j_a_p_out: (out): pointer to j_ell array at a_p
 * @j_b_p_out: (out): pointer to j_ell array at b_p
 *
 * Gets panel resources: j_ell arrays and appropriate ODE operator.
 * If a_p_idx == -1, computes j_a_p and configures temp_a if the panel changed.
 * If b_p_idx == -1, computes j_b_p and configures temp_b if the panel changed.
 * Otherwise, uses cached j_ell values and operators[a_p_idx] without reset.
 *
 * Returns: (transfer none): the appropriate ODE operator for this panel
 */
static NcmSBesselOdeOperator *
_ncm_sbessel_integrator_levin_get_panel_resources (NcmSBesselIntegratorLevin *sbilv,
                                                   gint a_p_idx, gint b_p_idx,
                                                   gdouble a_p, gdouble b_p,
                                                   guint ell_min, guint ell_max,
                                                   const gdouble **j_a_p_out,
                                                   const gdouble **j_b_p_out)
{
  const guint n_ell = sbilv->ell_cache_max + 1;
  NcmSBesselOdeOperator *op;

  /* Get j_a_p: compute at a moving endpoint or use a knot cache. */
  if (a_p_idx < 0)
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, a_p, sbilv->j_array_a);
    *j_a_p_out = sbilv->j_array_a;
  }
  else
  {
    *j_a_p_out = &sbilv->jl_knots[a_p_idx * n_ell];
  }

  /* Get j_b_p: compute at a moving endpoint or use a knot cache. */
  if (b_p_idx < 0)
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, b_p, sbilv->j_array_b);
    *j_b_p_out = sbilv->j_array_b;
  }
  else
  {
    *j_b_p_out = &sbilv->jl_knots[b_p_idx * n_ell];
  }

  /* Reuse an edge factorization until its interval or ell range changes. */
  if (a_p_idx < 0)
  {
    op = sbilv->ode_operator_temp_a;

    if (!sbilv->ode_operator_temp_a_valid ||
        (sbilv->ode_operator_temp_a_a != a_p) ||
        (sbilv->ode_operator_temp_a_b != b_p) ||
        (sbilv->ode_operator_temp_a_ell_min != ell_min) ||
        (sbilv->ode_operator_temp_a_ell_max != ell_max))
    {
      ncm_sbessel_ode_operator_reset (op, a_p, b_p, ell_min, ell_max);
      sbilv->ode_operator_temp_a_a       = a_p;
      sbilv->ode_operator_temp_a_b       = b_p;
      sbilv->ode_operator_temp_a_ell_min = ell_min;
      sbilv->ode_operator_temp_a_ell_max = ell_max;
      sbilv->ode_operator_temp_a_valid   = TRUE;
    }
  }
  else if (b_p_idx < 0)
  {
    op = sbilv->ode_operator_temp_b;

    if (!sbilv->ode_operator_temp_b_valid ||
        (sbilv->ode_operator_temp_b_a != a_p) ||
        (sbilv->ode_operator_temp_b_b != b_p) ||
        (sbilv->ode_operator_temp_b_ell_min != ell_min) ||
        (sbilv->ode_operator_temp_b_ell_max != ell_max))
    {
      ncm_sbessel_ode_operator_reset (op, a_p, b_p, ell_min, ell_max);
      sbilv->ode_operator_temp_b_a       = a_p;
      sbilv->ode_operator_temp_b_b       = b_p;
      sbilv->ode_operator_temp_b_ell_min = ell_min;
      sbilv->ode_operator_temp_b_ell_max = ell_max;
      sbilv->ode_operator_temp_b_valid   = TRUE;
    }
  }
  else
  {
    op = g_ptr_array_index (sbilv->operators, a_p_idx);
  }

  return op;
}

/* Relative level below which a panel cannot move the accumulated result. */
#define NCM_SBESSEL_LEVIN_PANEL_ABSTOL_EPS 1.0e-16

/*
 * Absolute floor for the next panel's Chebyshev RHS. Two independent scales
 * feed it, and the looser one wins:
 *
 *  - what the caller declared via ncm_sbessel_integrator_set_abstol(): an
 *    absolute error it can tolerate in this integral, because it knows the
 *    larger quantity the integral feeds into. This is the only scale available
 *    when the integral itself is numerically zero;
 *  - what the result already holds: a panel that cannot move any ell's
 *    accumulated value needs no further refinement.
 *
 * One RHS is shared by every ell in the batch, so the smallest accumulated
 * |result| over the batch sets the second scale. The panel contribution is
 * b_p * j_ell(b_p) * u'(b_p) - a_p * j_ell(a_p) * u'(a_p), so
 * max (b_p * |j_ell(b_p)|, a_p * |j_ell(a_p)|), maximized over the batch,
 * converts a bound on the result into a bound on the coefficients. Bounding
 * that scale by b_p alone (|j_ell| <= 1) costs many orders of magnitude where
 * the panel sits below the ell-th Bessel turning point: there j_ell is
 * evanescent, the panel cannot move the result whatever the RHS looks like,
 * and demanding a relative fit of the integrand there is both futile and
 * expensive. The scale used here is never larger than b_p, so the floor is
 * never tighter than the |j_ell| <= 1 one.
 *
 * Without either scale this returns 0.0 -- the pure relative criterion.
 */
static gdouble
_ncm_sbessel_integrator_levin_panel_abstol (NcmSBesselIntegratorLevin *sbilv,
                                            const gdouble *result_data,
                                            const gdouble *j_a_p, const gdouble *j_b_p,
                                            gdouble a_p, gdouble b_p,
                                            guint ell_min, guint ell_max)
{
  const gdouble caller_abstol = ncm_sbessel_integrator_get_abstol (NCM_SBESSEL_INTEGRATOR (sbilv));
  gdouble min_abs             = G_MAXDOUBLE;
  gdouble max_scale           = 0.0;
  gdouble abstol_result;
  guint ell;

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    const gdouble abs_ell   = fabs (result_data[ell - ell_min]);
    const gdouble scale_ell = MAX (b_p * fabs (j_b_p[ell]), a_p * fabs (j_a_p[ell]));

    min_abs   = MIN (min_abs, abs_ell);
    max_scale = MAX (max_scale, scale_ell);
  }

  abstol_result = MAX (NCM_SBESSEL_LEVIN_PANEL_ABSTOL_EPS * min_abs, caller_abstol);

  /* j_ell underflowed to zero for the whole batch: the panel contributes
   * exactly nothing, so any RHS will do. */
  if (max_scale == 0.0)
    return G_MAXDOUBLE;

  return abstol_result / max_scale;
}

/**
 * _ncm_sbessel_integrator_levin_solve_and_accumulate:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @spectral: spectral methods object
 * @operator: ODE operator for the panel
 * @F: integrand function
 * @a_p: lower bound of panel
 * @b_p: upper bound of panel
 * @j_a_p: array of j_ell(a_p) values
 * @j_b_p: array of j_ell(b_p) values
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @result_data: array to accumulate results
 * @user_data: user data for integrand
 *
 * Solves the Levin ODE for a panel [a_p, b_p] and accumulates (or assigns)
 * the contribution to the result array.
 */
static void
_ncm_sbessel_integrator_levin_solve_and_accumulate (NcmSBesselIntegratorLevin *sbilv,
                                                    NcmSpectral *spectral,
                                                    NcmSBesselOdeOperator *operator,
                                                    NcmSBesselIntegratorF F,
                                                    gdouble a_p, gdouble b_p,
                                                    const gdouble *j_a_p,
                                                    const gdouble *j_b_p,
                                                    gdouble k,
                                                    guint ell_min, guint ell_max,
                                                    gdouble *result_data,
                                                    gpointer user_data)
{
  guint ell;

  sbilv->panel_abstol = _ncm_sbessel_integrator_levin_panel_abstol (sbilv, result_data,
                                                                    j_a_p, j_b_p, a_p, b_p,
                                                                    ell_min, ell_max);

  _ncm_sbessel_integrator_levin_compute_rhs (sbilv, spectral, F, a_p, b_p, k, user_data);
  ncm_sbessel_ode_operator_solve_endpoints (operator, sbilv->rhs, &sbilv->endpoints_result);

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    const gint ell_idx      = ell - ell_min;
    const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, ell_idx * 3 + 0);
    const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, ell_idx * 3 + 1);
    const gdouble j_l_a     = j_a_p[ell];
    const gdouble j_l_b     = j_b_p[ell];
    const gdouble contrib   = b_p * j_l_b * y_prime_b - a_p * j_l_a * y_prime_a;

    result_data[ell_idx] += contrib;

    if (G_UNLIKELY (sbilv->record_panels))
    {
      const NcmSBesselIntegratorLevinPanelRec rec = {a_p, b_p, (gint) ell, contrib};

      g_array_append_val (sbilv->panel_records, rec);
    }
  }
}

static gdouble
_ncm_sbessel_integrator_levin_yj_deriv (guint ell, gdouble y, const gdouble *j)
{
  if (ell > 0)
    return y * j[ell - 1] - ell * j[ell];

  return cos (y);
}

static NcmSBesselOdeOperator *
_ncm_sbessel_integrator_levin_get_edge_operator (NcmSBesselIntegratorLevin *sbilv,
                                                 guint panel_idx, gboolean right_edge,
                                                 gdouble integral_a, gdouble integral_b,
                                                 guint ell_min, guint ell_max,
                                                 gdouble *panel_a, gdouble *panel_b)
{
  const gdouble coarse_a = g_array_index (sbilv->knots, gdouble, panel_idx);
  const gdouble coarse_b = g_array_index (sbilv->knots, gdouble, panel_idx + 1);
  const gdouble span     = integral_b - integral_a;
  gdouble width          = coarse_b - coarse_a;
  guint level            = 0;

  g_assert_cmpfloat (span, >, 0.0);

  /* Select the smallest dyadic cell containing the edge.  Consequently its
   * width is less than twice the requested span, avoiding the very high
   * solution orders of a complete coarse logarithmic panel. */
  while ((0.5 * width >= span) && (level < 52))
  {
    width *= 0.5;
    level++;
  }

  if (right_edge)
  {
    *panel_a = coarse_a;
    *panel_b = coarse_a + width;
  }
  else
  {
    *panel_a = coarse_b - width;
    *panel_b = coarse_b;
  }

  if (level == 0)
  {
    return g_ptr_array_index (sbilv->operators, panel_idx);
  }
  else
  {
    const guint64 key_value   = ((guint64) right_edge << 63) | ((guint64) panel_idx << 16) | level;
    NcmSBesselOdeOperator *op = g_hash_table_lookup (sbilv->edge_operators, &key_value);

    if (op == NULL)
    {
      guint64 *key = g_new (guint64, 1);

      *key = key_value;
      op   = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, *panel_a, *panel_b, ell_min, ell_max);
      g_hash_table_insert (sbilv->edge_operators, key, op);
    }

    return op;
  }
}

typedef struct _NcmSBesselIntegratorLevinEdgeExtension
{
  GArray *coeffs;
  gdouble a;
  gdouble b;
} NcmSBesselIntegratorLevinEdgeExtension;

static gdouble
_ncm_sbessel_integrator_levin_edge_extension_func (gpointer data, gdouble y)
{
  NcmSBesselIntegratorLevinEdgeExtension *extension = data;

  return ncm_spectral_chebyshev_eval_x (extension->coeffs, extension->a, extension->b, y);
}

static gboolean
_ncm_sbessel_integrator_levin_transform_edge_coeffs (NcmSBesselIntegratorLevin *sbilv,
                                                     gdouble panel_a, gdouble panel_b,
                                                     gdouble integral_a, gdouble integral_b,
                                                     gdouble reference_scale,
                                                     guint effective_len)
{
  const guint n          = effective_len;
  const guint output_len = sbilv->edge_cheb_coeffs->len;
  const gdouble alpha    = (panel_b - panel_a) / (integral_b - integral_a);
  const gdouble beta     = (panel_b + panel_a - integral_b - integral_a) / (integral_b - integral_a);
  const gdouble *a       = (gdouble *) sbilv->edge_cheb_coeffs->data;
  gdouble *previous;
  gdouble *current;
  gdouble *next;
  gdouble *b;
  gdouble transformed_scale = 0.0;
  guint degree, i;

  if (sbilv->edge_transform_work_len < 3 * n)
  {
    sbilv->edge_transform_work     = g_realloc_n (sbilv->edge_transform_work, 3 * n, sizeof (gdouble));
    sbilv->edge_transform_work_len = 3 * n;
  }

  previous = sbilv->edge_transform_work;
  current  = previous + n;
  next     = current + n;
  memset (sbilv->edge_transform_work, 0, 3 * n * sizeof (gdouble));
  g_array_set_size (sbilv->cheb_coeffs, output_len);
  b = (gdouble *) sbilv->cheb_coeffs->data;
  memset (b, 0, output_len * sizeof (gdouble));

  previous[0] = 1.0;
  b[0]        = a[0];

  if (n > 1)
  {
    current[0] = beta;
    current[1] = alpha;
    b[0]      += a[1] * beta;
    b[1]      += a[1] * alpha;
  }

  /* Recursively form T_degree(alpha t + beta) in the T_k(t) basis. */
  for (degree = 1; degree + 1 < n; degree++)
  {
    gdouble *tmp;

    memset (next, 0, n * sizeof (gdouble));

    for (i = 0; i <= degree; i++)
    {
      next[i] += 2.0 * beta * current[i] - previous[i];

      if (i == 0)
      {
        next[1] += 2.0 * alpha * current[0];
      }
      else
      {
        next[i - 1] += alpha * current[i];
        next[i + 1] += alpha * current[i];
      }
    }

    for (i = 0; i <= degree + 1; i++)
      b[i] += a[degree + 1] * next[i];

    tmp      = previous;
    previous = current;
    current  = next;
    next     = tmp;
  }

  for (i = 0; i < output_len; i++)
  {
    if (!isfinite (b[i]))
      return FALSE;

    transformed_scale += fabs (b[i]);
  }

  return transformed_scale <= 1.0e4 * MAX (reference_scale, G_MINDOUBLE);
}

static gboolean
_ncm_sbessel_integrator_levin_prepare_extended_rhs (NcmSBesselIntegratorLevin *sbilv,
                                                    NcmSpectral *spectral,
                                                    NcmSBesselIntegratorF F,
                                                    gdouble panel_a, gdouble panel_b,
                                                    gdouble integral_a, gdouble integral_b,
                                                    gdouble k, gpointer user_data)
{
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};
  NcmSBesselIntegratorLevinEdgeExtension extension;
  gdouble reference_scale = 0.0;
  gdouble discarded_scale = 0.0;
  gdouble discard_limit;
  guint effective_len;
  guint i;

  /* Fit only on the caller's interval.  The resulting polynomial supplies a
   * smooth extension over the fixed cell, so callbacks are never evaluated
   * outside their advertised integration domain. */
  ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, &_ncm_sbessel_integrator_levin_wrapper_func,
                                                       integral_a, integral_b,
                                                       sbilv->cheb_min_order, sbilv->cheb_reltol, sbilv->panel_abstol,
                                                       &sbilv->edge_cheb_coeffs, &wrapper);

  for (i = 0; i < sbilv->edge_cheb_coeffs->len; i++)
    reference_scale += fabs (g_array_index (sbilv->edge_cheb_coeffs, gdouble, i));

  /* Remove only roundoff-level tail coefficients before extrapolation.  Even
   * a 1e-16 coefficient can grow enormously under T_n(alpha t + beta). */
  discard_limit = 1.0e-4 * MAX (sbilv->cheb_reltol * reference_scale, sbilv->panel_abstol);
  effective_len = sbilv->edge_cheb_coeffs->len;

  while (effective_len > 1)
  {
    const gdouble next_scale = discarded_scale + fabs (g_array_index (sbilv->edge_cheb_coeffs,
                                                                      gdouble, effective_len - 1));

    if (next_scale > discard_limit)
      break;

    discarded_scale = next_scale;
    effective_len--;
  }

  for (i = effective_len; i < sbilv->edge_cheb_coeffs->len; i++)
    g_array_index (sbilv->edge_cheb_coeffs, gdouble, i) = 0.0;

  /* The affine coefficient transform is quadratic in the polynomial degree. */
  if (effective_len > 513)
    return FALSE;

  extension.coeffs = sbilv->edge_cheb_coeffs;
  extension.a      = integral_a;
  extension.b      = integral_b;

  /* Extrapolated high-order noise can grow exponentially outside [-1, 1].
   * Reject such a cell and use the moving-panel solver instead. */
  for (i = 0; i <= 8; i++)
  {
    const gdouble y     = panel_a + (panel_b - panel_a) * i / 8.0;
    const gdouble value = _ncm_sbessel_integrator_levin_edge_extension_func (&extension, y);

    if (!isfinite (value) || (fabs (value) > 1.0e4 * MAX (reference_scale, G_MINDOUBLE)))
      return FALSE;
  }

  if (!_ncm_sbessel_integrator_levin_transform_edge_coeffs (sbilv,
                                                            panel_a, panel_b, integral_a, integral_b,
                                                            reference_scale, effective_len))
    return FALSE;

  _ncm_sbessel_integrator_levin_build_rhs (sbilv);

  return TRUE;
}

/*
 * Let v_ell(y) = y j_ell(y) and W = v_ell u' - v_ell' u.  The forced ODE
 * gives W' = K(y / k, k) j_ell(y) / k.  We may therefore solve on a larger,
 * fixed cell with any smooth forcing extension: W(integral_b) - W(integral_a)
 * depends only on the forcing between the true bounds.  This is what makes the
 * fixed operator factorization reusable while a or b moves.
 */
static gboolean
_ncm_sbessel_integrator_levin_integrate_extended_panel (NcmSBesselIntegratorLevin *sbilv,
                                                        NcmSpectral *spectral,
                                                        NcmSBesselOdeOperator *operator,
                                                        NcmSBesselIntegratorF F,
                                                        gdouble panel_a, gdouble panel_b,
                                                        gdouble integral_a, gdouble integral_b,
                                                        const gdouble *j_integral_a, const gdouble *j_integral_b,
                                                        gdouble k, guint ell_min, guint ell_max,
                                                        gdouble *result_data, gpointer user_data)
{
  const gdouble *j_a;
  const gdouble *j_b;
  guint ell;

  if (j_integral_a != NULL)
  {
    j_a = j_integral_a;
  }
  else
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, integral_a, sbilv->j_array_a);
    j_a = sbilv->j_array_a;
  }

  if (j_integral_b != NULL)
  {
    j_b = j_integral_b;
  }
  else
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, integral_b, sbilv->j_array_b);
    j_b = sbilv->j_array_b;
  }

  sbilv->panel_abstol = _ncm_sbessel_integrator_levin_panel_abstol (sbilv, result_data,
                                                                    j_a, j_b, integral_a, integral_b,
                                                                    ell_min, ell_max);

  if (!_ncm_sbessel_integrator_levin_prepare_extended_rhs (sbilv, spectral, F,
                                                           panel_a, panel_b, integral_a, integral_b,
                                                           k, user_data))
    return FALSE;

  ncm_sbessel_ode_operator_solve_values (operator, sbilv->rhs,
                                         integral_a, integral_b, &sbilv->values_result);

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    const guint ell_idx   = ell - ell_min;
    const gdouble *values = &g_array_index (sbilv->values_result, gdouble, 4 * ell_idx);
    const gdouble u_a     = values[0];
    const gdouble du_a    = values[1];
    const gdouble u_b     = values[2];
    const gdouble du_b    = values[3];
    const gdouble yj_p_a  = _ncm_sbessel_integrator_levin_yj_deriv (ell, integral_a, j_a);
    const gdouble yj_p_b  = _ncm_sbessel_integrator_levin_yj_deriv (ell, integral_b, j_b);
    gdouble W_a, W_b;

    W_a                   = integral_a * j_a[ell] * du_a - yj_p_a * u_a;
    W_b                   = integral_b * j_b[ell] * du_b - yj_p_b * u_b;
    result_data[ell_idx] += W_b - W_a;

    if (G_UNLIKELY (sbilv->record_panels))
    {
      const NcmSBesselIntegratorLevinPanelRec rec = {integral_a, integral_b, (gint) ell, W_b - W_a};

      g_array_append_val (sbilv->panel_records, rec);
    }
  }

  return TRUE;
}

/**
 * _ncm_sbessel_integrator_levin_integrate_panel:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @a_p_idx: knot index for left endpoint (-1 if not on a knot)
 * @b_p_idx: knot index for right endpoint (-1 if not on a knot)
 * @a_p: lower bound of panel
 * @b_p: upper bound of panel
 * @spectral: spectral methods object
 * @F: integrand function
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @result_data: array to accumulate results
 * @user_data: user data for integrand
 *
 * High-level wrapper that integrates a single panel by:
 * 1. Acquiring panel resources (j_ell arrays and operator)
 * 2. Solving the Levin ODE and accumulating results
 *
 * This orchestrator function provides a clean interface for panel integration
 * while keeping the underlying implementation modular for testing and reuse.
 */
static void
_ncm_sbessel_integrator_levin_integrate_panel (NcmSBesselIntegratorLevin *sbilv,
                                               gint a_p_idx, gint b_p_idx,
                                               gdouble a_p, gdouble b_p,
                                               NcmSpectral *spectral,
                                               NcmSBesselIntegratorF F,
                                               gdouble k,
                                               guint ell_min, guint ell_max,
                                               gdouble *result_data,
                                               gpointer user_data)
{
  const gdouble *j_a_p, *j_b_p;
  NcmSBesselOdeOperator *op;

  op = _ncm_sbessel_integrator_levin_get_panel_resources (sbilv, a_p_idx, b_p_idx,
                                                          a_p, b_p, ell_min, ell_max,
                                                          &j_a_p, &j_b_p);

  _ncm_sbessel_integrator_levin_solve_and_accumulate (sbilv, spectral, op,
                                                      F, a_p, b_p, j_a_p, j_b_p, k,
                                                      ell_min, ell_max, result_data, user_data);
}

static void
_ncm_sbessel_integrator_levin_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

  /* Chain up : start */
  NCM_SBESSEL_INTEGRATOR_CLASS (ncm_sbessel_integrator_levin_parent_class)->set_ell_range (sbi, ell_min, ell_max);

  /* Skip operator preparation during construction - resources aren't ready yet */
  if (!sbilv->constructed)
    return;

  if ((ell_min != sbilv->alloc_ell_min) || (ell_max != sbilv->alloc_ell_max))
  {
    if (ell_max > sbilv->ell_cache_max)
      g_error ("Requested ell_max (%u) exceeds ell_cache_max (%u). "
               "Increase ell_cache_max to enable caching for the requested range.",
               ell_max, sbilv->ell_cache_max);

    _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, ell_min, ell_max);
    _ncm_sbessel_integrator_levin_prepare_knots_operators (sbilv, ell_min, ell_max);
  }
}

static guint
_ncm_sbessel_integrator_levin_get_ell_threshold (NcmSBesselIntegratorLevin *sbilv, gdouble a, gdouble b)
{
  return 1000000.0;

  /* The threshold is based on the upper bound */
  return (guint) floor (b);
}

static void
_ncm_sbessel_integrator_levin_integrate_direct (NcmSBesselIntegratorLevin *sbilv,
                                                const guint ell_min, guint ell_max,
                                                NcmSBesselIntegratorF F, const gdouble a, const gdouble b, gdouble k,
                                                NcmVector *result, gpointer user_data)
{
  const gdouble y_min                      = k * a; /* Transform to y-space */
  const gdouble y_max                      = k * b;
  const guint N                            = GSL_MAX (256, 4 * ell_max);
  const gdouble dy                         = (y_max - y_min) / N;
  gdouble * restrict result_ptr            = ncm_vector_data (result);
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};
  guint i, ell;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);
  /* Initialize direct results to zero */
  memset (result_ptr, 0, sizeof (gdouble) * (ell_max - ell_min + 1));
  ell_max = GSL_MIN (ell_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, y_max));

  /* First term */
  {
    const gdouble fa = _ncm_sbessel_integrator_levin_wrapper_func (&wrapper, y_min) / y_min;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, y_min, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] = fa * sbilv->jl_arr[ell];
    }
  }

  /* Interior terms with alternating weights 4 and 2 */
  for (i = 1; i < N; i++)
  {
    const gdouble y      = y_min + i * dy;
    const gdouble weight = (i % 2 == 1) ? 4.0 : 2.0;
    const gdouble fy     = _ncm_sbessel_integrator_levin_wrapper_func (&wrapper, y) / y;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, y, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] += weight * fy * sbilv->jl_arr[ell];
    }
  }

  /* Last term */
  {
    const gdouble fb = _ncm_sbessel_integrator_levin_wrapper_func (&wrapper, y_max) / y_max;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, y_max, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] += fb * sbilv->jl_arr[ell];
    }
  }

  /* Apply Simpson's rule factor */
  for (ell = ell_min; ell <= ell_max; ell++)
    result_ptr[ell - ell_min] *= dy / 3.0;
}

static void
_ncm_sbessel_integrator_levin_integrate_levin (NcmSBesselIntegratorLevin *sbilv,
                                               const guint ell_min, const guint ell_max,
                                               NcmSBesselIntegratorF F, const gdouble a, const gdouble b, gdouble k,
                                               NcmVector *result, gpointer user_data)
{
  NcmSBesselOdeSolver *solver = sbilv->ode_solver;
  NcmSpectral *spectral       = ncm_sbessel_ode_solver_peek_spectral (solver);
  gdouble *result_data        = ncm_vector_data (result);
  const gdouble y_min         = k * a; /* Transform to y-space */
  const gdouble y_max         = k * b;
  gint first_knot_idx         = -1;
  gint last_knot_idx          = -1;
  const guint n_ell           = ell_max - ell_min + 1;
  gdouble first_knot, last_knot;
  guint i;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);

  /* Initialize Levin results to zero */
  memset (result_data, 0, sizeof (gdouble) * n_ell);

  /* Find knots within [y_min, y_max] in y-space */
  for (i = 0; i < sbilv->knots->len; i++)
  {
    const gdouble knot = g_array_index (sbilv->knots, gdouble, i);

    if ((knot >= y_min) && (knot <= y_max))
    {
      if (first_knot_idx == -1)
        first_knot_idx = i;

      last_knot_idx = i;
    }
  }

  /* Use knots-based paneling if configured, otherwise single panel mode */
  if ((sbilv->knots->len > 1) && (first_knot_idx != -1) && (last_knot_idx != -1))
  {
    first_knot = g_array_index (sbilv->knots, gdouble, first_knot_idx);
    last_knot  = g_array_index (sbilv->knots, gdouble, last_knot_idx);

    /* Handle edge panel at the start [y_min, first_knot > y_min] if needed */
    if (y_min < first_knot)
    {
      const gdouble a_p = y_min;
      const gdouble b_p = g_array_index (sbilv->knots, gdouble, first_knot_idx);

      if (ell_min <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
      {
        if (first_knot_idx > 0)
        {
          const guint panel_idx = first_knot_idx - 1;
          const guint n_cache   = sbilv->ell_cache_max + 1;
          gdouble panel_a, panel_b;
          NcmSBesselOdeOperator *op;

          op = _ncm_sbessel_integrator_levin_get_edge_operator (sbilv, panel_idx, FALSE,
                                                                a_p, b_p, ell_min, ell_max,
                                                                &panel_a, &panel_b);

          if (!_ncm_sbessel_integrator_levin_integrate_extended_panel (sbilv, spectral, op, F,
                                                                       panel_a, panel_b, a_p, b_p,
                                                                       NULL,
                                                                       &sbilv->jl_knots[first_knot_idx * n_cache],
                                                                       k, ell_min, ell_max, result_data, user_data))
            _ncm_sbessel_integrator_levin_integrate_panel (sbilv, -1, first_knot_idx,
                                                           a_p, b_p, spectral, F, k,
                                                           ell_min, ell_max, result_data, user_data);
        }
        else
        {
          _ncm_sbessel_integrator_levin_integrate_panel (sbilv, -1, first_knot_idx,
                                                         a_p, b_p, spectral, F, k,
                                                         ell_min, ell_max, result_data, user_data);
        }
      }
    }

    for (i = first_knot_idx; i < (guint) last_knot_idx; i++)
    {
      const gdouble a_p = g_array_index (sbilv->knots, gdouble, i);
      const gdouble b_p = g_array_index (sbilv->knots, gdouble, i + 1);

      if (ell_min <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
        _ncm_sbessel_integrator_levin_integrate_panel (sbilv, i, i + 1,
                                                       a_p, b_p, spectral, F, k,
                                                       ell_min, ell_max, result_data, user_data);
    }

    if (y_max > last_knot)
    {
      const gdouble a_p = g_array_index (sbilv->knots, gdouble, last_knot_idx);
      const gdouble b_p = y_max;

      if (ell_min <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
      {
        if ((guint) last_knot_idx + 1 < sbilv->knots->len)
        {
          const guint panel_idx = last_knot_idx;
          const guint n_cache   = sbilv->ell_cache_max + 1;
          gdouble panel_a, panel_b;
          NcmSBesselOdeOperator *op;

          op = _ncm_sbessel_integrator_levin_get_edge_operator (sbilv, panel_idx, TRUE,
                                                                a_p, b_p, ell_min, ell_max,
                                                                &panel_a, &panel_b);

          if (!_ncm_sbessel_integrator_levin_integrate_extended_panel (sbilv, spectral, op, F,
                                                                       panel_a, panel_b, a_p, b_p,
                                                                       &sbilv->jl_knots[panel_idx * n_cache],
                                                                       NULL,
                                                                       k, ell_min, ell_max, result_data, user_data))
            _ncm_sbessel_integrator_levin_integrate_panel (sbilv, last_knot_idx, -1,
                                                           a_p, b_p, spectral, F, k,
                                                           ell_min, ell_max, result_data, user_data);
        }
        else
        {
          _ncm_sbessel_integrator_levin_integrate_panel (sbilv, last_knot_idx, -1,
                                                         a_p, b_p, spectral, F, k,
                                                         ell_min, ell_max, result_data, user_data);
        }
      }
    }
  }
  else
  {
    /* No paneling: integrate over full range [y_min, y_max]
     * Note: For single panel mode without knots, we use ode_operator directly
     * rather than temp operators. This is handled by passing -1 for both indices
     * but requires special handling in get_panel_resources. */
    const gdouble *j_a_p, *j_b_p;
    NcmSBesselOdeOperator *op;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, y_min, sbilv->j_array_a);
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, y_max, sbilv->j_array_b);
    j_a_p = sbilv->j_array_a;
    j_b_p = sbilv->j_array_b;

    op = sbilv->ode_operator;
    ncm_sbessel_ode_operator_reset (op, y_min, y_max, ell_min, ell_max);

    _ncm_sbessel_integrator_levin_solve_and_accumulate (sbilv, spectral, op,
                                                        F, y_min, y_max, j_a_p, j_b_p, k,
                                                        ell_min, ell_max, result_data, user_data);
  }
}

static void
_ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi,
                                         NcmSBesselIntegratorF F,
                                         gdouble a, gdouble b,
                                         gdouble k,
                                         NcmVector *result,
                                         gpointer user_data)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);
  const gdouble y_min              = k * a; /* Transform to y-space */
  const gdouble y_max              = k * b;
  guint ell_min, ell_max;
  guint n_ell, ell_threshold;

  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  n_ell         = ell_max - ell_min + 1;
  ell_threshold = _ncm_sbessel_integrator_levin_get_ell_threshold (sbilv, y_min, y_max);

  if (G_UNLIKELY (sbilv->record_panels))
    g_array_set_size (sbilv->panel_records, 0);

  g_assert_cmpuint (ncm_vector_len (result), >=, n_ell);

  /* Ensure resources are allocated */
  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, ell_min, ell_max);

  /* Use direct cubature integration for high ell values */
  if (ell_threshold <= ell_max)
    _ncm_sbessel_integrator_levin_integrate_direct (sbilv, ell_min, ell_max, F, a, b, k, result, user_data);
  else
    _ncm_sbessel_integrator_levin_integrate_levin (sbilv, ell_min, ell_max, F, a, b, k, result, user_data);
}

/**
 * ncm_sbessel_integrator_levin_set_record_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @record: whether to record per-panel diagnostics
 *
 * Enables or disables recording of the individual panel contributions.
 *
 * The integral is assembled as a sum of per-panel boundary terms which cancel
 * heavily, so the attainable relative accuracy is bounded by the cancellation
 * ratio $\sum_p |I_p| / |\sum_p I_p|$ times the machine epsilon. That ratio cannot
 * be recovered from the result alone. With recording enabled the contributions are
 * kept and can be read back with
 * ncm_sbessel_integrator_levin_get_panel_contrib().
 *
 * Recording is off by default and costs nothing when off. The records are cleared
 * at the start of every ncm_sbessel_integrator_integrate() call, so they always
 * describe the most recent one.
 */
void
ncm_sbessel_integrator_levin_set_record_panels (NcmSBesselIntegratorLevin *sbilv, gboolean record)
{
  sbilv->record_panels = record;

  if (!record)
    g_array_set_size (sbilv->panel_records, 0);
}

/**
 * ncm_sbessel_integrator_levin_get_record_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets whether per-panel diagnostics are being recorded.
 *
 * Returns: TRUE if recording is enabled.
 */
gboolean
ncm_sbessel_integrator_levin_get_record_panels (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->record_panels;
}

/**
 * ncm_sbessel_integrator_levin_get_n_panel_records:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the number of panel records from the most recent integration. One record is
 * produced per (panel, multipole) pair.
 *
 * Returns: the number of records.
 */
guint
ncm_sbessel_integrator_levin_get_n_panel_records (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->panel_records->len;
}

#define _NCM_SBILV_REC(sbilv, i) (&g_array_index ((sbilv)->panel_records, NcmSBesselIntegratorLevinPanelRec, (i)))

/**
 * ncm_sbessel_integrator_levin_get_panel_a:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the lower bound, in $y = kx$, of the panel of record @i.
 *
 * Returns: the panel lower bound.
 */
gdouble
ncm_sbessel_integrator_levin_get_panel_a (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->a;
}

/**
 * ncm_sbessel_integrator_levin_get_panel_b:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the upper bound, in $y = kx$, of the panel of record @i.
 *
 * Returns: the panel upper bound.
 */
gdouble
ncm_sbessel_integrator_levin_get_panel_b (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->b;
}

/**
 * ncm_sbessel_integrator_levin_get_panel_ell:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the multipole of record @i.
 *
 * Returns: the multipole.
 */
gint
ncm_sbessel_integrator_levin_get_panel_ell (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->ell;
}

/**
 * ncm_sbessel_integrator_levin_get_panel_contrib:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the contribution of record @i to the total integral. The sum of the
 * contributions over all records with the same multipole is that multipole's result.
 *
 * Returns: the panel contribution.
 */
gdouble
ncm_sbessel_integrator_levin_get_panel_contrib (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->contrib;
}

/**
 * ncm_sbessel_integrator_levin_new:
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 *
 * Creates a new #NcmSBesselIntegratorLevin with default parameters:
 * - y_knots_min = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_Y_KNOTS_MIN
 * - y_knots_max = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_Y_KNOTS_MAX
 * - n_knots = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS
 * - ell_cache_max = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX
 * - reltol = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_RELTOL
 * - cheb_min_order = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER
 * - cheb_reltol = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_RELTOL
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorLevin
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_new (guint ell_min, guint ell_max)
{
  return ncm_sbessel_integrator_levin_new_full (ell_min, ell_max,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_Y_KNOTS_MIN,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_Y_KNOTS_MAX,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_RELTOL,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_RELTOL);
}

/**
 * ncm_sbessel_integrator_levin_new_full:
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @y_knots_min: minimum value for knots in log-spaced grid (set to 0 to disable knots-based paneling)
 * @y_knots_max: maximum value for knots in log-spaced grid (set to 0 to disable knots-based paneling)
 * @n_knots: number of knots in the log-spaced grid (set to 0 to disable knots-based paneling)
 * @ell_cache_max: maximum ell value for precomputed spherical Bessel functions at knots
 * @reltol: relative tolerance for integration
 * @cheb_min_order: minimum order of Chebyshev decomposition for RHS computation
 * @cheb_reltol: relative tolerance for Chebyshev decomposition of integrand
 *
 * Creates a new #NcmSBesselIntegratorLevin with optional knots-based paneling. To
 * disable knots-based paneling and use single panel mode, set @y_knots_min,
 * @y_knots_max, or @n_knots to 0.
 *
 * The @ell_cache_max parameter controls the maximum ell value for which spherical
 * Bessel functions will be precomputed at all knots. For ell values beyond this,
 * spherical Bessel functions will be computed on-the-fly during integration.
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorLevin
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_new_full (guint ell_min, guint ell_max, gdouble y_knots_min, gdouble y_knots_max, guint n_knots, guint ell_cache_max, gdouble reltol, guint cheb_min_order, gdouble cheb_reltol)
{
  NcmDTuple2 ell_range             = NCM_DTUPLE2_STATIC_INIT ((gdouble) ell_min, (gdouble) ell_max);
  NcmSBesselIntegratorLevin *sbilv = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN,
                                                   "ell-range", &ell_range,
                                                   "y-knots-min", y_knots_min,
                                                   "y-knots-max", y_knots_max,
                                                   "n-knots", n_knots,
                                                   "ell-cache-max", ell_cache_max,
                                                   "reltol", reltol,
                                                   "cheb-min-order", cheb_min_order,
                                                   "cheb-reltol", cheb_reltol,
                                                   NULL);

  return sbilv;
}

/**
 * ncm_sbessel_integrator_levin_ref:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Increases the reference count of @sbilv by one.
 *
 * Returns: (transfer full): @sbilv
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_ref (NcmSBesselIntegratorLevin *sbilv)
{
  return g_object_ref (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_free:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Decreases the reference count of @sbilv by one.
 */
void
ncm_sbessel_integrator_levin_free (NcmSBesselIntegratorLevin *sbilv)
{
  g_object_unref (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_clear:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * If @sbilv is different from NULL, decreases the reference count of
 * @sbilv by one and sets @sbilv to NULL.
 */
void
ncm_sbessel_integrator_levin_clear (NcmSBesselIntegratorLevin **sbilv)
{
  g_clear_object (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_set_max_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @max_order: maximum order
 *
 * Sets the maximum order of Clenshaw-Curtis quadrature.
 */
void
ncm_sbessel_integrator_levin_set_max_order (NcmSBesselIntegratorLevin *sbilv, guint max_order)
{
  sbilv->max_order = max_order;
}

/**
 * ncm_sbessel_integrator_levin_get_max_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the maximum order.
 *
 * Returns: the maximum order
 */
guint
ncm_sbessel_integrator_levin_get_max_order (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->max_order;
}

/**
 * ncm_sbessel_integrator_levin_set_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance used to build ODE operators from here on.
 *
 * Operators take the tolerance in force when they are created and keep it for
 * their whole life, so this does not reach operators that already exist. Since
 * the panel operators are built as soon as the multipole range is known, an
 * integrator constructed and then given a tolerance keeps the tolerance it was
 * constructed with. Pass @reltol at construction --
 * ncm_sbessel_integrator_levin_new_full() or the `reltol` property -- to have
 * it apply to every operator.
 *
 * This is deliberate: an operator's accuracy is fixed when it is built, so it
 * cannot be changed underneath a computation that is already using it.
 */
void
ncm_sbessel_integrator_levin_set_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble reltol)
{
  if (sbilv->reltol == reltol)
    return;

  g_assert_cmpfloat (reltol, >, 0.0);

  sbilv->reltol = reltol;
  ncm_sbessel_ode_solver_set_tolerance (sbilv->ode_solver, reltol);
}

/**
 * ncm_sbessel_integrator_levin_get_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the relative tolerance.
 *
 * Returns: the relative tolerance
 */
gdouble
ncm_sbessel_integrator_levin_get_reltol (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->reltol;
}

/**
 * ncm_sbessel_integrator_levin_set_cheb_min_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @cheb_min_order: minimum Chebyshev order
 *
 * Sets the minimum order of Chebyshev decomposition for RHS computation.
 */
void
ncm_sbessel_integrator_levin_set_cheb_min_order (NcmSBesselIntegratorLevin *sbilv, guint cheb_min_order)
{
  sbilv->cheb_min_order = cheb_min_order;
}

/**
 * ncm_sbessel_integrator_levin_get_cheb_min_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the minimum Chebyshev order for RHS computation.
 *
 * Returns: the minimum Chebyshev order
 */
guint
ncm_sbessel_integrator_levin_get_cheb_min_order (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->cheb_min_order;
}

/**
 * ncm_sbessel_integrator_levin_set_cheb_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @cheb_reltol: Chebyshev decomposition relative tolerance
 *
 * Sets the relative tolerance for Chebyshev decomposition of the integrand.
 */
void
ncm_sbessel_integrator_levin_set_cheb_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble cheb_reltol)
{
  sbilv->cheb_reltol = cheb_reltol;
}

/**
 * ncm_sbessel_integrator_levin_get_cheb_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the relative tolerance for Chebyshev decomposition.
 *
 * Returns: the Chebyshev decomposition relative tolerance
 */
gdouble
ncm_sbessel_integrator_levin_get_cheb_reltol (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->cheb_reltol;
}

/**
 * ncm_sbessel_integrator_levin_get_y_knots_min:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the minimum knot value.
 *
 * Returns: the minimum knot value
 */
gdouble
ncm_sbessel_integrator_levin_get_y_knots_min (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->y_knots_min;
}

/**
 * ncm_sbessel_integrator_levin_get_y_knots_max:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the maximum knot value.
 *
 * Returns: the maximum knot value
 */
gdouble
ncm_sbessel_integrator_levin_get_y_knots_max (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->y_knots_max;
}

/**
 * ncm_sbessel_integrator_levin_get_n_knots:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the number of knots.
 *
 * Returns: the number of knots
 */
guint
ncm_sbessel_integrator_levin_get_n_knots (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_knots;
}

/**
 * ncm_sbessel_integrator_levin_get_ell_cache_max:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the maximum ell value for the precomputed spherical Bessel cache.
 *
 * Returns: the maximum ell value for cache
 */
guint
ncm_sbessel_integrator_levin_get_ell_cache_max (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->ell_cache_max;
}

