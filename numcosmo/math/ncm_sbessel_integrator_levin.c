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
 * For low ell values, the method solves the differential equation $x^2 y''(x) + 2x y'(x) + (x^2 -
 * \ell(\ell+1)) y(x) = f(x)$ with boundary conditions $y(x_n) = 0 = y(x_{n+1})$. The
 * integral is then given by $I_n = j_\ell(x_{n+1}) y'(x_{n+1}) - j_\ell(x_n) y'(x_n)$.
 *
 * For high ell values, it uses vector cubature integration where the integrand evaluates
 * $f(x)$ and all $j_\ell(x)$ values simultaneously for efficiency.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_integrator_levin.h"
#include "math/ncm_sbessel_ode_solver.h"
#include "math/ncm_spectral.h"
#include "math/ncm_sf_sbessel.h"
#include "math/ncm_lapack.h"
#include "math/ncm_c.h"
#include "math/ncm_integral_nd.h"
#include "ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_integration.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmSBesselIntegratorLevinArg
{
  NcmSBesselIntegratorLevin *sbilv;
  NcmSBesselIntegratorF F;
  gpointer user_data;
  guint ell_min;
  guint ell_max;
  gdouble *jl_arr;
} NcmSBesselIntegratorLevinArg;

struct _NcmSBesselIntegratorLevin
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  guint max_order;
  gdouble reltol;
  NcmSBesselOdeSolver *ode_solver;
  NcmSBesselOdeOperator *ode_operator;
  NcmSFSBesselArray *sba; /* Allocation tracking */
  guint alloc_max_order;
  guint alloc_ell_min;
  guint alloc_ell_max;
  /* Pre-allocated working arrays */
  GArray *cheb_coeffs;
  GArray *gegen_coeffs;
  GArray *rhs;
  gdouble *j_array_a;
  gdouble *j_array_b;
  GArray *endpoints_result;
  gdouble *jl_arr;
  /* Knots-based paneling */
  gdouble y_knots_min;
  gdouble y_knots_max;
  guint n_knots;
  guint ell_cache_max;                        /* Maximum ell for precomputed j_l at knots */
  GArray *knots;                              /* Log-spaced knots array */
  gdouble *jl_knots;                          /* Precomputed j_l at knots: [n_knots * (ell_cache_max + 1)] */
  GPtrArray *operators;                       /* Operators for each panel between consecutive knots */
  NcmSBesselOdeOperator *ode_operator_temp_a; /* Temporary operator for [a, smallest_knot > a] */
  NcmSBesselOdeOperator *ode_operator_temp_b; /* Temporary operator for [largest_knot < b, b] */
};

enum
{
  PROP_0,
  PROP_MAX_ORDER,
  PROP_RELTOL,
  PROP_Y_KNOTS_MIN,
  PROP_Y_KNOTS_MAX,
  PROP_N_KNOTS,
  PROP_ELL_CACHE_MAX,
};

static void _ncm_sbessel_integrator_levin_prepare_knots_array (NcmSBesselIntegratorLevin *sbilv);
static void _ncm_sbessel_integrator_levin_prepare_ell_cache (NcmSBesselIntegratorLevin *sbilv);
static void _ncm_sbessel_integrator_levin_compute_rhs (NcmSBesselIntegratorLevin *sbilv, NcmSpectral *spectral, NcmSBesselIntegratorF F, gdouble a, gdouble b, gpointer user_data);
static void _ncm_sbessel_integrator_levin_solve_and_accumulate (NcmSBesselIntegratorLevin *sbilv, NcmSpectral *spectral, NcmSBesselOdeOperator *operator, NcmSBesselIntegratorF F, gdouble a_p, gdouble b_p, const gdouble *j_a_p, const gdouble *j_b_p, guint ell_min, guint ell_max, gdouble *result_data, gpointer user_data);
static void _ncm_sbessel_integrator_levin_integrate_panel (NcmSBesselIntegratorLevin *sbilv, gint a_p_idx, gint b_p_idx, gdouble a_p, gdouble b_p, NcmSpectral *spectral, NcmSBesselIntegratorF F, guint ell_min, guint ell_max, gdouble *result_data, gpointer user_data);
static void _ncm_sbessel_integrator_levin_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, NcmVector *result, gpointer user_data);

G_DEFINE_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_levin_init (NcmSBesselIntegratorLevin *sbilv)
{
  sbilv->max_order        = 0;
  sbilv->reltol           = 0.0;
  sbilv->ode_solver       = ncm_sbessel_ode_solver_new ();
  sbilv->ode_operator     = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, 2, 2);
  sbilv->sba              = ncm_sf_sbessel_array_new ();
  sbilv->alloc_max_order  = 0;
  sbilv->alloc_ell_min    = -1;
  sbilv->alloc_ell_max    = -1;
  sbilv->cheb_coeffs      = NULL;
  sbilv->gegen_coeffs     = NULL;
  sbilv->rhs              = NULL;
  sbilv->j_array_a        = NULL;
  sbilv->j_array_b        = NULL;
  sbilv->endpoints_result = NULL;
  sbilv->jl_arr           = NULL;
  /* Knots-based paneling */
  sbilv->y_knots_min         = 0.0;
  sbilv->y_knots_max         = 0.0;
  sbilv->n_knots             = 0;
  sbilv->ell_cache_max       = 0;
  sbilv->knots               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->jl_knots            = NULL;
  sbilv->operators           = NULL;
  sbilv->ode_operator_temp_a = NULL;
  sbilv->ode_operator_temp_b = NULL;

  ncm_sbessel_ode_solver_set_tolerance (sbilv->ode_solver, 1.0e-14);
}

static void
_ncm_sbessel_integrator_levin_dispose (GObject *object)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  ncm_sbessel_ode_solver_clear (&sbilv->ode_solver);
  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator);
  ncm_sf_sbessel_array_clear (&sbilv->sba);
  g_clear_pointer (&sbilv->cheb_coeffs, g_array_unref);
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

  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  /* Prepare knots array and precompute spherical Bessel functions at construction time
   * since knots parameters and ell_cache_max are CONSTRUCT_ONLY */
  _ncm_sbessel_integrator_levin_prepare_knots_array (sbilv);
  _ncm_sbessel_integrator_levin_prepare_ell_cache (sbilv);
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
    guint i;

    for (i = 0; i < sbilv->n_knots; i++)
    {
      const gdouble ln_y = ln_y_min + L * i / (sbilv->n_knots - 1);
      const gdouble y    = exp (ln_y);

      g_array_index (sbilv->knots, gdouble, i) = y;
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

      sbilv->ode_operator_temp_a = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_b = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, ell_min, ell_max);
    }
    else if (need_reset)
    {
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
    }

    sbilv->alloc_ell_min = ell_min;
    sbilv->alloc_ell_max = ell_max;
  }
}

/**
 * _ncm_sbessel_integrator_levin_compute_rhs:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @spectral: spectral methods object
 * @F: integrand function
 * @a: lower integration bound
 * @b: upper integration bound
 * @user_data: user data for integrand
 *
 * Computes the RHS for the Levin ODE by:
 * 1. Computing Chebyshev coefficients for f(x)
 * 2. Converting to Gegenbauer C^(2) basis
 * 3. Setting up RHS with homogeneous boundary conditions
 */
static void
_ncm_sbessel_integrator_levin_compute_rhs (NcmSBesselIntegratorLevin *sbilv,
                                           NcmSpectral *spectral,
                                           NcmSBesselIntegratorF F,
                                           gdouble a, gdouble b,
                                           gpointer user_data)
{
  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a, b, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);
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
 * If a_p_idx == -1, computes j_a_p and uses temp_a operator with reset.
 * If b_p_idx == -1, computes j_b_p and uses temp_b operator with reset.
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

  /* Get j_a_p: compute if a_p_idx == -1, else use cache */
  if (a_p_idx < 0)
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, a_p, sbilv->j_array_a);
    *j_a_p_out = sbilv->j_array_a;
  }
  else
  {
    *j_a_p_out = &sbilv->jl_knots[a_p_idx * n_ell];
  }

  /* Get j_b_p: compute if b_p_idx == -1, else use cache */
  if (b_p_idx < 0)
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, b_p, sbilv->j_array_b);
    *j_b_p_out = sbilv->j_array_b;
  }
  else
  {
    *j_b_p_out = &sbilv->jl_knots[b_p_idx * n_ell];
  }

  /* Select operator and reset if needed */
  if (a_p_idx < 0)
  {
    op = sbilv->ode_operator_temp_a;
    ncm_sbessel_ode_operator_reset (op, a_p, b_p, ell_min, ell_max);
  }
  else if (b_p_idx < 0)
  {
    op = sbilv->ode_operator_temp_b;
    ncm_sbessel_ode_operator_reset (op, a_p, b_p, ell_min, ell_max);
  }
  else
  {
    op = g_ptr_array_index (sbilv->operators, a_p_idx);
  }

  return op;
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
                                                    guint ell_min, guint ell_max,
                                                    gdouble *result_data,
                                                    gpointer user_data)
{
  guint ell;

  _ncm_sbessel_integrator_levin_compute_rhs (sbilv, spectral, F, a_p, b_p, user_data);
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
  }
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
                                                      F, a_p, b_p, j_a_p, j_b_p,
                                                      ell_min, ell_max, result_data, user_data);
}

static void
_ncm_sbessel_integrator_levin_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

  /* Chain up : start */
  NCM_SBESSEL_INTEGRATOR_CLASS (ncm_sbessel_integrator_levin_parent_class)->set_ell_range (sbi, ell_min, ell_max);

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
  /* The threshold is based on the upper bound */
  return (guint) floor (b);
}

static void
_ncm_sbessel_integrator_levin_integrate_direct (NcmSBesselIntegratorLevin *sbilv,
                                                const guint ell_min, guint ell_max,
                                                NcmSBesselIntegratorF F, const gdouble a, const gdouble b,
                                                NcmVector *result, gpointer user_data)
{
  const guint N                 = GSL_MAX (256, 4 * ell_max);
  const gdouble dx              = (b - a) / N;
  gdouble * restrict result_ptr = ncm_vector_data (result);
  guint i, ell;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);
  /* Initialize direct results to zero */
  memset (result_ptr, 0, sizeof (gdouble) * (ell_max - ell_min + 1));
  ell_max = GSL_MIN (ell_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b));

  /* First term */
  {
    const gdouble fa = F (user_data, a) / a;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, a, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] = fa * sbilv->jl_arr[ell];
    }
  }

  /* Interior terms with alternating weights 4 and 2 */
  for (i = 1; i < N; i++)
  {
    const gdouble x      = a + i * dx;
    const gdouble weight = (i % 2 == 1) ? 4.0 : 2.0;
    const gdouble fx     = F (user_data, x) / x;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, x, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] += weight * fx * sbilv->jl_arr[ell];
    }
  }

  /* Last term */
  {
    const gdouble fb = F (user_data, b) / b;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, b, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] += fb * sbilv->jl_arr[ell];
    }
  }

  /* Apply Simpson's rule factor */
  for (ell = ell_min; ell <= ell_max; ell++)
    result_ptr[ell - ell_min] *= dx / 3.0;
}

static void
_ncm_sbessel_integrator_levin_integrate_levin (NcmSBesselIntegratorLevin *sbilv,
                                               const guint ell_min, const guint ell_max,
                                               NcmSBesselIntegratorF F, const gdouble a, const gdouble b,
                                               NcmVector *result, gpointer user_data)
{
  NcmSBesselOdeSolver *solver = sbilv->ode_solver;
  NcmSpectral *spectral       = ncm_sbessel_ode_solver_peek_spectral (solver);
  gdouble *result_data        = ncm_vector_data (result);
  gint first_knot_idx         = -1;
  gint last_knot_idx          = -1;
  gdouble first_knot, last_knot;
  guint i;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);

  /* Initialize Levin results to zero */
  memset (result_data, 0, sizeof (gdouble) * (ell_max - ell_min + 1));

  /* Find knots within [a, b] */
  for (i = 0; i < sbilv->knots->len; i++)
  {
    const gdouble knot = g_array_index (sbilv->knots, gdouble, i);

    if ((knot >= a) && (knot <= b))
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

    /* Handle edge panel at the start [a, first_knot > a] if needed */
    if (a < first_knot)
    {
      const gdouble a_p = a;
      const gdouble b_p = g_array_index (sbilv->knots, gdouble, first_knot_idx);

      if (ell_min <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
        _ncm_sbessel_integrator_levin_integrate_panel (sbilv, -1, first_knot_idx,
                                                       a_p, b_p, spectral, F,
                                                       ell_min, ell_max, result_data, user_data);
    }

    for (i = first_knot_idx; i < (guint) last_knot_idx; i++)
    {
      const gdouble a_p = g_array_index (sbilv->knots, gdouble, i);
      const gdouble b_p = g_array_index (sbilv->knots, gdouble, i + 1);

      if (ell_min <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
        _ncm_sbessel_integrator_levin_integrate_panel (sbilv, i, i + 1,
                                                       a_p, b_p, spectral, F,
                                                       ell_min, ell_max, result_data, user_data);
    }

    if (b > last_knot)
    {
      const gdouble a_p = g_array_index (sbilv->knots, gdouble, last_knot_idx);
      const gdouble b_p = b;

      if (ell_min <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
        _ncm_sbessel_integrator_levin_integrate_panel (sbilv, last_knot_idx, -1,
                                                       a_p, b_p, spectral, F,
                                                       ell_min, ell_max, result_data, user_data);
    }
  }
  else
  {
    /* No paneling: integrate over full range [a, b]
     * Note: For single panel mode without knots, we use ode_operator directly
     * rather than temp operators. This is handled by passing -1 for both indices
     * but requires special handling in get_panel_resources. */
    const gdouble *j_a_p, *j_b_p;
    NcmSBesselOdeOperator *op;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, a, sbilv->j_array_a);
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, b, sbilv->j_array_b);
    j_a_p = sbilv->j_array_a;
    j_b_p = sbilv->j_array_b;

    op = sbilv->ode_operator;
    ncm_sbessel_ode_operator_reset (op, a, b, ell_min, ell_max);

    _ncm_sbessel_integrator_levin_solve_and_accumulate (sbilv, spectral, op,
                                                        F, a, b, j_a_p, j_b_p,
                                                        ell_min, ell_max, result_data, user_data);
  }
}

static void
_ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi,
                                         NcmSBesselIntegratorF F,
                                         gdouble a, gdouble b,
                                         NcmVector *result,
                                         gpointer user_data)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);
  guint ell_min, ell_max;
  guint n_ell, ell_threshold;

  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  n_ell         = ell_max - ell_min + 1;
  ell_threshold = _ncm_sbessel_integrator_levin_get_ell_threshold (sbilv, a, b);

  g_assert_cmpuint (ncm_vector_len (result), ==, n_ell);

  /* Ensure resources are allocated */
  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, ell_min, ell_max);

  /* Use direct cubature integration for high ell values */
  if (ell_threshold <= ell_max)
    _ncm_sbessel_integrator_levin_integrate_direct (sbilv, ell_min, ell_max, F, a, b, result, user_data);
  else
    _ncm_sbessel_integrator_levin_integrate_levin (sbilv, ell_min, ell_max, F, a, b, result, user_data);
}

/**
 * ncm_sbessel_integrator_levin_new:
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @y_knots_min: minimum value for knots in log-spaced grid (set to 0 to disable knots-based paneling)
 * @y_knots_max: maximum value for knots in log-spaced grid (set to 0 to disable knots-based paneling)
 * @n_knots: number of knots in the log-spaced grid (set to 0 to disable knots-based paneling)
 * @ell_cache_max: maximum ell value for precomputed spherical Bessel functions at knots
 *
 * Creates a new #NcmSBesselIntegratorLevin with optional knots-based paneling.
 * To disable knots-based paneling and use single panel mode, set @y_knots_min,
 * @y_knots_max, or @n_knots to 0.
 *
 * The @ell_cache_max parameter controls the maximum ell value for which spherical
 * Bessel functions will be precomputed at all knots. For ell values beyond this,
 * spherical Bessel functions will be computed on-the-fly during integration.
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorLevin
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_new (guint ell_min, guint ell_max, gdouble y_knots_min, gdouble y_knots_max, guint n_knots, guint ell_cache_max)
{
  NcmSBesselIntegratorLevin *sbilv = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN,
                                                   "ell_min", ell_min,
                                                   "ell_max", ell_max,
                                                   "y-knots-min", y_knots_min,
                                                   "y-knots-max", y_knots_max,
                                                   "n-knots", n_knots,
                                                   "ell-cache-max", ell_cache_max,
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
 * Sets the relative tolerance for convergence.
 */
void
ncm_sbessel_integrator_levin_set_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble reltol)
{
  sbilv->reltol = reltol;
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

