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
  guint lmin;
  guint lmax;
  gdouble *jl_arr;
} NcmSBesselIntegratorLevinArg;

struct _NcmSBesselIntegratorLevin
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  guint max_order;
  gdouble reltol;
  guint n_panels;
  NcmSBesselOdeSolver *ode_solver;
  NcmSBesselOdeOperator *ode_operator;
  NcmSFSBesselArray *sba; /* Allocation tracking */
  guint alloc_max_order;
  gint alloc_lmin;
  gint alloc_lmax;
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
  GArray *knots;                              /* Log-spaced knots array */
  GPtrArray *operators;                       /* Operators for each panel between consecutive knots */
  NcmSBesselOdeOperator *ode_operator_temp_a; /* Temporary operator for [a, smallest_knot > a] */
  NcmSBesselOdeOperator *ode_operator_temp_b; /* Temporary operator for [largest_knot < b, b] */
};

enum
{
  PROP_0,
  PROP_MAX_ORDER,
  PROP_RELTOL,
  PROP_N_PANELS,
  PROP_Y_KNOTS_MIN,
  PROP_Y_KNOTS_MAX,
  PROP_N_KNOTS,
};

static void ncm_sbessel_integrator_levin_direct_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void ncm_sbessel_integrator_levin_direct_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

NCM_INTEGRAL_ND_DEFINE_TYPE (NCM, SBESSEL_INTEGRATOR_LEVIN_DIRECT, NcmSBesselIntegratorLevinDirect, ncm_sbessel_integrator_levin_direct, ncm_sbessel_integrator_levin_direct_dim, ncm_sbessel_integrator_levin_direct_integ, NcmSBesselIntegratorLevinArg);

G_DEFINE_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_levin_init (NcmSBesselIntegratorLevin *sbilv)
{
  sbilv->max_order        = 0;
  sbilv->reltol           = 0.0;
  sbilv->n_panels         = 1;
  sbilv->ode_solver       = ncm_sbessel_ode_solver_new ();
  sbilv->ode_operator     = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, 2, 2);
  sbilv->sba              = ncm_sf_sbessel_array_new ();
  sbilv->alloc_max_order  = 0;
  sbilv->alloc_lmin       = -1;
  sbilv->alloc_lmax       = -1;
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
  sbilv->knots               = NULL;
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

static void _ncm_sbessel_integrator_levin_prepare (NcmSBesselIntegrator *sbi);
static gdouble _ncm_sbessel_integrator_levin_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gint ell, gpointer user_data);
static void _ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, NcmVector *result, gpointer user_data);

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
    case PROP_N_PANELS:
      ncm_sbessel_integrator_levin_set_n_panels (sbilv, g_value_get_uint (value));
      break;
    case PROP_Y_KNOTS_MIN:
      ncm_sbessel_integrator_levin_set_y_knots_min (sbilv, g_value_get_double (value));
      break;
    case PROP_Y_KNOTS_MAX:
      ncm_sbessel_integrator_levin_set_y_knots_max (sbilv, g_value_get_double (value));
      break;
    case PROP_N_KNOTS:
      ncm_sbessel_integrator_levin_set_n_knots (sbilv, g_value_get_uint (value));
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
    case PROP_N_PANELS:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_n_panels (sbilv));
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
   * NcmSBesselIntegratorLevin:n-panels:
   *
   * Number of panels to divide the integration range in log space. Use 1 to disable
   * paneling (default). Values greater than 1 will divide the integration range [a, b]
   * into n_panels in logarithmic space, integrating each panel separately and summing the results.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_PANELS,
                                   g_param_spec_uint ("n-panels",
                                                      NULL,
                                                      "Number of integration panels",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:y-knots-min:
   *
   * Minimum value for knots in log-spaced grid. Set to 0 to disable knots-based paneling.
   */
  g_object_class_install_property (object_class,
                                   PROP_Y_KNOTS_MIN,
                                   g_param_spec_double ("y-knots-min",
                                                        NULL,
                                                        "Minimum knot value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:y-knots-max:
   *
   * Maximum value for knots in log-spaced grid. Set to 0 to disable knots-based paneling.
   */
  g_object_class_install_property (object_class,
                                   PROP_Y_KNOTS_MAX,
                                   g_param_spec_double ("y-knots-max",
                                                        NULL,
                                                        "Maximum knot value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:n-knots:
   *
   * Number of knots in the log-spaced grid. The knots will be equally spaced in
   * log space between y-knots-min and y-knots-max. Set to 0 to disable knots-based paneling.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_KNOTS,
                                   g_param_spec_uint ("n-knots",
                                                      NULL,
                                                      "Number of knots",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->prepare       = &_ncm_sbessel_integrator_levin_prepare;
  parent_class->integrate_ell = &_ncm_sbessel_integrator_levin_integrate_ell;
  parent_class->integrate     = &_ncm_sbessel_integrator_levin_integrate;
}

static void
_ncm_sbessel_integrator_levin_ensure_prepared (NcmSBesselIntegratorLevin *sbilv, guint max_order, gint lmin, gint lmax)
{
  gboolean need_realloc = FALSE;

  /* Check if reallocation is needed */
  if (sbilv->alloc_max_order != max_order)
    need_realloc = TRUE;

  if ((sbilv->alloc_lmin != lmin) || (sbilv->alloc_lmax != lmax))
    need_realloc = TRUE;

  if (!need_realloc)
    return;

  /* Free existing allocations */
  g_clear_pointer (&sbilv->cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->gegen_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->rhs, g_array_unref);
  g_clear_pointer (&sbilv->endpoints_result, g_array_unref);

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

  /* Allocate arrays for spectral coefficients */
  sbilv->cheb_coeffs  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order);
  sbilv->gegen_coeffs = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order);
  sbilv->rhs          = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order + 2);

  /* Allocate arrays for spherical Bessel functions */
  if (lmax >= 0)
  {
    sbilv->j_array_a = g_new0 (gdouble, lmax + 1);
    sbilv->j_array_b = g_new0 (gdouble, lmax + 1);
    sbilv->jl_arr    = g_new (gdouble, lmax + 1);
  }

  /* Allocate result matrix for batched endpoint computation (max block size is 8) */
  sbilv->endpoints_result = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 8 * 3);

  /* Update allocation tracking */
  sbilv->alloc_max_order = max_order;
  sbilv->alloc_lmin      = lmin;
  sbilv->alloc_lmax      = lmax;
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
  /* Only prepare if knots are properly configured and not already prepared */
  if ((sbilv->n_knots < 2) || (sbilv->y_knots_min <= 0.0) || (sbilv->y_knots_max <= sbilv->y_knots_min))
    return;

  if (sbilv->knots != NULL)
    return;  /* Already prepared */

  /* Create log-spaced knots array */
  sbilv->knots = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), sbilv->n_knots);

  {
    const gdouble ln_y_min = log (sbilv->y_knots_min);
    const gdouble ln_y_max = log (sbilv->y_knots_max);
    const gdouble L        = ln_y_max - ln_y_min;
    guint i;

    for (i = 0; i < sbilv->n_knots; i++)
    {
      const gdouble ln_y = ln_y_min + L * i / (sbilv->n_knots - 1);
      const gdouble y    = exp (ln_y);

      g_array_append_val (sbilv->knots, y);
    }

    printf ("Prepared %u knots from %e to %e\n", sbilv->n_knots, sbilv->y_knots_min, sbilv->y_knots_max);
  }
}

/**
 * _ncm_sbessel_integrator_levin_prepare_knots_operators:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 *
 * Prepares ODE operators for the knots-based paneling system:
 * - Pre-allocated ODE operators for each panel between consecutive knots
 * - Two temporary operators for edge panels [a, smallest_knot > a] and [largest_knot < b, b]
 *
 * Uses ncm_sbessel_ode_operator_reset() to efficiently update operators when ell range changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_knots_operators (NcmSBesselIntegratorLevin *sbilv, gint lmin, gint lmax)
{
  guint i;

  /* Only prepare if knots array exists */
  if (sbilv->knots == NULL)
    return;

  /* Check if operators need to be created or reset */
  const gboolean need_create = (sbilv->operators == NULL);
  const gboolean need_reset  = !need_create && ((sbilv->alloc_lmin != lmin) || (sbilv->alloc_lmax != lmax));

  if (!need_create && !need_reset)
    return;  /* Operators already prepared for this ell range */

  if (need_create)
  {
    /* First time: create all operators */
    sbilv->operators = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_sbessel_ode_operator_unref);

    for (i = 0; i < sbilv->n_knots - 1; i++)
    {
      const gdouble y_a = g_array_index (sbilv->knots, gdouble, i);
      const gdouble y_b = g_array_index (sbilv->knots, gdouble, i + 1);
      NcmSBesselOdeOperator *op;

      op = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, y_a, y_b, lmin, lmax);
      g_ptr_array_add (sbilv->operators, op);
    }

    /* Create temporary operators for edge panels */
    sbilv->ode_operator_temp_a = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, lmin, lmax);
    sbilv->ode_operator_temp_b = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, lmin, lmax);

    printf ("Created %u knots-based operators for ell=[%d, %d]\n", sbilv->n_knots - 1, lmin, lmax);
  }
  else if (need_reset)
  {
    /* Operators exist but ell range changed: reset them efficiently */
    for (i = 0; i < sbilv->operators->len; i++)
    {
      NcmSBesselOdeOperator *op = g_ptr_array_index (sbilv->operators, i);
      const gdouble y_a         = g_array_index (sbilv->knots, gdouble, i);
      const gdouble y_b         = g_array_index (sbilv->knots, gdouble, i + 1);

      ncm_sbessel_ode_operator_reset (op, y_a, y_b, lmin, lmax);
    }

    /* Reset temporary operators */
    ncm_sbessel_ode_operator_reset (sbilv->ode_operator_temp_a, 0.0, 1.0, lmin, lmax);
    ncm_sbessel_ode_operator_reset (sbilv->ode_operator_temp_b, 0.0, 1.0, lmin, lmax);

    printf ("Reset %u knots-based operators for ell=[%d, %d]\n", sbilv->n_knots - 1, lmin, lmax);
  }
}

static void
_ncm_sbessel_integrator_levin_prepare (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);
  const guint lmin                 = ncm_sbessel_integrator_get_lmin (sbi);
  const guint lmax                 = ncm_sbessel_integrator_get_lmax (sbi);

  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, lmin, lmax);
  _ncm_sbessel_integrator_levin_prepare_knots_array (sbilv);
  _ncm_sbessel_integrator_levin_prepare_knots_operators (sbilv, lmin, lmax);
}

static gdouble
_ncm_sbessel_integrator_levin_integrate_ell (NcmSBesselIntegrator *sbi,
                                             NcmSBesselIntegratorF F,
                                             gdouble a, gdouble b,
                                             gint ell,
                                             gpointer user_data)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);
  NcmSBesselOdeSolver *solver      = sbilv->ode_solver;
  NcmSpectral *spectral            = ncm_sbessel_ode_solver_peek_spectral (solver);

  printf ("Integrating ell = %u single panel\n", ell);

  /* Ensure resources are allocated */
  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, ell, ell);

  /* Set the interval and l value for this integration */
  ncm_sbessel_ode_operator_reset (sbilv->ode_operator, a, b, ell, ell);

  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a, b, 2, 1.0e-15, &sbilv->cheb_coeffs, user_data);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

  g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
  {
    gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
    const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

    rhs_data[0] = 0.0; /* BC at x=a (t=-1) */
    rhs_data[1] = 0.0; /* BC at x=b (t=+1) */

    /* Copy Gegenbauer coefficients to RHS starting from index 2 */
    memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
  }

  ncm_sbessel_ode_operator_solve_endpoints (sbilv->ode_operator, sbilv->rhs, &sbilv->endpoints_result);

  {
    const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, 0);
    const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, 1);
    const gdouble j_l_a     = gsl_sf_bessel_jl (ell, a);
    const gdouble j_l_b     = gsl_sf_bessel_jl (ell, b);
    const gdouble integral  = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a;

    return integral;
  }
}

static void
ncm_sbessel_integrator_levin_direct_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcmSBesselIntegratorLevinDirect *direct = NCM_SBESSEL_INTEGRATOR_LEVIN_DIRECT (intnd);
  NcmSBesselIntegratorLevinArg *arg       = &direct->data;

  *dim  = 1;
  *fdim = arg->lmax - arg->lmin + 1;
}

static void
ncm_sbessel_integrator_levin_direct_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcmSBesselIntegratorLevinDirect *direct = NCM_SBESSEL_INTEGRATOR_LEVIN_DIRECT (intnd);
  NcmSBesselIntegratorLevinArg *arg       = &direct->data;
  NcmSBesselIntegratorLevin *sbilv        = arg->sbilv;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble xi   = ncm_vector_fast_get (x, i);
    const gdouble f_xi = arg->F (arg->user_data, xi) / xi;

    /* Evaluate all spherical Bessel functions at once */
    ncm_sf_sbessel_array_eval (sbilv->sba, arg->lmax, xi, arg->jl_arr);

    /* Compute integrand for all ells */
    for (j = 0; j < fdim; j++)
    {
      const guint ell   = arg->lmin + j;
      const gdouble res = f_xi * arg->jl_arr[ell];

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
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
                                                const guint lmin, const guint ell_direct_min, guint ell_direct_max,
                                                NcmSBesselIntegratorF F, const gdouble a, const gdouble b,
                                                NcmVector *result, gpointer user_data)
{
  NcmSBesselIntegratorLevinDirect *direct = g_object_new (ncm_sbessel_integrator_levin_direct_get_type (), NULL);
  NcmSBesselIntegratorLevinArg *arg       = &direct->data;
  NcmIntegralND *intnd                    = NCM_INTEGRAL_ND (direct);
  gdouble x_min_val                       = a;
  gdouble x_max_val                       = b;
  NcmVector *x_min                        = ncm_vector_new_data_static (&x_min_val, 1, 1);
  NcmVector *x_max                        = ncm_vector_new_data_static (&x_max_val, 1, 1);
  const guint size                        = ell_direct_max - ell_direct_min + 1;
  NcmVector *result_sub                   = ncm_vector_get_subvector (result, ell_direct_min - lmin, size);
  NcmVector *err                          = ncm_vector_new (size);

  /* Setup argument structure */
  arg->sbilv     = sbilv;
  arg->F         = F;
  arg->user_data = user_data;
  arg->lmin      = ell_direct_min;
  arg->lmax      = ell_direct_max;
  arg->jl_arr    = sbilv->jl_arr;

  /* Configure integrator */
  ncm_integral_nd_set_reltol (intnd, sbilv->reltol);
  ncm_integral_nd_set_abstol (intnd, 1.0e-50);
  ncm_integral_nd_set_method (intnd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P);

  ncm_integral_nd_eval (intnd, x_min, x_max, result_sub, err);

  /* Cleanup */
  ncm_vector_free (x_min);
  ncm_vector_free (x_max);
  ncm_vector_free (result_sub);
  ncm_vector_free (err);
  g_object_unref (direct);
}

static void
_ncm_sbessel_integrator_levin_integrate_levin (NcmSBesselIntegratorLevin *sbilv,
                                               const guint lmin, const guint ell_levin_min, const guint ell_levin_max,
                                               NcmSBesselIntegratorF F, const gdouble a, const gdouble b,
                                               NcmVector *result, gpointer user_data)
{
  NcmSBesselOdeSolver *solver = sbilv->ode_solver;
  NcmSpectral *spectral       = ncm_sbessel_ode_solver_peek_spectral (solver);
  gdouble *result_data        = ncm_vector_data (result);

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);

  /* Initialize Levin results to zero */
  memset (&result_data[ell_levin_min - lmin], 0, sizeof (gdouble) * (ell_levin_max - ell_levin_min + 1));

  /* Panel dynamics: divide integration range if n_panels > 1 */
  if (sbilv->n_panels == 1)
  {
    /* No paneling: integrate over full range [a, b] */
    const gint my_lmax = GSL_MIN (ell_levin_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b));

    ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, a, sbilv->j_array_a);
    ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, b, sbilv->j_array_b);

    /* Step 1: Compute Chebyshev coefficients for f(x) */
    ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a, b, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);

    /* Step 2: Convert to Gegenbauer C^(2) basis */
    ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

    printf ("Computed Gegenbauer with len %u\n", sbilv->gegen_coeffs->len);
    fflush (stdout);

    /* Step 3: Set up RHS with homogeneous boundary conditions */
    g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
    {
      gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
      const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

      rhs_data[0] = 0.0; /* BC at x=a (t=-1) */
      rhs_data[1] = 0.0; /* BC at x=b (t=+1) */

      /* Copy Gegenbauer coefficients to RHS starting from index 2 */
      memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
    }

    /* Step 4-6: Process l values in blocks for better cache locality */
    {
      const guint block_size = 8;
      gint l_start;

      for (l_start = ell_levin_min; l_start <= (gint) ell_levin_max; l_start += block_size)
      {
        const gint l_end = GSL_MIN (l_start + block_size - 1, ell_levin_max);
        gint l;

        ncm_sbessel_ode_operator_reset (sbilv->ode_operator, a, b, ell_levin_min, ell_levin_max);

        /* Solve for all l values in this block using batched solver */
        printf ("Solving ODE for l = %d to %d\n", l_start, l_end);
        fflush (stdout);
        ncm_sbessel_ode_operator_solve_endpoints (sbilv->ode_operator, sbilv->rhs, &sbilv->endpoints_result);
        fflush (stderr);
        fflush (stdout);
        printf ("Done solving ODE for l = %d to %d\n", l_start, l_end);
        fflush (stdout);

        /* Extract derivatives and compute integrals */
        for (l = l_start; l <= l_end; l++)
        {
          const gint l_idx        = l - lmin;
          const gint block_idx    = l - l_start;
          const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 0);
          const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 1);

          if (l <= my_lmax)
          {
            const gdouble j_l_a = sbilv->j_array_a[l];
            const gdouble j_l_b = sbilv->j_array_b[l];

            result_data[l_idx] = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a;
          }
        }
      }
    }
  }
  else if ((sbilv->knots != NULL) && (sbilv->knots->len > 0))
  {
    /* Knots-based panel mode: use pre-computed knots and operators */
    guint i;
    gint first_knot_idx = -1;
    gint last_knot_idx  = -1;

    /* Find knots within [a, b] */
    for (i = 0; i < sbilv->knots->len; i++)
    {
      const gdouble knot = g_array_index (sbilv->knots, gdouble, i);

      if ((knot > a) && (first_knot_idx == -1))
        first_knot_idx = i;

      if (knot < b)
        last_knot_idx = i;
    }

    printf ("Knots-based panel mode: range [%e, %e], found knots from %d to %d\n",
            a, b, first_knot_idx, last_knot_idx);
    fflush (stdout);

    /* Handle edge panel at the start [a, first_knot > a] if needed */
    if (first_knot_idx > 0)
    {
      const gdouble a_p  = a;
      const gdouble b_p  = g_array_index (sbilv->knots, gdouble, first_knot_idx);
      const gint my_lmax = GSL_MIN (ell_levin_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p));

      printf ("Edge panel (start): integrating [%e, %e]\n", a_p, b_p);
      fflush (stdout);

      ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, a_p, sbilv->j_array_a);
      ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, b_p, sbilv->j_array_b);

      /* Compute Chebyshev coefficients */
      ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a_p, b_p, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);

      /* Convert to Gegenbauer C^(2) basis */
      ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

      printf ("Edge panel (start): Computed Gegenbauer with len %u\n", sbilv->gegen_coeffs->len);
      fflush (stdout);

      /* Set up RHS */
      g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
      {
        gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
        const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

        rhs_data[0] = 0.0;
        rhs_data[1] = 0.0;
        memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
      }

      /* Process l values in blocks using temp operator */
      {
        const guint block_size = 8;
        gint l_start;

        for (l_start = ell_levin_min; l_start <= (gint) ell_levin_max; l_start += block_size)
        {
          const gint l_end = GSL_MIN (l_start + block_size - 1, ell_levin_max);
          gint l;

          ncm_sbessel_ode_operator_reset (sbilv->ode_operator_temp_a, a_p, b_p, ell_levin_min, ell_levin_max);
          ncm_sbessel_ode_operator_solve_endpoints (sbilv->ode_operator_temp_a, sbilv->rhs, &sbilv->endpoints_result);
          printf ("Edge panel (start): Solved ODE for l = %d to %d, solution order %ld\n", l_start, l_end,
                  ncm_sbessel_ode_operator_get_n_cols (sbilv->ode_operator_temp_a));

          for (l = l_start; l <= l_end; l++)
          {
            const gint l_idx        = l - lmin;
            const gint block_idx    = l - l_start;
            const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 0);
            const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 1);

            if (l <= my_lmax)
            {
              const gdouble j_l_a = sbilv->j_array_a[l];
              const gdouble j_l_b = sbilv->j_array_b[l];

              result_data[l_idx] += b_p * j_l_b * y_prime_b - a_p * j_l_a * y_prime_a;
            }
          }
        }
      }

      printf ("Edge panel (start): Complete\n");
      fflush (stdout);
    }

    /* Process full panels between consecutive knots */
    if ((first_knot_idx >= 0) && (last_knot_idx >= 0))
    {
      for (i = first_knot_idx; i <= (guint) last_knot_idx; i++)
      {
        const gdouble a_p        = g_array_index (sbilv->knots, gdouble, i);
        const gdouble b_p        = (i + 1 < sbilv->knots->len) ? g_array_index (sbilv->knots, gdouble, i + 1) : b;
        const gint my_lmax       = GSL_MIN (ell_levin_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p));
        const guint operator_idx = i - first_knot_idx;

        /* Skip if this panel is beyond our integration range */
        if (a_p >= b)
          break;

        /* Clip to integration range */
        const gdouble a_panel = GSL_MAX (a_p, a);
        const gdouble b_panel = GSL_MIN (b_p, b);

        printf ("Full panel %u: integrating [%e, %e] using operator %u\n", i, a_panel, b_panel, operator_idx);
        fflush (stdout);

        ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, a_panel, sbilv->j_array_a);
        ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, b_panel, sbilv->j_array_b);

        /* Compute Chebyshev coefficients */
        ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a_panel, b_panel, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);

        /* Convert to Gegenbauer C^(2) basis */
        ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

        printf ("Full panel %u: Computed Gegenbauer with len %u\n", i, sbilv->gegen_coeffs->len);
        fflush (stdout);

        /* Set up RHS */
        g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
        {
          gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
          const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

          rhs_data[0] = 0.0;
          rhs_data[1] = 0.0;
          memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
        }

        /* Process l values in blocks using pre-allocated operator */
        {
          const guint block_size    = 8;
          NcmSBesselOdeOperator *op = g_ptr_array_index (sbilv->operators, operator_idx);
          gint l_start;

          for (l_start = ell_levin_min; l_start <= (gint) ell_levin_max; l_start += block_size)
          {
            const gint l_end = GSL_MIN (l_start + block_size - 1, ell_levin_max);
            gint l;

            ncm_sbessel_ode_operator_reset (op, a_panel, b_panel, ell_levin_min, ell_levin_max);
            ncm_sbessel_ode_operator_solve_endpoints (op, sbilv->rhs, &sbilv->endpoints_result);
            printf ("Full panel %u: Solved ODE for l = %d to %d, solution order %ld\n", i, l_start, l_end,
                    ncm_sbessel_ode_operator_get_n_cols (op));

            for (l = l_start; l <= l_end; l++)
            {
              const gint l_idx        = l - lmin;
              const gint block_idx    = l - l_start;
              const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 0);
              const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 1);

              if (l <= my_lmax)
              {
                const gdouble j_l_a = sbilv->j_array_a[l];
                const gdouble j_l_b = sbilv->j_array_b[l];

                result_data[l_idx] += b_panel * j_l_b * y_prime_b - a_panel * j_l_a * y_prime_a;
              }
            }
          }
        }

        printf ("Full panel %u: Complete\n", i);
        fflush (stdout);
      }
    }

    /* Handle edge panel at the end [last_knot < b, b] if needed */
    if ((last_knot_idx >= 0) && (last_knot_idx < (gint) (sbilv->knots->len - 1)))
    {
      const gdouble a_p  = g_array_index (sbilv->knots, gdouble, last_knot_idx + 1);
      const gdouble b_p  = b;
      const gint my_lmax = GSL_MIN (ell_levin_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p));

      /* Only process if there's an interval beyond the last full knot */
      if (a_p < b_p)
      {
        printf ("Edge panel (end): integrating [%e, %e]\n", a_p, b_p);
        fflush (stdout);

        ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, a_p, sbilv->j_array_a);
        ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, b_p, sbilv->j_array_b);

        /* Compute Chebyshev coefficients */
        ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a_p, b_p, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);

        /* Convert to Gegenbauer C^(2) basis */
        ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

        printf ("Edge panel (end): Computed Gegenbauer with len %u\n", sbilv->gegen_coeffs->len);
        fflush (stdout);

        /* Set up RHS */
        g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
        {
          gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
          const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

          rhs_data[0] = 0.0;
          rhs_data[1] = 0.0;
          memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
        }

        /* Process l values in blocks using temp operator */
        {
          const guint block_size = 8;
          gint l_start;

          for (l_start = ell_levin_min; l_start <= (gint) ell_levin_max; l_start += block_size)
          {
            const gint l_end = GSL_MIN (l_start + block_size - 1, ell_levin_max);
            gint l;

            ncm_sbessel_ode_operator_reset (sbilv->ode_operator_temp_b, a_p, b_p, ell_levin_min, ell_levin_max);
            ncm_sbessel_ode_operator_solve_endpoints (sbilv->ode_operator_temp_b, sbilv->rhs, &sbilv->endpoints_result);
            printf ("Edge panel (end): Solved ODE for l = %d to %d, solution order %ld\n", l_start, l_end,
                    ncm_sbessel_ode_operator_get_n_cols (sbilv->ode_operator_temp_b));

            for (l = l_start; l <= l_end; l++)
            {
              const gint l_idx        = l - lmin;
              const gint block_idx    = l - l_start;
              const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 0);
              const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 1);

              if (l <= my_lmax)
              {
                const gdouble j_l_a = sbilv->j_array_a[l];
                const gdouble j_l_b = sbilv->j_array_b[l];

                result_data[l_idx] += b_p * j_l_b * y_prime_b - a_p * j_l_a * y_prime_a;
              }
            }
          }
        }

        printf ("Edge panel (end): Complete\n");
        fflush (stdout);
      }
    }
  }
  else
  {
    /* Panel mode: divide [a, b] into n_panels in logarithmic space */
    const gdouble ln_a = log (a);
    const gdouble ln_b = log (b);
    const gdouble L    = ln_b - ln_a;
    const gdouble dln  = L / sbilv->n_panels;
    guint panel_idx;

    printf ("Panel mode: dividing range [%e, %e] into %u panels in log space\n", a, b, sbilv->n_panels);
    fflush (stdout);

    for (panel_idx = 0; panel_idx < sbilv->n_panels; panel_idx++)
    {
      const gdouble ln_a_p = ln_a + panel_idx * dln;
      const gdouble ln_b_p = ln_a + (panel_idx + 1) * dln;
      const gdouble a_p    = exp (ln_a_p);
      const gdouble b_p    = exp (ln_b_p);
      const gint my_lmax   = GSL_MIN (ell_levin_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p));

      printf ("Panel %u/%u: integrating [%e, %e]\n", panel_idx + 1, sbilv->n_panels, a_p, b_p);
      fflush (stdout);

      ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, a_p, sbilv->j_array_a);
      ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, b_p, sbilv->j_array_b);

      /* Step 1: Compute Chebyshev coefficients for f(x) for this panel */
      ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a_p, b_p, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);

      /* Step 2: Convert to Gegenbauer C^(2) basis */
      ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

      printf ("Panel %u/%u: Computed Gegenbauer with len %u\n", panel_idx + 1, sbilv->n_panels, sbilv->gegen_coeffs->len);
      fflush (stdout);

      /* Step 3: Set up RHS with homogeneous boundary conditions */
      g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
      {
        gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
        const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

        rhs_data[0] = 0.0; /* BC at x=a_p (t=-1) */
        rhs_data[1] = 0.0; /* BC at x=b_p (t=+1) */

        /* Copy Gegenbauer coefficients to RHS starting from index 2 */
        memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
      }

      /* Step 4-6: Process l values in blocks for better cache locality */
      {
        const guint block_size = 8;
        gint l_start;

        for (l_start = ell_levin_min; l_start <= (gint) ell_levin_max; l_start += block_size)
        {
          const gint l_end = GSL_MIN (l_start + block_size - 1, ell_levin_max);
          gint l;

          ncm_sbessel_ode_operator_reset (sbilv->ode_operator, a_p, b_p, ell_levin_min, ell_levin_max);

          /* Solve for all l values in this block using batched solver */
          printf ("Panel %u/%u: Solving ODE for l = %d to %d\n", panel_idx + 1, sbilv->n_panels, l_start, l_end);
          fflush (stdout);
          ncm_sbessel_ode_operator_solve_endpoints (sbilv->ode_operator, sbilv->rhs, &sbilv->endpoints_result);
          fflush (stderr);
          fflush (stdout);
          printf ("Panel %u/%u: Done solving ODE for l = %d to %d\n", panel_idx + 1, sbilv->n_panels, l_start, l_end);
          fflush (stdout);

          /* Extract derivatives and accumulate integrals */
          for (l = l_start; l <= l_end; l++)
          {
            const gint l_idx        = l - lmin;
            const gint block_idx    = l - l_start;
            const gdouble y_prime_a = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 0);
            const gdouble y_prime_b = g_array_index (sbilv->endpoints_result, gdouble, block_idx * 3 + 1);

            if (l <= my_lmax)
            {
              const gdouble j_l_a = sbilv->j_array_a[l];
              const gdouble j_l_b = sbilv->j_array_b[l];

              /* Accumulate contribution from this panel */
              result_data[l_idx] += b_p * j_l_b * y_prime_b - a_p * j_l_a * y_prime_a;
            }
          }
        }
      }

      printf ("Panel %u/%u: Complete\n", panel_idx + 1, sbilv->n_panels);
      fflush (stdout);
    }
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
  const guint lmin                 = ncm_sbessel_integrator_get_lmin (sbi);
  const guint lmax                 = ncm_sbessel_integrator_get_lmax (sbi);
  const guint n_ell                = lmax - lmin + 1;
  const guint ell_threshold        = _ncm_sbessel_integrator_levin_get_ell_threshold (sbilv, a, b);

  g_assert_cmpuint (ncm_vector_len (result), ==, n_ell);

  /* Ensure resources are allocated */
  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, lmin, lmax);

  /* Use direct cubature integration for high ell values */
  if (lmax >= ell_threshold)
  {
    const guint ell_direct_max = lmax;
    const guint ell_direct_min = GSL_MAX (lmin, ell_threshold);

    if (ell_direct_min <= ell_direct_max)
      _ncm_sbessel_integrator_levin_integrate_direct (sbilv, lmin, ell_direct_min, ell_direct_max, F, a, b, result, user_data);
  }

  /* Use Levin method for low ell values */
  if (lmin < ell_threshold)
  {
    const guint ell_levin_min = lmin;
    const guint ell_levin_max = GSL_MIN (ell_threshold - 1, lmax);

    if (ell_levin_min <= ell_levin_max)
      _ncm_sbessel_integrator_levin_integrate_levin (sbilv, lmin, ell_levin_min, ell_levin_max, F, a, b, result, user_data);
  }
}

/**
 * ncm_sbessel_integrator_levin_new:
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 *
 * Creates a new #NcmSBesselIntegratorLevin.
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorLevin
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_new (guint lmin, guint lmax)
{
  NcmSBesselIntegratorLevin *sbilv = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN,
                                                   "lmin", lmin,
                                                   "lmax", lmax,
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
 * ncm_sbessel_integrator_levin_prepare_for_ell_range:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 *
 * Prepares the integrator for a specific ell range [lmin, lmax].
 * This should be called before integrating for a block of ells.
 *
 * The typical usage pattern is:
 * - Prepare for ell_min=2, ell_max=9
 * - Integrate all kernels for those ells
 * - Prepare for ell_min=10, ell_max=17
 * - Integrate for those ells
 * - And so on...
 *
 * If knots-based paneling is configured (y_knots_min, y_knots_max, n_knots set),
 * this will create or update the ODE operators for the knot panels.
 */
void
ncm_sbessel_integrator_levin_prepare_for_ell_range (NcmSBesselIntegratorLevin *sbilv, gint lmin, gint lmax)
{
  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, lmin, lmax);
  _ncm_sbessel_integrator_levin_prepare_knots_array (sbilv);
  _ncm_sbessel_integrator_levin_prepare_knots_operators (sbilv, lmin, lmax);
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
 * ncm_sbessel_integrator_levin_set_n_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @n_panels: number of panels
 *
 * Sets the number of panels for logarithmic division of the integration range.
 * Use 1 to disable paneling (integrate full range at once). Values greater
 * than 1 will divide the range [a, b] into n_panels in logarithmic space.
 */
void
ncm_sbessel_integrator_levin_set_n_panels (NcmSBesselIntegratorLevin *sbilv, guint n_panels)
{
  g_assert_cmpuint (n_panels, >, 0);
  sbilv->n_panels = n_panels;
}

/**
 * ncm_sbessel_integrator_levin_get_n_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the number of panels.
 *
 * Returns: the number of panels
 */
guint
ncm_sbessel_integrator_levin_get_n_panels (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_panels;
}

/**
 * ncm_sbessel_integrator_levin_set_y_knots_min:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @y_knots_min: minimum knot value
 *
 * Sets the minimum value for the knots grid. Use with y_knots_max and n_knots
 * to define a pre-computed log-spaced knots grid for paneling.
 */
void
ncm_sbessel_integrator_levin_set_y_knots_min (NcmSBesselIntegratorLevin *sbilv, gdouble y_knots_min)
{
  sbilv->y_knots_min = y_knots_min;
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
 * ncm_sbessel_integrator_levin_set_y_knots_max:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @y_knots_max: maximum knot value
 *
 * Sets the maximum value for the knots grid. Use with y_knots_min and n_knots
 * to define a pre-computed log-spaced knots grid for paneling.
 */
void
ncm_sbessel_integrator_levin_set_y_knots_max (NcmSBesselIntegratorLevin *sbilv, gdouble y_knots_max)
{
  sbilv->y_knots_max = y_knots_max;
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
 * ncm_sbessel_integrator_levin_set_n_knots:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @n_knots: number of knots
 *
 * Sets the number of knots in the log-spaced grid. The knots will be
 * equally spaced in log space between y_knots_min and y_knots_max.
 * An ODE operator will be pre-allocated for each panel between consecutive knots.
 */
void
ncm_sbessel_integrator_levin_set_n_knots (NcmSBesselIntegratorLevin *sbilv, guint n_knots)
{
  sbilv->n_knots = n_knots;
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

