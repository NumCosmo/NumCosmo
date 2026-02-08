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
 * functions using a Levin-type method.
 *
 * The method solves the differential equation $x^2 y''(x) + 2x y'(x) + (x^2 -
 * \ell(\ell+1)) y(x) = f(x)$ with boundary conditions $y(x_n) = 0 = y(x_{n+1})$. The
 * integral is then given by $I_n = j_\ell(x_{n+1}) y'(x_{n+1}) - j_\ell(x_n) y'(x_n)$.
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
#include "ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmSBesselIntegratorLevin
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  guint max_order;
  gdouble reltol;
  NcmSBesselOdeSolver *ode_solver;
};

enum
{
  PROP_0,
  PROP_MAX_ORDER,
  PROP_RELTOL,
};

G_DEFINE_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_levin_init (NcmSBesselIntegratorLevin *sbilv)
{
  sbilv->max_order  = 0;
  sbilv->reltol     = 0.0;
  sbilv->ode_solver = ncm_sbessel_ode_solver_new (0, -1.0, 1.0);
}

static void
_ncm_sbessel_integrator_levin_dispose (GObject *object)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  ncm_sbessel_ode_solver_clear (&sbilv->ode_solver);

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

  parent_class->prepare       = &_ncm_sbessel_integrator_levin_prepare;
  parent_class->integrate_ell = &_ncm_sbessel_integrator_levin_integrate_ell;
}

static void
_ncm_sbessel_integrator_levin_prepare (NcmSBesselIntegrator *sbi)
{
  /* Preparation could involve pre-computing panel boundaries */
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
  const guint N                    = sbilv->max_order;
  NcmSpectral *spectral            = ncm_sbessel_ode_solver_peek_spectral (solver);
  NcmVector *cheb_coeffs           = ncm_vector_new (N);
  NcmVector *gegen_coeffs          = ncm_vector_new (N);
  NcmVector *rhs                   = ncm_vector_new (N + 2);
  gdouble y_prime_a, y_prime_b, error;

  /* Set the interval and l value for this integration */
  ncm_sbessel_ode_solver_set_l (solver, ell);
  ncm_sbessel_ode_solver_set_interval (solver, a, b);

  ncm_spectral_compute_chebyshev_coeffs (spectral, F, a, b, cheb_coeffs, user_data);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (cheb_coeffs, gegen_coeffs);

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (gegen_coeffs), ==, 1);
  {
    gdouble *rhs_data                = ncm_vector_data (rhs);
    const gdouble *gegen_coeffs_data = ncm_vector_data (gegen_coeffs);

    rhs_data[0] = 0.0; /* BC at x=a (t=-1) */
    rhs_data[1] = 0.0; /* BC at x=b (t=+1) */

    /* Copy Gegenbauer coefficients to RHS starting from index 2 */
    memcpy (&rhs_data[2], gegen_coeffs_data, N * sizeof (gdouble));
  }

  ncm_sbessel_ode_solver_solve_endpoints (solver, rhs, &y_prime_a, &y_prime_b, &error);

  {
    const gdouble j_l_a    = gsl_sf_bessel_jl (ell, a);
    const gdouble j_l_b    = gsl_sf_bessel_jl (ell, b);
    const gdouble integral = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a;

    /* Clean up */
    ncm_vector_free (cheb_coeffs);
    ncm_vector_free (gegen_coeffs);
    ncm_vector_free (rhs);

    return integral;
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

