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
};

enum
{
  PROP_0,
  PROP_MAX_ORDER,
  PROP_RELTOL,
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

static void
_ncm_sbessel_integrator_levin_prepare (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);
  const guint lmin                 = ncm_sbessel_integrator_get_lmin (sbi);
  const guint lmax                 = ncm_sbessel_integrator_get_lmax (sbi);

  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, lmin, lmax);
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
  const gint my_lmax          = GSL_MIN (ell_levin_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b));
  gdouble *result_data        = ncm_vector_data (result);

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);

  ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, a, sbilv->j_array_a);
  ncm_sf_sbessel_array_eval (sbilv->sba, my_lmax, b, sbilv->j_array_b);

  /* Initialize Levin results to zero */
  memset (&result_data[ell_levin_min - lmin], 0, sizeof (gdouble) * (ell_levin_max - ell_levin_min + 1));

  /* Step 1: Compute Chebyshev coefficients for f(x) - done once */
  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, F, a, b, 3, 1.0e-11, &sbilv->cheb_coeffs, user_data);

  /* Step 2: Convert to Gegenbauer C^(2) basis - done once */
  ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);

  /* Step 3: Set up RHS with homogeneous boundary conditions - done once */
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
      ncm_sbessel_ode_operator_solve_endpoints (sbilv->ode_operator, sbilv->rhs, &sbilv->endpoints_result);

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

