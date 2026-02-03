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
 * This class implements integration of functions multiplied by spherical
 * Bessel functions using a Levin-type method. The interval is divided into
 * panels, and in each panel the function is discretized on Clenshaw-Curtis
 * knots of increasing order. The method solves the differential equation
 * $x^2 y''(x) + 2x y'(x) + (x^2 - \ell(\ell+1)) y(x) = f(x)$
 * with boundary conditions $y(x_n) = 0 = y(x_{n+1})$.
 * The integral in each panel is then given by
 * $I_n = j_\ell(x_{n+1}) y'(x_{n+1}) - j_\ell(x_n) y'(x_n)$.
 * Note: The formulation uses $x^2 D^2 + 2x D + x^2 I$ as the ell-independent
 * part of the matrix, with only the diagonal term $-\ell(\ell+1)$ varying.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_integrator_levin.h"
#include "math/ncm_sbessel_ode_solver.h"
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

/* Matrix layout configuration */
#define NCM_LEVIN_USE_COLMAJOR 1

#if NCM_LEVIN_USE_COLMAJOR
#define NCM_LEVIN_CBLAS_LAYOUT CblasColMajor
#define NCM_LEVIN_MAT(m, i, j, N) ((m)[(i) + (j) * (N)])
#else
#define NCM_LEVIN_CBLAS_LAYOUT CblasRowMajor
#define NCM_LEVIN_MAT(m, i, j, N) ((m)[(i) * (N) + (j)])
#endif

struct _NcmSBesselIntegratorLevin
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  guint n_panels;
  guint min_order;
  guint max_order;
  gdouble reltol;
  NcmLapackWS *ws;
  gdouble *A;     /* Tridiagonal system matrix */
  gdouble *B;     /* Right-hand side vector */
  gdouble *y_sol; /* Solution vector */
  gdouble *x_cc;  /* Clenshaw-Curtis knots */
  gint *ipiv;     /* Pivot indices for LAPACK */
  guint alloc_order;
  NcmSBesselOdeSolver *ode_solver; /* Spherical Bessel ODE solver */
};

enum
{
  PROP_0,
  PROP_N_PANELS,
  PROP_MIN_ORDER,
  PROP_MAX_ORDER,
  PROP_RELTOL,
};

G_DEFINE_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_levin_init (NcmSBesselIntegratorLevin *sbilv)
{
  sbilv->n_panels    = 10;
  sbilv->min_order   = 4;
  sbilv->max_order   = 32;
  sbilv->reltol      = 1.0e-7;
  sbilv->ws          = ncm_lapack_ws_new ();
  sbilv->A           = NULL;
  sbilv->B           = NULL;
  sbilv->y_sol       = NULL;
  sbilv->x_cc        = NULL;
  sbilv->ipiv        = NULL;
  sbilv->alloc_order = 0;
  sbilv->ode_solver  = ncm_sbessel_ode_solver_new (0, -1.0, 1.0);
}

static void
_ncm_sbessel_integrator_levin_dispose (GObject *object)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  ncm_lapack_ws_clear (&sbilv->ws);

  g_clear_pointer (&sbilv->A, g_free);
  g_clear_pointer (&sbilv->B, g_free);
  g_clear_pointer (&sbilv->y_sol, g_free);
  g_clear_pointer (&sbilv->x_cc, g_free);
  g_clear_pointer (&sbilv->ipiv, g_free);

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
    case PROP_N_PANELS:
      ncm_sbessel_integrator_levin_set_n_panels (sbilv, g_value_get_uint (value));
      break;
    case PROP_MIN_ORDER:
      ncm_sbessel_integrator_levin_set_min_order (sbilv, g_value_get_uint (value));
      break;
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
    case PROP_N_PANELS:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_n_panels (sbilv));
      break;
    case PROP_MIN_ORDER:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_min_order (sbilv));
      break;
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
   * NcmSBesselIntegratorLevin:n-panels:
   *
   * Number of panels to divide the integration interval.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_PANELS,
                                   g_param_spec_uint ("n-panels",
                                                      NULL,
                                                      "Number of panels",
                                                      1, G_MAXUINT, 20,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:min-order:
   *
   * Minimum order of Clenshaw-Curtis quadrature.
   */
  g_object_class_install_property (object_class,
                                   PROP_MIN_ORDER,
                                   g_param_spec_uint ("min-order",
                                                      NULL,
                                                      "Minimum Clenshaw-Curtis order",
                                                      2, G_MAXUINT, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:max-order:
   *
   * Maximum order of Clenshaw-Curtis quadrature.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum Clenshaw-Curtis order",
                                                      2, G_MAXUINT, 512,
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

/*
 * _ncm_sbessel_integrator_levin_alloc_workspace:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @order: order of the Clenshaw-Curtis quadrature
 *
 * Allocates workspace arrays for the given order.
 */
static void
_ncm_sbessel_integrator_levin_alloc_workspace (NcmSBesselIntegratorLevin *sbilv, guint order)
{
  if (sbilv->alloc_order >= order)
    return;

  g_clear_pointer (&sbilv->A, g_free);
  g_clear_pointer (&sbilv->B, g_free);
  g_clear_pointer (&sbilv->y_sol, g_free);
  g_clear_pointer (&sbilv->x_cc, g_free);
  g_clear_pointer (&sbilv->ipiv, g_free);

  sbilv->A           = g_new (gdouble, order * order);
  sbilv->B           = g_new (gdouble, order);
  sbilv->y_sol       = g_new (gdouble, order);
  sbilv->x_cc        = g_new (gdouble, order);
  sbilv->ipiv        = g_new (gint, order);
  sbilv->alloc_order = order;
}

/*
 * _ncm_sbessel_integrator_levin_chebyshev_nodes:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @a: left endpoint
 * @b: right endpoint
 * @order: number of nodes
 * @x: (out): array to store nodes
 *
 * Computes Chebyshev nodes on the interval [a, b].
 * The nodes are given by x_k = (a+b)/2 + (b-a)/2 * cos(k*pi/(N-1))
 * for k = 0, ..., N-1.
 */
static void
_ncm_sbessel_integrator_levin_chebyshev_nodes (NcmSBesselIntegratorLevin *sbilv,
                                               gdouble a, gdouble b,
                                               guint N, gdouble *x)
{
  const gdouble mid    = 0.5 * (a + b);
  const gdouble half_h = 0.5 * (b - a);
  guint k;

  for (k = 0; k < N; k++)
  {
    gdouble theta = M_PI * k / (N - 1);

    x[k] = mid - half_h * cos (theta);
  }
}

/*
 * _ncm_sbessel_integrator_levin_chebyshev_diff_matrix:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @N: number of nodes
 * @a: left endpoint
 * @b: right endpoint
 * @D: (out): differentiation matrix (N x N)
 *
 * Computes the Chebyshev differentiation matrix for the interval [a, b].
 */
static void
_ncm_sbessel_integrator_levin_chebyshev_diff_matrix (NcmSBesselIntegratorLevin *sbilv,
                                                     guint N, gdouble a, gdouble b,
                                                     gdouble *D)
{
  guint i, j;
  gdouble *x = g_new (gdouble, N);
  gdouble *c = g_new (gdouble, N);

  /* Compute Chebyshev nodes in [-1, 1] */
  for (i = 0; i < N; i++)
  {
    x[i] = -cos (M_PI * i / (N - 1));
    c[i] = 1.0;
  }

  c[0]     = 2.0;
  c[N - 1] = 2.0;

  /* Build differentiation matrix */
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      if (i != j)
      {
        gint sign = ((i + j) % 2 == 0) ? 1 : -1;

        NCM_LEVIN_MAT (D, i, j, N) = (c[i] / c[j]) * sign / (x[i] - x[j]);
      }
      else
      {
        NCM_LEVIN_MAT (D, i, i, N) = 0.0;
      }
    }
  }

  /* Diagonal entries: negative sum of row */
  for (i = 0; i < N; i++)
  {
    gdouble row_sum = 0.0;

    for (j = 0; j < N; j++)
    {
      if (i != j)
        row_sum += NCM_LEVIN_MAT (D, i, j, N);
    }

    NCM_LEVIN_MAT (D, i, i, N) = -row_sum;
  }

  /* Scale for interval [a, b] */
  gdouble scale = 2.0 / (b - a);

  for (i = 0; i < N * N; i++)
    D[i] *= scale;

  g_free (x);
  g_free (c);
}

/*
 * _ncm_sbessel_integrator_levin_solve_ode:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @F: function to integrate
 * @user_data: user data for @F
 * @a: left endpoint of panel
 * @b: right endpoint of panel
 * @ell: multipole order
 * @N: number of Chebyshev nodes
 * @u: (out): solution vector
 * @up_a: (out): derivative at left endpoint (first node)
 * @up_b: (out): derivative at right endpoint (last node)
 *
 * Solves the Levin ODE on a panel using Chebyshev differentiation matrices.
 * The ODE is: x^2*y''(x) + 2x*y'(x) + (x^2 - ell(ell+1))*y(x) = f(x)
 * with boundary conditions y(a) = 0 = y(b).
 *
 * Returns: %TRUE on success, %FALSE on failure
 */
static gboolean
_ncm_sbessel_integrator_levin_solve_ode (NcmSBesselIntegratorLevin *sbilv,
                                         NcmSBesselIntegratorF F,
                                         gpointer user_data,
                                         gdouble a, gdouble b,
                                         gint ell,
                                         guint N,
                                         gdouble iu_a,
                                         gdouble iup_a,
                                         gdouble *u,
                                         gdouble *up_a,
                                         gdouble *up_b)
{
  guint i, j;
  gint info;
  const gdouble ell_factor = ell * (ell + 1.0);
  gdouble *D               = NULL;
  gdouble *D2              = NULL;
  gdouble *up              = NULL;

  _ncm_sbessel_integrator_levin_alloc_workspace (sbilv, N);
  _ncm_sbessel_integrator_levin_chebyshev_nodes (sbilv, a, b, N, sbilv->x_cc);

  D  = g_new (gdouble, N * N);
  D2 = g_new (gdouble, N * N);
  up = g_new (gdouble, N);

  /* Compute Chebyshev differentiation matrix */
  _ncm_sbessel_integrator_levin_chebyshev_diff_matrix (sbilv, N, a, b, D);

  /* Compute D2 = D @ D (matrix multiplication) using BLAS */
  cblas_dgemm (NCM_LEVIN_CBLAS_LAYOUT, CblasNoTrans, CblasNoTrans, N, N, N,
               1.0, D, N, D, N, 0.0, D2, N);

  /* Build the system matrix A = x^2*D2 + 2x*D + (x^2 - ell(ell+1))*I */
  /* and right-hand side rhs = f(x) */
  for (i = 0; i < N; i++)
  {
    const gdouble x  = sbilv->x_cc[i];
    const gdouble x2 = x * x;

    sbilv->B[i] = F (user_data, x);

    for (j = 0; j < N; j++)
    {
      NCM_LEVIN_MAT (sbilv->A, i, j, N) = x2 * NCM_LEVIN_MAT (D2, i, j, N) + 2.0 * x * NCM_LEVIN_MAT (D, i, j, N);

      if (i == j)
        NCM_LEVIN_MAT (sbilv->A, i, j, N) += x2 - ell_factor;
    }
  }

  /* Apply boundary conditions: u[0] = 0 and u[N-1] = 0 */
  for (j = 0; j < N; j++)
  {
    NCM_LEVIN_MAT (sbilv->A, 0, j, N)     = 0.0;
    NCM_LEVIN_MAT (sbilv->A, N - 1, j, N) = NCM_LEVIN_MAT (D, 0, j, N);
  }

  NCM_LEVIN_MAT (sbilv->A, 0, 0, N) = 1.0;

  sbilv->B[0]     = iu_a;
  sbilv->B[N - 1] = iup_a;

  /* Solve the linear system A * u = rhs using LAPACK */
  /* Use DGESV (general dense solver) for non-symmetric matrix */
  info = ncm_lapack_dgesv (N, 1, sbilv->A, N, sbilv->ipiv, sbilv->B, N);

  if (info != 0)
  {
    g_free (D);
    g_free (D2);
    g_free (up);

    return FALSE;
  }

  /* Copy solution */
  for (i = 0; i < N; i++)
  {
    u[i] = sbilv->B[i];
  }

  /* Compute derivative up = D @ u using BLAS */
  cblas_dgemv (NCM_LEVIN_CBLAS_LAYOUT, CblasNoTrans, N, N, 1.0, D, N, u, 1, 0.0, up, 1);

  *up_a = up[0];
  *up_b = up[N - 1];

  g_free (D);
  g_free (D2);
  g_free (up);

  return TRUE;
}

/*
 * _ncm_sbessel_integrator_levin_panel_integral:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @F: function to integrate
 * @user_data: user data for @F
 * @a: left endpoint of panel
 * @b: right endpoint of panel
 * @ell: multipole order
 * @integral: (out): computed integral value
 *
 * Computes the integral over a single panel using the Levin method
 * with adaptive order refinement.
 *
 * Returns: whether the computation was successful
 */
static gboolean
_ncm_sbessel_integrator_levin_panel_integral (NcmSBesselIntegratorLevin *sbilv,
                                              NcmSBesselIntegratorF F,
                                              gpointer user_data,
                                              gdouble a, gdouble b,
                                              gint ell, gdouble *integral)
{
  {
    ncm_sbessel_ode_solver_set_interval (sbilv->ode_solver, a, b);
    ncm_sbessel_ode_solver_set_l (sbilv->ode_solver, ell);
    *integral = ncm_sbessel_ode_solver_integrate (sbilv->ode_solver, F, 4096 * 16, user_data);

    return TRUE;
  }


  {
    guint order    = sbilv->min_order * 4;
    gdouble result = 0.0;
    gdouble prev_result;
    gboolean converged = FALSE;

    while (order <= sbilv->max_order && !converged)
    {
      gdouble *u = g_new (gdouble, order);
      gdouble up_a, up_b;
      gboolean success;

      success = _ncm_sbessel_integrator_levin_solve_ode (sbilv, F, user_data, a, b, ell, order, 0.0, 0.0, u, &up_a, &up_b);

      if (!success)
      {
        g_free (u);
        order *= 2;
        continue;
      }

      /* Compute integral: I = j_l(b) * u'(b) * b^2 - j_l(a) * u'(a) * a^2 */
      const gdouble jl_a   = gsl_sf_bessel_jl (ell, a);
      const gdouble jlp1_a = gsl_sf_bessel_jl (ell + 1, a);
      const gdouble jl_b   = gsl_sf_bessel_jl (ell, b);
      const gdouble jlp1_b = gsl_sf_bessel_jl (ell + 1, b);
      const gdouble djl_a  = (ell * jl_a) / a - jlp1_a;
      const gdouble djl_b  = (ell * jl_b) / b - jlp1_b;
      const gdouble u_a    = u[0];
      const gdouble u_b    = u[order - 1];
      const gdouble w_a    = (jl_a * up_a - djl_a * u_a) * a * a;
      const gdouble w_b    = (jl_b * up_b - djl_b * u_b) * b * b;

      prev_result = result;
      result      = w_b - w_a;

      g_free (u);

      /* Check convergence */
      if (order > sbilv->min_order)
      {
        const gdouble abs_diff = fabs (result - prev_result);
        const gdouble abs_val  = fabs (result);

        if (abs_diff < sbilv->reltol * abs_val)
          converged = TRUE;
      }

      order *= 2;
    }

    *integral = result;

    return converged;
  }
}

/*
 * _ncm_sbessel_integrator_levin_integrate_panel_recursive:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @F: function to integrate
 * @user_data: user data for @F
 * @a: left endpoint of panel
 * @b: right endpoint of panel
 * @ell: multipole order
 * @depth: current recursion depth
 *
 * Recursively integrates a panel, splitting it in half if convergence fails.
 *
 * Returns: the integral value
 */
static gdouble
_ncm_sbessel_integrator_levin_integrate_panel_recursive (NcmSBesselIntegratorLevin *sbilv,
                                                         NcmSBesselIntegratorF F,
                                                         gpointer user_data,
                                                         gdouble a, gdouble b,
                                                         gint ell,
                                                         guint depth)
{
  gdouble result = 0.0;
  gboolean success;

  /* Maximum recursion depth to avoid infinite loops */
  const guint max_depth = 10000000;

  success = _ncm_sbessel_integrator_levin_panel_integral (sbilv, F, user_data, a, b, ell, &result);

  if (!success && (depth < max_depth))
  {
    const gdouble L     = b - a;
    const gint n_pi     = floor (L / M_PI);
    const gdouble width = n_pi > 1 ? ((1) * M_PI) : (L / 2);
    gdouble left_result, right_result;

    left_result  = _ncm_sbessel_integrator_levin_integrate_panel_recursive (sbilv, F, user_data, a, a + width, ell, depth + 1);
    right_result = _ncm_sbessel_integrator_levin_integrate_panel_recursive (sbilv, F, user_data, a + width, b, ell, depth + 1);

    result = left_result + right_result;
  }

  return result;
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
  const gdouble L                  = b - a;
  const gdouble lll                = 230;
  const guint n_pi                 = floor (L / (lll * M_PI));
  gdouble result                   = 0.0;
  gdouble x0                       = a;
  gdouble x1                       = a;
  gdouble r0;
  guint i;

  for (i = 0; i < n_pi; i++)
  {
    x1      = x0 + lll * M_PI;
    r0      = _ncm_sbessel_integrator_levin_integrate_panel_recursive (sbilv, F, user_data, x0, x1, ell, 0);
    result += r0;

    printf ("Intermediate result for ell %d (% 22.15g, % 22.15g): % 22.15g | % 22.15g\n", ell, x0, x1, result, r0);
    x0 = x1;
  }

  /* Left over panel */
  {
    if (x0 < b)
    {
      r0      = _ncm_sbessel_integrator_levin_integrate_panel_recursive (sbilv, F, user_data, x0, b, ell, 0);
      result += r0;

      printf ("Intermediate result for ell %d (% 22.15g, % 22.15g): % 22.15g | % 22.15g\n", ell, x0, b, result, r0);
    }
  }

  printf ("Result for ell %d (% 22.15g, % 22.15g): % 22.15g\n", ell, a, b, result);


  return result;
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
 * ncm_sbessel_integrator_levin_set_n_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @n_panels: number of panels
 *
 * Sets the number of panels to divide the integration interval.
 */
void
ncm_sbessel_integrator_levin_set_n_panels (NcmSBesselIntegratorLevin *sbilv, guint n_panels)
{
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
 * ncm_sbessel_integrator_levin_set_min_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @min_order: minimum order
 *
 * Sets the minimum order of Clenshaw-Curtis quadrature.
 */
void
ncm_sbessel_integrator_levin_set_min_order (NcmSBesselIntegratorLevin *sbilv, guint min_order)
{
  sbilv->min_order = min_order;
}

/**
 * ncm_sbessel_integrator_levin_get_min_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the minimum order.
 *
 * Returns: the minimum order
 */
guint
ncm_sbessel_integrator_levin_get_min_order (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->min_order;
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

