/***************************************************************************
 *            test_nc_wl_ellipticity_series.c
 *
 *  Wed Jul 15 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

/* Checks #NcWLEllipticitySeriesTrace / #NcWLEllipticitySeriesTraceDet's
 * truncated g-Taylor series directly against nc_wl_ellipticity.h's own
 * already-shipped, independently-tested closed forms
 * (apply_shear_inv_*_c, det_jac_*_c): Horner-evaluating chi/jac at a real
 * g inside the truncation radius must match the exact closed form to
 * O(g^(order+1)) -- no finite differences, no dependence on this project's
 * own internal derivation being self-consistent. */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <complex.h>
#include <math.h>
#include <glib.h>
#include <glib-object.h>

/* nc_wl_ellipticity_series.c's own series objects work in the formal
 * variable w; theta -> w=exp(i*theta) is this test's own physics-specific
 * reading of that variable, so it stays local to this file rather than
 * living inside ncm_laurent_series_tps_eval_c() itself. */
static complex double
_tps_eval (NcmLaurentSeriesTPS *tps, gdouble theta, gdouble g)
{
  return ncm_laurent_series_tps_eval_c (tps, cexp (I * theta), g);
}

/* A spread of (rho,theta) shear-map points and small real g values well
 * inside the truncation radius (order below is high enough that the
 * truncation error at these g is far below the tolerance). */
static const gdouble test_rhos[]   = { 0.05, 0.25, 0.60, 0.85 };
static const gdouble test_thetas[] = { 0.0, 0.7, 2.1, -1.3 };
static const gdouble test_gs[]     = { 0.0, 0.02, -0.03, 0.05 };

#define TEST_ORDER 10
#define TEST_TOL   1.0e-9

void test_nc_wl_ellipticity_series_trace_matches_closed_form (void);
void test_nc_wl_ellipticity_series_trace_det_matches_closed_form (void);
void test_nc_wl_ellipticity_series_abs_sq_matches_chi_squared (void);
void test_nc_wl_ellipticity_series_trace_matches_python_golden_value (void);
void test_nc_wl_ellipticity_series_eval_is_reusable (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/wl_ellipticity_series/trace/matches_closed_form",
                   &test_nc_wl_ellipticity_series_trace_matches_closed_form);
  g_test_add_func ("/nc/wl_ellipticity_series/trace_det/matches_closed_form",
                   &test_nc_wl_ellipticity_series_trace_det_matches_closed_form);
  g_test_add_func ("/nc/wl_ellipticity_series/abs_sq_matches_chi_squared",
                   &test_nc_wl_ellipticity_series_abs_sq_matches_chi_squared);
  g_test_add_func ("/nc/wl_ellipticity_series/trace/matches_python_golden_value",
                   &test_nc_wl_ellipticity_series_trace_matches_python_golden_value);
  g_test_add_func ("/nc/wl_ellipticity_series/eval_is_reusable",
                   &test_nc_wl_ellipticity_series_eval_is_reusable);

  return g_test_run ();
}

void
test_nc_wl_ellipticity_series_trace_matches_closed_form (void)
{
  NcWLEllipticitySeriesTrace *ser = nc_wl_ellipticity_series_trace_new (TEST_ORDER);
  guint ir, it, ig;

  g_assert_cmpuint (nc_wl_ellipticity_series_trace_get_order (ser), ==, TEST_ORDER);

  for (ir = 0; ir < G_N_ELEMENTS (test_rhos); ir++)
  {
    const gdouble rho = test_rhos[ir];

    nc_wl_ellipticity_series_trace_eval (ser, rho);

    for (it = 0; it < G_N_ELEMENTS (test_thetas); it++)
    {
      const gdouble theta       = test_thetas[it];
      const complex double chiL = rho * cexp (I * theta);

      for (ig = 0; ig < G_N_ELEMENTS (test_gs); ig++)
      {
        const gdouble g = test_gs[ig];

        const complex double chi_series = _tps_eval (nc_wl_ellipticity_series_trace_get_chi (ser), theta, g);
        const gdouble jac_series        = creal (_tps_eval (nc_wl_ellipticity_series_trace_get_jac (ser), theta, g));

        const complex double chi_exact = nc_wl_ellipticity_apply_shear_inv_trace_c (g, chiL);
        const gdouble jac_exact        = nc_wl_ellipticity_det_jac_trace_c (g, chiL);

        g_assert_cmpfloat (cabs (chi_series - chi_exact), <, TEST_TOL);
        g_assert_cmpfloat (fabs (jac_series - jac_exact), <, TEST_TOL);
      }
    }
  }

  nc_wl_ellipticity_series_trace_free (ser);
}

void
test_nc_wl_ellipticity_series_trace_det_matches_closed_form (void)
{
  NcWLEllipticitySeriesTraceDet *ser = nc_wl_ellipticity_series_trace_det_new (TEST_ORDER);
  guint ir, it, ig;

  g_assert_cmpuint (nc_wl_ellipticity_series_trace_det_get_order (ser), ==, TEST_ORDER);

  for (ir = 0; ir < G_N_ELEMENTS (test_rhos); ir++)
  {
    const gdouble rho = test_rhos[ir];

    nc_wl_ellipticity_series_trace_det_eval (ser, rho);

    for (it = 0; it < G_N_ELEMENTS (test_thetas); it++)
    {
      const gdouble theta       = test_thetas[it];
      const complex double chiL = rho * cexp (I * theta);

      for (ig = 0; ig < G_N_ELEMENTS (test_gs); ig++)
      {
        const gdouble g = test_gs[ig];

        const complex double chi_series = _tps_eval (nc_wl_ellipticity_series_trace_det_get_chi (ser), theta, g);
        const gdouble jac_series        = creal (_tps_eval (nc_wl_ellipticity_series_trace_det_get_jac (ser), theta, g));

        const complex double chi_exact = nc_wl_ellipticity_apply_shear_inv_trace_det_c (g, chiL);
        const gdouble jac_exact        = nc_wl_ellipticity_det_jac_trace_det_c (g, chiL);

        g_assert_cmpfloat (cabs (chi_series - chi_exact), <, TEST_TOL);
        g_assert_cmpfloat (fabs (jac_series - jac_exact), <, TEST_TOL);
      }
    }
  }

  nc_wl_ellipticity_series_trace_det_free (ser);
}

/* abs_sq must be exactly |chi|^2, evaluated the same way as chi itself. */
void
test_nc_wl_ellipticity_series_abs_sq_matches_chi_squared (void)
{
  NcWLEllipticitySeriesTrace *trace   = nc_wl_ellipticity_series_trace_new (TEST_ORDER);
  NcWLEllipticitySeriesTraceDet *tdet = nc_wl_ellipticity_series_trace_det_new (TEST_ORDER);
  const gdouble rho                   = 0.4;
  const gdouble theta                 = 0.9;
  const gdouble g                     = 0.03;
  complex double chi, abs_sq;

  nc_wl_ellipticity_series_trace_eval (trace, rho);
  chi    = _tps_eval (nc_wl_ellipticity_series_trace_get_chi (trace), theta, g);
  abs_sq = _tps_eval (nc_wl_ellipticity_series_trace_get_abs_sq (trace), theta, g);
  g_assert_cmpfloat (cabs (abs_sq - chi * conj (chi)), <, TEST_TOL);

  nc_wl_ellipticity_series_trace_det_eval (tdet, rho);
  chi    = _tps_eval (nc_wl_ellipticity_series_trace_det_get_chi (tdet), theta, g);
  abs_sq = _tps_eval (nc_wl_ellipticity_series_trace_det_get_abs_sq (tdet), theta, g);
  g_assert_cmpfloat (cabs (abs_sq - chi * conj (chi)), <, TEST_TOL);

  nc_wl_ellipticity_series_trace_free (trace);
  nc_wl_ellipticity_series_trace_det_free (tdet);
}

/* Cross-language regression: the exact same golden value
 * test_chi_taylor_matches_python_reference() (test_ncm_laurent_series.c)
 * checks against its own hand-rolled recursion -- checked here directly
 * against the real, shipped object instead. rho=0.25, theta=0.7, g=0.09,
 * N=12: chi_I ~ 0.0131042280+0.1640678007j. */
void
test_nc_wl_ellipticity_series_trace_matches_python_golden_value (void)
{
  NcWLEllipticitySeriesTrace *ser = nc_wl_ellipticity_series_trace_new (12);
  complex double chi;

  nc_wl_ellipticity_series_trace_eval (ser, 0.25);
  chi = _tps_eval (nc_wl_ellipticity_series_trace_get_chi (ser), 0.7, 0.09);

  g_assert_cmpfloat (fabs (creal (chi) - 0.0131042280), <, 1.0e-6);
  g_assert_cmpfloat (fabs (cimag (chi) - 0.1640678007), <, 1.0e-6);

  nc_wl_ellipticity_series_trace_free (ser);
}

/* eval() must be a pure refill: calling it again at a different rho must
 * not leave any trace of the previous call's content (private conv/fold
 * scratch reused correctly, no stale-state leakage). */
void
test_nc_wl_ellipticity_series_eval_is_reusable (void)
{
  NcWLEllipticitySeriesTrace *trace   = nc_wl_ellipticity_series_trace_new (TEST_ORDER);
  NcWLEllipticitySeriesTraceDet *tdet = nc_wl_ellipticity_series_trace_det_new (TEST_ORDER);
  const gdouble theta                 = 0.4;
  const gdouble g                     = 0.02;
  guint pass;

  for (pass = 0; pass < 3; pass++)
  {
    const gdouble rho         = 0.1 + 0.2 * pass;
    const complex double chiL = rho * cexp (I * theta);

    nc_wl_ellipticity_series_trace_eval (trace, rho);
    g_assert_cmpfloat (cabs (_tps_eval (nc_wl_ellipticity_series_trace_get_chi (trace), theta, g) -
                             nc_wl_ellipticity_apply_shear_inv_trace_c (g, chiL)), <, TEST_TOL);

    nc_wl_ellipticity_series_trace_det_eval (tdet, rho);
    g_assert_cmpfloat (cabs (_tps_eval (nc_wl_ellipticity_series_trace_det_get_chi (tdet), theta, g) -
                             nc_wl_ellipticity_apply_shear_inv_trace_det_c (g, chiL)), <, TEST_TOL);
  }

  nc_wl_ellipticity_series_trace_free (trace);
  nc_wl_ellipticity_series_trace_det_free (tdet);
}

