/***************************************************************************
 *            test_ncm_stats_dist2d_spline.c
 *
 *  Tue September 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_stats_dist2d_spline.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/* m2lnp = x^2 + 2 y^2 + x y on [-2, 3] x [-1, 1], which the bicubic spline reproduces */
static gdouble
_test_m2lnp (gdouble x, gdouble y)
{
  return x * x + 2.0 * y * y + x * y;
}

static NcmStatsDist2d *
_test_dist2d_new (void)
{
  const guint nx      = 21;
  const guint ny      = 11;
  NcmVector *xv       = ncm_vector_new (nx);
  NcmVector *yv       = ncm_vector_new (ny);
  NcmMatrix *zm       = ncm_matrix_new (ny, nx);
  NcmSpline2d *m2lnp  = ncm_spline2d_bicubic_notaknot_new ();
  NcmStatsDist2d *sd2 = NULL;
  guint i, j;

  for (i = 0; i < nx; i++)
    ncm_vector_set (xv, i, -2.0 + 5.0 * i / (nx - 1.0));

  for (j = 0; j < ny; j++)
    ncm_vector_set (yv, j, -1.0 + 2.0 * j / (ny - 1.0));

  for (j = 0; j < ny; j++)
    for (i = 0; i < nx; i++)
      ncm_matrix_set (zm, j, i, _test_m2lnp (ncm_vector_get (xv, i), ncm_vector_get (yv, j)));

  ncm_spline2d_set (m2lnp, xv, yv, zm, FALSE);
  sd2 = NCM_STATS_DIST2D (ncm_stats_dist2d_spline_new (m2lnp));

  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_matrix_free (zm);
  ncm_spline2d_free (m2lnp);

  return sd2;
}

static void
test_ncm_stats_dist2d_spline_eval (void)
{
  NcmStatsDist2d *sd2 = _test_dist2d_new ();
  gdouble xi, xf, yi, yf;
  guint i, j;

  ncm_stats_dist2d_prepare (sd2);

  ncm_stats_dist2d_xbounds (sd2, &xi, &xf);
  ncm_stats_dist2d_ybounds (sd2, &yi, &yf);
  g_assert_cmpfloat (xi, ==, -2.0);
  g_assert_cmpfloat (xf, ==, 3.0);
  g_assert_cmpfloat (yi, ==, -1.0);
  g_assert_cmpfloat (yf, ==, 1.0);

  for (i = 0; i <= 10; i++)
  {
    for (j = 0; j <= 10; j++)
    {
      const gdouble x = -2.0 + 0.5 * i + 0.013;
      const gdouble y = -1.0 + 0.2 * j - 0.007 * (j == 10);
      const gdouble m = ncm_stats_dist2d_eval_m2lnp (sd2, x, y);

      ncm_assert_cmpdouble_e (m, ==, _test_m2lnp (x, y), 1.0e-12, 1.0e-12);
      ncm_assert_cmpdouble_e (ncm_stats_dist2d_eval_pdf (sd2, x, y), ==, exp (-0.5 * m), 1.0e-15, 0.0);
    }
  }

  ncm_stats_dist2d_free (sd2);
}

static void
test_ncm_stats_dist2d_spline_unimplemented_subprocess (void)
{
  NcmStatsDist2d *sd2 = _test_dist2d_new ();
  const gchar *which  = g_getenv ("TEST_NCM_STATS_DIST2D_METHOD");

  ncm_stats_dist2d_prepare (sd2);

  if (g_str_equal (which, "cdf"))
    ncm_stats_dist2d_eval_cdf (sd2, 0.0, 0.0);
  else if (g_str_equal (which, "marginal_pdf"))
    ncm_stats_dist2d_eval_marginal_pdf (sd2, 0.0);
  else if (g_str_equal (which, "marginal_cdf"))
    ncm_stats_dist2d_eval_marginal_cdf (sd2, 0.0);
  else if (g_str_equal (which, "marginal_inv_cdf"))
    ncm_stats_dist2d_eval_marginal_inv_cdf (sd2, 0.5);
  else
    ncm_stats_dist2d_eval_inv_cond (sd2, 0.5, 0.0);
}

static void
test_ncm_stats_dist2d_spline_no_m2lnp_subprocess (void)
{
  NcmStatsDist2d *sd2 = g_object_new (NCM_TYPE_STATS_DIST2D_SPLINE, NULL);

  ncm_stats_dist2d_prepare (sd2);
}

static void
test_ncm_stats_dist2d_spline_traps (void)
{
  const gchar *methods[] = {"cdf", "marginal_pdf", "marginal_cdf", "marginal_inv_cdf", "inv_cond"};
  guint k;

  for (k = 0; k < G_N_ELEMENTS (methods); k++)
  {
    gchar *pattern = g_strdup_printf ("*`NcmStatsDist2dSpline' does not implement %s*", methods[k]);

    g_setenv ("TEST_NCM_STATS_DIST2D_METHOD", methods[k], TRUE);
    g_test_trap_subprocess ("/ncm/stats/dist2d_spline/unimplemented/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (pattern);
    g_free (pattern);
  }

  g_unsetenv ("TEST_NCM_STATS_DIST2D_METHOD");

  g_test_trap_subprocess ("/ncm/stats/dist2d_spline/no_m2lnp/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*no m2lnp spline set*");
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/stats/dist2d_spline/eval", &test_ncm_stats_dist2d_spline_eval);
  g_test_add_func ("/ncm/stats/dist2d_spline/traps", &test_ncm_stats_dist2d_spline_traps);
  g_test_add_func ("/ncm/stats/dist2d_spline/unimplemented/subprocess", &test_ncm_stats_dist2d_spline_unimplemented_subprocess);
  g_test_add_func ("/ncm/stats/dist2d_spline/no_m2lnp/subprocess", &test_ncm_stats_dist2d_spline_no_m2lnp_subprocess);

  g_test_run ();
}

