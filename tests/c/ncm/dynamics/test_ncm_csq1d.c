/***************************************************************************
 *            test_ncm_csq1d.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_csq1d.c
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

#include "test_ncm_csq1d_bessel.h"

static void
test_ncm_csq1d_defaults (void)
{
  /*
   * A subclass implementing only xi, nu and F1 gets nu^2 = nu nu, m = e^xi / nu and
   * F2 = F1' / (2 nu) from NcmCSQ1D; the Bessel system has them in closed form.
   */
  const gdouble a          = 2.0;
  const gdouble k          = 1.3;
  TestCSQ1DBesselMin *bmin = test_csq1d_bessel_min_new (a, k, TRUE);
  NcmCSQ1D *csq1d          = NCM_CSQ1D (bmin);
  gdouble max_F2_err       = 0.0;
  guint i;

  for (i = 0; i < 50; i++)
  {
    const gdouble t = -pow (10.0, -1.0 + 5.0 * i / 49.0);

    ncm_assert_cmpdouble_e (ncm_csq1d_eval_nu2 (csq1d, NULL, t), ==, k * k, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_csq1d_eval_m (csq1d, NULL, t), ==, test_csq1d_bessel_m (a, k, -1.0, t), 1.0e-13, 0.0);

    max_F2_err = GSL_MAX (max_F2_err, fabs (ncm_csq1d_eval_F2 (csq1d, NULL, t) / test_csq1d_bessel_F2 (a, k, -1.0, t) - 1.0));
  }

  /* Measured 3.0e-13 over t in [-1e4, -0.1]. */
  g_assert_cmpfloat (max_F2_err, <, 3.0e-12);

  ncm_csq1d_free (csq1d);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/csq1d/defaults", &test_ncm_csq1d_defaults);

  g_test_run ();
}

