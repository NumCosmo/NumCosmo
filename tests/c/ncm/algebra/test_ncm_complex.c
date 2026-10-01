/***************************************************************************
 *            test_ncm_complex.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_complex.c
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

static void
test_ncm_complex_ops (void)
{
  {
    NcmComplex *c1 = ncm_complex_new ();

    ncm_complex_set (c1, 1.0, 2.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 2.0);

    ncm_complex_set_zero (c1);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 0.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 0.0);

    ncm_complex_free (c1);
  }

  {
    NcmComplex c1 = NCM_COMPLEX_INIT (1.0 + 2.0 * I);

    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);
  }

  {
    NcmComplex *c1 = ncm_complex_new ();
    NcmComplex *c2 = ncm_complex_new ();
    NcmComplex *c3 = ncm_complex_new ();

    ncm_complex_set (c1, 1.0, 2.0);
    ncm_complex_set (c2, 3.0, 4.0);
    ncm_complex_set (c3, 5.0, 6.0);

    ncm_complex_mul_real (c1, 3.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 3.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 6.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (3.0, 6.0));

    ncm_complex_res_add_mul (c1, c2, c3);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -6.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 44.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (-6.0, 44.0));

    ncm_complex_res_add_mul_real (c1, c2, 2.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 0.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 52.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (0.0, 52.0));

    ncm_complex_res_mul (c1, c2);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -208.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 156.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (-208.0, 156.0));

    ncm_complex_free (c1);
    ncm_complex_free (c2);
    ncm_complex_free (c3);
  }

  {
    NcmComplex c1;
    complex double z;

    ncm_complex_set_c (&c1, 1.0 + 2.0 * I);
    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);
    g_assert_cmpfloat (ncm_complex_Abs (&c1), ==, hypot (1.0, 2.0));

    z = ncm_complex_c (&c1);
    g_assert_cmpfloat (creal (z), ==, 1.0);
    g_assert_cmpfloat (cimag (z), ==, 2.0);
  }

  {
    NcmComplex c1  = NCM_COMPLEX_INIT (1.0 + 2.0 * I);
    NcmComplex *c2 = ncm_complex_dup (&c1);

    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);
    g_assert_cmpfloat (ncm_complex_Abs (&c1), ==, hypot (1.0, 2.0));

    g_assert_cmpfloat (ncm_complex_Re (c2), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (c2), ==, 2.0);

    ncm_complex_free (c2);
  }

  {
    NcmComplex *c1 = ncm_complex_new ();

    ncm_complex_clear (&c1);
    g_assert_null (c1);
  }
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/complex/ops", &test_ncm_complex_ops);

  g_test_run ();
}

