/***************************************************************************
 *            test_ncm_util.c
 *
 *  Tue July 30 19:54:22 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiooliveiraCode@proton.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira <caiooliveiraCode@proton.me>
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

void test_ncm_util_projected_radius (void);
void test_ncm_util_complex (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/util/projected_radius", test_ncm_util_projected_radius);
  g_test_add_func ("/ncm/util/complex", test_ncm_util_complex);

  g_test_run ();
}

void
test_ncm_util_projected_radius (void)
{
  NcmRNG *rng    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  gdouble ntests = 1000;
  guint i;

  for (i = 0; i < ntests; i++)
  {
    gdouble d     = ncm_rng_uniform_gen (rng, 0.0, 100.0);
    gdouble theta = ncm_rng_uniform_gen (rng, 0.0, M_PI);

    g_assert_cmpfloat (ncm_util_projected_radius (theta, d), >=, 0);
  }
}

void
test_ncm_util_complex (void)
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

    ncm_complex_res_add_mul (c1, c2, c3);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -6.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 44.0);

    ncm_complex_res_add_mul_real (c1, c2, 2.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 0.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 52.0);

    ncm_complex_free (c1);
    ncm_complex_free (c2);
    ncm_complex_free (c3);
  }
}

