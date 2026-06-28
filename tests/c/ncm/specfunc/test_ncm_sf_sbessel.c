/***************************************************************************
 *            test_ncm_sf_sbessel.c
 *
 *  Tue July 03 13:35:29 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#include <gsl/gsl_sf_bessel.h>

typedef struct _TestNcmSFSBessel
{
  guint ntests;
} TestNcmSFSBessel;

void test_ncm_sf_sbessel_new (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_free (TestNcmSFSBessel *test, gconstpointer pdata);

void test_ncm_sf_sbessel_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_taylor_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_spline_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata);

void test_ncm_sf_sbessel_traps (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_invalid_st (TestNcmSFSBessel *test, gconstpointer pdata);

#define NTOT 100
#define XMAX 2.0
#define L 80

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/sf/sbessel/cmp/gsl", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_cmp_gsl,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/taylor/cmp/gsl", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_taylor_cmp_gsl,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/spline/cmp/gsl", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_spline_cmp_gsl,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/traps", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_traps,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/invalid/st/subprocess", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_invalid_st,
              &test_ncm_sf_sbessel_free);

  g_test_run ();
}

void
test_ncm_sf_sbessel_new (TestNcmSFSBessel *test, gconstpointer pdata)
{
}

void
test_ncm_sf_sbessel_free (TestNcmSFSBessel *test, gconstpointer pdata)
{
  ncm_mpsf_sbessel_free_cache ();
}

void
test_ncm_sf_sbessel_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata)
{
  guint i, j;

  for (j = 0; j <= L; j++)
  {
    if (j == 0)
      ncm_assert_cmpdouble_e (ncm_sf_sbessel (j, 0.0), ==, 1.0, 1.0e-7, 0.0);
    else
      ncm_assert_cmpdouble_e (ncm_sf_sbessel (j, 0.0), ==, 0.0, 1.0e-7, 0.0);

    for (i = 0; i < NTOT; i++)
    {
      const gdouble x      = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * i);
      const gdouble ncm_jl = ncm_sf_sbessel (j, x);
      const gdouble gsl_jl = gsl_sf_bessel_jl (j, x);

      ncm_assert_cmpdouble_e (ncm_jl, ==, gsl_jl, 1.0e-7, 0.0);
    }
  }
}

static gdouble
_gsl_sf_bessel_jl (const gdouble x, gpointer user_data)
{
  guint *j = (guint *) user_data;

  /*printf ("%u % 22.15g % 22.15g % 22.15g\n", j[0], x, gsl_sf_bessel_jl (j[0], x), gsl_sf_bessel_jl (j[0], x) / x);*/

  return gsl_sf_bessel_jl (j[0], x);
}

void
test_ncm_sf_sbessel_taylor_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata)
{
  NcmDiff *diff = ncm_diff_new ();
  gdouble jla[4];
  guint i, j;

  for (j = 0; j <= L; j++)
  {
    for (i = 0; i < NTOT; i++)
    {
      const gdouble x        = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * i) + 1.0e-2;
      const gdouble gsl_jl   = gsl_sf_bessel_jl (j, x);
      const gdouble gsl_djl  = ncm_diff_rf_d1_1_to_1 (diff, x, _gsl_sf_bessel_jl, &j, NULL);
      const gdouble gsl_d2jl = ncm_diff_rc_d2_1_to_1 (diff, x, _gsl_sf_bessel_jl, &j, NULL) / 2.0;

      ncm_sf_sbessel_taylor (j, x, jla);

      ncm_assert_cmpdouble_e (jla[0], ==, gsl_jl,   1.0e-7, 0.0);
      ncm_assert_cmpdouble_e (jla[1], ==, gsl_djl,  1.0e-3, 0.0);
      ncm_assert_cmpdouble_e (jla[2], ==, gsl_d2jl, 1.0e-3, 0.0);
    }
  }

  ncm_diff_free (diff);
}

void
test_ncm_sf_sbessel_spline_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata)
{
  guint i, j;

  for (j = 0; j <= 10; j++)
  {
    const gdouble xi = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * 0) + 1.0e-2;
    const gdouble xf = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * (NTOT - 1.0)) + 2.0e-2;
    NcmSpline *s     = ncm_sf_sbessel_spline (j, xi, xf, 1.0e-5);

    for (i = 0; i < NTOT; i++)
    {
      const gdouble x      = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * i) + 1.0e-2;
      const gdouble ncm_jl = ncm_spline_eval (s, x);
      const gdouble gsl_jl = gsl_sf_bessel_jl (j, x);

      ncm_assert_cmpdouble_e (ncm_jl, ==, gsl_jl, 1.0e-3, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_sf_sbessel_traps (TestNcmSFSBessel *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/sf/sbessel/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_sf_sbessel_invalid_st (TestNcmSFSBessel *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

