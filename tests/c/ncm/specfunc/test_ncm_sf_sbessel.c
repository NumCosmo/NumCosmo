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
void test_ncm_sf_sbessel_deriv_from_array_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_exact_argument (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_array_cmp_mp (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_array_cutoff (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_array_large_x (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_array_negative (TestNcmSFSBessel *test, gconstpointer pdata);
void test_ncm_sf_sbessel_array_table (TestNcmSFSBessel *test, gconstpointer pdata);

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

  g_test_add ("/ncm/sf/sbessel/deriv_from_array/cmp/gsl", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_deriv_from_array_cmp_gsl,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/exact_argument", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_exact_argument,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/array/cmp/mp", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_array_cmp_mp,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/array/cutoff", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_array_cutoff,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/array/large_x", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_array_large_x,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/array/negative", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_array_negative,
              &test_ncm_sf_sbessel_free);

  g_test_add ("/ncm/sf/sbessel/array/table", TestNcmSFSBessel, NULL,
              &test_ncm_sf_sbessel_new,
              &test_ncm_sf_sbessel_array_table,
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
    g_assert_cmpfloat (ncm_sf_sbessel (j, 0.0), ==, (j == 0) ? 1.0 : 0.0);

    for (i = 0; i < NTOT; i++)
    {
      const gdouble x      = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * i);
      const gdouble ncm_jl = ncm_sf_sbessel (j, x);
      const gdouble gsl_jl = gsl_sf_bessel_jl (j, x);
      const gdouble env    = (x < j) ? fabs (gsl_jl) : GSL_MAX (fabs (gsl_jl), 1.0 / x);

      /* Scaled by the envelope, not the value, which vanishes at the zeros; measured 4.3e-13 */
      if (x > 0.0)
        g_assert_cmpfloat (fabs (ncm_jl - gsl_jl), <=, 5.0e-12 * env);
    }
  }
}

static gdouble
_gsl_sf_bessel_jl (const gdouble x, gpointer user_data)
{
  guint *j = (guint *) user_data;


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
test_ncm_sf_sbessel_deriv_from_array_cmp_gsl (TestNcmSFSBessel *test, gconstpointer pdata)
{
  NcmSFSBesselArray *sba = ncm_sf_sbessel_array_new ();
  gdouble *jl_x          = g_new (gdouble, L + 1);
  guint i, j;

  for (i = 0; i < NTOT; i++)
  {
    const gdouble x = 120.0 / NTOT * (i + 1);

    ncm_sf_sbessel_array_eval (sba, L, x, jl_x);

    for (j = 0; j <= L; j++)
    {
      /* j_l'(x) = (l j_{l-1} - (l+1) j_{l+1}) / (2 l + 1) via GSL values */
      const gdouble gsl_jlm1 = (j > 0) ? gsl_sf_bessel_jl (j - 1, x) : cos (x) / x;
      const gdouble gsl_jlp1 = gsl_sf_bessel_jl (j + 1, x);
      const gdouble gsl_djl  = (j * gsl_jlm1 - (j + 1.0) * gsl_jlp1) / (2.0 * j + 1.0);
      const gdouble gsl_jl   = gsl_sf_bessel_jl (j, x);
      const gdouble gsl_dxjl = gsl_jl + x * gsl_djl;
      const gdouble ncm_djl  = ncm_sf_sbessel_jl_deriv_from_array (j, x, jl_x);
      const gdouble ncm_dxjl = ncm_sf_sbessel_xjl_deriv_from_array (j, x, jl_x);
      const gdouble scale_d  = GSL_MAX (fabs (gsl_jlm1), fabs (gsl_jlp1)) + GSL_DBL_MIN;
      const gdouble scale_dx = x * scale_d + GSL_DBL_MIN;

      /* deep in the evanescent tail the array underflows to zero by design */
      if (jl_x[j] == 0.0)
        continue;

      ncm_assert_cmpdouble_e (ncm_djl, ==, gsl_djl, 1.0e-11, 1.0e-11 * scale_d);
      ncm_assert_cmpdouble_e (ncm_dxjl, ==, gsl_dxjl, 1.0e-11, 1.0e-11 * scale_dx);
    }
  }

  /* l = 0 and l = 1 at x = 0 */
  ncm_sf_sbessel_array_eval (sba, L, 0.0, jl_x);
  ncm_assert_cmpdouble_e (ncm_sf_sbessel_jl_deriv_from_array (0, 0.0, jl_x), ==, 0.0, 0.0, 1.0e-15);
  ncm_assert_cmpdouble_e (ncm_sf_sbessel_jl_deriv_from_array (1, 0.0, jl_x), ==, 1.0 / 3.0, 1.0e-15, 0.0);

  g_free (jl_x);
  ncm_sf_sbessel_array_free (sba);
}

/*
 * The argument is taken exactly. Truth from mpmath at 60 digits; errors scaled by the
 * oscillation amplitude 1/x, the second point sits at a zero. The first is where the old
 * continued-fraction conversion truncated x, off by 1.7e-12 of the amplitude.
 */
void
test_ncm_sf_sbessel_exact_argument (TestNcmSFSBessel *test, gconstpointer pdata)
{
  const guint l_a[]       = {78, 67, 1000};
  const gdouble x_a[]     = {7445.477961962307, 569.3348020587918, 1234.56789};
  const gdouble truth_a[] = {-4.189145352754247670805325e-5, -6.266742084230278044156579e-9, -4.367534182255071180371313e-4};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (l_a); i++)
    g_assert_cmpfloat (fabs (ncm_sf_sbessel (l_a[i], x_a[i]) - truth_a[i]), <=, 2.0 * GSL_DBL_EPSILON / x_a[i]);
}

/*
 * Error of an array value against the multiple precision j_l(x), scaled by the envelope:
 * |j_l(x)| below the turning point, where it decays monotonically, and the oscillation
 * amplitude 1/x above it.
 */
static gdouble
_test_array_err (guint l, gdouble x, gdouble jl)
{
  const gdouble truth = ncm_sf_sbessel (l, x);
  const gdouble env   = (fabs (x) < l) ? fabs (truth) : GSL_MAX (fabs (truth), 1.0 / fabs (x));

  return fabs (jl - truth) / env;
}

void
test_ncm_sf_sbessel_array_cmp_mp (TestNcmSFSBessel *test, gconstpointer pdata)
{
  /* Every branch: Taylor below 2.4e-4, Steed/Barnett up to lmax + 1, upward above */
  const gdouble x_a[] = {
    1.0e-5, 2.3e-4, 2.5e-4, 1.0e-3, 0.1, 1.0, 7.3, 30.0, 99.5, 101.5,
    300.0, 999.0, 1001.5, 1500.0, 1999.0, 2001.5, 5.0e3, 2.0e4, 1.0e6
  };
  const guint l_a[]      = {0, 1, 2, 3, 5, 10, 50, 100, 200, 500, 1000, 1500, 1990, 2000};
  const guint lmax       = 2000;
  NcmSFSBesselArray *sba = ncm_sf_sbessel_array_new_full (lmax, 1.0e-100);
  gdouble *jl_x          = g_new (gdouble, lmax + 1);
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (x_a); i++)
  {
    const gdouble x = x_a[i];
    const guint cut = ncm_sf_sbessel_array_eval_ell_cutoff (sba, x);

    ncm_sf_sbessel_array_eval (sba, lmax, x, jl_x);

    for (j = 0; j < G_N_ELEMENTS (l_a); j++)
    {
      /* Measured at most 1.5e-13 */
      if (l_a[j] <= cut)
        g_assert_cmpfloat (_test_array_err (l_a[j], x, jl_x[l_a[j]]), <, 1.0e-12);
    }
  }

  ncm_sf_sbessel_array_eval (sba, lmax, 0.0, jl_x);
  g_assert_cmpfloat (jl_x[0], ==, 1.0);

  for (j = 1; j <= lmax; j++)
    g_assert_cmpfloat (jl_x[j], ==, 0.0);

  g_free (jl_x);
  ncm_sf_sbessel_array_free (sba);
}

void
test_ncm_sf_sbessel_array_cutoff (TestNcmSFSBessel *test, gconstpointer pdata)
{
  const gdouble thr_a[] = {1.0e-300, 1.0e-100, 1.0e-30, 1.0e-10};
  const gdouble x_a[]   = {0.01, 0.3, 1.0, 3.7, 10.0, 55.5, 200.0, 777.0, 1500.0, 2500.0};
  const guint lmax      = 3000;
  gdouble *jl_x         = g_new (gdouble, lmax + 1);
  guint i, j, l;

  for (i = 0; i < G_N_ELEMENTS (thr_a); i++)
  {
    NcmSFSBesselArray *sba = ncm_sf_sbessel_array_new_full (lmax, thr_a[i]);

    for (j = 0; j < G_N_ELEMENTS (x_a); j++)
    {
      const gdouble x = x_a[j];
      const guint cut = ncm_sf_sbessel_array_eval_ell_cutoff (sba, x);

      /* Nothing is cut at this x */
      if (cut >= lmax)
        continue;

      ncm_sf_sbessel_array_eval (sba, lmax, x, jl_x);

      /* The orders above the cutoff are zero and below the threshold */
      for (l = cut + 1; l <= GSL_MIN (cut + 5, lmax); l++)
      {
        g_assert_cmpfloat (jl_x[l], ==, 0.0);
        g_assert_cmpfloat (fabs (ncm_sf_sbessel (l, x)), <=, thr_a[i]);
      }

      /*
       * The cutoff is not far above the crossing: the old estimate started the recursion
       * where |j_l| was 1e-211 times the threshold and overflowed. Measured at least 1e-4.
       */
      g_assert_cmpfloat (fabs (ncm_sf_sbessel (cut, x)), >, 1.0e-6 * thr_a[i]);
    }

    ncm_sf_sbessel_array_free (sba);
  }

  g_free (jl_x);
}

void
test_ncm_sf_sbessel_array_large_x (TestNcmSFSBessel *test, gconstpointer pdata)
{
  /* Default lmax and threshold: these returned inf and NaN before the Debye cutoff */
  const gdouble x_a[]    = {2640.0, 3000.0, 5000.0, 9000.0};
  NcmSFSBesselArray *sba = ncm_sf_sbessel_array_new ();
  const guint lmax       = ncm_sf_sbessel_array_get_lmax (sba);
  gdouble *jl_x          = g_new (gdouble, lmax + 1);
  guint i, l;

  for (i = 0; i < G_N_ELEMENTS (x_a); i++)
  {
    const gdouble x = x_a[i];

    ncm_sf_sbessel_array_eval (sba, lmax, x, jl_x);

    for (l = 0; l <= lmax; l++)
      g_assert_true (gsl_finite (jl_x[l]));

    g_assert_cmpfloat (_test_array_err (0, x, jl_x[0]), <, 1.0e-12);
    g_assert_cmpfloat (_test_array_err ((guint) x, x, jl_x[(guint) x]), <, 1.0e-12);
  }

  g_free (jl_x);
  ncm_sf_sbessel_array_free (sba);
}

void
test_ncm_sf_sbessel_array_negative (TestNcmSFSBessel *test, gconstpointer pdata)
{
  const gdouble x_a[]    = {1.0e-4, 0.5, 3.0, 30.0, 1500.0};
  const guint lmax       = 2000;
  NcmSFSBesselArray *sba = ncm_sf_sbessel_array_new_full (lmax, 1.0e-100);
  gdouble *jl_p          = g_new (gdouble, lmax + 1);
  gdouble *jl_m          = g_new (gdouble, lmax + 1);
  guint i, l;

  for (i = 0; i < G_N_ELEMENTS (x_a); i++)
  {
    const gdouble x = x_a[i];

    g_assert_cmpuint (ncm_sf_sbessel_array_eval_ell_cutoff (sba, -x), ==, ncm_sf_sbessel_array_eval_ell_cutoff (sba, x));

    ncm_sf_sbessel_array_eval (sba, lmax, x, jl_p);
    ncm_sf_sbessel_array_eval (sba, lmax, -x, jl_m);

    for (l = 0; l <= lmax; l++)
      g_assert_cmpfloat (jl_m[l], ==, (l % 2 == 0) ? jl_p[l] : -jl_p[l]);

    g_assert_cmpfloat (_test_array_err (1, -x, jl_m[1]), <, 1.0e-12);
  }

  g_free (jl_p);
  g_free (jl_m);
  ncm_sf_sbessel_array_free (sba);
}

void
test_ncm_sf_sbessel_array_table (TestNcmSFSBessel *test, gconstpointer pdata)
{
  const gdouble x_a[]      = {0.5, 3.0, 40.0, 700.0};
  const guint ell_max      = 300;
  NcmSFSBesselArray *sba   = ncm_sf_sbessel_array_new_full (1000, 1.0e-100);
  NcmSFSBesselArray *sba_t = ncm_sf_sbessel_array_new_full (1000, 1.0e-10);
  GArray *x                = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *x_copy           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble *jl_x            = g_new (gdouble, ell_max + 1);
  NcmMatrix *table, *same, *other_thr, *other_ell;
  guint i, l;

  g_array_append_vals (x, x_a, G_N_ELEMENTS (x_a));
  g_array_append_vals (x_copy, x_a, G_N_ELEMENTS (x_a));

  table = ncm_sf_sbessel_array_ref_table (sba, x, ell_max);
  g_assert_cmpuint (ncm_matrix_nrows (table), ==, x->len);
  g_assert_cmpuint (ncm_matrix_ncols (table), ==, ell_max + 1);

  for (i = 0; i < x->len; i++)
  {
    ncm_sf_sbessel_array_eval (sba, ell_max, x_a[i], jl_x);

    for (l = 0; l <= ell_max; l++)
      g_assert_cmpfloat (ncm_matrix_get (table, i, l), ==, jl_x[l]);
  }

  /* Keyed on the values of the abscissae, not on the array */
  same = ncm_sf_sbessel_array_ref_table (sba, x_copy, ell_max);
  g_assert_true (same == table);

  /* The threshold changes the rows, so it is part of the key */
  other_thr = ncm_sf_sbessel_array_ref_table (sba_t, x, ell_max);
  g_assert_true (other_thr != table);

  ncm_sf_sbessel_array_eval (sba_t, ell_max, x_a[0], jl_x);

  for (l = 0; l <= ell_max; l++)
    g_assert_cmpfloat (ncm_matrix_get (other_thr, 0, l), ==, jl_x[l]);

  other_ell = ncm_sf_sbessel_array_ref_table (sba, x, ell_max - 1);
  g_assert_true (other_ell != table);

  ncm_matrix_free (table);
  ncm_matrix_free (same);
  ncm_matrix_free (other_thr);
  ncm_matrix_free (other_ell);
  g_array_unref (x);
  g_array_unref (x_copy);
  g_free (jl_x);
  ncm_sf_sbessel_array_free (sba);
  ncm_sf_sbessel_array_free (sba_t);
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

