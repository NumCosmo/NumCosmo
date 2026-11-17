/***************************************************************************
 *            test_nc_galaxy_redshift_factor_spline.c
 *
 *  Fri Nov 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_redshift_factor_spline.c
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

/* sqrt(p) is the polynomial Q_COEF[0] + Q_COEF[1] z + ... of degree 5, positive on the
 * table; the order-6 B-spline reproduces it exactly, so p = q^2 is known in closed form. */
#define TEST_Q_DEGREE (5)
#define TEST_POLY_N (40)
#define TEST_POLY_Z_MIN (0.2)
#define TEST_POLY_Z_MAX (1.4)

static const gdouble Q_COEF[TEST_Q_DEGREE + 1] = {1.0, 0.5, -0.8, 0.3, 0.1, -0.05};

typedef struct _TestNcGalaxyRedshiftFactorSpline
{
  NcGalaxyRedshiftFactor *gsdr;
  NcGalaxyRedshiftFactorData *data;
  NcmMSet *mset;
  NcmVector *zv;
  NcmVector *pv;
  NcGalaxyPz *pz;
} TestNcGalaxyRedshiftFactorSpline;

void test_nc_galaxy_redshift_factor_spline_new_poly (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_new_zeros (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_free (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);

void test_nc_galaxy_redshift_factor_spline_poly_eval (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_poly_norm (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_poly_gen (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_zeros (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_obs (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/galaxy_redshift_factor_spline/poly/eval", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_poly,
              &test_nc_galaxy_redshift_factor_spline_poly_eval,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_add ("/nc/galaxy_redshift_factor_spline/poly/norm", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_poly,
              &test_nc_galaxy_redshift_factor_spline_poly_norm,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_add ("/nc/galaxy_redshift_factor_spline/poly/gen", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_poly,
              &test_nc_galaxy_redshift_factor_spline_poly_gen,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_add ("/nc/galaxy_redshift_factor_spline/zeros", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_zeros,
              &test_nc_galaxy_redshift_factor_spline_zeros,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_add ("/nc/galaxy_redshift_factor_spline/obs", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_poly,
              &test_nc_galaxy_redshift_factor_spline_obs,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_run ();

  return 0;
}

static gdouble
_test_q (const gdouble z)
{
  gdouble q = 0.0;
  gint i;

  for (i = TEST_Q_DEGREE; i >= 0; i--)
    q = q * z + Q_COEF[i];

  return q;
}

/* Antiderivative of q^2: the coefficients of q^2 are the self-convolution of Q_COEF. */
static gdouble
_test_q2_antideriv (const gdouble z)
{
  gdouble res = 0.0;
  gint m;

  for (m = 2 * TEST_Q_DEGREE; m >= 0; m--)
  {
    gdouble c_m = 0.0;
    gint i;

    for (i = MAX (0, m - TEST_Q_DEGREE); i <= MIN (m, TEST_Q_DEGREE); i++)
      c_m += Q_COEF[i] * Q_COEF[m - i];

    res = res * z + c_m / (m + 1.0);
  }

  return res * z;
}

static gint
_test_cmp_double (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

/* Creates the scheme and a galaxy whose density is that of the table (z_j, p_j) */
static void
_test_set_table (TestNcGalaxyRedshiftFactorSpline *test, const gdouble *z, const gdouble *p, const guint n)
{
  test->gsdr = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());
  test->mset = ncm_mset_empty_new ();
  test->data = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);
  test->zv   = ncm_vector_new_data_dup ((gdouble *) z, n, 1);
  test->pv   = ncm_vector_new_data_dup ((gdouble *) p, n, 1);
  test->pz   = NC_GALAXY_PZ (nc_galaxy_pz_bspline_sqrt_new_from_table (test->zv, test->pv));

  nc_galaxy_redshift_factor_spline_data_set (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (test->gsdr), test->data, test->pz);
}

/* Non-uniform knots, so the breakpoints of the B-spline are not the knots. */
void
test_nc_galaxy_redshift_factor_spline_new_poly (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  gdouble z[TEST_POLY_N], p[TEST_POLY_N];
  guint j;

  for (j = 0; j < TEST_POLY_N; j++)
  {
    z[j] = TEST_POLY_Z_MIN + (TEST_POLY_Z_MAX - TEST_POLY_Z_MIN) * pow (j / (TEST_POLY_N - 1.0), 1.3);
    p[j] = gsl_pow_2 (_test_q (z[j]));
  }

  _test_set_table (test, z, p, TEST_POLY_N);
}

/* A cut table as written by a photo-z code: zeros at interior knots and a narrow peak. */
void
test_nc_galaxy_redshift_factor_spline_new_zeros (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  const gdouble p[] = {0.02, 0.1, 0.6, 3.0, 0.6, 0.0, 0.0, 0.0, 0.03, 1.0, 2.0, 1.0, 0.0, 0.0, 0.05, 0.02};
  const guint n     = G_N_ELEMENTS (p);
  gdouble z[G_N_ELEMENTS (p)];
  guint j;

  for (j = 0; j < n; j++)
    z[j] = 0.5 + 0.01 * j;

  _test_set_table (test, z, p, n);
}

void
test_nc_galaxy_redshift_factor_spline_free (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  nc_galaxy_redshift_factor_data_unref (test->data);
  nc_galaxy_pz_free (test->pz);
  ncm_vector_free (test->zv);
  ncm_vector_free (test->pv);
  ncm_mset_free (test->mset);
  nc_galaxy_redshift_factor_free (test->gsdr);
}

/* The density is q^2 and its logarithm 2 ln q everywhere on the table, and the
 * integration limits are the first and last knots. */
void
test_nc_galaxy_redshift_factor_spline_poly_eval (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  NcGalaxyRedshiftFactorIntegrand *integ    = nc_galaxy_redshift_factor_integ (test->gsdr, test->mset, FALSE);
  NcGalaxyRedshiftFactorIntegrand *ln_integ = nc_galaxy_redshift_factor_integ (test->gsdr, test->mset, TRUE);
  gdouble z_min, z_max;
  guint i;

  nc_galaxy_redshift_factor_get_integ_lim (test->gsdr, test->mset, test->data, &z_min, &z_max);
  g_assert_cmpfloat (z_min, ==, TEST_POLY_Z_MIN);
  g_assert_cmpfloat (z_max, ==, TEST_POLY_Z_MAX);

  for (i = 0; i <= 1000; i++)
  {
    const gdouble z = z_min + (z_max - z_min) * i / 1000.0;
    const gdouble q = _test_q (z);

    g_assert_cmpfloat (q, >, 0.0);
    ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_integrand_eval (integ, z, test->data), ==, q * q, 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_integrand_eval (ln_integ, z, test->data), ==, 2.0 * log (q), 0.0, 1.0e-13);
  }

  nc_galaxy_redshift_factor_integrand_free (integ);
  nc_galaxy_redshift_factor_integrand_free (ln_integ);
}

/* The normalization and the fixed-node integral of the degree-10 polynomial q^2 are
 * exact. */
void
test_nc_galaxy_redshift_factor_spline_poly_norm (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  const gdouble exact = _test_q2_antideriv (TEST_POLY_Z_MAX) - _test_q2_antideriv (TEST_POLY_Z_MIN);
  NcmIntegralFixed *intf;

  ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, test->data), ==, exact, 1.0e-14, 0.0);

  intf = nc_galaxy_redshift_factor_make_fixed_nodes (test->gsdr, test->mset, test->data, TEST_POLY_Z_MIN, TEST_POLY_Z_MAX, 3, 6);
  ncm_assert_cmpdouble_e (ncm_integral_fixed_nodes_eval (intf), ==, exact, 1.0e-14, 0.0);
  ncm_integral_fixed_free (intf);
}

/* The draws follow q^2: Kolmogorov-Smirnov distance to the exact CDF below the 0.1%
 * critical value 1.95 / sqrt(N). */
void
test_nc_galaxy_redshift_factor_spline_poly_gen (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  const guint n_draws = 20000;
  const gdouble F_min = _test_q2_antideriv (TEST_POLY_Z_MIN);
  const gdouble total = _test_q2_antideriv (TEST_POLY_Z_MAX) - F_min;
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 1);
  GArray *draws       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_draws);
  gdouble ks          = 0.0;
  guint i;

  for (i = 0; i < n_draws; i++)
  {
    nc_galaxy_redshift_factor_gen (test->gsdr, test->mset, test->data, rng);
    g_assert_cmpfloat (test->data->z, >=, TEST_POLY_Z_MIN);
    g_assert_cmpfloat (test->data->z, <=, TEST_POLY_Z_MAX);
    g_array_append_val (draws, test->data->z);
  }

  g_array_sort (draws, &_test_cmp_double);

  for (i = 0; i < n_draws; i++)
  {
    const gdouble F = (_test_q2_antideriv (g_array_index (draws, gdouble, i)) - F_min) / total;

    ks = MAX (ks, MAX (fabs (F - (gdouble) i / n_draws), fabs (F - (i + 1.0) / n_draws)));
  }

  g_assert_cmpfloat (ks, <, 1.95 / sqrt (n_draws));

  g_array_unref (draws);
  ncm_rng_free (rng);
}

/* On a cut table the density is non-negative everywhere and equals the table at the
 * knots, to round-off (1e-30) at the zero knots. */
void
test_nc_galaxy_redshift_factor_spline_zeros (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  NcGalaxyRedshiftFactorIntegrand *integ = nc_galaxy_redshift_factor_integ (test->gsdr, test->mset, FALSE);
  NcmVector *xv                          = test->zv;
  NcmVector *yv                          = test->pv;
  const guint n                          = ncm_vector_len (xv);
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble z = ncm_vector_get (xv, i);
    const gdouble p = ncm_vector_get (yv, i);

    ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_integrand_eval (integ, z, test->data), ==, p, 1.0e-13, 1.0e-30);
  }

  for (i = 0; i <= 3000; i++)
  {
    const gdouble z = ncm_vector_get (xv, 0) + (ncm_vector_get (xv, n - 1) - ncm_vector_get (xv, 0)) * i / 3000.0;

    g_assert_cmpfloat (nc_galaxy_redshift_factor_integrand_eval (integ, z, test->data), >=, 0.0);
  }

  nc_galaxy_redshift_factor_integrand_free (integ);
}

/* A row is read from and written to the pz slot of the NcGalaxyWLObs as the same object. */
void
test_nc_galaxy_redshift_factor_spline_obs (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  GStrv columns                     = g_strsplit ("z", " ", -1);
  NcGalaxyWLObs *obs                = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, 2, columns);
  NcGalaxyRedshiftFactorData *data2 = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);

  nc_galaxy_wl_obs_set_pz (obs, 1, test->pz);
  nc_galaxy_redshift_factor_data_read_row (data2, obs, 1);
  g_assert_true (nc_galaxy_redshift_factor_spline_data_peek (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (test->gsdr), data2) == test->pz);
  ncm_assert_cmpdouble (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, data2), ==,
                        nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, test->data));

  nc_galaxy_redshift_factor_data_write_row (data2, obs, 0);
  g_assert_true (nc_galaxy_wl_obs_peek_pz (obs, 0) == test->pz);

  nc_galaxy_redshift_factor_data_unref (data2);
  nc_galaxy_wl_obs_free (obs);
  g_strfreev (columns);
}

