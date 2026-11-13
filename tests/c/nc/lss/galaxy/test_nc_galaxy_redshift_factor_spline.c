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
  NcmSpline *pz;
} TestNcGalaxyRedshiftFactorSpline;

void test_nc_galaxy_redshift_factor_spline_new_poly (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_new_zeros (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_free (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);

void test_nc_galaxy_redshift_factor_spline_poly_eval (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_poly_norm (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_poly_gen (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_zeros (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_obs (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_table_changed (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_threads (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_spline_errors (void);
void test_nc_galaxy_redshift_factor_spline_negative_subprocess (void);
void test_nc_galaxy_redshift_factor_spline_few_knots_subprocess (void);

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

  g_test_add ("/nc/galaxy_redshift_factor_spline/table_changed", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_poly,
              &test_nc_galaxy_redshift_factor_spline_table_changed,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_add ("/nc/galaxy_redshift_factor_spline/threads", TestNcGalaxyRedshiftFactorSpline, NULL,
              &test_nc_galaxy_redshift_factor_spline_new_poly,
              &test_nc_galaxy_redshift_factor_spline_threads,
              &test_nc_galaxy_redshift_factor_spline_free);

  g_test_add_func ("/nc/galaxy_redshift_factor_spline/errors", &test_nc_galaxy_redshift_factor_spline_errors);
  g_test_add_func ("/nc/galaxy_redshift_factor_spline/errors/negative/subprocess",
                   &test_nc_galaxy_redshift_factor_spline_negative_subprocess);
  g_test_add_func ("/nc/galaxy_redshift_factor_spline/errors/few_knots/subprocess",
                   &test_nc_galaxy_redshift_factor_spline_few_knots_subprocess);

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

static NcmSpline *
_test_table_spline (const gdouble *z, const gdouble *p, const guint n)
{
  NcmVector *xv = ncm_vector_new_data_dup ((gdouble *) z, n, 1);
  NcmVector *yv = ncm_vector_new_data_dup ((gdouble *) p, n, 1);
  NcmSpline *pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return pz;
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

  test->gsdr = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());
  test->mset = ncm_mset_empty_new ();
  test->data = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);
  test->pz   = _test_table_spline (z, p, TEST_POLY_N);

  nc_galaxy_redshift_factor_spline_data_set (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (test->gsdr), test->data, test->pz);
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

  test->gsdr = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());
  test->mset = ncm_mset_empty_new ();
  test->data = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);
  test->pz   = _test_table_spline (z, p, n);

  nc_galaxy_redshift_factor_spline_data_set (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (test->gsdr), test->data, test->pz);
}

void
test_nc_galaxy_redshift_factor_spline_free (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  nc_galaxy_redshift_factor_data_unref (test->data);
  ncm_spline_free (test->pz);
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
  NcmVector *xv                          = ncm_spline_peek_xv (test->pz);
  NcmVector *yv                          = ncm_spline_peek_yv (test->pz);
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

/* A row is read from and written to the NcGalaxyWLObs P(z) slot as the same spline. */
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

/* The data derived from a table is shared by every galaxy that reads the same spline,
 * and is rebuilt, not reused, after the table's values are changed in place: scaling the
 * values by 4 scales the normalization by 4, for a new galaxy and for one that reads the
 * table again. The galaxy that keeps the old data keeps the old normalization. */
void
test_nc_galaxy_redshift_factor_spline_table_changed (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  const gdouble exact               = _test_q2_antideriv (TEST_POLY_Z_MAX) - _test_q2_antideriv (TEST_POLY_Z_MIN);
  NcGalaxyRedshiftFactorSpline *gs  = NC_GALAXY_REDSHIFT_FACTOR_SPLINE (test->gsdr);
  NcGalaxyRedshiftFactorData *data2 = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);
  NcGalaxyRedshiftFactorData *data3 = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);

  nc_galaxy_redshift_factor_spline_data_set (gs, data2, test->pz);
  ncm_assert_cmpdouble (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, data2), ==,
                        nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, test->data));

  ncm_vector_scale (ncm_spline_peek_yv (test->pz), 4.0);

  nc_galaxy_redshift_factor_spline_data_set (gs, data2, test->pz);
  nc_galaxy_redshift_factor_spline_data_set (gs, data3, test->pz);

  ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, data2), ==, 4.0 * exact, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, data3), ==, 4.0 * exact, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, test->data), ==, exact, 1.0e-14, 0.0);

  nc_galaxy_redshift_factor_data_unref (data2);
  nc_galaxy_redshift_factor_data_unref (data3);
}

void
test_nc_galaxy_redshift_factor_spline_errors (void)
{
  g_test_trap_subprocess ("/nc/galaxy_redshift_factor_spline/errors/negative/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the P(z) table value 2 is -0.1 at z = 0.52; values must be non-negative*");

  g_test_trap_subprocess ("/nc/galaxy_redshift_factor_spline/errors/few_knots/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the P(z) table has 5 knots, at least 6 are required*");
}

void
test_nc_galaxy_redshift_factor_spline_negative_subprocess (void)
{
  const gdouble z[]                = {0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56};
  const gdouble p[]                = {0.1, 0.5, -0.1, 0.5, 1.0, 0.5, 0.1};
  NcGalaxyRedshiftFactor *gsdr     = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());
  NcmMSet *mset                    = ncm_mset_empty_new ();
  NcGalaxyRedshiftFactorData *data = nc_galaxy_redshift_factor_data_new (gsdr, mset);
  NcmSpline *pz                    = _test_table_spline (z, p, G_N_ELEMENTS (p));

  nc_galaxy_redshift_factor_spline_data_set (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (gsdr), data, pz);
}

void
test_nc_galaxy_redshift_factor_spline_few_knots_subprocess (void)
{
  const gdouble z[]                = {0.5, 0.51, 0.52, 0.53, 0.54};
  const gdouble p[]                = {0.1, 0.5, 1.0, 0.5, 0.1};
  NcGalaxyRedshiftFactor *gsdr     = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());
  NcmMSet *mset                    = ncm_mset_empty_new ();
  NcGalaxyRedshiftFactorData *data = nc_galaxy_redshift_factor_data_new (gsdr, mset);
  NcmSpline *pz                    = NCM_SPLINE (ncm_spline_gsl_new_by_id (NCM_SPLINE_GSL_LINEAR));
  NcmVector *xv                    = ncm_vector_new_data_dup ((gdouble *) z, G_N_ELEMENTS (z), 1);
  NcmVector *yv                    = ncm_vector_new_data_dup ((gdouble *) p, G_N_ELEMENTS (p), 1);

  ncm_spline_set (pz, xv, yv, TRUE);

  nc_galaxy_redshift_factor_spline_data_set (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (gsdr), data, pz);
}

#define TEST_THREADS_N (8)
#define TEST_THREADS_ROUNDS (40)
#define TEST_THREADS_KNOTS (4000)

typedef struct _TestThreadsArg
{
  NcGalaxyRedshiftFactorSpline *gs;
  NcGalaxyRedshiftFactorData *data;
  NcmSpline *pz;
  gint *ready;
  gint *go;
} TestThreadsArg;

/* Waits until every thread is ready and the start is given, then reads the table. */
static gpointer
_test_threads_read (gpointer userdata)
{
  TestThreadsArg *arg = (TestThreadsArg *) userdata;

  g_atomic_int_inc (arg->ready);

  while (!g_atomic_int_get (arg->go))
    ;

  nc_galaxy_redshift_factor_spline_data_set (arg->gs, arg->data, arg->pz);

  return NULL;
}

/* Threads that read the same new table at the same time all get the density of that
 * table: each round starts TEST_THREADS_N threads on a new table of q^2 at once, so
 * several of them build its data together and all but one drop their build. */
void
test_nc_galaxy_redshift_factor_spline_threads (TestNcGalaxyRedshiftFactorSpline *test, gconstpointer pdata)
{
  const gdouble exact              = _test_q2_antideriv (TEST_POLY_Z_MAX) - _test_q2_antideriv (TEST_POLY_Z_MIN);
  NcGalaxyRedshiftFactorSpline *gs = NC_GALAXY_REDSHIFT_FACTOR_SPLINE (test->gsdr);
  gdouble *z                       = g_new (gdouble, TEST_THREADS_KNOTS);
  gdouble *p                       = g_new (gdouble, TEST_THREADS_KNOTS);
  guint round, t, j;

  for (j = 0; j < TEST_THREADS_KNOTS; j++)
  {
    z[j] = TEST_POLY_Z_MIN + (TEST_POLY_Z_MAX - TEST_POLY_Z_MIN) * pow (j / (TEST_THREADS_KNOTS - 1.0), 1.3);
    p[j] = gsl_pow_2 (_test_q (z[j]));
  }

  for (round = 0; round < TEST_THREADS_ROUNDS; round++)
  {
    NcmSpline *pz = _test_table_spline (z, p, TEST_THREADS_KNOTS);
    TestThreadsArg arg[TEST_THREADS_N];
    GThread *th[TEST_THREADS_N];
    gint ready = 0;
    gint go    = 0;

    for (t = 0; t < TEST_THREADS_N; t++)
    {
      arg[t].gs    = gs;
      arg[t].data  = nc_galaxy_redshift_factor_data_new (test->gsdr, test->mset);
      arg[t].pz    = pz;
      arg[t].ready = &ready;
      arg[t].go    = &go;
      th[t]        = g_thread_new ("read", &_test_threads_read, &arg[t]);
    }

    while (g_atomic_int_get (&ready) < TEST_THREADS_N)
      ;

    g_atomic_int_set (&go, 1);

    for (t = 0; t < TEST_THREADS_N; t++)
      g_thread_join (th[t]);

    {
      const gdouble norm0 = nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, arg[0].data);

      ncm_assert_cmpdouble_e (norm0, ==, exact, 1.0e-13, 0.0);

      for (t = 0; t < TEST_THREADS_N; t++)
      {
        g_assert_true (nc_galaxy_redshift_factor_spline_data_peek (gs, arg[t].data) == pz);
        ncm_assert_cmpdouble (nc_galaxy_redshift_factor_norm (test->gsdr, test->mset, arg[t].data), ==, norm0);
        nc_galaxy_redshift_factor_data_unref (arg[t].data);
      }
    }

    ncm_spline_free (pz);
  }

  g_free (z);
  g_free (p);
}

