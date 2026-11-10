/***************************************************************************
 *            test_ncm_stats_rows.c
 *
 *  Fri Oct 24 12:00:00 2026
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

/*
 * The six analyses of NcmStatsRows, against values captured from the NcmStatsVec
 * implementation they replaced (commit 701a5314), on an AR(1) sample of 1500 rows and 3
 * columns drawn from NcmRNG seed 12345: columns with coefficients 0.9, 0.5 and 0, started
 * at 30, 15 and 6, shifted by 1, 2 and 3. Doubles are compared to TEST_RELTOL (libm, BLAS
 * and LAPACK differ between platforms in the last bits); indices, AR orders and break counts
 * exactly. The analyses run on the rows given to a NcmStatsRows, on a serialized copy of it
 * given the same rows, and on the deprecated NcmStatsVec saved rows.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_NROWS 1500
#define TEST_NCOLS 3
#define TEST_SEED 12345
#define TEST_RELTOL 1.0e-10

typedef struct _TestNcmStatsRows
{
  GPtrArray *sample;
  NcmStatsRows *rows;
  NcmStatsVec *svec;
} TestNcmStatsRows;

/* Counts the rows freed by the arrays of the reset test. */
static guint _test_ncm_stats_rows_nfreed = 0;

/* Reference values of make_stats_rows_cases.py, printed at 701a5314 (the original
 * NcmStatsVec implementation). */
/* BEGIN make_stats_rows_cases.py */
static const gdouble ref_ess[TEST_NCOLS]                 = {80.40645533906867, 532.8605949564081, 1500.0};
static const gint ref_ess_bindex                         = 0;
static const guint ref_ess_wp                            = 0;
static const guint ref_ess_wp_order                      = 0;
static const gdouble ref_pvals[TEST_NCOLS]               = {0.36474076801478694, 0.879391555732373, 0.7936728780752365};
static const gint ref_hd_bindex                          = 150;
static const guint ref_hd_wp                             = 1;
static const guint ref_hd_wp_order                       = 2;
static const gdouble ref_cs_first[3]                     = {-1.0662337610479606, -2.4493441870296446, -4.719078339740692};
static const gdouble ref_cs_last[3]                      = {1169.9410931368232, 1181.5741052595858, 1193.4054846628392};
static const gdouble ref_visual_mean                     = 0.7993338812209221;
static const gdouble ref_visual_var                      = 172804.55373951263;
static const gdouble ref_const_break[TEST_NCOLS]         = {9.0, 1.0, 0.0};
static const gdouble ref_cov_diag[TEST_NCOLS]            = {5.687214613813791, 1.197942771240321, 0.9577769831182713};
static const gdouble ref_cov_ogk[TEST_NCOLS][TEST_NCOLS] = {
  {5.7821562810551645, 0.24316839146534633, -0.013580397437948503},
  {0.24316839146534633, 1.2195039935012058, -0.030518580734787224},
  {-0.013580397437948503, -0.030518580734787224, 0.9614220826107162},
};

/* END make_stats_rows_cases.py */

static void test_ncm_stats_rows_new (TestNcmStatsRows *test, gconstpointer pdata);
static void test_ncm_stats_rows_free (TestNcmStatsRows *test, gconstpointer pdata);

static void _test_ncm_stats_rows_check (NcmStatsRows *rows);

static void test_ncm_stats_rows_reference (TestNcmStatsRows *test, gconstpointer pdata);
static void test_ncm_stats_rows_reference_serialized (TestNcmStatsRows *test, gconstpointer pdata);
static void test_ncm_stats_rows_reference_stats_vec (TestNcmStatsRows *test, gconstpointer pdata);
static void test_ncm_stats_rows_set_nrows (TestNcmStatsRows *test, gconstpointer pdata);
static void test_ncm_stats_rows_properties (TestNcmStatsRows *test, gconstpointer pdata);
static void test_ncm_stats_rows_reset (TestNcmStatsRows *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

  /* The first NcmStatsVec with save-x warns once per process; later ones in the tests do not. */
  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*save-x is deprecated*");
  ncm_stats_vec_free (ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, TRUE));
  g_test_assert_expected_messages ();

  g_test_add ("/ncm/stats_rows/reference", TestNcmStatsRows, NULL,
              &test_ncm_stats_rows_new,
              &test_ncm_stats_rows_reference,
              &test_ncm_stats_rows_free);
  g_test_add ("/ncm/stats_rows/reference/serialized", TestNcmStatsRows, NULL,
              &test_ncm_stats_rows_new,
              &test_ncm_stats_rows_reference_serialized,
              &test_ncm_stats_rows_free);
  g_test_add ("/ncm/stats_rows/reference/stats_vec", TestNcmStatsRows, NULL,
              &test_ncm_stats_rows_new,
              &test_ncm_stats_rows_reference_stats_vec,
              &test_ncm_stats_rows_free);
  g_test_add ("/ncm/stats_rows/set_nrows", TestNcmStatsRows, NULL,
              &test_ncm_stats_rows_new,
              &test_ncm_stats_rows_set_nrows,
              &test_ncm_stats_rows_free);
  g_test_add ("/ncm/stats_rows/properties", TestNcmStatsRows, NULL,
              &test_ncm_stats_rows_new,
              &test_ncm_stats_rows_properties,
              &test_ncm_stats_rows_free);
  g_test_add ("/ncm/stats_rows/reset", TestNcmStatsRows, NULL,
              &test_ncm_stats_rows_new,
              &test_ncm_stats_rows_reset,
              &test_ncm_stats_rows_free);

  g_test_run ();

  return 0;
}

static void
_test_ncm_stats_rows_count_free (gpointer row)
{
  _test_ncm_stats_rows_nfreed++;
  ncm_vector_free (row);
}

/* Runs the six analyses on @rows and compares them with the reference table. */
static void
_test_ncm_stats_rows_check (NcmStatsRows *rows)
{
  gint ess_bindex, hd_bindex, hd5_bindex;
  guint ess_wp, ess_wp_order, hd_wp, hd_wp_order, hd5_wp, hd5_wp_order;
  gdouble ess_wp_ess, hd_wp_pvalue, hd5_wp_pvalue, cs_mean, cs_var;
  NcmVector *ess      = ncm_stats_rows_max_ess_time (rows, 0, &ess_bindex, &ess_wp, &ess_wp_order, &ess_wp_ess);
  NcmVector *pvals    = ncm_stats_rows_heidel_diag (rows, 0, 0.0, &hd_bindex, &hd_wp, &hd_wp_order, &hd_wp_pvalue);
  NcmVector *pvals5   = ncm_stats_rows_heidel_diag (rows, 5, 0.1, &hd5_bindex, &hd5_wp, &hd5_wp_order, &hd5_wp_pvalue);
  NcmVector *cs       = ncm_stats_rows_visual_heidel_diag (rows, 0, 7, &cs_mean, &cs_var);
  NcmMatrix *cov_diag = ncm_matrix_new (TEST_NCOLS, TEST_NCOLS);
  NcmMatrix *cov_ogk  = ncm_matrix_new (TEST_NCOLS, TEST_NCOLS);
  guint i, j;

  ncm_stats_rows_cov_robust_diag (rows, cov_diag);
  ncm_stats_rows_cov_robust_ogk (rows, cov_ogk);

  g_assert_cmpuint (ncm_stats_rows_get_nrows (rows), ==, TEST_NROWS);
  g_assert_cmpuint (ncm_stats_rows_get_ncols (rows), ==, TEST_NCOLS);

  for (i = 0; i < TEST_NCOLS; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (ess, i), ==, ref_ess[i], TEST_RELTOL, 0.0);

  g_assert_cmpint (ess_bindex, ==, ref_ess_bindex);
  g_assert_cmpuint (ess_wp, ==, ref_ess_wp);
  g_assert_cmpuint (ess_wp_order, ==, ref_ess_wp_order);
  ncm_assert_cmpdouble_e (ess_wp_ess, ==, ref_ess[ref_ess_wp], TEST_RELTOL, 0.0);

  for (i = 0; i < TEST_NCOLS; i++)
  {
    ncm_assert_cmpdouble_e (ncm_vector_get (pvals, i), ==, ref_pvals[i], TEST_RELTOL, 0.0);
    /* Five tests at p-value 0.1 give the same index and values as the defaults here. */
    ncm_assert_cmpdouble_e (ncm_vector_get (pvals5, i), ==, ref_pvals[i], TEST_RELTOL, 0.0);
  }

  g_assert_cmpint (hd_bindex, ==, ref_hd_bindex);
  g_assert_cmpuint (hd_wp, ==, ref_hd_wp);
  g_assert_cmpuint (hd_wp_order, ==, ref_hd_wp_order);
  ncm_assert_cmpdouble_e (hd_wp_pvalue, ==, ref_pvals[ref_hd_wp], TEST_RELTOL, 0.0);
  g_assert_cmpint (hd5_bindex, ==, ref_hd_bindex);
  g_assert_cmpuint (hd5_wp, ==, ref_hd_wp);
  g_assert_cmpuint (hd5_wp_order, ==, ref_hd_wp_order);
  ncm_assert_cmpdouble_e (hd5_wp_pvalue, ==, ref_pvals[ref_hd_wp], TEST_RELTOL, 0.0);

  g_assert_cmpuint (ncm_vector_len (cs), ==, TEST_NROWS - 7);

  for (i = 0; i < 3; i++)
  {
    ncm_assert_cmpdouble_e (ncm_vector_get (cs, i), ==, ref_cs_first[i], TEST_RELTOL, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (cs, TEST_NROWS - 7 - 3 + i), ==, ref_cs_last[i], TEST_RELTOL, 0.0);
  }

  ncm_assert_cmpdouble_e (cs_mean, ==, ref_visual_mean, TEST_RELTOL, 0.0);
  ncm_assert_cmpdouble_e (cs_var, ==, ref_visual_var, TEST_RELTOL, 0.0);

  for (i = 0; i < TEST_NCOLS; i++)
    g_assert_cmpfloat (ncm_stats_rows_estimate_const_break (rows, i), ==, ref_const_break[i]);

  for (i = 0; i < TEST_NCOLS; i++)
  {
    for (j = 0; j < TEST_NCOLS; j++)
    {
      if (i == j)
        ncm_assert_cmpdouble_e (ncm_matrix_get (cov_diag, i, j), ==, ref_cov_diag[i], TEST_RELTOL, 0.0);
      else
        g_assert_cmpfloat (ncm_matrix_get (cov_diag, i, j), ==, 0.0);

      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_ogk, i, j), ==, ref_cov_ogk[i][j], TEST_RELTOL, 0.0);
    }
  }

  ncm_vector_free (ess);
  ncm_vector_free (pvals);
  ncm_vector_free (pvals5);
  ncm_vector_free (cs);
  ncm_matrix_free (cov_diag);
  ncm_matrix_free (cov_ogk);
}

static void
test_ncm_stats_rows_new (TestNcmStatsRows *test, gconstpointer pdata)
{
  const gdouble rho[TEST_NCOLS] = {0.9, 0.5, 0.0};
  gdouble x[TEST_NCOLS]         = {30.0, 15.0, 6.0};
  NcmRNG *rng                   = ncm_rng_seeded_new (NULL, TEST_SEED);
  guint t, p;

  test->sample = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  test->svec   = ncm_stats_vec_new (TEST_NCOLS, NCM_STATS_VEC_COV, TRUE);
  test->rows   = ncm_stats_rows_new (TEST_NCOLS);

  for (t = 0; t < TEST_NROWS; t++)
  {
    NcmVector *row = ncm_vector_new (TEST_NCOLS);

    for (p = 0; p < TEST_NCOLS; p++)
    {
      x[p] = rho[p] * x[p] + ncm_rng_gaussian_gen (rng, 0.0, 1.0);
      ncm_vector_set (row, p, x[p] + (p + 1.0));
    }

    ncm_stats_vec_append (test->svec, row, TRUE);
    g_ptr_array_add (test->sample, row);
  }

  ncm_stats_rows_set_rows (test->rows, test->sample);

  ncm_rng_free (rng);
}

static void
test_ncm_stats_rows_free (TestNcmStatsRows *test, gconstpointer pdata)
{
  g_ptr_array_unref (test->sample);
  ncm_stats_vec_clear (&test->svec);
  NCM_TEST_FREE (ncm_stats_rows_free, test->rows);
}

static void
test_ncm_stats_rows_reference (TestNcmStatsRows *test, gconstpointer pdata)
{
  _test_ncm_stats_rows_check (test->rows);
}

/* A serialized copy carries the number of columns, and analyzes the rows it is given. */
static void
test_ncm_stats_rows_reference_serialized (TestNcmStatsRows *test, gconstpointer pdata)
{
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmStatsRows *rows = NCM_STATS_ROWS (ncm_serialize_dup_obj (ser, G_OBJECT (test->rows)));

  g_assert_cmpuint (ncm_stats_rows_get_ncols (rows), ==, TEST_NCOLS);
  g_assert_cmpuint (ncm_stats_rows_get_nrows (rows), ==, 0);

  ncm_stats_rows_set_rows (rows, test->sample);
  _test_ncm_stats_rows_check (rows);

  ncm_stats_rows_free (rows);
  ncm_serialize_free (ser);
}

static void
test_ncm_stats_rows_reference_stats_vec (TestNcmStatsRows *test, gconstpointer pdata)
{
  gint ess_bindex, hd_bindex;
  guint ess_wp, ess_wp_order, hd_wp, hd_wp_order;
  gdouble ess_wp_ess, hd_wp_pvalue, cs_mean, cs_var;
  NcmVector *ess      = ncm_stats_vec_max_ess_time (test->svec, 0, &ess_bindex, &ess_wp, &ess_wp_order, &ess_wp_ess);
  NcmVector *pvals    = ncm_stats_vec_heidel_diag (test->svec, 0, 0.0, &hd_bindex, &hd_wp, &hd_wp_order, &hd_wp_pvalue);
  NcmVector *cs       = ncm_stats_vec_visual_heidel_diag (test->svec, 0, 7, &cs_mean, &cs_var);
  NcmMatrix *cov_diag = ncm_stats_vec_compute_cov_robust_diag (test->svec);
  NcmMatrix *cov_ogk  = ncm_stats_vec_compute_cov_robust_ogk (test->svec);
  guint i, j;

  for (i = 0; i < TEST_NCOLS; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (ess, i), ==, ref_ess[i], TEST_RELTOL, 0.0);

  g_assert_cmpint (ess_bindex, ==, ref_ess_bindex);
  g_assert_cmpuint (ess_wp, ==, ref_ess_wp);
  g_assert_cmpuint (ess_wp_order, ==, ref_ess_wp_order);
  ncm_assert_cmpdouble_e (ess_wp_ess, ==, ref_ess[ref_ess_wp], TEST_RELTOL, 0.0);

  for (i = 0; i < TEST_NCOLS; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (pvals, i), ==, ref_pvals[i], TEST_RELTOL, 0.0);

  g_assert_cmpint (hd_bindex, ==, ref_hd_bindex);
  g_assert_cmpuint (hd_wp, ==, ref_hd_wp);
  g_assert_cmpuint (hd_wp_order, ==, ref_hd_wp_order);
  ncm_assert_cmpdouble_e (hd_wp_pvalue, ==, ref_pvals[ref_hd_wp], TEST_RELTOL, 0.0);

  g_assert_cmpuint (ncm_vector_len (cs), ==, TEST_NROWS - 7);
  ncm_assert_cmpdouble_e (cs_mean, ==, ref_visual_mean, TEST_RELTOL, 0.0);
  ncm_assert_cmpdouble_e (cs_var, ==, ref_visual_var, TEST_RELTOL, 0.0);

  for (i = 0; i < TEST_NCOLS; i++)
    g_assert_cmpfloat (ncm_stats_vec_estimate_const_break (test->svec, i), ==, ref_const_break[i]);

  for (i = 0; i < TEST_NCOLS; i++)
  {
    ncm_assert_cmpdouble_e (ncm_matrix_get (cov_diag, i, i), ==, ref_cov_diag[i], TEST_RELTOL, 0.0);

    for (j = 0; j < TEST_NCOLS; j++)
      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_ogk, i, j), ==, ref_cov_ogk[i][j], TEST_RELTOL, 0.0);
  }

  ncm_vector_free (ess);
  ncm_vector_free (pvals);
  ncm_vector_free (cs);
  ncm_matrix_free (cov_diag);
  ncm_matrix_free (cov_ogk);
}

/* A shorter sample gives the results of its own rows: the robust variances of the first
 * 100 rows, which the transient dominates, exceed those of the whole sample. */
static void
test_ncm_stats_rows_set_nrows (TestNcmStatsRows *test, gconstpointer pdata)
{
  NcmStatsRows *rows  = test->rows;
  NcmMatrix *cov_full = ncm_matrix_new (TEST_NCOLS, TEST_NCOLS);
  NcmMatrix *cov_head = ncm_matrix_new (TEST_NCOLS, TEST_NCOLS);
  guint i;

  ncm_stats_rows_cov_robust_diag (rows, cov_full);

  ncm_stats_rows_set_nrows (rows, 100);
  g_assert_cmpuint (ncm_stats_rows_get_nrows (rows), ==, 100);
  g_assert_true (ncm_stats_rows_peek_row (rows, 99) == g_ptr_array_index (test->sample, 99));

  ncm_stats_rows_cov_robust_diag (rows, cov_head);

  for (i = 0; i < 2; i++)
    g_assert_cmpfloat (ncm_matrix_get (cov_head, i, i), >, ncm_matrix_get (cov_full, i, i));

  ncm_matrix_free (cov_full);
  ncm_matrix_free (cov_head);
}

/* The property reads back what the constructor set. */
static void
test_ncm_stats_rows_properties (TestNcmStatsRows *test, gconstpointer pdata)
{
  guint ncols;

  g_object_get (test->rows, "ncols", &ncols, NULL);

  g_assert_cmpuint (ncols, ==, TEST_NCOLS);
}

/* The intended cycle: set_rows, the analyses, reset. The object is given an array only it
 * references, whose element free function counts: reset releases every row, and the next
 * set_rows gives the reference results again. */
static void
test_ncm_stats_rows_reset (TestNcmStatsRows *test, gconstpointer pdata)
{
  GPtrArray *own = g_ptr_array_new_full (TEST_NROWS, &_test_ncm_stats_rows_count_free);
  NcmMatrix *cov = ncm_matrix_new (TEST_NCOLS, TEST_NCOLS);
  guint t;

  for (t = 0; t < TEST_NROWS; t++)
    g_ptr_array_add (own, ncm_vector_ref (g_ptr_array_index (test->sample, t)));

  ncm_stats_rows_set_rows (test->rows, own);
  g_ptr_array_unref (own);

  ncm_stats_rows_cov_robust_ogk (test->rows, cov);
  ncm_matrix_free (cov);

  _test_ncm_stats_rows_nfreed = 0;
  ncm_stats_rows_reset (test->rows);
  g_assert_cmpuint (_test_ncm_stats_rows_nfreed, ==, TEST_NROWS);
  g_assert_cmpuint (ncm_stats_rows_get_nrows (test->rows), ==, 0);
  g_assert_cmpuint (ncm_stats_rows_get_ncols (test->rows), ==, TEST_NCOLS);

  ncm_stats_rows_set_rows (test->rows, test->sample);
  _test_ncm_stats_rows_check (test->rows);
  ncm_stats_rows_reset (test->rows);
}

