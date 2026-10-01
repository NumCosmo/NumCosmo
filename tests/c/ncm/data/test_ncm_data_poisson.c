/***************************************************************************
 *            test_ncm_data_poisson.c
 *
 *  Tue Sep 29 02:00:00 2026
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/*
 * A toy NcmDataPoisson: the mean of bin i is A times its width, with A the
 * NcmModelMVND mean mu_0, so the model reads the bin edges.
 */

#define TEST_TYPE_POISSON (test_poisson_get_type ())
G_DECLARE_FINAL_TYPE (TestPoisson, test_poisson, TEST, POISSON, NcmDataPoisson)

struct _TestPoisson
{
  NcmDataPoisson parent_instance;
};

G_DEFINE_TYPE (TestPoisson, test_poisson, NCM_TYPE_DATA_POISSON)

static void
test_poisson_init (TestPoisson *tp)
{
}

static gdouble
_test_poisson_mean_func (NcmDataPoisson *poisson, NcmMSet *mset, guint n)
{
  const gdouble A = ncm_model_orig_vparam_get (ncm_mset_peek (mset, ncm_model_mvnd_id ()), NCM_MODEL_MVND_MEAN, 0);
  gdouble lower, upper;

  ncm_data_poisson_get_bin_range (poisson, n, &lower, &upper);

  return A * (upper - lower);
}

static void
test_poisson_class_init (TestPoissonClass *klass)
{
  NCM_DATA_CLASS (klass)->name              = "Test Poisson data";
  NCM_DATA_POISSON_CLASS (klass)->mean_func = &_test_poisson_mean_func;
}

typedef struct _TestNcmDataPoisson
{
  NcmMSet *mset;
  NcmData *data;
  NcmDataPoisson *poisson;
} TestNcmDataPoisson;

void test_ncm_data_poisson_new (TestNcmDataPoisson *test, gconstpointer pdata);
void test_ncm_data_poisson_free (TestNcmDataPoisson *test, gconstpointer pdata);

void test_ncm_data_poisson_edges (TestNcmDataPoisson *test, gconstpointer pdata);
void test_ncm_data_poisson_binning (TestNcmDataPoisson *test, gconstpointer pdata);
void test_ncm_data_poisson_m2lnL (TestNcmDataPoisson *test, gconstpointer pdata);
void test_ncm_data_poisson_serialize (TestNcmDataPoisson *test, gconstpointer pdata);
void test_ncm_data_poisson_errors (TestNcmDataPoisson *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/data_poisson/edges", TestNcmDataPoisson, NULL, &test_ncm_data_poisson_new, &test_ncm_data_poisson_edges, &test_ncm_data_poisson_free);
  g_test_add ("/ncm/data_poisson/binning", TestNcmDataPoisson, NULL, &test_ncm_data_poisson_new, &test_ncm_data_poisson_binning, &test_ncm_data_poisson_free);
  g_test_add ("/ncm/data_poisson/m2lnL", TestNcmDataPoisson, NULL, &test_ncm_data_poisson_new, &test_ncm_data_poisson_m2lnL, &test_ncm_data_poisson_free);
  g_test_add ("/ncm/data_poisson/serialize", TestNcmDataPoisson, NULL, &test_ncm_data_poisson_new, &test_ncm_data_poisson_serialize, &test_ncm_data_poisson_free);
  g_test_add ("/ncm/data_poisson/errors", TestNcmDataPoisson, NULL, &test_ncm_data_poisson_new, &test_ncm_data_poisson_errors, &test_ncm_data_poisson_free);

  g_test_run ();
}

/* Bins [0, 1], [1, 3], [3, 6], of widths 1, 2 and 3, with counts 1, 2 and 0. */
static void
_test_init_three_bins (TestNcmDataPoisson *test)
{
  const gdouble edges_data[]  = {0.0, 1.0, 3.0, 6.0};
  const gdouble counts_data[] = {1.0, 2.0, 0.0};
  NcmVector *edges            = ncm_vector_new_data_dup ((gdouble *) edges_data, 4, 1);
  NcmVector *counts           = ncm_vector_new_data_dup ((gdouble *) counts_data, 3, 1);

  ncm_data_poisson_init_from_vector (test->poisson, edges, counts);

  ncm_vector_free (edges);
  ncm_vector_free (counts);
}

void
test_ncm_data_poisson_new (TestNcmDataPoisson *test, gconstpointer pdata)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (1);

  test->mset    = ncm_mset_new (mvnd, NULL, NULL);
  test->poisson = g_object_new (TEST_TYPE_POISSON, "n-bins", 3, NULL);
  test->data    = NCM_DATA (test->poisson);

  ncm_model_orig_vparam_set (NCM_MODEL (mvnd), NCM_MODEL_MVND_MEAN, 0, 1.0);
  ncm_model_mvnd_free (mvnd);
}

void
test_ncm_data_poisson_free (TestNcmDataPoisson *test, gconstpointer pdata)
{
  ncm_data_free (test->data);
  ncm_mset_free (test->mset);
}

void
test_ncm_data_poisson_edges (TestNcmDataPoisson *test, gconstpointer pdata)
{
  NcmVector *edges = ncm_data_poisson_get_bin_edges (test->poisson);
  gdouble lower, upper;
  guint i;

  /* A new size has the edges 0, 1, ..., n. */
  g_assert_cmpuint (ncm_vector_len (edges), ==, 4);

  for (i = 0; i < 4; i++)
    g_assert_cmpfloat (ncm_vector_get (edges, i), ==, i);

  ncm_vector_free (edges);

  _test_init_three_bins (test);

  ncm_data_poisson_get_bin_range (test->poisson, 1, &lower, &upper);
  g_assert_cmpfloat (lower, ==, 1.0);
  g_assert_cmpfloat (upper, ==, 3.0);

  g_object_get (test->poisson, "bin-edges", &edges, NULL);
  g_assert_cmpfloat (ncm_vector_get (edges, 3), ==, 6.0);
  ncm_vector_free (edges);
}

void
test_ncm_data_poisson_binning (TestNcmDataPoisson *test, gconstpointer pdata)
{
  const gdouble edges_data[] = {0.0, 1.0, 3.0, 6.0};
  const gdouble x_data[]     = {0.5, 1.5, 2.5, 2.9, 5.0, 7.0, -1.0};
  NcmVector *edges           = ncm_vector_new_data_dup ((gdouble *) edges_data, 4, 1);
  NcmVector *x               = ncm_vector_new_data_dup ((gdouble *) x_data, 7, 1);
  NcmVector *counts;

  ncm_data_poisson_init_from_binning (test->poisson, edges, x);
  counts = ncm_data_poisson_get_hist_vals (test->poisson);

  /* Values outside the edges are ignored. */
  g_assert_cmpfloat (ncm_vector_get (counts, 0), ==, 1.0);
  g_assert_cmpfloat (ncm_vector_get (counts, 1), ==, 3.0);
  g_assert_cmpfloat (ncm_vector_get (counts, 2), ==, 1.0);
  g_assert_cmpfloat (ncm_data_poisson_get_sum (test->poisson), ==, 5.0);
  g_assert_true (ncm_data_is_init (test->data));

  ncm_data_poisson_init_zero (test->poisson, edges);
  g_assert_cmpfloat (ncm_data_poisson_get_sum (test->poisson), ==, 0.0);

  ncm_vector_free (counts);
  ncm_vector_free (edges);
  ncm_vector_free (x);
}

void
test_ncm_data_poisson_m2lnL (TestNcmDataPoisson *test, gconstpointer pdata)
{
  NcmVector *f = ncm_vector_new (3);
  NcmVector *means;
  gdouble m2lnL;
  guint i;

  _test_init_three_bins (test);

  /* The means are the widths 1, 2, 3; the first two bins hold their mean and add
   * nothing to the deviance, the empty bin adds 2 lambda = 6. */
  means = ncm_data_poisson_get_hist_means (test->poisson, test->mset);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (ncm_vector_get (means, i), ==, i + 1.0);

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 6.0);

  /* The least-squares vector holds the square roots of the terms. */
  ncm_data_leastsquares_f (test->data, test->mset, f);
  g_assert_cmpfloat (ncm_vector_get (f, 0), ==, 0.0);
  g_assert_cmpfloat (ncm_vector_get (f, 1), ==, 0.0);
  g_assert_cmpfloat (ncm_vector_get (f, 2), ==, sqrt (6.0));

  /* The Fisher whitening divides by sqrt (lambda). */
  ncm_vector_set_all (f, 1.0);
  ncm_data_inv_cov_Uf (test->data, test->mset, f);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (ncm_vector_get (f, i), ==, 1.0 / sqrt (i + 1.0));

  ncm_vector_free (means);
  ncm_vector_free (f);
}

void
test_ncm_data_poisson_serialize (TestNcmDataPoisson *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmDataPoisson *dup;
  NcmVector *edges;
  gdouble m2lnL, m2lnL_dup;

  _test_init_three_bins (test);

  dup   = NCM_DATA_POISSON (ncm_serialize_dup_obj (ser, G_OBJECT (test->poisson)));
  edges = ncm_data_poisson_get_bin_edges (dup);

  /* The edges survive, so the copy's means and likelihood are the same. */
  g_assert_cmpfloat (ncm_vector_get (edges, 1), ==, 1.0);
  g_assert_cmpfloat (ncm_vector_get (edges, 2), ==, 3.0);
  g_assert_cmpfloat (ncm_vector_get (edges, 3), ==, 6.0);

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  ncm_data_set_init (NCM_DATA (dup), TRUE);
  ncm_data_m2lnL_val (NCM_DATA (dup), test->mset, &m2lnL_dup);
  g_assert_cmpfloat (m2lnL_dup, ==, m2lnL);

  ncm_vector_free (edges);
  ncm_data_free (NCM_DATA (dup));
  ncm_serialize_free (ser);
}

void
test_ncm_data_poisson_errors (TestNcmDataPoisson *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    const gchar *which = g_getenv ("TEST_NCM_DATA_POISSON_ERROR");
    NcmVector *v2      = ncm_vector_new (2);
    NcmVector *v4      = ncm_vector_new (4);

    if (g_strcmp0 (which, "lengths") == 0)
    {
      ncm_data_poisson_init_from_vector (test->poisson, v2, v2);
    }
    else if (g_strcmp0 (which, "order") == 0)
    {
      ncm_vector_set_all (v4, 1.0);
      ncm_data_poisson_init_zero (test->poisson, v4);
    }
    else if (g_strcmp0 (which, "empty") == 0)
    {
      NcmVector *v1 = ncm_vector_new (1);

      ncm_data_poisson_init_zero (test->poisson, v1);
    }
    else if (g_strcmp0 (which, "counts") == 0)
    {
      g_object_set (test->poisson, "mean", v2, NULL);
    }
    else
    {
      gdouble lower, upper;

      ncm_data_poisson_get_bin_range (test->poisson, 3, &lower, &upper);
    }

    return; /* LCOV_EXCL_LINE */
  }

  {
    const gchar *cases[][2] = {
      {"lengths", "*2 bin edges for 2 counts, expected 3*"},
      {"order",   "*the bin edges must be increasing, but edge 0 is 1 and edge 1 is 1*"},
      {"empty",   "*at least two bin edges are needed, but 1 were given*"},
      {"counts",  "*has 3 bins, but the counts have 2 components*"},
      {"range",   "*has 3 bins, bin 3 requested*"},
    };
    guint i;

    for (i = 0; i < G_N_ELEMENTS (cases); i++)
    {
      g_setenv ("TEST_NCM_DATA_POISSON_ERROR", cases[i][0], TRUE);
      g_test_trap_subprocess (NULL, 0, 0);
      g_test_trap_assert_failed ();
      g_test_trap_assert_stderr (cases[i][1]);
    }

    g_unsetenv ("TEST_NCM_DATA_POISSON_ERROR");
  }
}

