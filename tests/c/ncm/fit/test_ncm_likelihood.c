/***************************************************************************
 *            test_ncm_likelihood.c
 *
 *  Wed September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_likelihood.c
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

/* An m2lnL prior returning a constant. */
#define TEST_TYPE_PRIOR_CONST (test_prior_const_get_type ())

G_DECLARE_FINAL_TYPE (TestPriorConst, test_prior_const, TEST, PRIOR_CONST, NcmPrior)

struct _TestPriorConst
{
  NcmPrior parent_instance;
};

G_DEFINE_TYPE (TestPriorConst, test_prior_const, NCM_TYPE_PRIOR)

#define TEST_PRIOR_CONST_VALUE (3.25)

static void
test_prior_const_init (TestPriorConst *p)
{
}

static void
_test_prior_const_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  res[0] = TEST_PRIOR_CONST_VALUE;
}

static void
test_prior_const_class_init (TestPriorConstClass *klass)
{
  NCM_PRIOR_CLASS (klass)->is_m2lnL = TRUE;
  NCM_MSET_FUNC_CLASS (klass)->eval = &_test_prior_const_eval;
}

typedef struct _TestNcmLikelihood
{
  NcmModel *model;
  NcmMSet *mset;
  NcmDataset *dset;
  NcmLikelihood *lh;
  NcmPrior *prior_g;
  NcmPrior *prior_f;
  guint dim;
} TestNcmLikelihood;

static void
test_ncm_likelihood_new (TestNcmLikelihood *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20260930);
  NcmDataGaussCovMVND *data_mvnd = NULL;
  guint i;

  test->dim = 3;
  data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 1.0, 50.0, -1.0, 1.0, rng);

  /* Without the normalization the data -2 ln L is the sum of squares of its least-squares f. */
  ncm_data_gauss_cov_use_norma (NCM_DATA_GAUSS_COV (data_mvnd), FALSE);
  test->model = NCM_MODEL (ncm_model_mvnd_new (test->dim));
  test->mset  = ncm_mset_new (test->model, NULL, NULL);
  test->dset  = ncm_dataset_new_list (data_mvnd, NULL);
  test->lh    = ncm_likelihood_new (test->dset);

  for (i = 0; i < test->dim; i++)
    ncm_model_orig_param_set (test->model, i, 0.1 * (i + 1.0));

  test->prior_g = NCM_PRIOR (ncm_prior_gauss_param_new (test->model, 0, 0.3, 0.2));
  test->prior_f = NCM_PRIOR (ncm_prior_flat_param_new (test->model, 1, -1.0, 0.2005, 1.0e-3));

  ncm_likelihood_priors_add (test->lh, test->prior_g);
  ncm_likelihood_priors_add (test->lh, test->prior_f);

  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
}

static void
test_ncm_likelihood_free (TestNcmLikelihood *test, gconstpointer pdata)
{
  ncm_prior_free (test->prior_g);
  ncm_prior_free (test->prior_f);
  ncm_likelihood_free (test->lh);
  ncm_dataset_free (test->dset);
  ncm_mset_free (test->mset);
  ncm_model_free (test->model);
}

static void
test_ncm_likelihood_m2lnL (TestNcmLikelihood *test, gconstpointer pdata)
{
  const gdouble f_g = ncm_mset_func_eval0 (NCM_MSET_FUNC (test->prior_g), test->mset);
  const gdouble f_f = ncm_mset_func_eval0 (NCM_MSET_FUNC (test->prior_f), test->mset);
  NcmPrior *prior_c = g_object_new (TEST_TYPE_PRIOR_CONST, NULL);
  gdouble m2lnL_data, m2lnL, m2lnL_priors;
  NcmVector *v;

  ncm_likelihood_priors_take (test->lh, prior_c);

  g_assert_cmpuint (ncm_likelihood_priors_length_f (test->lh), ==, 2);
  g_assert_cmpuint (ncm_likelihood_priors_length_m2lnL (test->lh), ==, 1);
  g_assert_true (ncm_likelihood_priors_peek_f (test->lh, 1) == test->prior_f);
  g_assert_true (ncm_likelihood_priors_peek_m2lnL (test->lh, 0) == prior_c);

  /* The flat prior sits a quarter width inside its upper limit: f^2 is not negligible. */
  g_assert_cmpfloat (f_f * f_f, >, 1.0e-3);

  ncm_dataset_m2lnL_val (test->dset, test->mset, &m2lnL_data);
  ncm_likelihood_m2lnL_val (test->lh, test->mset, &m2lnL);
  ncm_likelihood_priors_m2lnL_val (test->lh, test->mset, &m2lnL_priors);

  ncm_assert_cmpdouble_e (m2lnL_priors, ==, f_g * f_g + f_f * f_f + TEST_PRIOR_CONST_VALUE, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (m2lnL, ==, m2lnL_data + m2lnL_priors, 1.0e-14, 0.0);

  /* Terms: the data, then f^2 of the least-squares priors, then the m2lnL priors. */
  v = ncm_likelihood_peek_m2lnL_v (test->lh);
  g_assert_cmpuint (ncm_vector_len (v), ==, 4);
  ncm_assert_cmpdouble_e (ncm_vector_get (v, 0), ==, m2lnL_data, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_vector_get (v, 1), ==, f_g * f_g, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_vector_get (v, 2), ==, f_f * f_f, 1.0e-15, 0.0);
  g_assert_cmpfloat (ncm_vector_get (v, 3), ==, TEST_PRIOR_CONST_VALUE);
}

static void
test_ncm_likelihood_leastsquares (TestNcmLikelihood *test, gconstpointer pdata)
{
  const guint n    = ncm_dataset_get_n (test->dset);
  NcmVector *f     = ncm_vector_new (n + 2);
  NcmVector *f_ref = ncm_vector_new (n);
  gdouble m2lnL, sum_f2 = 0.0;
  guint i;

  ncm_likelihood_leastsquares_f (test->lh, test->mset, f);
  ncm_dataset_leastsquares_f (test->dset, test->mset, f_ref);

  for (i = 0; i < n; i++)
    g_assert_cmpfloat (ncm_vector_get (f, i), ==, ncm_vector_get (f_ref, i));

  g_assert_cmpfloat (ncm_vector_get (f, n), ==, ncm_mset_func_eval0 (NCM_MSET_FUNC (test->prior_g), test->mset));
  g_assert_cmpfloat (ncm_vector_get (f, n + 1), ==, ncm_mset_func_eval0 (NCM_MSET_FUNC (test->prior_f), test->mset));

  /* With least-squares priors only, the posterior is the sum of squares. */
  for (i = 0; i < n + 2; i++)
    sum_f2 += gsl_pow_2 (ncm_vector_get (f, i));

  ncm_likelihood_m2lnL_val (test->lh, test->mset, &m2lnL);
  ncm_assert_cmpdouble_e (m2lnL, ==, sum_f2, 1.0e-12, 0.0);

  ncm_vector_free (f);
  ncm_vector_free (f_ref);
}

static void
test_ncm_likelihood_leastsquares_m2lnL_prior (TestNcmLikelihood *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    NcmVector *f = ncm_vector_new (ncm_dataset_get_n (test->dset) + 2);

    ncm_likelihood_priors_take (test->lh, g_object_new (TEST_TYPE_PRIOR_CONST, NULL));
    ncm_likelihood_leastsquares_f (test->lh, test->mset, f);

    return;
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_likelihood_leastsquares_f: cannot calculate least-squares f*");
}

static void
test_ncm_likelihood_peek_range (TestNcmLikelihood *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    ncm_likelihood_priors_peek_m2lnL (test->lh, 0);

    return;
  }

  /* The test runs with nonfatal assertions, so the child reports and returns. */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_stderr ("*assertion failed (i < lh->priors_m2lnL->len)*");
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/likelihood/m2lnL", TestNcmLikelihood, NULL, &test_ncm_likelihood_new, &test_ncm_likelihood_m2lnL, &test_ncm_likelihood_free);
  g_test_add ("/ncm/likelihood/leastsquares", TestNcmLikelihood, NULL, &test_ncm_likelihood_new, &test_ncm_likelihood_leastsquares, &test_ncm_likelihood_free);
  g_test_add ("/ncm/likelihood/leastsquares/m2lnL_prior", TestNcmLikelihood, NULL, &test_ncm_likelihood_new, &test_ncm_likelihood_leastsquares_m2lnL_prior, &test_ncm_likelihood_free);
  g_test_add ("/ncm/likelihood/peek_range", TestNcmLikelihood, NULL, &test_ncm_likelihood_new, &test_ncm_likelihood_peek_range, &test_ncm_likelihood_free);

  g_test_run ();
}

