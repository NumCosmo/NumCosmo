/***************************************************************************
 *            test_ncm_data_gauss.c
 *
 *  Tue Sep 29 00:00:00 2026
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
 * A NcmDataGauss with a fixed inverse covariance: the mean is mu_1 in every
 * component of the NcmModelMVND mean (mu_0, mu_1).
 */

#define TEST_TYPE_GAUSS (test_gauss_get_type ())
G_DECLARE_DERIVABLE_TYPE (TestGauss, test_gauss, TEST, GAUSS, NcmDataGauss)

struct _TestGaussClass
{
  NcmDataGaussClass parent_class;
};

G_DEFINE_TYPE (TestGauss, test_gauss, NCM_TYPE_DATA_GAUSS)

static void
test_gauss_init (TestGauss *tg)
{
}

static void
_test_gauss_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcmModel *mvnd = ncm_mset_peek (mset, ncm_model_mvnd_id ());

  ncm_vector_set_all (vp, ncm_model_orig_vparam_get (mvnd, NCM_MODEL_MVND_MEAN, 1));
}

static void
test_gauss_class_init (TestGaussClass *klass)
{
  NCM_DATA_CLASS (klass)->name            = "Test Gaussian data";
  NCM_DATA_GAUSS_CLASS (klass)->mean_func = &_test_gauss_mean_func;
}

/*
 * The same data with the inverse covariance I / s^2, s = mu_0; inv_cov_func returns
 * TRUE only when s changed.
 */

#define TEST_TYPE_GAUSS_COV (test_gauss_cov_get_type ())
G_DECLARE_FINAL_TYPE (TestGaussCov, test_gauss_cov, TEST, GAUSS_COV, TestGauss)

struct _TestGaussCov
{
  TestGauss parent_instance;
  gdouble last_s;
};

G_DEFINE_TYPE (TestGaussCov, test_gauss_cov, TEST_TYPE_GAUSS)

static void
test_gauss_cov_init (TestGaussCov *tgc)
{
  tgc->last_s = GSL_NAN;
}

static gboolean
_test_gauss_cov_inv_cov_func (NcmDataGauss *gauss, NcmMSet *mset, NcmMatrix *inv_cov)
{
  TestGaussCov *tgc = TEST_GAUSS_COV (gauss);
  const gdouble s   = ncm_model_orig_vparam_get (ncm_mset_peek (mset, ncm_model_mvnd_id ()), NCM_MODEL_MVND_MEAN, 0);

  if (s == tgc->last_s)
    return FALSE;

  tgc->last_s = s;
  ncm_matrix_set_identity (inv_cov);
  ncm_matrix_scale (inv_cov, 1.0 / (s * s));

  return TRUE;
}

static void
test_gauss_cov_class_init (TestGaussCovClass *klass)
{
  NCM_DATA_GAUSS_CLASS (klass)->inv_cov_func = &_test_gauss_cov_inv_cov_func;
}

typedef struct _TestNcmDataGauss
{
  NcmMSet *mset;
  NcmDataGauss *gauss;
  NcmRNG *rng;
} TestNcmDataGauss;

static void
test_ncm_data_gauss_set_s (TestNcmDataGauss *test, const gdouble s)
{
  ncm_model_orig_vparam_set (ncm_mset_peek (test->mset, ncm_model_mvnd_id ()), NCM_MODEL_MVND_MEAN, 0, s);
}

void test_ncm_data_gauss_new (TestNcmDataGauss *test, gconstpointer pdata);
void test_ncm_data_gauss_free (TestNcmDataGauss *test, gconstpointer pdata);

void test_ncm_data_gauss_inv_cov_prop (TestNcmDataGauss *test, gconstpointer pdata);
void test_ncm_data_gauss_inv_cov_func (TestNcmDataGauss *test, gconstpointer pdata);
void test_ncm_data_gauss_resize (TestNcmDataGauss *test, gconstpointer pdata);
void test_ncm_data_gauss_errors (void);
void test_ncm_data_gauss_not_posdef_subprocess (void);
void test_ncm_data_gauss_ls_bootstrap_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/data_gauss/inv_cov/prop", TestNcmDataGauss, NULL, &test_ncm_data_gauss_new, &test_ncm_data_gauss_inv_cov_prop, &test_ncm_data_gauss_free);
  g_test_add ("/ncm/data_gauss/inv_cov/func", TestNcmDataGauss, NULL, &test_ncm_data_gauss_new, &test_ncm_data_gauss_inv_cov_func, &test_ncm_data_gauss_free);
  g_test_add ("/ncm/data_gauss/resize", TestNcmDataGauss, NULL, &test_ncm_data_gauss_new, &test_ncm_data_gauss_resize, &test_ncm_data_gauss_free);
  g_test_add_func ("/ncm/data_gauss/errors", &test_ncm_data_gauss_errors);
  g_test_add_func ("/ncm/data_gauss/errors/not_posdef/subprocess", &test_ncm_data_gauss_not_posdef_subprocess);
  g_test_add_func ("/ncm/data_gauss/errors/ls_bootstrap/subprocess", &test_ncm_data_gauss_ls_bootstrap_subprocess);

  g_test_run ();
}

void
test_ncm_data_gauss_new (TestNcmDataGauss *test, gconstpointer pdata)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (2);

  test->mset  = ncm_mset_new (mvnd, NULL, NULL);
  test->gauss = g_object_new (TEST_TYPE_GAUSS, "n-points", 2, NULL);
  test->rng   = ncm_rng_seeded_new (NULL, 1);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  ncm_matrix_set_identity (ncm_data_gauss_peek_inv_cov (test->gauss));
  ncm_vector_set (ncm_data_gauss_peek_mean (test->gauss), 0, 1.0);
  ncm_vector_set (ncm_data_gauss_peek_mean (test->gauss), 1, 2.0);
  ncm_data_set_init (NCM_DATA (test->gauss), TRUE);

  ncm_model_mvnd_free (mvnd);
}

void
test_ncm_data_gauss_free (TestNcmDataGauss *test, gconstpointer pdata)
{
  ncm_data_free (NCM_DATA (test->gauss));
  ncm_mset_free (test->mset);
  ncm_rng_free (test->rng);
}

void
test_ncm_data_gauss_inv_cov_prop (TestNcmDataGauss *test, gconstpointer pdata)
{
  NcmData *data  = NCM_DATA (test->gauss);
  NcmVector *f   = ncm_vector_new (2);
  NcmMatrix *inv = ncm_matrix_new (2, 2);
  gdouble m2lnL;

  /* The least-squares vector builds the Cholesky factor of the identity. */
  ncm_data_leastsquares_f (data, test->mset, f);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 5.0);
  g_assert_cmpfloat (ncm_vector_dot (f, f), ==, 5.0);

  /* Setting inv-cov invalidates the factor. */
  ncm_matrix_set_identity (inv);
  ncm_matrix_scale (inv, 4.0);
  g_object_set (test->gauss, "inv-cov", inv, NULL);

  ncm_data_leastsquares_f (data, test->mset, f);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 20.0);
  g_assert_cmpfloat (ncm_vector_dot (f, f), ==, 20.0);

  ncm_vector_free (f);
  ncm_matrix_free (inv);
}

void
test_ncm_data_gauss_inv_cov_func (TestNcmDataGauss *test, gconstpointer pdata)
{
  NcmData *data;
  NcmVector *f = ncm_vector_new (2);
  gdouble m2lnL;

  ncm_data_free (NCM_DATA (test->gauss));
  test->gauss = g_object_new (TEST_TYPE_GAUSS_COV, "n-points", 2, NULL);
  data        = NCM_DATA (test->gauss);
  ncm_vector_set_all (ncm_data_gauss_peek_mean (test->gauss), 1.0);
  ncm_data_set_init (data, TRUE);

  /* The factor is built at s = 1. */
  test_ncm_data_gauss_set_s (test, 1.0);
  ncm_data_leastsquares_f (data, test->mset, f);
  g_assert_cmpfloat (ncm_vector_dot (f, f), ==, 2.0);

  /* m2lnL sees the change to s = 2; the paths that use the factor must see it too. */
  test_ncm_data_gauss_set_s (test, 2.0);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 0.5);

  ncm_data_leastsquares_f (data, test->mset, f);
  g_assert_cmpfloat (ncm_vector_dot (f, f), ==, 0.5);

  ncm_vector_set_all (f, 1.0);
  ncm_data_inv_cov_Uf (data, test->mset, f);
  g_assert_cmpfloat (ncm_vector_get (f, 0), ==, 0.5);
  g_assert_cmpfloat (ncm_vector_get (f, 1), ==, 0.5);

  ncm_vector_free (f);
}

void
test_ncm_data_gauss_resize (TestNcmDataGauss *test, gconstpointer pdata)
{
  NcmData *data = NCM_DATA (test->gauss);
  gdouble m2lnL;

  ncm_data_resample (data, test->mset, test->rng);

  /* A new size frees the factor; resampling rebuilds it. */
  ncm_data_gauss_set_size (test->gauss, 3);
  g_assert_false (ncm_data_is_init (data));
  g_assert_cmpuint (ncm_data_get_length (data), ==, 3);

  ncm_matrix_set_identity (ncm_data_gauss_peek_inv_cov (test->gauss));
  ncm_data_resample (data, test->mset, test->rng);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));
}

void
test_ncm_data_gauss_errors (void)
{
  g_test_trap_subprocess ("/ncm/data_gauss/errors/not_posdef/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*data `Test Gaussian data': the inverse covariance is not positive definite*");

  g_test_trap_subprocess ("/ncm/data_gauss/errors/ls_bootstrap/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*bootstrap is not supported with least squares*");
}

void
test_ncm_data_gauss_not_posdef_subprocess (void)
{
  TestNcmDataGauss test;
  NcmVector *f = ncm_vector_new (2);

  test_ncm_data_gauss_new (&test, NULL);
  ncm_matrix_set (ncm_data_gauss_peek_inv_cov (test.gauss), 0, 0, -1.0);
  ncm_data_leastsquares_f (NCM_DATA (test.gauss), test.mset, f);

  ncm_vector_free (f);
  test_ncm_data_gauss_free (&test, NULL);
}

void
test_ncm_data_gauss_ls_bootstrap_subprocess (void)
{
  TestNcmDataGauss test;
  NcmVector *f = ncm_vector_new (2);

  test_ncm_data_gauss_new (&test, NULL);
  ncm_data_bootstrap_create (NCM_DATA (test.gauss));
  ncm_data_leastsquares_f (NCM_DATA (test.gauss), test.mset, f);

  ncm_vector_free (f);
  test_ncm_data_gauss_free (&test, NULL);
}

