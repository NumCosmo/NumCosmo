/***************************************************************************
 *            test_ncm_data.c
 *
 *  Mon Sep 28 22:00:00 2026
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
 * A minimal NcmData: its mean is the NcmModelMVND mean, its covariance the identity.
 * It counts begin calls and records the free parameters at which the covariance
 * methods are called.
 */

#define TEST_DATA_DIM 2

#define TEST_TYPE_DATA (test_data_get_type ())
G_DECLARE_DERIVABLE_TYPE (TestData, test_data, TEST, DATA, NcmData)

struct _TestDataClass
{
  NcmDataClass parent_class;
};

typedef struct _TestDataPrivate
{
  guint begin_count;
  gdouble y[TEST_DATA_DIM];
  gdouble cov_point[TEST_DATA_DIM];
} TestDataPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (TestData, test_data, NCM_TYPE_DATA)

static void
test_data_init (TestData *td)
{
  TestDataPrivate * const self = test_data_get_instance_private (td);
  guint i;

  self->begin_count = 0;

  for (i = 0; i < TEST_DATA_DIM; i++)
  {
    self->y[i]         = 0.0;
    self->cov_point[i] = GSL_NAN;
  }
}

static guint
_test_data_get_length (NcmData *data)
{
  return TEST_DATA_DIM;
}

static void
_test_data_begin (NcmData *data)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (data));

  self->begin_count++;
}

static void
_test_data_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  ncm_model_mvnd_mean (NCM_MODEL_MVND (ncm_mset_peek (mset, ncm_model_mvnd_id ())), mu);
}

static void
_test_data_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (data));
  NcmVector *mu                = ncm_vector_new (TEST_DATA_DIM);
  guint i;

  _test_data_mean_vector (data, mset, mu);

  for (i = 0; i < TEST_DATA_DIM; i++)
    self->y[i] = ncm_vector_get (mu, i) + 1.0;

  ncm_vector_free (mu);
}

static void
_test_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (data));
  NcmVector *mu                = ncm_vector_new (TEST_DATA_DIM);
  guint i;

  _test_data_mean_vector (data, mset, mu);

  *m2lnL = 0.0;

  for (i = 0; i < TEST_DATA_DIM; i++)
    *m2lnL += gsl_pow_2 (self->y[i] - ncm_vector_get (mu, i));

  ncm_vector_free (mu);
}

static void
_test_data_record_point (NcmData *data, NcmMSet *mset)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (data));
  guint i;

  for (i = 0; i < TEST_DATA_DIM; i++)
    self->cov_point[i] = ncm_mset_fparam_get (mset, i);
}

static void
_test_data_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H)
{
  _test_data_record_point (data, mset);
}

static void
_test_data_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  _test_data_record_point (data, mset);
}

static void
test_data_class_init (TestDataClass *klass)
{
  NcmDataClass *data_class = NCM_DATA_CLASS (klass);

  data_class->name        = "Test data";
  data_class->bootstrap   = TRUE;
  data_class->get_length  = &_test_data_get_length;
  data_class->begin       = &_test_data_begin;
  data_class->resample    = &_test_data_resample;
  data_class->m2lnL_val   = &_test_data_m2lnL_val;
  data_class->mean_vector = &_test_data_mean_vector;
  data_class->inv_cov_UH  = &_test_data_inv_cov_UH;
  data_class->inv_cov_Uf  = &_test_data_inv_cov_Uf;
}

/* The same data without resample. */

#define TEST_TYPE_DATA_NO_RESAMPLE (test_data_no_resample_get_type ())
G_DECLARE_FINAL_TYPE (TestDataNoResample, test_data_no_resample, TEST, DATA_NO_RESAMPLE, TestData)

struct _TestDataNoResample
{
  TestData parent_instance;
};

G_DEFINE_TYPE (TestDataNoResample, test_data_no_resample, TEST_TYPE_DATA)

static void
test_data_no_resample_init (TestDataNoResample *td)
{
}

static void
test_data_no_resample_class_init (TestDataNoResampleClass *klass)
{
  NCM_DATA_CLASS (klass)->resample = NULL;
}

typedef struct _TestNcmData
{
  NcmMSet *mset;
  NcmData *data;
  NcmRNG *rng;
  gdouble p[TEST_DATA_DIM];
} TestNcmData;

void test_ncm_data_new (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_free (TestNcmData *test, gconstpointer pdata);

void test_ncm_data_begin (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_fisher (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_fisher_bias (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_fisher_no_fparams (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_bootstrap (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_desc (TestNcmData *test, gconstpointer pdata);
void test_ncm_data_errors (void);
void test_ncm_data_prepare_uninit_subprocess (void);
void test_ncm_data_resample_missing_subprocess (void);
void test_ncm_data_fisher_bad_size_subprocess (void);
void test_ncm_data_fisher_bias_bad_f_true_subprocess (void);
void test_ncm_data_bootstrap_uninit_subprocess (void);
void test_ncm_data_bootstrap_set_null_subprocess (void);
void test_ncm_data_bootstrap_set_bad_size_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/data/begin", TestNcmData, NULL, &test_ncm_data_new, &test_ncm_data_begin, &test_ncm_data_free);
  g_test_add ("/ncm/data/fisher", TestNcmData, NULL, &test_ncm_data_new, &test_ncm_data_fisher, &test_ncm_data_free);
  g_test_add ("/ncm/data/fisher_bias", TestNcmData, NULL, &test_ncm_data_new, &test_ncm_data_fisher_bias, &test_ncm_data_free);
  g_test_add ("/ncm/data/fisher_no_fparams", TestNcmData, NULL, &test_ncm_data_new, &test_ncm_data_fisher_no_fparams, &test_ncm_data_free);
  g_test_add ("/ncm/data/bootstrap", TestNcmData, NULL, &test_ncm_data_new, &test_ncm_data_bootstrap, &test_ncm_data_free);
  g_test_add ("/ncm/data/desc", TestNcmData, NULL, &test_ncm_data_new, &test_ncm_data_desc, &test_ncm_data_free);

  g_test_add_func ("/ncm/data/errors", &test_ncm_data_errors);
  g_test_add_func ("/ncm/data/errors/prepare_uninit/subprocess", &test_ncm_data_prepare_uninit_subprocess);
  g_test_add_func ("/ncm/data/errors/resample_missing/subprocess", &test_ncm_data_resample_missing_subprocess);
  g_test_add_func ("/ncm/data/errors/fisher_bad_size/subprocess", &test_ncm_data_fisher_bad_size_subprocess);
  g_test_add_func ("/ncm/data/errors/fisher_bias_bad_f_true/subprocess", &test_ncm_data_fisher_bias_bad_f_true_subprocess);
  g_test_add_func ("/ncm/data/errors/bootstrap_uninit/subprocess", &test_ncm_data_bootstrap_uninit_subprocess);
  g_test_add_func ("/ncm/data/errors/bootstrap_set_null/subprocess", &test_ncm_data_bootstrap_set_null_subprocess);
  g_test_add_func ("/ncm/data/errors/bootstrap_set_bad_size/subprocess", &test_ncm_data_bootstrap_set_bad_size_subprocess);

  g_test_run ();
}

void
test_ncm_data_new (TestNcmData *test, gconstpointer pdata)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (TEST_DATA_DIM);
  guint i;

  test->mset = ncm_mset_new (mvnd, NULL, NULL);
  test->data = g_object_new (TEST_TYPE_DATA, NULL);
  test->rng  = ncm_rng_seeded_new (NULL, 1);
  test->p[0] = 0.75;
  test->p[1] = -1.25;

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  for (i = 0; i < TEST_DATA_DIM; i++)
    ncm_mset_fparam_set (test->mset, i, test->p[i]);

  ncm_model_mvnd_free (mvnd);
}

void
test_ncm_data_free (TestNcmData *test, gconstpointer pdata)
{
  ncm_data_free (test->data);
  ncm_mset_free (test->mset);
  ncm_rng_free (test->rng);
}

void
test_ncm_data_begin (TestNcmData *test, gconstpointer pdata)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (test->data));
  gdouble m2lnL;

  g_assert_false (ncm_data_is_init (test->data));

  /* resample prepares, so begin runs once on the data about to change. */
  ncm_data_resample (test->data, test->mset, test->rng);
  g_assert_true (ncm_data_is_init (test->data));
  g_assert_cmpuint (self->begin_count, ==, 1);

  /* The data changed, so begin runs again before the next evaluation, once. */
  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  g_assert_cmpuint (self->begin_count, ==, 2);
  g_assert_cmpfloat (m2lnL, ==, 2.0);

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  g_assert_cmpuint (self->begin_count, ==, 2);

  /* A resample of initialized data also runs begin again afterwards. */
  ncm_data_resample (test->data, test->mset, test->rng);
  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  g_assert_cmpuint (self->begin_count, ==, 3);

  g_assert_false (ncm_data_is_resampling (test->data));
}

void
test_ncm_data_fisher (TestNcmData *test, gconstpointer pdata)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (test->data));
  NcmMatrix *IM                = NULL;
  NcmMatrix *IM_in;
  guint i;

  ncm_data_set_init (test->data, TRUE);

  /* A NULL *IM allocates. */
  ncm_data_fisher_matrix (test->data, test->mset, &IM);
  g_assert_nonnull (IM);
  g_assert_cmpuint (ncm_matrix_nrows (IM), ==, TEST_DATA_DIM);
  g_assert_cmpuint (ncm_matrix_ncols (IM), ==, TEST_DATA_DIM);

  /* The covariance is evaluated at the point, and the point is restored. */
  for (i = 0; i < TEST_DATA_DIM; i++)
  {
    g_assert_cmpfloat (self->cov_point[i], ==, test->p[i]);
    g_assert_cmpfloat (ncm_mset_fparam_get (test->mset, i), ==, test->p[i]);
  }

  /* A matrix passed in is overwritten in place. */
  IM_in = IM;
  ncm_matrix_set_all (IM, -7.0);
  ncm_data_fisher_matrix (test->data, test->mset, &IM);
  g_assert_true (IM == IM_in);
  /* Symmetric up to the summation order of the BLAS product. */
  ncm_assert_cmpdouble_e (ncm_matrix_get (IM, 0, 1), ==, ncm_matrix_get (IM, 1, 0), 1.0e-14, 0.0);
  g_assert_cmpfloat (ncm_matrix_get (IM, 0, 0), >, 0.0);

  ncm_matrix_free (IM);
}

void
test_ncm_data_fisher_bias (TestNcmData *test, gconstpointer pdata)
{
  TestDataPrivate * const self = test_data_get_instance_private (TEST_DATA (test->data));
  NcmVector *f_true            = ncm_vector_new (TEST_DATA_DIM);
  NcmMatrix *IM                = NULL;
  NcmVector *delta_theta       = NULL;
  NcmVector *delta_theta_in;
  guint i;

  ncm_data_set_init (test->data, TRUE);

  /* With the true mean equal to the mean at the point the shift is zero. */
  ncm_data_mean_vector (test->data, test->mset, f_true);
  ncm_data_fisher_matrix_bias (test->data, test->mset, f_true, &IM, &delta_theta);

  for (i = 0; i < TEST_DATA_DIM; i++)
  {
    g_assert_cmpfloat (ncm_vector_get (delta_theta, i), ==, 0.0);
    g_assert_cmpfloat (self->cov_point[i], ==, test->p[i]);
    g_assert_cmpfloat (ncm_mset_fparam_get (test->mset, i), ==, test->p[i]);
  }

  /* A vector passed in is overwritten in place. */
  delta_theta_in = delta_theta;
  ncm_vector_set_all (delta_theta, -7.0);
  ncm_data_fisher_matrix_bias (test->data, test->mset, f_true, &IM, &delta_theta);
  g_assert_true (delta_theta == delta_theta_in);
  g_assert_cmpfloat (ncm_vector_get (delta_theta, 0), ==, 0.0);

  ncm_vector_free (f_true);
  ncm_vector_free (delta_theta);
  ncm_matrix_free (IM);
}

void
test_ncm_data_fisher_no_fparams (TestNcmData *test, gconstpointer pdata)
{
  NcmVector *f_true      = ncm_vector_new (TEST_DATA_DIM);
  NcmMatrix *IM          = ncm_matrix_new (TEST_DATA_DIM, TEST_DATA_DIM);
  NcmVector *delta_theta = ncm_vector_new (TEST_DATA_DIM);

  ncm_data_set_init (test->data, TRUE);
  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
  ncm_vector_set_zero (f_true);

  /* Without free parameters the matrix passed in is freed (valgrind) and cleared. */
  ncm_data_fisher_matrix (test->data, test->mset, &IM);
  g_assert_null (IM);

  IM = ncm_matrix_new (TEST_DATA_DIM, TEST_DATA_DIM);
  ncm_data_fisher_matrix_bias (test->data, test->mset, f_true, &IM, &delta_theta);
  g_assert_null (IM);
  g_assert_null (delta_theta);

  ncm_vector_free (f_true);
}

void
test_ncm_data_bootstrap (TestNcmData *test, gconstpointer pdata)
{
  NcmBootstrap *bstrap;

  ncm_data_set_init (test->data, TRUE);
  g_assert_false (ncm_data_bootstrap_enabled (test->data));

  ncm_data_bootstrap_create (test->data);
  g_assert_true (ncm_data_bootstrap_enabled (test->data));
  g_assert_cmpuint (ncm_bootstrap_get_fsize (ncm_data_peek_bootstrap (test->data)), ==, TEST_DATA_DIM);

  bstrap = ncm_bootstrap_sized_new (TEST_DATA_DIM);
  ncm_data_bootstrap_set (test->data, bstrap);
  g_assert_true (ncm_data_peek_bootstrap (test->data) == bstrap);
  ncm_bootstrap_free (bstrap);

  ncm_data_bootstrap_resample (test->data, test->rng);

  ncm_data_bootstrap_remove (test->data);
  g_assert_false (ncm_data_bootstrap_enabled (test->data));
}

void
test_ncm_data_desc (TestNcmData *test, gconstpointer pdata)
{
  gchar *desc;

  /* Without a description the class name is used. */
  g_assert_cmpstr (ncm_data_peek_desc (test->data), ==, "Test data");

  ncm_data_set_desc (test->data, "Set description");
  g_assert_cmpstr (ncm_data_peek_desc (test->data), ==, "Set description");

  ncm_data_take_desc (test->data, g_strdup ("Taken description"));
  desc = ncm_data_get_desc (test->data);
  g_assert_cmpstr (desc, ==, "Taken description");
  g_free (desc);

  g_object_set (test->data, "long-desc", "A longer description", NULL);
  g_object_get (test->data, "long-desc", &desc, NULL);
  g_assert_cmpstr (desc, ==, "A longer description");
  g_free (desc);
}

void
test_ncm_data_errors (void)
{
  g_test_trap_subprocess ("/ncm/data/errors/prepare_uninit/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*data `Test data' is not initialized; set its data or resample it first*");

  g_test_trap_subprocess ("/ncm/data/errors/resample_missing/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*The data (Test data) does not implement resample*");

  g_test_trap_subprocess ("/ncm/data/errors/fisher_bad_size/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the Fisher matrix passed in is 3 x 3, but there are 2 free parameters*");

  g_test_trap_subprocess ("/ncm/data/errors/fisher_bias_bad_f_true/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*data `Test data' has 2 points, but f_true has 3*");

  g_test_trap_subprocess ("/ncm/data/errors/bootstrap_uninit/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_data_bootstrap_create: data `Test data' is not initialized*");

  g_test_trap_subprocess ("/ncm/data/errors/bootstrap_set_null/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*use ncm_data_bootstrap_remove() to remove the bootstrap*");

  g_test_trap_subprocess ("/ncm/data/errors/bootstrap_set_bad_size/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*data `Test data' has 2 points, but the bootstrap has a full size of 5*");
}

void
test_ncm_data_prepare_uninit_subprocess (void)
{
  TestNcmData test;

  test_ncm_data_new (&test, NULL);
  ncm_data_prepare (test.data, test.mset);
  test_ncm_data_free (&test, NULL);
}

void
test_ncm_data_resample_missing_subprocess (void)
{
  TestNcmData test;

  test_ncm_data_new (&test, NULL);
  ncm_data_free (test.data);
  test.data = g_object_new (TEST_TYPE_DATA_NO_RESAMPLE, NULL);
  ncm_data_resample (test.data, test.mset, test.rng);
  test_ncm_data_free (&test, NULL);
}

void
test_ncm_data_fisher_bad_size_subprocess (void)
{
  TestNcmData test;
  NcmMatrix *IM = ncm_matrix_new (3, 3);

  test_ncm_data_new (&test, NULL);
  ncm_data_set_init (test.data, TRUE);
  ncm_data_fisher_matrix (test.data, test.mset, &IM);
  ncm_matrix_free (IM);
  test_ncm_data_free (&test, NULL);
}

void
test_ncm_data_fisher_bias_bad_f_true_subprocess (void)
{
  TestNcmData test;
  NcmVector *f_true      = ncm_vector_new (3);
  NcmMatrix *IM          = NULL;
  NcmVector *delta_theta = NULL;

  test_ncm_data_new (&test, NULL);
  ncm_data_set_init (test.data, TRUE);
  ncm_data_fisher_matrix_bias (test.data, test.mset, f_true, &IM, &delta_theta);
  ncm_vector_free (f_true);
  test_ncm_data_free (&test, NULL);
}

void
test_ncm_data_bootstrap_uninit_subprocess (void)
{
  TestNcmData test;

  test_ncm_data_new (&test, NULL);
  ncm_data_bootstrap_create (test.data);
  test_ncm_data_free (&test, NULL);
}

void
test_ncm_data_bootstrap_set_null_subprocess (void)
{
  TestNcmData test;

  test_ncm_data_new (&test, NULL);
  ncm_data_set_init (test.data, TRUE);
  ncm_data_bootstrap_set (test.data, NULL);
  test_ncm_data_free (&test, NULL);
}

void
test_ncm_data_bootstrap_set_bad_size_subprocess (void)
{
  TestNcmData test;
  NcmBootstrap *bstrap = ncm_bootstrap_sized_new (5);

  test_ncm_data_new (&test, NULL);
  ncm_data_set_init (test.data, TRUE);
  ncm_data_bootstrap_set (test.data, bstrap);
  ncm_bootstrap_free (bstrap);
  test_ncm_data_free (&test, NULL);
}

