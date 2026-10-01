/***************************************************************************
 *            test_ncm_dataset.c
 *
 *  Mon Sep 28 23:00:00 2026
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
 * Data blocks sharing a resource, as the CMB likelihoods share a Boltzmann solver:
 * each block raises the shared requirement in its prepare, and its mean is the
 * requirement in force when it is evaluated. Evaluated after every block has been
 * prepared, every mean is the largest requirement.
 */

typedef struct _TestShared
{
  guint req;
} TestShared;

#define TEST_TYPE_BLOCK (test_block_get_type ())
G_DECLARE_FINAL_TYPE (TestBlock, test_block, TEST, BLOCK, NcmData)

struct _TestBlock
{
  NcmData parent_instance;
  TestShared *shared;
  guint req;
  guint req_at_cov;
};

G_DEFINE_TYPE (TestBlock, test_block, NCM_TYPE_DATA)

static void
test_block_init (TestBlock *tb)
{
  tb->shared     = NULL;
  tb->req        = 0;
  tb->req_at_cov = 0;
}

static guint
_test_block_get_length (NcmData *data)
{
  return 1;
}

static void
_test_block_prepare (NcmData *data, NcmMSet *mset)
{
  TestBlock *tb = TEST_BLOCK (data);

  tb->shared->req = MAX (tb->shared->req, tb->req);
}

static void
_test_block_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  TestBlock *tb = TEST_BLOCK (data);

  ncm_vector_set (mu, 0, tb->shared->req + ncm_mset_fparam_get (mset, 0));
}

static void
_test_block_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  TestBlock *tb = TEST_BLOCK (data);

  *m2lnL = tb->shared->req;
}

static void
_test_block_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H)
{
  TestBlock *tb = TEST_BLOCK (data);

  tb->req_at_cov = tb->shared->req;
}

static void
_test_block_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  TestBlock *tb = TEST_BLOCK (data);

  tb->req_at_cov = tb->shared->req;
}

static void
test_block_class_init (TestBlockClass *klass)
{
  NcmDataClass *data_class = NCM_DATA_CLASS (klass);

  data_class->name        = "Test block";
  data_class->get_length  = &_test_block_get_length;
  data_class->prepare     = &_test_block_prepare;
  data_class->mean_vector = &_test_block_mean_vector;
  data_class->m2lnL_val   = &_test_block_m2lnL_val;
  data_class->inv_cov_UH  = &_test_block_inv_cov_UH;
  data_class->inv_cov_Uf  = &_test_block_inv_cov_Uf;
}

static TestBlock *
test_block_new (TestShared *shared, const guint req)
{
  TestBlock *tb = g_object_new (TEST_TYPE_BLOCK, NULL);

  tb->shared = shared;
  tb->req    = req;
  ncm_data_set_init (NCM_DATA (tb), TRUE);

  return tb;
}

typedef struct _TestNcmDataset
{
  TestShared shared;
  TestBlock *b1;
  TestBlock *b2;
  NcmDataset *dset;
  NcmMSet *mset;
} TestNcmDataset;

void test_ncm_dataset_new (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_free (TestNcmDataset *test, gconstpointer pdata);

void test_ncm_dataset_shared_m2lnL (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_shared_mean_vector (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_shared_fisher (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_shared_fisher_bias (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_fisher_inout (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_fisher_no_fparams (TestNcmDataset *test, gconstpointer pdata);
void test_ncm_dataset_copy_bootstrap (void);
void test_ncm_dataset_bootstrap_total_empty (void);
void test_ncm_dataset_errors (void);
void test_ncm_dataset_fisher_bad_size_subprocess (void);
void test_ncm_dataset_no_realization_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/dataset/shared/m2lnL", TestNcmDataset, NULL, &test_ncm_dataset_new, &test_ncm_dataset_shared_m2lnL, &test_ncm_dataset_free);
  g_test_add ("/ncm/dataset/shared/mean_vector", TestNcmDataset, NULL, &test_ncm_dataset_new, &test_ncm_dataset_shared_mean_vector, &test_ncm_dataset_free);
  g_test_add ("/ncm/dataset/shared/fisher", TestNcmDataset, NULL, &test_ncm_dataset_new, &test_ncm_dataset_shared_fisher, &test_ncm_dataset_free);
  g_test_add ("/ncm/dataset/shared/fisher_bias", TestNcmDataset, NULL, &test_ncm_dataset_new, &test_ncm_dataset_shared_fisher_bias, &test_ncm_dataset_free);
  g_test_add ("/ncm/dataset/fisher/inout", TestNcmDataset, NULL, &test_ncm_dataset_new, &test_ncm_dataset_fisher_inout, &test_ncm_dataset_free);
  g_test_add ("/ncm/dataset/fisher/no_fparams", TestNcmDataset, NULL, &test_ncm_dataset_new, &test_ncm_dataset_fisher_no_fparams, &test_ncm_dataset_free);
  g_test_add_func ("/ncm/dataset/copy/bootstrap", &test_ncm_dataset_copy_bootstrap);
  g_test_add_func ("/ncm/dataset/bootstrap/total_empty", &test_ncm_dataset_bootstrap_total_empty);
  g_test_add_func ("/ncm/dataset/errors", &test_ncm_dataset_errors);
  g_test_add_func ("/ncm/dataset/errors/fisher_bad_size/subprocess", &test_ncm_dataset_fisher_bad_size_subprocess);
  g_test_add_func ("/ncm/dataset/errors/no_realization/subprocess", &test_ncm_dataset_no_realization_subprocess);

  g_test_run ();
}

void
test_ncm_dataset_new (TestNcmDataset *test, gconstpointer pdata)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (1);

  test->shared.req = 0;
  test->b1         = test_block_new (&test->shared, 1);
  test->b2         = test_block_new (&test->shared, 2);
  test->dset       = ncm_dataset_new_list (test->b1, test->b2, NULL);
  test->mset       = ncm_mset_new (mvnd, NULL, NULL);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  ncm_model_mvnd_free (mvnd);
}

void
test_ncm_dataset_free (TestNcmDataset *test, gconstpointer pdata)
{
  ncm_dataset_free (test->dset);
  ncm_data_free (NCM_DATA (test->b1));
  ncm_data_free (NCM_DATA (test->b2));
  ncm_mset_free (test->mset);
}

void
test_ncm_dataset_shared_m2lnL (TestNcmDataset *test, gconstpointer pdata)
{
  gdouble m2lnL;

  ncm_dataset_m2lnL_val (test->dset, test->mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 4.0);
}

void
test_ncm_dataset_shared_mean_vector (TestNcmDataset *test, gconstpointer pdata)
{
  NcmVector *mu = ncm_vector_new (2);

  /* The first call already sees every block's requirement. */
  ncm_dataset_mean_vector (test->dset, test->mset, mu);
  g_assert_cmpfloat (ncm_vector_get (mu, 0), ==, 2.0);
  g_assert_cmpfloat (ncm_vector_get (mu, 1), ==, 2.0);

  ncm_vector_free (mu);
}

void
test_ncm_dataset_shared_fisher (TestNcmDataset *test, gconstpointer pdata)
{
  NcmMatrix *IM = NULL;

  ncm_dataset_fisher_matrix (test->dset, test->mset, &IM);
  g_assert_cmpuint (test->b1->req_at_cov, ==, 2);
  g_assert_cmpuint (test->b2->req_at_cov, ==, 2);

  ncm_matrix_free (IM);
}

void
test_ncm_dataset_shared_fisher_bias (TestNcmDataset *test, gconstpointer pdata)
{
  NcmVector *f_true      = ncm_vector_new (2);
  NcmMatrix *IM          = NULL;
  NcmVector *delta_theta = NULL;

  ncm_vector_set_all (f_true, 2.0);
  ncm_dataset_fisher_matrix_bias (test->dset, test->mset, f_true, &IM, &delta_theta);
  g_assert_cmpuint (test->b1->req_at_cov, ==, 2);
  g_assert_cmpuint (test->b2->req_at_cov, ==, 2);

  ncm_vector_free (f_true);
  ncm_vector_free (delta_theta);
  ncm_matrix_free (IM);
}

void
test_ncm_dataset_fisher_inout (TestNcmDataset *test, gconstpointer pdata)
{
  NcmVector *f_true      = ncm_vector_new (2);
  NcmMatrix *IM          = NULL;
  NcmVector *delta_theta = NULL;
  NcmMatrix *IM_in;
  NcmVector *delta_theta_in;

  /* NULL pointers allocate. */
  ncm_dataset_fisher_matrix (test->dset, test->mset, &IM);
  g_assert_nonnull (IM);
  g_assert_cmpuint (ncm_matrix_nrows (IM), ==, 1);

  /* A matrix passed in is overwritten in place (valgrind: nothing is lost). */
  IM_in = IM;
  ncm_dataset_fisher_matrix (test->dset, test->mset, &IM);
  g_assert_true (IM == IM_in);

  ncm_vector_set_all (f_true, 2.0);
  ncm_dataset_fisher_matrix_bias (test->dset, test->mset, f_true, &IM, &delta_theta);
  g_assert_true (IM == IM_in);
  g_assert_nonnull (delta_theta);

  delta_theta_in = delta_theta;
  ncm_dataset_fisher_matrix_bias (test->dset, test->mset, f_true, &IM, &delta_theta);
  g_assert_true (IM == IM_in);
  g_assert_true (delta_theta == delta_theta_in);

  ncm_vector_free (f_true);
  ncm_vector_free (delta_theta);
  ncm_matrix_free (IM);
}

void
test_ncm_dataset_fisher_no_fparams (TestNcmDataset *test, gconstpointer pdata)
{
  NcmVector *f_true      = ncm_vector_new (2);
  NcmMatrix *IM          = ncm_matrix_new (1, 1);
  NcmVector *delta_theta = ncm_vector_new (1);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
  ncm_vector_set_all (f_true, 2.0);

  /* Without free parameters the objects passed in are freed and cleared. */
  ncm_dataset_fisher_matrix (test->dset, test->mset, &IM);
  g_assert_null (IM);

  IM = ncm_matrix_new (1, 1);
  ncm_dataset_fisher_matrix_bias (test->dset, test->mset, f_true, &IM, &delta_theta);
  g_assert_null (IM);
  g_assert_null (delta_theta);

  ncm_vector_free (f_true);
}

void
test_ncm_dataset_copy_bootstrap (void)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 1);
  NcmDataGaussCovMVND *data      = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-2, 5.0e-2, 20.0, 1.0, 2.0, rng);
  NcmDataset *dset               = ncm_dataset_new_list (data, NULL);
  const NcmDatasetBStrapType t[] = {NCM_DATASET_BSTRAP_PARTIAL, NCM_DATASET_BSTRAP_TOTAL};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (t); i++)
  {
    NcmDataset *copy;
    NcmDatasetBStrapType copy_type;

    ncm_dataset_bootstrap_set (dset, t[i]);
    copy = ncm_dataset_copy (dset);

    g_object_get (copy, "bootstrap-type", &copy_type, NULL);
    g_assert_cmpint (copy_type, ==, t[i]);
    g_assert_true (ncm_dataset_peek_data (copy, 0) == NCM_DATA (data));

    /* The copy resamples like the original. */
    ncm_dataset_bootstrap_resample (copy, rng);

    ncm_dataset_free (copy);
  }

  ncm_dataset_free (dset);
  ncm_data_free (NCM_DATA (data));
  ncm_rng_free (rng);
}

/* A total bootstrap can give a block no draws; that block contributes zero */
void
test_ncm_dataset_bootstrap_total_empty (void)
{
  NcmRNG *rng                = ncm_rng_seeded_new (NULL, 1);
  NcmDataGaussCovMVND *data1 = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-2, 5.0e-2, 20.0, 1.0, 2.0, rng);
  NcmDataGaussCovMVND *data2 = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-2, 5.0e-2, 20.0, 1.0, 2.0, rng);
  NcmDataset *dset           = ncm_dataset_new_list (data1, data2, NULL);
  NcmModelMVND *mvnd         = ncm_model_mvnd_new (2);
  NcmMSet *mset              = ncm_mset_new (mvnd, NULL, NULL);
  guint n_empty              = 0;
  guint i;

  ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_TOTAL);

  for (i = 0; i < 100; i++)
  {
    NcmVector *m2lnL_v = ncm_vector_new (2);
    guint j;

    ncm_dataset_bootstrap_resample (dset, rng);
    ncm_dataset_m2lnL_vec (dset, mset, m2lnL_v);

    for (j = 0; j < 2; j++)
    {
      if (ncm_bootstrap_get_bsize (ncm_data_peek_bootstrap (ncm_dataset_peek_data (dset, j))) == 0)
      {
        g_assert_cmpfloat (ncm_vector_get (m2lnL_v, j), ==, 0.0);
        n_empty++;
      }
    }

    ncm_vector_free (m2lnL_v);
  }

  g_assert_cmpuint (n_empty, >, 0);

  ncm_mset_free (mset);
  ncm_model_mvnd_free (mvnd);
  ncm_dataset_free (dset);
  ncm_data_free (NCM_DATA (data1));
  ncm_data_free (NCM_DATA (data2));
  ncm_rng_free (rng);
}

void
test_ncm_dataset_errors (void)
{
  g_test_trap_subprocess ("/ncm/dataset/errors/no_realization/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the bootstrap has no realization, call ncm_dataset_bootstrap_resample() first*");

  g_test_trap_subprocess ("/ncm/dataset/errors/fisher_bad_size/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the Fisher matrix passed in is 3 x 3, but there are 1 free parameters*");
}

void
test_ncm_dataset_fisher_bad_size_subprocess (void)
{
  TestNcmDataset test;
  NcmMatrix *IM = ncm_matrix_new (3, 3);

  test_ncm_dataset_new (&test, NULL);
  ncm_dataset_fisher_matrix (test.dset, test.mset, &IM);
  ncm_matrix_free (IM);
  test_ncm_dataset_free (&test, NULL);
}

void
test_ncm_dataset_no_realization_subprocess (void)
{
  NcmRNG *rng               = ncm_rng_seeded_new (NULL, 1);
  NcmDataGaussCovMVND *data = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-2, 5.0e-2, 20.0, 1.0, 2.0, rng);
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmModelMVND *mvnd        = ncm_model_mvnd_new (2);
  NcmMSet *mset             = ncm_mset_new (mvnd, NULL, NULL);
  gdouble m2lnL;

  ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_PARTIAL);
  ncm_dataset_m2lnL_val (dset, mset, &m2lnL);
}

