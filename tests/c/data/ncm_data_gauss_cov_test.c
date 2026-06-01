/***************************************************************************
 *            ncm_data_gauss_cov_test.c
 *
 *  Wed June 19 14:25:23 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isotfware.com.br>
 ****************************************************************************/
/*
 * ncm_data_gauss_cov_test.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isotfware.com.br>
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
#include "ncm_data_gauss_cov_test.h"

G_DEFINE_TYPE (NcmDataGaussCovTest, ncm_data_gauss_cov_test, NCM_TYPE_DATA_GAUSS_COV);

static void
ncm_data_gauss_cov_test_init (NcmDataGaussCovTest *gcov_test)
{
  gcov_test->a = 0.0;
  gcov_test->b = 0.0;
  gcov_test->c = 0.0;
  gcov_test->d = 0.0;
}

static void
ncm_data_gauss_cov_test_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_cov_test_parent_class)->finalize (object);
}

static void _ncm_data_gauss_cov_test_prepare (NcmData *data, NcmMSet *mset);
static gboolean _ncm_data_gauss_cov_test_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov);
static gboolean _ncm_data_gauss_cov_test_cov_func0 (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov_pass);

static void
ncm_data_gauss_cov_test_class_init (NcmDataGaussCovTestClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->finalize = &ncm_data_gauss_cov_test_finalize;

  data_class->prepare    = &_ncm_data_gauss_cov_test_prepare;
  gauss_class->mean_func = &ncm_data_gauss_cov_test_mean_func;
  gauss_class->cov_func  = &_ncm_data_gauss_cov_test_cov_func0;
}

static void
_ncm_data_gauss_cov_test_prepare (NcmData *data, NcmMSet *mset)
{
}

void
ncm_data_gauss_cov_test_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcmDataGaussCovTest *gcov_test = NCM_DATA_GAUSS_COV_TEST (gauss);
  guint np                       = ncm_data_gauss_cov_get_size (gauss);
  guint i;

  for (i = 0; i < np; i++)
  {
    gdouble x = 1.0 / (np - 1.0) * i;

    ncm_vector_set (vp, i, gcov_test->a + gcov_test->b * cos (gcov_test->c * x + gcov_test->d));
  }
}

static gboolean
_ncm_data_gauss_cov_test_cov_func0 (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov_pass)
{
  NcmMatrix *cov = ncm_data_gauss_cov_peek_cov (gauss);

  if (cov != cov_pass)
    ncm_matrix_memcpy (cov_pass, cov);

  return TRUE;
}

static gboolean
_ncm_data_gauss_cov_test_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov_pass)
{
  NcmMatrix *cov = ncm_data_gauss_cov_peek_cov (gauss);
  NcmVector *y   = ncm_data_gauss_cov_peek_mean (gauss);
  guint np       = ncm_data_gauss_cov_get_size (gauss);
  guint i;

  g_assert_true (np == ncm_matrix_nrows (cov));
  g_assert_true (np == ncm_matrix_ncols (cov));

  for (i = 0; i < np; i++)
  {
    gdouble var_i = g_test_rand_double_range (1e-3, 1.0e-2) * fabs (ncm_vector_get (y, i)) * 1.0e-3;

    ncm_matrix_set (cov, i, i, var_i);
  }

  for (i = 0; i < np; i++)
  {
    guint j;

    for (j = i + 1; j < np; j++)
    {
      gdouble cor_ij = g_test_rand_double_range (1.0e-1, 1.0) * (g_test_rand_double_range (0.0, 1.0) > 0.5 ? -1.0 : 1.0);
      gdouble cov_ij = cor_ij * sqrt (ncm_matrix_get (cov, i, i) * ncm_matrix_get (cov, j, j));

      ncm_matrix_set (cov, i, j, cov_ij);
      ncm_matrix_set (cov, j, i, cov_ij);
    }
  }

  {
    NcmMatrix *t1 = ncm_matrix_dup (cov);
    gint ret      = gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, ncm_matrix_gsl (t1), ncm_matrix_gsl (t1), 0.0, ncm_matrix_gsl (cov));

    NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_test_cov_func", ret);
    ncm_matrix_copy_triangle (cov, 'U');
    ncm_matrix_free (t1);
  }

  return FALSE;
}

#define _TEST_NCM_DATA_GAUSS_COV_MIN_SIZE 10
#define _TEST_NCM_DATA_GAUSS_COV_MAX_SIZE 20

void
ncm_data_gauss_cov_test_gen_cov (NcmDataGaussCovTest *gcov_test)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (gcov_test);
  guint size             = g_test_rand_int_range (_TEST_NCM_DATA_GAUSS_COV_MIN_SIZE, _TEST_NCM_DATA_GAUSS_COV_MAX_SIZE);
  NcmVector *y;

  ncm_data_gauss_cov_set_size (gauss, size);

  y = ncm_data_gauss_cov_peek_mean (gauss);

  gcov_test->a = g_test_rand_double_range (1.0, 3.0);
  gcov_test->b = g_test_rand_double_range (1.0, 2.0) * 1e-1;
  gcov_test->c = g_test_rand_double_range (-1.0, 1.0) * 2.0 * M_PI;
  gcov_test->d = g_test_rand_double_range (-1.0, 1.0) * 2.0 * M_PI;

  ncm_data_gauss_cov_test_mean_func (gauss, NULL, y);
  _ncm_data_gauss_cov_test_cov_func (gauss, NULL, NULL);
  ncm_data_set_init (NCM_DATA (gauss), TRUE);
}

NcmData *
ncm_data_gauss_cov_test_new ()
{
  return g_object_new (NCM_TYPE_DATA_GAUSS_COV_TEST, NULL);
}

