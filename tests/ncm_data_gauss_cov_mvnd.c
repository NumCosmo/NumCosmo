/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_gauss_cov_mvnd.c
 *
 *  Sun February 04 15:07:28 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_data_gauss_cov_mvnd.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#include "ncm_data_gauss_cov_mvnd.h"
#include "ncm_model_mvnd.h"

struct _NcmDataGaussCovMVNDPrivate
{
  gint unused;
};

G_DEFINE_TYPE (NcmDataGaussCovMVND, ncm_data_gauss_cov_mvnd, NCM_TYPE_DATA_GAUSS_COV);

static void
ncm_data_gauss_cov_mvnd_init (NcmDataGaussCovMVND *gauss_mvnd)
{
  gauss_mvnd->priv = G_TYPE_INSTANCE_GET_PRIVATE (gauss_mvnd, NCM_TYPE_DATA_GAUSS_COV_MVND, NcmDataGaussCovMVNDPrivate);

}

static void
ncm_data_gauss_cov_mvnd_finalize (GObject *object)
{

  G_OBJECT_CLASS (ncm_data_gauss_cov_mvnd_parent_class)->finalize (object);
}

static void _ncm_data_gauss_cov_mvnd_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);

static void
ncm_data_gauss_cov_mvnd_class_init (NcmDataGaussCovMVNDClass *klass)
{
  GObjectClass* object_class        = G_OBJECT_CLASS (klass);
  NcmDataGaussCovClass* gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcmDataGaussCovMVNDPrivate));

  object_class->finalize = ncm_data_gauss_cov_mvnd_finalize;

  /*data_class->prepare    = &_ncm_data_gauss_cov_test_prepare;*/
  gauss_class->mean_func = &_ncm_data_gauss_cov_mvnd_mean_func;
  gauss_class->cov_func  = NULL;
}

static void 
_ncm_data_gauss_cov_mvnd_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcmModelMVND *model_mvnd = NCM_MODEL_MVND (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
  ncm_model_mvnd_mean (model_mvnd, vp);
}

NcmDataGaussCovMVND *
ncm_data_gauss_cov_mvnd_new (const guint dim)
{
  NcmDataGaussCovMVND *gauss_mvnd = g_object_new (NCM_TYPE_DATA_GAUSS_COV_MVND,
                                                  "n-points", dim,
                                                  "use-norma", TRUE,
                                                  NULL);
  return gauss_mvnd;
}

NcmDataGaussCovMVND *
ncm_data_gauss_cov_mvnd_new_full (const guint dim, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, const gdouble mean_min, const gdouble mean_max, NcmRNG *rng)
{
  NcmDataGaussCovMVND *gauss_mvnd = ncm_data_gauss_cov_mvnd_new (dim);
  ncm_data_gauss_cov_mvnd_gen_cov_mean (gauss_mvnd, sigma_min, sigma_max, cor_level, mean_min, mean_max, rng);

  return gauss_mvnd;
}

NcmDataGaussCovMVND *
ncm_data_gauss_cov_mvnd_ref (NcmDataGaussCovMVND *data_mvnd)
{
  return g_object_ref (data_mvnd);
}

void 
ncm_data_gauss_cov_mvnd_free (NcmDataGaussCovMVND *data_mvnd)
{
  g_object_unref (data_mvnd);
}

void 
ncm_data_gauss_cov_mvnd_clear (NcmDataGaussCovMVND **data_mvnd)
{
  g_clear_object (data_mvnd);
}

void 
ncm_data_gauss_cov_mvnd_gen_cov_mean (NcmDataGaussCovMVND *data_mvnd, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, const gdouble mean_min, const gdouble mean_max, NcmRNG *rng)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
  gint i;

  g_assert_cmpfloat (mean_min, <=, mean_max);

  ncm_matrix_fill_rand_cov (gcov->cov, sigma_min, sigma_max, cor_level, rng);

  if (mean_min == mean_max)
  {
    ncm_vector_set_all (gcov->y, mean_min);
  }
  else
  {
    for (i = 0; i < gcov->np; i++)
    {
      const gdouble mean = ncm_rng_uniform_gen (rng, mean_min, mean_max);
      ncm_vector_set (gcov->y, i, mean);
    }
  }

/*
  g_message ("\n");
  ncm_matrix_log_vals (gcov->cov, "# COV : ", "% 12.5g");
  ncm_vector_log_vals (gcov->y, "# MEAN:", "% 12.5g", TRUE); 
*/  
  ncm_data_set_init (NCM_DATA (gcov), TRUE);
}

NcmVector *
ncm_data_gauss_cov_mvnd_peek_mean (NcmDataGaussCovMVND *data_mvnd)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

  return gcov->y;
}

