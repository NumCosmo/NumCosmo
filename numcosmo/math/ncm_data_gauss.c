/***************************************************************************
 *            ncm_data_gauss.c
 *
 *  Fri Mar 19 14:57:35 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

/**
 * SECTION:ncm_data_gauss
 * @title: Gaussian Data - InvCov
 * @short_description: Gaussian data object, inverse covariance
 *
 * Gaussian distribution which uses the inverse covariance matrix as input.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_gauss.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_MEAN,
  PROP_INV_COV,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataGauss, ncm_data_gauss, NCM_TYPE_DATA);

static void
ncm_data_gauss_init (NcmDataGauss *gauss)
{
  gauss->np           = 0;
  gauss->y            = NULL;
  gauss->v            = NULL;
  gauss->inv_cov      = NULL;
  gauss->LLT          = NULL;
  gauss->prepared_LLT = FALSE;
}

static void
_ncm_data_gauss_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_gauss_parent_class)->constructed (object);
}

static void
ncm_data_gauss_dispose (GObject *object)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (object);

  ncm_vector_clear (&gauss->y);
  ncm_vector_clear (&gauss->v);
  ncm_matrix_clear (&gauss->inv_cov);
  ncm_matrix_clear (&gauss->LLT);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_parent_class)->dispose (object);
}

static void
ncm_data_gauss_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_parent_class)->finalize (object);
}

static void
ncm_data_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_gauss_set_size (gauss, g_value_get_uint (value));
      break;
    case PROP_MEAN:
      ncm_vector_set_from_variant (gauss->y, g_value_get_variant (value));
      break;
    case PROP_INV_COV:
      ncm_matrix_set_from_variant (gauss->inv_cov, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, gauss->np);
      break;
    case PROP_MEAN:
      g_value_take_variant (value, ncm_vector_get_variant (gauss->y));
      break;
    case PROP_INV_COV:
      g_value_take_variant (value, ncm_matrix_get_variant (gauss->inv_cov));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static guint _ncm_data_gauss_get_length (NcmData *data); 
/* static void _ncm_data_gauss_begin (NcmData *data); */
static void _ncm_data_gauss_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_gauss_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_gauss_set_size (NcmDataGauss *gauss, guint np);
static guint _ncm_data_gauss_get_size (NcmDataGauss *gauss);

static void
ncm_data_gauss_class_init (NcmDataGaussClass *klass)
{
  GObjectClass* object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_CLASS (klass);

  object_class->constructed  = &_ncm_data_gauss_constructed;
  object_class->set_property = &ncm_data_gauss_set_property;
  object_class->get_property = &ncm_data_gauss_get_property;
  object_class->dispose      = &ncm_data_gauss_dispose;
  object_class->finalize     = &ncm_data_gauss_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MEAN,
                                   g_param_spec_variant ("mean",
                                                         NULL,
                                                         "Data mean",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_INV_COV,
                                   g_param_spec_variant ("inv-cov",
                                                         NULL,
                                                         "Data covariance inverse",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap      = TRUE;
  data_class->get_length     = &_ncm_data_gauss_get_length;
  data_class->begin          = NULL;

  data_class->resample       = &_ncm_data_gauss_resample;
  data_class->m2lnL_val      = &_ncm_data_gauss_m2lnL_val;
  data_class->leastsquares_f = &_ncm_data_gauss_leastsquares_f;

  gauss_class->mean_func     = NULL;
  gauss_class->inv_cov_func  = NULL;
  gauss_class->set_size      = &_ncm_data_gauss_set_size;
  gauss_class->get_size      = &_ncm_data_gauss_get_size;
}

static guint 
_ncm_data_gauss_get_length (NcmData *data) 
{ 
  return NCM_DATA_GAUSS (data)->np; 
}

static void
_ncm_data_gauss_prepare_LLT (NcmData *data)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (data);

  if (gauss->LLT == NULL)
    gauss->LLT = ncm_matrix_dup (gauss->inv_cov);
  else
    ncm_matrix_memcpy (gauss->LLT, gauss->inv_cov);

  ncm_matrix_cholesky_decomp (gauss->LLT);

  gauss->prepared_LLT = TRUE;
}

static void
_ncm_data_gauss_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (data);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_GET_CLASS (gauss);
  gboolean inv_cov_update = FALSE;
  gint ret;
  guint i;

  if (gauss_class->inv_cov_func != NULL)
    inv_cov_update = gauss_class->inv_cov_func (gauss, mset, gauss->inv_cov);

  if (inv_cov_update || !gauss->prepared_LLT)
    _ncm_data_gauss_prepare_LLT (data);

  ncm_rng_lock (rng);
  for (i = 0; i < gauss->np; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    ncm_vector_set (gauss->v, i, u_i);
  }
  ncm_rng_unlock (rng);
  
  ret = gsl_blas_dtrsv (CblasLower, CblasTrans, CblasNonUnit, 
                        ncm_matrix_gsl (gauss->LLT), ncm_vector_gsl (gauss->v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_resample", ret);
  
  gauss_class->mean_func (gauss, mset, gauss->y);

  ncm_vector_sub (gauss->y, gauss->v);
}

static void
_ncm_data_gauss_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (data);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_GET_CLASS (gauss);
  gboolean inv_cov_update = FALSE;
  guint i, j;

  *m2lnL = 0.0;

  gauss_class->mean_func (gauss, mset, gauss->v);
  ncm_vector_sub (gauss->v, gauss->y);

  if (gauss_class->inv_cov_func != NULL)
    inv_cov_update = gauss_class->inv_cov_func (gauss, mset, gauss->inv_cov);

  if (!ncm_data_bootstrap_enabled (data))
  {

    for (i = 0; i < gauss->np; i++)
    {
      const gdouble f_i = ncm_vector_get (gauss->v, i);
      gdouble u_i = 0.0;
      for (j = 0; j < gauss->np; j++)
        u_i += ncm_matrix_get (gauss->inv_cov, i, j) * ncm_vector_get (gauss->v, j);
      *m2lnL += u_i * f_i;
    }
  }
  else
  {
    const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
    gint ret;

    if (inv_cov_update || !gauss->prepared_LLT)
      _ncm_data_gauss_prepare_LLT (data);

    ret = gsl_blas_dtrmv (CblasLower, CblasTrans, CblasNonUnit, 
                          ncm_matrix_gsl (gauss->LLT), ncm_vector_gsl (gauss->v));
    NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_resample", ret);
    
    for (i = 0; i < bsize; i++)
    {
      guint k = ncm_bootstrap_get (data->bstrap, i);
      gdouble u_i = ncm_vector_get (gauss->v, k);
      *m2lnL += u_i * u_i;
    }
  }
}

static void
_ncm_data_gauss_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (data);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_GET_CLASS (gauss);
  gboolean inv_cov_update = FALSE;
  gint ret;

  gauss_class->mean_func (gauss, mset, v);
  ncm_vector_sub (v, gauss->y);

  if (ncm_data_bootstrap_enabled (data))
    g_error ("NcmDataGauss: does not support bootstrap with least squares");
  
  if (gauss_class->inv_cov_func != NULL)
    inv_cov_update = gauss_class->inv_cov_func (gauss, mset, gauss->inv_cov);

  if (inv_cov_update || !gauss->prepared_LLT)
    _ncm_data_gauss_prepare_LLT (data);

  ret = gsl_blas_dtrmv (CblasLower, CblasTrans, CblasNonUnit, 
                        ncm_matrix_gsl (gauss->LLT), ncm_vector_gsl (v));
  NCM_TEST_GSL_RESULT("_ncm_data_gauss_leastsquares_f", ret);
}

void 
_ncm_data_gauss_set_size (NcmDataGauss *gauss, guint np)
{
  NcmData *data = NCM_DATA (gauss);
  if ((np == 0) || (np != gauss->np))
  {
    gauss->np = 0;
    ncm_vector_clear (&gauss->y);
    ncm_vector_clear (&gauss->v);
    ncm_matrix_clear (&gauss->inv_cov);
    ncm_matrix_clear (&gauss->LLT);
    data->init = FALSE;
  }
  if ((np != 0) && (np != gauss->np))
  {
    gauss->np      = np;
    gauss->y       = ncm_vector_new (gauss->np);
    gauss->v       = ncm_vector_new (gauss->np);
    gauss->inv_cov = ncm_matrix_new (gauss->np, gauss->np);
    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (data->bstrap, np);
      ncm_bootstrap_set_bsize (data->bstrap, np);
    }
    data->init = FALSE;
  }
}

static guint 
_ncm_data_gauss_get_size (NcmDataGauss *gauss)
{
  return gauss->np;
}

/**
 * ncm_data_gauss_set_size: (virtual set_size)
 * @gauss: a #NcmDataGauss
 * @np: data size.
 *
 * Sets the data size to @np.
 * 
 */
void 
ncm_data_gauss_set_size (NcmDataGauss *gauss, guint np)
{
  NCM_DATA_GAUSS_GET_CLASS (gauss)->set_size (gauss, np);
}

/**
 * ncm_data_gauss_get_size: (virtual get_size)
 * @gauss: a #NcmDataGauss
 *
 * Gets the data size.
 * 
 * Returns: Data size.
 * 
 */
guint 
ncm_data_gauss_get_size (NcmDataGauss *gauss)
{
  return NCM_DATA_GAUSS_GET_CLASS (gauss)->get_size (gauss);
}
