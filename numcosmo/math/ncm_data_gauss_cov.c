/***************************************************************************
 *            ncm_data_gauss_cov.c
 *
 *  Tue November 20 18:46:02 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_data_gauss_cov
 * @title: Gaussian Data - Cov
 * @short_description: Gaussian data object, covariance
 *
 * Generic gaussian distribution which uses the covariance matrix as input.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_gauss_cov.h"
#include "math/ncm_cfg.h"
#include "math/ncm_c.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_USE_DET,
  PROP_MEAN,
  PROP_COV,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataGaussCov, ncm_data_gauss_cov, NCM_TYPE_DATA);

static void
ncm_data_gauss_cov_init (NcmDataGaussCov *gauss)
{
  gauss->np               = 0;
  gauss->y                = NULL;
  gauss->v                = NULL;
  gauss->cov              = NULL;
  gauss->LLT              = NULL;
  gauss->prepared_LLT     = FALSE;
  gauss->use_det          = FALSE;
}

static void
_ncm_data_gauss_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_gauss_cov_parent_class)->constructed (object);
  {
  }
}

static void
_ncm_data_gauss_cov_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS_COV (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_gauss_cov_set_size (gauss, g_value_get_uint (value));
      break;
    case PROP_USE_DET:
      gauss->use_det = g_value_get_boolean (value);
      break;
    case PROP_MEAN:
      ncm_vector_set_from_variant (gauss->y, g_value_get_variant (value));
      break;
    case PROP_COV:
      ncm_matrix_set_from_variant (gauss->cov, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_data_gauss_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS_COV (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, gauss->np);
      break;
    case PROP_USE_DET:
      g_value_set_boolean (value, gauss->use_det);
      break;
    case PROP_MEAN:
      g_value_take_variant (value, ncm_vector_get_variant (gauss->y));
      break;
    case PROP_COV:
      g_value_take_variant (value, ncm_matrix_get_variant (gauss->cov));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_data_gauss_cov_dispose (GObject *object)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (object);

  ncm_vector_clear (&gauss->y);
  ncm_vector_clear (&gauss->v);
  ncm_matrix_clear (&gauss->cov);
  ncm_matrix_clear (&gauss->LLT);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_cov_parent_class)->dispose (object);
}

static void
_ncm_data_gauss_cov_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_cov_parent_class)->finalize (object);
}

static guint _ncm_data_gauss_cov_get_length (NcmData *data); 
/* static void _ncm_data_gauss_cov_begin (NcmData *data); */
static void _ncm_data_gauss_cov_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_gauss_cov_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_cov_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np);
static guint _ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss);

static void
ncm_data_gauss_cov_class_init (NcmDataGaussCovClass *klass)
{
  GObjectClass* object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->constructed  = &_ncm_data_gauss_cov_constructed;
  object_class->set_property = &_ncm_data_gauss_cov_set_property;
  object_class->get_property = &_ncm_data_gauss_cov_get_property;
  object_class->dispose      = &_ncm_data_gauss_cov_dispose;
  object_class->finalize     = &_ncm_data_gauss_cov_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_USE_DET,
                                   g_param_spec_boolean ("use-det",
                                                         NULL,
                                                         "Use determinant to calculate -2lnL",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MEAN,
                                   g_param_spec_variant ("mean",
                                                         NULL,
                                                         "Data mean",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_COV,
                                   g_param_spec_variant ("cov",
                                                         NULL,
                                                         "Data covariance",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap          = TRUE;
  
  data_class->get_length         = &_ncm_data_gauss_cov_get_length;
  data_class->begin              = NULL;

  data_class->resample           = &_ncm_data_gauss_cov_resample;
  data_class->m2lnL_val          = &_ncm_data_gauss_cov_m2lnL_val;
  data_class->leastsquares_f     = &_ncm_data_gauss_cov_leastsquares_f;

  gauss_cov_class->mean_func = NULL;
  gauss_cov_class->cov_func  = NULL;
  gauss_cov_class->set_size  = &_ncm_data_gauss_cov_set_size;
  gauss_cov_class->get_size  = &_ncm_data_gauss_cov_get_size;
}

static guint 
_ncm_data_gauss_cov_get_length (NcmData *data) 
{ 
  return NCM_DATA_GAUSS_COV (data)->np; 
}

static void
_ncm_data_gauss_cov_prepare_LLT (NcmData *data)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (data);

  if (gauss->LLT == NULL)
    gauss->LLT = ncm_matrix_dup (gauss->cov);
  else
    ncm_matrix_memcpy (gauss->LLT, gauss->cov);

  ncm_matrix_cholesky_decomp (gauss->LLT);

  gauss->prepared_LLT = TRUE;
}

static void
_ncm_data_gauss_cov_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gboolean cov_update = FALSE;
  gint ret;
  guint i;

  if (gauss_cov_class->cov_func != NULL)
    cov_update = gauss_cov_class->cov_func (gauss, mset, gauss->cov);

  if (cov_update || !gauss->prepared_LLT)
    _ncm_data_gauss_cov_prepare_LLT (data);

  ncm_rng_lock (rng);
  for (i = 0; i < gauss->np; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    ncm_vector_set (gauss->v, i, u_i);
  }
  ncm_rng_unlock (rng);
  
  ret = gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, 
                        ncm_matrix_gsl (gauss->LLT), ncm_vector_gsl (gauss->v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_resample", ret);

  gauss_cov_class->mean_func (gauss, mset, gauss->y);
  ncm_vector_sub (gauss->y, gauss->v);
}

static void
_ncm_data_gauss_cov_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gboolean cov_update = FALSE;
  gint ret;

  *m2lnL = 0.0;

  gauss_cov_class->mean_func (gauss, mset, gauss->v);

  ncm_vector_sub (gauss->v, gauss->y);
  
  if (gauss_cov_class->cov_func != NULL)
    cov_update = gauss_cov_class->cov_func (gauss, mset, gauss->cov);

  if (cov_update || !gauss->prepared_LLT)
    _ncm_data_gauss_cov_prepare_LLT (data);

  ret = gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, 
                        ncm_matrix_gsl (gauss->LLT), ncm_vector_gsl (gauss->v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_m2lnL_val", ret);

  if (!ncm_data_bootstrap_enabled (data))
  {
    ret = gsl_blas_ddot (ncm_vector_gsl (gauss->v),
                         ncm_vector_gsl (gauss->v),
                         m2lnL);
    NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_m2lnL_val", ret);

    if (gauss->use_det)
    {
      guint i;
      gdouble lndetL = 0.0;
      *m2lnL += gauss->np * ncm_c_ln2pi ();

      for (i = 0; i < gauss->np; i++)
      {
        lndetL += log (ncm_matrix_get (gauss->LLT, i, i));
      }
      *m2lnL += 2.0 * lndetL;
    }
  }
  else
  {
    const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
    guint i;
    g_assert (ncm_bootstrap_is_init (data->bstrap));
    
    if (gauss->use_det)
    {
      gdouble lndetL = 0.0;
      *m2lnL += bsize * ncm_c_ln2pi ();
      for (i = 0; i < bsize; i++)
      {
        guint k = ncm_bootstrap_get (data->bstrap, i);
        const gdouble u_i = ncm_vector_get (gauss->v, k); 
        lndetL += log (ncm_matrix_get (gauss->LLT, k, k));
        *m2lnL += u_i * u_i;
      }
      *m2lnL += 2.0 * lndetL;
    }
    else
    {
      for (i = 0; i < bsize; i++)
      {
        guint k = ncm_bootstrap_get (data->bstrap, i);
        const gdouble u_i = ncm_vector_get (gauss->v, k); 
        *m2lnL += u_i * u_i;
      }
    }
  }
}

static void
_ncm_data_gauss_cov_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gboolean cov_update = FALSE;
  gint ret;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("NcmDataGaussCov: does not support bootstrap with least squares");
  
  gauss_cov_class->mean_func (gauss, mset, v);
  ncm_vector_sub (v, gauss->y);
  
  if (gauss_cov_class->cov_func != NULL)
    cov_update = gauss_cov_class->cov_func (gauss, mset, gauss->cov);

  if (cov_update || !gauss->prepared_LLT)
    _ncm_data_gauss_cov_prepare_LLT (data);

  ret = gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, 
                        ncm_matrix_gsl (gauss->LLT), ncm_vector_gsl (v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_leastsquares_f", ret);
}

static void 
_ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcmData *data = NCM_DATA (gauss);
  if ((np == 0) || (np != gauss->np))
  {
    gauss->np = 0;
    ncm_vector_clear (&gauss->y);
    ncm_vector_clear (&gauss->v);
    ncm_matrix_clear (&gauss->cov);
    ncm_matrix_clear (&gauss->LLT);
    data->init = FALSE;
  }
  if ((np != 0) && (np != gauss->np))
  {
    NcmData *data = NCM_DATA (gauss);
    gauss->np  = np;
    gauss->y   = ncm_vector_new (gauss->np);
    gauss->v   = ncm_vector_new (gauss->np);
    gauss->cov = ncm_matrix_new (gauss->np, gauss->np);
    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (data->bstrap, np);
      ncm_bootstrap_set_bsize (data->bstrap, np);
    }
    data->init = FALSE;
  }
}

static guint 
_ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss)
{
  return gauss->np;
}

/**
 * ncm_data_gauss_cov_set_size: (virtual set_size)
 * @gauss: a #NcmDataGaussCov
 * @np: data size.
 *
 * Sets the data size to @np.
 * 
 */
void 
ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np)
{
  NCM_DATA_GAUSS_COV_GET_CLASS (gauss)->set_size (gauss, np);
}

/**
 * ncm_data_gauss_cov_get_size: (virtual get_size)
 * @gauss: a #NcmDataGaussCov
 *
 * Gets the data size.
 * 
 * Returns: Data size.
 * 
 */
guint 
ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss)
{
  return NCM_DATA_GAUSS_COV_GET_CLASS (gauss)->get_size (gauss);
}
