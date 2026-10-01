/***************************************************************************
 *            ncm_data_gauss.c
 *
 *  Fri Mar 19 14:57:35 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmDataGauss:
 *
 * Gaussian data given by its inverse covariance.
 *
 * Abstract class for data $y$ of #NcmDataGauss:n-points points drawn from a normal
 * distribution with mean $\mu$ and covariance $C$, given by its inverse $C^{-1}$
 * (#NcmDataGauss:inv-cov); $-2\ln L = (\mu - y)^T C^{-1} (\mu - y)$, without the
 * normalization. The Cholesky factor $C^{-1} = U^T U$ is computed when needed for
 * resampling, least squares, bootstrap and the Fisher matrix.
 *
 * Subclasses implement the mean_func virtual method, which fills the mean $\mu$
 * given the models in a #NcmMSet. When $C^{-1}$ depends on the models, they also
 * implement inv_cov_func, which updates the matrix passed to it and returns %TRUE
 * when it changed it; the Cholesky factor is then recomputed before its next use.
 * Setting #NcmDataGauss:inv-cov also invalidates the factor, while changing the
 * matrix returned by ncm_data_gauss_peek_inv_cov() in place after the first
 * evaluation does not: set the property instead.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/data/ncm_data_gauss.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/algebra/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_MEAN,
  PROP_INV_COV,
  PROP_SIZE,
};

typedef struct _NcmDataGaussPrivate
{
  guint np;
  NcmVector *y;
  NcmVector *v;
  NcmMatrix *inv_cov;
  NcmMatrix *LLT;
  gboolean prepared_LLT;
} NcmDataGaussPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmDataGauss, ncm_data_gauss, NCM_TYPE_DATA)

static void
ncm_data_gauss_init (NcmDataGauss *gauss)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  self->np           = 0;
  self->y            = NULL;
  self->v            = NULL;
  self->inv_cov      = NULL;
  self->LLT          = NULL;
  self->prepared_LLT = FALSE;
}

static void
_ncm_data_gauss_constructed (GObject *object)
{
  NcmDataGauss *gauss            = NCM_DATA_GAUSS (object);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_GET_CLASS (gauss);

  if (gauss_class->mean_func == NULL)
    g_error ("NcmDataGauss: mean_func not set, class %s cannot be used.", G_OBJECT_CLASS_NAME (gauss_class));

  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_gauss_parent_class)->constructed (object);
}

static void
ncm_data_gauss_dispose (GObject *object)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (object);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  ncm_vector_clear (&self->y);
  ncm_vector_clear (&self->v);
  ncm_matrix_clear (&self->inv_cov);
  ncm_matrix_clear (&self->LLT);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_parent_class)->dispose (object);
}

static void
ncm_data_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (object);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  g_return_if_fail (NCM_IS_DATA_GAUSS (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_gauss_set_size (gauss, g_value_get_uint (value));
      break;
    case PROP_MEAN:
      ncm_vector_substitute (&self->y, g_value_get_object (value), TRUE);
      break;
    case PROP_INV_COV:
      ncm_matrix_substitute (&self->inv_cov, g_value_get_object (value), TRUE);
      self->prepared_LLT = FALSE;
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_data_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (object);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  g_return_if_fail (NCM_IS_DATA_GAUSS (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, self->np);
      break;
    case PROP_MEAN:
      g_value_set_object (value, self->y);
      break;
    case PROP_INV_COV:
      g_value_set_object (value, self->inv_cov);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static guint _ncm_data_gauss_get_length (NcmData *data);

/* static void _ncm_data_gauss_begin (NcmData *data); */
static void _ncm_data_gauss_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_gauss_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_gauss_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu);
static void _ncm_data_gauss_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H);
static void _ncm_data_gauss_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f);

static void _ncm_data_gauss_set_size (NcmDataGauss *gauss, guint np);
static guint _ncm_data_gauss_get_size (NcmDataGauss *gauss);

static void
ncm_data_gauss_class_init (NcmDataGaussClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_CLASS (klass);

  object_class->constructed  = &_ncm_data_gauss_constructed;
  object_class->set_property = &ncm_data_gauss_set_property;
  object_class->get_property = &ncm_data_gauss_get_property;
  object_class->dispose      = &ncm_data_gauss_dispose;

  /**
   * NcmDataGauss:n-points:
   *
   * The number of data points; changing it reallocates the data and marks it not
   * initialized.
   */
  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmDataGauss:mean:
   *
   * The data vector $y$.
   */
  g_object_class_install_property (object_class,
                                   PROP_MEAN,
                                   g_param_spec_object ("mean",
                                                        NULL,
                                                        "Data mean",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmDataGauss:inv-cov:
   *
   * The inverse covariance $C^{-1}$; setting it invalidates its Cholesky factor.
   */
  g_object_class_install_property (object_class,
                                   PROP_INV_COV,
                                   g_param_spec_object ("inv-cov",
                                                        NULL,
                                                        "Data covariance inverse",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = TRUE;
  data_class->get_length = &_ncm_data_gauss_get_length;
  data_class->begin      = NULL;

  data_class->resample       = &_ncm_data_gauss_resample;
  data_class->m2lnL_val      = &_ncm_data_gauss_m2lnL_val;
  data_class->leastsquares_f = &_ncm_data_gauss_leastsquares_f;

  data_class->mean_vector = &_ncm_data_gauss_mean_vector;
  data_class->inv_cov_UH  = &_ncm_data_gauss_inv_cov_UH;
  data_class->inv_cov_Uf  = &_ncm_data_gauss_inv_cov_Uf;

  gauss_class->mean_func    = NULL;
  gauss_class->inv_cov_func = NULL;
  gauss_class->set_size     = &_ncm_data_gauss_set_size;
  gauss_class->get_size     = &_ncm_data_gauss_get_size;
}

static guint
_ncm_data_gauss_get_length (NcmData *data)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (data);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  return self->np;
}

/*
 * Calls inv_cov_func, if any; when it changed the inverse covariance, the Cholesky
 * factor is marked stale. Every path that may update the matrix goes through here, so
 * no path consumes the change without the factor noticing.
 */
static void
_ncm_data_gauss_update_inv_cov (NcmDataGauss *gauss, NcmMSet *mset)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);
  NcmDataGaussClass *gauss_class   = NCM_DATA_GAUSS_GET_CLASS (gauss);

  if ((gauss_class->inv_cov_func != NULL) && gauss_class->inv_cov_func (gauss, mset, self->inv_cov))
    self->prepared_LLT = FALSE;
}

/* The upper Cholesky factor U of C^{-1} = U^T U, recomputed when stale. */
static NcmMatrix *
_ncm_data_gauss_peek_LLT (NcmDataGauss *gauss)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  if (!self->prepared_LLT)
  {
    gint ret;

    if (self->LLT == NULL)
      self->LLT = ncm_matrix_dup (self->inv_cov);
    else
      ncm_matrix_memcpy (self->LLT, self->inv_cov);

    ret = ncm_matrix_cholesky_decomp (self->LLT, 'U');

    if (ret != 0)
      g_error ("_ncm_data_gauss_peek_LLT: data `%s': the inverse covariance is not positive "
               "definite (Cholesky decomposition returned %d).",
               ncm_data_peek_desc (NCM_DATA (gauss)), ret);

    self->prepared_LLT = TRUE;
  }

  return self->LLT;
}

static void
_ncm_data_gauss_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (data);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);
  NcmDataGaussClass *gauss_class   = NCM_DATA_GAUSS_GET_CLASS (gauss);
  NcmMatrix *LLT;
  gint ret;
  guint i;

  _ncm_data_gauss_update_inv_cov (gauss, mset);
  LLT = _ncm_data_gauss_peek_LLT (gauss);

  ncm_rng_lock (rng);

  for (i = 0; i < self->np; i++)
  {
    const gdouble u_i = ncm_rng_ugaussian_gen (rng);

    ncm_vector_set (self->v, i, u_i);
  }

  ncm_rng_unlock (rng);

  /* CblasLower, CblasTrans => CblasUpper, CblasNoTrans */
  ret = gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit,
                        ncm_matrix_gsl (LLT), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_resample", ret);

  gauss_class->mean_func (gauss, mset, self->y);

  ncm_vector_sub (self->y, self->v);
}

static void
_ncm_data_gauss_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (data);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);
  NcmDataGaussClass *gauss_class   = NCM_DATA_GAUSS_GET_CLASS (gauss);
  guint i, j;

  *m2lnL = 0.0;

  gauss_class->mean_func (gauss, mset, self->v);
  ncm_vector_sub (self->v, self->y);

  _ncm_data_gauss_update_inv_cov (gauss, mset);

  if (!ncm_data_bootstrap_enabled (data))
  {
    for (i = 0; i < self->np; i++)
    {
      const gdouble f_i = ncm_vector_get (self->v, i);
      gdouble u_i       = 0.0;

      for (j = 0; j < self->np; j++)
        u_i += ncm_matrix_get (self->inv_cov, i, j) * ncm_vector_get (self->v, j);

      *m2lnL += u_i * f_i;
    }
  }
  else
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
    const guint bsize    = ncm_bootstrap_get_bsize (bstrap);
    NcmMatrix *LLT       = _ncm_data_gauss_peek_LLT (gauss);
    gint ret;

    /* CblasLower, CblasTrans => CblasUpper, CblasNoTrans */
    ret = gsl_blas_dtrmv (CblasUpper, CblasNoTrans, CblasNonUnit,
                          ncm_matrix_gsl (LLT), ncm_vector_gsl (self->v));
    NCM_TEST_GSL_RESULT ("_ncm_data_gauss_m2lnL_val", ret);

    for (i = 0; i < bsize; i++)
    {
      guint k     = ncm_bootstrap_get (bstrap, i);
      gdouble u_i = ncm_vector_get (self->v, k);

      *m2lnL += u_i * u_i;
    }
  }
}

static void
_ncm_data_gauss_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataGauss *gauss              = NCM_DATA_GAUSS (data);
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);
  NcmDataGaussClass *gauss_class   = NCM_DATA_GAUSS_GET_CLASS (gauss);
  gint ret;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("_ncm_data_gauss_leastsquares_f: data `%s': bootstrap is not supported with least squares.",
             ncm_data_peek_desc (data));

  gauss_class->mean_func (gauss, mset, v);
  ncm_vector_sub (v, self->y);

  _ncm_data_gauss_update_inv_cov (gauss, mset);

  /* CblasLower, CblasTrans => CblasUpper, CblasNoTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasNoTrans, CblasNonUnit,
                        ncm_matrix_gsl (_ncm_data_gauss_peek_LLT (gauss)), ncm_vector_gsl (v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_leastsquares_f", ret);
}

static void
_ncm_data_gauss_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  NcmDataGauss *gauss            = NCM_DATA_GAUSS (data);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_GET_CLASS (gauss);

  gauss_class->mean_func (gauss, mset, mu);
}

static void
_ncm_data_gauss_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (data);
  gint ret;

  _ncm_data_gauss_update_inv_cov (gauss, mset);

  ret = gsl_blas_dtrmm (CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
                        1.0, ncm_matrix_gsl (_ncm_data_gauss_peek_LLT (gauss)), ncm_matrix_gsl (H));

  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_inv_cov_UH", ret);
}

static void
_ncm_data_gauss_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (data);
  gint ret;

  _ncm_data_gauss_update_inv_cov (gauss, mset);

  ret = gsl_blas_dtrmv (CblasUpper, CblasNoTrans, CblasNonUnit,
                        ncm_matrix_gsl (_ncm_data_gauss_peek_LLT (gauss)), ncm_vector_gsl (f));

  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_inv_cov_Uf", ret);
}

void
_ncm_data_gauss_set_size (NcmDataGauss *gauss, guint np)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);
  NcmData *data                    = NCM_DATA (gauss);

  if ((np == 0) || (np != self->np))
  {
    self->np = 0;
    ncm_vector_clear (&self->y);
    ncm_vector_clear (&self->v);
    ncm_matrix_clear (&self->inv_cov);
    ncm_matrix_clear (&self->LLT);
    self->prepared_LLT = FALSE;

    ncm_data_set_init (data, FALSE);
  }

  if ((np != 0) && (np != self->np))
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

    self->np      = np;
    self->y       = ncm_vector_new (self->np);
    self->v       = ncm_vector_new (self->np);
    self->inv_cov = ncm_matrix_new (self->np, self->np);

    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (bstrap, np);
      ncm_bootstrap_set_bsize (bstrap, np);
    }

    ncm_data_set_init (data, FALSE);
  }
}

static guint
_ncm_data_gauss_get_size (NcmDataGauss *gauss)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  return self->np;
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

/**
 * ncm_data_gauss_peek_inv_cov:
 * @gauss: a #NcmDataGauss
 *
 * Gets the inverse covariance matrix.
 *
 * Returns: (transfer none): Inverse covariance matrix.
 */
NcmMatrix *
ncm_data_gauss_peek_inv_cov (NcmDataGauss *gauss)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  return self->inv_cov;
}

/**
 * ncm_data_gauss_peek_mean:
 * @gauss: a #NcmDataGauss
 *
 * Gets the mean vector.
 *
 * Returns: (transfer none): Mean vector.
 */
NcmVector *
ncm_data_gauss_peek_mean (NcmDataGauss *gauss)
{
  NcmDataGaussPrivate * const self = ncm_data_gauss_get_instance_private (gauss);

  return self->y;
}

