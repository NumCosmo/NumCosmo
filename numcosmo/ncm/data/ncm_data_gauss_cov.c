/***************************************************************************
 *            ncm_data_gauss_cov.c
 *
 *  Tue November 20 18:46:02 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * NcmDataGaussCov:
 *
 * Gaussian data given by its covariance.
 *
 * Abstract class for data $y$ of #NcmDataGaussCov:n-points points drawn from a
 * normal distribution with mean $\mu$ and covariance $C$ (#NcmDataGaussCov:cov);
 * $-2\ln L = (\mu - y)^T C^{-1} (\mu - y)$, plus $n\ln 2\pi + \ln\det C$ when
 * #NcmDataGaussCov:use-norma is set. Every evaluation uses the Cholesky factor
 * $C = U^T U$, computed when needed.
 *
 * Subclasses implement the mean_func virtual method, which fills the mean $\mu$
 * given the models in a #NcmMSet. When $C$ depends on the models, they also
 * implement cov_func, which fills the matrix passed to it and returns %TRUE when it
 * changed that matrix; the Cholesky factor is then recomputed before its next use.
 * ncm_data_gauss_cov_compute_cov() calls cov_func on a copy, so cov_func must decide
 * on the matrix it is given, not on what it returned before.
 *
 * Replacing the covariance with ncm_data_gauss_cov_set_cov() or
 * #NcmDataGaussCov:cov invalidates the factor, while changing the matrix returned by
 * ncm_data_gauss_cov_peek_cov() in place after the first evaluation does not:
 * replace it instead.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/data/ncm_data_gauss_cov.h"
#include "ncm/algebra/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_USE_NORMA,
  PROP_MEAN,
  PROP_COV,
  PROP_SIZE,
};

typedef struct _NcmDataGaussCovPrivate
{
  guint np;
  NcmVector *y;
  NcmVector *v;
  NcmMatrix *cov;
  NcmMatrix *LLT;
  gboolean prepared_LLT;
  gboolean use_norma;
} NcmDataGaussCovPrivate;


G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmDataGaussCov, ncm_data_gauss_cov, NCM_TYPE_DATA) /* LCOV_EXCL_BR_LINE */

static void
ncm_data_gauss_cov_init (NcmDataGaussCov *gauss)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  self->np           = 0;
  self->y            = NULL;
  self->v            = NULL;
  self->cov          = NULL;
  self->LLT          = NULL;
  self->prepared_LLT = FALSE;
  self->use_norma    = FALSE;
}

static void
_ncm_data_gauss_cov_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataGaussCov *gauss              = NCM_DATA_GAUSS_COV (object);
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  g_return_if_fail (NCM_IS_DATA_GAUSS_COV (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_gauss_cov_set_size (gauss, g_value_get_uint (value));
      break;
    case PROP_USE_NORMA:
      ncm_data_gauss_cov_use_norma (gauss, g_value_get_boolean (value));
      break;
    case PROP_MEAN:
    {
      NcmVector *v = g_value_get_object (value);

      ncm_vector_substitute (&self->y, v, TRUE);
      break;
    }
    case PROP_COV:
      ncm_matrix_substitute (&self->cov, g_value_get_object (value), TRUE);
      self->prepared_LLT = FALSE;
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_data_gauss_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataGaussCov *gauss              = NCM_DATA_GAUSS_COV (object);
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  g_return_if_fail (NCM_IS_DATA_GAUSS_COV (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, self->np);
      break;
    case PROP_USE_NORMA:
      g_value_set_boolean (value, self->use_norma);
      break;
    case PROP_MEAN:
      g_value_set_object (value, self->y);
      break;
    case PROP_COV:
      g_value_set_object (value, self->cov);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_data_gauss_cov_dispose (GObject *object)
{
  NcmDataGaussCov *gauss              = NCM_DATA_GAUSS_COV (object);
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  ncm_vector_clear (&self->y);
  ncm_vector_clear (&self->v);
  ncm_matrix_clear (&self->cov);
  ncm_matrix_clear (&self->LLT);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_cov_parent_class)->dispose (object);
}

static guint _ncm_data_gauss_cov_get_length (NcmData *data);

/* static void _ncm_data_gauss_cov_begin (NcmData *data); */
static void _ncm_data_gauss_cov_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_gauss_cov_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_cov_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_gauss_cov_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu);
static void _ncm_data_gauss_cov_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H);
static void _ncm_data_gauss_cov_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f);

static void _ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np);
static guint _ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss);
static void _ncm_data_gauss_cov_lnNorma2 (NcmDataGaussCov *gauss, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_cov_lnNorma2_bs (NcmDataGaussCov *gauss, NcmMSet *mset, NcmBootstrap *bstrap, gdouble *m2lnL);

static void
ncm_data_gauss_cov_class_init (NcmDataGaussCovClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class              = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &_ncm_data_gauss_cov_set_property;
  object_class->get_property = &_ncm_data_gauss_cov_get_property;
  object_class->dispose      = &_ncm_data_gauss_cov_dispose;

  /**
   * NcmDataGaussCov:n-points:
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
   * NcmDataGaussCov:use-norma:
   *
   * Whether $-2\ln L$ includes the normalization $n\ln 2\pi + \ln\det C$.
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_NORMA,
                                   g_param_spec_boolean ("use-norma",
                                                         NULL,
                                                         "Use the likelihood normalization to calculate -2lnL",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmDataGaussCov:mean:
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
   * NcmDataGaussCov:cov:
   *
   * The covariance $C$; setting it invalidates its Cholesky factor.
   */
  g_object_class_install_property (object_class,
                                   PROP_COV,
                                   g_param_spec_object ("cov",
                                                        NULL,
                                                        "Data covariance",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap = TRUE;

  data_class->get_length = &_ncm_data_gauss_cov_get_length;
  data_class->get_dof    = NULL;
  data_class->begin      = NULL;

  data_class->resample       = &_ncm_data_gauss_cov_resample;
  data_class->m2lnL_val      = &_ncm_data_gauss_cov_m2lnL_val;
  data_class->leastsquares_f = &_ncm_data_gauss_cov_leastsquares_f;

  data_class->mean_vector = &_ncm_data_gauss_cov_mean_vector;
  data_class->inv_cov_UH  = &_ncm_data_gauss_cov_inv_cov_UH;
  data_class->inv_cov_Uf  = &_ncm_data_gauss_cov_inv_cov_Uf;

  gauss_cov_class->mean_func   = NULL;
  gauss_cov_class->cov_func    = NULL;
  gauss_cov_class->lnNorma2    = &_ncm_data_gauss_cov_lnNorma2;
  gauss_cov_class->lnNorma2_bs = &_ncm_data_gauss_cov_lnNorma2_bs;
  gauss_cov_class->set_size    = &_ncm_data_gauss_cov_set_size;
  gauss_cov_class->get_size    = &_ncm_data_gauss_cov_get_size;
}

static guint
_ncm_data_gauss_cov_get_length (NcmData *data)
{
  NcmDataGaussCov *gauss              = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  return self->np;
}

/*
 * Calls cov_func, if any; when it changed the covariance, the Cholesky factor is
 * marked stale. Every path that may update the matrix goes through here.
 */
static void
_ncm_data_gauss_cov_update_cov (NcmDataGaussCov *gauss, NcmMSet *mset)
{
  NcmDataGaussCovPrivate * const self   = ncm_data_gauss_cov_get_instance_private (gauss);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);

  if ((gauss_cov_class->cov_func != NULL) && gauss_cov_class->cov_func (gauss, mset, self->cov))
    self->prepared_LLT = FALSE;
}

/* The upper Cholesky factor U of C = U^T U, recomputed when stale. */
static NcmMatrix *
_ncm_data_gauss_cov_peek_LLT (NcmDataGaussCov *gauss)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  if (!self->prepared_LLT)
  {
    gint ret;

    if (self->LLT == NULL)
      self->LLT = ncm_matrix_dup (self->cov);
    else
      ncm_matrix_memcpy (self->LLT, self->cov);

    ret = ncm_matrix_cholesky_decomp (self->LLT, 'U');

    if (ret != 0)
      g_error ("_ncm_data_gauss_cov_peek_LLT: data `%s': the covariance is not positive "
               "definite (Cholesky decomposition returned %d).",
               ncm_data_peek_desc (NCM_DATA (gauss)), ret);

    self->prepared_LLT = TRUE;
  }

  return self->LLT;
}

static void
_ncm_data_gauss_cov_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataGaussCov *gauss                = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovPrivate * const self   = ncm_data_gauss_cov_get_instance_private (gauss);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gint ret;
  guint i;

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  ncm_rng_lock (rng);

  for (i = 0; i < self->np; i++)
  {
    const gdouble u_i = ncm_rng_ugaussian_gen (rng);

    ncm_vector_set (self->v, i, u_i);
  }

  ncm_rng_unlock (rng);

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (_ncm_data_gauss_cov_peek_LLT (gauss)), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_resample", ret);

  gauss_cov_class->mean_func (gauss, mset, self->y);
  ncm_vector_sub (self->y, self->v);
}

static void
_ncm_data_gauss_cov_lnNorma2 (NcmDataGaussCov *gauss, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  *m2lnL += self->np * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (_ncm_data_gauss_cov_peek_LLT (gauss));
}

static void
_ncm_data_gauss_cov_lnNorma2_bs (NcmDataGaussCov *gauss, NcmMSet *mset, NcmBootstrap *bstrap, gdouble *m2lnL)
{
  NcmMatrix *LLT    = _ncm_data_gauss_cov_peek_LLT (gauss);
  const gdouble lb  = 1.0e-200;
  const gdouble ub  = 1.0e+200;
  gdouble detL      = 1.0;
  glong exponent    = 0;
  const guint bsize = ncm_bootstrap_get_bsize (bstrap);
  guint i;

  for (i = 0; i < bsize; i++)
  {
    const guint k       = ncm_bootstrap_get (bstrap, i);
    const gdouble Lii   = ncm_matrix_get (LLT, k, k);
    const gdouble ndetL = detL * Lii;

    if (G_UNLIKELY ((ndetL < lb) || (ndetL > ub)))
    {
      gint exponent_i = 0;

      detL      = frexp (ndetL, &exponent_i);
      exponent += exponent_i;
    }
    else
    {
      detL = ndetL;
    }
  }

  /* One Gaussian term per bootstrap draw, as in the sum of squares. */
  *m2lnL += bsize * ncm_c_ln2pi () + 2.0 * (log (detL) + exponent * M_LN2);
}

static void
_ncm_data_gauss_cov_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGaussCov *gauss                = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovPrivate * const self   = ncm_data_gauss_cov_get_instance_private (gauss);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gint ret;

  *m2lnL = 0.0;

  gauss_cov_class->mean_func (gauss, mset, self->v);

  ncm_vector_sub (self->v, self->y);

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (_ncm_data_gauss_cov_peek_LLT (gauss)), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_m2lnL_val", ret);

  if (!ncm_data_bootstrap_enabled (data))
  {
    ret = gsl_blas_ddot (ncm_vector_gsl (self->v),
                         ncm_vector_gsl (self->v),
                         m2lnL);
    NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_m2lnL_val", ret);

    if (self->use_norma)
      gauss_cov_class->lnNorma2 (gauss, mset, m2lnL);
  }
  else
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
    const guint bsize    = ncm_bootstrap_get_bsize (bstrap);
    guint i;

    if (bsize == 0)
    {
      *m2lnL = 0.0;

      return;
    }

    g_assert (ncm_bootstrap_is_init (bstrap));

    for (i = 0; i < bsize; i++)
    {
      guint k           = ncm_bootstrap_get (bstrap, i);
      const gdouble u_i = ncm_vector_get (self->v, k);

      *m2lnL += u_i * u_i;
    }

    if (self->use_norma)
      gauss_cov_class->lnNorma2_bs (gauss, mset, bstrap, m2lnL);
  }
}

static void
_ncm_data_gauss_cov_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataGaussCov *gauss                = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovPrivate * const self   = ncm_data_gauss_cov_get_instance_private (gauss);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gint ret;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("_ncm_data_gauss_cov_leastsquares_f: data `%s': bootstrap is not supported with least squares.",
             ncm_data_peek_desc (data));

  gauss_cov_class->mean_func (gauss, mset, v);
  ncm_vector_sub (v, self->y);

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (_ncm_data_gauss_cov_peek_LLT (gauss)), ncm_vector_gsl (v));
  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_leastsquares_f", ret);
}

static void
_ncm_data_gauss_cov_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  NcmDataGaussCov *gauss                = NCM_DATA_GAUSS_COV (data);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);

  gauss_cov_class->mean_func (gauss, mset, mu);
}

static void
_ncm_data_gauss_cov_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (data);
  gint ret;

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  ret = gsl_blas_dtrsm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                        1.0, ncm_matrix_gsl (_ncm_data_gauss_cov_peek_LLT (gauss)), ncm_matrix_gsl (H));

  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_inv_cov_UH", ret);
}

static void
_ncm_data_gauss_cov_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (data);
  gint ret;

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (_ncm_data_gauss_cov_peek_LLT (gauss)), ncm_vector_gsl (f));

  NCM_TEST_GSL_RESULT ("_ncm_data_gauss_cov_inv_cov_Uf", ret);
}

static void
_ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);
  NcmData *data                       = NCM_DATA (gauss);

  if ((np == 0) || (np != self->np))
  {
    self->np = 0;
    ncm_vector_clear (&self->y);
    ncm_vector_clear (&self->v);
    ncm_matrix_clear (&self->cov);
    ncm_matrix_clear (&self->LLT);
    self->prepared_LLT = FALSE;

    ncm_data_set_init (data, FALSE);
  }

  if ((np != 0) && (np != self->np))
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

    self->np  = np;
    self->y   = ncm_vector_new (self->np);
    self->v   = ncm_vector_new (self->np);
    self->cov = ncm_matrix_new (self->np, self->np);

    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (bstrap, np);
      ncm_bootstrap_set_bsize (bstrap, np);
    }

    ncm_data_set_init (data, FALSE);
  }
}

static guint
_ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  return self->np;
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

/**
 * ncm_data_gauss_cov_use_norma:
 * @gauss: a #NcmDataGaussCov
 * @use_norma: a boolean
 *
 * Sets whether the value of $-2\ln(L)$ will be properly normalized.
 *
 */
void
ncm_data_gauss_cov_use_norma (NcmDataGaussCov *gauss, gboolean use_norma)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  self->use_norma = use_norma;
}

/**
 * ncm_data_gauss_cov_replace_mean:
 * @gauss: a #NcmDataGaussCov
 * @mean: new mean #NcmVector
 *
 * Replaces the current mean vector for @mean.
 */
void
ncm_data_gauss_cov_replace_mean (NcmDataGaussCov *gauss, NcmVector *mean)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  ncm_vector_substitute (&self->y, mean, TRUE);
}

/**
 * ncm_data_gauss_cov_peek_mean:
 * @gauss: a #NcmDataGaussCov
 *
 * Returns: (transfer none): the current data mean #NcmVector.
 */
NcmVector *
ncm_data_gauss_cov_peek_mean (NcmDataGaussCov *gauss)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  return self->y;
}

/**
 * ncm_data_gauss_cov_set_cov:
 * @gauss: a #NcmDataGaussCov
 * @cov: a #NcmMatrix
 *
 * Replaces the covariance matrix with @cov, which must be n-points x n-points, and
 * invalidates its Cholesky factor.
 *
 */
void
ncm_data_gauss_cov_set_cov (NcmDataGaussCov *gauss, NcmMatrix *cov)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  if ((ncm_matrix_nrows (cov) != self->np) || (ncm_matrix_ncols (cov) != self->np))
    g_error ("ncm_data_gauss_cov_set_cov: data `%s' has %u points, but the covariance is %u x %u.",
             ncm_data_peek_desc (NCM_DATA (gauss)), self->np, ncm_matrix_nrows (cov), ncm_matrix_ncols (cov));

  ncm_matrix_ref (cov);
  ncm_matrix_clear (&self->cov);
  self->cov          = cov;
  self->prepared_LLT = FALSE;
}

/**
 * ncm_data_gauss_cov_peek_cov:
 * @gauss: a #NcmDataGaussCov
 *
 * Returns: (transfer none): the current data covariance #NcmMatrix.
 */
NcmMatrix *
ncm_data_gauss_cov_peek_cov (NcmDataGaussCov *gauss)
{
  NcmDataGaussCovPrivate * const self = ncm_data_gauss_cov_get_instance_private (gauss);

  return self->cov;
}

/**
 * ncm_data_gauss_cov_compute_cov:
 * @gauss: a #NcmDataGaussCov
 * @mset: a #NcmMSet
 * @updated: (out) (allow-none): a #gboolean
 *
 * Computes the covariance matrix based on the models in @mset. If
 * the ncm_data_gauss_cov_func is not set, returns %NULL.
 *
 * Returns: (transfer full): the current data covariance #NcmMatrix.
 */
NcmMatrix *
ncm_data_gauss_cov_compute_cov (NcmDataGaussCov *gauss, NcmMSet *mset, gboolean *updated)
{
  NcmDataGaussCovClass * const gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  NcmDataGaussCovPrivate * const self          = ncm_data_gauss_cov_get_instance_private (gauss);

  if (gauss_cov_class->cov_func != NULL)
  {
    NcmMatrix *cov = ncm_matrix_dup (self->cov);
    gboolean l_updated;

    ncm_data_prepare (NCM_DATA (gauss), mset);
    l_updated = gauss_cov_class->cov_func (gauss, mset, cov);

    if (updated != NULL)
      *updated = l_updated;

    return cov;
  }

  return NULL;
}

/**
 * ncm_data_gauss_cov_get_log_norma:
 * @gauss: a #NcmDataGaussCov
 * @mset: a #NcmMSet
 *
 * Returns: the log-normalization factor for $-2\ln(L)$.
 */
gdouble
ncm_data_gauss_cov_get_log_norma (NcmDataGaussCov *gauss, NcmMSet *mset)
{
  NcmDataGaussCovClass * const gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  gdouble m2lnL                                = 0.0;

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  gauss_cov_class->lnNorma2 (gauss, mset, &m2lnL);

  return m2lnL;
}

/**
 * ncm_data_gauss_cov_bulk_resample:
 * @gauss: a #NcmDataGaussCov
 * @mset: a #NcmMSet
 * @resample: a #NcmMatrix
 * @rng: a #NcmRNG
 *
 * Resamples the data based on the models in @mset according to the current
 * data distribution. The resampled data is stored in @resample. The resampling
 * is done in bulk, i.e., all the data is resampled at once.
 *
 */
void
ncm_data_gauss_cov_bulk_resample (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *resample, NcmRNG *rng)
{
  NcmData *data                         = NCM_DATA (gauss);
  NcmDataGaussCovPrivate * const self   = ncm_data_gauss_cov_get_instance_private (gauss);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  guint nrealizations                   = ncm_matrix_nrows (resample);
  gint ret;
  guint i, j;

  if ((ncm_matrix_ncols (resample) != self->np) || (nrealizations < 1))
    g_error ("ncm_data_gauss_cov_bulk_resample: data `%s' has %u points, but the matrix is %u x %u.",
             ncm_data_peek_desc (data), self->np, nrealizations, ncm_matrix_ncols (resample));

  ncm_data_prepare (data, mset);

  _ncm_data_gauss_cov_update_cov (gauss, mset);

  ncm_rng_lock (rng);

  for (j = 0; j < nrealizations; j++)
  {
    for (i = 0; i < self->np; i++)
    {
      const gdouble u_i = ncm_rng_ugaussian_gen (rng);

      ncm_matrix_set (resample, j, i, u_i);
    }
  }

  ncm_rng_unlock (rng);

  ret = gsl_blas_dtrmm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                        1.0, ncm_matrix_gsl (_ncm_data_gauss_cov_peek_LLT (gauss)),
                        ncm_matrix_gsl (resample));
  NCM_TEST_GSL_RESULT ("ncm_data_gauss_cov_bulk_resample", ret); /* LCOV_EXCL_BR_LINE */

  gauss_cov_class->mean_func (gauss, mset, self->v);

  for (j = 0; j < nrealizations; j++)
  {
    gdouble *theta_j = ncm_matrix_ptr (resample, j, 0);

    cblas_daxpy (self->np,
                 -1.0,
                 ncm_vector_const_data (self->v), ncm_vector_stride (self->v),
                 theta_j, 1);
  }

  ncm_matrix_scale (resample, -1.0);
}

/**
 * ncm_data_gauss_cov_compute_mean:
 * @gauss: a #NcmDataGaussCov
 * @mset: a #NcmMSet
 *
 * Computes the mean vector based on the models in @mset. If
 * the ncm_data_gauss_cov_mean_func is not set, returns %NULL.
 *
 * Returns: (transfer full): the current data mean #NcmVector.
 */
NcmVector *
ncm_data_gauss_cov_compute_mean (NcmDataGaussCov *gauss, NcmMSet *mset)
{
  NcmDataGaussCovClass * const gauss_cov_class = NCM_DATA_GAUSS_COV_GET_CLASS (gauss);
  NcmDataGaussCovPrivate * const self          = ncm_data_gauss_cov_get_instance_private (gauss);

  if (gauss_cov_class->mean_func != NULL)
  {
    NcmVector *mean = ncm_vector_dup (self->y);

    ncm_data_prepare (NCM_DATA (gauss), mset);
    gauss_cov_class->mean_func (gauss, mset, mean);

    return mean;
  }

  return NULL;
}

