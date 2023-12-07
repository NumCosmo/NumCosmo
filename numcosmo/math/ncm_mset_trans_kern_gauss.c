/***************************************************************************
 *            ncm_mset_trans_kern_gauss.c
 *
 *  Wed September 03 14:55:40 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_gauss.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_mset_trans_kern_gauss
 * @title: NcmMSetTransKernGauss
 * @short_description: A multivariate gaussian sampler.
 *
 * This object subclasses NcmMSetTransKern and implements a multivariate gaussian
 * sampler.
 *
 * Implementation of a multivariate Gaussian sampler, providing a straightforward
 * method for generating random parameter vectors with multivariate parameters. This
 * sampler generates vectors with a normal distribution. The covariance of parameters
 * can be configured directly using ncm_mset_trans_kern_gauss_set_cov() or by
 * specifying individual standard deviations as parameter scales, assuming zero
 * correlation.
 *
 * **Key Functionality:**
 *
 * - Generates random parameter vectors with multivariate parameters.
 * - Utilizes a multivariate Gaussian distribution for sampling.
 * - Allows direct setting of covariance using ncm_mset_trans_kern_gauss_set_cov().
 * - Supports alternative methods:
 *    - Using ncm_mset_trans_kern_gauss_set_cov_from_scale() sets covariance using
 *      the scale property of parameters as standard deviation with zero correlation.
 *    - Using ncm_mset_trans_kern_gauss_set_cov_from_rescale() sets covariance using
 *      the scale property of parameters times @epsilon as standard deviation with zero
 *      correlation.
 *
 * This implementation is particularly useful when a Gaussian sampling approach is
 * required for generating random parameter vectors with multivariate parameters,
 * offering flexibility in specifying covariance through direct settings or individual
 * standard deviations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_trans_kern_gauss.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_LEN,
  PROP_COV,
  PROP_SIZE
};

G_DEFINE_TYPE (NcmMSetTransKernGauss, ncm_mset_trans_kern_gauss, NCM_TYPE_MSET_TRANS_KERN)

static void
ncm_mset_trans_kern_gauss_init (NcmMSetTransKernGauss *tkerng)
{
  tkerng->len  = 0;
  tkerng->cov  = NULL;
  tkerng->LLT  = NULL;
  tkerng->v    = NULL;
  tkerng->init = FALSE;
}

static void
ncm_mset_trans_kern_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernGauss *tkerng = NCM_MSET_TRANS_KERN_GAUSS (object);

  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_GAUSS (object));

  switch (prop_id)
  {
    case PROP_LEN:
      ncm_mset_trans_kern_gauss_set_size (tkerng, g_value_get_uint (value));
      break;
    case PROP_COV:
      ncm_mset_trans_kern_gauss_set_cov (tkerng, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mset_trans_kern_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernGauss *tkerng = NCM_MSET_TRANS_KERN_GAUSS (object);

  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_GAUSS (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, tkerng->len);
      break;
    case PROP_COV:
      g_value_set_object (value, tkerng->cov);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mset_trans_kern_gauss_dispose (GObject *object)
{
  NcmMSetTransKernGauss *tkerng = NCM_MSET_TRANS_KERN_GAUSS (object);

  ncm_matrix_clear (&tkerng->cov);
  ncm_matrix_clear (&tkerng->LLT);
  ncm_vector_clear (&tkerng->v);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_gauss_parent_class)->dispose (object);
}

static void
ncm_mset_trans_kern_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_gauss_parent_class)->finalize (object);
}

static void _ncm_mset_trans_kern_gauss_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset);
static void _ncm_mset_trans_kern_gauss_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static gdouble _ncm_mset_trans_kern_gauss_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
static const gchar *_ncm_mset_trans_kern_gauss_get_name (NcmMSetTransKern *tkern);

static void
ncm_mset_trans_kern_gauss_class_init (NcmMSetTransKernGaussClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmMSetTransKernClass *tkern_class = NCM_MSET_TRANS_KERN_CLASS (klass);

  object_class->set_property = ncm_mset_trans_kern_gauss_set_property;
  object_class->get_property = ncm_mset_trans_kern_gauss_get_property;
  object_class->dispose      = ncm_mset_trans_kern_gauss_dispose;
  object_class->finalize     = ncm_mset_trans_kern_gauss_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "length",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_COV,
                                   g_param_spec_object ("cov",
                                                        NULL,
                                                        "covariance",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  tkern_class->set_mset = &_ncm_mset_trans_kern_gauss_set_mset;
  tkern_class->generate = &_ncm_mset_trans_kern_gauss_generate;
  tkern_class->pdf      = &_ncm_mset_trans_kern_gauss_pdf;
  tkern_class->get_name = &_ncm_mset_trans_kern_gauss_get_name;
}

static void
_ncm_mset_trans_kern_gauss_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  NCM_UNUSED (mset);

  guint fparam_len = ncm_mset_fparam_len (tkern->mset);

  ncm_mset_trans_kern_gauss_set_size (NCM_MSET_TRANS_KERN_GAUSS (tkern), fparam_len);
}

static void
_ncm_mset_trans_kern_gauss_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernGauss *tkerng = NCM_MSET_TRANS_KERN_GAUSS (tkern);
  gint ret;
  guint i;

  g_assert (tkerng->init);

  while (TRUE)
  {
    ncm_rng_lock (rng);

    for (i = 0; i < tkerng->len; i++)
    {
      const gdouble u_i = gsl_ran_ugaussian (rng->r);

      ncm_vector_set (thetastar, i, u_i);
    }

    ncm_rng_unlock (rng);

    ret = gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit,
                          ncm_matrix_gsl (tkerng->LLT), ncm_vector_gsl (thetastar));
    NCM_TEST_GSL_RESULT ("ncm_mset_trans_kern_gauss_sample", ret);

    ncm_vector_add (thetastar, theta);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }
}

static gdouble
_ncm_mset_trans_kern_gauss_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar)
{
  NcmMSetTransKernGauss *tkerng = NCM_MSET_TRANS_KERN_GAUSS (tkern);
  gdouble m2lnP                 = 0.0;
  gint ret;
  guint i;

  g_assert (tkerng->init);

  ncm_vector_memcpy (tkerng->v, theta);
  ncm_vector_sub (tkerng->v, thetastar);

  ret = gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit,
                        ncm_matrix_gsl (tkerng->LLT), ncm_vector_gsl (tkerng->v));
  NCM_TEST_GSL_RESULT ("_ncm_mset_trans_kern_gauss_pdf", ret);

  ret = gsl_blas_ddot (ncm_vector_gsl (tkerng->v),
                       ncm_vector_gsl (tkerng->v),
                       &m2lnP);
  NCM_TEST_GSL_RESULT ("_ncm_mset_trans_kern_gauss_pdf", ret);

  m2lnP += tkerng->len * ncm_c_ln2pi ();

  for (i = 0; i < tkerng->len; i++)
  {
    m2lnP += 2.0 * log (ncm_matrix_get (tkerng->LLT, i, i));
  }

  return exp (-0.5 * m2lnP);
}

static const gchar *
_ncm_mset_trans_kern_gauss_get_name (NcmMSetTransKern *tkern)
{
  return "Multivariate Gaussian Sampler";
}

/**
 * ncm_mset_trans_kern_gauss_new:
 * @len: Number of variables
 *
 * New NcmMSetTransKern gauss for @len multivariate gaussian.
 *
 * Returns: (transfer full): a new #NcmMSetTransKernGauss.
 *
 */
NcmMSetTransKernGauss *
ncm_mset_trans_kern_gauss_new (guint len)
{
  NcmMSetTransKernGauss *tkerng = g_object_new (NCM_TYPE_MSET_TRANS_KERN_GAUSS,
                                                "length", len,
                                                NULL);

  return tkerng;
}

/**
 * ncm_mset_trans_kern_gauss_set_size:
 * @tkerng: a #NcmMSetTransKernGauss.
 * @len: Number of variables.
 *
 * Sets size of #NcmMSetTransKernGauss.
 *
 */
void
ncm_mset_trans_kern_gauss_set_size (NcmMSetTransKernGauss *tkerng, guint len)
{
  if ((len == 0) || (len != tkerng->len))
  {
    tkerng->len = 0;
    ncm_matrix_clear (&tkerng->cov);
    ncm_matrix_clear (&tkerng->LLT);
    ncm_vector_clear (&tkerng->v);
  }

  if ((len != 0) && (len != tkerng->len))
  {
    tkerng->len = len;
    tkerng->cov = ncm_matrix_new (tkerng->len, tkerng->len);
    tkerng->LLT = ncm_matrix_new (tkerng->len, tkerng->len);
    tkerng->v   = ncm_vector_new (tkerng->len);
  }
}

/**
 * ncm_mset_trans_kern_gauss_get_size:
 * @tkerng: a #NcmMSetTransKernGauss.
 *
 * Gets size of #NcmMSetTransKernGauss.
 *
 * Returns: size of the gaussian multivariate.
 */
guint
ncm_mset_trans_kern_gauss_get_size (NcmMSetTransKernGauss *tkerng)
{
  return tkerng->len;
}

/**
 * ncm_mset_trans_kern_gauss_set_cov:
 * @tkerng: a #NcmMSetTransKernGauss.
 * @cov: a #NcmMatrix.
 *
 * Sets the covariance given by the #NcmMatrix @cov.
 *
 */
void
ncm_mset_trans_kern_gauss_set_cov (NcmMSetTransKernGauss *tkerng, const NcmMatrix *cov)
{
  gint ret;

  g_assert_cmpuint (ncm_matrix_ncols (tkerng->cov), ==, ncm_matrix_ncols (cov));
  g_assert_cmpuint (ncm_matrix_nrows (tkerng->cov), ==, ncm_matrix_nrows (cov));
  ncm_matrix_memcpy (tkerng->cov, cov);
  ncm_matrix_memcpy (tkerng->LLT, cov);

  ret = ncm_matrix_cholesky_decomp (tkerng->LLT, 'L');

  if (ret != 0)
    g_error ("ncm_mset_trans_kern_gauss_set_cov[ncm_matrix_cholesky_decomp]: %d.", ret);

  tkerng->init = TRUE;
}

/**
 * ncm_mset_trans_kern_gauss_set_cov_variant:
 * @tkerng: a #NcmMSetTransKernGauss.
 * @cov: a #GVariant.
 *
 * Sets the covariance given by the #GVariant @cov.
 *
 */
void
ncm_mset_trans_kern_gauss_set_cov_variant (NcmMSetTransKernGauss *tkerng, GVariant *cov)
{
  gint ret;

  ncm_matrix_set_from_variant (tkerng->cov, cov);
  ncm_matrix_memcpy (tkerng->LLT, tkerng->cov);

  ret = ncm_matrix_cholesky_decomp (tkerng->LLT, 'L');

  if (ret != 0)
    g_error ("ncm_mset_trans_kern_gauss_set_cov_variant[ncm_matrix_cholesky_decomp]: %d.", ret);

  tkerng->init = TRUE;
}

/**
 * ncm_mset_trans_kern_gauss_set_cov_data:
 * @tkerng: a #NcmMSetTransKernGauss.
 * @cov: a #GVariant.
 *
 * Sets the covariance given by the double array @cov.
 *
 */
void
ncm_mset_trans_kern_gauss_set_cov_data (NcmMSetTransKernGauss *tkerng, gdouble *cov)
{
  gint ret;

  ncm_matrix_set_from_data (tkerng->cov, cov);
  ncm_matrix_memcpy (tkerng->LLT, tkerng->cov);

  ret = ncm_matrix_cholesky_decomp (tkerng->LLT, 'L');

  if (ret != 0)
    g_error ("ncm_mset_trans_kern_gauss_set_cov_variant[ncm_matrix_cholesky_decomp]: %d.", ret);

  tkerng->init = TRUE;
}

/**
 * ncm_mset_trans_kern_gauss_get_cov:
 * @tkerng: a #NcmMSetTransKernGauss.
 *
 * Gets the covariance.
 *
 * Returns: (transfer full): the covariance.
 */
NcmMatrix *
ncm_mset_trans_kern_gauss_get_cov (NcmMSetTransKernGauss *tkerng)
{
  return ncm_matrix_ref (tkerng->cov);
}

/**
 * ncm_mset_trans_kern_gauss_set_cov_from_scale:
 * @tkerng: a #NcmMSetTransKernGauss
 *
 * Sets the covariance using the scale property of the parameters as
 * standard deviation and zero correlation.
 *
 */
void
ncm_mset_trans_kern_gauss_set_cov_from_scale (NcmMSetTransKernGauss *tkerng)
{
  NcmMSetTransKern *tkern = NCM_MSET_TRANS_KERN (tkerng);
  guint i;
  gint ret;

  g_assert (tkern->mset != NULL);

  ncm_matrix_set_identity (tkerng->cov);

  for (i = 0; i < tkerng->len; i++)
  {
    const gdouble scale = ncm_mset_fparam_get_scale (tkern->mset, i);

    ncm_matrix_set (tkerng->cov, i, i, scale * scale);
  }

  ncm_matrix_memcpy (tkerng->LLT, tkerng->cov);

  ret = ncm_matrix_cholesky_decomp (tkerng->LLT, 'L');

  if (ret != 0)
    g_error ("ncm_mset_trans_kern_gauss_set_cov_from_scale[ncm_matrix_cholesky_decomp]: %d.", ret);

  tkerng->init = TRUE;
}

/**
 * ncm_mset_trans_kern_gauss_set_cov_from_rescale:
 * @tkerng: a #NcmMSetTransKernGauss
 * @epsilon: the overall rescale
 *
 * Sets the covariance using the scale property of the parameters times
 * @epsilon as standard deviation and zero correlation.
 *
 */
void
ncm_mset_trans_kern_gauss_set_cov_from_rescale (NcmMSetTransKernGauss *tkerng, const gdouble epsilon)
{
  NcmMSetTransKern *tkern = NCM_MSET_TRANS_KERN (tkerng);
  guint i;
  gint ret;

  g_assert (tkern->mset != NULL);

  ncm_matrix_set_identity (tkerng->cov);

  for (i = 0; i < tkerng->len; i++)
  {
    const gdouble scale = ncm_mset_fparam_get_scale (tkern->mset, i) * epsilon;

    ncm_matrix_set (tkerng->cov, i, i, scale * scale);
  }

  ncm_matrix_memcpy (tkerng->LLT, tkerng->cov);

  ret = ncm_matrix_cholesky_decomp (tkerng->LLT, 'L');

  if (ret != 0)
    g_error ("ncm_mset_trans_kern_gauss_set_cov_from_scale[ncm_matrix_cholesky_decomp]: %d.", ret);

  tkerng->init = TRUE;
}

