/***************************************************************************
 *            ncm_mc_sampler_gauss.c
 *
 *  Wed September 03 14:55:40 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mc_sampler_gauss.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_mc_sampler_gauss
 * @title: Markov Chain Multivariate Gaussian Sampler 
 * @short_description: Object implementing a multivariate gaussian sampler.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm_mc_sampler_gauss.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

enum
{
  PROP_0,
  PROP_LEN,
  PROP_COV,
  PROP_SIZE
};

G_DEFINE_TYPE (NcmMCSamplerGauss, ncm_mc_sampler_gauss, NCM_TYPE_MC_SAMPLER);

static void
ncm_mc_sampler_gauss_init (NcmMCSamplerGauss *mcsg)
{
  mcsg->len  = 0;
  mcsg->cov  = NULL;
  mcsg->LLT  = NULL;
  mcsg->init = FALSE;
}

static void
ncm_mc_sampler_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMCSamplerGauss *mcsg = NCM_MC_SAMPLER_GAUSS (object);
  g_return_if_fail (NCM_IS_MC_SAMPLER_GAUSS (object));

  switch (prop_id)
  {
    case PROP_LEN:
      ncm_mc_sampler_gauss_set_size (mcsg, g_value_get_uint (value));
      break;
    case PROP_COV:
      ncm_mc_sampler_gauss_set_cov_variant (mcsg, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mc_sampler_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMCSamplerGauss *mcsg = NCM_MC_SAMPLER_GAUSS (object);
  g_return_if_fail (NCM_IS_MC_SAMPLER_GAUSS (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, mcsg->len);
      break;
    case PROP_COV:
      g_value_take_variant (value, ncm_matrix_get_variant (mcsg->cov));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mc_sampler_gauss_dispose (GObject *object)
{
  NcmMCSamplerGauss *mcsg = NCM_MC_SAMPLER_GAUSS (object);
  
  ncm_matrix_clear (&mcsg->cov);
  ncm_matrix_clear (&mcsg->LLT);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mc_sampler_gauss_parent_class)->dispose (object);
}

static void
ncm_mc_sampler_gauss_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mc_sampler_gauss_parent_class)->finalize (object);
}

static void _ncm_mc_sampler_gauss_set_mset (NcmMCSampler *mcs, NcmMSet *mset);
static void _ncm_mc_sampler_gauss_generate (NcmMCSampler *mcs, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static const gchar *_ncm_mc_sampler_gauss_get_name (NcmMCSampler *mcs);

static void
ncm_mc_sampler_gauss_class_init (NcmMCSamplerGaussClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmMCSamplerClass *mcs_class = NCM_MC_SAMPLER_CLASS (klass);

  object_class->set_property = ncm_mc_sampler_gauss_set_property;
  object_class->get_property = ncm_mc_sampler_gauss_get_property;
  object_class->dispose      = ncm_mc_sampler_gauss_dispose;
  object_class->finalize     = ncm_mc_sampler_gauss_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "length",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
    
  g_object_class_install_property (object_class,
                                   PROP_COV,
                                   g_param_spec_variant ("cov",
                                                         NULL,
                                                         "covariance",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  mcs_class->set_mset = &_ncm_mc_sampler_gauss_set_mset;
  mcs_class->generate = &_ncm_mc_sampler_gauss_generate;
  mcs_class->get_name = &_ncm_mc_sampler_gauss_get_name;
}

static void 
_ncm_mc_sampler_gauss_set_mset (NcmMCSampler *mcs, NcmMSet *mset)
{
  NCM_UNUSED (mset);
  guint fparam_len = ncm_mset_fparam_len (mcs->mset);
  ncm_mc_sampler_gauss_set_size (NCM_MC_SAMPLER_GAUSS (mcs), fparam_len);
}

static void 
_ncm_mc_sampler_gauss_generate (NcmMCSampler *mcs, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMCSamplerGauss *mcsg = NCM_MC_SAMPLER_GAUSS (mcs);
  gint ret;
  guint i, j;
  gboolean valid = FALSE;

  g_assert (mcsg->init);

  while (!valid)
  {
    ncm_rng_lock (rng);
    for (i = 0; i < mcsg->len; i++)
    {
      const gdouble u_i = gsl_ran_ugaussian (rng->r);
      ncm_vector_set (thetastar, i, u_i);
    }
    ncm_rng_unlock (rng);

    ret = gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, 
                          ncm_matrix_gsl (mcsg->LLT), ncm_vector_gsl (thetastar));
    NCM_TEST_GSL_RESULT ("ncm_mc_sampler_gauss_sample", ret);

    ncm_vector_add (thetastar, theta);

    valid = TRUE;
    for (j = 0; j < mcsg->len; j++)
    {
      const gdouble lb  = ncm_mset_fparam_get_lower_bound (mcs->mset, j);
      const gdouble ub  = ncm_mset_fparam_get_upper_bound (mcs->mset, j);
      const gdouble val = ncm_vector_get (thetastar, j);
      if (val < lb || val > ub)
      {
        valid = FALSE;
        break;
      }
    }
  }
}

static const gchar *
_ncm_mc_sampler_gauss_get_name (NcmMCSampler *mcs)
{
  return "Multivariate Gaussian Sampler";
}

/**
 * ncm_mc_sampler_gauss_new:
 * @len: Number of variables
 *
 * New NcmMCSampler gauss for @len multivariate gaussian.
 * 
 * Returns: (transfer full): a new #NcmMCSamplerGauss.
 * 
 */
NcmMCSamplerGauss *
ncm_mc_sampler_gauss_new (guint len)
{
  NcmMCSamplerGauss *mcsg = g_object_new (NCM_TYPE_MC_SAMPLER_GAUSS, 
                                          "length", len,
                                          NULL);
  return mcsg;
}

/**
 * ncm_mc_sampler_gauss_set_size:
 * @mcsg: a #NcmMCSamplerGauss.
 * @len: Number of variables.
 *
 * Sets size of #NcmMCSamplerGauss.
 * 
 */
void 
ncm_mc_sampler_gauss_set_size (NcmMCSamplerGauss *mcsg, guint len)
{
  if ((len == 0) || (len != mcsg->len))
  {
    mcsg->len = 0;
    ncm_matrix_clear (&mcsg->cov);
    ncm_matrix_clear (&mcsg->LLT);
  }
  if ((len != 0) && (len != mcsg->len))
  {
    mcsg->len = len;
    mcsg->cov = ncm_matrix_new (mcsg->len, mcsg->len);
    mcsg->LLT = ncm_matrix_new (mcsg->len, mcsg->len);
  }
}

/**
 * ncm_mc_sampler_gauss_get_size:
 * @mcsg: a #NcmMCSamplerGauss.
 *
 * Gets size of #NcmMCSamplerGauss.
 * 
 * Returns: size of the gaussian multivariate.
 */
guint 
ncm_mc_sampler_gauss_get_size (NcmMCSamplerGauss *mcsg)
{
  return mcsg->len;
}

/**
 * ncm_mc_sampler_gauss_set_cov:
 * @mcsg: a #NcmMCSamplerGauss.
 * @cov: a #NcmMatrix.
 *
 * Sets the covariance given by the #NcmMatrix @cov.
 * 
 */
void 
ncm_mc_sampler_gauss_set_cov (NcmMCSamplerGauss *mcsg, const NcmMatrix *cov)
{
  g_assert_cmpuint (ncm_matrix_ncols (mcsg->cov), ==, ncm_matrix_ncols (cov));
  g_assert_cmpuint (ncm_matrix_nrows (mcsg->cov), ==, ncm_matrix_nrows (cov));
  ncm_matrix_memcpy (mcsg->cov, cov);
  ncm_matrix_memcpy (mcsg->LLT, cov);
  ncm_matrix_cholesky_decomp (mcsg->LLT);
  mcsg->init = TRUE;
}

/**
 * ncm_mc_sampler_gauss_set_cov_variant:
 * @mcsg: a #NcmMCSamplerGauss.
 * @cov: a #GVariant.
 *
 * Sets the covariance given by the #GVariant @cov.
 * 
 */
void 
ncm_mc_sampler_gauss_set_cov_variant (NcmMCSamplerGauss *mcsg, GVariant *cov)
{
  ncm_matrix_set_from_variant (mcsg->cov, cov);
  ncm_matrix_memcpy (mcsg->LLT, mcsg->cov);
  ncm_matrix_cholesky_decomp (mcsg->LLT);
  mcsg->init = TRUE;
}

/**
 * ncm_mc_sampler_gauss_set_cov_data:
 * @mcsg: a #NcmMCSamplerGauss.
 * @cov: a #GVariant.
 *
 * Sets the covariance given by the double array @cov.
 * 
 */
void 
ncm_mc_sampler_gauss_set_cov_data (NcmMCSamplerGauss *mcsg, gdouble *cov)
{
  ncm_matrix_set_from_data (mcsg->cov, cov);
  ncm_matrix_memcpy (mcsg->LLT, mcsg->cov);
  ncm_matrix_cholesky_decomp (mcsg->LLT);
  mcsg->init = TRUE;
}

/**
 * ncm_mc_sampler_gauss_get_cov:
 * @mcsg: a #NcmMCSamplerGauss.
 *
 * Gets the covariance.
 * 
 * Returns: (transfer full): the covariance.
 */
NcmMatrix *
ncm_mc_sampler_gauss_get_cov (NcmMCSamplerGauss *mcsg)
{
  return ncm_matrix_ref (mcsg->cov);
}

/**
 * ncm_mc_sampler_gauss_set_cov_from_scale:
 * @mcsg: a #NcmMCSamplerGauss.
 *
 * Sets the covariance using the scale property of the parameters.
 * 
 */
void 
ncm_mc_sampler_gauss_set_cov_from_scale (NcmMCSamplerGauss *mcsg)
{
  NcmMCSampler *mcs = NCM_MC_SAMPLER (mcsg);
  guint i;
  g_assert (mcs->mset != NULL);

  ncm_matrix_set_identity (mcsg->cov);
  for (i = 0; i < mcsg->len; i++)
  {
    const gdouble scale = ncm_mset_fparam_get_scale (mcs->mset, i);
    ncm_matrix_set (mcsg->cov, i, i, scale * scale);
  }
  ncm_matrix_memcpy (mcsg->LLT, mcsg->cov);
  ncm_matrix_cholesky_decomp (mcsg->LLT);
  mcsg->init = TRUE;
}
