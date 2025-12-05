/***************************************************************************
 *            nc_cluster_mass_richness.c
 *
 *  Tue May 27 14:05:00 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2025 <vitenti@uel.br>
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
 * NcClusterMassRichness:
 *
 * Abstract class for cluster mass-richness observable relations.
 *
 * This abstract class implements the truncated lognormal distribution for
 * cluster mass-richness relations. Subclasses must implement the virtual
 * functions mu() and sigma() that define the mean and standard deviation
 * of the log-richness as a function of mass and redshift.
 *
 * The probability distribution is given by:
 * $$
 * P(\ln \lambda | M, z) = \frac{1}{\sqrt{2\pi}\sigma} \exp\left[-\frac{(\ln \lambda - \mu)^2}{2\sigma^2}\right]
 * $$
 * where $\mu = \mu(M, z)$ and $\sigma = \sigma(M, z)$ are model-specific functions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_richness.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_c.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"
#include "math/ncm_matrix.h"

typedef struct _NcClusterMassRichnessPrivate
{
  gdouble M0;
  gdouble z0;
  gdouble lnM0;
  gdouble ln1pz0;
  gdouble lnR_min;
  gdouble lnR_max;
  gboolean sample_full_dist;
} NcClusterMassRichnessPrivate;

enum
{
  PROP_0,
  PROP_M0,
  PROP_Z0,
  PROP_LNRICHNESS_MIN,
  PROP_LNRICHNESS_MAX,
  PROP_SAMPLE_FULL_DIST,
  PROP_LEN,
};

#define VECTOR (NCM_MODEL (mr))
#define CUT    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_CUT))

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcClusterMassRichness, nc_cluster_mass_richness, NC_TYPE_CLUSTER_MASS)

static void
nc_cluster_mass_richness_init (NcClusterMassRichness *mr)
{
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  self->M0               = 0.0;
  self->z0               = 0.0;
  self->lnM0             = 0.0;
  self->ln1pz0           = 0.0;
  self->lnR_min          = 0.0;
  self->lnR_max          = 0.0;
  self->sample_full_dist = TRUE;
}

static void
nc_cluster_mass_richness_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassRichness *mr                 = NC_CLUSTER_MASS_RICHNESS (object);
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  g_return_if_fail (NC_IS_CLUSTER_MASS_RICHNESS (object));

  switch (prop_id)
  {
    case PROP_M0:
      self->M0   = g_value_get_double (value);
      self->lnM0 = log (self->M0);
      break;
    case PROP_Z0:
      self->z0     = g_value_get_double (value);
      self->ln1pz0 = log1p (self->z0);
      break;
    case PROP_LNRICHNESS_MIN:
      self->lnR_min = g_value_get_double (value);
      break;
    case PROP_LNRICHNESS_MAX:
      self->lnR_max = g_value_get_double (value);
      break;
    case PROP_SAMPLE_FULL_DIST:
      self->sample_full_dist = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_richness_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassRichness *mr                 = NC_CLUSTER_MASS_RICHNESS (object);
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  g_return_if_fail (NC_IS_CLUSTER_MASS_RICHNESS (object));

  switch (prop_id)
  {
    case PROP_M0:
      g_value_set_double (value, self->M0);
      break;
    case PROP_Z0:
      g_value_set_double (value, self->z0);
      break;
    case PROP_LNRICHNESS_MIN:
      g_value_set_double (value, self->lnR_min);
      break;
    case PROP_LNRICHNESS_MAX:
      g_value_set_double (value, self->lnR_max);
      break;
    case PROP_SAMPLE_FULL_DIST:
      g_value_set_boolean (value, self->sample_full_dist);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_richness_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_richness_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_richness_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_richness_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
static void _nc_cluster_mass_richness_mu_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z, gdouble *mu, gdouble *sigma);

static gdouble _nc_cluster_mass_richness_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_richness_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_richness_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_richness_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_richness_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_richness_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_richness_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
static gdouble _nc_cluster_mass_richness_volume (NcClusterMass *clusterm);
static void _nc_cluster_mass_richness_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, NcmVector *res);

static void
nc_cluster_mass_richness_class_init (NcClusterMassRichnessClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcClusterMassClass *cm_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  model_class->set_property = &nc_cluster_mass_richness_set_property;
  model_class->get_property = &nc_cluster_mass_richness_get_property;
  object_class->finalize    = &nc_cluster_mass_richness_finalize;

  ncm_model_class_set_name_nick (model_class, "Cluster mass-richness relation abstract class", "NcClusterMassRichness");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcClusterMassRichness:M0:
   *
   * Pivot mass in the richness-mass scaling relation.
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Pivot mass",
                                                        11.0 * M_LN10, G_MAXDOUBLE, 3.0e14 / 0.71,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassRichness:z0:
   *
   * Pivot redshift in the richness-mass scaling relation.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z0,
                                   g_param_spec_double ("z0",
                                                        NULL,
                                                        "Pivot redshift",
                                                        0.0, G_MAXDOUBLE, 0.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassRichness:lnRichness-min:
   *
   * Minimum logarithm (base e) of richness for cluster selection.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MIN,
                                   g_param_spec_double ("lnRichness-min",
                                                        NULL,
                                                        "Minimum LnRichness",
                                                        0.0, G_MAXDOUBLE, M_LN10 * 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassRichness:lnRichness-max:
   *
   * Maximum logarithm (base e) of richness for cluster selection.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MAX,
                                   g_param_spec_double ("lnRichness-max",
                                                        NULL,
                                                        "Maximum LnRichness",
                                                        0.0, G_MAXDOUBLE, M_LN10 * 2.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassRichness:sample-full-dist:
   *
   * Controls the sampling strategy for richness values:
   *
   * When TRUE (default): Samples are drawn from the full (untruncated) Gaussian
   * distribution. This means some samples will fall outside the allowed richness
   * range [CUT, lnR_max] and should be discarded by the caller. This approach
   * correctly models the selection effect where clusters with low richness are
   * not observed.
   *
   * When FALSE: Samples are drawn from the truncated Gaussian distribution using
   * tail sampling. All samples are guaranteed to fall within the allowed range.
   * This is useful when you want to ensure every input cluster produces exactly
   * one output cluster, but does not model the selection effect.
   */
  g_object_class_install_property (object_class,
                                   PROP_SAMPLE_FULL_DIST,
                                   g_param_spec_boolean ("sample-full-dist",
                                                         NULL,
                                                         "Whether to sample from the full (untruncated) distribution",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassRichness:CUT:
   *
   * Cut in ln-richness, i.e., clusters with $\ln(\lambda) < $ CUT are excluded.
   * Default is 0.0, corresponding to a richness cut of 1.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_RICHNESS_CUT, "CUT", "cut",
                              -2.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_RICHNESS_DEFAULT_CUT,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  /* Set virtual methods for NcClusterMass */
  cm_class->P               = &_nc_cluster_mass_richness_p;
  cm_class->intP            = &_nc_cluster_mass_richness_intp;
  cm_class->intP_bin        = &_nc_cluster_mass_richness_intp_bin;
  cm_class->resample        = &_nc_cluster_mass_richness_resample;
  cm_class->P_limits        = &_nc_cluster_mass_richness_p_limits;
  cm_class->P_bin_limits    = &_nc_cluster_mass_richness_p_bin_limits;
  cm_class->N_limits        = &_nc_cluster_mass_richness_n_limits;
  cm_class->volume          = &_nc_cluster_mass_richness_volume;
  cm_class->P_vec_z_lnMobs  = &_nc_cluster_mass_richness_p_vec_z_lnMobs;
  cm_class->_obs_len        = 1;
  cm_class->_obs_params_len = 1;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);

  /* Set virtual methods for this class */
  klass->mu       = &_nc_cluster_mass_richness_mu;
  klass->sigma    = &_nc_cluster_mass_richness_sigma;
  klass->mu_sigma = &_nc_cluster_mass_richness_mu_sigma;
}

static gdouble
_nc_cluster_mass_richness_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  g_error ("_nc_cluster_mass_richness_mu: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (mr));

  return 0.0;
}

static gdouble
_nc_cluster_mass_richness_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  g_error ("_nc_cluster_mass_richness_sigma: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (mr));

  return 0.0;
}

static void
_nc_cluster_mass_richness_mu_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z, gdouble *mu, gdouble *sigma)
{
  *mu    = nc_cluster_mass_richness_mu (mr, lnM, z);
  *sigma = nc_cluster_mass_richness_sigma (mr, lnM, z);
}

static gdouble
_nc_cluster_mass_richness_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  NcClusterMassRichness *mr = NC_CLUSTER_MASS_RICHNESS (clusterm);
  gdouble mu, sigma;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);

  if (lnM_obs[0] < CUT)
  {
    return 0.0;
  }
  else
  {
    const gdouble sigma_lnR_cat = lnM_obs_params[0];
    const gdouble sigma_total   = sqrt (sigma * sigma + sigma_lnR_cat * sigma_lnR_cat);
    const gdouble x             = (lnM_obs[0] - mu) / sigma_total;

    return 1.0 / (ncm_c_sqrt_2pi () * sigma_total) * exp (-0.5 * x * x);
  }
}

static gdouble
_nc_cluster_mass_richness_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassRichness *mr                 = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);
  gdouble mu, sigma;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);
  {
    /* Note: intP does not use catalog sigma - it represents the selection function */
    const gdouble x_cut = (mu - CUT) / (M_SQRT2 * sigma);
    const gdouble x_max = (mu - self->lnR_max) / (M_SQRT2 * sigma);

    if ((fabs (x_max) > 4.0) || (fabs (x_cut) > 4.0))
      return -(erfc (x_cut) - erfc (x_max)) / 2.0;
    else
      return (erf (x_cut) - erf (x_max)) / 2.0;
  }
}

static gdouble
_nc_cluster_mass_richness_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  NcClusterMassRichness *mr = NC_CLUSTER_MASS_RICHNESS (clusterm);

  if ((lnM_obs_lower[0] < CUT) && (lnM_obs_upper[0] < 0.0))
  {
    return 0.0;
  }
  else
  {
    gdouble mu, sigma;
    const gdouble sigma_lnR_cat = lnM_obs_params[0];

    nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);
    sigma = sqrt (sigma * sigma + sigma_lnR_cat * sigma_lnR_cat);
    {
      const gdouble x_min = (mu - lnM_obs_lower[0]) / (M_SQRT2 * sigma);
      const gdouble x_max = (mu - lnM_obs_upper[0]) / (M_SQRT2 * sigma);
      const gdouble x_cut = (mu - CUT) / (M_SQRT2 * sigma);
      gdouble intp_bin    = 0.0;

      if ((lnM_obs_lower[0] < 0.0) && (lnM_obs_upper[0] >= 0.0))
      {
        if ((fabs (x_max) > 4.0) || (fabs (x_cut) > 4.0))
          return -(erfc (x_cut) - erfc (x_max)) / 2.0;
        else
          return (erf (x_cut) - erf (x_max)) / 2.0;
      }
      else
      {
        if ((fabs (x_max) > 4.0) || (fabs (x_min) > 4.0))
          intp_bin = -(erfc (x_min) - erfc (x_max)) / 2.0;
        else
          intp_bin = (erf (x_min) - erf (x_max)) / 2.0;

        if (intp_bin < 0.0)
          return 0.0;
        else
          return intp_bin;
      }
    }
  }
}

static gboolean
_nc_cluster_mass_richness_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassRichness *mr                 = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);
  gdouble mu, sigma, sigma_total;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);

  sigma_total = sqrt (sigma * sigma + lnM_obs_params[0] * lnM_obs_params[0]);

  ncm_rng_lock (rng);

  if (self->sample_full_dist)
  {
    lnM_obs[0] = ncm_rng_gaussian_gen (rng, mu, sigma_total);
  }
  else
  {
    lnM_obs[0]  = ncm_rng_gaussian_tail_gen (rng, CUT - mu, sigma_total);
    lnM_obs[0] += mu;
  }

  ncm_rng_unlock (rng);

  return (lnM_obs[0] <= self->lnR_max) && (lnM_obs[0] >= CUT);
}

static void
_nc_cluster_mass_richness_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  const gdouble lnMl = M_LN10 * 13.0;
  const gdouble lnMu = M_LN10 * 16.0;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_richness_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  const gdouble lnMl = M_LN10 * 13.0;
  const gdouble lnMu = M_LN10 * 16.0;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;
}

static void
_nc_cluster_mass_richness_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  const gdouble lnMl = M_LN10 * 13.0;
  const gdouble lnMu = M_LN10 * 16.0;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static gdouble
_nc_cluster_mass_richness_volume (NcClusterMass *clusterm)
{
  NcClusterMassRichness *mr                 = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  return (self->lnR_max - CUT);
}

static void
_nc_cluster_mass_richness_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, NcmVector *res)
{
  NcClusterMassRichness *mr  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  const gdouble *lnM_obs_ptr = ncm_matrix_const_data (lnM_obs);
  const gdouble *z_ptr       = ncm_vector_const_data (z);
  const guint tda            = ncm_matrix_tda (lnM_obs);
  const guint sz             = ncm_vector_stride (z);
  const guint len            = ncm_vector_len (z);
  const gdouble sqrt_2pi     = ncm_c_sqrt_2pi ();
  const gdouble cut          = CUT;
  gdouble *res_ptr           = ncm_vector_ptr (res, 0);
  guint i;

  g_assert_cmpuint (ncm_vector_stride (res), ==, 1);

  if ((tda == 1) && (sz == 1))
  {
    const gdouble *params_ptr = ncm_matrix_const_data (lnM_obs_params);
    const guint params_tda    = ncm_matrix_tda (lnM_obs_params);

    for (i = 0; i < len; i++)
    {
      gdouble lnR, sigma;
      const gdouble sigma_lnR_cat = params_ptr[i * params_tda];

      nc_cluster_mass_richness_mu_sigma (mr, lnM, z_ptr[i], &lnR, &sigma);
      {
        const gdouble sigma_total = sqrt (sigma * sigma + sigma_lnR_cat * sigma_lnR_cat);
        const gdouble x           = (lnM_obs_ptr[i] - lnR) / sigma_total;

        if (lnM_obs_ptr[i] < cut)
          res_ptr[i] = 0.0;
        else
          res_ptr[i] = 1.0 * exp (-0.5 * x * x) / (sqrt_2pi * sigma_total);
      }
    }
  }
  else
  {
    const gdouble *params_ptr = ncm_matrix_const_data (lnM_obs_params);
    const guint params_tda    = ncm_matrix_tda (lnM_obs_params);

    for (i = 0; i < len; i++)
    {
      gdouble lnR, sigma;
      const gdouble sigma_lnR_cat = params_ptr[i * params_tda];

      nc_cluster_mass_richness_mu_sigma (mr, lnM, z_ptr[i * sz], &lnR, &sigma);
      {
        const gdouble sigma_total = sqrt (sigma * sigma + sigma_lnR_cat * sigma_lnR_cat);
        const gdouble x           = (lnM_obs_ptr[i * tda] - lnR) / sigma_total;

        if (lnM_obs_ptr[i * tda] < cut)
          res_ptr[i] = 0.0;
        else
          res_ptr[i] = 1.0 * exp (-0.5 * x * x) / (sqrt_2pi * sigma_total);
      }
    }
  }
}

/**
 * nc_cluster_mass_richness_mu:
 * @mr: a #NcClusterMassRichness
 * @lnM: natural logarithm of the mass
 * @z: redshift
 *
 * Computes the mean of the log-richness distribution as a function of
 * mass and redshift.
 *
 * Returns: the mean log-richness $\mu(\ln M, z)$
 */
gdouble
nc_cluster_mass_richness_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  return NC_CLUSTER_MASS_RICHNESS_GET_CLASS (mr)->mu (mr, lnM, z);
}

/**
 * nc_cluster_mass_richness_sigma:
 * @mr: a #NcClusterMassRichness
 * @lnM: natural logarithm of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the log-richness distribution as a
 * function of mass and redshift.
 *
 * Returns: the standard deviation $\sigma(\ln M, z)$
 */
gdouble
nc_cluster_mass_richness_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  return NC_CLUSTER_MASS_RICHNESS_GET_CLASS (mr)->sigma (mr, lnM, z);
}

/**
 * nc_cluster_mass_richness_mu_sigma:
 * @mr: a #NcClusterMassRichness
 * @lnM: natural logarithm of the mass
 * @z: redshift
 * @mu: (out): location to store the mean log-richness
 * @sigma: (out): location to store the standard deviation
 *
 * Computes both the mean and standard deviation of the log-richness distribution
 * simultaneously. This can be more efficient than calling mu() and sigma()
 * separately when subclasses share intermediate computations.
 */
void
nc_cluster_mass_richness_mu_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z, gdouble *mu, gdouble *sigma)
{
  NC_CLUSTER_MASS_RICHNESS_GET_CLASS (mr)->mu_sigma (mr, lnM, z, mu, sigma);
}

/**
 * nc_cluster_mass_richness_get_cut:
 * @mr: a #NcClusterMassRichness
 *
 * Gets the cut in richness.
 *
 * Returns: the cut in richness.
 */
gdouble
nc_cluster_mass_richness_get_cut (NcClusterMassRichness *mr)
{
  return CUT;
}

/**
 * nc_cluster_mass_richness_set_sample_full_dist:
 * @mr: a #NcClusterMassRichness
 * @on: whether to sample from the full distribution
 *
 * Sets the sampling strategy for richness values.
 *
 * When @on is TRUE: Samples are drawn from the full Gaussian distribution.
 * Out-of-range samples (below CUT or above lnR_max) should be discarded,
 * modeling the selection effect where low-richness clusters are not observed.
 *
 * When @on is FALSE: Samples are drawn from the truncated distribution,
 * guaranteeing all values fall within [CUT, lnR_max].
 */
void
nc_cluster_mass_richness_set_sample_full_dist (NcClusterMassRichness *mr, gboolean on)
{
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  self->sample_full_dist = on;
}

/**
 * nc_cluster_mass_richness_get_sample_full_dist:
 * @mr: a #NcClusterMassRichness
 *
 * Gets the current sampling strategy.
 *
 * Returns: TRUE if sampling from the full distribution (may produce
 *          out-of-range values), FALSE if sampling from the truncated
 *          distribution (all values within allowed range).
 */
gboolean
nc_cluster_mass_richness_get_sample_full_dist (NcClusterMassRichness *mr)
{
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  return self->sample_full_dist;
}

/**
 * nc_cluster_mass_richness_get_mean:
 * @mr: a #NcClusterMassRichness
 * @lnM: natural logarithm of the mass
 * @z: redshift
 *
 * Computes the mean of the truncated richness distribution.
 *
 * Returns: the mean of the truncated distribution.
 */
gdouble
nc_cluster_mass_richness_get_mean (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  gdouble lnR_mean, lnR_sigma;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &lnR_mean, &lnR_sigma);

  return nc_cluster_mass_richness_compute_truncated_mean (mr, lnR_mean, lnR_sigma);
}

/**
 * nc_cluster_mass_richness_get_std:
 * @mr: a #NcClusterMassRichness
 * @lnM: natural logarithm of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the truncated richness distribution.
 *
 * Returns: the standard deviation of the truncated distribution.
 */
gdouble
nc_cluster_mass_richness_get_std (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  gdouble lnR_mean, lnR_sigma;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &lnR_mean, &lnR_sigma);

  return nc_cluster_mass_richness_compute_truncated_std (mr, lnR_mean, lnR_sigma);
}

/**
 * nc_cluster_mass_richness_compute_truncated_mean:
 * @mr: a #NcClusterMassRichness
 * @lnR_mean: mean of the untruncated log-richness distribution
 * @lnR_sigma: standard deviation of the untruncated log-richness distribution
 *
 * Computes the mean of the truncated log-richness distribution given the
 * untruncated mean and standard deviation. The distribution is truncated
 * at the cut threshold.
 *
 * Returns: the mean of the truncated distribution.
 */
gdouble
nc_cluster_mass_richness_compute_truncated_mean (NcClusterMassRichness *mr, gdouble lnR_mean, gdouble lnR_sigma)
{
  const gdouble A               = (CUT - lnR_mean) / lnR_sigma;
  const gdouble B               = (1.0 / ncm_c_sqrt_2pi ()) * exp (-0.5 * A * A);
  const gdouble C               = 0.5 * erfc (A / M_SQRT2);
  const gdouble mean_correction = (lnR_sigma * B / C);

  return lnR_mean + mean_correction;
}

/**
 * nc_cluster_mass_richness_compute_truncated_std:
 * @mr: a #NcClusterMassRichness
 * @lnR_mean: mean of the untruncated log-richness distribution
 * @lnR_sigma: standard deviation of the untruncated log-richness distribution
 *
 * Computes the standard deviation of the truncated log-richness distribution
 * given the untruncated mean and standard deviation. The distribution is
 * truncated at the cut threshold.
 *
 * Returns: the standard deviation of the truncated distribution.
 */
gdouble
nc_cluster_mass_richness_compute_truncated_std (NcClusterMassRichness *mr, gdouble lnR_mean, gdouble lnR_sigma)
{
  const gdouble A              = (CUT - lnR_mean) / lnR_sigma;
  const gdouble B              = (1.0 / ncm_c_sqrt_2pi ()) * exp (-0.5 * A * A);
  const gdouble C              = 0.5 * erfc (A / M_SQRT2);
  const gdouble std_correction = pow (1.0 + (A * B / C) - (B / C) * (B / C), 0.5);

  return lnR_sigma * std_correction;
}

/**
 * nc_cluster_mass_richness_lnM0:
 * @mr: a #NcClusterMassRichness
 *
 * Gets the natural logarithm of the pivot mass.
 *
 * Returns: $\ln M_0$
 */
gdouble
nc_cluster_mass_richness_lnM0 (NcClusterMassRichness *mr)
{
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  return self->lnM0;
}

/**
 * nc_cluster_mass_richness_ln1pz0:
 * @mr: a #NcClusterMassRichness
 *
 * Gets the natural logarithm of (1 + pivot redshift).
 *
 * Returns: $\ln(1 + z_0)$
 */
gdouble
nc_cluster_mass_richness_ln1pz0 (NcClusterMassRichness *mr)
{
  NcClusterMassRichnessPrivate * const self = nc_cluster_mass_richness_get_instance_private (mr);

  return self->ln1pz0;
}

