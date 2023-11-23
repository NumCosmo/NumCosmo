/***************************************************************************
 *            nc_cluster_mass_ascaso.c
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Mariana Penna Lima and Begoña Ascaso
 *  <pennalima@gmail.com>, <bego.ascaso.work@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima and Begoña Ascaso 2017 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_mass_ascaso
 * @title: NcClusterMassAscaso
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_ascaso.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcClusterMassAscasoPrivate
{
  gdouble M0;
  gdouble z0;
  gdouble lnM0;
  gdouble ln1pz0;
  gdouble lnR_max;
  gdouble lnR_min;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassAscaso, nc_cluster_mass_ascaso, NC_TYPE_CLUSTER_MASS)

#define VECTOR (NCM_MODEL (ascaso)->params)
#define MU_P0 (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P0))
#define MU_P1 (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P1))
#define MU_P2 (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P2))
#define SIGMA_P0  (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_ASCASO_SIGMA_P0))
#define SIGMA_P1  (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_ASCASO_SIGMA_P1))
#define SIGMA_P2  (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_ASCASO_SIGMA_P2))

enum
{
  PROP_0,
  PROP_M0,
  PROP_Z0,
  PROP_LNRICHNESS_MIN,
  PROP_LNRICHNESS_MAX,
  PROP_SIZE,
};

static void
nc_cluster_mass_ascaso_init (NcClusterMassAscaso *ascaso)
{
  NcClusterMassAscasoPrivate * const self = ascaso->priv = nc_cluster_mass_ascaso_get_instance_private (ascaso);

  self->M0      = 0.0;
  self->z0      = 0.0;
  self->lnM0    = 0.0;
  self->ln1pz0  = 0.0;
  self->lnR_min = GSL_NEGINF;
  self->lnR_max = GSL_POSINF;
}

static void
_nc_cluster_mass_ascaso_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (object);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;

  g_return_if_fail (NC_IS_CLUSTER_MASS_ASCASO (object));

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
      g_assert (self->lnR_min < self->lnR_max);
      break;
    case PROP_LNRICHNESS_MAX:
      self->lnR_max = g_value_get_double (value);
      g_assert (self->lnR_min < self->lnR_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_ascaso_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (object);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;

  g_return_if_fail (NC_IS_CLUSTER_MASS_ASCASO (object));

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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_ascaso_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_ascaso_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_ascaso_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_ascaso_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_ascaso_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_ascaso_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_ascaso_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_ascaso_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_ascaso_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
static gdouble _nc_cluster_mass_ascaso_volume (NcClusterMass *clusterm);

static void _nc_cluster_mass_ascaso_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, GArray *res);

static void
nc_cluster_mass_ascaso_class_init (NcClusterMassAscasoClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_cluster_mass_ascaso_set_property;
  model_class->get_property = &_nc_cluster_mass_ascaso_get_property;
  object_class->finalize    = &_nc_cluster_mass_ascaso_finalize;

  ncm_model_class_set_name_nick (model_class, "Ascaso Ln-normal richness distribution", "Ascaso");
  ncm_model_class_add_params (model_class, 6, 0, PROP_SIZE);


  /**
   * NcClusterMassAscaso:M0:
   *
   * Pivot mass FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Pivot mass",
                                                        11.0 * M_LN10, G_MAXDOUBLE, 3.0e14 / 0.71,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

/**
 * NcClusterMassAscaso:Z0:
 *
 * Pivot redshift FIXME Set correct values (limits)
 */
  g_object_class_install_property (object_class,
                                   PROP_Z0,
                                   g_param_spec_double ("z0",
                                                        NULL,
                                                        "Pivot redshift",
                                                        0.0, G_MAXDOUBLE, 0.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



  /**
   * NcClusterMassAscaso:lnRichness_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MIN,
                                   g_param_spec_double ("lnRichness-min",
                                                        NULL,
                                                        "Minimum LnRichness",
                                                        0.0, G_MAXDOUBLE, M_LN10 * 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassAscaso:lnRichness_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MAX,
                                   g_param_spec_double ("lnRichness-max",
                                                        NULL,
                                                        "Maximum LnRichness",
                                                        0.0, G_MAXDOUBLE,  M_LN10 * 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassAscaso:MU_P0:
   *
   * Distribution's  bias in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P0, "mu_p0", "mup0",
                              0.0,  6.0, 1.0e-1,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:MU_P1:
   *
   * Distribution's slope with respect to the mass in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P1, "mu_p1", "mup1",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:MU_P2:
   *
   * Distribution's slope with respect to the redshift in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P2, "mu_p2", "mup2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:sigma_P0:
   *
   * Distribution's bias in the standard deviation, $\sigma \in [10^{-4}, 10]$.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P0, "\\sigma_p0", "sigmap0",
                              1.0e-4, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:sigma_P1:
   *
   * Distribution's slope with respect to the mass in the standard deviation.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P1, "\\sigma_p1", "sigmap1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P1,
                              NCM_PARAM_TYPE_FIXED);


/**
 * NcClusterMassAscaso:sigma_P2:
 *
 * Distribution's slope with respect to the redshift in the standard deviation.
 *
 */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P2, "\\sigma_p2", "sigmap2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P2,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->P              = &_nc_cluster_mass_ascaso_p;
  parent_class->intP           = &_nc_cluster_mass_ascaso_intp;
  parent_class->intP_bin       = &_nc_cluster_mass_ascaso_intp_bin;
  parent_class->resample       = &_nc_cluster_mass_ascaso_resample;
  parent_class->P_limits       = &_nc_cluster_mass_ascaso_p_limits;
  parent_class->P_bin_limits   = &_nc_cluster_mass_ascaso_p_bin_limits;
  parent_class->N_limits       = &_nc_cluster_mass_ascaso_n_limits;
  parent_class->volume         = &_nc_cluster_mass_ascaso_volume;
  parent_class->P_vec_z_lnMobs = &_nc_cluster_mass_ascaso_p_vec_z_lnMobs;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 0;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static void
_nc_cluster_mass_ascaso_lnR_sigma (NcClusterMass *clusterm, const gdouble lnM, const gdouble z, gdouble *lnR, gdouble *sigma)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (clusterm);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;
  const gdouble DlnM                      = lnM - self->lnM0;
  const gdouble Dln1pz                    = log1p (z) - self->ln1pz0;

  lnR[0]   = MU_P0    + MU_P1    * DlnM + MU_P2    * Dln1pz;
  sigma[0] = SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;
}

static gdouble
_nc_cluster_mass_ascaso_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  /*NcClusterMassAscaso *ascaso    = NC_CLUSTER_MASS_ASCASO (clusterm);*/
  gdouble lnR_true, sigma;

  _nc_cluster_mass_ascaso_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x = (lnM_obs[0] - lnR_true) / sigma;

    return 1.0 / (ncm_c_sqrt_2pi () * sigma) * exp (-0.5 * x * x);
  }
}

static gdouble
_nc_cluster_mass_ascaso_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (clusterm);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;
  gdouble lnR_true, sigma;

  _nc_cluster_mass_ascaso_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x_min = (lnR_true - self->lnR_min) / (M_SQRT2 * sigma);
    const gdouble x_max = (lnR_true - self->lnR_max) / (M_SQRT2 * sigma);

    if (x_max > 4.0)
      return -(erfc (x_min) - erfc (x_max)) / 2.0;
    else
      return (erf (x_min) - erf (x_max)) / 2.0;
  }
}

static gdouble
_nc_cluster_mass_ascaso_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  /*NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);*/
  gdouble lnR_true, sigma;

  _nc_cluster_mass_ascaso_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x_min = (lnR_true - lnM_obs_lower[0]) / (M_SQRT2 * sigma);
    const gdouble x_max = (lnR_true - lnM_obs_upper[0]) / (M_SQRT2 * sigma);

    if (x_max > 4.0)
      return -(erfc (x_min) - erfc (x_max)) / 2.0;
    else
      return (erf (x_min) - erf (x_max)) / 2.0;
  }
}

static gboolean
_nc_cluster_mass_ascaso_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (clusterm);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;
  gdouble lnR_true, sigma;

  _nc_cluster_mass_ascaso_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  ncm_rng_lock (rng);
  lnM_obs[0] = lnR_true + gsl_ran_gaussian (rng->r, sigma);
  ncm_rng_unlock (rng);

  return (lnM_obs[0] <= self->lnR_max) && (lnM_obs[0] >= self->lnR_min);
}

static void
_nc_cluster_mass_ascaso_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble mean          = lnM_obs[0] - MU_P0; /* - P2 * log10(1.0 + z);  FIX This!!!! What is the mean richeness? */
  const gdouble logRichnessl  = M_LN10 * log10 (1e13);
  const gdouble logRichnessu  = M_LN10 * log10 (1e15);

  *lnM_lower = logRichnessl;
  *lnM_upper = logRichnessu;

  return;
}

static void
_nc_cluster_mass_ascaso_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble lnMl          = M_LN10 * log10 (1e13);
  const gdouble lnMu          = M_LN10 * log10 (1e15);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;
}

static void
_nc_cluster_mass_ascaso_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble lnMl          =  M_LN10 * log10 (1e13);
  const gdouble lnMu          =  M_LN10 * log10 (1e15);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static gdouble
_nc_cluster_mass_ascaso_volume (NcClusterMass *clusterm)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (clusterm);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;

  return (self->lnR_max - self->lnR_min);
}

static void
_nc_cluster_mass_ascaso_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, GArray *res)
{
  NcClusterMassAscaso *ascaso             = NC_CLUSTER_MASS_ASCASO (clusterm);
  NcClusterMassAscasoPrivate * const self = ascaso->priv;

  const gdouble *lnM_obs_ptr = ncm_matrix_const_data (lnM_obs);
  const gdouble *z_ptr       = ncm_vector_const_data (z);
  const guint tda            = ncm_matrix_tda (lnM_obs);
  const guint sz             = ncm_vector_stride (z);
  const guint len            = ncm_vector_len (z);
  const gdouble DlnM         = lnM - self->lnM0;
  const gdouble sqrt_2pi     = ncm_c_sqrt_2pi ();
  const gdouble lnR_pre      = MU_P0    + MU_P1    * DlnM;
  const gdouble sigma_pre    = SIGMA_P0 + SIGMA_P1 * DlnM;
  const gdouble mu_p2        = MU_P2;
  const gdouble sigma_p2     = SIGMA_P2;
  gdouble *res_ptr           = &g_array_index (res, gdouble, 0);
  guint i;

  if ((tda == 1) && (sz == 1))
  {
    for (i = 0; i < len; i++)
    {
      const gdouble Dln1pz = log1p (z_ptr[i]) - self->ln1pz0;
      const gdouble lnR    = lnR_pre + mu_p2 * Dln1pz;
      const gdouble sigma  = sigma_pre + sigma_p2 * Dln1pz;

      const gdouble x = (lnM_obs_ptr[i] - lnR) / sigma;

      res_ptr[i] = exp (-0.5 * x * x) / (sqrt_2pi * sigma);
    }
  }
  else
  {
    for (i = 0; i < len; i++)
    {
      const gdouble Dln1pz = log1p (z_ptr[i * sz]) - self->ln1pz0;
      const gdouble lnR    = lnR_pre + mu_p2 * Dln1pz;
      const gdouble sigma  = sigma_pre + sigma_p2 * Dln1pz;

      const gdouble x = (lnM_obs_ptr[i * tda] - lnR) / sigma;

      res_ptr[i] = exp (-0.5 * x * x) / (sqrt_2pi * sigma);
    }
  }
}

