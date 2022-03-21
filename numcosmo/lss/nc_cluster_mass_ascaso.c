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
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcClusterMassAscaso, nc_cluster_mass_ascaso, NC_TYPE_CLUSTER_MASS);

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
  ascaso->M0             = 5.0e14;
  ascaso->Z0             = 0.0;
  ascaso->lnRichness_min = GSL_NEGINF;
  ascaso->lnRichness_max = GSL_POSINF;
}

static void
_nc_cluster_mass_ascaso_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (object);
  
  g_return_if_fail (NC_IS_CLUSTER_MASS_ASCASO (object));
  
  switch (prop_id)
  {
    case PROP_M0:
      ascaso->M0 = g_value_get_double (value);
      break;
    case PROP_Z0:
      ascaso->Z0 = g_value_get_double (value);
      break;
    case PROP_LNRICHNESS_MIN:
      ascaso->lnRichness_min = g_value_get_double (value);
      g_assert (ascaso->lnRichness_min < ascaso->lnRichness_max);
      break;
    case PROP_LNRICHNESS_MAX:
      ascaso->lnRichness_max = g_value_get_double (value);
      g_assert (ascaso->lnRichness_min < ascaso->lnRichness_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_ascaso_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (object);
  
  g_return_if_fail (NC_IS_CLUSTER_MASS_ASCASO (object));
  
  switch (prop_id)
  {
    case PROP_M0:
      g_value_set_double (value, ascaso->M0);
      break;
    case PROP_Z0:
      g_value_set_double (value, ascaso->Z0);
      break; 
    case PROP_LNRICHNESS_MIN:
      g_value_set_double (value, ascaso->lnRichness_min);
      break;
    case PROP_LNRICHNESS_MAX:
      g_value_set_double (value, ascaso->lnRichness_max);
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
static gboolean _nc_cluster_mass_ascaso_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_ascaso_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_ascaso_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);

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
                                                        11.0 * M_LN10, G_MAXDOUBLE, log (5.0) + 13.0 * M_LN10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

/**
   * NcClusterMassAscaso:Z0:
   *
   * Pivot redshift FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Z0,
                                   g_param_spec_double ("Z0",
                                                        NULL,
                                                        "Pivot redshift",
                                                        0.0, G_MAXDOUBLE, 10.0,
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
                                                        11.0 * M_LN10, G_MAXDOUBLE, log (5.0) + 13.0 * M_LN10,
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
                                                        11.0 * M_LN10, G_MAXDOUBLE, 16.0 * M_LN10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcClusterMassAscaso:MU_P0:
   *
   * Distribution's  bias in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P0, "mu_p0", "mup0",
                              0.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P0,
                              NCM_PARAM_TYPE_FIXED);
  
  /**
   * NcClusterMassAscaso:MU_P1:
   *
   * Distribution's slope with respect to the mass in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P1, "mu_p1", "mup1",
                              0.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P1,
                              NCM_PARAM_TYPE_FIXED);
  
  /**
   * NcClusterMassAscaso:MU_P2:
   *
   * Distribution's slope with respect to the redshift in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P2, "mu_p2", "mup2",
                              0.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P2,
                              NCM_PARAM_TYPE_FIXED);
  
  /**
   * NcClusterMassAscaso:sigma_P0:
   *
   * Distribution's bias in the standard deviation, $\sigma \in [10^{-4}, 10]$.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P0, "\\sigma_p0", "sigmap0",
                              1.0e-4,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:sigma_P1:
   *
   * Distribution's slope with respect to the mass in the standard deviation.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P1, "\\sigma_p1", "sigmap1",
                              0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P1,
                              NCM_PARAM_TYPE_FIXED);


/**
   * NcClusterMassAscaso:sigma_P2:
   *
   * Distribution's slope with respect to the redshift in the standard deviation.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P2, "\\sigma_p2", "sigmap2",
                              0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P2,
                              NCM_PARAM_TYPE_FIXED);
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->P              = &_nc_cluster_mass_ascaso_p;
  parent_class->intP           = &_nc_cluster_mass_ascaso_intp;
  parent_class->resample       = &_nc_cluster_mass_ascaso_resample;
  parent_class->P_limits       = &_nc_cluster_mass_ascaso_p_limits;
  parent_class->N_limits       = &_nc_cluster_mass_ascaso_n_limits;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 0;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static gdouble
_nc_cluster_mass_ascaso_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble logM, gdouble z, const gdouble *logM_obs, const gdouble *logM_obs_params)
{
  NcClusterMassAscaso *ascaso    = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble logM0            = log10 (ascaso->M0);
  const gdouble  z0              = ascaso->Z0;
  const gdouble logRichness_true = MU_P0 + MU_P1 * (logM - logM0) + MU_P2 * (log10 (1.0 + z) - log10 (1.0 + z0));
  const gdouble SIGMA            = SIGMA_P0 + SIGMA_P1 * (logM - logM0) + SIGMA_P2 * (log10 (1.0 + z) - log10 (1.0 + z0));
  const gdouble sqrt2_sigma      = M_SQRT2 * SIGMA;
  const gdouble x                = (logM_obs[0] - logRichness_true) / sqrt2_sigma;
  
  NCM_UNUSED (cosmo);
  
  return M_2_SQRTPI / (2.0 * M_SQRT2) * exp (-x * x) / (SIGMA);
}

static gdouble
_nc_cluster_mass_ascaso_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble logM, gdouble z)
{
  NcClusterMassAscaso *ascaso    = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble logM0            = log10 (ascaso->M0);
  const gdouble z0               = ascaso->Z0;
  const gdouble SIGMA            = SIGMA_P0 + SIGMA_P1 * (logM - logM0) + SIGMA_P2 * (log10 (1.0 + z) - log10 (1.0 + z0));
  const gdouble sqrt2_sigma      = M_SQRT2 * SIGMA;
  const gdouble logRichness_true = MU_P0 + MU_P1 * (logM - logM0) + MU_P2 * (log10 (1.0 + z) - log10 (1.0 + z0) );
  const gdouble x_min            = (logRichness_true - ascaso->lnRichness_min) / sqrt2_sigma;
  const gdouble x_max            = (logRichness_true - ascaso->lnRichness_max) / sqrt2_sigma;
  
  NCM_UNUSED (cosmo);
  
  if (x_max > 4.0)
    return -(erfc (x_min) - erfc (x_max)) / 2.0;
  else
    return (erf (x_min) - erf (x_max)) / 2.0;
}

static gboolean
_nc_cluster_mass_ascaso_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble logM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);
  gdouble logRichness_true;
  gdouble SIGMA;
  const gdouble logM0 = log10 (ascaso->M0);
  const gdouble z0    = ascaso->Z0;

  NCM_UNUSED (cosmo);
  
  ncm_rng_lock (rng);
  logRichness_true = MU_P0 + MU_P1 * (logM - logM0) + MU_P2 * (log10 (1.0 + z) - log10(1.0+z0));
  SIGMA            = SIGMA_P0 + SIGMA_P1 * (logM - logM0) + SIGMA_P2 * (log10 (1.0 + z) - log10 (1.0 + z0));
  lnM_obs[0]       =  logRichness_true + gsl_ran_gaussian (rng->r, SIGMA);
  ncm_rng_unlock (rng);
  
  return (lnM_obs[0] <= ascaso->lnRichness_max) && (lnM_obs[0] >= ascaso->lnRichness_min);
}

static void
_nc_cluster_mass_ascaso_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble mean          = lnM_obs[0] - MU_P0; /* - P2 * log10(1.0 + z);  FIX This!!!! What is the mean richeness? */
  const gdouble logRichnessl  = mean - 7.0 * SIGMA_P0;
  const gdouble logRichnessu  = mean + 7.0 * SIGMA_P0;
  
  NCM_UNUSED (cosmo);
  
  *lnM_lower = logRichnessl;
  *lnM_upper = logRichnessu;
  
  return;
}

static void
_nc_cluster_mass_ascaso_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (clusterm);
  const gdouble lnMl          = ascaso->lnRichness_min - 7.0 * SIGMA_P0;
  const gdouble lnMu          = ascaso->lnRichness_max + 7.0 * SIGMA_P0;
  
  NCM_UNUSED (cosmo);
  
  *lnM_lower = lnMl;
  *lnM_upper = lnMu;
  
  return;
}

