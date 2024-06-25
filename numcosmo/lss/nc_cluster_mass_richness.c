/***************************************************************************
 *            nc_cluster_mass_richness.c
 *
 *  Fri June 22 17:48:55 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_mass_richness
 * @title: NcClusterMassRichness
 * @short_description: Cluster mass-richenss ln-normal distribution with a mass-redshift dependent bias and richness measurement erros.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_richness.h"
#include "math/ncm_cfg.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcClusterMassRichnessPrivate
{
  gdouble M0;
  gdouble lnM0;
  gdouble richnessobs_min;
  gdouble richnessobs_max;
  
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassRichness, nc_cluster_mass_richness, NC_TYPE_CLUSTER_MASS)

#define VECTOR   (NCM_MODEL (richness))
#define MU       (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_MU))
#define MU_M     (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_MU_M))
#define MU_Z     (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_MU_Z))
#define MU_MZ    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_MU_MZ))
#define SIGMA_0  (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_SIGMA_0))

enum
{
  PROP_0,
  PROP_M0,
  PROP_LNM0,
  PROP_RICHNESSOBS_MIN,
  PROP_RICHNESSOBS_MAX,
  PROP_SIZE,
};

#define _NC_CLUSTER_MASS_RICHNESS_DEFAULT_INT_KEY 6

static void
nc_cluster_mass_richness_init (NcClusterMassRichness *richness)
{
  NcClusterMassRichnessPrivate * const self = richness->priv = nc_cluster_mass_richness_get_instance_private (richness);


  self->richnessobs_min = 0;
  self->richnessobs_max = GSL_POSINF;
  self->M0      = 0.0;
  self->lnM0    = 0.0;
}

static void
_nc_cluster_mass_richness_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassRichness *richness             = NC_CLUSTER_MASS_RICHNESS (object);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  g_return_if_fail (NC_IS_CLUSTER_MASS_RICHNESS (object));

  switch (prop_id)
  {
    case PROP_RICHNESSOBS_MIN:
      self->richnessobs_min = g_value_get_double (value);
      g_assert (self->richnessobs_min < self->richnessobs_max);
      g_assert (self->richnessobs_min > 0);
      break;
    case PROP_RICHNESSOBS_MAX:
      self->richnessobs_max = g_value_get_double (value);
      g_assert (self->richnessobs_min < self->richnessobs_max);
      break;
    case PROP_M0:
      self->M0   = g_value_get_double (value);
      self->lnM0 = log (self->M0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_richness_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (object);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  g_return_if_fail (NC_IS_CLUSTER_MASS_RICHNESS (object));

  switch (prop_id)
  {
    case PROP_M0:
      g_value_set_double (value, self->M0);
      break;
    case PROP_RICHNESSOBS_MIN:
      g_value_set_double (value, self->richnessobs_min);
      break;
    case PROP_RICHNESSOBS_MAX:
      g_value_set_double (value, self->richnessobs_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
  }

static void
_nc_cluster_mass_richness_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_richness_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_richness_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_richness_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_richness_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_richness_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_richness_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_richness_p_limits_bin (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_richness_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);

static void
nc_cluster_mass_richness_class_init (NcClusterMassRichnessClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_cluster_mass_richness_finalize;

  model_class->set_property = &_nc_cluster_mass_richness_set_property;
  model_class->get_property = &_nc_cluster_mass_richness_get_property;

  ncm_model_class_set_name_nick (model_class, "richness-mass distribution", "Richness_Mass");
  ncm_model_class_add_params (model_class, 5, 0, PROP_SIZE);

  /**
   * NcClusterMassRichness:richnessobs_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_RICHNESSOBS_MIN,
                                   g_param_spec_double ("richnessobs-min",
                                                        NULL,
                                                        "Minimum observed richness",
                                                        1, G_MAXDOUBLE, 50,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassRichness:richnessobs_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_RICHNESSOBS_MAX,
                                   g_param_spec_double ("richnessobs-max",
                                                        NULL,
                                                        "Maximum observed richness",
                                                        80, G_MAXDOUBLE, 1000,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

 /**
   * NcClusterMassRichness:MU :
   *
   * Distribution's  bias in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_RICHNESS_MU, "mu", "mu",
                              0.0,  6.0, 1.0e-1,
                              NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassRichness:MU_M:
   *
   * Distribution's slope with respect to the mass in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_RICHNESS_MU_M, "mu_m", "mum",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU_M,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassRichness:MU_Z:
   *
   * Distribution's slope with respect to the redshift in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_RICHNESS_MU_Z, "mu_z", "muz",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU_Z,
                              NCM_PARAM_TYPE_FIXED);

     /**
   * NcClusterMassRichness:MU_MZ:
   *
   * Distribution's slope with respect to the redshift  and mass in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_RICHNESS_MU_MZ, "mu_mz", "mumz",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU_MZ,
                              NCM_PARAM_TYPE_FIXED);

      /**
   * NcClusterMassRichness:SIGMA_0:
   *
   * Distribution's mass-richness scatter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_RICHNESS_SIGMA_0, "sigma_0", "sigma0",
                              0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_RICHNESS_DEFAULT_SIGMA_0,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->P               = &_nc_cluster_mass_richness_p;
  parent_class->intP            = &_nc_cluster_mass_richness_intp;
  parent_class->intP_bin        = &_nc_cluster_mass_richness_intp_bin;
  parent_class->resample        = &_nc_cluster_mass_richness_resample;
  parent_class->P_limits        = &_nc_cluster_mass_richness_p_limits;
  parent_class->P_bin_limits    = &_nc_cluster_mass_richness_p_limits_bin;
  parent_class->N_limits        = &_nc_cluster_mass_richness_n_limits;
  parent_class->_obs_len        = 1;
  parent_class->_obs_params_len = 1;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

typedef struct _NcClusterMassRichnessInt
{ 
  NcHICosmo *cosmo;
  NcClusterMassRichness *richness;
  const gdouble *richness_obs_params;
  const gdouble *richness_obs;
  const gdouble *richness_obs_lower;
  const gdouble *richness_obs_upper;
  gdouble lnM;
  gdouble z;
} NcClusterMassRichnessInt;


static void
_nc_cluster_mass_richness_mean (NcClusterMass *clusterm, const gdouble lnM, const gdouble z, gdouble *mean, gdouble *sigma)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  const gdouble DlnM                      = lnM - self->lnM0;
  const gdouble Dln1pz                    = log1p (z);

  mean[0]   = MU    + MU_M * DlnM + MU_Z * Dln1pz + MU_MZ * (DlnM) * (Dln1pz);
  sigma[0]  = SIGMA_0;
}

static gdouble
_nc_cluster_mass_richness_p_integrand (gdouble richness_true ,gpointer userdata)
{
  NcClusterMassRichnessInt *obs_data = (NcClusterMassRichnessInt *) userdata;
  NcClusterMass *clusterm = NC_CLUSTER_MASS ( obs_data->richness);

  gdouble mean;
  gdouble sigma;
  const gdouble sigma_obs       = obs_data->richness_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];
  const gdouble sqrt2_sigma_obs   = M_SQRT2 * sigma_obs;
  const gdouble x_obs             = (obs_data->richness_obs[0] - richness_true) / sqrt2_sigma_obs;
  
  const gdouble richness_obs_true = M_2_SQRTPI / (2.0 * M_SQRT2) * exp (-x_obs * x_obs) / (sigma_obs);
  

  _nc_cluster_mass_richness_mean (clusterm, obs_data->lnM, obs_data->z, &mean , &sigma);

  const gdouble sqrt2_sigma   = M_SQRT2 * sigma;
  const gdouble x             = (log(richness_true) - mean) / sqrt2_sigma;

  const gdouble mass_richness = M_2_SQRTPI / (2.0 * M_SQRT2) * exp (-x * x) / (sigma);

  return richness_obs_true * mass_richness;
}

static gdouble
_nc_cluster_mass_richness_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{  
  
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;
  NcClusterMassRichnessInt obs_data;

  gdouble clusterm_p, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  obs_data.cosmo               = cosmo;
  obs_data.richness            = richness;
  obs_data.richness_obs        = lnM_obs;
  obs_data.richness_obs_params = lnM_obs_params;
  obs_data.z                   = z;
  obs_data.lnM                 = lnM;

  const gdouble sigma_obs      = obs_data.richness_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];

  F.function = &_nc_cluster_mass_richness_p_integrand;
  F.params   = &obs_data;

  const gdouble richness_true_lb = obs_data.richness_obs[0]  - 7 * sigma_obs;
  const gdouble richness_true_ub = obs_data.richness_obs[0]  + 7 * sigma_obs;

  gsl_integration_qag (&F, richness_true_lb, richness_true_ub, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_RICHNESS_DEFAULT_INT_KEY, *w, &clusterm_p, &err);
  ncm_memory_pool_return (w);

  return clusterm_p;
}

static gdouble
_nc_cluster_mass_richness_intp_integrand (gdouble richness_true ,gpointer userdata)
{
  NcClusterMassRichnessInt *obs_data = (NcClusterMassRichnessInt *) userdata;
  NcClusterMassRichnessPrivate * const self = obs_data->richness->priv;
  NcClusterMass *clusterm = NC_CLUSTER_MASS ( obs_data->richness);
  gdouble richness_obs_true = 0;
  

  gdouble mean;
  gdouble sigma;
  const gdouble sigma_obs       = 25; //
  const gdouble sqrt2_sigma_obs   = M_SQRT2 * sigma_obs;
  const gdouble x_obs_min             = (self->richnessobs_min - richness_true) / sqrt2_sigma_obs;
  const gdouble x_obs_max             = (self->richnessobs_max - richness_true) / sqrt2_sigma_obs;
  
  if (x_obs_max > 4.0)
  {
    richness_obs_true =  -(erfc (x_obs_min) - erfc (x_obs_max)) / 2.0;}
  else
  {
    richness_obs_true =  (erf (x_obs_min) - erf (x_obs_max)) / 2.0;}
  

  _nc_cluster_mass_richness_mean (clusterm, obs_data->lnM, obs_data->z, &mean , &sigma);

  const gdouble sqrt2_sigma   = M_SQRT2 * sigma;
  const gdouble x             = (log(richness_true) - mean) / sqrt2_sigma;

  const gdouble mass_richness = M_2_SQRTPI / (2.0 * M_SQRT2) * exp (-x * x) / (sigma);

  return richness_obs_true * mass_richness;
}

static gdouble
_nc_cluster_mass_richness_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessInt obs_data;

  gdouble clusterm_intp, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  obs_data.cosmo          = cosmo;
  obs_data.richness       = richness;
  obs_data.z              = z;
  obs_data.lnM              = lnM;


  F.function = &_nc_cluster_mass_richness_p_integrand;
  F.params   = &obs_data;
  
  NcClusterMassRichnessPrivate * const self = richness->priv;

  const gdouble sigma_obs       = 25; //arrumar isso
  const gdouble richness_true_lb = self->richnessobs_min - 7 * sigma_obs;
  const gdouble richness_true_ub = self->richnessobs_max + 7 * sigma_obs;

  gsl_integration_qag (&F, richness_true_lb, richness_true_ub, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_RICHNESS_DEFAULT_INT_KEY, *w, &clusterm_intp, &err);
  ncm_memory_pool_return (w);

  return clusterm_intp;
}

static gdouble
_nc_cluster_mass_richness_intp_bin_integrand (gdouble richness_true ,gpointer userdata)
{
  NcClusterMassRichnessInt *obs_data = (NcClusterMassRichnessInt *) userdata;
  NcClusterMassRichnessPrivate * const self = obs_data->richness->priv;
  NcClusterMass *clusterm = NC_CLUSTER_MASS ( obs_data->richness);
  gdouble richness_obs_true = 0;
  

  gdouble mean;
  gdouble sigma;
  const gdouble sigma_obs       = obs_data->richness_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];
  const gdouble sqrt2_sigma_obs   = M_SQRT2 * sigma_obs;
  const gdouble x_obs_min             = (obs_data->richness_obs_lower[0] - richness_true) / sqrt2_sigma_obs;
  const gdouble x_obs_max             = (obs_data->richness_obs_upper[0] - richness_true) / sqrt2_sigma_obs;
  
  if (x_obs_max > 4.0)
  {
    richness_obs_true =  -(erfc (x_obs_min) - erfc (x_obs_max)) / 2.0;}
  else
  {
    richness_obs_true =  (erf (x_obs_min) - erf (x_obs_max)) / 2.0;}
  

  _nc_cluster_mass_richness_mean (clusterm, obs_data->lnM, obs_data->z, &mean , &sigma);

  const gdouble sqrt2_sigma   = M_SQRT2 * sigma;
  const gdouble x             = (log(richness_true) - mean) / sqrt2_sigma;

  const gdouble mass_richness = M_2_SQRTPI / (2.0 * M_SQRT2) * exp (-x * x) / (sigma);

  return richness_obs_true * mass_richness;
}

static gdouble
_nc_cluster_mass_richness_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessInt obs_data;

  gdouble clusterm_intp_bin, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  obs_data.cosmo          = cosmo;
  obs_data.richness       = richness;
  obs_data.richness_obs_lower = lnM_obs_lower;
  obs_data.richness_obs_upper = lnM_obs_upper;
  obs_data.richness_obs_params = lnM_obs_params;
  obs_data.z              = z;
  obs_data.lnM              = lnM;


  F.function = &_nc_cluster_mass_richness_p_integrand;
  F.params   = &obs_data;
  
  NcClusterMassRichnessPrivate * const self = richness->priv;

  const gdouble sigma_obs        = obs_data.richness_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];
  const gdouble richness_true_lb = obs_data.richness_obs_lower[0] - 7 * sigma_obs;
  const gdouble richness_true_ub = obs_data.richness_obs_upper[0] + 7 * sigma_obs;

  gsl_integration_qag (&F, richness_true_lb, richness_true_ub, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_RICHNESS_DEFAULT_INT_KEY, *w, &clusterm_intp_bin, &err);
  ncm_memory_pool_return (w);

  return clusterm_intp_bin;
}

static gboolean
_nc_cluster_mass_richness_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  gdouble mean;
  gdouble sigma;
  gdouble sigma_obs =  lnM_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];

  _nc_cluster_mass_richness_mean (clusterm, lnM, z, &mean , &sigma);


  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  ncm_rng_lock (rng);

  const gdouble richness_true = ncm_rng_gaussian_gen (rng, mean, sigma);
  lnM_obs[0] = ncm_rng_gaussian_gen (rng, richness_true, sigma_obs);
  ncm_rng_unlock (rng);

  return (lnM_obs[0] <= self->richnessobs_max) && (lnM_obs[0] >= self->richnessobs_min);
}

static void
_nc_cluster_mass_richness_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  gdouble x_lb , x_ub;
  const gdouble sigma_obs       = lnM_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];
  gdouble richness_true_lb = lnM_obs[0] - 7 * sigma_obs;
  gdouble richness_true_ub = lnM_obs[0] + 7 * sigma_obs;

  const gdouble m_lb =  (log(richness_true_lb) - MU - MU_Z) / (MU_M + MU_MZ);
  const gdouble m_ub =  (log(richness_true_ub) - MU - MU_Z) / (MU_M + MU_MZ);
  if (MU_MZ >= 0)
  {
    if (m_lb >=1)
    {
       x_lb = (log(richness_true_lb) - MU ) / (MU_M) + self->lnM0;
    }
    else
    {
       x_lb = -MU_Z / MU_MZ + self->lnM0;
    }
  }
  else
  {
     x_lb = M_LN10 * log10 (1e12);
  }
  
  const gdouble lnMl          = GSL_MAX(x_lb - 7.0 * SIGMA_0 , 0.0);

  if (MU_MZ >= 0)
  {
    if (m_ub >= 1)
    {
       x_ub = (log(richness_true_ub) - MU ) / (MU_M) + self->lnM0;
    }
    else
    {
       x_ub = -MU_Z / MU_MZ + self->lnM0;
    }
  }
  else
  {
     x_ub = M_LN10 * log10 (1e17);
  }
  const gdouble lnMu          = x_ub + 7.0 * SIGMA_0;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_richness_p_limits_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  gdouble x_lb , x_ub;
  gdouble sigma;
  const gdouble sigma_obs       = lnM_obs_params[NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS];
  gdouble richness_true_lb = lnM_obs_lower[0] - 7 * sigma_obs;
  gdouble richness_true_ub = lnM_obs_upper[0] + 7 * sigma_obs;

  const gdouble m_lb =  (log(richness_true_lb) - MU - MU_Z) / (MU_M + MU_MZ);
  const gdouble m_ub =  (log(richness_true_ub) - MU - MU_Z) / (MU_M + MU_MZ);
  if (MU_MZ >= 0)
  {
    if (m_lb >=1)
    {
       x_lb = (log(richness_true_lb) - MU ) / (MU_M) + self->lnM0;
    }
    else
    {
       x_lb = -MU_Z / MU_MZ + self->lnM0;
    }
  }
  else
  {
     x_lb = M_LN10 * log10 (1e12);
  }
  
  const gdouble lnMl          = GSL_MAX(x_lb - 7.0 * SIGMA_0 , 0.0);

  if (MU_MZ >= 0)
  {
    if (m_ub >= 1)
    {
       x_ub = (log(richness_true_ub) - MU ) / (MU_M) + self->lnM0;
    }
    else
    {
       x_ub = -MU_Z / MU_MZ + self->lnM0;
    }
  }
  else
  {
     x_ub = M_LN10 * log10 (1e17);
  }
  const gdouble lnMu          = x_ub + 7.0 * SIGMA_0;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;
}

static void
_nc_cluster_mass_richness_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  gdouble x_lb , x_ub;
  const gdouble sigma_obs       = 1; //arrumar
  gdouble richness_true_lb = self->richnessobs_min - 7 * sigma_obs;
  gdouble richness_true_ub = self->richnessobs_min + 7 * sigma_obs;

  const gdouble m_lb =  (log(richness_true_lb) - MU - MU_Z) / (MU_M + MU_MZ);
  const gdouble m_ub =  (log(richness_true_ub) - MU - MU_Z) / (MU_M + MU_MZ);
  if (MU_MZ >= 0)
  {
    if (m_lb >=1)
    {
       x_lb = (log(richness_true_lb) - MU ) / (MU_M) + self->lnM0;
    }
    else
    {
       x_lb = -MU_Z / MU_MZ + self->lnM0;
    }
  }
  else
  {
     x_lb = M_LN10 * log10 (1e12);
  }
  
  const gdouble lnMl          = GSL_MAX(x_lb - 7.0 * SIGMA_0 , 0.0);

  if (MU_MZ >= 0)
  {
    if (m_ub >= 1)
    {
       x_ub = (log(richness_true_ub) - MU ) / (MU_M) + self->lnM0;
    }
    else
    {
       x_ub = -MU_Z / MU_MZ + self->lnM0;
    }
  }
  else
  {
     x_ub = M_LN10 * log10 (1e17);
  }
  const gdouble lnMu          = x_ub + 7.0 * SIGMA_0;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static gdouble
_nc_cluster_mass_richness_volume (NcClusterMass *clusterm)
{
  NcClusterMassRichness *richness  = NC_CLUSTER_MASS_RICHNESS (clusterm);
  NcClusterMassRichnessPrivate * const self = richness->priv;

  return (self->richnessobs_min - self->richnessobs_max); //pensar nisso 
}