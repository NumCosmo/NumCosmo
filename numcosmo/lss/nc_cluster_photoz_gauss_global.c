/***************************************************************************
 *            nc_cluster_photoz_gauss_global.c
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_photoz_gauss_global
 * @title: NcClusterPhotozGaussGlobal
 * @short_description: Global gaussian photometric distribution for clusters.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_photoz_gauss_global.h"
#include "math/ncm_data.h"
#include "math/ncm_rng.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define VECTOR (model->params)
#define Z_BIAS (ncm_vector_get (VECTOR, NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_Z_BIAS))
#define SIGMA0 (ncm_vector_get (VECTOR, NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SIGMA0))

G_DEFINE_TYPE (NcClusterPhotozGaussGlobal, nc_cluster_photoz_gauss_global, NC_TYPE_CLUSTER_REDSHIFT);

enum
{
  PROP_0,
  PROP_PZ_MIN,
  PROP_PZ_MAX,
  PROP_SIZE
};

static void
nc_cluster_photoz_gauss_global_init (NcClusterPhotozGaussGlobal *pzg_global)
{
}

static void
_nc_cluster_photoz_gauss_global_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_cluster_photoz_gauss_global_parent_class)->constructed (object);
  {
    NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object);
    
    g_assert_cmpfloat (pzg_global->pz_min, <, pzg_global->pz_max);
  }
}

static void
_nc_cluster_photoz_gauss_global_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object);
  
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object));
  
  switch (prop_id)
  {
    case PROP_PZ_MIN:
      pzg_global->pz_min = g_value_get_double (value);
      break;
    case PROP_PZ_MAX:
      pzg_global->pz_max = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_photoz_gauss_global_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object);
  
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object));
  
  switch (prop_id)
  {
    case PROP_PZ_MIN:
      g_value_set_double (value, pzg_global->pz_min);
      break;
    case PROP_PZ_MAX:
      g_value_set_double (value, pzg_global->pz_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_photoz_gauss_global_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_photoz_gauss_global_parent_class)->finalize (object);
}

static gdouble _nc_cluster_photoz_gauss_global_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params);
static gdouble _nc_cluster_photoz_gauss_global_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z);
static gdouble _nc_cluster_photoz_gauss_global_intp_bin(NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params);
static gboolean _nc_cluster_photoz_gauss_global_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng);
static void _nc_cluster_photoz_gauss_global_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
static void _nc_cluster_photoz_gauss_global_p_bin_limits(NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
static void _nc_cluster_photoz_gauss_global_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper);

static void
nc_cluster_photoz_gauss_global_class_init (NcClusterPhotozGaussGlobalClass *klass)
{
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcClusterRedshiftClass *parent_class = NC_CLUSTER_REDSHIFT_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);
  
  object_class->constructed = &_nc_cluster_photoz_gauss_global_constructed;
  model_class->set_property = &_nc_cluster_photoz_gauss_global_set_property;
  model_class->get_property = &_nc_cluster_photoz_gauss_global_get_property;
  object_class->finalize    = &_nc_cluster_photoz_gauss_global_finalize;
  
  ncm_model_class_set_name_nick (model_class, "Global Gaussian distribution", "GaussianGlobal");
  ncm_model_class_add_params (model_class, 2, 0, PROP_SIZE);
  
  /**
   * NcClusterPhotozGaussGlobal:pz_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ_MIN,
                                   g_param_spec_double ("pz-min",
                                                        NULL,
                                                        "Minimum photoz",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcClusterPhotozGaussGlobal:pz_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ_MAX,
                                   g_param_spec_double ("pz-max",
                                                        NULL,
                                                        "Maximum photoz",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcClusterPhotozGaussGlobal:z_bias:
   *
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_Z_BIAS, "z-bias", "z-bias",
                              -G_MAXDOUBLE, G_MAXDOUBLE, 1.0e-2,
                              NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_BIAS,
                              NCM_PARAM_TYPE_FIXED);
  
  /**
   * NcClusterPhotozGaussGlobal:sigma0:
   *
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SIGMA0, "sigma0", "sigma0",
                              0.0, G_MAXDOUBLE, 1.0e-2,
                              NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_SIGMA0,
                              NCM_PARAM_TYPE_FIXED);
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
  
  parent_class->P              = &_nc_cluster_photoz_gauss_global_p;
  parent_class->intP           = &_nc_cluster_photoz_gauss_global_intp;
  parent_class->intP_bin       = &_nc_cluster_photoz_gauss_global_intp_bin;
  parent_class->resample       = &_nc_cluster_photoz_gauss_global_resample;
  parent_class->P_limits       = &_nc_cluster_photoz_gauss_global_p_limits;
  parent_class->P_bin_limits   = &_nc_cluster_photoz_gauss_global_p_bin_limits;
  parent_class->N_limits       = &_nc_cluster_photoz_gauss_global_n_limits;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 0;
  
  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_REDSHIFT_IMPL_ALL);
}

static gdouble
_nc_cluster_photoz_gauss_global_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params)
{
  NcmModel *model           = NCM_MODEL (clusterz);
  const gdouble z_eff       = z + Z_BIAS;
  const gdouble sqrt2_sigma = M_SQRT2 * SIGMA0 * (1.0 + z);
  const gdouble y1          = (z_obs[0] - z_eff) / sqrt2_sigma;
  
  NCM_UNUSED (lnM);
  NCM_UNUSED (z_obs_params);
  
  return M_2_SQRTPI / M_SQRT2 * exp (-y1 * y1) / (SIGMA0 * (1.0 + z) * (1.0 + erf (z_eff / sqrt2_sigma)));
}

static gdouble
_nc_cluster_photoz_gauss_global_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z)
{
  NcmModel *model                        = NCM_MODEL (clusterz);
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble z_eff                    = z + Z_BIAS;
  const gdouble sqrt2_sigma              = M_SQRT2 * SIGMA0 * (1.0 + z);
  const gdouble x_min                    = (z_eff - pzg_global->pz_min) / sqrt2_sigma;
  const gdouble x_max                    = (z_eff - pzg_global->pz_max) / sqrt2_sigma;
  
  NCM_UNUSED (lnM);
  
  if (x_max > 4.0)
    return -(erfc (x_min) - erfc (x_max)) /
           (1.0 + erf (z_eff / sqrt2_sigma));
  else
    return (erf (x_min) - erf (x_max)) /
           (1.0 + erf (z_eff / sqrt2_sigma));
}

static gdouble _nc_cluster_photoz_gauss_global_intp_bin(NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params)
{
NcmModel *model                        = NCM_MODEL (clusterz);
/*NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);*/
const gdouble z_eff                    = z + Z_BIAS;
const gdouble sqrt2_sigma              = M_SQRT2 * SIGMA0 *(1.0 + z);
const gdouble x_min                    = (z_eff - z_obs_lower[0]) / sqrt2_sigma;
const gdouble x_max                    = (z_eff - z_obs_upper[0]) / sqrt2_sigma;
NCM_UNUSED (lnM);

if (x_max > 4.0)
    return -(erfc (x_min) - erfc (x_max)) /
           (1.0 + erf (z_eff / sqrt2_sigma));
  else
    return (erf (x_min) - erf (x_max)) /
           (1.0 + erf (z_eff / sqrt2_sigma)); 
}




static gboolean
_nc_cluster_photoz_gauss_global_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng)
{
  NcmModel *model                        = NCM_MODEL (clusterz);
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble sigma_z                  = SIGMA0 * (1.0 + z);
  
  ncm_rng_lock (rng);
  
  do {
    z_obs[0] = z + Z_BIAS + gsl_ran_gaussian (rng->r, sigma_z);
  } while (z_obs[0] < 0.0);
  
  ncm_rng_unlock (rng);
  
  return (z_obs[0] <= pzg_global->pz_max) && (z_obs[0] >= pzg_global->pz_min);
}

static void
_nc_cluster_photoz_gauss_global_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  NcmModel *model    = NCM_MODEL (clusterz);
  const gdouble mean = z_obs[0] - Z_BIAS;
  const gdouble zl   = GSL_MAX (mean - 10.0 * SIGMA0 * (1.0 + z_obs[0]), 0.0);
  const gdouble zu   = mean + 10.0 * SIGMA0 * (1.0 + z_obs[0]);
  
  NCM_UNUSED (z_obs_params);
  
  *z_lower = zl;
  *z_upper = zu;
  
  return;
}

static void
_nc_cluster_photoz_gauss_global_p_bin_limits(NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
NcmModel *model                        = NCM_MODEL (clusterz);
/*NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);*/
const gdouble zl                       = GSL_MAX (z_obs_lower[0] - Z_BIAS - 10.0 * SIGMA0 * (1.0 + z_obs_lower[0]),0.0);
const gdouble zu                       = GSL_MAX (z_obs_upper[0] - Z_BIAS + 10.0 * SIGMA0 * (1.0 + z_obs_upper[0]),0.0);
*z_lower = zl;
*z_upper = zu;

return;
}


static void
_nc_cluster_photoz_gauss_global_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper)
{
  NcmModel *model                        = NCM_MODEL (clusterz);
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble zl                       = GSL_MAX (pzg_global->pz_min - Z_BIAS - 10.0 * SIGMA0 * (1.0 + pzg_global->pz_min), 0.0);
  const gdouble zu                       =          pzg_global->pz_max - Z_BIAS + 10.0 * SIGMA0 * (1.0 + pzg_global->pz_max);
  
  *z_lower = zl;
  *z_upper = zu;
  
  return;
}

/**
 * nc_cluster_photoz_gauss_global_new:
 * @pz_min: FIXME
 * @pz_max: FIXME
 * @z_bias: FIXME
 * @sigma0: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcClusterRedshift.
 */
NcClusterRedshift *
nc_cluster_photoz_gauss_global_new (gdouble pz_min, const gdouble pz_max, const gdouble z_bias, const gdouble sigma0)
{
  return g_object_new (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL,
                       "pz-min", pz_min,
                       "pz-max", pz_max,
                       "z-bias", z_bias,
                       "sigma0", sigma0,
                       NULL);
}

/**
 * nc_cluster_photoz_gauss_global_set_z_bias:
 * @pzg_global: a #NcClusterPhotozGaussGlobal.
 * @z_bias: value of #NcClusterPhotozGaussGlobal:z-bias.
 *
 * Sets the value @z_bias to the #NcClusterPhotozGaussGlobal:z-bias property.
 *
 */
void
nc_cluster_photoz_gauss_global_set_z_bias (NcClusterPhotozGaussGlobal *pzg_global, const gdouble z_bias)
{
  NcmModel *model = NCM_MODEL (pzg_global);
  
  ncm_model_param_set (model, NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_Z_BIAS, z_bias);
}

/**
 * nc_cluster_photoz_gauss_global_get_z_bias:
 * @pzg_global: a #NcClusterPhotozGaussGlobal.
 *
 * Returns: the value of #NcClusterPhotozGaussGlobal:z-bias property.
 */
gdouble
nc_cluster_photoz_gauss_global_get_z_bias (const NcClusterPhotozGaussGlobal *pzg_global)
{
  NcmModel *model = NCM_MODEL (pzg_global);
  
  return Z_BIAS;
}

/**
 * nc_cluster_photoz_gauss_global_set_sigma0:
 * @pzg_global: a #NcClusterPhotozGaussGlobal.
 * @sigma0: value of #NcClusterPhotozGaussGlobal:sigma0.
 *
 * Sets the value @sigma0 to the #NcClusterPhotozGaussGlobal:sigma0 property.
 *
 */
void
nc_cluster_photoz_gauss_global_set_sigma0 (NcClusterPhotozGaussGlobal *pzg_global, const gdouble sigma0)
{
  NcmModel *model = NCM_MODEL (pzg_global);
  
  ncm_model_param_set (model, NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SIGMA0, sigma0);
}

/**
 * nc_cluster_photoz_gauss_global_get_sigma0:
 * @pzg_global: a #NcClusterPhotozGaussGlobal.
 *
 * Returns: the value of #NcClusterPhotozGaussGlobal:sigma0 property.
 */
gdouble
nc_cluster_photoz_gauss_global_get_sigma0 (const NcClusterPhotozGaussGlobal *pzg_global)
{
  NcmModel *model = NCM_MODEL (pzg_global);
  
  return SIGMA0;
}

