/***************************************************************************
 *            nc_cluster_photoz_gauss.c
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
 * SECTION:nc_cluster_photoz_gauss
 * @title: NcClusterPhotozGauss
 * @short_description: Individual gaussian photometric distribution for clusters.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_photoz_gauss.h"
#include "math/ncm_data.h"
#include "math/ncm_rng.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcClusterPhotozGauss, nc_cluster_photoz_gauss, NC_TYPE_CLUSTER_REDSHIFT);

enum
{
  PROP_0,
  PROP_PZ_MIN,
  PROP_PZ_MAX,
  PROP_SIZE
};

static void
nc_cluster_photoz_gauss_init (NcClusterPhotozGauss *pzg)
{
  pzg->pz_max = 0.0;
  pzg->pz_min = 0.0;
}

static void
_nc_cluster_photoz_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_photoz_gauss_parent_class)->finalize (object);
}

static void
_nc_cluster_photoz_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (object);
  
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS (object));
  
  switch (prop_id)
  {
    case PROP_PZ_MIN:
      pzg->pz_min = g_value_get_double (value);
      break;
    case PROP_PZ_MAX:
      pzg->pz_max = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_photoz_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (object);
  
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS (object));
  
  switch (prop_id)
  {
    case PROP_PZ_MIN:
      g_value_set_double (value, pzg->pz_min);
      break;
    case PROP_PZ_MAX:
      g_value_set_double (value, pzg->pz_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble _nc_cluster_photoz_gauss_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params);
static gdouble _nc_cluster_photoz_gauss_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z);
static gboolean _nc_cluster_photoz_gauss_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng);
static void _nc_cluster_photoz_gauss_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
static void _nc_cluster_photoz_gauss_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper);

static void
nc_cluster_photoz_gauss_class_init (NcClusterPhotozGaussClass *klass)
{
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcClusterRedshiftClass *parent_class = NC_CLUSTER_REDSHIFT_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);
  
  object_class->finalize = _nc_cluster_photoz_gauss_finalize;
  
  model_class->set_property = _nc_cluster_photoz_gauss_set_property;
  model_class->get_property = _nc_cluster_photoz_gauss_get_property;
  
  ncm_model_class_set_name_nick (model_class, "Gaussian distribution", "Gaussian");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);
  
  /**
   * NcClusterPhotozGauss:pz_min:
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
   * NcClusterPhotozGauss:pz_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ_MAX,
                                   g_param_spec_double ("pz-max",
                                                        NULL,
                                                        "Maximum photoz",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
  
  parent_class->P              = &_nc_cluster_photoz_gauss_p;
  parent_class->intP           = &_nc_cluster_photoz_gauss_intp;
  parent_class->resample       = &_nc_cluster_photoz_gauss_resample;
  parent_class->P_limits       = &_nc_cluster_photoz_gauss_p_limits;
  parent_class->N_limits       = &_nc_cluster_photoz_gauss_n_limits;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 2;
  
  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_REDSHIFT_IMPL_ALL);
}

static gdouble
_nc_cluster_photoz_gauss_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params)
{
  const gdouble z_bias      = z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS];
  const gdouble z_eff       = z_bias + z;
  const gdouble sigma       = z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA];
  const gdouble sqrt2_sigma = M_SQRT2 * sigma;
  const gdouble y1          = (z_obs[0] - z_eff) / sqrt2_sigma;
  
  return M_2_SQRTPI / M_SQRT2 * exp (-y1 * y1) / (sigma * (1.0 + erf (z_eff / sqrt2_sigma)));
}

static gdouble
_nc_cluster_photoz_gauss_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (clusterz);
  const gdouble sigma       = 0.03 * (1.0 + z);
  const gdouble sqrt2_sigma = M_SQRT2 * sigma;
  const gdouble x_min       = (z - pzg->pz_min) / sqrt2_sigma;
  const gdouble x_max       = (z - pzg->pz_max) / sqrt2_sigma;
  
  if (x_max > 4.0)
    return -(erfc (x_min) - erfc (x_max)) /
           (1.0 + erf (z / sqrt2_sigma));
  else
    return (erf (x_min) - erf (x_max)) /
           (1.0 + erf (z / sqrt2_sigma));
}

static gboolean
_nc_cluster_photoz_gauss_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (clusterz);
  gdouble sigma_z;
  
  sigma_z = z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA];
  
  ncm_rng_lock (rng);
  
  do {
    z_obs[0] = z + z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS] + gsl_ran_gaussian (rng->r, sigma_z);
  } while (z_obs[0] < 0);
  
  ncm_rng_unlock (rng);
  
  return (z_obs[0] <= pzg->pz_max) && (z_obs[0] >= pzg->pz_min);
}

static void
_nc_cluster_photoz_gauss_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  const gdouble mean = z_obs[0] - z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS];
  const gdouble zl   = GSL_MAX (mean - 10.0 * z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA], 0.0);
  const gdouble zu   = mean + 10.0 * z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA];
  
  *z_lower = zl;
  *z_upper = zu;
  
  return;
}

static void
_nc_cluster_photoz_gauss_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (clusterz);
  const gdouble zl          = GSL_MAX (pzg->pz_min - 10.0 * 0.03 * (1.0 + pzg->pz_min), 0.0);
  const gdouble zu          = pzg->pz_max + 10.0 * 0.03 * (1.0 + pzg->pz_max);
  
  *z_lower = zl;
  *z_upper = zu;
  
  return;
}

/**
 * nc_cluster_photoz_gauss_new:
 *
 * FIXME
 *
 * Returns: A new #NcClusterRedshift.
 */
NcClusterRedshift *
nc_cluster_photoz_gauss_new ()
{
  return g_object_new (NC_TYPE_CLUSTER_PHOTOZ_GAUSS, NULL);
}

