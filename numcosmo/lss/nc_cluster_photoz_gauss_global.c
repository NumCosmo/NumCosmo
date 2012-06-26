/***************************************************************************
 *            nc_cluster_photoz_gauss_global.c
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
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
 * @title: Global Gaussian Photoz Cluster
 * @short_description: Gaussian photometric redshift
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>

G_DEFINE_TYPE (NcClusterPhotozGaussGlobal, nc_cluster_photoz_gauss_global, NC_TYPE_CLUSTER_REDSHIFT);

enum
{
  PROP_0,
  PROP_PZ_MIN,
  PROP_PZ_MAX,
  PROP_Z_BIAS,
  PROP_SIGMA0,
};

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
nc_cluster_photoz_gauss_global_new (gdouble pz_min, gdouble pz_max, gdouble z_bias, gdouble sigma0)
{
  return g_object_new (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL,
                       "pz-min", pz_min,
                       "pz-max", pz_max,
                       "z-bias", z_bias,
                       "sigma0", sigma0,
                       NULL);
}

static gdouble
_nc_cluster_photoz_gauss_global_p (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble z_eff = z + pzg_global->z_bias;
  const gdouble sqrt2_sigma = M_SQRT2 * pzg_global->sigma0 * (1.0 + z);
  const gdouble y1 = (z_obs[0] - z_eff) / sqrt2_sigma;

  return M_2_SQRTPI / M_SQRT2 * exp (- y1 * y1) / (pzg_global->sigma0 * (1.0 + z) * (1.0 + erf (z_eff / sqrt2_sigma)));
}

static gdouble
_nc_cluster_photoz_gauss_global_intp (NcClusterRedshift *clusterz, gdouble lnM, gdouble z)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble z_eff = z + pzg_global->z_bias;
  const gdouble sqrt2_sigma = M_SQRT2 * pzg_global->sigma0 * (1.0 + z);
  const gdouble x_min = (z_eff - pzg_global->pz_min) / sqrt2_sigma;
  const gdouble x_max = (z_eff - pzg_global->pz_max) / sqrt2_sigma;

  if (x_max > 4.0)
  {
	return -(erfc (x_min) - erfc (x_max)) / 
	  (1.0 + erf (z_eff / sqrt2_sigma));
  }
  else
  {
	return (erf (x_min) - erf (x_max)) / 
	  (1.0 + erf (z_eff / sqrt2_sigma));
  }
}

static gboolean
_nc_cluster_photoz_gauss_global_resample (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble sigma_z = pzg_global->sigma0 * (1.0 + z);
  gsl_rng *rng = ncm_get_rng ();

  do {
	z_obs[0] = z + pzg_global->z_bias + gsl_ran_gaussian (rng, sigma_z);
  } while (z_obs[0] < 0.0);

  return (z_obs[0] <= pzg_global->pz_max) && (z_obs[0] >= pzg_global->pz_min);
}

static void
_nc_cluster_photoz_gauss_global_p_limits (NcClusterRedshift *clusterz, gdouble *z_obs, gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble mean = z_obs[0] - pzg_global->z_bias;
  const gdouble zl = GSL_MAX (mean - 10.0 * pzg_global->sigma0 * (1.0 + z_obs[0]), 0.0);
  const gdouble zu = mean + 10.0 * pzg_global->sigma0 * (1.0 + z_obs[0]);

  *z_lower = zl;
  *z_upper = zu;

  return;
}

static void
_nc_cluster_photoz_gauss_global_n_limits (NcClusterRedshift *clusterz, gdouble *z_lower, gdouble *z_upper)
{
  NcClusterPhotozGaussGlobal *pzg_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (clusterz);
  const gdouble zl = GSL_MAX (pzg_global->pz_min - pzg_global->z_bias - 10.0 * pzg_global->sigma0 * (1.0 + pzg_global->pz_min), 0.0);
  const gdouble zu = pzg_global->pz_max - pzg_global->z_bias + 10.0 * pzg_global->sigma0 * (1.0 + pzg_global->pz_max);

  *z_lower = zl;
  *z_upper = zu;

  return;
}

guint _nc_cluster_photoz_gauss_global_obs_len (NcClusterRedshift *clusterz) { return 1; }
guint _nc_cluster_photoz_gauss_global_obs_params_len (NcClusterRedshift *clusterz) { return 0; }

/**
 * nc_cluster_photoz_gauss_global_set_z_bias:
 * @pzg_global: a #NcClusterPhotozGaussGlobal.
 * @z_bias: value of #NcClusterPhotozGaussGlobal:z-bias.
 *
 * Sets the value @z_bias to the #NcClusterPhotozGaussGlobal:z-bias property.
 *
 */
void
nc_cluster_photoz_gauss_global_set_z_bias (NcClusterPhotozGaussGlobal *pzg_global, gdouble z_bias)
{
  pzg_global->z_bias = z_bias;
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
  return pzg_global->z_bias;
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
nc_cluster_photoz_gauss_global_set_sigma0 (NcClusterPhotozGaussGlobal *pzg_global, gdouble sigma0)
{
  pzg_global->sigma0 = sigma0;
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
  return pzg_global->sigma0;
}

static void
nc_cluster_photoz_gauss_global_init (NcClusterPhotozGaussGlobal *pzg_global)
{
  pzg_global->z_bias = 0.0;
  pzg_global->sigma0 = 0.0;
}

static void
_nc_cluster_photoz_gauss_global_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_photoz_gauss_global_parent_class)->finalize (object);
}

static void
_nc_cluster_photoz_gauss_global_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
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
    case PROP_Z_BIAS:
      pzg_global->z_bias = g_value_get_double (value);
      break;
	case PROP_SIGMA0:
      pzg_global->sigma0 = g_value_get_double (value);
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
    case PROP_Z_BIAS:
      g_value_set_double (value, pzg_global->z_bias);
      break;
	case PROP_SIGMA0:
      g_value_set_double (value, pzg_global->sigma0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_photoz_gauss_global_class_init (NcClusterPhotozGaussGlobalClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterRedshiftClass* parent_class = NC_CLUSTER_REDSHIFT_CLASS (klass);

  parent_class->P              = &_nc_cluster_photoz_gauss_global_p;
  parent_class->intP           = &_nc_cluster_photoz_gauss_global_intp;
  parent_class->resample       = &_nc_cluster_photoz_gauss_global_resample;
  parent_class->P_limits       = &_nc_cluster_photoz_gauss_global_p_limits;
  parent_class->N_limits       = &_nc_cluster_photoz_gauss_global_n_limits;
  parent_class->obs_len        = &_nc_cluster_photoz_gauss_global_obs_len;
  parent_class->obs_params_len = &_nc_cluster_photoz_gauss_global_obs_params_len;

  parent_class->impl = NC_CLUSTER_REDSHIFT_IMPL_ALL;

  object_class->finalize = _nc_cluster_photoz_gauss_global_finalize;
  object_class->set_property = _nc_cluster_photoz_gauss_global_set_property;
  object_class->get_property = _nc_cluster_photoz_gauss_global_get_property;

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterPhotozGaussGlobal:z_bias:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_BIAS,
                                   g_param_spec_double ("z-bias",
                                                        NULL,
                                                        "z-bias",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterPhotozGaussGlobal:sigma0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA0,
                                   g_param_spec_double ("sigma0",
                                                        NULL,
                                                        "sigma0",
                                                        0.0, G_MAXDOUBLE, 0.03,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

