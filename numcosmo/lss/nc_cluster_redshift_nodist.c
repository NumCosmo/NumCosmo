/***************************************************************************
 *            nc_cluster_redshift_nodist.c
 *
 *  Fri June 22 13:44:51 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_cluster_redshift_nodist
 * @title: Cluster Abundance Redshift No Distribution
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_redshift_nodist.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#include <gsl/gsl_math.h>

G_DEFINE_TYPE (NcClusterRedshiftNodist, nc_cluster_redshift_nodist, NC_TYPE_CLUSTER_REDSHIFT);

enum
{
  PROP_0,
  PROP_Z_MIN,
  PROP_Z_MAX,
  PROP_SIZE
};

static gdouble
_nc_cluster_redshift_nodist_p (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params)
{
  g_error ("This object don't implement p.");
  NCM_UNUSED (clusterz);
  NCM_UNUSED (lnM);
  NCM_UNUSED (z);
  NCM_UNUSED (z_obs);
  NCM_UNUSED (z_obs_params);
  return GSL_NAN;
}

static gdouble
_nc_cluster_redshift_nodist_intp (NcClusterRedshift *clusterz, gdouble lnM, gdouble z)
{
  g_error ("This object don't implement n_z_lnm.");
  NCM_UNUSED (clusterz);
  NCM_UNUSED (lnM);
  NCM_UNUSED (z);
  return GSL_NAN;
}

static gboolean
_nc_cluster_redshift_nodist_resample (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params, NcmRNG *rng)
{
  NcClusterRedshiftNodist *zn = NC_CLUSTER_REDSHIFT_NODIST (clusterz);

  NCM_UNUSED (lnM);
  NCM_UNUSED (z);
  NCM_UNUSED (z_obs_params);

  z_obs[0] = z;
  return (z_obs[0] <= zn->z_max) && (z_obs[0] >= zn->z_min);
}

static void
_nc_cluster_redshift_nodist_p_limits (NcClusterRedshift *clusterz, gdouble *z_obs, gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  g_error ("This object don't implement p_limits.");

  NCM_UNUSED (clusterz);
  NCM_UNUSED (z_obs);
  NCM_UNUSED (z_obs_params);
  NCM_UNUSED (z_lower);
  NCM_UNUSED (z_upper);
  
  return;
}

static void
_nc_cluster_redshift_nodist_n_limits (NcClusterRedshift *clusterz, gdouble *z_lower, gdouble *z_upper)
{
  NcClusterRedshiftNodist *zn = NC_CLUSTER_REDSHIFT_NODIST (clusterz);

  *z_lower = zn->z_min;
  *z_upper = zn->z_max;

  return;
}

guint _nc_cluster_redshift_nodist_obs_len (NcClusterRedshift *clusterz) { NCM_UNUSED (clusterz); return 1; }
guint _nc_cluster_redshift_nodist_obs_params_len (NcClusterRedshift *clusterz) { NCM_UNUSED (clusterz); return 0; }

static void
_nc_cluster_redshift_nodist_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterRedshiftNodist *zn = NC_CLUSTER_REDSHIFT_NODIST (object);
  g_return_if_fail (NC_IS_CLUSTER_REDSHIFT_NODIST (object));

  switch (prop_id)
  {
    case PROP_Z_MIN:
      zn->z_min = g_value_get_double (value);
      break;
	case PROP_Z_MAX:
      zn->z_max = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_redshift_nodist_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterRedshiftNodist *zn = NC_CLUSTER_REDSHIFT_NODIST (object);
  g_return_if_fail (NC_IS_CLUSTER_REDSHIFT_NODIST (object));

  switch (prop_id)
  {
    case PROP_Z_MIN:
      g_value_set_double (value, zn->z_min);
      break;
	case PROP_Z_MAX:
      g_value_set_double (value, zn->z_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_redshift_nodist_init (NcClusterRedshiftNodist *zn)
{
  zn->z_min = 0.0;
  zn->z_max = 0.0;
}

static void
nc_cluster_redshift_nodist_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_redshift_nodist_parent_class)->finalize (object);
}

static void
nc_cluster_redshift_nodist_class_init (NcClusterRedshiftNodistClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterRedshiftClass* parent_class = NC_CLUSTER_REDSHIFT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  parent_class->P              = &_nc_cluster_redshift_nodist_p;
  parent_class->intP           = &_nc_cluster_redshift_nodist_intp;
  parent_class->resample       = &_nc_cluster_redshift_nodist_resample;
  parent_class->P_limits       = &_nc_cluster_redshift_nodist_p_limits;
  parent_class->N_limits       = &_nc_cluster_redshift_nodist_n_limits;
  parent_class->obs_len        = &_nc_cluster_redshift_nodist_obs_len;
  parent_class->obs_params_len = &_nc_cluster_redshift_nodist_obs_params_len;

  parent_class->impl           = NC_CLUSTER_REDSHIFT_N_LIMTS | NC_CLUSTER_REDSHIFT_RESAMPLE;

  object_class->finalize    =  &nc_cluster_redshift_nodist_finalize;
  
  model_class->set_property = &_nc_cluster_redshift_nodist_set_property;
  model_class->get_property = &_nc_cluster_redshift_nodist_get_property;

  ncm_model_class_set_name_nick (model_class, "No redshift distribution", "No_distribution");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);
  
  /**
   * NcClusterRedshiftNodist:z_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_MIN,
                                   g_param_spec_double ("z-min",
                                                        NULL,
                                                        "Minimum z",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterRedshiftNodist:z_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_MAX,
                                   g_param_spec_double ("z-max",
                                                        NULL,
                                                        "Maximum z",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}
