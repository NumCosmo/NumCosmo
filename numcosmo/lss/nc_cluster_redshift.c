/***************************************************************************
 *            nc_cluster_redshift.c
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
 * SECTION:nc_cluster_redshift
 * @title: Abstract Cluster Redshift Object
 * @short_description: Observed redshift distribution
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_redshift.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_ABSTRACT_TYPE (NcClusterRedshift, nc_cluster_redshift, G_TYPE_OBJECT);

/**
 * nc_cluster_redshift_new_from_name:
 * @redshift_name: string which specifies the type of the redshift distribution.
 *
 * This function returns a new #NcClusterRedshift whose type is defined by @redshift_name.
 *
 * Returns: A new #NcClusterRedshift.
 */
NcClusterRedshift *
nc_cluster_redshift_new_from_name (gchar *redshift_name)
{
  GObject *obj = ncm_serialize_global_from_string (redshift_name);
  GType redshift_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (redshift_type, NC_TYPE_CLUSTER_REDSHIFT))
	g_error ("nc_cluster_redshift_new_from_name: NcClusterRedshift %s do not descend from %s\n", redshift_name, g_type_name (NC_TYPE_CLUSTER_REDSHIFT));
  return NC_CLUSTER_REDSHIFT (obj);
}

/**
 * nc_cluster_redshift_ref:
 * @clusterz: FIXME.
 *
 * FIXME
 *
 * Returns: (transfer full): @clusterz.
 */
NcClusterRedshift *
nc_cluster_redshift_ref (NcClusterRedshift *clusterz)
{
  return g_object_ref (clusterz);
}

/**
 * nc_cluster_redshift_free:
 * @clusterz: a #NcClusterRedshift.
 *
 * Atomically decrements the reference count of @clusterz by one. If the reference count drops to 0,
 * all memory allocated by @clusterz is released.
 *
 */
void
nc_cluster_redshift_free (NcClusterRedshift *clusterz)
{
  g_object_unref (clusterz);
}

/**
 * nc_cluster_redshift_clear:
 * @clusterz: a #NcClusterRedshift.
 *
 * Atomically decrements the reference count of @clusterz by one. If the reference count drops to 0,
 * all memory allocated by @clusterz is released. Set pointer to NULL.
 *
 */
void
nc_cluster_redshift_clear (NcClusterRedshift **clusterz)
{
  g_clear_object (clusterz);
}

/**
 * nc_cluster_redshift_impl:
 * @clusterz: FIXME.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcClusterRedshiftImpl 
nc_cluster_redshift_impl (NcClusterRedshift *clusterz)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->impl;
}

/**
 * nc_cluster_redshift_obs_len:
 * @clusterz: FIXME.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
nc_cluster_redshift_obs_len (NcClusterRedshift *clusterz)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->obs_len (clusterz);
}

/**
 * nc_cluster_redshift_obs_params_len:
 * @clusterz: FIXME.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
nc_cluster_redshift_obs_params_len (NcClusterRedshift *clusterz)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->obs_params_len (clusterz);
}

/**
 * nc_cluster_redshift_p:
 * @clusterz: a #NcClusterRedshift.
 * @z: true redshift.
 * @lnM: true mass.
 * @z_obs: observed redshift.
 * @z_obs_params: observed redshift params.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_redshift_p (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->P (clusterz, lnM, z, z_obs, z_obs_params);
}

/**
 * nc_cluster_redshift_intp:
 * @clusterz: a #NcClusterRedshift.
 * @z: true redshift.
 * @lnM: true mass.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_cluster_redshift_intp (NcClusterRedshift *clusterz, gdouble lnM, gdouble z)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->intP (clusterz, lnM, z);
}

/**
 * nc_cluster_redshift_resample:
 * @clusterz: a #NcClusterRedshift.
 * @z: true redshift.
 * @lnM: true mass.
 * @z_obs: (out): observed redshift.
 * @z_obs_params: (out): observed redshift params.
   @rng: a #NcmRNG.
 *
 * FIXME
 * The function which will call this one is responsible to allocate enough memory for @z_lower and @z_upper.
 *
 * Returns: FIXME
 */
gboolean
nc_cluster_redshift_resample (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params, NcmRNG *rng)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->resample (clusterz, lnM, z, z_obs, z_obs_params, rng);
}

/**
 * nc_cluster_redshift_p_limits:
 * @clusterz: a #NcClusterRedshift.
 * @z_obs: observed redshift.
 * @z_obs_params: observed redshift params.
 * @z_lower: (out): pointer to the lower limit of the true redshift integration.
 * @z_upper: (out): pointer to the upper limit of the true redshift integration.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @z_lower and @z_upper.
*/
void
nc_cluster_redshift_p_limits (NcClusterRedshift *clusterz, gdouble *z_obs, gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->P_limits (clusterz, z_obs, z_obs_params, z_lower, z_upper);
}

/**
 * nc_cluster_redshift_n_limits:
 * @clusterz: a #NcClusterRedshift.
 * @z_lower: (out): pointer to the lower limit of the true redshift.
 * @z_upper: (out): pointer to the upper limit of the true redshift.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @z_lower and @z_upper.
 */
void
nc_cluster_redshift_n_limits (NcClusterRedshift *clusterz, gdouble *z_lower, gdouble *z_upper)
{
  NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->N_limits (clusterz, z_lower, z_upper);
}

static void
nc_cluster_redshift_init (NcClusterRedshift *nc_cluster_redshift)
{
  NCM_UNUSED (nc_cluster_redshift);
}

static void
_nc_cluster_redshift_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_redshift_parent_class)->finalize (object);
}

static void
nc_cluster_redshift_class_init (NcClusterRedshiftClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->finalize = _nc_cluster_redshift_finalize;
}

