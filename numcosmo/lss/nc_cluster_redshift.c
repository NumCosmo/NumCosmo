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
 * @title: NcClusterRedshift
 * @short_description: Abstract class for cluster redshift distributions.
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

struct _NcClusterRedshiftPrivate
{
  guint place_holder;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcClusterRedshift, nc_cluster_redshift, NCM_TYPE_MODEL);

static void
nc_cluster_redshift_init (NcClusterRedshift *clusterz)
{
  NcClusterRedshiftPrivate * const self = clusterz->priv = nc_cluster_redshift_get_instance_private (clusterz);
  
  self->place_holder = 0;
}

static void
_nc_cluster_redshift_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_redshift_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_cluster_redshift, NC_TYPE_CLUSTER_REDSHIFT);

static gdouble _nc_cluster_redshift_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params);
static gdouble _nc_cluster_redshift_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z);
static gdouble _nc_cluster_redshift_intp_bin (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params);
static gboolean _nc_cluster_redshift_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng);
static void _nc_cluster_redshift_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
static void _nc_cluster_redshift_p_bin_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
static void _nc_cluster_redshift_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper);
static gdouble _nc_cluster_redshift_volume (NcClusterRedshift *clusterz);

static void
nc_cluster_redshift_class_init (NcClusterRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  
  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */
  
  object_class->finalize = _nc_cluster_redshift_finalize;
  
  ncm_model_class_set_name_nick (model_class, "Cluster redshift abstract class", "NcClusterRedshift");
  ncm_model_class_add_params (model_class, 0, 0, 1);
  
  ncm_mset_model_register_id (model_class,
                              "NcClusterRedshift",
                              "Cluster redshift observable models.",
                              NULL,
                              TRUE,
                              NCM_MSET_MODEL_MAIN);
  
  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));
  
  klass->P              = &_nc_cluster_redshift_p;
  klass->intP           = &_nc_cluster_redshift_intp;
  klass->intP_bin       = &_nc_cluster_redshift_intp_bin;
  klass->resample       = &_nc_cluster_redshift_resample;
  klass->P_limits       = &_nc_cluster_redshift_p_limits;
  klass->P_bin_limits   = &_nc_cluster_redshift_p_bin_limits;
  klass->N_limits       = &_nc_cluster_redshift_n_limits;
  klass->volume         = &_nc_cluster_redshift_volume;
  klass->obs_len        = 0;
  klass->obs_params_len = 0;
}

/* LCOV_EXCL_START */

static gdouble
_nc_cluster_redshift_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params)
{
  g_error ("_nc_cluster_redshift_p: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
  
  return 0.0;
}

static gdouble
_nc_cluster_redshift_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z)
{
  g_error ("_nc_cluster_redshift_intp: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
  
  return 0.0;
}

static gdouble
_nc_cluster_redshift_intp_bin (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params)
{
  g_error ("_nc_cluster_redshift_intp_bin: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
  
  return 0.0;
}

static gboolean
_nc_cluster_redshift_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng)
{
  g_error ("_nc_cluster_redshift_resample: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
  
  return FALSE;
}

static void
_nc_cluster_redshift_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  g_error ("_nc_cluster_redshift_p_limits: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
}

static void
_nc_cluster_redshift_p_bin_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  g_error ("_nc_cluster_redshift_p_limits: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
}

static void
_nc_cluster_redshift_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper)
{
  g_error ("_nc_cluster_redshift_n_limits: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
}

static gdouble
_nc_cluster_redshift_volume (NcClusterRedshift *clusterz)
{
  g_error ("_nc_cluster_redshift_volume: not implemented by `%s'\n", G_OBJECT_TYPE_NAME (clusterz));
  return 0.0;
}

/* LCOV_EXCL_STOP */

/**
 * nc_cluster_redshift_class_obs_len:
 * @clusterz_class: a #NcClusterRedshiftClass
 *
 * The number of observable redshifts of each cluster, e.g.,
 * 1 - only photometric redshift,
 * 1 - only spectroscopic redshift,
 * 2 - both photometric and spectroscopic redshifts.
 *
 * Returns: The number of observable redshifts.
 */
guint
nc_cluster_redshift_class_obs_len (NcClusterRedshiftClass *clusterz_class)
{
  return clusterz_class->obs_len;
}

/**
 * nc_cluster_redshift_class_obs_params_len:
 * @clusterz_class: a #NcClusterRedshiftClass
 *
 * The number of parameters related to the observable redshifts of each cluster, e.g.,
 * 1 - measured error of the photometric redshift.
 *
 * Returns: The number of parameters related to the observable redshifts.
 */
guint
nc_cluster_redshift_class_obs_params_len (NcClusterRedshiftClass *clusterz_class)
{
  return clusterz_class->obs_params_len;
}

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
  GObject *obj        = ncm_serialize_global_from_string (redshift_name);
  GType redshift_type = G_OBJECT_TYPE (obj);
  
  if (!g_type_is_a (redshift_type, NC_TYPE_CLUSTER_REDSHIFT))
    g_error ("nc_cluster_redshift_new_from_name: NcClusterRedshift %s do not descend from %s.", redshift_name, g_type_name (NC_TYPE_CLUSTER_REDSHIFT));
  
  return NC_CLUSTER_REDSHIFT (obj);
}

/**
 * nc_cluster_redshift_ref:
 * @clusterz: a #NcClusterRedshift
 *
 * Increases the reference count of @clusterz by one.
 *
 * Returns: (transfer full): @clusterz
 */
NcClusterRedshift *
nc_cluster_redshift_ref (NcClusterRedshift *clusterz)
{
  return g_object_ref (clusterz);
}

/**
 * nc_cluster_redshift_free:
 * @clusterz: a #NcClusterRedshift
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
 * @clusterz: a #NcClusterRedshift
 *
 * Atomically decrements the reference count of @clusterz by one. Set pointer to NULL.
 *
 */
void
nc_cluster_redshift_clear (NcClusterRedshift **clusterz)
{
  g_clear_object (clusterz);
}

/**
 * nc_cluster_redshift_obs_len:
 * @clusterz: a #NcClusterRedshift
 *
 * nc_cluster_redshift_obs_len().
 *
 * Returns: The number of observable redshifts.
 */
guint
nc_cluster_redshift_obs_len (NcClusterRedshift *clusterz)
{
  return nc_cluster_redshift_class_obs_len (NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz));
}

/**
 * nc_cluster_redshift_obs_params_len:
 * @clusterz: a #NcClusterRedshift
 *
 * See nc_cluster_redshift_class_obs_len().
 *
 * Returns: The number of parameters related to the observable redshifts.
 */
guint
nc_cluster_redshift_obs_params_len (NcClusterRedshift *clusterz)
{
  return nc_cluster_redshift_class_obs_params_len (NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz));
}

/**
 * nc_cluster_redshift_p:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z: true redshift
 * @lnM: true mass
 * @z_obs: (array) (element-type gdouble): measured redshift
 * @z_obs_params: (array) (element-type gdouble): measured redshift params
 *
 * It computes the probability density function (pdf) of the cluster redshift distribution @clusterz
 * given @z, @lnM and the measured redshit @z_obs and its parameter(s) @z_obs_params.
 *
 * Returns: The pdf of @clusterz.
 */
gdouble
nc_cluster_redshift_p (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs, const gdouble *z_obs_params)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->P (clusterz, cosmo, lnM, z, z_obs, z_obs_params);
}

/**
 * nc_cluster_redshift_intp:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z: true redshift
 * @lnM: true mass
 *
 * It computes the @clusterz probability distribution of @z lying
 * in the range $[z^{obs}_{min}, z^{obs}_{max}]$, namely,
 * $$ intp = \int_{z^{obs}_{min}}^{z^{obs}_{max}} p \, dz^{obs},$$
 * where $p$ is [nc_cluster_redshift_p()].
 *
 * Returns: The probability distribution of @z lying within $[z^{obs}_{min}, z^{obs}_{max}]$.
 */
gdouble
nc_cluster_redshift_intp (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->intP (clusterz, cosmo, lnM, z);
}

/**
 * nc_cluster_redshift_intp_bin:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z: true redshift
 * @lnM: true mass
 * @z_obs_lower: (array) (element-type gdouble): FIXME
 * @z_obs_upper: (array) (element-type gdouble): FIXME
 * @z_obs_params:(array) (element-type gdouble) (allow-none): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_cluster_redshift_intp_bin (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->intP_bin (clusterz, cosmo, lnM, z, z_obs_lower, z_obs_upper, z_obs_params);
}

/**
 * nc_cluster_redshift_resample:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z: true redshift
 * @lnM: true mass
 * @z_obs: (out): observed redshift
 * @z_obs_params: (out): observed redshift params
 *  @rng: a #NcmRNG
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_cluster_redshift_resample (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *z_obs, const gdouble *z_obs_params, NcmRNG *rng)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->resample (clusterz, cosmo, lnM, z, z_obs, z_obs_params, rng);
}

/**
 * nc_cluster_redshift_p_limits:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z_obs: (array) (element-type gdouble): observed redshift
 * @z_obs_params: (array) (element-type gdouble): observed redshift params
 * @z_lower: (out): pointer to the lower limit of the true redshift integration
 * @z_upper: (out): pointer to the upper limit of the true redshift integration
 *
 * FIXME
 *
 */
void
nc_cluster_redshift_p_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->P_limits (clusterz, cosmo, z_obs, z_obs_params, z_lower, z_upper);
}

/**
 * nc_cluster_redshift_p_bin_limits:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z_obs_lower: (array) (element-type gdouble): observed redshift
 * @z_obs_upper: (array) (element-type gdouble): observed redshift
 * @z_obs_params: (array) (element-type gdouble): observed redshift params
 * @z_lower: (out): pointer to the lower limit of the true redshift integration
 * @z_upper: (out): pointer to the upper limit of the true redshift integration
 *
 * FIXME
 *
 */
void
nc_cluster_redshift_p_bin_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->P_bin_limits (clusterz, cosmo, z_obs_lower, z_obs_upper, z_obs_params, z_lower, z_upper);
}

/**
 * nc_cluster_redshift_n_limits:
 * @clusterz: a #NcClusterRedshift
 * @cosmo: a #NcHICosmo
 * @z_lower: (out): pointer to the lower limit of the true redshift
 * @z_upper: (out): pointer to the upper limit of the true redshift
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @z_lower and @z_upper.
 */
void
nc_cluster_redshift_n_limits (NcClusterRedshift *clusterz, NcHICosmo *cosmo, gdouble *z_lower, gdouble *z_upper)
{
  NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->N_limits (clusterz, cosmo, z_lower, z_upper);
}

/**
 * nc_cluster_redshift_volume:
 * @clusterz: a #NcClusterRedshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_cluster_redshift_volume (NcClusterRedshift *clusterz)
{
  return NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz)->volume (clusterz);
}

static void
_nc_cluster_redshift_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);
  
  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType *modelsc = g_type_children (models[i], &ncc);
    
    g_message ("#  ");
    
    for (j = 0; j < n; j++)
      g_message (" ");
    
    g_message ("%s\n", g_type_name (models[i]));
    
    if (ncc)
      _nc_cluster_redshift_log_all_models_go (models[i], n + 2);
    
    g_free (modelsc);
  }
  
  g_free (models);
}

/**
 * nc_cluster_redshift_log_all_models:
 *
 * This function lists all implemented models of cluster redshift distributions.
 *
 */
void
nc_cluster_redshift_log_all_models (void)
{
  g_message ("# Registred NcClusterRedshift:%s are:\n", g_type_name (NC_TYPE_CLUSTER_REDSHIFT));
  _nc_cluster_redshift_log_all_models_go (NC_TYPE_CLUSTER_REDSHIFT, 0);
}

