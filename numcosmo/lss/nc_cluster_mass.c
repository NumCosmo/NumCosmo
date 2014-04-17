/***************************************************************************
 *            nc_cluster_mass.c
 *
 *  Thu June 21 23:27:03 2012
 *  Copyright  2012  Mariana Penna Lima, Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com> 
 * Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_cluster_mass
 * @title: Cluster Mass Distribution
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcClusterMass, nc_cluster_mass, NCM_TYPE_MODEL);

/**
 * nc_cluster_mass_new_from_name: 
 * @mass_name: string which specifies the type of the mass distribution.
 *
 * This function returns a new #NcClusterMass whose type is defined by @mass_name.
 *
 * Returns: A new #NcClusterMass.
 */
NcClusterMass *
nc_cluster_mass_new_from_name (gchar *mass_name)
{
  GObject *obj = ncm_serialize_global_from_string (mass_name);
  GType mass_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (mass_type, NC_TYPE_CLUSTER_MASS))
	g_error ("nc_cluster_mass_new_from_name: NcClusterMass %s do not descend from %s\n", mass_name, g_type_name (NC_TYPE_CLUSTER_MASS));
  return NC_CLUSTER_MASS (obj);
}

/**
 * nc_cluster_mass_ref:
 * @clusterm: a #NcClusterMass.
 *
 * FIXME
 *
 * Returns: (transfer full): @clusterm.
 */
NcClusterMass *
nc_cluster_mass_ref (NcClusterMass *clusterm)
{
  return g_object_ref (clusterm);
}

/**
 * nc_cluster_mass_free:
 * @clusterm: a #NcClusterMass.
 *
 * FIXME
 *
 */
void
nc_cluster_mass_free (NcClusterMass *clusterm)
{
  g_object_unref (clusterm);
}

/**
 * nc_cluster_mass_clear:
 * @clusterm: a #NcClusterMass.
 *
 * FIXME
 *
 */
void
nc_cluster_mass_clear (NcClusterMass **clusterm)
{
  g_clear_object (clusterm);
}

/**
 * nc_cluster_mass_impl:
 * @clusterm: a #NcClusterMass.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcClusterMassImpl 
nc_cluster_mass_impl (NcClusterMass *clusterm)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->impl;
}

/**
 * nc_cluster_mass_obs_len:
 * @clusterm: a #NcClusterMass.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
nc_cluster_mass_obs_len (NcClusterMass *clusterm)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->obs_len (clusterm);
}

/**
 * nc_cluster_mass_obs_params_len:
 * @clusterm: a #NcClusterMass.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
nc_cluster_mass_obs_params_len (NcClusterMass *clusterm)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->obs_params_len (clusterm);
}

/**
 * nc_cluster_mass_p:
 * @clusterm: a #NcClusterMass.
 * @model: a #NcHICosmo. 
 * @z: true redshift.
 * @lnM: logarithm base e of the true mass.
 * @lnM_obs: logarithm base e of the observed mass.
 * @lnM_obs_params: observed mass params.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_mass_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->P (clusterm, model, lnM, z, lnM_obs, lnM_obs_params);
}

/**
 * nc_cluster_mass_intp:
 * @clusterm: a #NcClusterMass.
 * @model: a #NcHICosmo.
 * @z: true redshift.
 * @lnM: logarithm base e of the true mass.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_mass_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->intP (clusterm, model, lnM, z);
}

/**
 * nc_cluster_mass_resample:
 * @clusterm: a #NcClusterMass.
 * @model: a #NcHICosmo. 
 * @z: true redshift.
 * @lnM: logarithm base e of the true mass.
 * @lnM_obs: (out): logarithm base e of the observed mass.
 * @lnM_obs_params: (out): observed mass params.
 * @rng: a #NcmRNG.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gboolean
nc_cluster_mass_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params, NcmRNG *rng)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->resample (clusterm,model, lnM, z, lnM_obs, lnM_obs_params, rng);
}

/**
 * nc_cluster_mass_p_limits:
 * @clusterm: a #NcClusterMass.
 * @model: a #NcHICosmo. 
 * @lnM_obs: observed mass.
 * @lnM_obs_params: observed mass params.
 * @lnM_lower: (out): pointer to the lower limit of the real mass integration.
 * @lnM_upper: (out): pointer to the upper limit of the real mass integration.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @lnM_lower and @lnM_upper.
 */
void
nc_cluster_mass_p_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NC_CLUSTER_MASS_GET_CLASS (clusterm)->P_limits (clusterm, model, lnM_obs, lnM_obs_params, lnM_lower, lnM_upper);
}

/**
 * nc_cluster_mass_n_limits:
 * @clusterm: a #NcClusterMass.
 * @model: a #NcHICosmo. 
 * @lnM_lower: (out): lower limit of the logarithm base e of the true mass.
 * @lnM_upper: (out): upper limit of the lgarithm base e of the true mass.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @lnM_lower and @lnM_upper.
 */
void
nc_cluster_mass_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NC_CLUSTER_MASS_GET_CLASS (clusterm)->N_limits (clusterm, model, lnM_lower, lnM_upper);
}

static void
_nc_cluster_mass_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);
  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType *modelsc = g_type_children (models[i], &ncc);

    g_message ("#  ");
    for (j = 0; j < n; j++) g_message (" ");
    g_message ("%s\n", g_type_name (models[i]));
    if (ncc)
      _nc_cluster_mass_log_all_models_go (models[i], n + 2);

    g_free (modelsc);
  }
  g_free (models);
}

/**
 * nc_cluster_mass_log_all_models:
 *
 * FIXME
 *
 */
void
nc_cluster_mass_log_all_models ()
{
  g_message ("# Registred NcClusterMass:%s are:\n", g_type_name (NC_TYPE_CLUSTER_MASS));
  _nc_cluster_mass_log_all_models_go (NC_TYPE_CLUSTER_MASS, 0);
}

static void
nc_cluster_mass_init (NcClusterMass *nc_cluster_mass)
{
  NCM_UNUSED (nc_cluster_mass);
}

static void
nc_cluster_mass_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_cluster_mass, NC_TYPE_CLUSTER_MASS);

static void
nc_cluster_mass_class_init (NcClusterMassClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = nc_cluster_mass_finalize;

  ncm_model_class_add_params (NCM_MODEL_CLASS (klass), 0, 0, 1);
  
  ncm_mset_model_register_id (NCM_MODEL_CLASS (klass), 
                              "NcClusterMass",
                              "Cluster mass observable relation models.",
                              NULL);
  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));
}

