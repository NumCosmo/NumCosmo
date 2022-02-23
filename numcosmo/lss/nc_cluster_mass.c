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
 * @title: NcClusterMass
 * @short_description: Abstract class for cluster mass distributions.
 *
 * NcClusterMass is the abstract class designed to abridge the functions
 * that any cluster mass distribution should implement, see NcClusterMassImpl.
 * Its parent_class is NcmModel.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

struct _NcClusterMassPrivate
{
  guint nbins;
  GPtrArray *bin_knots;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcClusterMass, nc_cluster_mass, NCM_TYPE_MODEL);

static void
nc_cluster_mass_init (NcClusterMass *clusterm)
{
  NcClusterMassPrivate * const self = clusterm->priv = nc_cluster_mass_get_instance_private (clusterm);

  self->nbins = 0;
  self->bin_knots = g_ptr_array_new ();

  g_ptr_array_set_free_func (self->bin_knots, ncm_vector_free);
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
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = nc_cluster_mass_finalize;

  ncm_model_class_set_name_nick (model_class, "Cluster mass abstract class", "NcClusterMass");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_mset_model_register_id (model_class,
                              "NcClusterMass",
                              "Cluster mass observable relation models.",
                              NULL,
                              TRUE,
                              NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));
}

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
  GObject *obj    = ncm_serialize_global_from_string (mass_name);
  GType mass_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (mass_type, NC_TYPE_CLUSTER_MASS))
    g_error ("nc_cluster_mass_new_from_name: NcClusterMass %s do not descend from %s.", 
             mass_name, g_type_name (NC_TYPE_CLUSTER_MASS));

  return NC_CLUSTER_MASS (obj);
}

/**
 * nc_cluster_mass_ref:
 * @clusterm: a #NcClusterMass
 *
 * Increases the reference count of @clusterm by one.
 *
 * Returns: (transfer full): @clusterm
 */
NcClusterMass *
nc_cluster_mass_ref (NcClusterMass *clusterm)
{
  return g_object_ref (clusterm);
}

/**
 * nc_cluster_mass_free:
 * @clusterm: a #NcClusterMass
 *
 * Atomically decrements the reference count of @clusterm by one. If the reference count drops to 0,
 * all memory allocated by @clusterm is released.
 *
 */
void
nc_cluster_mass_free (NcClusterMass *clusterm)
{
  g_object_unref (clusterm);
}

/**
 * nc_cluster_mass_clear:
 * @clusterm: a #NcClusterMass
 *
 * The reference count of @clusterm is decreased and the pointer is set to NULL.
 *
 */
void
nc_cluster_mass_clear (NcClusterMass **clusterm)
{
  g_clear_object (clusterm);
}

/**
 * nc_cluster_mass_obs_len:
 * @clusterm: a #NcClusterMass
 *
 * The number of observable masses (or just the observable which is related to the cluster mass)
 * of each cluster, e.g.,
 * 1 - SZ mass,
 * 1 - X-ray mass,
 * 1 - Lensing mass,
 * 2 - SZ and X-ray masses,
 * 3 - SZ, X-ray and lensing masses.
 *
 * Returns: The number of observable masses.
 */
guint
nc_cluster_mass_obs_len (NcClusterMass *clusterm)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->obs_len (clusterm);
}

/**
 * nc_cluster_mass_obs_params_len:
 * @clusterm: a #NcClusterMass
 *
 * The number of parameters related to the observable masses of each cluster, e.g.,
 * 1 - error of the SZ mass,
 * 1 - error of the X-ray mass,
 * 2 - errors of SZ and X-ray masses.
 *
 * Returns: The number of parameters related to the observable masses.
 */
guint
nc_cluster_mass_obs_params_len (NcClusterMass *clusterm)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->obs_params_len (clusterm);
}

/**
 * nc_cluster_mass_set_bins:
 * @clusterm: a #NcClusterMass
 * @lnM_obs_knots: (array) (element-type NcmVector): An array of #NcmVector's containing the knots for each dimension.
 *
 * Receives an array of nc_cluster_mass_obs_len() dimension containing
 * the bin knots for each dimension.
 *
 */
void
nc_cluster_mass_set_bins (NcClusterMass *clusterm, GPtrArray *lnM_obs_knots)
{
  NcClusterMassPrivate * const self = clusterm->priv;
  const guint dim = nc_cluster_mass_obs_len (clusterm);
  gint i;
  guint nbins = 1;

  g_assert_cmpuint (lnM_obs_knots->len, ==, dim);

  g_ptr_array_set_size (self->bin_knots, 0);

  for (i = 0; i < dim; i++)
  {
    NcmVector *knots_i = g_ptr_array_index (lnM_obs_knots, i);
    const guint knots_i_len = ncm_vector_len (knots_i);

    g_assert_cmpuint (knots_i_len, >, 1);

    nbins *= knots_i_len - 1;

    g_ptr_array_append (self->bin_knots, ncm_vector_dup (knots_i));
  }

  self->nbins = nbins;
}

/**
 * nc_cluster_mass_peek_bins:
 * @clusterm: a #NcClusterMass
 *
 * Provides the current binning returning the array of
 * knots for each dimension.
 *
 * Returns: (transfer none): an array of #NcmVector's
 */
GPtrArray *
nc_cluster_mass_peek_bins (NcClusterMass *clusterm)
{
  NcClusterMassPrivate * const self = clusterm->priv;

  return self->bin_knots;
}


/**
 * nc_cluster_mass_get_bin:
 * @clusterm: a #NcClusterMass
 * @bin_index: bin index
 * @lnM_obs_lower: (out callee-allocates) (array) (element-type gdouble): bin lower bounds
 * @lnM_obs_upper: (out callee-allocates) (array) (element-type gdouble): bin upper bounds
 *
 * FIXME
 *
 */
void
nc_cluster_mass_get_bin (NcClusterMass *clusterm, const guint bin_index, gdouble **lnM_obs_lower, gdouble **lnM_obs_upper)
{
  NcClusterMassPrivate * const self = clusterm->priv;

  g_assert_cmpuint (bin_index, <, self->nbins);



}

gboolean nc_cluster_mass_has_bins (NcClusterMass *clusterm);


/**
 * nc_cluster_mass_p:
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo
 * @lnM:logarithm base e of the true mass
 * @z: true redshift
 * @lnM_obs: (array) (element-type double): logarithm base e of the observed mass
 * @lnM_obs_params: (array) (element-type double): observed mass paramaters
 *
 * It computes the probability density function (pdf) of the cluster mass distribution @clusterm
 * given @cosmo, @lnM, @z and the observable cluster mass (or just the observable).
 *
 * Returns: The pdf of @clusterm.
*/
gdouble
nc_cluster_mass_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->P (clusterm, cosmo, lnM, z, lnM_obs, lnM_obs_params);
}

/**
 * nc_cluster_mass_intp:
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo
 * @z: true redshift
 * @lnM: logarithm base e of the true mass
 *
 * It computes the @clusterm probability distribution of @lnM lying
 * in the range $[]$, namely,
 * $$ intp = \int_{\ln M^{obs}_{min}}^{\ln M^{obs}_{max}} p \, d\ln M^{obs},$$
 * where $p$ is [nc_cluster_mass_p()].
 *
 * Returns: The probability distribution of @lnM lying within $[\ln M^{obs}_{min}, \ln M^{obs}_{max}]$.
*/
gdouble
nc_cluster_mass_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->intP (clusterm, cosmo, lnM, z);
}

/**
 * nc_cluster_mass_resample:
 * @clusterm: a #NcClusterMass.
 * @cosmo: a #NcHICosmo.
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
nc_cluster_mass_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  return NC_CLUSTER_MASS_GET_CLASS (clusterm)->resample (clusterm, cosmo, lnM, z, lnM_obs, lnM_obs_params, rng);
}

/**
 * nc_cluster_mass_p_limits:
 * @clusterm: a #NcClusterMass.
 * @cosmo: a #NcHICosmo.
 * @lnM_obs: (array) (element-type gdouble): observed mass.
 * @lnM_obs_params: (array) (element-type gdouble): observed mass params.
 * @lnM_lower: (out): pointer to the lower limit of the real mass integration.
 * @lnM_upper: (out): pointer to the upper limit of the real mass integration.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @lnM_lower and @lnM_upper.
 */
void
nc_cluster_mass_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NC_CLUSTER_MASS_GET_CLASS (clusterm)->P_limits (clusterm, cosmo, lnM_obs, lnM_obs_params, lnM_lower, lnM_upper);
}

/**
 * nc_cluster_mass_n_limits:
 * @clusterm: a #NcClusterMass.
 * @cosmo: a #NcHICosmo.
 * @lnM_lower: (out): lower limit of the logarithm base e of the true mass.
 * @lnM_upper: (out): upper limit of the logarithm base e of the true mass.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @lnM_lower and @lnM_upper.
 */
void
nc_cluster_mass_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NC_CLUSTER_MASS_GET_CLASS (clusterm)->N_limits (clusterm, cosmo, lnM_lower, lnM_upper);
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
 * This function lists all implemented models of cluster mass distributions.
 *
 */
void
nc_cluster_mass_log_all_models (void)
{
  g_message ("# Registred NcClusterMass:%s are:\n", g_type_name (NC_TYPE_CLUSTER_MASS));
  _nc_cluster_mass_log_all_models_go (NC_TYPE_CLUSTER_MASS, 0);
}
