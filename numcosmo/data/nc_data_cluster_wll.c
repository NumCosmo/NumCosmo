/***************************************************************************
 *            nc_data_cluster_wll.c
 *
 *  Tue Jun 15 16:00:13 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wll.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_data_cluster_wll
 * @title: NcDataClusterWLL
 * @short_description: Cluster weak lensing likelihood.
 * @stability: Unstable
 *
 * This class implements the weak lensing likelihood for galaxy clusters.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_wll.h"

#include "nc_hicosmo.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "galaxy/nc_galaxy_wl_likelihood.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"

struct _NcDataClusterWLLPrivate
{
  NcmObjArray *galaxy_array;
  gdouble z_cluster;
  gboolean kde;
};

enum
{
  PROP_0,
  PROP_GALAXY_ARRAY,
  PROP_Z_CLUSTER,
  PROP_KDE,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterWLL, nc_data_cluster_wll, NCM_TYPE_DATA);

static void
nc_data_cluster_wll_init (NcDataClusterWLL *dcwll)
{
  NcDataClusterWLLPrivate * const self = dcwll->priv = nc_data_cluster_wll_get_instance_private (dcwll);

  self->galaxy_array = ncm_obj_array_new ();
  self->z_cluster    = 0.0;
  self->kde          = FALSE;
}

static void
nc_data_cluster_wll_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (object);
  NcDataClusterWLLPrivate * const self = dcwll->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_WLL (object));

  switch (prop_id)
  {
    case PROP_GALAXY_ARRAY:
    {
      NcmObjArray *galaxy_array = g_value_get_boxed (value);

      g_clear_pointer (&self->galaxy_array, ncm_obj_array_unref);
      self->galaxy_array = ncm_obj_array_ref (galaxy_array);

      break;
    }
    case PROP_Z_CLUSTER:
      self->z_cluster = g_value_get_double (value);
      break;
    case PROP_KDE:
      nc_data_cluster_wll_set_kde (dcwll, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_wll_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (object);
  NcDataClusterWLLPrivate * const self = dcwll->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_WLL (object));

  switch (prop_id)
  {
    case PROP_GALAXY_ARRAY:
      g_value_set_boxed (value, self->galaxy_array);
      break;
    case PROP_Z_CLUSTER:
      g_value_set_double (value, self->z_cluster);
      break;
    case PROP_KDE:
      g_value_set_boolean (value, self->kde);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_wll_dispose (GObject *object)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (object);
  NcDataClusterWLLPrivate * const self = dcwll->priv;

  ncm_obj_array_clear (&self->galaxy_array);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wll_parent_class)->dispose (object);
}

static void
nc_data_cluster_wll_finalize (GObject *object)
{
  /*NcDataClusterWLL *dcwll = NC_DATA_CLUSTER_WLL (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wll_parent_class)->finalize (object);
}

static void _nc_data_cluster_wll_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_cluster_wll_get_len (NcmData *data);
static void _nc_data_cluster_wll_prepare (NcmData *data, NcmMSet *mset);

static void
nc_data_cluster_wll_class_init (NcDataClusterWLLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_wll_set_property;
  object_class->get_property = nc_data_cluster_wll_get_property;
  object_class->dispose      = nc_data_cluster_wll_dispose;
  object_class->finalize     = nc_data_cluster_wll_finalize;

  g_object_class_install_property (object_class,
                                   PROP_GALAXY_ARRAY,
                                   g_param_spec_boxed ("galaxy-array",
                                                       NULL,
                                                       "Array of galaxy weak lensing objects",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_double ("z-cluster",
                                                        NULL,
                                                        "Cluster (halo) redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_KDE,
                                   g_param_spec_boolean ("kde",
                                                         NULL,
                                                         "Whether to use KDE method",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  data_class->m2lnL_val  = &_nc_data_cluster_wll_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_wll_get_len;
  data_class->prepare    = &_nc_data_cluster_wll_prepare;
}

static void
_nc_data_cluster_wll_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (data);
  NcDataClusterWLLPrivate * const self = dcwll->priv;
  NcHICosmo *cosmo                     = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd          = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp             = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const guint ngal                     = self->galaxy_array->len;
  guint i;

  m2lnL[0] = 0.0;

  for (i = 0; i < ngal; i++)
  {
    NcGalaxyWLLikelihood *gwll_i = NC_GALAXY_WL_LIKELIHOOD (ncm_obj_array_peek (self->galaxy_array, i));

    if (self->kde)
      m2lnL[0] += nc_galaxy_wl_likelihood_kde_eval_m2lnP (gwll_i, cosmo, dp, smd, self->z_cluster);
    else
      m2lnL[0] += nc_galaxy_wl_likelihood_eval_m2lnP (gwll_i, cosmo, dp, smd, self->z_cluster);
  }

  return;
}

static guint
_nc_data_cluster_wll_get_len (NcmData *data)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (data);
  NcDataClusterWLLPrivate * const self = dcwll->priv;

  if (self->galaxy_array != NULL)
  {
    const guint ngal = self->galaxy_array->len;
    guint len        = 0;
    guint i;

    for (i = 0; i < ngal; i++)
    {
      NcGalaxyWLLikelihood *gwll_i = NC_GALAXY_WL_LIKELIHOOD (ncm_obj_array_peek (self->galaxy_array, i));

      len += nc_galaxy_wl_likelihood_len (gwll_i);
    }

    return len;
  }
  else
  {
    return 0;
  }
}

static void
_nc_data_cluster_wll_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (data);
  NcDataClusterWLLPrivate * const self = dcwll->priv;
  NcHICosmo *cosmo                     = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd          = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp             = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));

  g_assert ((cosmo != NULL) && (smd != NULL) && (dp != NULL));

  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);

  if (self->kde)
    printf ("KDE\n");
  else
    printf ("No KDE\n");
}

/**
 * nc_data_cluster_wll_new:
 *
 * Creates a new #NcDataClusterWLL.
 *
 * Returns: the newly created #NcDataClusterWLL.
 */
NcDataClusterWLL *
nc_data_cluster_wll_new (void)
{
  NcDataClusterWLL *dcwll = g_object_new (NC_TYPE_DATA_CLUSTER_WLL,
                                          NULL);

  return dcwll;
}

/**
 * nc_data_cluster_wll_new_from_file:
 * @filename: file containing a serialized #NcDataClusterWLL
 *
 * Creates a new #NcDataClusterWLL from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataClusterWLL.
 */
NcDataClusterWLL *
nc_data_cluster_wll_new_from_file (const gchar *filename)
{
  NcDataClusterWLL *dcwll = NC_DATA_CLUSTER_WLL (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_CLUSTER_WLL (dcwll));

  return dcwll;
}

/**
 * nc_data_cluster_wll_ref:
 * @dcwll: a #NcDataClusterWLL
 *
 * Increases the reference count of @dcwll by one.
 *
 * Returns: (transfer full): @dcwll
 */
NcDataClusterWLL *
nc_data_cluster_wll_ref (NcDataClusterWLL *dcwll)
{
  return g_object_ref (dcwll);
}

/**
 * nc_data_cluster_wll_free:
 * @dcwll: a #NcDataClusterWLL
 *
 * Atomically decrements the reference count of @dcwll by one. If the reference count drops to 0,
 * all memory allocated by @dcwll is released.
 *
 */
void
nc_data_cluster_wll_free (NcDataClusterWLL *dcwll)
{
  g_object_unref (dcwll);
}

/**
 * nc_data_cluster_wll_clear:
 * @dcwll: a #NcDataClusterWLL
 *
 * The reference count of @dcwll is decreased and the pointer is set to NULL.
 *
 */
void
nc_data_cluster_wll_clear (NcDataClusterWLL **dcwll)
{
  g_clear_object (dcwll);
}

/**
 * nc_data_cluster_wll_set_kde:
 * @dcwll: a #NcDataClusterWLL
 * @kde: whether to use KDE method
 *
 * The reference count of @dcwll is decreased and the pointer is set to NULL.
 *
 */
void
nc_data_cluster_wll_set_kde (NcDataClusterWLL *data, gboolean kde)
{
  NcDataClusterWLL *dcwll              = NC_DATA_CLUSTER_WLL (data);
  NcDataClusterWLLPrivate * const self = dcwll->priv;

  self->kde = kde;
}

