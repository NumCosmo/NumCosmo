/***************************************************************************
 *            nc_data_cluster_wl.c
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.c
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_cluster_wl
 * @title: NcDataClusterWL
 * @short_description: Cluster weak lensing likelihood.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_wl.h"

#include "nc_hicosmo.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "lss/nc_galaxy_wl.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"

struct _NcDataClusterWLPrivate
{
  NcmObjArray *galaxy_array;
  gdouble psf_size;
  gdouble z_cluster;
  gdouble ra_cluster;
  gdouble dec_cluster;
};

enum
{
  PROP_0,
  PROP_GALAXY_ARRAY,
  PROP_PSF_SIZE,
  PROP_Z_CLUSTER,
  PROP_RA_CLUSTER,
  PROP_DEC_CLUSTER,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterWL, nc_data_cluster_wl, NCM_TYPE_DATA);

static void
nc_data_cluster_wl_init (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = dcwl->priv = nc_data_cluster_wl_get_instance_private (dcwl);
  
  self->galaxy_array = ncm_obj_array_new ();
  self->z_cluster    = 0.0;
  self->ra_cluster   = 0.0;
  self->dec_cluster  = 0.0;
  self->psf_size     = 0.0;
}

static void
nc_data_cluster_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  
  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));
  
  switch (prop_id)
  {
    case PROP_GALAXY_ARRAY:
    {
      NcmObjArray *galaxy_array = g_value_get_boxed (value);
      
      g_clear_pointer (&self->galaxy_array, ncm_obj_array_unref);
      self->galaxy_array = ncm_obj_array_ref (galaxy_array);
      
      break;
    }
    case PROP_PSF_SIZE:
      self->psf_size = g_value_get_double (value);
      break;
    case PROP_Z_CLUSTER:
      self->z_cluster = g_value_get_double (value);
      break;
    case PROP_RA_CLUSTER:
      self->ra_cluster = g_value_get_double (value);
      break;
    case PROP_DEC_CLUSTER:
      self->dec_cluster = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_wl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  
  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));
  
  switch (prop_id)
  {
    case PROP_GALAXY_ARRAY:
      g_value_set_boxed (value, self->galaxy_array);
      break;
    case PROP_PSF_SIZE:
      g_value_set_double (value, self->psf_size);
      break;
    case PROP_Z_CLUSTER:
      g_value_set_double (value, self->z_cluster);
      break;
    case PROP_RA_CLUSTER:
      g_value_set_double (value, self->ra_cluster);
      break;
    case PROP_DEC_CLUSTER:
      g_value_set_double (value, self->dec_cluster);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_wl_dispose (GObject *object)
{
  NcDataClusterWL *dcwl                = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  
  ncm_obj_array_clear (&self->galaxy_array);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->dispose (object);
}

static void
nc_data_cluster_wl_finalize (GObject *object)
{
  /*NcDataClusterWL *dcwl = NC_DATA_CLUSTER_WL (object);*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->finalize (object);
}

static void _nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_cluster_wl_get_len (NcmData *data);
static void _nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset);

static void
nc_data_cluster_wl_class_init (NcDataClusterWLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  
  object_class->set_property = nc_data_cluster_wl_set_property;
  object_class->get_property = nc_data_cluster_wl_get_property;
  object_class->dispose      = nc_data_cluster_wl_dispose;
  object_class->finalize     = nc_data_cluster_wl_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_GALAXY_ARRAY,
                                   g_param_spec_boxed ("galaxy-array",
                                                       NULL,
                                                       "Array of galaxy weak lensing objects",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PSF_SIZE,
                                   g_param_spec_double ("psf-size",
                                                        NULL,
                                                        "PSF size",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_double ("z-cluster",
                                                        NULL,
                                                        "Cluster (halo) redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RA_CLUSTER,
                                   g_param_spec_double ("ra-cluster",
                                                        NULL,
                                                        "Cluster (halo) RA",
                                                        -360.0, 360.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DEC_CLUSTER,
                                   g_param_spec_double ("dec-cluster",
                                                        NULL,
                                                        "Cluster (halo) DEC",
                                                        -360.0, 360.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->m2lnL_val  = &_nc_data_cluster_wl_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_wl_get_len;
  data_class->prepare    = &_nc_data_cluster_wl_prepare;
}

static void
_nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcHICosmo *cosmo                    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd         = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const guint ngal                    = self->galaxy_array->len;
  gint i;

  m2lnL[0] = 0.0;

  for (i = 0; i < ngal; i++)
  {
    NcGalaxyWL *gwl_i = NC_GALAXY_WL (ncm_obj_array_peek (self->galaxy_array, i));
    m2lnL[0] += nc_galaxy_wl_eval_m2lnP (gwl_i, cosmo, dp, smd, self->z_cluster);
  }
  
  return;
}

static guint
_nc_data_cluster_wl_get_len (NcmData *data)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  
  if (self->galaxy_array != NULL)
  {

    const guint ngal = self->galaxy_array->len;
    guint len = 0;
    gint i;

    for (i = 0; i < ngal; i++)
    {
      NcGalaxyWL *gwl_i = NC_GALAXY_WL (ncm_obj_array_peek (self->galaxy_array, i));
      len += nc_galaxy_wl_len (gwl_i);
    }
    return len;
  }
  else
    return 0;
}

static void
_nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset)
{
  /*NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);*/
  /*NcDataClusterWLPrivate * const self = dcwl->priv;*/
  NcHICosmo *cosmo                    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd         = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  
  g_assert ((cosmo != NULL) && (smd != NULL) && (dp != NULL));
  
  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);
}

/**
 * nc_data_cluster_wl_new:
 *
 * Creates a new #NcDataClusterWL.
 *
 * Returns: the newly created #NcDataClusterWL.
 */
NcDataClusterWL *
nc_data_cluster_wl_new (void)
{
  NcDataClusterWL *dcwl = g_object_new (NC_TYPE_DATA_CLUSTER_WL,
                                        NULL);
  
  return dcwl;
}

/**
 * nc_data_cluster_wl_new_from_file:
 * @filename: file containing a serialized #NcDataClusterWL
 *
 * Creates a new #NcDataClusterWL from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataClusterWL.
 */
NcDataClusterWL *
nc_data_cluster_wl_new_from_file (const gchar *filename)
{
  NcDataClusterWL *dcwl = NC_DATA_CLUSTER_WL (ncm_serialize_global_from_file (filename));
  
  g_assert (NC_IS_DATA_CLUSTER_WL (dcwl));
  
  return dcwl;
}

/**
 * nc_data_cluster_wl_ref:
 * @dcwl: a #NcDataClusterWL
 *
 * Increases the reference count of @dcwl by one.
 *
 * Returns: (transfer full): @dcwl
 */
NcDataClusterWL *
nc_data_cluster_wl_ref (NcDataClusterWL *dcwl)
{
  return g_object_ref (dcwl);
}

/**
 * nc_data_cluster_wl_free:
 * @dcwl: a #NcDataClusterWL
 *
 * Atomically decrements the reference count of @dcwl by one. If the reference count drops to 0,
 * all memory allocated by @dcwl is released.
 *
 */
void
nc_data_cluster_wl_free (NcDataClusterWL *dcwl)
{
  g_object_unref (dcwl);
}

/**
 * nc_data_cluster_wl_clear:
 * @dcwl: a #NcDataClusterWL
 *
 * The reference count of @dcwl is decreased and the pointer is set to NULL.
 *
 */
void
nc_data_cluster_wl_clear (NcDataClusterWL **dcwl)
{
  g_clear_object (dcwl);
}
