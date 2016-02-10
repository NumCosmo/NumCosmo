/***************************************************************************
 *            nc_data_cmb_dist_priors.c
 *
 *  Thu Apr 22 15:56:44 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_data_cmb_dist_priors
 * @title: NcDataCMBDistPriors
 * @short_description: Cosmic microwave background data -- distance priors.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cmb_dist_priors.h"
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataCMBDistPriors, nc_data_cmb_dist_priors, NCM_TYPE_DATA_GAUSS);

static void
nc_data_cmb_dist_priors_init (NcDataCMBDistPriors *cmb_dist_priors)
{
  cmb_dist_priors->dist = NULL;
}

static void
_nc_data_cmb_dist_priors_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_cmb_dist_priors_parent_class)->constructed (object);
}

static void
nc_data_cmb_dist_priors_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataCMBDistPriors *cmb_dist_priors = NC_DATA_CMB_DIST_PRIORS (object);
  g_return_if_fail (NC_IS_DATA_CMB_DIST_PRIORS (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&cmb_dist_priors->dist);
      cmb_dist_priors->dist = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cmb_dist_priors_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataCMBDistPriors *cmb_dist_priors = NC_DATA_CMB_DIST_PRIORS (object);
  g_return_if_fail (NC_IS_DATA_CMB_DIST_PRIORS (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, cmb_dist_priors->dist);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cmb_dist_priors_dispose (GObject *object)
{
  NcDataCMBDistPriors *cmb_dist_priors = NC_DATA_CMB_DIST_PRIORS (object);

  nc_distance_clear (&cmb_dist_priors->dist);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cmb_dist_priors_parent_class)->dispose (object);
}

static void
nc_data_cmb_dist_priors_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cmb_dist_priors_parent_class)->finalize (object);
}

static void _nc_data_cmb_dist_priors_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cmb_dist_priors_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp);

static void
nc_data_cmb_dist_priors_class_init (NcDataCMBDistPriorsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_CLASS (klass);

  object_class->constructed  = &_nc_data_cmb_dist_priors_constructed;
  object_class->set_property = &nc_data_cmb_dist_priors_set_property;
  object_class->get_property = &nc_data_cmb_dist_priors_get_property;
  object_class->dispose      = &nc_data_cmb_dist_priors_dispose;
  object_class->finalize     = &nc_data_cmb_dist_priors_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_cmb_dist_priors_prepare;
  gauss_class->mean_func = &_nc_data_cmb_dist_priors_mean_func;
}

static void
_nc_data_cmb_dist_priors_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataCMBDistPriors *cmb_dist_priors = NC_DATA_CMB_DIST_PRIORS (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (cmb_dist_priors->dist, cosmo);
}

static void 
_nc_data_cmb_dist_priors_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataCMBDistPriors *cmb_dist_priors = NC_DATA_CMB_DIST_PRIORS (gauss);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  
  const gdouble Ascale = nc_distance_acoustic_scale (cmb_dist_priors->dist, cosmo);
  const gdouble Rlss = nc_distance_shift_parameter_lss (cmb_dist_priors->dist, cosmo);
  const gdouble zdec = nc_distance_decoupling_redshift (cmb_dist_priors->dist, cosmo);

  ncm_vector_set (vp, 0, Ascale);
  ncm_vector_set (vp, 1, Rlss);
  ncm_vector_set (vp, 2, zdec);
}

/**
 * nc_data_cmb_dist_priors_new_empty:
 * @dist: a #NcDistance
 *
 * This function allocates an empty #NcDataCMBDistPriors object.
 *
 * Returns: A #NcDataCMBDistPriors.
 */
NcDataCMBDistPriors *
nc_data_cmb_dist_priors_new_empty (NcDistance *dist)
{
  NcDataCMBDistPriors *cmb_dist_priors = g_object_new (NC_TYPE_DATA_CMB_DIST_PRIORS,
                                                       "dist", dist,
                                                       NULL);
  return cmb_dist_priors;
}

/**
 * nc_data_cmb_dist_priors_new_from_file:
 * @filename: file containing a serialized #NcDataCMBDistPriors
 * 
 * Creates a new #NcDataCMBDistPriors from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataCMBDistPriors.
 */
NcDataCMBDistPriors *
nc_data_cmb_dist_priors_new_from_file (const gchar *filename)
{
  NcDataCMBDistPriors *cmb_dist_priors = NC_DATA_CMB_DIST_PRIORS (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_CMB_DIST_PRIORS (cmb_dist_priors));

  return cmb_dist_priors;
}


/**
 * nc_data_cmb_dist_priors_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataCMBId
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataCMBDistPriors *
nc_data_cmb_dist_priors_new_from_id (NcDistance *dist, NcDataCMBId id)
{
  NcDataCMBDistPriors *cmb_dist_priors;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_CMB_DIST_PRIORS_WMAP5:
      filename = ncm_cfg_get_data_filename ("nc_data_cmb_wmap5_dist_priors.obj", TRUE);
      break;
    case NC_DATA_CMB_DIST_PRIORS_WMAP7:
      filename = ncm_cfg_get_data_filename ("nc_data_cmb_wmap7_dist_priors.obj", TRUE);
      break;
    case NC_DATA_CMB_DIST_PRIORS_WMAP9:
      filename = ncm_cfg_get_data_filename ("nc_data_cmb_wmap9_dist_priors.obj", TRUE);
      break;
    default:
      g_error ("nc_data_cmb_dist_priors_new_from_id: id %d not recognized.", id);
      break;
  }

  cmb_dist_priors = nc_data_cmb_dist_priors_new_from_file (filename);
  nc_data_cmb_dist_priors_set_dist (cmb_dist_priors, dist);
  g_free (filename);

  return cmb_dist_priors;
}

/**
 * nc_data_cmb_dist_priors_set_dist:
 * @cmb_dist_priors: a #NcDataCMBDistPriors
 * @dist: a #NcDistance
 * 
 * Sets the distance object.
 * 
 */
void 
nc_data_cmb_dist_priors_set_dist (NcDataCMBDistPriors *cmb_dist_priors, NcDistance *dist)
{
  nc_distance_clear (&cmb_dist_priors->dist);
  cmb_dist_priors->dist = nc_distance_ref (dist);
}
