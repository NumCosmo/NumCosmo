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
 * SECTION:nc_data_cmb_shift_param
 * @title: NcDataCMBShiftParam
 * @short_description: Cosmic microwave background data -- shift parameter.
 * @stability: Stable
 * @include: numcosmo/data/nc_data_cmb_shift_param.h
 *
 * This object creates a CMB shift parameter $R$ data (#NcmData).
 * See #NcDataCMBId and #nc_distance_shift_parameter for more informations.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cmb_shift_param.h"
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataCMBShiftParam, nc_data_cmb_shift_param, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_cmb_shift_param_init (NcDataCMBShiftParam *cmb_shift_param)
{
  cmb_shift_param->dist = NULL;
  cmb_shift_param->x    = NULL;
}

static void
nc_data_cmb_shift_param_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (object);
  
  g_return_if_fail (NC_IS_DATA_CMB_SHIFT_PARAM (object));
  
  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_cmb_shift_param_set_dist (cmb_shift_param, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&cmb_shift_param->x, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cmb_shift_param_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (object);
  
  g_return_if_fail (NC_IS_DATA_CMB_SHIFT_PARAM (object));
  
  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, cmb_shift_param->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, cmb_shift_param->x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cmb_shift_param_dispose (GObject *object)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (object);
  
  nc_distance_clear (&cmb_shift_param->dist);
  ncm_vector_clear (&cmb_shift_param->x);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cmb_shift_param_parent_class)->dispose (object);
}

static void
nc_data_cmb_shift_param_finalize (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_cmb_shift_param_parent_class)->finalize (object);
}

static void _nc_data_cmb_shift_param_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cmb_shift_param_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_cmb_shift_param_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_cmb_shift_param_class_init (NcDataCMBShiftParamClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass *diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);
  
  object_class->set_property = &nc_data_cmb_shift_param_set_property;
  object_class->get_property = &nc_data_cmb_shift_param_get_property;
  object_class->dispose      = &nc_data_cmb_shift_param_dispose;
  object_class->finalize     = &nc_data_cmb_shift_param_finalize;
  
  /**
   * NcDataCMBShiftParam:dist:
   *
   * The #NcDistance object to be used.
   */
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcDataCMBShiftParam:z:
   *
   * The data redshift.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshift",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->prepare   = &_nc_data_cmb_shift_param_prepare;
  diag_class->mean_func = &_nc_data_cmb_shift_param_mean_func;
  diag_class->set_size  = &_nc_data_cmb_shift_param_set_size;
}

static void
_nc_data_cmb_shift_param_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (data);
  NcHICosmo *cosmo                     = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  
  nc_distance_prepare_if_needed (cmb_shift_param->dist, cosmo);
}

/***************************************************************************
 * Dilation scale ratio Dv(0.35)/Dv(0.2) = 1.858 +/- 0.051 (arXiv:0705.3323)
 *
 ****************************************************************************/

static void
_nc_data_cmb_shift_param_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (diag);
  NcHICosmo *cosmo                     = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble zrec                   = ncm_vector_get (cmb_shift_param->x, 0);
  const gdouble R                      = nc_distance_shift_parameter (cmb_shift_param->dist, cosmo, zrec);
  
  ncm_vector_set (vp, 0, R);
}

static void
_nc_data_cmb_shift_param_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (diag);
  
  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&cmb_shift_param->x);
  
  if ((np != 0) && (np != diag->np))
    cmb_shift_param->x = ncm_vector_new (np);
  
  /* Chain up : start */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_cmb_shift_param_parent_class)->set_size (diag, np);
}

/**
 * nc_data_cmb_shift_param_new_empty:
 * @dist: a #NcDistance
 *
 * This function allocates an empty #NcDataCMBShiftParam object. 
 *
 * Returns: a #NcDataCMBShiftParam.
 */
NcDataCMBShiftParam *
nc_data_cmb_shift_param_new_empty (NcDistance *dist)
{
  NcDataCMBShiftParam *cmb_shift_param = g_object_new (NC_TYPE_DATA_CMB_SHIFT_PARAM,
                                                       "dist", dist,
                                                       NULL);
  
  return cmb_shift_param;
}

/**
 * nc_data_cmb_shift_param_new_from_file:
 * @filename: file containing a serialized #NcDataCMBShiftParam.
 *
 * Creates a new #NcDataCMBShiftParam from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataCMBShiftParam.
 */
NcDataCMBShiftParam *
nc_data_cmb_shift_param_new_from_file (const gchar *filename)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (ncm_serialize_global_from_file (filename));
  
  g_assert (NC_IS_DATA_CMB_SHIFT_PARAM (cmb_shift_param));
  
  return cmb_shift_param;
}

/**
 * nc_data_cmb_shift_param_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataCMBId
 *
 * Creates a new #NcDataCMBShiftParam from @id.
 * See #NcDataCMBId for the available options.
 *
 * Returns: the newly created #NcDataCMBShiftParam.
 */
NcDataCMBShiftParam *
nc_data_cmb_shift_param_new_from_id (NcDistance *dist, NcDataCMBId id)
{
  NcDataCMBShiftParam *cmb_shift_param;
  gchar *filename;
  
  switch (id)
  {
    case NC_DATA_CMB_SHIFT_PARAM_WMAP3:
      filename = ncm_cfg_get_data_filename ("nc_data_cmb_wmap3_shift_param.obj", TRUE);
      break;
    case NC_DATA_CMB_SHIFT_PARAM_WMAP5:
      filename = ncm_cfg_get_data_filename ("nc_data_cmb_wmap5_shift_param.obj", TRUE);
      break;
    case NC_DATA_CMB_SHIFT_PARAM_WMAP7:
      filename = ncm_cfg_get_data_filename ("nc_data_cmb_wmap7_shift_param.obj", TRUE);
      break;
    default:
      g_error ("nc_data_cmb_shift_param_new_from_id: id %d not recognized.", id);
      break;
  }
  
  cmb_shift_param = nc_data_cmb_shift_param_new_from_file (filename);
  nc_data_cmb_shift_param_set_dist (cmb_shift_param, dist);
  g_free (filename);
  
  return cmb_shift_param;
}

/**
 * nc_data_cmb_shift_param_set_dist:
 * @cmb_shift_param: a #NcDataCMBShiftParam
 * @dist: a #NcDistance
 *
 * Sets the distance object @dist to @cmb_shift_param.
 *
 */
void
nc_data_cmb_shift_param_set_dist (NcDataCMBShiftParam *cmb_shift_param, NcDistance *dist)
{
  nc_distance_clear (&cmb_shift_param->dist);
  cmb_shift_param->dist = nc_distance_ref (dist);
}

