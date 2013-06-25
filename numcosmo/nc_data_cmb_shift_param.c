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
 * @title: Cosmic Microwave Background Data -- Shift Parameter
 * @short_description: CMB shift parameter implementation
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_cmb_shift_param.h"
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_ID,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataCMBShiftParam, nc_data_cmb_shift_param, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_cmb_shift_param_init (NcDataCMBShiftParam *cmb_shift_param)
{
  cmb_shift_param->dist = NULL;
  cmb_shift_param->id   = NC_DATA_CMB_NSAMPLES;
}

static void
nc_data_cmb_shift_param_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (object);
  g_return_if_fail (NC_IS_DATA_CMB_SHIFT_PARAM (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&cmb_shift_param->dist);
      cmb_shift_param->dist = g_value_dup_object (value);
      break;
    case PROP_ID:
      nc_data_cmb_shift_param_set_sample (cmb_shift_param, g_value_get_enum (value));
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
    case PROP_ID:
      g_value_set_enum (value, nc_data_cmb_shift_param_get_sample (cmb_shift_param));
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

static void
nc_data_cmb_shift_param_class_init (NcDataCMBShiftParamClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->set_property = &nc_data_cmb_shift_param_set_property;
  object_class->get_property = &nc_data_cmb_shift_param_get_property;
  object_class->dispose      = &nc_data_cmb_shift_param_dispose;
  object_class->finalize     = &nc_data_cmb_shift_param_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ID,
                                   g_param_spec_enum ("sample-id",
                                                      NULL,
                                                      "Sample id",
                                                      NC_TYPE_DATA_CMB_ID, NC_DATA_CMB_SHIFT_PARAM_WMAP7,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  

  data_class->prepare   = &_nc_data_cmb_shift_param_prepare;
  diag_class->mean_func = &_nc_data_cmb_shift_param_mean_func;
}

static void
_nc_data_cmb_shift_param_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataCMBShiftParam *cmb_shift_param = NC_DATA_CMB_SHIFT_PARAM (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
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
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble zrec = ncm_vector_get (cmb_shift_param->x, 0);
  const gdouble R    = nc_distance_shift_parameter (cmb_shift_param->dist, cosmo, zrec);

  ncm_vector_set (vp, 0, R);
}

/**
 * nc_data_cmb_shift_param_new:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cmb_shift_param_new (NcDistance *dist, NcDataCMBId id)
{
  return g_object_new (NC_TYPE_DATA_CMB_SHIFT_PARAM,
                       "sample-id", id,
                       "dist", dist,
                       NULL);
}

/**
 * nc_data_cmb_shift_param_set_sample:
 * @cmb_shift_param: a #NcDataCMBShiftParam.
 * @id: FIXME
 *
 * FIXME
 *
 */
void 
nc_data_cmb_shift_param_set_sample (NcDataCMBShiftParam *cmb_shift_param, NcDataCMBId id)
{
  NcmData *data = NCM_DATA (cmb_shift_param);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (cmb_shift_param);
  
  g_assert (id < NC_DATA_CMB_NSAMPLES);

  if (cmb_shift_param->x != NULL && ncm_vector_len (cmb_shift_param->x) != 1)
    ncm_vector_clear (&cmb_shift_param->x);

  if (cmb_shift_param->x == NULL)
    cmb_shift_param->x = ncm_vector_new (1);
  
  if (data->desc != NULL)
    g_free (data->desc);

  ncm_data_gauss_diag_set_size (diag, 1);
  cmb_shift_param->id = id;

  switch (id)
  {
    case NC_DATA_CMB_SHIFT_PARAM_WMAP3:
      data->desc = g_strdup ("WMAP3 shift parameter");
      ncm_vector_set (cmb_shift_param->x, 0, ncm_c_wmap3_cmb_z ());
      ncm_vector_set (diag->y,            0, ncm_c_wmap3_cmb_R ());
      ncm_vector_set (diag->sigma,        0, ncm_c_wmap3_cmb_sigma_R ());
      break;
    case NC_DATA_CMB_SHIFT_PARAM_WMAP5:
      data->desc = g_strdup ("WMAP5 shift parameter");
      ncm_vector_set (cmb_shift_param->x, 0, ncm_c_wmap5_cmb_z ());
      ncm_vector_set (diag->y,            0, ncm_c_wmap5_cmb_R ());
      ncm_vector_set (diag->sigma,        0, ncm_c_wmap5_cmb_sigma_R ());
      break;
    case NC_DATA_CMB_SHIFT_PARAM_WMAP7:
      data->desc = g_strdup ("WMAP7 shift parameter");
      ncm_vector_set (cmb_shift_param->x, 0, ncm_c_wmap7_cmb_z ());
      ncm_vector_set (diag->y,            0, ncm_c_wmap7_cmb_R ());
      ncm_vector_set (diag->sigma,        0, ncm_c_wmap7_cmb_sigma_R ());
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  ncm_data_set_init (data);
}

/**
 * nc_data_cmb_shift_param_get_sample:
 * @cmb_shift_param: a #NcDataCMBShiftParam
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcDataCMBId 
nc_data_cmb_shift_param_get_sample (NcDataCMBShiftParam *cmb_shift_param)
{
  return cmb_shift_param->id;
}
