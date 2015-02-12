/***************************************************************************
 *            nc_data_bao_dv.c
 *
 *  Thu Apr 22 15:31:32 2010
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
 * SECTION:nc_data_bao_dv
 * @title: Baryonic Oscillation Data -- Volume Mean
 * @short_description: BAO averaged volume $D_V$ estimator
 *
 * See <link linkend="XEisenstein2005">Eisenstein et al. (2005)</link>.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_dv.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDV, nc_data_bao_dv, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_bao_dv_init (NcDataBaoDV *bao_dv)
{
  bao_dv->dist = NULL;
  bao_dv->x    = NULL;
}

static void
nc_data_bao_dv_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (object);
  g_return_if_fail (NC_IS_DATA_BAO_DV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&bao_dv->dist);
      bao_dv->dist = g_value_dup_object (value);
      break;
    case PROP_Z:
      ncm_vector_set_from_variant (bao_dv->x, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dv_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (object);
  g_return_if_fail (NC_IS_DATA_BAO_DV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, bao_dv->dist);
      break;
    case PROP_Z:
      g_value_take_variant (value, ncm_vector_get_variant (bao_dv->x));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dv_dispose (GObject *object)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (object);

  nc_distance_clear (&bao_dv->dist);
  ncm_vector_clear (&bao_dv->x);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dv_parent_class)->dispose (object);
}


static void
nc_data_bao_dv_finalize (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_bao_dv_parent_class)->finalize (object);
}

static void _nc_data_bao_dv_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_dv_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_dv_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_bao_dv_class_init (NcDataBaoDVClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->set_property = &nc_data_bao_dv_set_property;
  object_class->get_property = &nc_data_bao_dv_get_property;
  object_class->dispose      = &nc_data_bao_dv_dispose;
  object_class->finalize     = &nc_data_bao_dv_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_variant ("z",
                                                         NULL,
                                                         "Data redshift",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare   = &_nc_data_bao_dv_prepare;
  diag_class->mean_func = &_nc_data_bao_dv_mean_func;
  diag_class->set_size  = &_nc_data_bao_dv_set_size;
}

static void
_nc_data_bao_dv_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (bao_dv->dist, cosmo);
}

static void 
_nc_data_bao_dv_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  guint i;

  for (i = 0; i < diag->np; i++)
  {
    const gdouble z = ncm_vector_get (bao_dv->x, i);
    const gdouble DV = nc_distance_dilation_scale (bao_dv->dist, cosmo, z);
    ncm_vector_set (vp, i, DV);
  }
}

static void
_nc_data_bao_dv_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (diag);

  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&bao_dv->x);

  if ((np != 0) && (np != diag->np))
    bao_dv->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_bao_dv_parent_class)->set_size (diag, np);
}

/**
 * nc_data_bao_dv_new:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: a #NcmData
 */
NcmData *
nc_data_bao_dv_new (NcDistance *dist, NcDataBaoId id)
{
  NcmData *data = g_object_new (NC_TYPE_DATA_BAO_DV,
                                "dist", dist,
                                NULL);
  nc_data_bao_dv_set_sample (NC_DATA_BAO_DV (data), id);
  return data;
}

/**
 * nc_data_bao_dv_set_sample:
 * @bao_dv: a #NcDataBaoDV.
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 */
void 
nc_data_bao_dv_set_sample (NcDataBaoDV *bao_dv, NcDataBaoId id)
{
  NcmData *data = NCM_DATA (bao_dv);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (bao_dv);
  
  g_assert (id == NC_DATA_BAO_DV_EISENSTEIN2005);

  ncm_data_set_desc (data, "Eisenstein 2005, BAO Sample Dv");

  ncm_data_gauss_diag_set_size (diag, 1);

  ncm_vector_set (bao_dv->x,   0, ncm_c_bao_eisenstein_z ());
  ncm_vector_set (diag->y,     0, ncm_c_bao_eisenstein_DV ());
  ncm_vector_set (diag->sigma, 0, ncm_c_bao_eisenstein_sigma_DV ());

  ncm_data_set_init (data, TRUE);
}
