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
 * @title: Baryonic Oscillation Data
 * @short_description: BAO Data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_bao_dv.h"

#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_ID,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDV, nc_data_bao_dv, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_bao_dv_init (NcDataBaoDV *bao_dv)
{
  bao_dv->dist = NULL;
  bao_dv->id   = NC_DATA_BAO_NSAMPLES;
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
    case PROP_ID:
      nc_data_bao_dv_set_sample (bao_dv, g_value_get_enum (value));
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
    case PROP_ID:
      g_value_set_enum (value, nc_data_bao_dv_get_sample (bao_dv));
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
                                   PROP_ID,
                                   g_param_spec_enum ("sample-id",
                                                      NULL,
                                                      "Sample id",
                                                      NC_TYPE_DATA_BAO_ID, NC_DATA_BAO_DV_EISENSTEIN,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  

  data_class->prepare   = &_nc_data_bao_dv_prepare;
  diag_class->mean_func = &_nc_data_bao_dv_mean_func;
}

static void
_nc_data_bao_dv_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  nc_distance_prepare_if_needed (bao_dv->dist, cosmo);
}

static void 
_nc_data_bao_dv_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDV *bao_dv = NC_DATA_BAO_DV (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  gint i;

  for (i = 0; i < diag->np; i++)
  {
    const gdouble z = ncm_vector_get (bao_dv->x, i);
    const gdouble DV = nc_distance_dilation_scale (bao_dv->dist, cosmo, z);
    ncm_vector_set (vp, i, DV);
  }
}

/**
 * nc_data_bao_dv_new:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_bao_dv_new (NcDistance *dist, NcDataBaoId id)
{
  return g_object_new (NC_TYPE_DATA_BAO_DV,
                       "sample-id", id,
                       "dist", dist,
                       NULL);
}

/**
 * nc_data_bao_dv_set_size:
 * @bao_dv: a #NcDataBaoDV
 * @np: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void 
nc_data_bao_dv_set_size (NcDataBaoDV *bao_dv, guint np)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (bao_dv);

  if (diag->np != 0)
    g_assert (bao_dv->x != NULL && ncm_vector_len (bao_dv->x) == diag->np);
  
  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&bao_dv->x);

  if ((np != 0) && (np != diag->np))
    bao_dv->x = ncm_vector_new_sunk (np);

  ncm_data_gauss_diag_set_size (NCM_DATA_GAUSS_DIAG (bao_dv), np);
}

/**
 * nc_data_bao_dv_get_size:
 * @bao_dv: a #NcDataBaoDV
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint 
nc_data_bao_dv_get_size (NcDataBaoDV *bao_dv)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (bao_dv);

  if (diag->np != 0)
    g_assert (bao_dv->x != NULL && ncm_vector_len (bao_dv->x) == diag->np);

  return ncm_data_gauss_diag_get_size (NCM_DATA_GAUSS_DIAG (bao_dv));
}

/**
 * nc_data_bao_dv_set_sample:
 * @bao_dv: a #NcDataBaoDV.
 * @id: FIXME
 *
 * FIXME
 *
 */
void 
nc_data_bao_dv_set_sample (NcDataBaoDV *bao_dv, NcDataBaoId id)
{
  NcmData *data = NCM_DATA (bao_dv);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (bao_dv);
  
  g_assert (id == NC_DATA_BAO_DV_EISENSTEIN);

  if (data->desc != NULL)
    g_free (data->desc);
  data->desc = g_strdup ("Eisenstein 2005, BAO Sample Dv");

  nc_data_bao_dv_set_size (bao_dv, 1);
  bao_dv->id = id;

  ncm_vector_set (bao_dv->x,   0, ncm_c_bao_eisenstein_z ());
  ncm_vector_set (diag->y,     0, ncm_c_bao_eisenstein_DV ());
  ncm_vector_set (diag->sigma, 0, ncm_c_bao_eisenstein_sigma_DV ());

  ncm_data_set_init (data);
}

/**
 * nc_data_bao_dv_get_sample:
 * @bao_dv: a #NcDataBaoDV
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcDataBaoId 
nc_data_bao_dv_get_sample (NcDataBaoDV *bao_dv)
{
  return bao_dv->id;
}
