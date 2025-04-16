/***************************************************************************
 *            nc_data_bao_dvr_dtdh.c
 *
 *  Tue Apr 15 17:27:40 2025
 *  Copyright  2025  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * nc_data_bao_dtr_dhr.c
 * Copyright (C) 2025 Mariana Penna-Lima <pennalima@unb.br>
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
 * NcDataBaoDvrDtdh:
 *
 * Baryon Oscillation Data -- $(D_V/r,\; D_t/D_H)$ data.
 *
 * The data is stored in a #NcDataBaoDvrDtDh object. The data is stored in a
 * #NcmDataGaussCov base class object, which is a subclass of #NcmData. The data
 * represents the mean values of the distance dilation scale $D_V$ divided by the
 * sound horizon at the last scattering surface $r_s$, and the mean values of the
 * transverse distance $D_t$ over the Hubble distance $D_H$ at the redshift $z$.
 *
 * DESI DR1 data - LRG and ELG 2024.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "data/nc_data_bao_dvr_dtdh.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

typedef struct _NcDataBaoDvrDtDhPrivate
{
  gint placeholder;
} NcDataBaoDvrDtDhPrivate;

struct _NcDataBaoDvrDtDh
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcDistance *dist;
  NcmVector *x;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDvrDtDh, nc_data_bao_dvr_dtdh, NCM_TYPE_DATA_GAUSS_COV)

static void
nc_data_bao_dvr_dtdh_init (NcDataBaoDvrDtDh *dvdtdh)
{
  dvdtdh->dist = nc_distance_new (3.0);
  dvdtdh->x    = NULL;
}

static void
nc_data_bao_dvr_dtdh_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (object);

  g_return_if_fail (NC_IS_DATA_BAO_DVR_DTDH (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_dvr_dtdh_set_dist (dvdtdh, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&dvdtdh->x, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dvr_dtdh_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (object);

  g_return_if_fail (NC_IS_DATA_BAO_DVR_DTDH (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dvdtdh->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, dvdtdh->x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dvr_dtdh_dispose (GObject *object)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (object);

  nc_distance_clear (&dvdtdh->dist);
  ncm_vector_clear (&dvdtdh->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dvr_dtdh_parent_class)->dispose (object);
}

static void
nc_data_bao_dvr_dtdh_finalize (GObject *object)
{
  /*NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dvr_dtdh_parent_class)->finalize (object);
}

static void _nc_data_bao_dvr_dtdh_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_dvr_dtdh_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_dvr_dtdh_set_size (NcmDataGaussCov *gauss, guint np);

static void
nc_data_bao_dvr_dtdh_class_init (NcDataBaoDvrDtDhClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_bao_dvr_dtdh_set_property;
  object_class->get_property = &nc_data_bao_dvr_dtdh_get_property;
  object_class->dispose      = &nc_data_bao_dvr_dtdh_dispose;
  object_class->finalize     = &nc_data_bao_dvr_dtdh_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshift",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_bao_dvr_dtdh_prepare;
  gauss_class->mean_func = &_nc_data_bao_dvr_dtdh_mean_func;
  gauss_class->set_size  = &_nc_data_bao_dvr_dtdh_set_size;
}

static void
_nc_data_bao_dvr_dtdh_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (data);
  NcHICosmo *cosmo         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  nc_distance_prepare_if_needed (dvdtdh->dist, cosmo);
}

static void
_nc_data_bao_dvr_dtdh_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (gauss);
  NcHICosmo *cosmo         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const guint np           = ncm_data_gauss_cov_get_size (gauss);
  guint i;

  for (i = 0; i < np; i++)
  {
    if (i % 2 == 0)
    {
      const gdouble z1     = ncm_vector_get (dvdtdh->x, i);
      const gdouble z2     = ncm_vector_get (dvdtdh->x, i + 1);
      const gdouble DV_r   = 1.0 / nc_distance_bao_r_Dv (dvdtdh->dist, cosmo, z1);
      const gdouble one_DH = nc_hicosmo_E (cosmo, z2);
      const gdouble Dt_DH  = nc_distance_transverse (dvdtdh->dist, cosmo, z2) * one_DH;

      ncm_vector_set (vp, i, DV_r);
      ncm_vector_set (vp, i + 1, Dt_DH);
    }
    else
    {
      continue;
    }
  }
}

static void
_nc_data_bao_dvr_dtdh_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (gauss);
  const guint cnp          = ncm_data_gauss_cov_get_size (gauss);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&dvdtdh->x);

  if ((np != 0) && (np != cnp))
    dvdtdh->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_bao_dvr_dtdh_parent_class)->set_size (gauss, np);
}

/**
 * nc_data_bao_dvr_dtdh_new_from_file:
 * @filename: file containing a serialized #NcDataBaoDvrDtDh.
 *
 * Creates a new #NcDataBaoDvrDtDh from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoDvrDtDh.
 */
NcDataBaoDvrDtDh *
nc_data_bao_dvr_dtdh_new_from_file (const gchar *filename)
{
  NcDataBaoDvrDtDh *dvdtdh = NC_DATA_BAO_DVR_DTDH (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_BAO_DVR_DTDH (dvdtdh));

  return dvdtdh;
}

/**
 * nc_data_bao_dvr_dtdh_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * Creates a new acustic scale data object #NcDataBaoDvrDtDh from @id.
 * This object requires a #NcDistance object to be set.
 *
 *
 * Returns: (transfer full): a #NcDataBaoDvrDtDh
 */
NcDataBaoDvrDtDh *
nc_data_bao_dvr_dtdh_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoDvrDtDh *dvdtdh;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dvr_dtdh_desi_dr1_lrg_elg_2024.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_dtr_dhr_new_from_id: id %d not recognized.", id);
      break;
  }

  dvdtdh = nc_data_bao_dvr_dtdh_new_from_file (filename);
  nc_data_bao_dvr_dtdh_set_dist (dvdtdh, dist);
  g_free (filename);

  return dvdtdh;
}

/**
 * nc_data_bao_dvr_dtdh_set_dist:
 * @dvdtdh: a #NcDataBaoDvrDtDh
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_dvr_dtdh_set_dist (NcDataBaoDvrDtDh *dvdtdh, NcDistance *dist)
{
  nc_distance_clear (&dvdtdh->dist);
  dvdtdh->dist = nc_distance_ref (dist);
}

