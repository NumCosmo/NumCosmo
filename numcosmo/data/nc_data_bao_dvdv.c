/***************************************************************************
 *            nc_data_bao_dvdv.c
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
 * SECTION:nc_data_bao_dvdv
 * @title: NcDataBaoDVDV
 * @short_description: Baryon oscillation data -- $D_V / D_V$ ratio.
 *
 * See [Percival et al. (2007)][XPercival2007].
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_dvdv.h"

#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDVDV, nc_data_bao_dvdv, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_bao_dvdv_init (NcDataBaoDVDV *bao_dvdv)
{
  bao_dvdv->dist = nc_distance_new (2.0);
}

static void
nc_data_bao_dvdv_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoDVDV *bao_dvdv = NC_DATA_BAO_DVDV (object);
  g_return_if_fail (NC_IS_DATA_BAO_DVDV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_dvdv_set_dist (bao_dvdv, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dvdv_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoDVDV *bao_dvdv = NC_DATA_BAO_DVDV (object);
  g_return_if_fail (NC_IS_DATA_BAO_DVDV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, bao_dvdv->dist);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dvdv_dispose (GObject *object)
{
  NcDataBaoDVDV *bao_dvdv = NC_DATA_BAO_DVDV (object);

  nc_distance_clear (&bao_dvdv->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dvdv_parent_class)->dispose (object);
}


static void
nc_data_bao_dvdv_finalize (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_bao_dvdv_parent_class)->finalize (object);
}

static void _nc_data_bao_dvdv_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_dvdv_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);

static void
nc_data_bao_dvdv_class_init (NcDataBaoDVDVClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->set_property = &nc_data_bao_dvdv_set_property;
  object_class->get_property = &nc_data_bao_dvdv_get_property;
  object_class->dispose      = &nc_data_bao_dvdv_dispose;
  object_class->finalize     = &nc_data_bao_dvdv_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare   = &_nc_data_bao_dvdv_prepare;
  diag_class->mean_func = &_nc_data_bao_dvdv_mean_func;
}

static void
_nc_data_bao_dvdv_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDVDV *bao_dvdv = NC_DATA_BAO_DVDV (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (bao_dvdv->dist, cosmo);
}

/***************************************************************************
 * Dilation scale ratio Dv(0.35)/Dv(0.2) = 1.858 +/- 0.051 (arXiv:0705.3323)
 *
 ****************************************************************************/

static void
_nc_data_bao_dvdv_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDVDV *bao_dvdv = NC_DATA_BAO_DVDV (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  const gdouble Dv_035 = nc_distance_dilation_scale (bao_dvdv->dist, cosmo, 0.35);
  const gdouble Dv_020 = nc_distance_dilation_scale (bao_dvdv->dist, cosmo, 0.20);

  ncm_vector_set (vp, 0, Dv_035 / Dv_020);
}

/**
 * nc_data_bao_dvdv_new_from_file:
 * @filename: file containing a serialized #NcDataBaoDVDV.
 *
 * Creates a new #NcDataBaoDVDV from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoDVDV.
 */
NcDataBaoDVDV *
nc_data_bao_dvdv_new_from_file (const gchar *filename)
{
  NcDataBaoDVDV *bao_dvdv = NC_DATA_BAO_DVDV (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_BAO_DVDV (bao_dvdv));

  return bao_dvdv;
}

/**
 * nc_data_bao_dvdv_new:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: a #NcDataBaoDVDV
 */
NcDataBaoDVDV *
nc_data_bao_dvdv_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoDVDV *bao_dvdv;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_DVDV_PERCIVAL2007:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dvdv_percival2007.obj", TRUE);
      break;
    case NC_DATA_BAO_DVDV_PERCIVAL2010:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dvdv_percival2010.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_dvdv_new_from_id: id %d not recognized.", id);
      break;
  }

  bao_dvdv = nc_data_bao_dvdv_new_from_file (filename);
  nc_data_bao_dvdv_set_dist (bao_dvdv, dist);
  g_free (filename);

  return bao_dvdv;
}

/**
 * nc_data_bao_dvdv_set_dist:
 * @bao_dvdv: a #NcDataBaoDVDV
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_dvdv_set_dist (NcDataBaoDVDV *bao_dvdv, NcDistance *dist)
{
  nc_distance_clear (&bao_dvdv->dist);
  bao_dvdv->dist = nc_distance_ref (dist);
}
