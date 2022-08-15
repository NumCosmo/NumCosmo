/***************************************************************************
 *            nc_data_bao_dtr_dhr.c
 *
 *  Mon Ago 15 14:38:48 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_bao_dtr_dhr.c
 * Copyright (C) 2022 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_data_bao_dtr_dhr
 * @title: NcDataBaoDtrDHr
 * @short_description: Baryon Oscillation Data -- $(D_H/r,\; D_t/r)$ data.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_dtr_dhr.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDtrDHr, nc_data_bao_dtr_dhr, NCM_TYPE_DATA_GAUSS_COV);

static void
nc_data_bao_dtr_dhr_init (NcDataBaoDtrDHr *dhdt)
{
  dhdt->dist = nc_distance_new (3.0);
  dhdt->x    = NULL;
}

static void
nc_data_bao_dtr_dhr_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (object);
  g_return_if_fail (NC_IS_DATA_BAO_DTR_DHR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_dtr_dhr_set_dist (dhdt, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&dhdt->x, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dtr_dhr_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (object);
  g_return_if_fail (NC_IS_DATA_BAO_DTR_DHR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dhdt->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, dhdt->x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dtr_dhr_dispose (GObject *object)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (object);

  nc_distance_clear (&dhdt->dist);
  ncm_vector_clear (&dhdt->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dtr_dhr_parent_class)->dispose (object);
}

static void
nc_data_bao_dtr_dhr_finalize (GObject *object)
{
  /*NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dtr_dhr_parent_class)->finalize (object);
}

static void _nc_data_bao_dtr_dhr_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_dtr_dhr_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_dtr_dhr_set_size (NcmDataGaussCov *gauss, guint np);

static void
nc_data_bao_dtr_dhr_class_init (NcDataBaoDtrDHrClass *klass)
{
  GObjectClass* object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_bao_dtr_dhr_set_property;
  object_class->get_property = &nc_data_bao_dtr_dhr_get_property;
  object_class->dispose      = &nc_data_bao_dtr_dhr_dispose;
  object_class->finalize     = &nc_data_bao_dtr_dhr_finalize;

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

  data_class->prepare    = &_nc_data_bao_dtr_dhr_prepare;
  gauss_class->mean_func = &_nc_data_bao_dtr_dhr_mean_func;
  gauss_class->set_size  = &_nc_data_bao_dtr_dhr_set_size;
}

static void
_nc_data_bao_dtr_dhr_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (data);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (dhdt->dist, cosmo);
}

static void
_nc_data_bao_dtr_dhr_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (gauss);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
    if (i % 2 == 0)
    {
      const gdouble z1      = ncm_vector_get (dhdt->x, i);
      const gdouble z2      = ncm_vector_get (dhdt->x, i + 1);
      const gdouble DH_r    = nc_distance_DH_r (dhdt->dist, cosmo, z1);
      const gdouble Dt_r    = nc_distance_DA_r (dhdt->dist, cosmo, z2);

      ncm_vector_set (vp, i, Dt_r);
      ncm_vector_set (vp, i + 1, DH_r);
    }
    else
      continue;
  }
}

static void
_nc_data_bao_dtr_dhr_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (gauss);

  if ((np == 0) || (np != gauss->np))
    ncm_vector_clear (&dhdt->x);

  if ((np != 0) && (np != gauss->np))
    dhdt->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_bao_dtr_dhr_parent_class)->set_size (gauss, np);
}

/**
 * nc_data_bao_dtr_dhr_new_from_file:
 * @filename: file containing a serialized #NcDataBaoDtrDHr.
 *
 * Creates a new #NcDataBaoDtrDHr from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoDtrDHr.
 */
NcDataBaoDtrDHr *
nc_data_bao_dtr_dhr_new_from_file (const gchar *filename)
{
  NcDataBaoDtrDHr *dhdt = NC_DATA_BAO_DTR_DHR (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_BAO_DTR_DHR (dhdt));
  return dhdt;
}

/**
 * nc_data_bao_dtr_dhr_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: (transfer full): a #NcDataBaoDtrDHr
 */
NcDataBaoDtrDHr *
nc_data_bao_dtr_dhr_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoDtrDHr *dhdt;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dtr_dhr_dr12_2016_dr16compatible.obj", TRUE);
      break;
    case NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dtr_dhr_dr16_lrg_2021.obj", TRUE);
      break;
    case NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dtr_dhr_dr16_qso_2021.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_dtr_dhr_new_from_id: id %d not recognized.", id);
      break;
  }
  dhdt = nc_data_bao_dtr_dhr_new_from_file (filename);
  nc_data_bao_dtr_dhr_set_dist (dhdt, dist);
  g_free (filename);

  return dhdt;
}

/**
 * nc_data_bao_dtr_dhr_set_dist:
 * @dhdt: a #NcDataBaoDtrDHr
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_dtr_dhr_set_dist (NcDataBaoDtrDHr *dhdt, NcDistance *dist)
{
  nc_distance_clear (&dhdt->dist);
  dhdt->dist = nc_distance_ref (dist);
}
