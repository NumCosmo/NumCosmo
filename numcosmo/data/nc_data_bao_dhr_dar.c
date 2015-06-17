/***************************************************************************
 *            nc_data_bao_dhr_dar.c
 *
 *  Sat May 23 14:38:48 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_bao_dhr_dar.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_data_bao_dhr_dar
 * @title: NcDataBaoDHrDAr
 * @short_description: Baryon Oscillation Data -- $(D_H/r,\; D_A/r)$ data.
 *
 * See: <link linkend="XFont-Ribera2014">Font-Ribera et al. (2014)</link>, 
 * <link linkend="XDelubac2015">Delubac et al. (2015)</link> and
 * <link linkend="XAubourg2014">Aubourg et al. (2014)</link>.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_dhr_dar.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDHrDAr, nc_data_bao_dhr_dar, NCM_TYPE_DATA_GAUSS_COV);

static void
nc_data_bao_dhr_dar_init (NcDataBaoDHrDAr *dhda)
{
  dhda->dist = nc_distance_new (2.0);
  dhda->x    = NULL;
}

static void
nc_data_bao_dhr_dar_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (object);
  g_return_if_fail (NC_IS_DATA_BAO_DHR_DAR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_dhr_dar_set_dist (dhda, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_set_from_variant (dhda->x, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dhr_dar_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (object);
  g_return_if_fail (NC_IS_DATA_BAO_DHR_DAR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dhda->dist);
      break;
    case PROP_Z:
      g_value_take_variant (value, ncm_vector_get_variant (dhda->x));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dhr_dar_dispose (GObject *object)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (object);

  nc_distance_clear (&dhda->dist);
  ncm_vector_clear (&dhda->x);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dhr_dar_parent_class)->dispose (object);
}

static void
nc_data_bao_dhr_dar_finalize (GObject *object)
{
  /*NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dhr_dar_parent_class)->finalize (object);
}

static void _nc_data_bao_dhr_dar_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_dhr_dar_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_dhr_dar_set_size (NcmDataGaussCov *gauss, guint np);

static void
nc_data_bao_dhr_dar_class_init (NcDataBaoDHrDArClass *klass)
{
  GObjectClass* object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_bao_dhr_dar_set_property;
  object_class->get_property = &nc_data_bao_dhr_dar_get_property;
  object_class->dispose      = &nc_data_bao_dhr_dar_dispose;
  object_class->finalize     = &nc_data_bao_dhr_dar_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_variant ("z",
                                                         NULL,
                                                         "Data redshift",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_bao_dhr_dar_prepare;
  gauss_class->mean_func = &_nc_data_bao_dhr_dar_mean_func;
  gauss_class->set_size  = &_nc_data_bao_dhr_dar_set_size;
}

static void
_nc_data_bao_dhr_dar_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (data);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (dhda->dist, cosmo);
}

static void 
_nc_data_bao_dhr_dar_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (gauss);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  
  const gdouble z1      = ncm_vector_get (dhda->x, 0);
  const gdouble z2      = ncm_vector_get (dhda->x, 1);
  const gdouble DH_r    = nc_distance_DH_r (dhda->dist, cosmo, z1);
  const gdouble DA_r    = nc_distance_DA_r (dhda->dist, cosmo, z2);

  ncm_vector_set (vp, 0, DH_r);
  ncm_vector_set (vp, 1, DA_r);
}

static void 
_nc_data_bao_dhr_dar_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (gauss);

  if ((np == 0) || (np != gauss->np))
    ncm_vector_clear (&dhda->x);

  if ((np != 0) && (np != gauss->np))
    dhda->x = ncm_vector_new (np);
  
  /* Chain up : end */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_bao_dhr_dar_parent_class)->set_size (gauss, np);
}

/**
 * nc_data_bao_dhr_dar_new_from_file:
 * @filename: file containing a serialized #NcDataBaoDHrDAr.
 * 
 * Creates a new #NcDataBaoDHrDAr from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataBaoDHrDAr.
 */
NcDataBaoDHrDAr *
nc_data_bao_dhr_dar_new_from_file (const gchar *filename)
{
  NcDataBaoDHrDAr *dhda = NC_DATA_BAO_DHR_DAR (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_BAO_DHR_DAR (dhda));
  return dhda;
}

/**
 * nc_data_bao_dhr_dar_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: (transfer full): a #NcDataBaoDHrDAr
 */
NcDataBaoDHrDAr *
nc_data_bao_dhr_dar_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoDHrDAr *dhda;
  gchar *filename;
  
  switch (id)
  {
    case NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dhr_dar_sdss_dr11_2015.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_dhr_dar_new_from_id: id %d not recognized.", id);
      break;
  }
  dhda = nc_data_bao_dhr_dar_new_from_file (filename);
  nc_data_bao_dhr_dar_set_dist (dhda, dist);
  g_free (filename);

  return dhda;
}

/**
 * nc_data_bao_dhr_dar_set_dist:
 * @dhda: a #NcDataBaoDHrDAr
 * @dist: a #NcDistance
 * 
 * Sets the distance object.
 * 
 */
void 
nc_data_bao_dhr_dar_set_dist (NcDataBaoDHrDAr *dhda, NcDistance *dist)
{
  nc_distance_clear (&dhda->dist);
  dhda->dist = nc_distance_ref (dist);
}
