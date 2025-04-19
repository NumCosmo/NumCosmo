/***************************************************************************
 *            nc_data_bao_rdv.c
 *
 *  Thu Apr 22 15:31:32 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcDataBaoRDV:
 *
 * Baryon Oscillation Data -- $r_s / D_V$ ratio.
 *
 * See [Percival et al. (2007)][XPercival2007].
 *
 * Kazin et al. (arXiv:1401.0358): our implementation of the inverse covariance matrix
 * is given by $$C^{-1}_{new} = \frac{1}{r_s^{\text{fid}}} C^{-1},$$ where
 * $r_s^{\text{fid}} = 148.6$ and $C^{-1}$ is given in table 4. This modification is due
 * the fact that we are using $D_V(z)/r_s(z_d)$ instead of $D_V(z)*
 * r_s^{\text{fid}}/r_s(z_d)$. Analogously, we implemented $D_V(z) / r_s(z_d) = [1716.4,
 * 2220.8, 2516.1] / 148.6 = [11.550, 14.945, 16.932]$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_rdv.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

typedef struct _NcDataBaoRDVPrivate
{
  gint placeholder;
} NcDataBaoRDVPrivate;

struct _NcDataBaoRDV
{
  NcmDataGauss parent_instance;
  NcDistance *dist;
  NcmVector *x;
  gboolean r_DV;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_DATA_FORM,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoRDV, nc_data_bao_rdv, NCM_TYPE_DATA_GAUSS)

static void
nc_data_bao_rdv_init (NcDataBaoRDV *bao_rdv)
{
  bao_rdv->x    = NULL;
  bao_rdv->dist = nc_distance_new (2.0);
  bao_rdv->r_DV = FALSE;
}

static void
_nc_data_bao_rdv_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_bao_rdv_parent_class)->constructed (object);
}

static void
nc_data_bao_rdv_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (object);

  g_return_if_fail (NC_IS_DATA_BAO_RDV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_rdv_set_dist (bao_rdv, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&bao_rdv->x, g_value_get_object (value), TRUE);
      break;
    case PROP_DATA_FORM:
      bao_rdv->r_DV = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_bao_rdv_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (object);

  g_return_if_fail (NC_IS_DATA_BAO_RDV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, bao_rdv->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, bao_rdv->x);
      break;
    case PROP_DATA_FORM:
      g_value_set_boolean (value, bao_rdv->r_DV);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_bao_rdv_dispose (GObject *object)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (object);

  ncm_vector_clear (&bao_rdv->x);
  nc_distance_clear (&bao_rdv->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_rdv_parent_class)->dispose (object);
}

static void
nc_data_bao_rdv_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_rdv_parent_class)->finalize (object);
}

static void _nc_data_bao_rdv_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_rdv_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_rdv_set_size (NcmDataGauss *gauss, guint np);

static void
nc_data_bao_rdv_class_init (NcDataBaoRDVClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_CLASS (klass);

  object_class->constructed  = &_nc_data_bao_rdv_constructed;
  object_class->set_property = &nc_data_bao_rdv_set_property;
  object_class->get_property = &nc_data_bao_rdv_get_property;
  object_class->dispose      = &nc_data_bao_rdv_dispose;
  object_class->finalize     = &nc_data_bao_rdv_finalize;

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
  g_object_class_install_property (object_class,
                                   PROP_DATA_FORM,
                                   g_param_spec_boolean ("is-rDV",
                                                         NULL,
                                                         "Whether the format is r/DV or DV/r",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  data_class->prepare    = &_nc_data_bao_rdv_prepare;
  gauss_class->mean_func = &_nc_data_bao_rdv_mean_func;
  gauss_class->set_size  = &_nc_data_bao_rdv_set_size;
}

static void
_nc_data_bao_rdv_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (data);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  nc_distance_prepare_if_needed (bao_rdv->dist, cosmo);
}

static void
_nc_data_bao_rdv_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (gauss);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const guint np        = ncm_data_gauss_get_size (gauss);
  guint i;

  if (bao_rdv->r_DV)
  {
    for (i = 0; i < np; i++)
    {
      const gdouble z    = ncm_vector_get (bao_rdv->x, i);
      const gdouble r_Dv = nc_distance_bao_r_Dv (bao_rdv->dist, cosmo, z);

      ncm_vector_set (vp, i, r_Dv);
    }
  }
  else
  {
    for (i = 0; i < np; i++)
    {
      const gdouble z    = ncm_vector_get (bao_rdv->x, i);
      const gdouble Dv_r = 1.0 / nc_distance_bao_r_Dv (bao_rdv->dist, cosmo, z);

      ncm_vector_set (vp, i, Dv_r);
    }
  }
}

static void
_nc_data_bao_rdv_set_size (NcmDataGauss *gauss, guint np)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (gauss);
  const guint cnp       = ncm_data_gauss_get_size (gauss);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&bao_rdv->x);

  if ((np != 0) && (np != cnp))
    bao_rdv->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_CLASS (nc_data_bao_rdv_parent_class)->set_size (gauss, np);
}

/**
 * nc_data_bao_rdv_new_from_file:
 * @filename: file containing a serialized #NcDataBaoRDV.
 *
 * Creates a new #NcDataBaoRDV from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoRDV.
 */
NcDataBaoRDV *
nc_data_bao_rdv_new_from_file (const gchar *filename)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_BAO_RDV (bao_rdv));

  return bao_rdv;
}

/**
 * nc_data_bao_rdv_new:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * Creates a new #NcDataBaoRDV from @id.
 * This object requires a #NcDistance object to be set.
 *
 * Returns: a #NcDataBaoRDV
 */
NcDataBaoRDV *
nc_data_bao_rdv_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoRDV *bao_rdv;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_RDV_PERCIVAL2007:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_percival2007.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_PERCIVAL2010:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_percival2010.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_BEUTLER2011:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_beutler2011.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_PADMANABHAN2012:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_padmanabhan2012.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_ANDERSON2012:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_anderson2012.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_BLAKE2012:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_blake2012.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_KAZIN2014:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_kazin2014.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_BOSS_QSO_ATA2017:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_boss_qso_ata2017.obj", TRUE);
      break;
    case NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_rdv_desi_dr1_bgs_qso_2024.obj", TRUE);
      break;
    default:                                                              /* LCOV_EXCL_LINE */
      g_error ("nc_data_bao_rdv_new_from_id: id %d not recognized.", id); /* LCOV_EXCL_LINE */
      break;                                                              /* LCOV_EXCL_LINE */
  }

  bao_rdv = nc_data_bao_rdv_new_from_file (filename);
  nc_data_bao_rdv_set_dist (bao_rdv, dist);
  g_free (filename);

  return bao_rdv;
}

/**
 * nc_data_bao_rdv_set_dist:
 * @bao_rdv: a #NcDataBaoRDV
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_rdv_set_dist (NcDataBaoRDV *bao_rdv, NcDistance *dist)
{
  nc_distance_clear (&bao_rdv->dist);
  bao_rdv->dist = nc_distance_ref (dist);
}

