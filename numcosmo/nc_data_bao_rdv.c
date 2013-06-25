/***************************************************************************
 *            nc_data_bao_rdv.c
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
 * SECTION:nc_data_bao_rdv
 * @title: Baryonic Oscillation Data -- rDv
 * @short_description: BAO $r/D_V$ ratio estimator
 *
 * See <link linkend="XPercival2007">Percival et al. (2007)</link>.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_bao_rdv.h"

#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_ID,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoRDV, nc_data_bao_rdv, NCM_TYPE_DATA_GAUSS);

static void
nc_data_bao_rdv_init (NcDataBaoRDV *bao_rdv)
{
  bao_rdv->x    = NULL;
  bao_rdv->dist = NULL;
  bao_rdv->id   = NC_DATA_BAO_NSAMPLES;
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
      nc_distance_clear (&bao_rdv->dist);
      bao_rdv->dist = g_value_dup_object (value);
      break;
    case PROP_ID:
      nc_data_bao_rdv_set_sample (bao_rdv, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    case PROP_ID:
      g_value_set_enum (value, nc_data_bao_rdv_get_sample (bao_rdv));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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

static void
nc_data_bao_rdv_class_init (NcDataBaoRDVClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ID,
                                   g_param_spec_enum ("sample-id",
                                                      NULL,
                                                      "Sample id",
                                                      NC_TYPE_DATA_BAO_ID, NC_DATA_BAO_RDV_PERCIVAL2010,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_bao_rdv_prepare;
  gauss_class->mean_func = &_nc_data_bao_rdv_mean_func;

}

static void
_nc_data_bao_rdv_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (bao_rdv->dist, cosmo);
}

static void 
_nc_data_bao_rdv_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (gauss);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
    const gdouble z  = ncm_vector_get (bao_rdv->x, i);
    const gdouble r_Dv = nc_distance_bao_r_Dv (bao_rdv->dist, cosmo, z);
    ncm_vector_set (vp, i, r_Dv);
  }
}

/**
 * nc_data_bao_rdv_new:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_bao_rdv_new (NcDistance *dist, NcDataBaoId id)
{
  return g_object_new (NC_TYPE_DATA_BAO_RDV,
                       "sample-id", id,
                       "dist", dist,
                       NULL);
}

/**
 * nc_data_bao_rdv_set_size:
 * @bao_rdv: a #NcDataBaoRDV
 * @np: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void 
nc_data_bao_rdv_set_size (NcDataBaoRDV *bao_rdv, guint np)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (bao_rdv);

  if (gauss->np != 0)
    g_assert (bao_rdv->x != NULL && ncm_vector_len (bao_rdv->x) == gauss->np);
  
  if ((np == 0) || (np != gauss->np))
    ncm_vector_clear (&bao_rdv->x);

  if ((np != 0) && (np != gauss->np))
    bao_rdv->x = ncm_vector_new (np);

  ncm_data_gauss_set_size (NCM_DATA_GAUSS (bao_rdv), np);
}

/**
 * nc_data_bao_rdv_get_size:
 * @bao_rdv: a #NcDataBaoRDV
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint 
nc_data_bao_rdv_get_size (NcDataBaoRDV *bao_rdv)
{
  NcmDataGauss *gauss = NCM_DATA_GAUSS (bao_rdv);

  if (gauss->np != 0)
    g_assert (bao_rdv->x != NULL && ncm_vector_len (bao_rdv->x) == gauss->np);

  return ncm_data_gauss_get_size (NCM_DATA_GAUSS (bao_rdv));
}


/***************************************************************************
 * BAO Percival priors data (arXiv:0705.3323)
 *
 ****************************************************************************/

static gdouble nc_bao_distance_priors_percival2007_z[]       = {    0.2,   0.35 };
static gdouble nc_bao_distance_priors_percival2007_bestfit[] = { 0.1980, 0.1094 };
static gdouble nc_bao_distance_priors_percival2007_inv_cov[][2] =
{ 
  { 35059.0, -24031.0},
  {-24031.0, 108300.0} 
};

static gdouble nc_bao_distance_priors_percival2010_z[]       = {    0.2,   0.35 };
static gdouble nc_bao_distance_priors_percival2010_bestfit[] = { 0.1905, 0.1097 };
static gdouble nc_bao_distance_priors_percival2010_inv_cov[][2] =
{ 
  { 30124.0, -17227.0},
  {-17227.0,  86977.0} 
};

/***************************************************************************
 * BAO 6dFGRS Beutler et al. (2011)
 ****************************************************************************/
static gdouble nc_bao_distance_priors_beutler2011_z[]       = {   0.1 };
static gdouble nc_bao_distance_priors_beutler2011_bestfit[] = { 0.336 };
static gdouble nc_bao_distance_priors_beutler2011_inv_cov[][1] =
{ 
  { 1.0 / (0.015 * 0.015) }, 
};

/***************************************************************************
 * BAO SDSS-DR7-rec Padmanabhan et al. (2012)
 ****************************************************************************/
static gdouble nc_bao_distance_priors_padmanabhan2012_z[]       = {  0.35 };
static gdouble nc_bao_distance_priors_padmanabhan2012_bestfit[] = { 1.0 / 8.88 };
static gdouble nc_bao_distance_priors_padmanabhan2012_inv_cov[][1] =
{ 
  { (8.88 * 8.88 * 8.88 * 8.88) / (0.17 * 0.17) },
};

/***************************************************************************
 * BAO SDSS-DR9-rec Anderson et al. (2012)
 ****************************************************************************/
static gdouble nc_bao_distance_priors_anderson2012_z[]       = {  0.57 };
static gdouble nc_bao_distance_priors_anderson2012_bestfit[] = { 1.0 / 13.67 };
static gdouble nc_bao_distance_priors_anderson2012_inv_cov[][1] =
{ 
  { (13.67 * 13.67 * 13.67 * 13.67) / (0.22 * 0.22) },
};

/***************************************************************************
 * BAO WiggleZ Blake et al. (2012)
 ****************************************************************************/
static gdouble nc_bao_distance_priors_blake2012_z[]       = {   0.44,   0.60,   0.73 };
static gdouble nc_bao_distance_priors_blake2012_bestfit[] = { 0.0916, 0.0726, 0.0592 };
static gdouble nc_bao_distance_priors_blake2012_inv_cov[][3] =
{ 
  {  24532.1, -25137.7,  12099.1 },
  { -25137.7, 134598.4, -64783.9 },
  {  12099.1, -64783.9, 128837.6 },
};

/**
 * nc_data_bao_rdv_set_sample:
 * @bao_rdv: a #NcDataBaoRDV.
 * @id: FIXME
 *
 * FIXME
 *
 */
void
nc_data_bao_rdv_set_sample (NcDataBaoRDV *bao_rdv, NcDataBaoId id)
{
  NcmData *data = NCM_DATA (bao_rdv);
  NcmDataGauss *gauss = NCM_DATA_GAUSS (bao_rdv);
  gint i, j;
  
  g_assert (id >= NC_DATA_BAO_RDV_PERCIVAL2007 && id <= NC_DATA_BAO_RDV_BLAKE2012);

  if (data->desc != NULL)
    g_free (data->desc);

  switch (id)
  {
    case NC_DATA_BAO_RDV_PERCIVAL2007:
    {
      data->desc = g_strdup ("Percival 2007, BAO Sample R-Dv");
      nc_data_bao_rdv_set_size (bao_rdv, 2);

      for (i = 0; i < 2; i++)
      {
        ncm_vector_set (bao_rdv->x, i, nc_bao_distance_priors_percival2007_z[i]);
        ncm_vector_set (gauss->y,   i, nc_bao_distance_priors_percival2007_bestfit[i]);
        for (j = 0; j < 2; j++)
          ncm_matrix_set (gauss->inv_cov, i, j, 
                          nc_bao_distance_priors_percival2007_inv_cov[i][j]);
      }
      break;
    }
    case NC_DATA_BAO_RDV_PERCIVAL2010:
    {
      data->desc = g_strdup ("Percival 2010, BAO Sample R-Dv");
      nc_data_bao_rdv_set_size (bao_rdv, 2);

      for (i = 0; i < 2; i++)
      {
        ncm_vector_set (bao_rdv->x, i, nc_bao_distance_priors_percival2010_z[i]);
        ncm_vector_set (gauss->y,   i, nc_bao_distance_priors_percival2010_bestfit[i]);
        for (j = 0; j < 2; j++)
          ncm_matrix_set (gauss->inv_cov, i, j, 
                          nc_bao_distance_priors_percival2010_inv_cov[i][j]);
      }
      break;
    }
    case NC_DATA_BAO_RDV_BEUTLER2011:
    {
      data->desc = g_strdup ("6dFGRS -- Beutler (2011), BAO Sample R-Dv");
      nc_data_bao_rdv_set_size (bao_rdv, 1);

      for (i = 0; i < 1; i++)
      {
        ncm_vector_set (bao_rdv->x, i, nc_bao_distance_priors_beutler2011_z[i]);
        ncm_vector_set (gauss->y,   i, nc_bao_distance_priors_beutler2011_bestfit[i]);
        for (j = 0; j < 1; j++)
          ncm_matrix_set (gauss->inv_cov, i, j, 
                          nc_bao_distance_priors_beutler2011_inv_cov[i][j]);
      }
      break;
    }
    case NC_DATA_BAO_RDV_PADMANABHAN2012:
    {
      data->desc = g_strdup ("SDSS-DR7-rec -- Padmanabhan (2012), BAO Sample R-Dv");
      nc_data_bao_rdv_set_size (bao_rdv, 1);

      for (i = 0; i < 1; i++)
      {
        ncm_vector_set (bao_rdv->x, i, nc_bao_distance_priors_padmanabhan2012_z[i]);
        ncm_vector_set (gauss->y,   i, nc_bao_distance_priors_padmanabhan2012_bestfit[i]);
        for (j = 0; j < 1; j++)
          ncm_matrix_set (gauss->inv_cov, i, j, 
                          nc_bao_distance_priors_padmanabhan2012_inv_cov[i][j]);
      }
      break;
    }
    case NC_DATA_BAO_RDV_ANDERSON2012:
    {
      data->desc = g_strdup ("SDSS-DR9-rec -- Anderson (2012), BAO Sample R-Dv");
      nc_data_bao_rdv_set_size (bao_rdv, 1);

      for (i = 0; i < 1; i++)
      {
        ncm_vector_set (bao_rdv->x, i, nc_bao_distance_priors_anderson2012_z[i]);
        ncm_vector_set (gauss->y,   i, nc_bao_distance_priors_anderson2012_bestfit[i]);
        for (j = 0; j < 1; j++)
          ncm_matrix_set (gauss->inv_cov, i, j, 
                          nc_bao_distance_priors_anderson2012_inv_cov[i][j]);
      }
      break;
    }
    case NC_DATA_BAO_RDV_BLAKE2012:
    {
      data->desc = g_strdup ("WiggleZ -- Blake (2012), BAO Sample R-Dv");
      nc_data_bao_rdv_set_size (bao_rdv, 3);

      for (i = 0; i < 3; i++)
      {
        ncm_vector_set (bao_rdv->x, i, nc_bao_distance_priors_blake2012_z[i]);
        ncm_vector_set (gauss->y,   i, nc_bao_distance_priors_blake2012_bestfit[i]);
        for (j = 0; j < 3; j++)
          ncm_matrix_set (gauss->inv_cov, i, j, 
                          nc_bao_distance_priors_blake2012_inv_cov[i][j]);
      }
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
  
  bao_rdv->id = id;
  ncm_data_set_init (data);
}

/**
 * nc_data_bao_rdv_get_sample:
 * @bao_rdv: a #NcDataBaoRDV
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcDataBaoId
nc_data_bao_rdv_get_sample (NcDataBaoRDV *bao_rdv)
{
  return bao_rdv->id;
}
