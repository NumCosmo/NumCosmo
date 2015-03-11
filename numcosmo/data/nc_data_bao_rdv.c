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
 * @title: NcDataBaoRDV
 * @short_description: Baryon Oscillation Data -- $r_s / D_V$ ratio.
 *
 * See <link linkend="XPercival2007">Percival et al. (2007)</link>.
 * 
 * Kazin et al. (arXiv:1401.0358): our implementation of the inverse covariance matrix is given by  
 * $$C^{-1}_{new} = \frac{1}{r_s^{\text{fid}}} C^{-1},$$
 * where $r_s^{\text{fid}} = 148.6$ and $C^{-1}$ is given in table 4. This modification is due the fact 
 * that we are using $D_V(z)/r_s(z_d)$ instead of $D_V(z)* r_s^{\text{fid}}/r_s(z_d)$. Analogously, we implemented 
 * $D_V(z) / r_s(z_d) = [1716.4, 2220.8, 2516.1] / 148.6 = [11.550, 14.945, 16.932]$. 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_rdv.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_DATA_FORM, 
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoRDV, nc_data_bao_rdv, NCM_TYPE_DATA_GAUSS);

static void
nc_data_bao_rdv_init (NcDataBaoRDV *bao_rdv)
{
  bao_rdv->x    = NULL;
  bao_rdv->dist = NULL;
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
      nc_distance_clear (&bao_rdv->dist);
      bao_rdv->dist = g_value_dup_object (value);
      break;
    case PROP_Z:
      ncm_vector_set_from_variant (bao_rdv->x, g_value_get_variant (value));
      break;
    case PROP_DATA_FORM:
      bao_rdv->r_DV = g_value_get_boolean (value);
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
    case PROP_Z:
      g_value_take_variant (value, ncm_vector_get_variant (bao_rdv->x));
      break;
    case PROP_DATA_FORM:
      g_value_set_boolean (value, bao_rdv->r_DV);
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
static void _nc_data_bao_rdv_set_size (NcmDataGauss *gauss, guint np);

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
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_variant ("z",
                                                         NULL,
                                                         "Data redshift",
                                                         G_VARIANT_TYPE ("ad"), NULL,
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
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (bao_rdv->dist, cosmo);
}

static void 
_nc_data_bao_rdv_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (gauss);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  guint i;

  if (bao_rdv->r_DV)
  {
    for (i = 0; i < gauss->np; i++)
    {
      const gdouble z  = ncm_vector_get (bao_rdv->x, i);
      const gdouble r_Dv = nc_distance_bao_r_Dv (bao_rdv->dist, cosmo, z);
      ncm_vector_set (vp, i, r_Dv);
    }
  }
  else
  {
    for (i = 0; i < gauss->np; i++)
    {
      const gdouble z  = ncm_vector_get (bao_rdv->x, i);
      const gdouble Dv_r = 1.0 / nc_distance_bao_r_Dv (bao_rdv->dist, cosmo, z);
      ncm_vector_set (vp, i, Dv_r);
    }
  }
}

/**
 * nc_data_bao_rdv_new:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: a #NcmData
 */
NcmData *
nc_data_bao_rdv_new (NcDistance *dist, NcDataBaoId id)
{
  NcmData *data = g_object_new (NC_TYPE_DATA_BAO_RDV,
                                "dist", dist,
                                NULL);
  nc_data_bao_rdv_set_sample (NC_DATA_BAO_RDV (data), id);
  return data;
}

static void 
_nc_data_bao_rdv_set_size (NcmDataGauss *gauss, guint np)
{
  NcDataBaoRDV *bao_rdv = NC_DATA_BAO_RDV (gauss);

  if ((np == 0) || (np != gauss->np))
    ncm_vector_clear (&bao_rdv->x);

  if ((np != 0) && (np != gauss->np))
    bao_rdv->x = ncm_vector_new (np);
  
  /* Chain up : end */
  NCM_DATA_GAUSS_CLASS (nc_data_bao_rdv_parent_class)->set_size (gauss, np);
}

#define _NC_DATA_BAO_RDV_MAX_LEN 5
#define _NC_DATA_BAO_RDV_MAX_LEN2 25

typedef struct _NcDataBaoRDVSample
{
  NcDataBaoId id;
  gboolean r_DV;
  guint len;
  gdouble z[_NC_DATA_BAO_RDV_MAX_LEN];
  gdouble bestfit[_NC_DATA_BAO_RDV_MAX_LEN];
  gdouble inv_cov[_NC_DATA_BAO_RDV_MAX_LEN2];
  gchar *desc;
} NcDataBaoRDVSample;

/***************************************************************************
 * Data of r_s/D_V and D_V/r_s
 ****************************************************************************/

/***************************************************************************
 * BAO Percival priors data (arXiv:0705.3323)
 * Idem (arXiv:0907.1660)
 * BAO 6dFGRS Beutler et al. (2011) (arXiv:1106.3366) - r_s computed using EH1998
 * BAO SDSS-DR7-rec Padmanabhan et al. (2012) (arXiv:1202.0090) - r_s computed using EH1998
 * BAO SDSS-DR9-rec Anderson et al. (2012) (arXiv:1203.6594) - r_s computed using EH1998
 * BAO WiggleZ Blake et al. (2012) arXiv:1108.2635 - r_s computed using EH1998 - Inverse covariance matrix arXiv:1212.5226
 * BAO WiggleZ Kazin et al. (2014) arXiv:1401.0358 - r_s computed using CAMB - r_s^{fid} = 148.6, D_V*(r_s^{fid}/r_s) = {1716.4, 2220.8, 2516.1}  
 ****************************************************************************/

NcDataBaoRDVSample nc_data_bao_rdv_samples[] = {
  { 
   NC_DATA_BAO_RDV_PERCIVAL2007, TRUE, 2,
   {    0.2,   0.35 }, 
   { 0.1980, 0.1094 }, 
   {
     35059.0, -24031.0,
    -24031.0, 108300.0 
   },
   "Percival 2007, BAO Sample R-Dv"},
  {
   NC_DATA_BAO_RDV_PERCIVAL2010, TRUE, 2,
   {    0.2,   0.35 }, 
   { 0.1905, 0.1097 }, 
   { 
     30124.0, -17227.0,
    -17227.0,  86977.0
   },
   "Percival 2010, BAO Sample R-Dv"},
   {
     NC_DATA_BAO_RDV_BEUTLER2011, TRUE, 1,
     { 0.106 },
     { 0.336 / 1.027}, /* Correction for r_{s,EH}/r_{sCAMB} = 1.027 see Kazin (2014) */
     { 1.0 * 1.027 * 1.027 / (0.015 * 0.015) },
     "6dFGRS -- Beutler (2011), BAO Sample R-Dv"},
   {
     NC_DATA_BAO_RDV_PADMANABHAN2012, FALSE, 1,
     { 0.35 },
     { 8.88 * 1.025}, /* Correction for r_{s,EH}/r_{sCAMB} = 1.025 see Kazin (2014) */
     { 1.0 / (0.17 * 0.17 * 1.025 * 1.025) },
     "SDSS-DR7-rec -- Padmanabhan (2012), BAO Sample R-Dv"},
   {
     NC_DATA_BAO_RDV_ANDERSON2012, FALSE, 1,
     { 0.57 },
     { 13.67 },
     { 1.0 / (0.22 * 0.22) },
     "SDSS-DR9-rec -- Anderson (2012), BAO Sample R-Dv"},
   {
     NC_DATA_BAO_RDV_BLAKE2012, TRUE, 3,
     {   0.44,   0.60,   0.73 },
     { 0.0916, 0.0726, 0.0592 },
     {  
       24532.1, -25137.7,  12099.1,
      -25137.7, 134598.4, -64783.9,
       12099.1, -64783.9, 128837.6
     },
     "WiggleZ -- Blake (2012), BAO Sample R-Dv"},
   {
     NC_DATA_BAO_RDV_KAZIN2014, FALSE, 3,
     {   0.44,   0.60,   0.73 },
     { 11.550, 14.945, 16.932 },
     {
       4.8116, -2.4651,  1.0375, /* e.g., first element: (148.6 * 148.6 * 2.179 * 1e-4) */
      -2.4651,  3.7697, -1.5865,
       1.0375, -1.5865,  3.6498
     },
     "WiggleZ -- Kazin (2014), BAO Sample R-Dv"},
};

/**
 * nc_data_bao_rdv_set_sample:
 * @bao_rdv: a #NcDataBaoRDV
 * @id: a #NcDataBaoId
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
  guint sindex = id - NC_DATA_BAO_RDV_FIRST;
  NcDataBaoRDVSample *sdata = &nc_data_bao_rdv_samples[sindex];
  
  g_assert (id >= NC_DATA_BAO_RDV_PERCIVAL2007 && id <= NC_DATA_BAO_RDV_KAZIN2014);

  if (data->desc != NULL)
    g_free (data->desc);

  data->desc = g_strdup (sdata->desc);  
  ncm_data_gauss_set_size (gauss, sdata->len);
  bao_rdv->r_DV = sdata->r_DV;

  for (i = 0; i < sdata->len; i++)
  {
    ncm_vector_set (bao_rdv->x, i, sdata->z[i]);
    ncm_vector_set (gauss->y,   i, sdata->bestfit[i]);
    for (j = 0; j < sdata->len; j++)
      ncm_matrix_set (gauss->inv_cov, i, j, sdata->inv_cov[i*sdata->len + j]);
  }
  g_assert_cmpuint (sdata->id, ==, id); 
  
  ncm_data_set_init (data, TRUE);
}
