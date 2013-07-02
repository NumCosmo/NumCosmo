/***************************************************************************
 *            nc_data_snia_cov.c
 *
 *  Sat December 08 15:58:15 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:nc_data_snia_cov
 * @title: Supernovae Ia Data -- Covariance
 * @short_description: SNIa data with covariance error matrix
 * 
 * See #NcSNIADistCov.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_snia_cov.h"
#include "math/ncm_model_ctrl.h"

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataSNIACov, nc_data_snia_cov, NCM_TYPE_DATA_GAUSS_COV);

static void
nc_data_snia_cov_init (NcDataSNIACov *snia_cov)
{
  snia_cov->dcov_ctrl = ncm_model_ctrl_new (NULL);
}

static void
nc_data_snia_cov_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
/*  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object); */
  g_return_if_fail (NC_IS_DATA_SNIA_COV (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_snia_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
/*  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object); */
  g_return_if_fail (NC_IS_DATA_SNIA_COV (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_snia_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->constructed (object);

  ncm_data_set_init (NCM_DATA (object));
}

static void
nc_data_snia_cov_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->dispose (object);
}

static void
nc_data_snia_cov_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->finalize (object);
}

static void _nc_data_snia_cov_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_snia_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static gboolean _nc_data_snia_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov);

static void
nc_data_snia_cov_class_init (NcDataSNIACovClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass* gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_snia_cov_set_property;
  object_class->get_property = &nc_data_snia_cov_get_property;
  object_class->constructed  = &nc_data_snia_cov_constructed;
  object_class->dispose      = &nc_data_snia_cov_dispose;
  object_class->finalize     = &nc_data_snia_cov_finalize;
  
  data_class->prepare    = &_nc_data_snia_cov_prepare;
  gauss_class->mean_func = &_nc_data_snia_cov_mean_func;
  gauss_class->cov_func  = &_nc_data_snia_cov_func;
  
}

static void
_nc_data_snia_cov_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (data);
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  nc_snia_dist_cov_prepare_if_needed (dcov, mset);
  if (ncm_model_ctrl_model_update (snia_cov->dcov_ctrl, NCM_MODEL (dcov)))
    nc_data_snia_cov_set_dcov (snia_cov, dcov);
}

static void 
_nc_data_snia_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  nc_snia_dist_cov_mean (dcov, cosmo, vp);
}

static gboolean 
_nc_data_snia_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  nc_snia_dist_cov_calc (dcov, cov);
  
  return TRUE;
}

/**
 * nc_data_snia_cov_new:
 * @use_det: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_snia_cov_new (gboolean use_det)
{
  return g_object_new (NC_TYPE_DATA_SNIA_COV,
                       "use-det", use_det,
                       NULL);
}

/**
 * nc_data_snia_cov_set_dcov:
 * @snia_cov: FIXME
 * @dcov: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void 
nc_data_snia_cov_set_dcov (NcDataSNIACov *snia_cov, NcSNIADistCov *dcov)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (snia_cov);

  ncm_data_gauss_cov_set_size (gauss, dcov->mu_len);
  if (dcov->mu_len > 0)
  {
    ncm_vector_memcpy (gauss->y, dcov->mag);
    ncm_data_set_init (NCM_DATA (snia_cov));
  }
}
