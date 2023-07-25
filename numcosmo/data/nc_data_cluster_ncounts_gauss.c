/***************************************************************************
 *            nc_data_cluster_ncounts_gauss.c
 *
 *  Tue Apr  6 01:11:23 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_cluster_ncounts_gauss
 * @title: NcDataClusterNCountsGauss
 * @short_description: Cluster number count data gaussian likelihood.
 *
 * #NcDataClusterNCountsGauss is a #NcmDataGaussCov that implements the
 * gaussian likelihood for the cluster number count data.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_ncounts_gauss.h"
#include "nc_hireion.h"

#include "math/ncm_func_eval.h"
#include "math/ncm_serialize.h"
#include "math/ncm_obj_array.h"
#include "math/ncm_cfg.h"

#include "misc/cubature.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <glib/gstdio.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_Z_OBS,
  PROP_Z_OBS_PARAMS,
  PROP_LNM_OBS,
  PROP_LNM_OBS_PARAMS,
  PROP_CAD,
  PROP_HAS_SSC,
  PROP_S_MATRIX,
  PROP_BIN_COUNT,
  PROP_SIZE,
};

struct _NcDataClusterNCountsGaussPrivate
{
  NcmVector *z_obs;
  NcmMatrix *z_obs_params;
  NcmVector *lnM_obs;
  NcmMatrix *lnM_obs_params;
  NcClusterAbundance *cad;
  gboolean has_ssc;
  NcmMatrix *s_matrix;
  NcmVector *bin_count;
  GArray *index_map;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterNCountsGauss, nc_data_cluster_ncounts_gauss, NCM_TYPE_DATA_GAUSS_COV);

typedef struct _NcDataClusterNCountsGaussIndex
{
  guint i_z;
  guint i_lnM;
  gdouble *z_obs_lb;
  gdouble *z_obs_ub;
  gdouble *lnM_obs_lb;
  gdouble *lnM_obs_ub;
} NcDataClusterNCountsGaussIndex;

static void
nc_data_cluster_ncounts_gauss_init (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv = nc_data_cluster_ncounts_gauss_get_instance_private (ncounts_gauss);

  self->z_obs          = NULL;
  self->z_obs_params   = NULL;
  self->lnM_obs        = NULL;
  self->lnM_obs_params = NULL;
  self->cad            = NULL;
  self->has_ssc        = FALSE;
  self->s_matrix       = NULL;
  self->bin_count      = NULL;
  self->index_map      = g_array_new (FALSE, FALSE, sizeof (NcDataClusterNCountsGaussIndex));
}

static void
nc_data_cluster_ncounts_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (object);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  g_return_if_fail (NC_DATA_CLUSTER_NCOUNTS_GAUSS (object));

  switch (prop_id)
  {
    case PROP_Z_OBS:
      nc_data_cluster_ncounts_gauss_set_z_obs (ncounts_gauss, g_value_get_object (value));
      break;
    case PROP_Z_OBS_PARAMS:
      nc_data_cluster_ncounts_gauss_set_z_obs_params (ncounts_gauss, g_value_get_object (value));
      break;
    case PROP_LNM_OBS:
      nc_data_cluster_ncounts_gauss_set_lnM_obs (ncounts_gauss, g_value_get_object (value));
      break;
    case PROP_LNM_OBS_PARAMS:
      nc_data_cluster_ncounts_gauss_set_lnM_obs_params (ncounts_gauss, g_value_get_object (value));
      break;
    case PROP_CAD:
      nc_cluster_abundance_clear (&self->cad);
      self->cad = g_value_dup_object (value);
      break;
    case PROP_HAS_SSC:
      nc_data_cluster_ncounts_gauss_set_has_ssc (ncounts_gauss, g_value_get_boolean (value));
      break;
    case PROP_S_MATRIX:
      nc_data_cluster_ncounts_gauss_set_s_matrix (ncounts_gauss,  g_value_get_object (value));
    case PROP_BIN_COUNT:
      nc_data_cluster_ncounts_gauss_set_bin_count (ncounts_gauss, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_ncounts_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (object);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  g_return_if_fail (NC_DATA_CLUSTER_NCOUNTS_GAUSS (object));

  switch (prop_id)
  {
    case PROP_Z_OBS:
      g_value_set_object (value, self->z_obs);
      break;
    case PROP_Z_OBS_PARAMS:
      g_value_set_object (value, self->z_obs_params);
      break;
    case PROP_LNM_OBS:
      g_value_set_object (value, self->lnM_obs);
      break;
    case PROP_LNM_OBS_PARAMS:
      g_value_set_object (value, self->lnM_obs_params);
      break;
    case PROP_CAD:
      g_value_set_object (value, self->cad);
      break;
    case PROP_HAS_SSC:
      g_value_set_boolean (value, self->has_ssc);
      break;
    case PROP_S_MATRIX:
      g_value_set_object (value, self->s_matrix);
      break;
    case PROP_BIN_COUNT:
      g_value_set_object (value, self->bin_count);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_ncounts_gauss_dispose (GObject *object)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (object);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_vector_clear (&self->z_obs);
  ncm_matrix_clear (&self->z_obs_params);
  ncm_vector_clear (&self->lnM_obs);
  ncm_matrix_clear (&self->lnM_obs_params);
  nc_cluster_abundance_clear (&self->cad);
  ncm_matrix_clear (&self->s_matrix);
  ncm_vector_clear (&self->bin_count);

  g_clear_pointer (&self->index_map, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncounts_gauss_parent_class)->dispose (object);
}

static void
nc_data_cluster_ncounts_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncounts_gauss_parent_class)->finalize (object);
}

static void _nc_data_cluster_ncounts_gauss_begin (NcmData *data);
static void _nc_data_cluster_ncounts_gauss_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_ncounts_gauss_set_size (NcmDataGaussCov *gauss_cov, guint np);
static void _nc_data_cluster_ncounts_gauss_mean_func (NcmDataGaussCov *gauss_cov, NcmMSet *mset, NcmVector *vp);
static gboolean _nc_data_cluster_ncounts_gauss_cov_func (NcmDataGaussCov *gauss_cov, NcmMSet *mset, NcmMatrix *cov);

static void
nc_data_cluster_ncounts_gauss_class_init (NcDataClusterNCountsGaussClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class              = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_cov_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_cluster_ncounts_gauss_set_property;
  object_class->get_property = &nc_data_cluster_ncounts_gauss_get_property;
  object_class->dispose      = &nc_data_cluster_ncounts_gauss_dispose;
  object_class->finalize     = &nc_data_cluster_ncounts_gauss_finalize;

  g_object_class_install_property (object_class,
                                   PROP_Z_OBS,
                                   g_param_spec_object ("z-obs",
                                                        NULL,
                                                        "Clusters redshift observables",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_OBS_PARAMS,
                                   g_param_spec_object ("z-obs-params",
                                                        NULL,
                                                        "Clusters redshift observables parameters",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS,
                                   g_param_spec_object ("lnM-obs",
                                                        NULL,
                                                        "Clusters mass observables",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS_PARAMS,
                                   g_param_spec_object ("lnM-obs-params",
                                                        NULL,
                                                        "Clusters mass observables parameters",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CAD,
                                   g_param_spec_object ("cluster-abundance",
                                                        NULL,
                                                        "Cluster abundance",
                                                        NC_TYPE_CLUSTER_ABUNDANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HAS_SSC,
                                   g_param_spec_boolean ("has-ssc",
                                                         NULL,
                                                         "Whether use super sample covariance",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_S_MATRIX,
                                   g_param_spec_object ("s-matrix",
                                                        NULL,
                                                        "Super sample covariance matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_COUNT,
                                   g_param_spec_object ("bin-count",
                                                        NULL,
                                                        "Bin count",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  data_class->begin          = &_nc_data_cluster_ncounts_gauss_begin;
  data_class->prepare        = &_nc_data_cluster_ncounts_gauss_prepare;
  gauss_cov_class->set_size  = &_nc_data_cluster_ncounts_gauss_set_size;
  gauss_cov_class->mean_func = &_nc_data_cluster_ncounts_gauss_mean_func;
  gauss_cov_class->cov_func  = &_nc_data_cluster_ncounts_gauss_cov_func;
}

static void
_nc_data_cluster_ncounts_gauss_begin (NcmData *data)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (data);
  NcmDataGaussCov *gauss_cov                    = NCM_DATA_GAUSS_COV (data);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;
  const guint np                                = ncm_data_gauss_cov_get_size (gauss_cov);
  const gint nz                                 = ncm_vector_len (self->z_obs) - 1;
  const gint nlnM                               = ncm_vector_len (self->lnM_obs) - 1;

  g_assert_cmpint (np, ==, nz * nlnM);

  g_array_set_size (self->index_map, 0);

  {
    gint i_z, i_lnM;

    for (i_z = 0; i_z < nz; i_z++)
    {
      for (i_lnM = 0; i_lnM < nlnM; i_lnM++)
      {
        NcDataClusterNCountsGaussIndex k;

        k.i_z        = i_z;
        k.i_lnM      = i_lnM;
        k.z_obs_lb   = ncm_vector_ptr (self->z_obs, i_z + 0);
        k.z_obs_ub   = ncm_vector_ptr (self->z_obs, i_z + 1);
        k.lnM_obs_lb = ncm_vector_ptr (self->lnM_obs, i_lnM + 0);
        k.lnM_obs_ub = ncm_vector_ptr (self->lnM_obs, i_lnM + 1);

        g_array_append_val (self->index_map, k);
      }
    }
  }
}

static void
_nc_data_cluster_ncounts_gauss_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (data);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;
  NcHICosmo *cosmo                              = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz                   = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm                       = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  g_assert ((cosmo != NULL) && (clusterz != NULL) && (clusterm != NULL));

  nc_cluster_abundance_prepare_if_needed (self->cad, cosmo, clusterz, clusterm);
}

static void
_nc_data_cluster_ncounts_gauss_set_size (NcmDataGaussCov *gauss_cov, guint np)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (gauss_cov);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;
  const guint cnp                               = ncm_data_gauss_cov_get_size (gauss_cov);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&self->bin_count);

  if ((np != 0) && (np != cnp))
    self->bin_count = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_cluster_ncounts_gauss_parent_class)->set_size (gauss_cov, np);
}

/**
 * _nc_data_cluster_ncounts_gauss_mean_func:
 * @filename: file containing a serialized #NcDataClusterNCountsGauss.
 *
 * Calculates the expected number of clusters inside each bin of mass and redshift
 */

static void
_nc_data_cluster_ncounts_gauss_mean_func (NcmDataGaussCov *gauss_cov, NcmMSet *mset, NcmVector *vp)
{
  NcDataClusterNCountsGauss *ncounts_gauss      = NC_DATA_CLUSTER_NCOUNTS_GAUSS (gauss_cov);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;
  NcHICosmo *cosmo                              = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz                   = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm                       = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterAbundance *cad                       = self->cad;
  gint i;


  for (i = 0; i < self->index_map->len; i++)
  {
    const NcDataClusterNCountsGaussIndex *k = &g_array_index (self->index_map, NcDataClusterNCountsGaussIndex, i);

    const gdouble mean = nc_cluster_abundance_intp_bin_d2n (cad, cosmo, clusterz, clusterm,
                                                            k->lnM_obs_lb,
                                                            k->lnM_obs_ub,
                                                            NULL,
                                                            k->z_obs_lb,
                                                            k->z_obs_ub,
                                                            NULL);

    ncm_vector_set (vp, i, mean);
  }

  return;
}

/**
 * _nc_data_cluster_ncounts_gauss_cov_func:
 * @filename: file containing a serialized #NcDataClusterNCountsGauss.
 *
 * Calculates the covariance of the number clusters between each bin of mass and redshift
 */

static gboolean
_nc_data_cluster_ncounts_gauss_cov_func (NcmDataGaussCov *gauss_cov, NcmMSet *mset, NcmMatrix *cov)
{
  NcDataClusterNCountsGauss *ncounts_gauss = NC_DATA_CLUSTER_NCOUNTS_GAUSS (gauss_cov);
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterAbundance *cad = self->cad;
  guint i, j;

  ncm_matrix_set_zero (cov);

  if (self->has_ssc)
  {
    for (i = 0; i < self->index_map->len; i++)
    {
      const NcDataClusterNCountsGaussIndex *k_i = &g_array_index (self->index_map, NcDataClusterNCountsGaussIndex, i);
      const gdouble bias_i                      = nc_cluster_abundance_intp_bin_d2n_bias (cad, cosmo, clusterz, clusterm,
                                                                                          k_i->lnM_obs_lb,
                                                                                          k_i->lnM_obs_ub,
                                                                                          NULL,
                                                                                          k_i->z_obs_lb,
                                                                                          k_i->z_obs_ub,
                                                                                          NULL);

      for (j = 0; j < self->index_map->len; j++)
      {
        const NcDataClusterNCountsGaussIndex *k_j = &g_array_index (self->index_map, NcDataClusterNCountsGaussIndex, j);
        const gdouble Sij                         = ncm_matrix_get (self->s_matrix, k_i->i_z, k_j->i_z);


        if (i == j)
        {
          const gdouble poisson_i = nc_cluster_abundance_intp_bin_d2n (cad, cosmo, clusterz, clusterm,
                                                                       k_i->lnM_obs_lb,
                                                                       k_i->lnM_obs_ub,
                                                                       NULL,
                                                                       k_i->z_obs_lb,
                                                                       k_i->z_obs_ub,
                                                                       NULL);

          ncm_matrix_set (cov, i, j, poisson_i + bias_i * bias_i * Sij);
        }
        else
        {
          const gdouble bias_j = nc_cluster_abundance_intp_bin_d2n_bias (cad, cosmo, clusterz, clusterm,
                                                                         k_j->lnM_obs_lb,
                                                                         k_j->lnM_obs_ub,
                                                                         NULL,
                                                                         k_j->z_obs_lb,
                                                                         k_j->z_obs_ub,
                                                                         NULL);

          ncm_matrix_set (cov, i, j, bias_i * bias_j * Sij);
        }
      }
    }
  }
  else
  {
    for (i = 0; i < self->index_map->len; i++)
    {
      const NcDataClusterNCountsGaussIndex *k_i = &g_array_index (self->index_map, NcDataClusterNCountsGaussIndex, i);
      const gdouble poisson_i                   = nc_cluster_abundance_intp_bin_d2n (cad, cosmo, clusterz, clusterm,
                                                                                     k_i->lnM_obs_lb,
                                                                                     k_i->lnM_obs_ub,
                                                                                     NULL,
                                                                                     k_i->z_obs_lb,
                                                                                     k_i->z_obs_ub,
                                                                                     NULL);

      ncm_matrix_set (cov, i, i, poisson_i);
    }
  }


  return TRUE;
}

/**
 * nc_data_cluster_ncounts_gauss_new:
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 * Returns: NcDataClusterNCountsGauss
 */
NcDataClusterNCountsGauss *
nc_data_cluster_ncounts_gauss_new (NcClusterAbundance *cad)
{
  NcDataClusterNCountsGauss *ncounts_gauss = g_object_new (NC_TYPE_DATA_CLUSTER_NCOUNTS_GAUSS,
                                                           "cluster-abundance", cad,
                                                           NULL);

  return ncounts_gauss;
}

/**
 * nc_data_cluster_ncounts_gauss_set_z_obs:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @z_obs: a #NcmVector
 *
 * Sets array of #NcmVector's representing the lower and upper bounds
 * of each bin.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_z_obs (NcDataClusterNCountsGauss *ncounts_gauss, NcmVector *z_obs)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_vector_clear (&self->z_obs);
  self->z_obs = ncm_vector_ref (z_obs);
}

/**
 * nc_data_cluster_ncounts_gauss_set_z_obs_params:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @z_obs_params: a #NcmMatrix
 *
 * Sets array of #NcmVector's representing the params
 * of each bin.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_z_obs_params (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *z_obs_params)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_matrix_clear (&self->z_obs_params);
  self->z_obs_params = ncm_matrix_ref (z_obs_params);
}

/**
 * nc_data_cluster_ncounts_gauss_set_lnM_obs:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @lnM_obs: a #NcmVector
 *
 * Sets array of #NcmVector's representing the lower and upper bounds
 * of each bin.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_lnM_obs (NcDataClusterNCountsGauss *ncounts_gauss, NcmVector *lnM_obs)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_vector_clear (&self->lnM_obs);
  self->lnM_obs = ncm_vector_ref (lnM_obs);
}

/**
 * nc_data_cluster_ncounts_gauss_set_lnM_obs_params:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @lnM_obs_params: a #NcmMatrix
 *
 * Sets array of #NcmVector's representing the params
 * of each bin.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_lnM_obs_params (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *lnM_obs_params)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_matrix_clear (&self->lnM_obs_params);
  self->lnM_obs_params = ncm_matrix_ref (lnM_obs_params);
}

/**
 * nc_data_cluster_ncounts_gauss_set_has_ssc:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @on: FIXME
 *
 * Sets array of #Set if the data has super sample covariance.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_has_ssc (NcDataClusterNCountsGauss *ncounts_gauss, gboolean on)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  self->has_ssc = on;
}

/**
 * nc_data_cluster_ncounts_gauss_set_s_matrix:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @s_matrix: a #NcmMatrix
 *
 * Sets array of #NcmVector's representing the super sample covariance effect in each bin.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_s_matrix (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *s_matrix)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_matrix_clear (&self->s_matrix);
  self->s_matrix = ncm_matrix_ref (s_matrix);
}

/**
 * nc_data_cluster_ncounts_gauss_set_bin_count:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 * @bin_count: a #NcmVector
 *
 * Sets array of #NcmVector's representing the observed number of clusters in each bin.
 *
 */
void
nc_data_cluster_ncounts_gauss_set_bin_count (NcDataClusterNCountsGauss *ncounts_gauss, NcmVector *bin_count)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  ncm_vector_clear (&self->bin_count);
  self->bin_count = ncm_vector_ref (bin_count);
}

/**
 * nc_data_cluster_ncounts_gauss_get_z_obs:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets the matrix containing the redshift observables.
 *
 * Returns: (transfer full): Redshift observable #NcmVector.
 */
NcmVector *
nc_data_cluster_ncounts_gauss_get_z_obs (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  if (self->z_obs != NULL)
    return ncm_vector_ref (self->z_obs);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncounts_gauss_get_z_obs_params:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets the matrix containing the redshift observables parameters.
 *
 * Returns: (transfer full): Redshift observable parameters #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncounts_gauss_get_z_obs_params (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  if (self->z_obs_params != NULL)
    return ncm_matrix_ref (self->z_obs_params);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncounts_gauss_get_lnM_obs:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets the matrix containing the mass observables.
 *
 * Returns: (transfer full): Mass observable #NcmMatrix.
 */
NcmVector *
nc_data_cluster_ncounts_gauss_get_lnM_obs (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  if (self->lnM_obs != NULL)
    return ncm_vector_ref (self->lnM_obs);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncounts_gauss_get_lnM_obs_params:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets the matrix containing the mass observables parameters.
 *
 * Returns: (transfer full): Mass observable parameters #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncounts_gauss_get_lnM_obs_params (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  if (self->lnM_obs_params != NULL)
    return ncm_matrix_ref (self->lnM_obs_params);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncounts_gauss_get_has_ssc:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets if the ssc option is on.
 *
 * Returns: TRUE or FALSE.
 */
gboolean
nc_data_cluster_ncounts_gauss_get_has_ssc (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  return self->has_ssc;
}

/**
 * nc_data_cluster_ncounts_gauss_get_s_matrix:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets the matrix containing the super sample covariance.
 *
 * Returns: (transfer full): Super sample covariance #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncounts_gauss_get_s_matrix (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  if (self->s_matrix != NULL)
    return ncm_matrix_ref (self->s_matrix);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncounts_gauss_get_bin_count:
 * @ncounts_gauss: a #NcDataClusterNCountsGauss
 *
 * Gets the vector containing the number of clusters in each bin
 *
 * Returns: (transfer full): The observed clusters in each bin #NcmVector.
 */
NcmVector *
nc_data_cluster_ncounts_gauss_get_bin_count (NcDataClusterNCountsGauss *ncounts_gauss)
{
  NcDataClusterNCountsGaussPrivate * const self = ncounts_gauss->priv;

  if (self->bin_count != NULL)
    return ncm_vector_ref (self->bin_count);
  else
    return NULL;
}

