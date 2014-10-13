/***************************************************************************
 *            nc_abc_cluster_ncount.c
 *
 *  Mon October 13 13:28:53 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_abc_cluster_ncount.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_abc_cluster_ncount
 * @title: Monte Carlo ABC analysis for cluster number counts
 * @short_description: Object implementing Approximate Bayesian Computation (ABC) for cluster number counts 
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "abc/nc_abc_cluster_ncount.h"
#include "nc_data_cluster_ncount.h"
#include "math/ncm_mset_trans_kern_gauss.h"

enum
{
  PROP_0,
  PROP_SCALE_COV,
  PROP_QUANTILES,
};


G_DEFINE_TYPE (NcABCClusterNCount, nc_abc_cluster_ncount, NCM_TYPE_ABC);

static void
nc_abc_cluster_ncount_init (NcABCClusterNCount *abcnc)
{
  abcnc->data_summary = NULL;
  abcnc->covar        = NULL;
  abcnc->scale_cov    = FALSE;
  abcnc->quantiles    = NULL;
}

static void
_nc_abc_cluster_ncount_constructed (GObject *object)
{
  NcmABC *abc = NCM_ABC (object);
  
  NcmMSetTransKernGauss *tkerng = ncm_mset_trans_kern_gauss_new (0);
  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (tkerng), abc->mcat->mset);
  ncm_mset_trans_kern_gauss_set_cov_from_scale (tkerng);
  ncm_abc_set_trans_kern (abc, NCM_MSET_TRANS_KERN (tkerng));
}

static void
_nc_abc_cluster_ncount_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (object);
  g_return_if_fail (NC_IS_ABC_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_SCALE_COV:
      nc_abc_cluster_ncount_set_scale_cov (abcnc, g_value_get_boolean (value));
      break;
    case PROP_QUANTILES:
    {
      GVariant *var = g_value_get_variant (value);
      NcmVector *v = ncm_vector_new_variant (var);
      ncm_vector_clear (&abcnc->quantiles);
      abcnc->quantiles = v;
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_abc_cluster_ncount_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (object);
  g_return_if_fail (NC_IS_ABC_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_SCALE_COV:
      g_value_set_boolean (value, abcnc->scale_cov);
      break;
    case PROP_QUANTILES:
    {
      if (abcnc->quantiles != NULL)
      {
        GVariant *var = ncm_vector_peek_variant (abcnc->quantiles); 
        g_value_take_variant (value, var);
      }
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_abc_cluster_ncount_dispose (GObject *object)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (object);
  
  ncm_matrix_clear (&abcnc->covar);
  ncm_vector_clear (&abcnc->quantiles);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_abc_cluster_ncount_parent_class)->dispose (object);
}

static void
_nc_abc_cluster_ncount_finalize (GObject *object)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (object);
  
  g_clear_pointer (&abcnc->data_summary, gsl_histogram2d_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_abc_cluster_ncount_parent_class)->finalize (object);
}

static gboolean _nc_abc_cluster_ncount_data_summary (NcmABC *abc);
static gdouble _nc_abc_cluster_ncount_mock_distance (NcmABC *abc, NcmDataset *dset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static gdouble _nc_abc_cluster_ncount_distance_prob (NcmABC *abc, gdouble distance);
static void _nc_abc_cluster_ncount_update_tkern (NcmABC *abc);
static const gchar *_nc_abc_cluster_ncount_get_desc (NcmABC *abc);
static const gchar *_nc_abc_cluster_ncount_log_info (NcmABC *abc);

static void
nc_abc_cluster_ncount_class_init (NcABCClusterNCountClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmABCClass *abc_class = NCM_ABC_CLASS (klass);

  object_class->constructed  = &_nc_abc_cluster_ncount_constructed;
  object_class->set_property = &_nc_abc_cluster_ncount_set_property;
  object_class->get_property = &_nc_abc_cluster_ncount_get_property;
  object_class->dispose  = &_nc_abc_cluster_ncount_dispose;
  object_class->finalize = &_nc_abc_cluster_ncount_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SCALE_COV,
                                   g_param_spec_boolean ("scale-cov",
                                                         NULL,
                                                         "Scaled covariance",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  abc_class->data_summary  = &_nc_abc_cluster_ncount_data_summary;
  abc_class->mock_distance = &_nc_abc_cluster_ncount_mock_distance;
  abc_class->distance_prob = &_nc_abc_cluster_ncount_distance_prob;
  abc_class->update_tkern  = &_nc_abc_cluster_ncount_update_tkern;
  abc_class->get_desc      = &_nc_abc_cluster_ncount_get_desc;
  abc_class->log_info      = &_nc_abc_cluster_ncount_log_info;
}

static gboolean 
_nc_abc_cluster_ncount_data_summary (NcmABC *abc)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (abc);
  g_assert (abc->dset != NULL);
  g_assert_cmpuint (ncm_dataset_get_ndata (abc->dset), ==, 1);
  {
    NcmData *data = ncm_dataset_get_data (abc->dset, 0);
    NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);

    if (abcnc->quantiles == NULL)
    {
      gdouble quantiles_data[7] = {0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98};
      NcmVector *quantiles_static = ncm_vector_new_data_static (quantiles_data, 7, 1);
      abcnc->quantiles = ncm_vector_dup (quantiles_static);
      ncm_vector_free (quantiles_static);
    }
    
    g_assert (NC_IS_DATA_CLUSTER_NCOUNT (data));

    nc_data_cluster_ncount_set_bin_by_quantile (ncount, abcnc->quantiles, abcnc->quantiles);

    g_clear_pointer (&abcnc->data_summary, gsl_histogram2d_free);
    abcnc->data_summary = gsl_histogram2d_clone (ncount->z_lnM);
  }
  

  return TRUE;
}

static gdouble 
_nc_abc_cluster_ncount_mock_distance (NcmABC *abc, NcmDataset *dset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (abc);
  NcmData *data = ncm_dataset_get_data (dset, 0);
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);

  const gsize nx = gsl_histogram2d_nx (abcnc->data_summary);
  const gsize ny = gsl_histogram2d_ny (abcnc->data_summary);
  const gsize total = nx * ny;
  gsize i; 

  gdouble res = 0.0;
  
  for (i = 0; i < total; i++)
  {
    res += gsl_pow_2 (abcnc->data_summary->bin[i] - ncount->z_lnM->bin[i]); 
  }
  
  return sqrt (res) / total;
}

static gdouble 
_nc_abc_cluster_ncount_distance_prob (NcmABC *abc, gdouble distance)
{
  if (distance < abc->epsilon)
    return 1.0;
  else
    return 0.0;
}

static void 
_nc_abc_cluster_ncount_update_tkern (NcmABC *abc)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (abc);
  const gdouble scale = abcnc->scale_cov ? 0.25 + 1.75 * ncm_abc_get_accept_rate (abc) : 2.0;
  
  ncm_mset_catalog_get_covar (abc->mcat, &abc->covar);
  
  ncm_matrix_scale (abc->covar, scale);
  ncm_mset_trans_kern_gauss_set_cov (NCM_MSET_TRANS_KERN_GAUSS (abc->tkern), abc->covar);

  if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("# NcABCClusterNCount: scale covariance by %f\n", scale);

  {
    gdouble dist_75 = ncm_abc_get_dist_quantile (abc, 0.75);
    ncm_abc_update_epsilon (abc, dist_75);
  }
}

static const gchar *
_nc_abc_cluster_ncount_get_desc (NcmABC *abc)
{
  return "NcABCClusterNCount";
}
static const gchar *
_nc_abc_cluster_ncount_log_info (NcmABC *abc)
{
  return ncm_dataset_get_info (abc->dset);
}

/**
 * nc_abc_cluster_ncount_new:
 * @mset: a #NcmMSet.
 * @prior: a #NcmMSetTransKern.
 * @dset: a #NcmDataset.
 * 
 * Creates a new #NcABCClusterNCount.
 * 
 * Returns: (transfer full): a new #NcABCClusterNCount.
 */
NcABCClusterNCount *
nc_abc_cluster_ncount_new (NcmMSet *mset, NcmMSetTransKern *prior, NcmDataset *dset)
{
  NcABCClusterNCount *abcnc = g_object_new (NC_TYPE_ABC_CLUSTER_NCOUNT, 
                                            "mset", mset,
                                            "prior", prior,
                                            "data-set", dset,
                                            NULL);
  return abcnc;
}

/**
 * nc_abc_cluster_ncount_set_scale_cov:
 * @abcnc: a #NcABCClusterNCount.
 * @on: whether sets on or off covariance scaling.
 * 
 * FIXME
 * 
 */
void 
nc_abc_cluster_ncount_set_scale_cov (NcABCClusterNCount *abcnc, gboolean on)
{
  abcnc->scale_cov = on;
}
