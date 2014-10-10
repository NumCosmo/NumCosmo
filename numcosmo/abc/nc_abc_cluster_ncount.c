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
  PROP_EPSILON
};


G_DEFINE_TYPE (NcABCClusterNCount, nc_abc_cluster_ncount, NCM_TYPE_ABC);

static void
nc_abc_cluster_ncount_init (NcABCClusterNCount *abcnc)
{
  abcnc->epsilon = 0.0;
  abcnc->data_summary = NULL;
  abcnc->covar = NULL;
}

static void
_nc_abc_cluster_ncount_constructed (GObject *object)
{
  NcmABC *abc = NCM_ABC (object);
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (object);
  
  NcmMSetTransKernGauss *tkerng = ncm_mset_trans_kern_gauss_new (0);
  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (tkerng), abc->mcat->mset);
  ncm_mset_trans_kern_gauss_set_cov_from_scale (tkerng);
  ncm_abc_set_trans_kern (abc, NCM_MSET_TRANS_KERN (tkerng));

  abcnc->epsilon = 1e20;
}

static void
_nc_abc_cluster_ncount_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (object);
  g_return_if_fail (NC_IS_ABC_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_EPSILON:
      abcnc->epsilon = g_value_get_double (value);
      break;
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
    case PROP_EPSILON:
      g_value_set_double (value, abcnc->epsilon);
      break;
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
                                   PROP_EPSILON,
                                   g_param_spec_double ("epsilon",
                                                        NULL,
                                                        "epsilon",
                                                        0.0, G_MAXDOUBLE, 1.0,
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
    gdouble quantiles[7] = {0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98};
    NcmData *data = ncm_dataset_get_data (abc->dset, 0);
    NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
    NcmVector *quantiles_v = ncm_vector_new_data_static (quantiles, 7, 1);
    g_assert (NC_IS_DATA_CLUSTER_NCOUNT (data));

    nc_data_cluster_ncount_set_bin_by_quantile (ncount, quantiles_v, quantiles_v);

    g_clear_pointer (&abcnc->data_summary, gsl_histogram2d_free);
    abcnc->data_summary = gsl_histogram2d_clone (ncount->z_lnM);
    
    ncm_vector_free (quantiles_v);
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
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (abc);

  if (distance < abcnc->epsilon)
    return 1.0;
  else
    return 0.0;
}

static void 
_nc_abc_cluster_ncount_update_tkern (NcmABC *abc)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (abc);
  ncm_mset_catalog_get_covar (abc->mcat, &abc->covar);
  
  ncm_matrix_scale (abc->covar, 2.0);
  ncm_mset_trans_kern_gauss_set_cov (NCM_MSET_TRANS_KERN_GAUSS (abc->tkern), abc->covar);

  {
    gdouble dist_75 = ncm_abc_get_dist_quantile (abc, 0.75);
    if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      guint i;
      g_message ("# NcABCClusterNCount: ");
      for (i = 0; i < 8; i++)
      {
        gdouble p = (15.0 * i) / 100.0;
        p = p > 1.0 ? 1.0 : p;
        g_message ("[%2.0f%% %4.2f] ", 100.0 * p, ncm_abc_get_dist_quantile (abc, p));
      }
      g_message ("\n");
      g_message ("# NcABCClusterNCount: epsilon_t   = %g.\n", 
                 abcnc->epsilon);
      g_message ("# NcABCClusterNCount: epsilon_t+1 = %g.\n", 
                 dist_75 );
      
    }
    abcnc->epsilon = dist_75;
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
