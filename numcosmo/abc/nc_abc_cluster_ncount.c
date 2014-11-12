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
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_SCALE_COV,
  PROP_BIN_TYPE,
  PROP_QUANTILES,
  PROP_Z_NODES,
  PROP_LNM_NODES,
  PROP_Z_BINS,
  PROP_LNM_BINS,
  PROP_EPSILON_UPDATE,
  PROP_EPSILON_UPDATE_TYPE, 
};


G_DEFINE_TYPE (NcABCClusterNCount, nc_abc_cluster_ncount, NCM_TYPE_ABC);

static void
nc_abc_cluster_ncount_init (NcABCClusterNCount *abcnc)
{
  abcnc->data_summary   = NULL;
  abcnc->ncount         = NULL;
  abcnc->data_total     = 0.0;
  abcnc->covar          = NULL;
  abcnc->scale_cov      = FALSE;
  abcnc->quantiles      = NULL;
  abcnc->z_nodes        = NULL;
  abcnc->lnM_nodes      = NULL;
  abcnc->z_bins         = 0;
  abcnc->lnM_bins       = 0;
  abcnc->epsilon_update = 0.0;
  abcnc->bin_type       = NC_ABC_CLUSTER_NCOUNT_BIN_NTYPES;
  abcnc->uptype         = NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_NTYPE;
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
    case PROP_BIN_TYPE:
    {
      NcABCClusterNCountBin bin_type = g_value_get_enum (value);
      switch (bin_type)
      {
        case NC_ABC_CLUSTER_NCOUNT_BIN_UNIFORM:
          break;
        case NC_ABC_CLUSTER_NCOUNT_BIN_QUANTILE:
          g_assert (abcnc->quantiles != NULL);
          break;
        case NC_ABC_CLUSTER_NCOUNT_BIN_NODES:
          g_assert (abcnc->z_nodes != NULL);
          g_assert (abcnc->lnM_nodes != NULL);
          break;
        default:
          g_assert_not_reached ();
          break;
      }
      abcnc->bin_type = bin_type;
      break;
    }
    case PROP_QUANTILES:
    {
      GVariant *var = g_value_get_variant (value);
      if (var != NULL)
      {
        NcmVector *v = ncm_vector_new_variant (var);
        ncm_vector_clear (&abcnc->quantiles);
        abcnc->quantiles = v;
      }
      break;
    }
    case PROP_Z_NODES:
    {
      GVariant *var = g_value_get_variant (value);
      if (var != NULL)
      {
        NcmVector *v = ncm_vector_new_variant (var);
        ncm_vector_clear (&abcnc->z_nodes);
        abcnc->z_nodes = v;
      }
      break;
    }
    case PROP_LNM_NODES:
    {
      GVariant *var = g_value_get_variant (value);
      if (var != NULL)
      {
        NcmVector *v = ncm_vector_new_variant (var);
        ncm_vector_clear (&abcnc->lnM_nodes);
        abcnc->lnM_nodes = v;
      }
      break;
    }
    case PROP_Z_BINS:
      abcnc->z_bins = g_value_get_uint (value);
      g_assert_cmpuint (abcnc->z_bins, >, 0);
      break;
    case PROP_LNM_BINS:
      abcnc->lnM_bins = g_value_get_uint (value);
      g_assert_cmpuint (abcnc->lnM_bins, >, 0);
      break;
    case PROP_EPSILON_UPDATE:
      nc_abc_cluster_ncount_set_epsilon_update (abcnc, g_value_get_double (value));
      break;
    case PROP_EPSILON_UPDATE_TYPE:
      abcnc->uptype = g_value_get_enum (value);
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
    case PROP_SCALE_COV:
      g_value_set_boolean (value, abcnc->scale_cov);
      break;
    case PROP_BIN_TYPE:
      g_value_set_enum (value, abcnc->bin_type);
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
    case PROP_Z_NODES:
    {
      if (abcnc->quantiles != NULL)
      {
        GVariant *var = ncm_vector_peek_variant (abcnc->z_nodes); 
        g_value_take_variant (value, var);
      }
      break;
    }
    case PROP_LNM_NODES:
    {
      if (abcnc->quantiles != NULL)
      {
        GVariant *var = ncm_vector_peek_variant (abcnc->lnM_nodes); 
        g_value_take_variant (value, var);
      }
      break;
    }
    case PROP_Z_BINS:
      g_value_set_uint (value, abcnc->z_bins);
      break;
    case PROP_LNM_BINS:
      g_value_set_uint (value, abcnc->lnM_bins);
      break;
    case PROP_EPSILON_UPDATE:
      g_value_set_double (value, abcnc->epsilon_update);
      break;
    case PROP_EPSILON_UPDATE_TYPE:
      g_value_set_enum (value, abcnc->uptype);
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
  ncm_vector_clear (&abcnc->quantiles);
  ncm_vector_clear (&abcnc->z_nodes);
  ncm_vector_clear (&abcnc->lnM_nodes);

  nc_data_cluster_ncount_clear (&abcnc->ncount);
  
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

  g_object_class_install_property (object_class,
                                   PROP_BIN_TYPE,
                                   g_param_spec_enum ("binning-type",
                                                      NULL,
                                                      "Binning type",
                                                      NC_TYPE_ABC_CLUSTER_NCOUNT_BIN, NC_ABC_CLUSTER_NCOUNT_BIN_UNIFORM,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); 

  g_object_class_install_property (object_class,
                                   PROP_QUANTILES,
                                   g_param_spec_variant ("quantiles",
                                                         NULL,
                                                         "Quantiles for binning",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_Z_NODES,
                                   g_param_spec_variant ("z-nodes",
                                                         NULL,
                                                         "Nodes for z",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_LNM_NODES,
                                   g_param_spec_variant ("lnM-nodes",
                                                         NULL,
                                                         "Nodes for lnM",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z_BINS,
                                   g_param_spec_uint ("z-bins",
                                                      NULL,
                                                      "Number of bins in z",
                                                      1, G_MAXUINT32, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_LNM_BINS,
                                   g_param_spec_uint ("lnM-bins",
                                                      NULL,
                                                      "Number of bins in lnM",
                                                      1, G_MAXUINT32, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_EPSILON_UPDATE,
                                   g_param_spec_double ("epsilon-update",
                                                      NULL,
                                                      "Value used to update epsilon",
                                                      0.0, 1.0, 0.75,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_EPSILON_UPDATE_TYPE,
                                   g_param_spec_enum ("epsilon-update-type",
                                                      NULL,
                                                      "Method used to update epsilon",
                                                      NC_TYPE_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE, NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_QUANTILE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  abc_class->data_summary  = &_nc_abc_cluster_ncount_data_summary;
  abc_class->mock_distance = &_nc_abc_cluster_ncount_mock_distance;
  abc_class->distance_prob = &_nc_abc_cluster_ncount_distance_prob;
  abc_class->update_tkern  = &_nc_abc_cluster_ncount_update_tkern;
  abc_class->get_desc      = &_nc_abc_cluster_ncount_get_desc;
  abc_class->log_info      = &_nc_abc_cluster_ncount_log_info;
}

static gdouble
_nc_abc_summary_smooth (NcABCClusterNCount *abcnc, NcDataClusterNCount *ncount, gdouble z, gdouble lnM)
{
  const gsize nx = gsl_histogram2d_nx (abcnc->data_summary);
  const gsize ny = gsl_histogram2d_ny (abcnc->data_summary);
  const gdouble dz      = (abcnc->data_summary->xrange[nx - 1] - abcnc->data_summary->xrange[0]);
  const gdouble dlnM    = (abcnc->data_summary->yrange[ny - 1] - abcnc->data_summary->yrange[0]);
  const gdouble sigma_z   = dz * 0.05;
  const gdouble sigma_lnM = dlnM * 0.05;
  const gdouble logsqrt2pi_sigma_z = log (sqrt (2.0 * M_PI) * sigma_z);
  const gdouble logsqrt2pi_sigma_lnM = log (sqrt (2.0 * M_PI) * sigma_lnM);
  
  gdouble res = 0.0;
  /*gsize i, j;*/
  gsize i;

  for (i = 0; i < ncount->np; i++)
  {
    const gdouble zi = ncm_matrix_get (ncount->z_obs, i, 0);
    const gdouble lnMi = ncm_matrix_get (ncount->lnM_obs, i, 0);
    const gdouble x2      = -pow ((z - zi) / sigma_z, 2.0) * 0.5 - logsqrt2pi_sigma_z;
    const gdouble y2      = -pow ((lnM - lnMi) / sigma_lnM, 2.0) * 0.5 - logsqrt2pi_sigma_lnM;
    res += exp (x2 + y2);
//printf ("# dist % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", x2 + y2, exp (x2 + y2), nij, z, lnM, z_bin, lnM_bin);
  }
/*
  for (i = 0; i < nx - 1; i++)
  {
    for (j = 0; j < ny - 1; j++)
    {
      const gdouble nij = abcnc->data_summary->bin[i * ny + j];
      if (nij > 0)
      {
        const gdouble z_bin   = (abcnc->data_summary->xrange[i + 1] + abcnc->data_summary->xrange[i]) * 0.5;
        const gdouble lnM_bin = (abcnc->data_summary->yrange[j + 1] + abcnc->data_summary->yrange[j]) * 0.5;
        const gdouble x2      = -pow ((z - z_bin) / sigma_z, 2.0) * 0.5 - logsqrt2pi_sigma_z;
        const gdouble y2      = -pow ((lnM - lnM_bin) / sigma_lnM, 2.0) * 0.5 - logsqrt2pi_sigma_lnM;
//printf ("# dist % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", x2 + y2, exp (x2 + y2), nij, z, lnM, z_bin, lnM_bin);
        res += exp (x2 + y2) * nij;
      }
    }
  }
*/
  return res;
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
    g_assert (NC_IS_DATA_CLUSTER_NCOUNT (data));

    switch (abcnc->bin_type)
    {
      case NC_ABC_CLUSTER_NCOUNT_BIN_UNIFORM:
      {
        nc_data_cluster_ncount_set_bin_by_minmax (ncount, abcnc->z_bins, abcnc->lnM_bins);    
        break;
      }
      case NC_ABC_CLUSTER_NCOUNT_BIN_QUANTILE:
      {
        nc_data_cluster_ncount_set_bin_by_quantile (ncount, abcnc->quantiles, abcnc->quantiles);
        break;
      }
      case NC_ABC_CLUSTER_NCOUNT_BIN_NODES:
      {
        nc_data_cluster_ncount_set_bin_by_nodes (ncount, abcnc->z_nodes, abcnc->lnM_nodes);
        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    g_clear_pointer (&abcnc->data_summary, gsl_histogram2d_free);

    nc_data_cluster_ncount_clear (&abcnc->ncount);
    abcnc->ncount = nc_data_cluster_ncount_ref (ncount);
    abcnc->data_summary = gsl_histogram2d_clone (ncount->z_lnM);
    abcnc->data_total   = gsl_histogram2d_sum (abcnc->data_summary);

    if (TRUE)
    {
      gsize i;
      gdouble res = 0.0;
      gdouble te = 0.0;

      for (i = 0; i < abcnc->ncount->np; i++)
      {
        const gdouble zi = ncm_matrix_get (abcnc->ncount->z_obs, i, 0);
        const gdouble lnMi = ncm_matrix_get (abcnc->ncount->lnM_obs, i, 0);
        const gdouble smooth  = _nc_abc_summary_smooth (abcnc, abcnc->ncount, zi, lnMi);
        res += -log (smooth);
        te++;
      }

      abcnc->data_total = 2.0 * (res + te);
      /*printf ("# Real data gives % 20.15g % 20.15g % 20.15g\n", res, te, abcnc->data_total);*/
    }
  }
  
  return TRUE;
}

static gdouble 
_nc_abc_cluster_ncount_mock_distance (NcmABC *abc, NcmDataset *dset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcABCClusterNCount *abcnc = NC_ABC_CLUSTER_NCOUNT (abc);
  NcmData *data = ncm_dataset_peek_data (dset, 0);
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);

  gsize i; 
  gdouble res = 0.0, te = 0.0;

  for (i = 0; i < abcnc->ncount->np; i++)
  {
    const gdouble zi = ncm_matrix_get (abcnc->ncount->z_obs, i, 0);
    const gdouble lnMi = ncm_matrix_get (abcnc->ncount->lnM_obs, i, 0);
    const gdouble smooth  = _nc_abc_summary_smooth (abcnc, ncount, zi, lnMi);
    res += -log (smooth);
  }
  te = ncount->np;
  /*
  for (i = 0; i < nx - 1; i++)
  {
    for (j = 0; j < ny - 1; j++)
    {
      const gdouble nij = ncount->z_lnM->bin[i * ny + j];
      if (nij > 0)
      {
        const gdouble z_bin   = (abcnc->data_summary->xrange[i + 1] + abcnc->data_summary->xrange[i]) * 0.5;
        const gdouble lnM_bin = (abcnc->data_summary->yrange[j + 1] + abcnc->data_summary->yrange[j]) * 0.5;
        const gdouble smooth  = _nc_abc_summary_smooth (abcnc, z_bin, lnM_bin);
        
        res += -log (nij * smooth);
        te  += nij;
        // printf ("% zd % zd % 20.15g % 20.15g % 20.15g % 20.15g\n", i, j, res, te, nij, smooth);
      }
    }
  }
*/
  /*printf ("try: % 20.15g % 20.15g % 20.15g\n", res, te, (2.0 * (res + te) - abcnc->data_total) / abcnc->ncount->np);*/
  return (2.0 * (res + te) - abcnc->data_total) / abcnc->ncount->np;

  /*
  gdouble pdf_data = 0.0;
  gdouble pdf_mock = 0.0;
  gdouble nebins = 0.0;
  gint e = 0;

  for (i = 0; i < total; i++)
  {
    pdf_data += abcnc->data_summary->bin[i];
    pdf_mock += ncount->z_lnM->bin[i];
    
    res += gsl_pow_2 ((pdf_data - pdf_mock) / abcnc->data_total); 
  }
  */
/*
  for (i = 0; i < total; i++)
  {
    pdf_data += abcnc->data_summary->bin[i];
    pdf_mock += ncount->z_lnM->bin[i];

    if (abcnc->data_total - pdf_data <= 5.0)
      break;
    
    if (pdf_data >= 5.0)
    {
      nebins++;
      res += -2.0 * ((pdf_mock - pdf_data) * log (pdf_data) + lgamma_r (pdf_data + 1.0, &e) - lgamma_r (pdf_mock + 1, &e));
      pdf_data = 0.0;
      pdf_mock = 0.0;
    }
  }

  for (; i < total; i++)
  {
    pdf_data += abcnc->data_summary->bin[i];
    pdf_mock += ncount->z_lnM->bin[i];
  }
  if (pdf_data != 0.0 || pdf_mock != 0.0)
  {
    nebins++;
    res += -2.0 * ((pdf_mock - pdf_data) * log (pdf_data) + lgamma_r (pdf_data + 1.0, &e) - lgamma_r (pdf_mock + 1, &e));
    pdf_data = 0.0;
    pdf_mock = 0.0;    
  }
  
  return res / nebins;
*/
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
  const gdouble epsilon = NCM_ABC (abcnc)->epsilon;
  
  ncm_mset_catalog_get_covar (abc->mcat, &abc->covar);
  
  ncm_matrix_scale (abc->covar, scale);
  ncm_mset_trans_kern_gauss_set_cov (NCM_MSET_TRANS_KERN_GAUSS (abc->tkern), abc->covar);

  if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("# NcABCClusterNCount: scale covariance by %f\n", scale);

  switch (abcnc->uptype)
  {
    case NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_UNIFORM:
    {
      if (epsilon - abcnc->epsilon_update < 0)
        ncm_abc_update_epsilon (abc, epsilon * abcnc->epsilon_update);
      else
        ncm_abc_update_epsilon (abc, epsilon - abcnc->epsilon_update);
      break;
    }
    case NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_QUANTILE:
    {
      gdouble dist = ncm_abc_get_dist_quantile (abc, abcnc->epsilon_update);
      ncm_abc_update_epsilon (abc, dist);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
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

/**
 * nc_abc_cluster_ncount_set_epsilon_update:
 * @abcnc: a #NcABCClusterNCount.
 * @q: the quantile $q \in (0, 1)$.
 * 
 * Sets the quantile used to update epsilon.
 * 
 */
void 
nc_abc_cluster_ncount_set_epsilon_update (NcABCClusterNCount *abcnc, gdouble q)
{
  g_assert_cmpfloat (q, <, 1.0);
  g_assert_cmpfloat (q, >, 0.0);
  abcnc->epsilon_update = q;
}

/**
 * nc_abc_cluster_ncount_set_bin_uniform:
 * @abcnc: a #NcABCClusterNCount.
 * @z_bins: number of bins in z.
 * @lnM_bins: number of bins in lnM.
 * 
 * Sets the binning type to #NC_ABC_CLUSTER_NCOUNT_BIN_UNIFORM.
 * 
 */
void 
nc_abc_cluster_ncount_set_bin_uniform (NcABCClusterNCount *abcnc, guint z_bins, guint lnM_bins)
{
  g_assert_cmpuint (z_bins, >, 0);
  g_assert_cmpuint (lnM_bins, >, 0);
  abcnc->z_bins = z_bins;
  abcnc->lnM_bins = lnM_bins;
  abcnc->bin_type = NC_ABC_CLUSTER_NCOUNT_BIN_UNIFORM;
}

/**
 * nc_abc_cluster_ncount_set_bin_quantile:
 * @abcnc: a #NcABCClusterNCount.
 * @quantiles: (allow-none): a #NcmVector or NULL.
 * 
 * Sets the binning type to #NC_ABC_CLUSTER_NCOUNT_BIN_QUANTILE and uses
 * @quantiles as the quantiles for both z and lnM. If @quantiles is NULL
 * uses the defaults: (0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98).
 * 
 */
void 
nc_abc_cluster_ncount_set_bin_quantile (NcABCClusterNCount *abcnc, NcmVector *quantiles)
{
  if (quantiles == NULL)
  {
    gdouble quantiles_data[7] = {0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98};
    quantiles = ncm_vector_new_data_dup (quantiles_data, 7, 1);
  }
  else
    ncm_vector_ref (quantiles);

  abcnc->bin_type = NC_ABC_CLUSTER_NCOUNT_BIN_QUANTILE;
  abcnc->quantiles = quantiles;
}

/**
 * nc_abc_cluster_ncount_set_bin_nodes:
 * @abcnc: a #NcABCClusterNCount.
 * @z_nodes: a #NcmVector.
 * @lnM_nodes: a #NcmVector.
 * 
 * Sets the binning type to #NC_ABC_CLUSTER_NCOUNT_BIN_NODES and uses
 * @z_nodes and @lnM_nodes as nodes for binning.
 * 
 */
void 
nc_abc_cluster_ncount_set_bin_nodes (NcABCClusterNCount *abcnc, NcmVector *z_nodes, NcmVector *lnM_nodes)
{
  g_assert (z_nodes != NULL);
  g_assert (lnM_nodes != NULL);

  abcnc->bin_type  = NC_ABC_CLUSTER_NCOUNT_BIN_NODES;
  abcnc->z_nodes   = ncm_vector_ref (z_nodes);
  abcnc->lnM_nodes = ncm_vector_ref (lnM_nodes);
}
