/***************************************************************************
 *            nc_data_cluster_pseudo_counts.c
 *
 *  Sun Apr 5 20:23:19 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_cluster_pseudo_counts.c
 * Copyright (C) 2015 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_cluster_pseudo_counts
 * @title: NcDataClusterPseudoCounts
 * @short_description: Galaxy clusters data -- pseudo number counts likelihood.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_pseudo_counts.h"
#include "nc_hicosmo.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_OBS, 
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataClusterPseudoCounts, nc_data_cluster_pseudo_counts, NCM_TYPE_DATA);

static void
nc_data_cluster_pseudo_counts_init (NcDataClusterPseudoCounts *dcpc)
{
  dcpc->np  = 0;
  dcpc->obs = NULL;
}

static void
nc_data_cluster_pseudo_counts_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_PSEUDO_COUNTS (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_data_cluster_pseudo_counts_set_obs (dcpc, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_pseudo_counts_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_PSEUDO_COUNTS (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, dcpc->obs);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_pseudo_counts_dispose (GObject *object)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);

  ncm_matrix_clear (&dcpc->obs);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_pseudo_counts_parent_class)->dispose (object);
}

static void
nc_data_cluster_pseudo_counts_finalize (GObject *object)
{
  /*NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_pseudo_counts_parent_class)->finalize (object);
}

static void _nc_data_cluster_pseudo_counts_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_cluster_pseudo_counts_get_len (NcmData *data) { return NC_DATA_CLUSTER_PSEUDO_COUNTS (data)->np; }

static void
nc_data_cluster_pseudo_counts_class_init (NcDataClusterPseudoCountsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_pseudo_counts_set_property;
  object_class->get_property = nc_data_cluster_pseudo_counts_get_property;
  object_class->dispose      = nc_data_cluster_pseudo_counts_dispose;
  object_class->finalize     = nc_data_cluster_pseudo_counts_finalize;

  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                         NULL,
                                                         "Cluster observables",
                                                         NCM_TYPE_MATRIX,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->m2lnL_val  = &_nc_data_cluster_pseudo_counts_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_pseudo_counts_get_len;
}

/**
 * nc_data_cluster_pseudo_counts_new:
 * 
 * Creates a new #NcDataClusterPseudoCounts.
 * 
 * Returns: the newly created #NcDataClusterPseudoCounts.
 */
NcDataClusterPseudoCounts *
nc_data_cluster_pseudo_counts_new (void)
{
  NcDataClusterPseudoCounts *dcpc;

  dcpc = g_object_new (NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS,
                       NULL);

  return dcpc;
}

/**
 * nc_data_cluster_pseudo_counts_new_from_file:
 * @filename: file containing a serialized #NcDataClusterPseudoCounts
 * 
 * Creates a new #NcDataClusterPseudoCounts from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataClusterPseudoCounts.
 */
NcDataClusterPseudoCounts *
nc_data_cluster_pseudo_counts_new_from_file (const gchar *filename)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_CLUSTER_PSEUDO_COUNTS (dcpc));

  return dcpc;
}

/**
 * nc_data_cluster_pseudo_counts_ref:
 * @dcpc: a #NcDataClusterPseudoCounts
 *
 * Increases the reference count of @dcpc by one.
 * 
 * Returns: (transfer full): @dcpc
 */
NcDataClusterPseudoCounts *
nc_data_cluster_pseudo_counts_ref (NcDataClusterPseudoCounts *dcpc)
{
  return g_object_ref (dcpc);
}

/**
 * nc_data_cluster_pseudo_counts_free:
 * @dcpc: a #NcDataClusterPseudoCounts
 *
 * Atomically decrements the reference count of @dcpc by one. If the reference count drops to 0,
 * all memory allocated by @dcpc is released.
 * 
 */
void
nc_data_cluster_pseudo_counts_free (NcDataClusterPseudoCounts *dcpc)
{
  g_object_unref (dcpc);
}

/**
 * nc_data_cluster_pseudo_counts_clear:
 * @dcpc: a #NcDataClusterPseudoCounts
 *
 * The reference count of @dcpc is decreased and the pointer is set to NULL.
 * 
 */
void
nc_data_cluster_pseudo_counts_clear (NcDataClusterPseudoCounts **dcpc)
{
  g_clear_object (dcpc);
}

/**
 * nc_data_cluster_pseudo_counts_set_obs:
 * @dcpc: a #NcDataClusterPseudoCounts
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster mass observables, e.g., redshift 
 * observed mass(es) and the parameters of the redshift and/or mass-observable distributions.
 * 
 */
void 
nc_data_cluster_pseudo_counts_set_obs (NcDataClusterPseudoCounts *dcpc, const NcmMatrix *m)
{
  guint np;
  g_assert (m != NULL);
  g_assert_cmpuint (ncm_matrix_ncols (m), ==, NC_DATA_CLUSTER_PSEUDO_COUNTS_LEN);
  np = ncm_matrix_nrows (m);

  if (dcpc->np > 0)
  {
    if (np != dcpc->np)
    {
      ncm_matrix_clear (&dcpc->obs);
      dcpc->obs = ncm_matrix_new (np, NC_DATA_CLUSTER_PSEUDO_COUNTS_LEN);
    }
    ncm_matrix_memcpy (dcpc->obs, m);
  }
  else
  {
    dcpc->np = np;
    dcpc->obs = ncm_matrix_dup (m);
  }

  ncm_data_set_init (NCM_DATA (dcpc), TRUE);
}

static void
_nc_data_cluster_pseudo_counts_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterMass *clusterm = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterPseudoCounts *cpc = NC_CLUSTER_PSEUDO_COUNTS (ncm_mset_peek (mset, nc_cluster_pseudo_counts_id ()));

  g_assert (cosmo != NULL);
  g_assert (clusterm != NULL);
  g_assert (cpc != NULL);
  
  gdouble Ndet = nc_cluster_pseudo_counts_ndet (cpc, cosmo);
  gdouble lnNdet = log (Ndet);
  gint i;
  *m2lnL = 0.0;
  if (Ndet < 1.0)
  {
    *m2lnL = GSL_POSINF;
    return;
  }
  else
    *m2lnL = 0.0;

  //printf ("# Ndet % 20.15g\n", Ndet);
  for (i = 0; i < dcpc->np; i++)
  {
    const gdouble z = ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_Z);
    const gdouble *M = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL);
    const gdouble *M_params = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL);
    //const gdouble m2lnL_i = log (nc_cluster_pseudo_counts_posterior_numerator (cpc, clusterm, cosmo, z, M, M_params));
    const gdouble m2lnL_i = log (nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, clusterm, cosmo, z, M[0], M[1], M_params[0], M_params[1]));

    //printf ("%d % 20.15g % 20.15g % 20.15g % 20.15g\n", i, m2lnL_i, m2lnL_i - log (Ndet), log (Ndet), Ndet);
    if (!gsl_finite (m2lnL_i))
    {
      *m2lnL += m2lnL_i;
      break;
    }
    else
    {
      *m2lnL += m2lnL_i;
    }
  }

  *m2lnL -= dcpc->np * lnNdet;

  *m2lnL = -2.0 * (*m2lnL);
  return;
}
