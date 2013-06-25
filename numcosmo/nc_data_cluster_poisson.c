/***************************************************************************
 *            nc_data_cluster_poisson.c
 *
 *  Tue Apr  6 01:11:23 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
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
 * SECTION:nc_data_cluster_poisson
 * @title: Cluster number count data
 * @short_description: FIXME
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_cluster_poisson.h"
#include "nc_data_cluster_ncount.h"
#include "math/util.h"

#include <gsl/gsl_randist.h>

G_DEFINE_TYPE (NcDataClusterPoisson, nc_data_cluster_poisson, NCM_TYPE_DATA_POISSON);

static void
nc_data_cluster_poisson_init (NcDataClusterPoisson *cluster_poisson)
{
  cluster_poisson->ncount = NULL;
}

static void
nc_data_cluster_poisson_dispose (GObject *object)
{
  NcDataClusterPoisson *cluster_poisson = NC_DATA_CLUSTER_POISSON (object);
  NcmData *data = NCM_DATA (cluster_poisson->ncount);
  
  ncm_data_clear (&data);
  cluster_poisson->ncount = NULL;
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_poisson_parent_class)->dispose (object);
}

static void
nc_data_cluster_poisson_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_poisson_parent_class)->finalize (object);
}

static void _nc_data_cluster_poisson_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_poisson_resample (NcmData *data, NcmMSet *mset);

static void
nc_data_cluster_poisson_class_init (NcDataClusterPoissonClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->dispose  = nc_data_cluster_poisson_dispose;
  object_class->finalize = nc_data_cluster_poisson_finalize;

  data_class->prepare    = &_nc_data_cluster_poisson_prepare;
  data_class->resample   = &_nc_data_cluster_poisson_resample;
}

static void
_nc_data_cluster_poisson_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterPoisson *cpoisson = NC_DATA_CLUSTER_POISSON (data);
  ncm_data_prepare (NCM_DATA (cpoisson->ncount), mset);
}

static void
_nc_data_cluster_poisson_resample (NcmData *data, NcmMSet *mset)
{
  NcDataClusterPoisson *cpoisson = NC_DATA_CLUSTER_POISSON (data);
  NcDataClusterNCount *ncount = cpoisson->ncount;
  NcClusterAbundance *cad = ncount->cad;
  NcmDataPoisson *poisson = NCM_DATA_POISSON (cpoisson);
  gint i;

  ncm_data_resample (NCM_DATA (cpoisson->ncount), mset);
  
  nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())));
  gsl_histogram_reset (poisson->h);

  for (i = 0; i < ncount->np; i++)
  {
    g_assert_not_reached ();
    /* FIXME */
  }
}


/**
 * nc_data_cluster_poisson_new:
 * @ncount: a #NcClusterAbundance.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_poisson_new (NcDataClusterNCount *ncount)
{
  NcDataClusterPoisson *cpoisson = g_object_new (NC_TYPE_DATA_CLUSTER_POISSON, 
                                                 "data-cluster-ncount", ncount,
                                                 NULL);  

  return NCM_DATA (cpoisson);
}

/**
 * nc_data_cluster_poisson_new_cad:
 * @cad: a #NcClusterAbundance.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_poisson_new_cad (NcClusterAbundance *cad)
{
  NcmData *data_ncount = nc_data_cluster_ncount_new (cad);
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data_ncount);
  NcmData *data = nc_data_cluster_poisson_new (ncount);  

  ncm_data_free (data_ncount);
  return data;
}

