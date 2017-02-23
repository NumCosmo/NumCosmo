/***************************************************************************
 *            nc_data_cluster_counts_box_poisson.c
 *
 *  Mon Feb 20 15:31:47 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_cluster_counts_box_poisson
 * @title: NcDataClusterCountsBoxPoisson
 * @short_description: Binned (mass) cluster number count data in a box.
 *
 * FIXME Box (simulation), and not redshift space.
 * 
 */

enum
{
  PROP_0,
  PROP_NBIN_MASS,
  PROP_VOLUME,
  PROP_SIZE,
};

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_counts_box_poisson.h"
#include "lss/nc_halo_mass_function.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcDataClusterCountsBoxPoisson, nc_data_cluster_counts_box_poisson, NCM_TYPE_DATA_POISSON);

static void
nc_data_cluster_counts_box_poisson_init (NcDataClusterCountsBoxPoisson *cluster_poisson)
{
  cluster_poisson->mfp = NULL;
}

static void
nc_data_cluster_counts_box_poisson_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterCountsBoxPoisson *poisson = C_DATA_CLUSTER_COUNTS_BOX_POISSON (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_POISSON (object));

  switch (prop_id)
  {
    case PROP_NBIN_MASS:
      g_value_set_int (&poisson->nbin_mass);
      poisson->nbin_mass = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_counts_box_poisson_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterCountsBoxPoisson *poisson = C_DATA_CLUSTER_COUNTS_BOX_POISSON (object);

  g_return_if_fail (NC_IS_DATA_CLUSTER_POISSON (object));

  switch (prop_id)
  {
    case PROP_NBIN_MASS:
	  g_value_set_int (value, poisson->nbin_mass);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_counts_box_poisson_dispose (GObject *object)
{
  NcDataClusterCountsBoxPoisson *cluster_poisson = C_DATA_CLUSTER_COUNTS_BOX_POISSON (object);
  
  nc_halo_mass_function_clear (&cluster_poisson->mfp);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_counts_box_poisson_parent_class)->dispose (object);
}

static void
nc_data_cluster_counts_box_poisson_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_counts_box_poisson_parent_class)->finalize (object);
}

static void _nc_data_cluster_counts_box_poisson_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_counts_box_poisson_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);

static void
nc_data_cluster_counts_box_poisson_class_init (NcDataClusterCountsBoxPoissonClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &nc_data_cluster_counts_box_poisson_set_property;
  object_class->get_property = &nc_data_cluster_counts_box_poisson_get_property;
  object_class->dispose      = &nc_data_cluster_counts_box_poisson_dispose;
  object_class->finalize     = &nc_data_cluster_counts_box_poisson_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NBIN_MASS,
                                   g_param_spec_object ("nbin-mass",
                                                        NULL,
                                                        "Number of mass bins",
                                                        0, 1000, 20,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->prepare    = &_nc_data_cluster_counts_box_poisson_prepare;
  data_class->resample   = &_nc_data_cluster_counts_box_poisson_resample;
}

static void
_nc_data_cluster_counts_box_poisson_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterCountsBoxPoisson *cpoisson = C_DATA_CLUSTER_COUNTS_BOX_POISSON (data);
  ncm_data_prepare (NCM_DATA (cpoisson->ncount), mset);
}

static void
_nc_data_cluster_counts_box_poisson_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterCountsBoxPoisson *cpoisson = C_DATA_CLUSTER_COUNTS_BOX_POISSON (data);
  NcDataClusterNCount *ncount = cpoisson->ncount;
  NcClusterAbundance *cad = ncount->cad;
  NcmDataPoisson *poisson = NCM_DATA_POISSON (cpoisson);
  guint i;

  ncm_data_resample (NCM_DATA (cpoisson->ncount), mset, rng);
  
  nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())),
                                         cad->lnMi);

  gsl_histogram_reset (poisson->h);

  for (i = 0; i < ncount->np; i++)
  {
    g_assert_not_reached ();
    /* FIXME */
  }
}


/**
 * nc_data_cluster_counts_box_poisson_new:
 * @ncount: a #NcClusterAbundance.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_counts_box_poisson_new (NcDataClusterNCount *ncount)
{
  NcDataClusterCountsBoxPoisson *cpoisson = g_object_new (NC_TYPE_DATA_CLUSTER_POISSON, 
                                                 "data-cluster-ncount", ncount,
                                                 NULL);  

  return NCM_DATA (cpoisson);
}

/**
 * nc_data_cluster_counts_box_poisson_new_cad:
 * @cad: a #NcClusterAbundance.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_counts_box_poisson_new_cad (NcClusterAbundance *cad)
{
  NcDataClusterNCount *ncount = nc_data_cluster_ncount_new (cad);
  NcmData *data = nc_data_cluster_counts_box_poisson_new (ncount);  

  ncm_data_free (NCM_DATA (ncount));
  return data;
}

