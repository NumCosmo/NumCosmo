/***************************************************************************
 *            data_cluster.c
 *
 *  Sat Apr 24 14:29:17 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include "de_options.h"
#include "data_cluster.h"
#include "savedata.h"

static void _nc_de_data_cluster_append (NcDEDataClusterEntries *de_data_cluster, NcData *dca_unbinned, NcDataSet *ds);

GPtrArray *
nc_de_data_cluster_new (NcDistance *dist, NcmMSet *mset, NcDEDataClusterEntries *de_data_cluster, NcDataSet *ds, gint id)
{
  GPtrArray *ca_array = g_ptr_array_new ();
  NcWindow *wp = nc_window_new_from_name (de_data_cluster->window_name);
  NcTransferFunc *tf = nc_transfer_func_new_from_name (de_data_cluster->transfer_name);
  NcMatterVar *vp = nc_matter_var_new (NC_MATTER_VAR_FFT, wp, tf);
  NcGrowthFunc *gf = nc_growth_func_new ();
  NcMultiplicityFunc *mulf = nc_multiplicity_func_new_from_name (de_data_cluster->multiplicity_name);
  NcMassFunction *mfp = nc_mass_function_new (dist, vp, gf, mulf);
  NcClusterMass *clusterm = nc_cluster_mass_new_from_name (de_data_cluster->clusterm_ser);
  NcClusterRedshift *clusterz = nc_cluster_redshift_new_from_name (de_data_cluster->clusterz_ser);
  NcClusterAbundanceOpt opt =
	(de_data_cluster->use_photoz ? NC_CLUSTER_ABUNDANCE_PHOTOZ  : NC_CLUSTER_ABUNDANCE_NONE) |
	(de_data_cluster->use_Mobs   ? NC_CLUSTER_ABUNDANCE_MOBS    : NC_CLUSTER_ABUNDANCE_NONE) |
	(de_data_cluster->use_Mobs_local   ? NC_CLUSTER_ABUNDANCE_MOBS_LOCAL    : NC_CLUSTER_ABUNDANCE_NONE) |
	(de_data_cluster->use_selection   ? (NC_CLUSTER_ABUNDANCE_COMPLETENESS | NC_CLUSTER_ABUNDANCE_PURITY)  : NC_CLUSTER_ABUNDANCE_NONE) |
	(de_data_cluster->binmass    ? NC_CLUSTER_ABUNDANCE_BINMASS : NC_CLUSTER_ABUNDANCE_NONE);

  switch (id)
  {
    case 0:
    {
      gint i = 0;
      if (de_data_cluster->cata_file == NULL)
        g_error ("For --cluster-id 0, you must specify a fit catalog via --catalog file.fit");
      while (de_data_cluster->cata_file[i] != NULL)
      {
        NcClusterAbundance *cad = nc_cluster_abundance_new (opt, mfp, NULL, clusterz, clusterm);
        NcData *dca_unbinned = nc_data_cluster_abundance_unbinned_new (cad);
        nc_cluster_abundance_catalog_load (dca_unbinned, de_data_cluster->cata_file[i], opt);
        _nc_de_data_cluster_append (de_data_cluster, dca_unbinned, ds);
				g_ptr_array_add (ca_array, dca_unbinned);
        if ((i == 0) && (de_data_cluster->save_cata != NULL))
          nc_cluster_abundance_catalog_save (dca_unbinned, de_data_cluster->save_cata, TRUE);
        i++;
			}
    }
      break;
    case 1:
    {
      gint i = 0;
      if (de_data_cluster->cata_file == NULL)
        g_error ("For --cluster-id 1, you must specify a fit catalog via --catalog filename");
      while (de_data_cluster->cata_file[i] != NULL)
      {
        NcClusterAbundance *cad = nc_cluster_abundance_new (opt, mfp, NULL, clusterz, clusterm);
        NcData *dca_unbinned = nc_data_cluster_abundance_unbinned_new (cad);
				nc_data_cluster_abundance_unbinned_init_from_text_file (dca_unbinned, de_data_cluster->cata_file[i], opt, de_data_cluster->area_survey * gsl_pow_2 (M_PI / 180.0), log(de_data_cluster->Mi), log(de_data_cluster->Mf), de_data_cluster->z_initial, de_data_cluster->z_final, de_data_cluster->photoz_sigma0, de_data_cluster->photoz_bias, de_data_cluster->lnM_sigma0, de_data_cluster->lnM_bias);
        _nc_de_data_cluster_append (de_data_cluster, dca_unbinned, ds);
				g_ptr_array_add (ca_array, dca_unbinned);
        if ((i == 0) && (de_data_cluster->save_cata != NULL))
          nc_cluster_abundance_catalog_save (dca_unbinned, de_data_cluster->save_cata, TRUE);
        i++;
      }
    }
      break;
    case 2:
    {
      NcClusterAbundance *cad = nc_cluster_abundance_new (opt, mfp, NULL, clusterz, clusterm);
      NcData *dca_unbinned = nc_data_cluster_abundance_unbinned_new (cad);

      nc_data_cluster_abundance_unbinned_init_from_sampling (dca_unbinned, mset, clusterz, clusterm, de_data_cluster->area_survey * gsl_pow_2 (M_PI / 180.0));

	  if (de_data_cluster->save_cata != NULL)
        nc_cluster_abundance_catalog_save (dca_unbinned, de_data_cluster->save_cata, TRUE);
      _nc_de_data_cluster_append (de_data_cluster, dca_unbinned, ds);
			g_ptr_array_add (ca_array, dca_unbinned);
    }
      break;
    default:
      g_error ("The option --catalog-id must be between (0,2).");
  }
	return ca_array;
}

static void
_nc_de_data_cluster_append (NcDEDataClusterEntries *de_data_cluster, NcData *dca_unbinned, NcDataSet *ds)
{
  if (de_data_cluster->binned)
  {
    NcData *dca_binned;
    gsl_vector *nodes = gsl_vector_alloc (de_data_cluster->n_bins + 1);
    gint i;
    for (i = 0; i <= de_data_cluster->n_bins; i++)
      gsl_vector_set (nodes, i, de_data_cluster->z_initial + (de_data_cluster->z_final - de_data_cluster->z_initial) / (de_data_cluster->n_bins) * i);

    dca_binned = nc_data_cluster_abundance_unbinned_bin_data (dca_unbinned, nodes);

    /* nc_data_free0 (dca_unbinned, TRUE); */
    gsl_vector_free (nodes);

		nc_dataset_append_data (ds, dca_binned);
  }
	else
		nc_dataset_append_data (ds, dca_unbinned);
}
