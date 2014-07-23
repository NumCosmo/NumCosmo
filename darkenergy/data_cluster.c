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

static void _nc_de_data_cluster_append (NcDEDataClusterEntries *de_data_cluster, NcmData *dca_unbinned, NcmDataset *dset);

GPtrArray *
nc_de_data_cluster_new (NcDistance *dist, NcmMSet *mset, NcDEDataClusterEntries *de_data_cluster, NcmDataset *dset, NcDataClusterAbundanceId id, NcmRNG *rng)
{
  GPtrArray *ca_array = g_ptr_array_new ();
  NcWindow *wp = nc_window_new_from_name (de_data_cluster->window_name);
  NcTransferFunc *tf = nc_transfer_func_new_from_name (de_data_cluster->transfer_name);
  NcMatterVar *vp = nc_matter_var_new (NC_MATTER_VAR_FFT, wp, tf);
  NcGrowthFunc *gf = nc_growth_func_new ();
  NcMultiplicityFunc *mulf = nc_multiplicity_func_new_from_name (de_data_cluster->multiplicity_name);
  NcMassFunction *mfp = nc_mass_function_new (dist, vp, gf, mulf);

  nc_window_free (wp);
  nc_transfer_func_free (tf);
  nc_matter_var_free (vp);
  nc_growth_func_free (gf);
  nc_multiplicity_func_free (mulf);

  switch (id)
  {
#ifdef NUMCOSMO_HAVE_CFITSIO
    case NC_DATA_CLUSTER_ABUNDANCE_FIT:
    {
    gint i = 0;
    if (de_data_cluster->cata_file == NULL)
      g_error ("For --cluster-id 0, you must specify a fit catalog via --catalog file.fit");
    while (de_data_cluster->cata_file[i] != NULL)
    {
      NcClusterAbundance *cad = nc_cluster_abundance_nodist_new (mfp, NULL);
      NcDataClusterNCount *ncount = nc_data_cluster_ncount_new (cad);
      nc_cluster_abundance_free (cad);

      nc_data_cluster_ncount_catalog_load (ncount, de_data_cluster->cata_file[i]);

      ncm_mset_set (mset, NCM_MODEL (ncount->m));

      nc_data_cluster_ncount_true_data (ncount, de_data_cluster->use_true_data);
      _nc_de_data_cluster_append (de_data_cluster, NCM_DATA (ncount), dset);
      g_ptr_array_add (ca_array, NCM_DATA (ncount));
      if ((i == 0) && (de_data_cluster->save_cata != NULL))
        nc_data_cluster_ncount_catalog_save (ncount, de_data_cluster->save_cata, TRUE);
      i++;
    }
  }
      break;
    case NC_DATA_CLUSTER_ABUNDANCE_TXT:
    {
      gint i = 0;
      if (de_data_cluster->cata_file == NULL)
        g_error ("For --cluster-id 1, you must specify a fit catalog via --catalog filename");
      while (de_data_cluster->cata_file[i] != NULL)
      {
        NcClusterAbundance *cad = nc_cluster_abundance_nodist_new (mfp, NULL);
        NcDataClusterNCount *ncount = nc_data_cluster_ncount_new (cad);
        nc_cluster_abundance_free (cad);

        //nc_data_cluster_abundance_unbinned_init_from_text_file (dca_unbinned, de_data_cluster->cata_file[i], opt, de_data_cluster->area_survey * gsl_pow_2 (M_PI / 180.0), log(de_data_cluster->Mi), log(de_data_cluster->Mf), de_data_cluster->z_initial, de_data_cluster->z_final, de_data_cluster->photoz_sigma0, de_data_cluster->photoz_bias, de_data_cluster->lnM_sigma0, de_data_cluster->lnM_bias);
        g_assert_not_reached ();
        nc_data_cluster_ncount_true_data (ncount, de_data_cluster->use_true_data);

        _nc_de_data_cluster_append (de_data_cluster, NCM_DATA (ncount), dset);
        g_ptr_array_add (ca_array, NCM_DATA (ncount));
        if ((i == 0) && (de_data_cluster->save_cata != NULL))
          nc_data_cluster_ncount_catalog_save (ncount, de_data_cluster->save_cata, TRUE);
        i++;
      }
    }
      break;
#endif /* HAVE_CONFIG_H */
    case NC_DATA_CLUSTER_ABUNDANCE_SAMPLING:
    {
      NcClusterMass *clusterm = nc_cluster_mass_new_from_name (de_data_cluster->clusterm_ser);
      NcClusterRedshift *clusterz = nc_cluster_redshift_new_from_name (de_data_cluster->clusterz_ser);
      NcClusterAbundance *cad = nc_cluster_abundance_new (mfp, NULL, clusterz, clusterm);
      NcDataClusterNCount *ncount = nc_data_cluster_ncount_new (cad);

      ncm_mset_set (mset, NCM_MODEL (clusterm));
      nc_cluster_abundance_free (cad);

      nc_data_cluster_ncount_init_from_sampling (ncount, mset, clusterz, clusterm, de_data_cluster->area_survey * gsl_pow_2 (M_PI / 180.0), rng);
      nc_data_cluster_ncount_true_data (ncount, de_data_cluster->use_true_data);

      if (de_data_cluster->save_cata != NULL)
#ifdef NUMCOSMO_HAVE_CFITSIO
        nc_data_cluster_ncount_catalog_save (ncount, de_data_cluster->save_cata, TRUE);
#else
      g_error ("darkenergy: cannot save file numcosmo built without support for fits files");
#endif /* HAVE_CONFIG_H */
      _nc_de_data_cluster_append (de_data_cluster, NCM_DATA (ncount), dset);
      g_ptr_array_add (ca_array, NCM_DATA (ncount));
      nc_cluster_mass_free (clusterm);
      nc_cluster_redshift_free (clusterz);
    }
      break;
    default:
      g_error ("The option --catalog-id must be between (0,2).");
  }

  nc_mass_function_free (mfp);

  return ca_array;
}

static void
_nc_de_data_cluster_append (NcDEDataClusterEntries *de_data_cluster, NcmData *data, NcmDataset *dset)
{
  if (de_data_cluster->binned)
  {
    g_assert_not_reached ();
  }
  else
  {
    ncm_dataset_append_data (dset, data);
    ncm_data_free (data);
  }
}
