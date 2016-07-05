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
  gint filter_type;

  if (de_data_cluster->filter_type != NULL)
  {
    const GEnumValue *filter_type_id = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_POWSPEC_FILTER_TYPE, de_data_cluster->filter_type);
    if (filter_type_id == NULL)
    {
      g_message ("DataCluster: Filter type `%s' not found. Use one from the following list:", de_data_cluster->filter_type);
      ncm_cfg_enum_print_all (NCM_TYPE_POWSPEC_FILTER_TYPE, "Powerspectrum filters");
      g_error ("DataCluster: Giving up");
    }
    filter_type = filter_type_id->value;
  }
  else
    filter_type = NCM_POWSPEC_FILTER_TYPE_TOPHAT;

  if (de_data_cluster->ps_type == NULL)
  {
    de_data_cluster->ps_type = g_strdup ("NcPowspecMLTransfer{'transfer' : <{'NcTransferFuncEH', @a{sv} {}}>}");
  }

  if (de_data_cluster->multiplicity_name == NULL)
  {
    de_data_cluster->multiplicity_name = g_strdup ("NcMultiplicityFuncTinkerMean");
  }

  if (de_data_cluster->clusterm_ser == NULL)
  {
    de_data_cluster->clusterm_ser = g_strdup ("NcClusterMassNodist");
  }

  if (de_data_cluster->clusterz_ser == NULL)
  {
    de_data_cluster->clusterz_ser = g_strdup ("NcClusterRedshiftNodist");
  }

  {
    NcClusterMass *clusterm     = nc_cluster_mass_new_from_name (de_data_cluster->clusterm_ser);
    NcClusterRedshift *clusterz = nc_cluster_redshift_new_from_name (de_data_cluster->clusterz_ser);
    
    ncm_mset_set (mset, NCM_MODEL (clusterm));
    ncm_mset_set (mset, NCM_MODEL (clusterz));
  }

  {
    NcmPowspec *ps              = NCM_POWSPEC (ncm_serialize_global_from_string (de_data_cluster->ps_type));
    NcmPowspecFilter *psf       = ncm_powspec_filter_new (ps, filter_type);
    NcMultiplicityFunc *mulf    = nc_multiplicity_func_new_from_name (de_data_cluster->multiplicity_name);
    NcHaloMassFunction *mfp     = nc_halo_mass_function_new (dist, psf, mulf);
    NcClusterAbundance *cad     = nc_cluster_abundance_nodist_new (mfp, NULL);
    NcDataClusterNCount *ncount = nc_data_cluster_ncount_new (cad);
    
    ncm_powspec_clear (&ps);
    ncm_powspec_filter_clear (&psf);
    nc_multiplicity_func_free (mulf);
    nc_cluster_abundance_free (cad);

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

          nc_data_cluster_ncount_catalog_load (ncount, de_data_cluster->cata_file[i]);
          nc_data_cluster_ncount_true_data (ncount, de_data_cluster->use_true_data);

          _nc_de_data_cluster_append (de_data_cluster, NCM_DATA (ncount), dset);
          g_ptr_array_add (ca_array, NCM_DATA (ncount));
          if ((i == 0) && (de_data_cluster->save_cata != NULL))
            nc_data_cluster_ncount_catalog_save (ncount, de_data_cluster->save_cata, TRUE);
          i++;
        }
        break;
      }
#endif /* HAVE_CONFIG_H */
      case NC_DATA_CLUSTER_ABUNDANCE_SAMPLING:
      {
        nc_data_cluster_ncount_init_from_sampling (ncount, mset, de_data_cluster->area_survey * gsl_pow_2 (M_PI / 180.0), rng);
        nc_data_cluster_ncount_true_data (ncount, de_data_cluster->use_true_data);

        if (de_data_cluster->save_cata != NULL)
#ifdef NUMCOSMO_HAVE_CFITSIO
          nc_data_cluster_ncount_catalog_save (ncount, de_data_cluster->save_cata, TRUE);
#else
          g_error ("darkenergy: cannot save file numcosmo built without support for fits files");
#endif /* HAVE_CONFIG_H */
        _nc_de_data_cluster_append (de_data_cluster, NCM_DATA (ncount), dset);
        g_ptr_array_add (ca_array, NCM_DATA (ncount));
      }
        break;
      default:
        g_error ("The option --catalog-id must be between (0,2).");
    }

    nc_halo_mass_function_free (mfp);

    return ca_array;
  }
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
