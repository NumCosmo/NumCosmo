/***************************************************************************
 *            test_ncm_data_cluster_ncounts_gauss.c
 *
 *  Tue April 03 16:02:26 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <gsl/gsl_statistics_double.h>

#include "ncm_data_gauss_cov_test.h"

typedef struct _TestNcClusterNCountsGauss
{
  NcmMSet *mset;
  NcDataClusterNCountsGauss *ncounts_gauss;
  NcClusterAbundance *cad;
  NcmVector *z_obs;
  NcmVector *lnM_obs;
  gdouble area;
  guint ntests;

} TestNcClusterNCountsGauss;

void test_nc_data_cluster_ncounts_gauss_new (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_free (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_sanity (TestNcClusterNCountsGauss *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_cluster_ncounts_gauss/sanity", TestNcClusterNCountsGauss, NULL,
              &test_nc_data_cluster_ncounts_gauss_new,
              &test_nc_data_cluster_ncounts_gauss_sanity,
              &test_nc_data_cluster_ncounts_gauss_free);


  g_test_run ();

}

void
test_nc_data_cluster_ncounts_gauss_new (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHIReion *reion            = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim              = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcDistance *dist            = nc_distance_new (3.0);
  NcTransferFunc *tf          = NC_TRANSFER_FUNC (ncm_serialize_global_from_string ("NcTransferFuncEH"));
  NcPowspecML *ps_ml          = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf       = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mulf    = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL, 500.0));
  NcHaloMassFunction *mfp     = nc_halo_mass_function_new (dist, psf, mulf);
  NcHaloBias *hbias           = NC_HALO_BIAS (nc_halo_bias_ps_new (mfp));
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassAscaso"));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_serialize_global_from_string ("NcClusterRedshiftNodist{'z-min':<0.1>, 'z-max':<1.0>}"));
  gdouble z_obs_array[8] = {0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8};
  gdouble lnM_obs_array[2] = {0.1 , 2.0};
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));


  test->z_obs  = ncm_vector_new_data_static (z_obs_array, 8, 1);
  test->lnM_obs = ncm_vector_new_data_static (lnM_obs_array, 2, 1);
  test->cad    = nc_cluster_abundance_new (mfp, hbias);
  test->mset   = ncm_mset_new (cosmo, clusterm, clusterz, NULL);
  test->ncounts_gauss = nc_data_cluster_ncounts_gauss_new (test->cad);
  test->area   = g_test_rand_double_range (0.21, 0.27) * 4.0 * ncm_c_pi () / 100.0;

  nc_data_cluster_ncounts_gauss_set_z_obs(test->ncounts_gauss , test->z_obs);
  nc_data_cluster_ncounts_gauss_set_lnM_obs(test->ncounts_gauss , test->lnM_obs);

  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
  ncm_model_free (NCM_MODEL (clusterm));
  ncm_model_free (NCM_MODEL (clusterz));

  nc_halo_mass_function_free (mfp);
  nc_halo_bias_free (hbias);
  nc_multiplicity_func_free (mulf);
  ncm_powspec_filter_free (psf);
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
  nc_distance_free (dist);
}

void
test_nc_data_cluster_ncounts_gauss_free (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{

  NCM_TEST_FREE (ncm_mset_free, test->mset);
}

void
test_nc_data_cluster_ncounts_gauss_sanity (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{ 
  NcmRNG *rng           = ncm_rng_seeded_new (NULL, g_test_rand_int ());


  {
    gchar *desc = ncm_data_get_desc (NCM_DATA (test->ncounts_gauss));

    g_assert_true (strlen (desc) > 0);
    g_free (desc);
  }

  {
    NcmMatrix *m;
    NcmVector *v;
  
    v = nc_data_cluster_ncounts_gauss_get_lnM_obs (test->ncounts_gauss);
    g_assert_cmpuint (ncm_vector_len(v), ==, ncm_vector_len(test->lnM_obs));
    g_assert_true (v == test->lnM_obs);
    ncm_vector_free (v);

    m = nc_data_cluster_ncounts_gauss_get_lnM_obs_params (test->ncounts_gauss);
    g_assert_null (m);

    v = nc_data_cluster_ncounts_gauss_get_z_obs (test->ncounts_gauss);
    g_assert_cmpuint (ncm_vector_len(v), ==, ncm_vector_len(test->z_obs));
    g_assert_true (v == test->z_obs);
    ncm_vector_free (v);

    m = nc_data_cluster_ncounts_gauss_get_z_obs_params (test->ncounts_gauss);
    g_assert_null (m);
  }

  {
    g_assert_true (nc_data_cluster_ncounts_gauss_get_has_ssc(test->ncounts_gauss) == FALSE);
    g_assert_true (nc_data_cluster_ncounts_gauss_get_fix_cov(test->ncounts_gauss) == FALSE);

    g_assert_null (nc_data_cluster_ncounts_gauss_get_s_matrix(test->ncounts_gauss));
    g_assert_null (nc_data_cluster_ncounts_gauss_get_resample_s_matrix(test->ncounts_gauss));

  }

}

