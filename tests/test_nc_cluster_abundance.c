/***************************************************************************
 *            test_nc_cluster_abundance.c
 *
 *  Tue April 19 14:34:27 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2022 <sandro@isoftware.com.br>
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcClusterAbundance
{
  NcmMSet *mset;
  NcClusterAbundance *cad;
  NcDataClusterNCount *ncdata;
  gdouble area;
  guint ntests;
  NcClusterMass *clusterm;
  NcClusterRedshift *clusterz;
  NcHICosmo *cosmo;
} TestNcClusterAbundance;

void test_nc_cluster_abundance_new (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_free (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_sanity (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_serialize (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_true_n (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_n (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_intp_d2n (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_intp_bin_d2n (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_intp_d2n_bias (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_intp_bin_d2n_bias (TestNcClusterAbundance *test, gconstpointer pdata);
void test_nc_cluster_abundance_mean_bias (TestNcClusterAbundance *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/nc/cluster_abundance/sanity", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_sanity,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/serialize", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_serialize,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_true_n", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_true_n,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_n", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_n,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_intp_d2n", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_intp_d2n,
              &test_nc_cluster_abundance_free);

  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_intp_bin_d2n", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_intp_bin_d2n,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_intp_d2n_bias", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_intp_d2n_bias,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_intp_bin_d2n_bias", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_intp_bin_d2n_bias,
              &test_nc_cluster_abundance_free);
  
  g_test_add ("/nc/cluster_abundance/nc_cluster_abundance_mean_bias", TestNcClusterAbundance, NULL,
              &test_nc_cluster_abundance_new,
              &test_nc_cluster_abundance_mean_bias,
              &test_nc_cluster_abundance_free);
  
  g_test_run ();
}

void
test_nc_cluster_abundance_new (TestNcClusterAbundance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcHIReion *reion            = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim              = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcDistance *dist            = nc_distance_new (3.0);
  NcTransferFunc *tf          = nc_transfer_func_new_from_name ("NcTransferFuncEH");
  NcPowspecML *ps_ml          = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf       = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mulf    = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL, 200));
  NcHaloMassFunction *mfp     = nc_halo_mass_function_new (dist, psf, mulf);
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (nc_cluster_mass_new_from_name ("NcClusterMassAscaso"));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (nc_cluster_redshift_new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<0.1>, 'pz-max':<1.0>, 'z-bias':<0.0>, 'sigma0':<1.0e-2>}"));
  NcHaloBiasType *hbias_ps    = nc_halo_bias_type_ps_new (200);
  NcHaloBiasFunc *mbiasf      = nc_halo_bias_func_new (mfp, hbias_ps);
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
  
  test->cosmo    = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  test->cad      = nc_cluster_abundance_new (mfp, mbiasf);
  test->mset     = ncm_mset_new (cosmo, clusterm, clusterz, NULL);
  test->ncdata   = nc_data_cluster_ncount_new (test->cad, "NcClusterPhotozGaussGlobal", "NcClusterMassAscaso");
  test->clusterm = clusterm;
  test->clusterz = clusterz;
  test->area     = g_test_rand_double_range (0.21, 0.27) * 4.0 * ncm_c_pi () / 100.0;
  
  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
  
  nc_halo_mass_function_free (mfp);
  nc_multiplicity_func_free (mulf);
  ncm_powspec_filter_free (psf);
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
  nc_halo_bias_func_free (mbiasf);
  nc_distance_free (dist);
}

void
test_nc_cluster_abundance_free (TestNcClusterAbundance *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_data_cluster_ncount_free, test->ncdata);
  NCM_TEST_FREE (nc_cluster_abundance_free, test->cad);
  NCM_TEST_FREE (ncm_mset_free, test->mset);
  NCM_TEST_FREE (ncm_model_free, NCM_MODEL (test->clusterm));
  NCM_TEST_FREE (ncm_model_free, NCM_MODEL (test->clusterz));
}

void
test_nc_cluster_abundance_sanity (TestNcClusterAbundance *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_new (NULL);
  
  nc_data_cluster_ncount_init_from_sampling (test->ncdata, test->mset, test->area, rng);
  
  {
    gchar *desc = ncm_data_get_desc (NCM_DATA (test->ncdata));
    
    g_assert (strlen (desc) > 0);
    g_free (desc);
  }
  
  g_assert (nc_data_cluster_ncount_has_lnM_true (test->ncdata));
  g_assert (nc_data_cluster_ncount_has_z_true (test->ncdata));
  
  g_assert_cmpuint (nc_data_cluster_ncount_get_len (test->ncdata), >, 0);
  g_assert_cmpuint (nc_data_cluster_ncount_lnM_obs_len (test->ncdata), ==, 1);
  g_assert_cmpuint (nc_data_cluster_ncount_lnM_obs_params_len (test->ncdata), ==, 0);
  g_assert_cmpuint (nc_data_cluster_ncount_z_obs_len (test->ncdata), ==, 1);
  g_assert_cmpuint (nc_data_cluster_ncount_z_obs_params_len (test->ncdata), ==, 0);
  
  {
    NcmVector *v;
    
    v = nc_data_cluster_ncount_get_lnM_true (test->ncdata);
    ncm_vector_free (v);
    
    v = nc_data_cluster_ncount_get_z_true (test->ncdata);
    ncm_vector_free (v);
  }
  
  {
    NcmMatrix *m;
    
    m = nc_data_cluster_ncount_get_lnM_obs (test->ncdata);
    ncm_matrix_free (m);
    
    m = nc_data_cluster_ncount_get_lnM_obs_params (test->ncdata);
    g_assert_null (m);
    
    m = nc_data_cluster_ncount_get_z_obs (test->ncdata);
    ncm_matrix_free (m);
    
    m = nc_data_cluster_ncount_get_z_obs_params (test->ncdata);
    g_assert_null (m);
  }
  
  nc_data_cluster_ncount_true_data (test->ncdata, TRUE);
  g_assert (nc_data_cluster_ncount_using_true_data (test->ncdata));
  
  nc_data_cluster_ncount_true_data (test->ncdata, FALSE);
  g_assert (!nc_data_cluster_ncount_using_true_data (test->ncdata));
  
  ncm_rng_free (rng);
}

void
test_nc_cluster_abundance_serialize (TestNcClusterAbundance *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GVariant *var     = ncm_serialize_to_variant (ser, G_OBJECT (test->cad));
  
  nc_cluster_abundance_free (test->cad);
  test->cad = NC_CLUSTER_ABUNDANCE (ncm_serialize_from_variant (ser, var));
  
  g_variant_unref (var);
  ncm_serialize_free (ser);
}

void
test_nc_cluster_abundance_true_n (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  {
    gdouble n1 = nc_cluster_abundance_true_n (test->cad, test->cosmo, test->clusterz, test->clusterm);
    gdouble n2 = nc_cluster_abundance_true_n (test->cad, test->cosmo, test->clusterz, test->clusterm);

    g_assert (gsl_finite (n1));
    g_assert (gsl_finite (n2));

    ncm_assert_cmpdouble_e (n1, ==, n2, 1.0e-15, 0.0);
  }
}

void
test_nc_cluster_abundance_n (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  {
    gdouble N1 = nc_cluster_abundance_n (test->cad, test->cosmo, test->clusterz, test->clusterm);
    gdouble N2 = nc_cluster_abundance_n (test->cad, test->cosmo, test->clusterz, test->clusterm);

    g_assert (gsl_finite (N1));
    g_assert (gsl_finite (N2));

    ncm_assert_cmpdouble_e (N1, ==, N2, 1.0e-15, 0.0);
  }
}

void
test_nc_cluster_abundance_intp_d2n (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  {
    gdouble d2n_1 = nc_cluster_abundance_intp_d2n (test->cad, test->cosmo, test->clusterz, test->clusterm, 14, 1);
    gdouble d2n_2 = nc_cluster_abundance_intp_d2n (test->cad, test->cosmo, test->clusterz, test->clusterm, 14, 1);

    g_assert (gsl_finite (d2n_1));
    g_assert (gsl_finite (d2n_2));

    ncm_assert_cmpdouble_e (d2n_1, ==, d2n_2, 1.0e-15, 0.0);
  }
}

void
test_nc_cluster_abundance_intp_bin_d2n (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  {
    const gdouble lnM_obs_bin_lower[1] = {0.0};
    const gdouble lnM_obs_bin_upper[1] = {0.5};
    const gdouble z_obs_bin_lower[1]   = {0.1};
    const gdouble z_obs_bin_upper[1]   = {0.3};
    gdouble d2n_bin1                   = nc_cluster_abundance_intp_bin_d2n (test->cad, test->cosmo, test->clusterz, test->clusterm, lnM_obs_bin_lower, lnM_obs_bin_upper, NULL, z_obs_bin_lower, z_obs_bin_upper, NULL);
    gdouble d2n_bin2                   = nc_cluster_abundance_intp_bin_d2n (test->cad, test->cosmo, test->clusterz, test->clusterm, lnM_obs_bin_lower, lnM_obs_bin_upper, NULL, z_obs_bin_lower, z_obs_bin_upper, NULL);

    g_assert (gsl_finite (d2n_bin1));
    g_assert (gsl_finite (d2n_bin2));

    ncm_assert_cmpdouble_e (d2n_bin1, ==, d2n_bin2, 1.0e-15, 0.0);
  }
}

void
test_nc_cluster_abundance_intp_d2n_bias (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  gdouble lnM_obs[1] = {1.0};
  gdouble z_obs[1]   = {0.3};
  gdouble d2n_bias1  = nc_cluster_abundance_intp_d2n_bias (test->cad, test->cosmo, test->clusterz, test->clusterm, lnM_obs, NULL, z_obs, NULL);
  gdouble d2n_bias2  = nc_cluster_abundance_intp_d2n_bias (test->cad, test->cosmo, test->clusterz, test->clusterm, lnM_obs, NULL, z_obs, NULL);

  g_assert (gsl_finite (d2n_bias1));
  g_assert (gsl_finite (d2n_bias2));
  
  ncm_assert_cmpdouble_e (d2n_bias1, ==, d2n_bias2, 1.0e-15, 0.0);
}

void
test_nc_cluster_abundance_intp_bin_d2n_bias (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  {
    const gdouble lnM_obs_bin_lower[1] = {0.4};
    const gdouble lnM_obs_bin_upper[1] = {1.0};
    const gdouble z_obs_bin_lower[1]   = {0.1};
    const gdouble z_obs_bin_upper[1]   = {0.3};
    gdouble d2n_bias_bin1 = nc_cluster_abundance_intp_bin_d2n_bias (test->cad, test->cosmo, test->clusterz, test->clusterm, lnM_obs_bin_lower, lnM_obs_bin_upper, NULL, z_obs_bin_lower, z_obs_bin_upper, NULL);
    gdouble d2n_bias_bin2 = nc_cluster_abundance_intp_bin_d2n_bias (test->cad, test->cosmo, test->clusterz, test->clusterm, lnM_obs_bin_lower, lnM_obs_bin_upper, NULL, z_obs_bin_lower, z_obs_bin_upper, NULL);

    g_assert (gsl_finite (d2n_bias_bin1));
    g_assert (gsl_finite (d2n_bias_bin2));

    ncm_assert_cmpdouble_e (d2n_bias_bin1, ==, d2n_bias_bin2, 1.0e-15, 0.0);
  }
}

void
test_nc_cluster_abundance_mean_bias (TestNcClusterAbundance *test, gconstpointer pdata)
{
  test_nc_cluster_abundance_sanity (test, pdata);
  {
    gdouble mean_bias1 = nc_cluster_abundance_mean_bias (test->cad, test->cosmo, test->clusterz, test->clusterm);
    gdouble mean_bias2 = nc_cluster_abundance_mean_bias (test->cad, test->cosmo, test->clusterz, test->clusterm);

    g_assert (gsl_finite (mean_bias1));
    g_assert (gsl_finite (mean_bias2));

    ncm_assert_cmpdouble_e (mean_bias1, ==, mean_bias2, 1.0e-15, 0.0);
  }
}
