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
  NcmMatrix *z_obs_params;
  NcmVector *lnM_obs;
  NcmMatrix *lnM_obs_params;
  NcmMatrix *s_matrix;
  NcmMatrix *resample_s_matrix;
  gdouble area;
  guint ntests;
} TestNcClusterNCountsGauss;

void test_nc_data_cluster_ncounts_gauss_new (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_free (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_sanity (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_serialize (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_s_matrix (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_mean_func (TestNcClusterNCountsGauss *test, gconstpointer pdata);
void test_nc_data_cluster_ncounts_gauss_cov (TestNcClusterNCountsGauss *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/data_cluster_ncounts_gauss/sanity", TestNcClusterNCountsGauss, NULL,
              &test_nc_data_cluster_ncounts_gauss_new,
              &test_nc_data_cluster_ncounts_gauss_sanity,
              &test_nc_data_cluster_ncounts_gauss_free);

  g_test_add ("/nc/data_cluster_ncounts_gauss/serialize", TestNcClusterNCountsGauss, NULL,
              &test_nc_data_cluster_ncounts_gauss_new,
              &test_nc_data_cluster_ncounts_gauss_serialize,
              &test_nc_data_cluster_ncounts_gauss_free);

  g_test_add ("/nc/data_cluster_ncounts_gauss/s_matrix", TestNcClusterNCountsGauss, NULL,
              &test_nc_data_cluster_ncounts_gauss_new,
              &test_nc_data_cluster_ncounts_gauss_s_matrix,
              &test_nc_data_cluster_ncounts_gauss_free);

  g_test_add ("/nc/data_cluster_ncounts_gauss/mean_func", TestNcClusterNCountsGauss, NULL,
              &test_nc_data_cluster_ncounts_gauss_new,
              &test_nc_data_cluster_ncounts_gauss_mean_func,
              &test_nc_data_cluster_ncounts_gauss_free);

  g_test_add ("/nc/data_cluster_ncounts_gauss/cov", TestNcClusterNCountsGauss, NULL,
              &test_nc_data_cluster_ncounts_gauss_new,
              &test_nc_data_cluster_ncounts_gauss_cov,
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
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassNodist"));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_serialize_global_from_string ("NcClusterRedshiftNodist{'z-min':<0.1>, 'z-max':<1.0>}"));
  gdouble z_obs_array[8]      = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  gdouble lnM_obs_array[2]    = {34, 35};
  NcmMatrix *s                = ncm_matrix_new (7, 7);
  NcmMatrix *s_resample       = ncm_matrix_new (7, 7);
  guint i, j;

  ncm_matrix_set_identity (s_resample);

  for (i = 0; i < ncm_matrix_row_len (s); i++)
  {
    for (j = 0; j < ncm_matrix_col_len (s); j++)
    {
      ncm_matrix_set (s, i, j, 2);
    }
  }

  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));


  test->z_obs             = ncm_vector_new_data_dup (z_obs_array, 8, 1);
  test->z_obs_params      = ncm_matrix_new0 (7, 7);
  test->lnM_obs           = ncm_vector_new_data_dup (lnM_obs_array, 2, 1);
  test->lnM_obs_params    = ncm_matrix_new0 (7, 7);
  test->cad               = nc_cluster_abundance_new (mfp, hbias);
  test->mset              = ncm_mset_new (cosmo, clusterm, clusterz, NULL);
  test->ncounts_gauss     = nc_data_cluster_ncounts_gauss_new (test->cad);
  test->area              = g_test_rand_double_range (0.21, 0.27) * 4.0 * ncm_c_pi () / 100.0;
  test->s_matrix          = s;
  test->resample_s_matrix = s_resample;

  ncm_data_set_init (NCM_DATA (test->ncounts_gauss), TRUE);

  nc_data_cluster_ncounts_gauss_set_z_obs (test->ncounts_gauss, test->z_obs);
  nc_data_cluster_ncounts_gauss_set_lnM_obs (test->ncounts_gauss, test->lnM_obs);

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
  const guint np = (ncm_vector_len (test->z_obs) - 1) * (ncm_vector_len (test->lnM_obs) - 1);
  {
    gchar *desc = ncm_data_get_desc (NCM_DATA (test->ncounts_gauss));

    g_assert_true (strlen (desc) > 0);
    g_free (desc);
  }

  {
    NcmMatrix *m;
    NcmVector *v;

    v = nc_data_cluster_ncounts_gauss_get_lnM_obs (test->ncounts_gauss);
    g_assert_cmpuint (ncm_vector_len (v), ==, ncm_vector_len (test->lnM_obs));
    g_assert_true (v == test->lnM_obs);
    ncm_vector_free (v);

    m = nc_data_cluster_ncounts_gauss_get_lnM_obs_params (test->ncounts_gauss);
    g_assert_null (m);
    nc_data_cluster_ncounts_gauss_set_lnM_obs_params (test->ncounts_gauss, test->lnM_obs_params);
    m = nc_data_cluster_ncounts_gauss_get_lnM_obs_params (test->ncounts_gauss);
    g_assert_true (m == test->lnM_obs_params);
    ncm_matrix_free (m);

    v = nc_data_cluster_ncounts_gauss_get_z_obs (test->ncounts_gauss);
    g_assert_cmpuint (ncm_vector_len (v), ==, ncm_vector_len (test->z_obs));
    g_assert_true (v == test->z_obs);
    ncm_vector_free (v);

    m = nc_data_cluster_ncounts_gauss_get_z_obs_params (test->ncounts_gauss);
    g_assert_null (m);
    nc_data_cluster_ncounts_gauss_set_z_obs_params (test->ncounts_gauss, test->z_obs_params);
    m = nc_data_cluster_ncounts_gauss_get_z_obs_params (test->ncounts_gauss);
    g_assert_true (m == test->z_obs_params);
    ncm_matrix_free (m);
  }

  {
    g_assert_true (nc_data_cluster_ncounts_gauss_get_has_ssc (test->ncounts_gauss) == FALSE);
    g_assert_true (nc_data_cluster_ncounts_gauss_get_fix_cov (test->ncounts_gauss) == FALSE);
  }

  {
    NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (test->ncounts_gauss);

    ncm_data_set_init (NCM_DATA (gauss_cov), TRUE);
    g_assert_true (ncm_data_gauss_cov_get_size (gauss_cov) == np);

    nc_data_cluster_ncounts_gauss_set_z_obs (test->ncounts_gauss, test->z_obs);
  }

  {
    NcmDataGaussCov *gauss_cov_2 = NCM_DATA_GAUSS_COV (test->ncounts_gauss);

    ncm_data_set_init (NCM_DATA (gauss_cov_2), TRUE);
    g_assert_true (ncm_data_gauss_cov_get_size (gauss_cov_2) == np);
  }
}

void
test_nc_data_cluster_ncounts_gauss_serialize (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{
  {
    nc_data_cluster_ncounts_gauss_set_z_obs_params (test->ncounts_gauss, test->z_obs_params);
    nc_data_cluster_ncounts_gauss_set_lnM_obs_params (test->ncounts_gauss, test->lnM_obs_params);
    nc_data_cluster_ncounts_gauss_set_s_matrix (test->ncounts_gauss, test->s_matrix);
    nc_data_cluster_ncounts_gauss_set_resample_s_matrix (test->ncounts_gauss, test->resample_s_matrix);
  }
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    GVariant *var     = ncm_serialize_to_variant (ser, G_OBJECT (test->ncounts_gauss));

    nc_data_cluster_ncounts_gauss_free (test->ncounts_gauss);
    test->ncounts_gauss = NC_DATA_CLUSTER_NCOUNTS_GAUSS (ncm_serialize_from_variant (ser, var));

    g_variant_unref (var);
    ncm_serialize_free (ser);
  }
}

void
test_nc_data_cluster_ncounts_gauss_s_matrix (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{
  nc_data_cluster_ncounts_gauss_set_has_ssc (test->ncounts_gauss, TRUE);
  g_assert_true (nc_data_cluster_ncounts_gauss_get_has_ssc (test->ncounts_gauss) == TRUE);

  g_assert_null (nc_data_cluster_ncounts_gauss_get_s_matrix (test->ncounts_gauss));
  g_assert_null (nc_data_cluster_ncounts_gauss_get_resample_s_matrix (test->ncounts_gauss));


  nc_data_cluster_ncounts_gauss_set_s_matrix (test->ncounts_gauss, test->s_matrix);
  nc_data_cluster_ncounts_gauss_set_resample_s_matrix (test->ncounts_gauss, test->resample_s_matrix);

  g_assert_true (nc_data_cluster_ncounts_gauss_get_s_matrix (test->ncounts_gauss) == test->s_matrix);
  g_assert_true (nc_data_cluster_ncounts_gauss_get_resample_s_matrix (test->ncounts_gauss) == test->resample_s_matrix);
}

void
test_nc_data_cluster_ncounts_gauss_mean_func (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (test->ncounts_gauss);
  NcmVector *mean;
  guint i;

  ncm_data_set_init (NCM_DATA (gauss_cov), TRUE);

  mean = ncm_vector_dup (ncm_data_gauss_cov_compute_mean (gauss_cov, test->mset));

  for (i = 0; i < ncm_vector_len (mean); i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (mean, i)));
  }
}

void
test_nc_data_cluster_ncounts_gauss_cov (TestNcClusterNCountsGauss *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (test->ncounts_gauss);
  NcmRNG *rng                = ncm_rng_new (NULL);
  NcmMatrix *cov_resample;
  NcmMatrix *cov_diag;
  NcmMatrix *cov_1;
  NcmMatrix *cov_2;
  NcmMatrix *cov_3;
  NcmMatrix *cov_4;
  guint i, j;

  nc_data_cluster_ncounts_gauss_set_s_matrix (test->ncounts_gauss, test->s_matrix);
  nc_data_cluster_ncounts_gauss_set_resample_s_matrix (test->ncounts_gauss, test->resample_s_matrix);

  ncm_data_set_init (NCM_DATA (gauss_cov), TRUE);

  g_assert_false (nc_data_cluster_ncounts_gauss_get_fix_cov (test->ncounts_gauss));

  cov_diag = ncm_matrix_dup (ncm_data_gauss_cov_compute_cov (gauss_cov, test->mset, NULL));
  ncm_data_gauss_cov_set_cov (gauss_cov, cov_diag);

  for (i = 0; i < ncm_matrix_row_len (cov_diag); i++)
  {
    for (j = 0; j < ncm_matrix_col_len (cov_diag); j++)
    {
      g_assert_true (gsl_finite (ncm_matrix_get (cov_diag, i, j)));

      if (i != j)
        ncm_assert_cmpdouble_e (ncm_matrix_get (cov_diag, i, j), ==, 0, 1.0e-15, 0.0);

      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_diag, i, j), ==, ncm_matrix_get (ncm_data_gauss_cov_peek_cov (gauss_cov), i, j), 1.0e-15, 0.0);
    }
  }

  nc_data_cluster_ncounts_gauss_set_has_ssc (test->ncounts_gauss, TRUE);

  cov_1 = ncm_matrix_dup (ncm_data_gauss_cov_compute_cov (gauss_cov, test->mset,  NULL));

  ncm_data_resample (NCM_DATA (test->ncounts_gauss), test->mset, rng);

  cov_resample = ncm_matrix_dup (ncm_data_gauss_cov_peek_cov (gauss_cov));

  for (i = 0; i < ncm_matrix_row_len (cov_1); i++)
  {
    for (j = 0; j < ncm_matrix_col_len (cov_1); j++)
    {
      g_assert_true (gsl_finite (ncm_matrix_get (cov_diag, i, j)));
      g_assert_true (gsl_finite (ncm_matrix_get (cov_1, i, j)));
      g_assert_true (gsl_finite (ncm_matrix_get (cov_resample, i, j)));
      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_1, i, j),  !=, ncm_matrix_get (cov_resample, i, j), 1.0e-15, 0.0);

      if (i != j)
        ncm_assert_cmpdouble_e (ncm_matrix_get (cov_diag, i, j), ==, 0, 1.0e-15, 0.0);
    }
  }


  nc_data_cluster_ncounts_gauss_set_fix_cov (test->ncounts_gauss, TRUE);

  g_assert_true (nc_data_cluster_ncounts_gauss_get_fix_cov (test->ncounts_gauss) == TRUE);


  cov_2 = ncm_matrix_dup (ncm_data_gauss_cov_compute_cov (gauss_cov, test->mset, NULL));

  for (i = 0; i < ncm_matrix_row_len (cov_2); i++)
  {
    for (j = 0; j < ncm_matrix_col_len (cov_2); j++)
    {
      g_assert_true (gsl_finite (ncm_matrix_get (cov_2, i, j)));
      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_1, i, j),  !=, ncm_matrix_get (cov_2, i, j), 1.0e-15, 0.0);
      g_assert_true (ncm_matrix_get (cov_2, i, j)  ==  ncm_matrix_get (cov_resample, i, j));
    }
  }

  nc_data_cluster_ncounts_gauss_set_fix_cov (test->ncounts_gauss, FALSE);
  cov_3 = ncm_matrix_dup (ncm_data_gauss_cov_compute_cov (gauss_cov, test->mset, NULL));

  for (i = 0; i < ncm_matrix_row_len (cov_1); i++)
  {
    for (j = 0; j < ncm_matrix_col_len (cov_1); j++)
    {
      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_1, i, j),  ==, ncm_matrix_get (cov_3, i, j), 1.0e-15, 0.0);
    }
  }

  nc_data_cluster_ncounts_gauss_set_has_ssc (test->ncounts_gauss, FALSE);

  nc_data_cluster_ncounts_gauss_set_fix_cov (test->ncounts_gauss, FALSE);
  cov_4 = ncm_matrix_dup (ncm_data_gauss_cov_compute_cov (gauss_cov, test->mset, NULL));

  for (i = 0; i < ncm_matrix_row_len (cov_1); i++)
  {
    for (j = 0; j < ncm_matrix_col_len (cov_1); j++)
    {
      ncm_assert_cmpdouble_e (ncm_matrix_get (cov_4, i, j),  !=, ncm_matrix_get (cov_3, i, j), 1.0e-15, 0.0);
    }
  }
}

