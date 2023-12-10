/***************************************************************************
 *            test_nc_data_cluster_ncount.c
 *
 *  Tue April 19 14:34:27 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2022 <vitenti@uel.br>
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

typedef struct _TestNcDataClusterNCount
{
  NcmMSet *mset;
  NcDataClusterNCount *ncdata;
  NcClusterAbundance *cad;
  gdouble area;
  guint ntests;
} TestNcDataClusterNCount;

void test_nc_data_cluster_ncount_new (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_free (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_sanity (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_bin (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_serialize (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_saveload (TestNcDataClusterNCount *test, gconstpointer pdata);

void test_nc_data_cluster_ncount_m2lnL_true_data_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_binned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_serialize_true_data_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_serialize_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_serialize_binned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_saveload_true_data_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_saveload_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_m2lnL_saveload_binned (TestNcDataClusterNCount *test, gconstpointer pdata);

void test_nc_data_cluster_ncount_traps (TestNcDataClusterNCount *test, gconstpointer pdata);
void test_nc_data_cluster_ncount_invalid_test (TestNcDataClusterNCount *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_cluster_ncount/sanity", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_sanity,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/bin", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_bin,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/serialize", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_serialize,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/saveload", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_saveload,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/true_data/unbinned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_true_data_unbinned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/unbinned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_unbinned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/binned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_binned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/serialize/true_data/unbinned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_serialize_true_data_unbinned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/serialize/unbinned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_serialize_unbinned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/serialize/binned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_serialize_binned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/saveload/true_data/unbinned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_saveload_true_data_unbinned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/saveload/unbinned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_saveload_unbinned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/m2lnL/saveload/binned", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_m2lnL_saveload_binned,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/traps", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_traps,
              &test_nc_data_cluster_ncount_free);

  g_test_add ("/nc/data_cluster_ncount/invalid/test/subprocess", TestNcDataClusterNCount, NULL,
              &test_nc_data_cluster_ncount_new,
              &test_nc_data_cluster_ncount_invalid_test,
              &test_nc_data_cluster_ncount_free);

  g_test_run ();
}

void
test_nc_data_cluster_ncount_new (TestNcDataClusterNCount *test, gconstpointer pdata)
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
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassAscaso"));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_serialize_global_from_string ("NcClusterRedshiftNodist{'z-min':<0.1>, 'z-max':<1.0>}"));

  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));

  test->cad    = nc_cluster_abundance_new (mfp, NULL);
  test->mset   = ncm_mset_new (cosmo, clusterm, clusterz, NULL);
  test->ncdata = nc_data_cluster_ncount_new (test->cad, "NcClusterRedshiftNodist", "NcClusterMassAscaso");
  test->area   = g_test_rand_double_range (0.21, 0.27) * 4.0 * ncm_c_pi () / 100.0;

  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
  ncm_model_free (NCM_MODEL (clusterm));
  ncm_model_free (NCM_MODEL (clusterz));

  nc_halo_mass_function_free (mfp);
  nc_multiplicity_func_free (mulf);
  ncm_powspec_filter_free (psf);
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
  nc_distance_free (dist);
}

void
test_nc_data_cluster_ncount_free (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_data_cluster_ncount_free, test->ncdata);
  NCM_TEST_FREE (nc_cluster_abundance_free, test->cad);
  NCM_TEST_FREE (ncm_mset_free, test->mset);
}

void
test_nc_data_cluster_ncount_sanity (TestNcDataClusterNCount *test, gconstpointer pdata)
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
test_nc_data_cluster_ncount_bin (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  const gdouble zl   = 0.1;
  const gdouble zu   = 1.0;
  const gdouble lnRl = 1.0;
  const gdouble lnRu = 2.0;
  gint i, j;

  nc_data_cluster_ncount_del_bins (test->ncdata);

  for (i = 0; i < 5; i++)
  {
    gdouble bin_zl = zl + (zu - zl) / 5.0 * i;
    gdouble bin_zu = zl + (zu - zl) / 5.0 * (i + 1.0);

    NcmVector *bin_zl_vec = ncm_vector_new_data_static (&bin_zl, 1, 1);
    NcmVector *bin_zu_vec = ncm_vector_new_data_static (&bin_zu, 1, 1);

    for (j = 0; j < 4; j++)
    {
      gdouble bin_lnRl = lnRl + (lnRu - lnRl) / 5.0 * j;
      gdouble bin_lnRu = lnRl + (lnRu - lnRl) / 5.0 * (j + 1.0);

      NcmVector *bin_lnRl_vec = ncm_vector_new_data_static (&bin_lnRl, 1, 1);
      NcmVector *bin_lnRu_vec = ncm_vector_new_data_static (&bin_lnRu, 1, 1);

      nc_data_cluster_ncount_add_bin (test->ncdata, bin_lnRl_vec, bin_lnRu_vec, bin_zl_vec, bin_zu_vec);

      ncm_vector_free (bin_lnRl_vec);
      ncm_vector_free (bin_lnRu_vec);
    }

    ncm_vector_free (bin_zl_vec);
    ncm_vector_free (bin_zu_vec);
  }

  nc_data_cluster_ncount_bin_data (test->ncdata);
}

void
test_nc_data_cluster_ncount_serialize (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GVariant *var     = ncm_serialize_to_variant (ser, G_OBJECT (test->ncdata));

  nc_data_cluster_ncount_free (test->ncdata);
  test->ncdata = NC_DATA_CLUSTER_NCOUNT (ncm_serialize_from_variant (ser, var));

  g_variant_unref (var);
  ncm_serialize_free (ser);
}

static void
_test_nc_data_cluster_ncount_saveload (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  nc_data_cluster_ncount_catalog_save (test->ncdata, "test_nc_data_cluster_ncount_tmp.fits", TRUE);
  nc_data_cluster_ncount_free (test->ncdata);

  test->ncdata = nc_data_cluster_ncount_new (test->cad, "NcClusterRedshiftNodist", "NcClusterMassAscaso");
  nc_data_cluster_ncount_catalog_load (test->ncdata, "test_nc_data_cluster_ncount_tmp.fits");
}

void
test_nc_data_cluster_ncount_saveload (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  test_nc_data_cluster_ncount_sanity (test, pdata);
  _test_nc_data_cluster_ncount_saveload (test, pdata);
}

void
test_nc_data_cluster_ncount_m2lnL_true_data_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);

  nc_data_cluster_ncount_true_data (test->ncdata, TRUE);
  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL);

  g_assert (gsl_finite (m2lnL));
}

void
test_nc_data_cluster_ncount_m2lnL_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL);

  g_assert (gsl_finite (m2lnL));
}

void
test_nc_data_cluster_ncount_m2lnL_binned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);
  test_nc_data_cluster_ncount_bin (test, pdata);

  nc_data_cluster_ncount_set_binned (test->ncdata, TRUE);
  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL);

  g_assert (gsl_finite (m2lnL));
}

void
test_nc_data_cluster_ncount_m2lnL_serialize_true_data_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL0 = 0.0;
  gdouble m2lnL1 = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);

  nc_data_cluster_ncount_true_data (test->ncdata, TRUE);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL0);
  g_assert (gsl_finite (m2lnL0));

  test_nc_data_cluster_ncount_serialize (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL1);
  g_assert (gsl_finite (m2lnL1));

  ncm_assert_cmpdouble_e (m2lnL0, ==, m2lnL1, 1.0e-15, 0.0);
}

void
test_nc_data_cluster_ncount_m2lnL_serialize_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL0 = 0.0;
  gdouble m2lnL1 = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL0);
  g_assert (gsl_finite (m2lnL0));

  test_nc_data_cluster_ncount_serialize (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL1);
  g_assert (gsl_finite (m2lnL1));

  ncm_assert_cmpdouble_e (m2lnL0, ==, m2lnL1, 1.0e-15, 0.0);
}

void
test_nc_data_cluster_ncount_m2lnL_serialize_binned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL0 = 0.0;
  gdouble m2lnL1 = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);
  test_nc_data_cluster_ncount_bin (test, pdata);

  nc_data_cluster_ncount_set_binned (test->ncdata, TRUE);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL0);
  g_assert (gsl_finite (m2lnL0));

  test_nc_data_cluster_ncount_serialize (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL1);
  g_assert (gsl_finite (m2lnL1));

  ncm_assert_cmpdouble_e (m2lnL0, ==, m2lnL1, 1.0e-15, 0.0);
}

void
test_nc_data_cluster_ncount_m2lnL_saveload_true_data_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL0 = 0.0;
  gdouble m2lnL1 = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);

  nc_data_cluster_ncount_true_data (test->ncdata, TRUE);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL0);
  g_assert (gsl_finite (m2lnL0));

  _test_nc_data_cluster_ncount_saveload (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL1);
  g_assert (gsl_finite (m2lnL1));

  ncm_assert_cmpdouble_e (m2lnL0, ==, m2lnL1, 1.0e-12, 0.0);
}

void
test_nc_data_cluster_ncount_m2lnL_saveload_unbinned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL0 = 0.0;
  gdouble m2lnL1 = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL0);
  g_assert (gsl_finite (m2lnL0));

  _test_nc_data_cluster_ncount_saveload (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL1);
  g_assert (gsl_finite (m2lnL1));

  ncm_assert_cmpdouble_e (m2lnL0, ==, m2lnL1, 1.0e-12, 0.0);
}

void
test_nc_data_cluster_ncount_m2lnL_saveload_binned (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  gdouble m2lnL0 = 0.0;
  gdouble m2lnL1 = 0.0;

  test_nc_data_cluster_ncount_sanity (test, pdata);
  test_nc_data_cluster_ncount_bin (test, pdata);

  nc_data_cluster_ncount_set_binned (test->ncdata, TRUE);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL0);
  g_assert (gsl_finite (m2lnL0));

  _test_nc_data_cluster_ncount_saveload (test, pdata);

  ncm_data_m2lnL_val (NCM_DATA (test->ncdata), test->mset, &m2lnL1);
  g_assert (gsl_finite (m2lnL1));

  ncm_assert_cmpdouble_e (m2lnL0, ==, m2lnL1, 1.0e-12, 0.0);
}

void
test_nc_data_cluster_ncount_traps (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/nc/data_cluster_ncount/invalid/test/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_nc_data_cluster_ncount_invalid_test (TestNcDataClusterNCount *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

