/***************************************************************************
 *            test_ncm_mset_catalog.c
 *
 *  Wed February 07 10:36:55 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <sandro@isoftware.com.br>
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

#include "ncm_data_gauss_cov_mvnd.h"
#include "ncm_model_mvnd.h"

typedef struct _TestNcmMSetCatalog
{
  gint dim;
  NcmRNG *rng;
  NcmDataGaussCovMVND *data_mvnd;
  NcmMSetCatalog *mcat;
  guint ntests;
} TestNcmMSetCatalog;

void test_ncm_mset_catalog_new (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_traps (TestNcmMSetCatalog *test, gconstpointer pdata);

void test_ncm_mset_catalog_free (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_mean (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_cov (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_norma (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_vol (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_invalid_run (TestNcmMSetCatalog *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/mset/catalog/mean", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_mean,
              &test_ncm_mset_catalog_free);

  g_test_add ("/ncm/mset/catalog/cov", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_cov,
              &test_ncm_mset_catalog_free);

  g_test_add ("/ncm/mset/catalog/norma", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_norma,
              &test_ncm_mset_catalog_free);
  
  g_test_add ("/ncm/mset/catalog/vol", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_vol,
              &test_ncm_mset_catalog_free);
  
  g_test_add ("/ncm/mset/catalog/traps", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_traps,
              &test_ncm_mset_catalog_free);
  
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_add ("/ncm/mset/catalog/invalid/run/subprocess", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_invalid_run,
              &test_ncm_mset_catalog_free);
#endif
  g_test_run ();
}

void
test_ncm_mset_catalog_new (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  const gint dim                 = test->dim = g_test_rand_int_range (2, 20);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 1.0e-1, 1.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmMSetCatalog *mcat;

  /*g_object_set (data_mvnd,  "use-norma", FALSE, NULL);*/
  /*ncm_matrix_set_identity (NCM_DATA_GAUSS_COV (data_mvnd)->cov);*/
  
  ncm_mset_param_set_vector (mset, NCM_DATA_GAUSS_COV (data_mvnd)->y);
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE, 
                               "m2lnL", "-2\\ln(L)", 
                               NULL);

  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  
  test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd);
  test->mcat      = ncm_mset_catalog_ref (mcat);
  test->rng       = rng;

  g_assert (NCM_IS_MSET_CATALOG (test->mcat));

  ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_mset_clear (&mset);
  ncm_mset_catalog_clear (&mcat);
}

void
test_ncm_mset_catalog_free (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_mset_catalog_free, test->mcat);
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (test->data_mvnd));
  NCM_TEST_FREE (ncm_rng_free, test->rng);
}

void
test_ncm_mset_catalog_mean (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt = g_test_rand_int_range (100000, 1000000); 
  gint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_data_resample (data, mset, test->rng);
    ncm_mset_catalog_add_from_vector_array (test->mcat, cov->y, &m2lnL);
  }

  if (FALSE)
  {
    ncm_cfg_msg_sepa ();
    ncm_mset_catalog_estimate_autocorrelation_tau (test->mcat, FALSE);
    ncm_mset_catalog_log_current_stats (test->mcat);
    ncm_mset_params_log_vals (mset);
  }

  {
    NcmVector *mean = NULL;
    NcmVector *p    = ncm_vector_new (test->dim);
    
    ncm_mset_catalog_get_mean (test->mcat, &mean);
    ncm_mset_param_get_vector (mset, p);

    ncm_vector_cmp (p, mean);

    g_assert_cmpfloat (ncm_vector_get_max (p), <, 1.0e-3);
    /*printf ("MAX % 22.15e\n", ncm_vector_get_max (p));*/
      
    ncm_vector_clear (&mean);
    ncm_vector_clear (&p);
  }
  
}

void
test_ncm_mset_catalog_cov (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  NcmMatrix *data_cov  = cov->cov;
  const guint nt = g_test_rand_int_range (100000, 1000000); 
  gint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_data_resample (data, mset, test->rng);
    ncm_mset_catalog_add_from_vector_array (test->mcat, cov->y, &m2lnL);
  }

  {
    NcmMatrix *cat_cov = NULL;
    ncm_mset_catalog_get_covar (test->mcat, &cat_cov);

    g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, 1.0e-1);

    ncm_matrix_norma_diag (data_cov, data_cov);
    ncm_matrix_norma_diag (cat_cov, cat_cov);

    g_assert_cmpfloat (ncm_matrix_cmp (cat_cov, data_cov, 1.0), <, 1.0e-1);

    if (FALSE)
    {
      ncm_cfg_msg_sepa ();
      ncm_matrix_log_vals (cat_cov,  "# CAT  COV: ", "% 12.5g");
      ncm_matrix_log_vals (data_cov, "# DATA COV: ", "% 12.5g");

      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cat_cov, data_cov, 0.0));
      printf ("# WDIFFD  : % 22.15e\n", ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0));

      ncm_matrix_sub (cat_cov, data_cov);
      ncm_matrix_div_elements (cat_cov, data_cov);

      ncm_matrix_log_vals (cat_cov,  "# CMP     : ", "% 12.5e");
    }
  }
}

void
test_ncm_mset_catalog_norma (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (100000, 1000000);
  gint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, cov->y, &m2lnL);
  }

  g_assert_cmpfloat (fabs (ncm_mset_catalog_get_post_lnnorm (test->mcat) / test->dim), <, 0.2);    
}

void
test_ncm_mset_catalog_vol (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (100000, 1000000);
  gdouble glnvol;
  gint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, cov->y, &m2lnL);
  }

  {
    const gdouble lnevol = ncm_mset_catalog_get_post_lnvol (test->mcat, 0.6827, &glnvol);
    ncm_assert_cmpdouble_e (lnevol, ==, glnvol, 0.2, 0.0);
  }
/*
  printf ("# DIM %d LNNORMA = % 22.15g\n", test->dim, ncm_mset_catalog_get_post_lnnorm (test->mcat));
  printf ("# DIM %d VOL1SIG = % 22.15g ", test->dim, ncm_mset_catalog_get_post_lnvol (test->mcat, 0.6827, &glnvol));
  printf ("% 22.15g\n", glnvol);  
*/
}
 

#if GLIB_CHECK_VERSION(2,38,0)
void
test_ncm_mset_catalog_traps (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/invalid/run/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}
#else
void
test_ncm_mset_catalog_traps (TestNcmMSetCatalog *test, gconstpointer pdata)
{
}
#endif

void
test_ncm_mset_catalog_invalid_run (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}


