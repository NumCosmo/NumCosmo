/***************************************************************************
 *            test_ncm_fit_mc.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_fit_mc.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

static void test_ncm_fit_mc_from_model (void);
static void test_ncm_fit_mc_incompatible_file (void);
static void test_ncm_fit_mc_data_file (void);
static void _test_ncm_fit_mc_set_file (NcmFitMC *mc, const gchar *filename);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/fit_mc/from_model", &test_ncm_fit_mc_from_model);
  g_test_add_func ("/ncm/fit_mc/incompatible_file", &test_ncm_fit_mc_incompatible_file);
  g_test_add_func ("/ncm/fit_mc/data_file", &test_ncm_fit_mc_data_file);

  g_test_run ();
}

/* Gives mc a FITS storage on filename. */
static void
_test_ncm_fit_mc_set_file (NcmFitMC *mc, const gchar *filename)
{
  NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (filename, FALSE, NULL);

  ncm_fit_mc_set_storage (mc, NCM_MSET_CATALOG_STORAGE (storage_fits));
  ncm_mset_catalog_storage_fits_free (storage_fits);
}

/*
 * A data file whose description is not the fit's (one free parameter against two):
 * start_run reports the storage's error and the run is not started.
 */
static void
test_ncm_fit_mc_incompatible_file (void)
{
  gchar *tmp_dir                 = g_dir_make_tmp ("tmp_test_ncm_fit_mc_XXXXXX", NULL);
  gchar *filename                = g_build_filename (tmp_dir, "mc.fits", NULL);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20261101);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (2);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  NcmModelMVND *model1           = ncm_model_mvnd_new (1);
  NcmMSet *mset1                 = ncm_mset_new (NCM_MODEL (model1), NULL, NULL);
  const gchar *names[]           = {NCM_MSET_CATALOG_M2LNL_COLNAME, NULL};
  const gchar *symbols[]         = {NCM_MSET_CATALOG_M2LNL_SYMBOL, NULL};
  GError *error                  = NULL;
  NcmFit *fit;
  NcmFitMC *mc;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_param_set_all_ftype (mset1, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset1);

  {
    NcmMSetCatalogDesc *desc                = ncm_mset_catalog_desc_new (mset1, 1, FALSE, names, symbols, 0, NCM_MSET_CATALOG_RTYPE_MC, &error);
    NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (filename, FALSE, &error);

    g_assert_no_error (error);
    ncm_mset_catalog_storage_set_desc (NCM_MSET_CATALOG_STORAGE (storage_fits), desc, &error);
    g_assert_no_error (error);

    ncm_mset_catalog_desc_free (desc);
    ncm_mset_catalog_storage_fits_free (storage_fits);
  }

  fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  mc  = ncm_fit_mc_new (fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
  _test_ncm_fit_mc_set_file (mc, filename);

  ncm_fit_mc_start_run (mc, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_INCOMPATIBLE);
  g_assert_false (ncm_fit_mc_is_running (mc));
  g_clear_error (&error);

  ncm_fit_mc_free (mc);
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_mset_free (mset1);
  ncm_model_mvnd_free (model);
  ncm_model_mvnd_free (model1);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);

  g_unlink (filename);
  g_rmdir (tmp_dir);
  g_free (filename);
  g_free (tmp_dir);
}

/*
 * A run on a data file, continued by a second NcmFitMC on the same file: the second
 * takes the file's RNG; one set with another seed is the storage's RNG error. reset
 * erases the catalog.
 */
static void
test_ncm_fit_mc_data_file (void)
{
  gchar *tmp_dir                 = g_dir_make_tmp ("tmp_test_ncm_fit_mc_XXXXXX", NULL);
  gchar *filename                = g_build_filename (tmp_dir, "mc.fits", NULL);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20261101);
  NcmRNG *rng2                   = ncm_rng_seeded_new (NULL, 20261102);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (2);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  GError *error                  = NULL;
  NcmFit *fit;
  NcmFitMC *mc;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);

  mc = ncm_fit_mc_new (fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_mc_set_rng (mc, rng);
  _test_ncm_fit_mc_set_file (mc, filename);

  ncm_fit_mc_start_run (mc, &error);
  g_assert_no_error (error);
  ncm_fit_mc_run (mc, 10);
  ncm_fit_mc_end_run (mc);
  g_assert_cmpuint (ncm_mset_catalog_len (ncm_fit_mc_peek_catalog (mc)), ==, 10);
  ncm_fit_mc_free (mc);

  mc = ncm_fit_mc_new (fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_mc_set_rng (mc, rng2);
  _test_ncm_fit_mc_set_file (mc, filename);
  ncm_fit_mc_start_run (mc, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG);
  g_assert_false (ncm_fit_mc_is_running (mc));
  g_clear_error (&error);
  ncm_fit_mc_free (mc);

  mc = ncm_fit_mc_new (fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
  _test_ncm_fit_mc_set_file (mc, filename);
  ncm_fit_mc_start_run (mc, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_len (ncm_fit_mc_peek_catalog (mc)), ==, 10);
  g_assert_cmpuint (ncm_rng_get_seed (ncm_mset_catalog_peek_rng (ncm_fit_mc_peek_catalog (mc))), ==, 20261101);
  ncm_fit_mc_run (mc, 20);
  ncm_fit_mc_end_run (mc);
  g_assert_cmpuint (ncm_mset_catalog_len (ncm_fit_mc_peek_catalog (mc)), ==, 20);

  ncm_fit_mc_reset (mc);
  g_assert_cmpuint (ncm_mset_catalog_len (ncm_fit_mc_peek_catalog (mc)), ==, 0);
  ncm_fit_mc_free (mc);

  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_rng_free (rng2);

  {
    GDir *dir = g_dir_open (tmp_dir, 0, NULL);
    const gchar *name;

    g_assert_nonnull (dir);

    while ((name = g_dir_read_name (dir)) != NULL)
    {
      gchar *path = g_build_filename (tmp_dir, name, NULL);

      g_unlink (path);
      g_free (path);
    }

    g_dir_close (dir);
  }

  g_rmdir (tmp_dir);
  g_free (filename);
  g_free (tmp_dir);
}

static void
test_ncm_fit_mc_from_model (void)
{
  /*
   * The best fit of each realization resampled from the fiducial theta is the resampled
   * data vector, so the best fits have mean theta and covariance C, the data covariance:
   * the sample mean within sqrt (C_ii / n) and the sample covariance within
   * sqrt ((C_ii C_jj + C_ij^2) / n), checked at five of those.
   */
  const guint n                  = 4000;
  const gdouble theta[2]         = {0.3, -0.2};
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20260930);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (2);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  NcmMatrix *C                   = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmFit *fit;
  NcmFitMC *mc;
  NcmMSetCatalog *mcat;
  NcmVector *mean = NULL;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_model_orig_param_set (NCM_MODEL (model), 0, theta[0]);
  ncm_model_orig_param_set (NCM_MODEL (model), 1, theta[1]);
  fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);

  mc = ncm_fit_mc_new (fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_mc_set_rng (mc, rng);

  ncm_fit_mc_start_run (mc, NULL);
  ncm_fit_mc_run (mc, n);
  ncm_fit_mc_end_run (mc);

  mcat = ncm_fit_mc_peek_catalog (mc);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, n);

  /* A linear model fits every realization exactly. */
  for (i = 0; i < n; i++)
    g_assert_cmpfloat (fabs (ncm_vector_get (ncm_mset_catalog_peek_row (mcat, i), 0)), <, 1.0e-8);

  ncm_fit_mc_mean_covar (mc);
  ncm_mset_catalog_get_mean (mcat, &mean);

  {
    NcmFitState *fstate = ncm_fit_peek_state (fit);
    NcmVector *fparams  = ncm_fit_state_peek_fparams (fstate);
    NcmMatrix *covar    = ncm_fit_state_peek_covar (fstate);

    g_assert_true (ncm_fit_state_has_covar (fstate));

    for (i = 0; i < 2; i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (fparams, i), ==, ncm_vector_get (mean, i));
      g_assert_cmpfloat (fabs (ncm_vector_get (fparams, i) - theta[i]), <, 5.0 * sqrt (ncm_matrix_get (C, i, i) / n));

      for (j = 0; j < 2; j++)
      {
        const gdouble sd_ij = sqrt ((ncm_matrix_get (C, i, i) * ncm_matrix_get (C, j, j) + gsl_pow_2 (ncm_matrix_get (C, i, j))) / n);

        g_assert_cmpfloat (fabs (ncm_matrix_get (covar, i, j) - ncm_matrix_get (C, i, j)), <, 5.0 * sd_ij);
      }
    }
  }

  ncm_vector_free (mean);
  ncm_fit_mc_free (mc);
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
}

