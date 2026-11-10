/***************************************************************************
 *            test_ncm_fit_mcbs.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_fit_mcbs.c
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

static void test_ncm_fit_mcbs_run (void);
static void test_ncm_fit_mcbs_filename (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/fit_mcbs/run", &test_ncm_fit_mcbs_run);
  g_test_add_func ("/ncm/fit_mcbs/filename", &test_ncm_fit_mcbs_filename);

  g_test_run ();
}

static void
test_ncm_fit_mcbs_run (void)
{
  /* One catalog row, the bootstrap mean, per realization, and the fit gets their
   * mean and covariance. */
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20260930);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (2);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  NcmFit *fit;
  NcmFitMCBS *mcbs;
  NcmMSetCatalog *mcat;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  fit  = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  mcbs = ncm_fit_mcbs_new (fit);
  ncm_fit_mcbs_set_rng (mcbs, rng);

  ncm_fit_mcbs_run (mcbs, NULL, 4, 20, NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, NCM_FIT_RUN_MSGS_NONE, FALSE, NULL);

  mcat = ncm_fit_mcbs_get_catalog (mcbs);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 4);
  g_assert_true (ncm_fit_state_has_covar (ncm_fit_peek_state (fit)));

  ncm_mset_catalog_free (mcat);
  ncm_fit_mcbs_free (mcbs);
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
}

static void
test_ncm_fit_mcbs_filename (void)
{
  /*
   * The catalog on a FITS storage: the bootstrap catalogs are written next to its file
   * under the name without the .fits extension. A second run adds its realizations to
   * the same catalog.
   */
  gchar *tmp_dir                 = g_dir_make_tmp ("tmp_test_ncm_fit_mcbs_XXXXXX", NULL);
  gchar *base                    = g_build_filename (tmp_dir, "mcbs", NULL);
  gchar *filename                = g_strdup_printf ("%s.fits", base);
  gchar *bstrap0                 = g_strdup_printf ("%s-bstrap-%06d.fits", base, 0);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20261101);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (2);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  GError *error                  = NULL;
  NcmFit *fit;
  NcmFitMCBS *mcbs;
  NcmMSetCatalog *mcat;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  fit  = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  mcbs = ncm_fit_mcbs_new (fit);
  ncm_fit_mcbs_set_rng (mcbs, rng);

  {
    NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (filename, FALSE, NULL);
    NcmMSetCatalogStorage *storage_prop     = NULL;

    ncm_fit_mcbs_set_storage (mcbs, NCM_MSET_CATALOG_STORAGE (storage_fits));
    g_object_get (mcbs, "storage", &storage_prop, NULL);
    g_assert_true (storage_prop == NCM_MSET_CATALOG_STORAGE (storage_fits));
    ncm_mset_catalog_storage_free (storage_prop);
    ncm_mset_catalog_storage_fits_free (storage_fits);
  }

  ncm_fit_mcbs_run (mcbs, NULL, 2, 5, NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, NCM_FIT_RUN_MSGS_NONE, FALSE, &error);
  g_assert_no_error (error);

  mcat = ncm_fit_mcbs_get_catalog (mcbs);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 2);
  g_assert_cmpstr (ncm_mset_catalog_peek_filename (mcat), ==, filename);
  ncm_mset_catalog_free (mcat);
  g_assert_true (g_file_test (bstrap0, G_FILE_TEST_EXISTS));

  ncm_fit_mcbs_run (mcbs, NULL, 2, 5, NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, NCM_FIT_RUN_MSGS_NONE, FALSE, &error);
  g_assert_no_error (error);

  mcat = ncm_fit_mcbs_get_catalog (mcbs);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 4);
  ncm_mset_catalog_free (mcat);

  ncm_fit_mcbs_free (mcbs);
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);

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
  g_free (bstrap0);
  g_free (filename);
  g_free (base);
  g_free (tmp_dir);
}

