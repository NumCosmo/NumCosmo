/***************************************************************************
 *            test_ncm_lh_ratio1d.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_lh_ratio1d.c
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
#include <gsl/gsl_cdf.h>

typedef struct _TestNcmLHRatio1d
{
  NcmDataGaussCovMVND *data_mvnd;
  NcmModelMVND *model;
  NcmFit *fit;
} TestNcmLHRatio1d;

static void
test_ncm_lh_ratio1d_new (TestNcmLHRatio1d *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 20260930);
  NcmMSet *mset     = NULL;
  NcmDataset *dset  = NULL;
  NcmLikelihood *lh = NULL;

  test->data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-2, 1.0, 50.0, -1.0, 1.0, rng);
  test->model     = ncm_model_mvnd_new (2);
  mset            = ncm_mset_new (NCM_MODEL (test->model), NULL, NULL);
  dset            = ncm_dataset_new_list (test->data_mvnd, NULL);
  lh              = ncm_likelihood_new (dset);

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  test->fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);

  g_assert_true (ncm_fit_run (test->fit, NCM_FIT_RUN_MSGS_NONE));
  ncm_fit_obs_fisher (test->fit);

  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_rng_free (rng);
}

static void
test_ncm_lh_ratio1d_free (TestNcmLHRatio1d *test, gconstpointer pdata)
{
  ncm_fit_free (test->fit);
  ncm_model_mvnd_free (test->model);
  ncm_data_gauss_cov_mvnd_free (test->data_mvnd);
}

/*
 * The profile of -2 ln L over the other parameter is (theta_i - y_i)^2 / C_ii, so the
 * bounds at confidence level cl are the offsets -/+ sqrt (chi2_1 (cl) C_ii).
 */
static void
_test_ncm_lh_ratio1d_assert_bounds (NcmLHRatio1d *lhr1d, TestNcmLHRatio1d *test, guint i, gdouble clevel)
{
  NcmMatrix *C        = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (test->data_mvnd));
  const gdouble delta = sqrt (gsl_cdf_chisq_Qinv (1.0 - clevel, 1.0) * ncm_matrix_get (C, i, i));
  gdouble lb, ub;

  ncm_lh_ratio1d_find_bounds (lhr1d, clevel, NCM_FIT_RUN_MSGS_NONE, &lb, &ub);

  ncm_assert_cmpdouble_e (lb, ==, -delta, 1.0e-8, 0.0);
  ncm_assert_cmpdouble_e (ub, ==, delta, 1.0e-8, 0.0);
}

static void
test_ncm_lh_ratio1d_bounds (TestNcmLHRatio1d *test, gconstpointer pdata)
{
  NcmMSetPIndex *pi0  = ncm_mset_pindex_new (ncm_model_mvnd_id (), 0);
  NcmMSetPIndex *pi1  = ncm_mset_pindex_new (ncm_model_mvnd_id (), 1);
  NcmLHRatio1d *lhr1d = ncm_lh_ratio1d_new (test->fit, pi0);

  _test_ncm_lh_ratio1d_assert_bounds (lhr1d, test, 0, ncm_c_stats_1sigma ());
  _test_ncm_lh_ratio1d_assert_bounds (lhr1d, test, 0, ncm_c_stats_2sigma ());

  ncm_lh_ratio1d_set_pindex (lhr1d, pi1);
  _test_ncm_lh_ratio1d_assert_bounds (lhr1d, test, 1, ncm_c_stats_1sigma ());

  ncm_lh_ratio1d_free (lhr1d);
  ncm_mset_pindex_free (pi0);
  ncm_mset_pindex_free (pi1);
}

static void
test_ncm_lh_ratio1d_upper_bound (TestNcmLHRatio1d *test, gconstpointer pdata)
{
  /* An upper parameter bound 0.5 sigma above the best fit ends the interval there. */
  NcmMSetPIndex *pi0  = ncm_mset_pindex_new (ncm_model_mvnd_id (), 0);
  NcmMSet *mset       = ncm_fit_peek_mset (test->fit);
  const gdouble bf    = ncm_mset_param_get (mset, ncm_model_mvnd_id (), 0);
  const gdouble sd    = ncm_fit_covar_sd (test->fit, ncm_model_mvnd_id (), 0);
  NcmLHRatio1d *lhr1d = NULL;
  gdouble lb, ub;

  ncm_model_param_set_upper_bound (NCM_MODEL (test->model), 0, bf + 0.5 * sd);
  lhr1d = ncm_lh_ratio1d_new (test->fit, pi0);

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*reaches the upper bound*");
  ncm_lh_ratio1d_find_bounds (lhr1d, ncm_c_stats_1sigma (), NCM_FIT_RUN_MSGS_NONE, &lb, &ub);
  g_test_assert_expected_messages ();

  ncm_assert_cmpdouble_e (ub, ==, 0.5 * sd, 1.0e-12, 0.0);
  g_assert_cmpfloat (lb, <, 0.0);

  ncm_lh_ratio1d_free (lhr1d);
  ncm_mset_pindex_free (pi0);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/lh_ratio1d/bounds", TestNcmLHRatio1d, NULL,
              &test_ncm_lh_ratio1d_new,
              &test_ncm_lh_ratio1d_bounds,
              &test_ncm_lh_ratio1d_free);
  g_test_add ("/ncm/lh_ratio1d/upper_bound", TestNcmLHRatio1d, NULL,
              &test_ncm_lh_ratio1d_new,
              &test_ncm_lh_ratio1d_upper_bound,
              &test_ncm_lh_ratio1d_free);

  g_test_run ();
}

