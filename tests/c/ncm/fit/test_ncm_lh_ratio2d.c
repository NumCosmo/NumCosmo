/***************************************************************************
 *            test_ncm_lh_ratio2d.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_lh_ratio2d.c
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

typedef struct _TestNcmLHRatio2d
{
  NcmDataGaussCovMVND *data_mvnd;
  NcmModelMVND *model;
  NcmFit *fit;
  NcmLHRatio2d *lhr2d;
} TestNcmLHRatio2d;

static void
test_ncm_lh_ratio2d_new (TestNcmLHRatio2d *test, gconstpointer pdata)
{
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, 20260930);
  NcmMSetPIndex *pi0 = ncm_mset_pindex_new (ncm_model_mvnd_id (), 0);
  NcmMSetPIndex *pi1 = ncm_mset_pindex_new (ncm_model_mvnd_id (), 1);
  NcmMSet *mset      = NULL;
  NcmDataset *dset   = NULL;
  NcmLikelihood *lh  = NULL;

  test->data_mvnd = ncm_data_gauss_cov_mvnd_new_full (3, 1.0e-2, 1.0, 50.0, -1.0, 1.0, rng);
  test->model     = ncm_model_mvnd_new (3);
  mset            = ncm_mset_new (NCM_MODEL (test->model), NULL, NULL);
  dset            = ncm_dataset_new_list (test->data_mvnd, NULL);
  lh              = ncm_likelihood_new (dset);

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  test->fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);

  g_assert_true (ncm_fit_run (test->fit, NCM_FIT_RUN_MSGS_NONE));
  ncm_fit_obs_fisher (test->fit);
  test->lhr2d = ncm_lh_ratio2d_new (test->fit, pi0, pi1, 1.0e-5);

  ncm_mset_pindex_free (pi0);
  ncm_mset_pindex_free (pi1);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_rng_free (rng);
}

static void
test_ncm_lh_ratio2d_free (TestNcmLHRatio2d *test, gconstpointer pdata)
{
  ncm_lh_ratio2d_free (test->lhr2d);
  ncm_fit_free (test->fit);
  ncm_model_mvnd_free (test->model);
  ncm_data_gauss_cov_mvnd_free (test->data_mvnd);
}

/*
 * The profile of -2 ln L over the third parameter is d^T M^-1 d, d the offset from the
 * best fit and M the leading 2x2 block of the data covariance, so every border point at
 * confidence level cl has d^T M^-1 d = chi2_2 (cl). Returns the largest relative
 * deviation.
 */
static gdouble
_test_ncm_lh_ratio2d_border_dev (TestNcmLHRatio2d *test, NcmLHRatio2dRegion *rg)
{
  NcmMatrix *C       = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (test->data_mvnd));
  NcmMSet *mset      = ncm_fit_peek_mset (test->fit);
  const gdouble m00  = ncm_matrix_get (C, 0, 0);
  const gdouble m01  = ncm_matrix_get (C, 0, 1);
  const gdouble m11  = ncm_matrix_get (C, 1, 1);
  const gdouble det  = m00 * m11 - m01 * m01;
  const gdouble bf0  = ncm_mset_param_get (mset, ncm_model_mvnd_id (), 0);
  const gdouble bf1  = ncm_mset_param_get (mset, ncm_model_mvnd_id (), 1);
  const gdouble chi2 = gsl_cdf_chisq_Qinv (1.0 - rg->clevel, 2.0);
  gdouble dev        = 0.0;
  guint i;

  for (i = 0; i < rg->np; i++)
  {
    const gdouble d0 = ncm_vector_get (rg->p1, i) - bf0;
    const gdouble d1 = ncm_vector_get (rg->p2, i) - bf1;
    const gdouble q  = (m11 * d0 * d0 - 2.0 * m01 * d0 * d1 + m00 * d1 * d1) / det;

    dev = GSL_MAX (dev, fabs (q / chi2 - 1.0));
  }

  return dev;
}

static void
test_ncm_lh_ratio2d_fisher_border (TestNcmLHRatio2d *test, gconstpointer pdata)
{
  NcmLHRatio2dRegion *rg = ncm_lh_ratio2d_fisher_border (test->lhr2d, ncm_c_stats_1sigma (), 0.0, NCM_FIT_RUN_MSGS_NONE);

  g_assert_cmpuint (rg->np, ==, 601);
  g_assert_cmpfloat (_test_ncm_lh_ratio2d_border_dev (test, rg), <, 1.0e-8);

  ncm_lh_ratio2d_region_free (rg);
}

static void
test_ncm_lh_ratio2d_conf_region (TestNcmLHRatio2d *test, gconstpointer pdata)
{
  NcmLHRatio2dRegion *rg = ncm_lh_ratio2d_conf_region (test->lhr2d, ncm_c_stats_1sigma (), 40.0, NCM_FIT_RUN_MSGS_NONE);

  /* A closed polygon of about 40 points, every one on the border. */
  g_assert_cmpuint (rg->np, >=, 41);
  g_assert_cmpuint (rg->np, <=, 44);
  g_assert_cmpfloat (ncm_vector_get (rg->p1, 0), ==, ncm_vector_get (rg->p1, rg->np - 1));
  g_assert_cmpfloat (ncm_vector_get (rg->p2, 0), ==, ncm_vector_get (rg->p2, rg->np - 1));
  g_assert_cmpfloat (_test_ncm_lh_ratio2d_border_dev (test, rg), <, 1.0e-8);

  ncm_lh_ratio2d_region_free (rg);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/lh_ratio2d/fisher_border", TestNcmLHRatio2d, NULL,
              &test_ncm_lh_ratio2d_new,
              &test_ncm_lh_ratio2d_fisher_border,
              &test_ncm_lh_ratio2d_free);
  g_test_add ("/ncm/lh_ratio2d/conf_region", TestNcmLHRatio2d, NULL,
              &test_ncm_lh_ratio2d_new,
              &test_ncm_lh_ratio2d_conf_region,
              &test_ncm_lh_ratio2d_free);

  g_test_run ();
}

