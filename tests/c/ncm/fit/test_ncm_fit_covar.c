/***************************************************************************
 *            test_ncm_fit_covar.c
 *
 *  Wed September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_fit_covar.c
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

static void
test_ncm_fit_covar_indefinite (void)
{
  /*
   * -2 ln L = 0.1 [100 (x2 - x1^2)^2 + (1 - x1)^2] has, at (0.3, 1), the Hessian
   * H = [[-29, -12], [-12, 20]]: indefinite, so the covariance takes the LU fallback.
   * It must still be (H / 2)^-1, and the log-determinant -ln |det (H / 2)|.
   */
  NcmModelRosenbrock *mrb = ncm_model_rosenbrock_new ();
  NcmDataRosenbrock *drb  = ncm_data_rosenbrock_new ();
  NcmMSet *mset           = ncm_mset_new (NCM_MODEL (mrb), NULL, NULL);
  NcmDataset *dset        = ncm_dataset_new_list (drb, NULL);
  NcmLikelihood *lh       = ncm_likelihood_new (dset);
  NcmFit *fit;
  const gdouble F[2][2] = {
    {
      -14.5, -6.0
    }, {
      -6.0, 10.0
    }
  };
  const gdouble det     = F[0][0] * F[1][1] - F[0][1] * F[1][0];
  const gdouble C[2][2] = {
    {
      F[1][1] / det, -F[0][1] / det
    }, {
      -F[1][0] / det, F[0][0] / det
    }
  };
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  fit = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);

  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, 0.3);
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2, 1.0);

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*covariance matrix not positive definite*");
  ncm_fit_obs_fisher (fit);
  g_test_assert_expected_messages ();

  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      ncm_assert_cmpdouble_e (ncm_matrix_get (ncm_fit_state_peek_covar (ncm_fit_peek_state (fit)), i, j), ==, C[i][j], 1.0e-8, 0.0);

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*covariance matrix not positive definite*");
  ncm_assert_cmpdouble_e (ncm_fit_numdiff_m2lnL_lndet_covar (fit), ==, -log (fabs (det)), 1.0e-8, 0.0);
  g_test_assert_expected_messages ();

  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_data_rosenbrock_free (drb);
  ncm_model_rosenbrock_free (mrb);
}

static void
test_ncm_fit_covar_ls_J_accurate (void)
{
  /* One prior on top of d data points: J is (d + 1) x d, not square. The accurate
   * Jacobian must agree with the forward one to the latter's error. */
  const guint d                  = 3;
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20260930);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (d, 1.0e-2, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (d);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  NcmMatrix *J_fo                = ncm_matrix_new (d + 1, d);
  NcmMatrix *J_ac                = ncm_matrix_new (d + 1, d);
  NcmFit *fit;
  guint i, j;

  ncm_likelihood_priors_take (lh, NCM_PRIOR (ncm_prior_gauss_param_new (NCM_MODEL (model), 0, 0.3, 0.2)));
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);

  for (i = 0; i < d; i++)
    ncm_model_orig_param_set (NCM_MODEL (model), i, 0.1 * (i + 1.0));

  fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  ncm_fit_ls_J (fit, J_fo);
  ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_ACCURATE);
  ncm_fit_ls_J (fit, J_ac);

  for (i = 0; i < d + 1; i++)
    for (j = 0; j < d; j++)
      ncm_assert_cmpdouble_e (ncm_matrix_get (J_ac, i, j), ==, ncm_matrix_get (J_fo, i, j), 1.0e-6, 1.0e-6);

  ncm_fit_free (fit);
  ncm_matrix_free (J_fo);
  ncm_matrix_free (J_ac);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
}

static void
test_ncm_fit_covar_sub_fit_accurate (void)
{
  /*
   * The fit varies theta0 and its sub-fit profiles theta1 of a Gaussian with
   * covariance C. The profile of -2 ln L is (theta0 - y0)^2 / C00, so its gradient is
   * 2 (theta0 - y0) / C00 and its Fisher matrix 1 / C00. The accurate derivatives must
   * rerun the sub-fit at every point, as the forward and central ones do.
   */
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 20260930);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-2, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (2);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  NcmVector *y                   = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmMatrix *C                   = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmVector *grad                = ncm_vector_new (1);
  const gdouble theta0           = 0.1;
  const gdouble C00              = ncm_matrix_get (C, 0, 0);
  NcmMSet *mset_sub;
  NcmFit *fit, *sub_fit;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_param_set_ftype (mset, ncm_model_mvnd_id (), 1, NCM_PARAM_TYPE_FIXED);
  fit = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_ACCURATE);

  /* The models are shared, the free-parameter maps are not. */
  mset_sub = ncm_mset_shallow_copy (mset, NULL);
  ncm_mset_param_set_ftype (mset_sub, ncm_model_mvnd_id (), 0, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (mset_sub, ncm_model_mvnd_id (), 1, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset_sub);
  sub_fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset_sub, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  ncm_fit_set_sub_fit (fit, sub_fit);

  g_assert_cmpuint (ncm_mset_fparams_len (mset), ==, 1);
  g_assert_cmpuint (ncm_mset_fparam_get_pi (mset, 0)->pid, ==, 0);
  g_assert_cmpuint (ncm_mset_fparams_len (mset_sub), ==, 1);
  g_assert_cmpuint (ncm_mset_fparam_get_pi (mset_sub, 0)->pid, ==, 1);

  /* theta1 away from its profile value. */
  ncm_model_orig_param_set (NCM_MODEL (model), 0, theta0);
  ncm_model_orig_param_set (NCM_MODEL (model), 1, ncm_vector_get (y, 1) + 0.5);

  ncm_fit_m2lnL_grad (fit, grad);
  ncm_assert_cmpdouble_e (ncm_vector_get (grad, 0), ==, 2.0 * (theta0 - ncm_vector_get (y, 0)) / C00, 1.0e-6, 0.0);

  ncm_fit_obs_fisher (fit);
  ncm_assert_cmpdouble_e (ncm_matrix_get (ncm_fit_state_peek_covar (ncm_fit_peek_state (fit)), 0, 0), ==, C00, 1.0e-6, 0.0);

  ncm_fit_free (fit);
  ncm_fit_free (sub_fit);
  ncm_vector_free (grad);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset_sub);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/fit/covar/indefinite", &test_ncm_fit_covar_indefinite);
  g_test_add_func ("/ncm/fit/covar/ls_J_accurate", &test_ncm_fit_covar_ls_J_accurate);
  g_test_add_func ("/ncm/fit/covar/sub_fit_accurate", &test_ncm_fit_covar_sub_fit_accurate);

  g_test_run ();
}

