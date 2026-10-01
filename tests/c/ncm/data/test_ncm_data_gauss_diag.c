/***************************************************************************
 *            test_ncm_data_gauss_diag.c
 *
 *  Tue Sep 29 01:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#define TEST_NP 5

static const gdouble test_x[TEST_NP]      = {0.0, 0.25, 0.5, 0.75, 1.0};
static const gdouble test_sigma0[TEST_NP] = {0.1, 0.2, 0.3, 0.2, 0.1};

/*
 * A NcmDataGaussDiag with mean a + b x_i and standard deviations s sigma0_i, where
 * (a, b, s) is the NcmModelMVND mean; sigma_func returns TRUE only when s changed.
 */

#define TEST_TYPE_DIAG (test_diag_get_type ())
G_DECLARE_FINAL_TYPE (TestDiag, test_diag, TEST, DIAG, NcmDataGaussDiag)

struct _TestDiag
{
  NcmDataGaussDiag parent_instance;
  gdouble last_s;
};

G_DEFINE_TYPE (TestDiag, test_diag, NCM_TYPE_DATA_GAUSS_DIAG)

static void
test_diag_init (TestDiag *td)
{
  td->last_s = GSL_NAN;
}

static gdouble
_test_diag_param (NcmMSet *mset, const guint i)
{
  return ncm_model_orig_vparam_get (ncm_mset_peek (mset, ncm_model_mvnd_id ()), NCM_MODEL_MVND_MEAN, i);
}

static void
_test_diag_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  const gdouble a = _test_diag_param (mset, 0);
  const gdouble b = _test_diag_param (mset, 1);
  guint i;

  for (i = 0; i < TEST_NP; i++)
    ncm_vector_set (vp, i, a + b * test_x[i]);
}

static gboolean
_test_diag_sigma_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *sigma)
{
  TestDiag *td    = TEST_DIAG (diag);
  const gdouble s = _test_diag_param (mset, 2);
  guint i;

  if (s == td->last_s)
    return FALSE;

  td->last_s = s;

  for (i = 0; i < TEST_NP; i++)
    ncm_vector_set (sigma, i, s * test_sigma0[i]);

  return TRUE;
}

static void
test_diag_class_init (TestDiagClass *klass)
{
  NCM_DATA_CLASS (klass)->name                  = "Test diagonal data";
  NCM_DATA_GAUSS_DIAG_CLASS (klass)->mean_func  = &_test_diag_mean_func;
  NCM_DATA_GAUSS_DIAG_CLASS (klass)->sigma_func = &_test_diag_sigma_func;
}

typedef struct _TestNcmDataGaussDiag
{
  NcmMSet *mset;
  NcmData *data;
  NcmDataGaussDiag *diag;
  NcmRNG *rng;
} TestNcmDataGaussDiag;

static void
_test_set_params (TestNcmDataGaussDiag *test, const gdouble a, const gdouble b, const gdouble s)
{
  NcmModel *mvnd = ncm_mset_peek (test->mset, ncm_model_mvnd_id ());

  ncm_model_orig_vparam_set (mvnd, NCM_MODEL_MVND_MEAN, 0, a);
  ncm_model_orig_vparam_set (mvnd, NCM_MODEL_MVND_MEAN, 1, b);
  ncm_model_orig_vparam_set (mvnd, NCM_MODEL_MVND_MEAN, 2, s);
}

/* The profiled chi^2 of the current data, computed directly. */
static gdouble
_test_profiled_chi2 (TestNcmDataGaussDiag *test)
{
  NcmVector *mu    = ncm_vector_new (TEST_NP);
  NcmVector *y     = ncm_data_gauss_diag_peek_mean (test->diag);
  NcmVector *sigma = ncm_data_gauss_diag_peek_std (test->diag);
  gdouble wr2      = 0.0;
  gdouble wr       = 0.0;
  gdouble wt       = 0.0;
  guint i;

  _test_diag_mean_func (test->diag, test->mset, mu);

  for (i = 0; i < TEST_NP; i++)
  {
    const gdouble w_i = 1.0 / gsl_pow_2 (ncm_vector_get (sigma, i));
    const gdouble r_i = ncm_vector_get (mu, i) - ncm_vector_get (y, i);

    wr2 += w_i * r_i * r_i;
    wr  += w_i * r_i;
    wt  += w_i;
  }

  ncm_vector_free (mu);

  return wr2 - wr * wr / wt;
}

void test_ncm_data_gauss_diag_new (TestNcmDataGaussDiag *test, gconstpointer pdata);
void test_ncm_data_gauss_diag_free (TestNcmDataGaussDiag *test, gconstpointer pdata);

void test_ncm_data_gauss_diag_profiled_uh (TestNcmDataGaussDiag *test, gconstpointer pdata);
void test_ncm_data_gauss_diag_profiled_offset (TestNcmDataGaussDiag *test, gconstpointer pdata);
void test_ncm_data_gauss_diag_sigma_func (TestNcmDataGaussDiag *test, gconstpointer pdata);
void test_ncm_data_gauss_diag_sigma_prop (TestNcmDataGaussDiag *test, gconstpointer pdata);
void test_ncm_data_gauss_diag_resize (TestNcmDataGaussDiag *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/data_gauss_diag/profiled/uh", TestNcmDataGaussDiag, NULL, &test_ncm_data_gauss_diag_new, &test_ncm_data_gauss_diag_profiled_uh, &test_ncm_data_gauss_diag_free);
  g_test_add ("/ncm/data_gauss_diag/profiled/offset", TestNcmDataGaussDiag, NULL, &test_ncm_data_gauss_diag_new, &test_ncm_data_gauss_diag_profiled_offset, &test_ncm_data_gauss_diag_free);
  g_test_add ("/ncm/data_gauss_diag/sigma_func", TestNcmDataGaussDiag, NULL, &test_ncm_data_gauss_diag_new, &test_ncm_data_gauss_diag_sigma_func, &test_ncm_data_gauss_diag_free);
  g_test_add ("/ncm/data_gauss_diag/sigma_prop", TestNcmDataGaussDiag, NULL, &test_ncm_data_gauss_diag_new, &test_ncm_data_gauss_diag_sigma_prop, &test_ncm_data_gauss_diag_free);
  g_test_add ("/ncm/data_gauss_diag/resize", TestNcmDataGaussDiag, NULL, &test_ncm_data_gauss_diag_new, &test_ncm_data_gauss_diag_resize, &test_ncm_data_gauss_diag_free);

  g_test_run ();
}

void
test_ncm_data_gauss_diag_new (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (3);

  test->mset = ncm_mset_new (mvnd, NULL, NULL);
  test->diag = g_object_new (TEST_TYPE_DIAG, "n-points", TEST_NP, "w-mean", TRUE, NULL);
  test->data = NCM_DATA (test->diag);
  test->rng  = ncm_rng_seeded_new (NULL, 1);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_param_set_ftype (test->mset, ncm_model_mvnd_id (), 2, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
  _test_set_params (test, 0.3, 1.7, 1.0);

  ncm_data_resample (test->data, test->mset, test->rng);

  ncm_model_mvnd_free (mvnd);
}

void
test_ncm_data_gauss_diag_free (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  ncm_data_free (test->data);
  ncm_mset_free (test->mset);
  ncm_rng_free (test->rng);
}

void
test_ncm_data_gauss_diag_profiled_uh (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  NcmMatrix *H = ncm_matrix_new (2, TEST_NP);
  gdouble wt   = 0.0;
  gdouble w[TEST_NP];
  guint i, a, b;

  /* The Jacobian of the mean with respect to (a, b). */
  for (i = 0; i < TEST_NP; i++)
  {
    ncm_matrix_set (H, 0, i, 1.0);
    ncm_matrix_set (H, 1, i, test_x[i]);
    w[i] = 1.0 / gsl_pow_2 (test_sigma0[i]);
    wt  += w[i];
  }

  ncm_data_inv_cov_UH (test->data, test->mset, H);

  /* H H^T is the Fisher matrix of the profiled likelihood, J (W - w w^T / wt) J^T;
   * the two computations differ only by rounding: at most 5.4e-17 wt measured. */
  for (a = 0; a < 2; a++)
  {
    for (b = 0; b < 2; b++)
    {
      gdouble F_ab  = 0.0;
      gdouble E_ab  = 0.0;
      gdouble J_a_w = 0.0;
      gdouble J_b_w = 0.0;

      for (i = 0; i < TEST_NP; i++)
      {
        const gdouble J_ai = (a == 0) ? 1.0 : test_x[i];
        const gdouble J_bi = (b == 0) ? 1.0 : test_x[i];

        F_ab  += ncm_matrix_get (H, a, i) * ncm_matrix_get (H, b, i);
        E_ab  += J_ai * w[i] * J_bi;
        J_a_w += J_ai * w[i];
        J_b_w += J_bi * w[i];
      }

      E_ab -= J_a_w * J_b_w / wt;

      g_assert_cmpfloat (fabs (F_ab - E_ab), <=, 1.0e-14 * wt);
    }
  }

  ncm_matrix_free (H);
}

void
test_ncm_data_gauss_diag_profiled_offset (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  NcmVector *f_true      = ncm_vector_new (TEST_NP);
  NcmMatrix *IM          = NULL;
  NcmVector *delta_theta = NULL;
  guint i;

  /* A true mean that differs from the model mean by a constant: under the profiled
   * likelihood there is no bias. */
  _test_diag_mean_func (test->diag, test->mset, f_true);
  ncm_vector_add_constant (f_true, 0.7);

  ncm_data_fisher_matrix_bias (test->data, test->mset, f_true, &IM, &delta_theta);

  /* Rounding of the projection: at most 1.4e-14 measured. */
  for (i = 0; i < 2; i++)
    g_assert_cmpfloat (fabs (ncm_vector_get (delta_theta, i)), <=, 1.0e-12);

  /* The constant direction carries no information (2.4e-22 measured, against
   * F_bb = 53.1). */
  g_assert_cmpfloat (fabs (ncm_matrix_get (IM, 0, 0)), <=, 1.0e-12);

  ncm_vector_free (f_true);
  ncm_vector_free (delta_theta);
  ncm_matrix_free (IM);
}

void
test_ncm_data_gauss_diag_sigma_func (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  gdouble m2lnL;
  gdouble chi2;
  gdouble norma = 0.0;
  guint i;

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);

  /* The change of s is seen first by resample; m2lnL must use the new weights. */
  _test_set_params (test, 0.3, 1.7, 2.0);
  ncm_data_resample (test->data, test->mset, test->rng);
  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);

  chi2 = _test_profiled_chi2 (test);

  for (i = 0; i < TEST_NP; i++)
    norma += ncm_c_ln2pi () + 2.0 * log (2.0 * test_sigma0[i]);

  /* The same sums in another order: rounding only, at most 1.9e-16 relative measured. */
  g_assert_cmpfloat (fabs (m2lnL - (chi2 + norma)), <=, 1.0e-14 * fabs (m2lnL));
}

void
test_ncm_data_gauss_diag_sigma_prop (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  NcmVector *sigma = ncm_vector_new (TEST_NP);
  NcmVector *f     = ncm_vector_new (TEST_NP);
  gdouble m2lnL;
  gdouble norma = 0.0;
  guint i;

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);

  /* Setting sigma updates the weights (the parameter s is unchanged). */
  for (i = 0; i < TEST_NP; i++)
  {
    ncm_vector_set (sigma, i, 3.0 * test_sigma0[i]);
    norma += ncm_c_ln2pi () + 2.0 * log (3.0 * test_sigma0[i]);
  }

  g_object_set (test->diag, "sigma", sigma, NULL);

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  ncm_data_leastsquares_f (test->data, test->mset, f);

  g_assert_cmpfloat (fabs (m2lnL - (_test_profiled_chi2 (test) + norma)), <=, 1.0e-14 * fabs (m2lnL));
  g_assert_cmpfloat (fabs (ncm_vector_dot (f, f) - _test_profiled_chi2 (test)), <=, 1.0e-14 * fabs (m2lnL));

  ncm_vector_free (sigma);
  ncm_vector_free (f);
}

void
test_ncm_data_gauss_diag_resize (TestNcmDataGaussDiag *test, gconstpointer pdata)
{
  gdouble m2lnL;

  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);

  /* A new size frees the weights; the next evaluation rebuilds them. */
  ncm_data_gauss_diag_set_size (test->diag, TEST_NP + 1);
  ncm_data_gauss_diag_set_size (test->diag, TEST_NP);
  TEST_DIAG (test->diag)->last_s = GSL_NAN;
  ncm_data_resample (test->data, test->mset, test->rng);
  ncm_data_m2lnL_val (test->data, test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));
}

