/***************************************************************************
 *            test_ncm_fftlog.c
 *
 *  Sun September 03 11:45:13 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_fftlog.c
 *
 * Copyright (C) 2017 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_sf_gamma.h>
#include <complex.h>

/* A kernel that declares no bias range: the Gaussian window's coefficients at b = 0. */
#define TEST_TYPE_NCM_FFTLOG_NO_BIAS (test_ncm_fftlog_no_bias_get_type ())
G_DECLARE_FINAL_TYPE (TestNcmFftlogNoBias, test_ncm_fftlog_no_bias, TEST, NCM_FFTLOG_NO_BIAS, NcmFftlog)

struct _TestNcmFftlogNoBias
{
  NcmFftlog parent_instance;
};

G_DEFINE_TYPE (TestNcmFftlogNoBias, test_ncm_fftlog_no_bias, NCM_TYPE_FFTLOG)

static void
test_ncm_fftlog_no_bias_init (TestNcmFftlogNoBias *fftlog)
{
}

static void
_test_ncm_fftlog_no_bias_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
  const gdouble twopi_Lt = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  complex double *Ym     = (complex double *) Ym_0;
  gint i;

  for (i = 0; i < ncm_fftlog_get_full_size (fftlog); i++)
  {
    const complex double x = 0.5 * (1.0 + twopi_Lt * ncm_fftlog_get_mode_index (fftlog, i) * I);
    gsl_sf_result rho, theta;

    gsl_sf_lngamma_complex_e (creal (x), cimag (x), &rho, &theta);
    Ym[i] = 0.5 * cexp (rho.val + I * theta.val);
  }
}

static void
test_ncm_fftlog_no_bias_class_init (TestNcmFftlogNoBiasClass *klass)
{
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  fftlog_class->name       = "test_no_bias";
  fftlog_class->compute_Ym = &_test_ncm_fftlog_no_bias_compute_Ym;
}

static NcmFftlog *
test_ncm_fftlog_no_bias_new (void)
{
  return g_object_new (TEST_TYPE_NCM_FFTLOG_NO_BIAS, "lnr0", 0.0, "lnk0", 0.0, "Lk", 20.0, "N", 250, NULL);
}

typedef struct _TestNcmFftlogK
{
  gsl_function Fk;
  gdouble lnr;
  guint ell;
  guint ntests;
} TestNcmFftlogK;

typedef struct _TestNcmFftlog
{
  NcmFftlog *fftlog;
  gsl_function Fk;
  gsl_function KFk;
  gdouble lnk_i, lnk_f;
  guint ntests;
  TestNcmFftlogK *argK;
} TestNcmFftlog;

void test_ncm_fftlog_tophatwin2_new (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_gausswin2_new (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_sbessel_j_new (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_sbessel_j_bias0_5_new (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_free (TestNcmFftlog *test, gconstpointer pdata);

void test_ncm_fftlog_setget (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_vector (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_calibrate (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_calibrate_fail (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_serialized (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_deriv (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_use_eval_int (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_eval_smooth_padding (TestNcmFftlog *test, gconstpointer pdata);

void test_ncm_fftlog_tophatwin2_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_gausswin2_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_sbessel_j_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_invalid_st (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_invalid_length (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_tophatwin2_truth (void);
void test_ncm_fftlog_smooth_padding_negative_slope (void);
void test_ncm_fftlog_smooth_padding_refine_stable (void);
void test_ncm_fftlog_calibrate_zero (void);
void test_ncm_fftlog_calibrate_max_n_traps (void);
void test_ncm_fftlog_calibrate_max_n_subprocess (void);
void test_ncm_fftlog_calibrate_short_padding_subprocess (void);
void test_ncm_fftlog_get_Ym_keeps_eval (void);
void test_ncm_fftlog_smooth_padding_traps (void);
void test_ncm_fftlog_invalid_smooth_padding (void);
void test_ncm_fftlog_smooth_padding_power_law (void);
void test_ncm_fftlog_odd_full_size (void);
void test_ncm_fftlog_noring_needs_nyquist (void);
void test_ncm_fftlog_bias_tophatwin2_truth (void);
void test_ncm_fftlog_bias_gausswin2_truth (void);
void test_ncm_fftlog_bias_power_law_converges (void);
void test_ncm_fftlog_bias_sbessel_j_truth (void);
void test_ncm_fftlog_sbessel_j_best_lnr0 (void);
void test_ncm_fftlog_bias_best (void);
void test_ncm_fftlog_bias_traps (void);
void test_ncm_fftlog_bias_invalid_range (void);
void test_ncm_fftlog_bias_invalid_unsupported (void);
void test_ncm_fftlog_bias_invalid_end_slopes (void);

typedef struct _TestCases
{
  gchar *name;

  void (*func) (TestNcmFftlog *test, gconstpointer pdata);
} TestCases;

TestCases tests[] = {
  {"setget", &test_ncm_fftlog_setget},
  {"eval", &test_ncm_fftlog_eval},
  {"eval/vector", &test_ncm_fftlog_eval_vector},
  {"eval/calibrate", &test_ncm_fftlog_eval_calibrate},
  {"eval/calibrate/fail", &test_ncm_fftlog_eval_calibrate_fail},
  {"eval/serialized", &test_ncm_fftlog_eval_serialized},
  {"eval/deriv", &test_ncm_fftlog_eval_deriv},
  {"eval/use_eval_int", &test_ncm_fftlog_eval_use_eval_int},
  {"eval/smooth_padding", &test_ncm_fftlog_eval_smooth_padding},
};

TestCases fixtures[] = {
  {"tophatwin2", &test_ncm_fftlog_tophatwin2_new},
  {"gausswin2", &test_ncm_fftlog_gausswin2_new},
  {"sbessel_j", &test_ncm_fftlog_sbessel_j_new},
  {"sbessel_j_bias0_5", &test_ncm_fftlog_sbessel_j_bias0_5_new},
};

#define NUMINT_RELTOL1 1.0e-3
#define NUMINT_RELTOL2 1.0e-8

gint
main (gint argc, gchar *argv[])
{
  const gint nfixtures = sizeof (fixtures) / sizeof (TestCases);
  const gint ntests    = sizeof (tests) / sizeof (TestCases);
  gint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < nfixtures; i++)
  {
    for (j = 0; j < ntests; j++)
    {
      gchar *path = g_strdup_printf ("/ncm/fftlog/%s/%s", fixtures[i].name, tests[j].name);

      g_test_add (path, TestNcmFftlog, NULL,
                  fixtures[i].func,
                  tests[j].func,
                  &test_ncm_fftlog_free);

      g_free (path);
    }
  }

  g_test_add ("/ncm/fftlog/tophatwin2/traps", TestNcmFftlog, NULL,
              &test_ncm_fftlog_tophatwin2_new,
              &test_ncm_fftlog_tophatwin2_traps,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/gausswin2/traps", TestNcmFftlog, NULL,
              &test_ncm_fftlog_gausswin2_new,
              &test_ncm_fftlog_gausswin2_traps,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/sbessel_j/traps", TestNcmFftlog, NULL,
              &test_ncm_fftlog_sbessel_j_new,
              &test_ncm_fftlog_sbessel_j_traps,
              &test_ncm_fftlog_free);


  g_test_add ("/ncm/fftlog/tophatwin2/invalid/length/subprocess", TestNcmFftlog, NULL,
              &test_ncm_fftlog_tophatwin2_new,
              &test_ncm_fftlog_invalid_length,
              &test_ncm_fftlog_free);
  g_test_add ("/ncm/fftlog/tophatwin2/invalid/st/subprocess", TestNcmFftlog, NULL,
              &test_ncm_fftlog_tophatwin2_new,
              &test_ncm_fftlog_invalid_st,
              &test_ncm_fftlog_free);
  g_test_add ("/ncm/fftlog/gausswin2/invalid/st/subprocess", TestNcmFftlog, NULL,
              &test_ncm_fftlog_gausswin2_new,
              &test_ncm_fftlog_invalid_st,
              &test_ncm_fftlog_free);
  g_test_add ("/ncm/fftlog/sbessel_j/invalid/st/subprocess", TestNcmFftlog, NULL,
              &test_ncm_fftlog_sbessel_j_new,
              &test_ncm_fftlog_invalid_st,
              &test_ncm_fftlog_free);

  g_test_add_func ("/ncm/fftlog/tophatwin2/truth", test_ncm_fftlog_tophatwin2_truth);
  g_test_add_func ("/ncm/fftlog/smooth_padding/negative_slope", test_ncm_fftlog_smooth_padding_negative_slope);
  g_test_add_func ("/ncm/fftlog/smooth_padding/refine_stable", test_ncm_fftlog_smooth_padding_refine_stable);
  g_test_add_func ("/ncm/fftlog/calibrate/zero", test_ncm_fftlog_calibrate_zero);
  g_test_add_func ("/ncm/fftlog/calibrate/max_n/traps", test_ncm_fftlog_calibrate_max_n_traps);
  g_test_add_func ("/ncm/fftlog/calibrate/max_n/subprocess", test_ncm_fftlog_calibrate_max_n_subprocess);
  g_test_add_func ("/ncm/fftlog/calibrate/short_padding/subprocess", test_ncm_fftlog_calibrate_short_padding_subprocess);
  g_test_add_func ("/ncm/fftlog/get_Ym_keeps_eval", test_ncm_fftlog_get_Ym_keeps_eval);
  g_test_add_func ("/ncm/fftlog/smooth_padding/traps", test_ncm_fftlog_smooth_padding_traps);
  g_test_add_func ("/ncm/fftlog/smooth_padding/invalid/subprocess", test_ncm_fftlog_invalid_smooth_padding);
  g_test_add_func ("/ncm/fftlog/smooth_padding/power_law", test_ncm_fftlog_smooth_padding_power_law);
  g_test_add_func ("/ncm/fftlog/odd_full_size", test_ncm_fftlog_odd_full_size);
  g_test_add_func ("/ncm/fftlog/noring_needs_nyquist", test_ncm_fftlog_noring_needs_nyquist);
  g_test_add_func ("/ncm/fftlog/bias/tophatwin2/truth", test_ncm_fftlog_bias_tophatwin2_truth);
  g_test_add_func ("/ncm/fftlog/bias/gausswin2/truth", test_ncm_fftlog_bias_gausswin2_truth);
  g_test_add_func ("/ncm/fftlog/bias/power_law_converges", test_ncm_fftlog_bias_power_law_converges);
  g_test_add_func ("/ncm/fftlog/bias/sbessel_j/truth", test_ncm_fftlog_bias_sbessel_j_truth);
  g_test_add_func ("/ncm/fftlog/sbessel_j/best_lnr0", test_ncm_fftlog_sbessel_j_best_lnr0);
  g_test_add_func ("/ncm/fftlog/bias/best", test_ncm_fftlog_bias_best);
  g_test_add_func ("/ncm/fftlog/bias/traps", test_ncm_fftlog_bias_traps);
  g_test_add_func ("/ncm/fftlog/bias/invalid/range/subprocess", test_ncm_fftlog_bias_invalid_range);
  g_test_add_func ("/ncm/fftlog/bias/invalid/unsupported/subprocess", test_ncm_fftlog_bias_invalid_unsupported);
  g_test_add_func ("/ncm/fftlog/bias/invalid/end_slopes/subprocess", test_ncm_fftlog_bias_invalid_end_slopes);

  g_test_run ();
}

#define NTESTS 20

typedef struct _TestNcmFftlogPlaw
{
  gdouble lnA;
  gdouble ns;
} TestNcmFftlogPlaw;

static gdouble
_test_ncm_fftlog_plaw (gdouble k, gpointer user_data)
{
  TestNcmFftlogPlaw *args = (TestNcmFftlogPlaw *) user_data;

  return exp (args->lnA + log (k) * (args->ns - 1.0));
}

static gdouble
_test_ncm_fftlog_tophatwin2 (gdouble lnk, gpointer user_data)
{
  TestNcmFftlogK *args = (TestNcmFftlogK *) user_data;
  const gdouble kr     = exp (lnk + args->lnr);
  const gdouble k      = exp (lnk);

  return GSL_FN_EVAL (&args->Fk, k) * k * gsl_pow_2 (3.0 * ncm_sf_sbessel (1, kr) / kr);
}

static gdouble
_test_ncm_fftlog_gausswin2 (gdouble lnk, gpointer user_data)
{
  TestNcmFftlogK *args = (TestNcmFftlogK *) user_data;
  const gdouble kr     = exp (lnk + args->lnr);
  const gdouble k      = exp (lnk);

  return GSL_FN_EVAL (&args->Fk, k) * k * exp (-kr * kr);
}

static gdouble
_test_ncm_fftlog_sbessel_j (gdouble lnk, gpointer user_data)
{
  TestNcmFftlogK *args = (TestNcmFftlogK *) user_data;
  const gdouble kr     = exp (lnk + args->lnr);
  const gdouble k      = exp (lnk);

  return GSL_FN_EVAL (&args->Fk, k) * k * ncm_sf_sbessel (args->ell, kr);
}

void
test_ncm_fftlog_tophatwin2_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (10000, 20000);
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (1.0e+3), log (1.0e+6));

  test->fftlog       = fftlog;
  test->Fk.function  = &_test_ncm_fftlog_plaw;
  test->Fk.params    = arg;
  test->KFk.function = &_test_ncm_fftlog_tophatwin2;
  test->KFk.params   = argK;
  test->argK         = argK;

  test->lnk_i = g_test_rand_double_range (log (1.0e-4), log (1.0e0));
  test->lnk_f = test->lnk_i + Lk;

  test->ntests = NTESTS;

  arg->lnA = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns  = g_test_rand_double_range (0.5, 1.5);

  argK->lnr = 0.0;
  argK->Fk  = test->Fk;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_lnr0 (fftlog, -0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));
  g_assert_true (NCM_IS_FFTLOG_TOPHATWIN2 (fftlog));
}

void
test_ncm_fftlog_gausswin2_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (1000, 2000);
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (1.0e+3), log (1.0e+6));

  test->fftlog       = fftlog;
  test->Fk.function  = &_test_ncm_fftlog_plaw;
  test->Fk.params    = arg;
  test->KFk.function = &_test_ncm_fftlog_gausswin2;
  test->KFk.params   = argK;
  test->argK         = argK;

  test->lnk_i = g_test_rand_double_range (log (1.0e-4), log (1.0e0));
  test->lnk_f = test->lnk_i + Lk;

  test->ntests = NTESTS;

  arg->lnA = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns  = g_test_rand_double_range (0.5, 1.5);

  argK->lnr = 0.0;
  argK->Fk  = test->Fk;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_lnr0 (fftlog, -0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));
  g_assert_true (NCM_IS_FFTLOG_GAUSSWIN2 (fftlog));
}

void
test_ncm_fftlog_sbessel_j_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (7800, 8000);
  const guint ell        = g_test_rand_int_range  (0, 5);
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (ell, 0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (9.0e+3), log (1.0e+4));

  test->fftlog       = fftlog;
  test->Fk.function  = &_test_ncm_fftlog_plaw;
  test->Fk.params    = arg;
  test->KFk.function = &_test_ncm_fftlog_sbessel_j;
  test->KFk.params   = argK;
  test->argK         = argK;

  test->lnk_i = g_test_rand_double_range (log (1.0e-6), log (4.0e-6));
  test->lnk_f = test->lnk_i + Lk;

  test->ntests = NTESTS;

  arg->lnA = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns  = g_test_rand_double_range (0.5, 0.6);

  argK->lnr = 0.0;
  argK->Fk  = test->Fk;
  argK->ell = ell;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);

  ncm_fftlog_sbessel_j_set_best_lnr0 (NCM_FFTLOG_SBESSEL_J (fftlog));
  ncm_fftlog_sbessel_j_set_best_lnk0 (NCM_FFTLOG_SBESSEL_J (fftlog));

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));
  g_assert_true (NCM_IS_FFTLOG_SBESSEL_J (fftlog));
}

/* The bias leaves G unchanged, so the reference is the plain j_l integral; b = 1/2 also
 * takes the one-lngamma branch of the coefficients. */
void
test_ncm_fftlog_sbessel_j_bias0_5_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (7800, 8000);
  const guint ell        = g_test_rand_int_range  (0, 5);
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (ell, 0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (9.0e+3), log (1.0e+4));

  test->fftlog       = fftlog;
  test->Fk.function  = &_test_ncm_fftlog_plaw;
  test->Fk.params    = arg;
  test->KFk.function = &_test_ncm_fftlog_sbessel_j;
  test->KFk.params   = argK;
  test->argK         = argK;

  test->lnk_i = g_test_rand_double_range (log (1.0e-6), log (1.0e-4));
  test->lnk_f = test->lnk_i + Lk;

  test->ntests = NTESTS;

  arg->lnA = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns  = g_test_rand_double_range (0.5, 0.6);

  argK->lnr = 0.0;
  argK->Fk  = test->Fk;
  argK->ell = ell;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);
  ncm_fftlog_set_bias (fftlog, 0.5);

  ncm_fftlog_sbessel_j_set_best_lnr0 (NCM_FFTLOG_SBESSEL_J (fftlog));
  ncm_fftlog_sbessel_j_set_best_lnk0 (NCM_FFTLOG_SBESSEL_J (fftlog));

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));
  g_assert_true (NCM_IS_FFTLOG_SBESSEL_J (fftlog));
}

void
test_ncm_fftlog_free (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;

  NCM_TEST_FREE (ncm_fftlog_free, fftlog);

  g_free (test->Fk.params);
  g_free (test->KFk.params);
}

void
test_ncm_fftlog_setget (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  gdouble lnr0      = g_test_rand_double_range (log (1.0e-4), log (1.0e-2));
  gdouble lnk0      = g_test_rand_double_range (log (1.0e-4), log (1.0e-2));
  gdouble Lk        = g_test_rand_double_range (log (1.0e+2), log (1.0e+4));

  ncm_fftlog_set_lnr0 (fftlog, lnr0);
  ncm_fftlog_set_lnk0 (fftlog, lnk0);
  ncm_fftlog_set_length (fftlog, Lk);
  ncm_fftlog_set_padding (fftlog, 0.1);

  ncm_assert_cmpdouble_e (ncm_fftlog_get_lnr0 (fftlog), ==, lnr0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_fftlog_get_lnk0 (fftlog), ==, lnk0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_fftlog_get_length (fftlog), ==, Lk, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_fftlog_get_padding (fftlog), ==, 0.1, 1.0e-15, 0.0);

  ncm_fftlog_set_nderivs (fftlog, 1);
  g_assert_cmpint (ncm_fftlog_get_nderivs (fftlog), ==, 1);
  ncm_fftlog_set_nderivs (fftlog, 2);
  g_assert_cmpint (ncm_fftlog_get_nderivs (fftlog), ==, 2);
  ncm_fftlog_set_nderivs (fftlog, 3);
  g_assert_cmpint (ncm_fftlog_get_nderivs (fftlog), ==, 3);
  ncm_fftlog_set_nderivs (fftlog, 1);
  g_assert_cmpint (ncm_fftlog_get_nderivs (fftlog), ==, 1);

  ncm_fftlog_set_size (fftlog, 1000);
  g_assert_cmpint (ncm_fftlog_get_size (fftlog), >=, 1000);

  ncm_fftlog_set_noring (fftlog, TRUE);
  g_assert_true (ncm_fftlog_get_noring (fftlog));
  ncm_fftlog_set_noring (fftlog, FALSE);
  g_assert_false (ncm_fftlog_get_noring (fftlog));

  ncm_fftlog_set_bias (fftlog, 0.25);
  g_assert_cmpfloat (ncm_fftlog_get_bias (fftlog), ==, 0.25);
  {
    gdouble bias;

    g_object_get (G_OBJECT (fftlog), "bias", &bias, NULL);
    g_assert_cmpfloat (bias, ==, 0.25);
  }
  ncm_fftlog_set_bias (fftlog, 0.0);

  ncm_fftlog_set_eval_r_min (fftlog, 1.0e-3);
  ncm_fftlog_set_eval_r_max (fftlog, 1.0e+3);
  ncm_assert_cmpdouble_e (ncm_fftlog_get_eval_r_min (fftlog), ==, 1.0e-3, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_fftlog_get_eval_r_max (fftlog), ==, 1.0e+3, 1.0e-15, 0.0);
  ncm_fftlog_use_eval_interval (fftlog, FALSE);
  {
    gboolean use_eval_interval;

    g_object_get (G_OBJECT (fftlog), "use-eval-int", &use_eval_interval, NULL);
    g_assert_false (use_eval_interval);
  }
  ncm_fftlog_use_eval_interval (fftlog, TRUE);
  {
    gboolean use_eval_interval;

    g_object_get (G_OBJECT (fftlog), "use-eval-int", &use_eval_interval, NULL);
    g_assert_true (use_eval_interval);
  }

  ncm_fftlog_use_smooth_padding (fftlog, TRUE);
  {
    gboolean use_smooth_padding;

    g_object_get (G_OBJECT (fftlog), "use-smooth-padding", &use_smooth_padding, NULL);
    g_assert_true (use_smooth_padding);
  }
  ncm_fftlog_use_smooth_padding (fftlog, FALSE);
  {
    gboolean use_smooth_padding;

    g_object_get (G_OBJECT (fftlog), "use-smooth-padding", &use_smooth_padding, NULL);
    g_assert_false (use_smooth_padding);
  }
}

void
test_ncm_fftlog_eval (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  gdouble reltol    = 1.0e-1;
  NcmVector *lnr;
  guint i, len;

  ncm_fftlog_eval_by_function (fftlog, test->Fk.function, test->Fk.params);
  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  {
    guint size = 0;

    g_assert_nonnull (ncm_fftlog_get_Ym (fftlog, &size));
    g_assert_cmpuint (size, >, 0);
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint l                  = g_test_rand_int_range ((len / 3), 2 * (len / 3));
    const gdouble lnr_l      = ncm_vector_get (lnr, l);
    const gdouble fftlog_res = ncm_fftlog_eval_output (NCM_FFTLOG (fftlog), 0, lnr_l);
    gdouble res, err;

    test->argK->lnr = lnr_l;
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL1, &res, &err);

    if (fabs (fftlog_res / res - 1.0) > reltol)
      ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL2, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }
}

void
test_ncm_fftlog_eval_vector (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  gdouble reltol    = 1.0e-1;
  guint len         = ncm_fftlog_get_size (fftlog);
  NcmVector *lnk    = ncm_vector_new (len);
  NcmVector *Fk     = ncm_vector_new (len);
  NcmVector *lnr;
  guint i;

  ncm_fftlog_get_lnk_vector (fftlog, lnk);
  {
    for (i = 0; i < len; i++)
    {
      const gdouble k = exp (ncm_vector_get (lnk, i));

      ncm_vector_set (Fk, i, test->Fk.function (k, test->Fk.params));
    }
  }

  ncm_fftlog_eval_by_vector (fftlog, Fk);

  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  {
    guint size = 0;

    g_assert_nonnull (ncm_fftlog_get_Ym (fftlog, &size));
    g_assert_cmpuint (size, >, 0);
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint l                  = g_test_rand_int_range ((len / 3), 2 * (len / 3));
    const gdouble lnr_l      = ncm_vector_get (lnr, l);
    const gdouble fftlog_res = ncm_fftlog_eval_output (NCM_FFTLOG (fftlog), 0, lnr_l);
    gdouble res, err;

    test->argK->lnr = lnr_l;
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL1, &res, &err);

    if (fabs (fftlog_res / res - 1.0) > reltol)
      ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL2, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }
}

void
test_ncm_fftlog_eval_calibrate (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  gdouble reltol    = 1.0e-1;
  NcmVector *lnr;
  guint i, len;

  ncm_fftlog_calibrate_size (fftlog, test->Fk.function, test->Fk.params, 1.0e-1);
  ncm_fftlog_eval_by_gsl_function (fftlog, &test->Fk);
  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  for (i = 0; i < test->ntests; i++)
  {
    guint l                  = g_test_rand_int_range ((len / 3), 2 * (len / 3));
    const gdouble lnr_l      = ncm_vector_get (lnr, l);
    const gdouble fftlog_res = ncm_fftlog_eval_output (NCM_FFTLOG (fftlog), 0, lnr_l);
    gdouble res, err;

    test->argK->lnr = lnr_l;
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL1, &res, &err);

    if (fabs (fftlog_res / res - 1.0) > reltol)
      ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL2, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }
}

void
test_ncm_fftlog_eval_calibrate_fail (TestNcmFftlog *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    NcmFftlog *fftlog = test->fftlog;

    /* 1e-14 is out of reach within 100 knots: the calibration must abort, not return
     * a result short of the requested accuracy. */
    ncm_fftlog_set_max_size (fftlog, 100);

    ncm_fftlog_calibrate_size (fftlog, test->Fk.function, test->Fk.params, 1.0e-14);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exceeds the maximum (100)*");
}

void
test_ncm_fftlog_eval_serialized (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_serialize_dup_obj (ser, G_OBJECT (test->fftlog)));
  gdouble reltol    = 1.0e-1;
  NcmVector *lnr;
  guint i, len;

  ncm_fftlog_eval_by_gsl_function (fftlog, &test->Fk);
  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  for (i = 0; i < test->ntests; i++)
  {
    guint l                  = g_test_rand_int_range ((len / 3), 2 * (len / 3));
    const gdouble lnr_l      = ncm_vector_get (lnr, l);
    const gdouble fftlog_res = ncm_fftlog_eval_output (NCM_FFTLOG (fftlog), 0, lnr_l);
    gdouble res, err;

    test->argK->lnr = lnr_l;
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL1, &res, &err);

    if (fabs (fftlog_res / res - 1.0) > reltol)
      ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL2, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }

  ncm_serialize_free (ser);
  ncm_fftlog_free (fftlog);
}

void
test_ncm_fftlog_eval_deriv (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  gdouble reltol    = 1.0e-1;
  NcmVector *lnr;
  guint i, len;

  ncm_fftlog_set_nderivs (fftlog, 1);
  ncm_fftlog_eval_by_function (fftlog, test->Fk.function, test->Fk.params);
  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  {
    NcmVector *Gr;

    Gr = ncm_fftlog_get_vector_Gr (fftlog, 0);
    g_assert_nonnull (Gr);
    ncm_vector_free (Gr);

    Gr = ncm_fftlog_get_vector_Gr (fftlog, 1);
    g_assert_nonnull (Gr);
    ncm_vector_free (Gr);
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint l                  = g_test_rand_int_range ((len / 3), 2 * (len / 3));
    const gdouble lnr_l      = ncm_vector_get (lnr, l);
    const gdouble fftlog_res = ncm_fftlog_eval_output (NCM_FFTLOG (fftlog), 0, lnr_l);
    gdouble res, err;

    test->argK->lnr = lnr_l;
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL1, &res, &err);

    if (fabs (fftlog_res / res - 1.0) > reltol)
      ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL2, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }
}

void
test_ncm_fftlog_eval_use_eval_int (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  gdouble reltol    = 1.0e-1;
  NcmVector *lnr;
  guint i, len;

  ncm_fftlog_eval_by_function (fftlog, test->Fk.function, test->Fk.params);
  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  g_assert_cmpuint (len, >=, 100);
  {
    const gdouble lnr_l = ncm_vector_get (lnr, len / 10);
    const gdouble lnr_u = ncm_vector_get (lnr, 9 * len / 10);

    ncm_fftlog_set_eval_r_min (fftlog, exp (lnr_l));
    ncm_fftlog_set_eval_r_max (fftlog, exp (lnr_u));

    ncm_fftlog_use_eval_interval (fftlog, TRUE);
    ncm_fftlog_eval_by_gsl_function (fftlog, &test->Fk);
    ncm_fftlog_prepare_splines (fftlog);
    ncm_vector_free (lnr);

    lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
    len = ncm_vector_len (lnr);
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint l                  = g_test_rand_int_range ((len / 3), 2 * (len / 3));
    const gdouble lnr_l      = ncm_vector_get (lnr, l);
    const gdouble fftlog_res = ncm_fftlog_eval_output (NCM_FFTLOG (fftlog), 0, lnr_l);
    gdouble res, err;

    test->argK->lnr = lnr_l;
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL1, &res, &err);

    if (fabs (fftlog_res / res - 1.0) > reltol)
      ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, NUMINT_RELTOL2, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }
}

void
test_ncm_fftlog_eval_smooth_padding (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  NcmVector *lnr;
  guint len;

  ncm_fftlog_use_smooth_padding (fftlog, TRUE);
  ncm_fftlog_eval_by_function (fftlog, test->Fk.function, test->Fk.params);
  ncm_fftlog_prepare_splines (fftlog);
  lnr = ncm_fftlog_get_vector_lnr (NCM_FFTLOG (fftlog));
  len = ncm_vector_len (lnr);

  g_assert_cmpint (len, >, 0);
}

void
test_ncm_fftlog_tophatwin2_traps (TestNcmFftlog *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/fftlog/tophatwin2/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  /* A non-positive period makes the knot spacing L / N meaningless */
  g_test_trap_subprocess ("/ncm/fftlog/tophatwin2/invalid/length/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*period must be positive*");
}

void
test_ncm_fftlog_gausswin2_traps (TestNcmFftlog *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/fftlog/gausswin2/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fftlog_sbessel_j_traps (TestNcmFftlog *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/fftlog/sbessel_j/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fftlog_invalid_st (TestNcmFftlog *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

void
test_ncm_fftlog_invalid_length (TestNcmFftlog *test, gconstpointer pdata)
{
  ncm_fftlog_set_length (test->fftlog, 0.0);
}

static gdouble
_test_gauss_k3 (const gdouble k, gpointer user_data)
{
  return gsl_pow_3 (k) * exp (-k * k);
}

/*
 * G(r) = int F(k) W(kr)^2 dk for F = k^3 exp(-k^2), which vanishes at both ends of the
 * grid, at exact knots (no-ringing off, so ln r = n L / N'). Truth from mpmath at 30
 * digits; error measured at most 5.8e-15 of the peak.
 */
void
test_ncm_fftlog_tophatwin2_truth (void)
{
  const gint pos[]      = {-40, -10, 0, 10, 30};
  const gdouble truth[] = {
    0.49966783048136725739, 0.46163609991915580706, 0.34271556221491577223,
    0.11109960959342109508, 0.00071932026187923236137
  };
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));
  NcmVector *Gr;
  gdouble peak;
  guint i, N_2;

  ncm_fftlog_set_noring (fftlog, FALSE);
  ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);

  g_assert_cmpuint (ncm_fftlog_get_size (fftlog), ==, 250);
  N_2  = ncm_fftlog_get_size (fftlog) / 2;
  Gr   = ncm_fftlog_get_vector_Gr (fftlog, 0);
  peak = ncm_vector_get_max (Gr);

  for (i = 0; i < G_N_ELEMENTS (pos); i++)
    g_assert_cmpfloat (fabs (ncm_vector_get (Gr, N_2 + pos[i]) - truth[i]), <, 5.0e-14 * peak);

  ncm_vector_free (Gr);
  ncm_fftlog_free (fftlog);
}

/* ncm_fftlog_get_Ym() writes into the transform's own buffer; the next evaluation must
 * not use it as if it were still prepared. */
void
test_ncm_fftlog_get_Ym_keeps_eval (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (1.0, -0.5, 20.0, 250));
  NcmVector *before, *after;
  guint size, i;

  ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);
  before = ncm_vector_dup (ncm_fftlog_peek_output_vector (fftlog, 0));

  ncm_fftlog_get_Ym (fftlog, &size);
  g_assert_cmpuint (size, ==, 2 * ncm_fftlog_get_full_size (fftlog));

  ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);
  after = ncm_fftlog_peek_output_vector (fftlog, 0);

  for (i = 0; i < ncm_vector_len (before); i++)
    g_assert_cmpfloat (ncm_vector_get (after, i), ==, ncm_vector_get (before, i));

  ncm_vector_free (before);
  ncm_fftlog_free (fftlog);
}

void
test_ncm_fftlog_smooth_padding_traps (void)
{
  g_test_trap_subprocess ("/ncm/fftlog/smooth_padding/invalid/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*smooth padding needs F > 0*");
}

static gdouble
_test_negative_edge (const gdouble k, gpointer user_data)
{
  return -gsl_pow_3 (k) * exp (-k * k);
}

/* The power-law continuation is fitted to log F: F <= 0 at an end used to give NaN. */
void
test_ncm_fftlog_invalid_smooth_padding (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));

  ncm_fftlog_use_smooth_padding (fftlog, TRUE);
  ncm_fftlog_eval_by_function (fftlog, &_test_negative_edge, NULL);
}

static gdouble
_test_sqrt_k (const gdouble k, gpointer user_data)
{
  return sqrt (k);
}

/* For F = k^s the smooth padding continues the table with the exact power law at the
 * end where F decays, so G(r) follows the infinite-range law G ~ r^-(s+1) towards the
 * inverse of that end, where zeros in the padding are off by the truncation of the
 * table (3e-4 for the tophat at k_min r = 1e-3). Checked at k_min r from 1e-5 to 1e-3
 * through the log-derivative, which needs no normalisation; measured 6.2e-10 (tophat) and
 * 4.4e-10 (Gaussian) with the tapered padding. */
void
test_ncm_fftlog_smooth_padding_power_law (void)
{
  const gdouble lnk_min = log (1.0e-8);
  const gdouble lnk_max = log (1.0e3);
  const gdouble lnk0    = 0.5 * (lnk_min + lnk_max);
  const gdouble s       = 0.5;
  NcmFftlog *fftlogs[2] = {
    NCM_FFTLOG (ncm_fftlog_tophatwin2_new (-lnk0, lnk0, lnk_max - lnk_min, 1000)),
    NCM_FFTLOG (ncm_fftlog_gausswin2_new (-lnk0, lnk0, lnk_max - lnk_min, 1000)),
  };
  guint w;

  for (w = 0; w < 2; w++)
  {
    NcmFftlog *fftlog = fftlogs[w];
    gint i;

    ncm_fftlog_set_padding (fftlog, 1.0);
    ncm_fftlog_set_nderivs (fftlog, 1);
    ncm_fftlog_use_smooth_padding (fftlog, TRUE);
    ncm_fftlog_eval_by_function (fftlog, &_test_sqrt_k, NULL);
    ncm_fftlog_prepare_splines (fftlog);

    for (i = 0; i <= 20; i++)
    {
      const gdouble lnr    = log (1.0e3) + i * log (1.0e2) / 20.0;
      const gdouble G0     = ncm_fftlog_eval_output (fftlog, 0, lnr);
      const gdouble G1     = ncm_fftlog_eval_output (fftlog, 1, lnr);
      const gdouble dlnGdr = G1 / G0;

      ncm_assert_cmpdouble_e (dlnGdr, ==, -(s + 1.0), 1.0e-8, 0.0);
    }

    ncm_fftlog_free (fftlog);
  }
}

static gdouble
_test_k3_gauss (const gdouble k, gpointer user_data)
{
  return gsl_pow_3 (k) * exp (-k * k);
}

/* The knots sit at (i - Nf_2) L / N on both sides of the transform, which cancels in the
 * phase only for an even full size; an odd one used to shift the output by one knot
 * (5e-3 at N = 1000). Compare an odd full size (padding 1.187 gives 2187) with an even
 * one (1.25 gives 2250) for an input that vanishes at both ends; both periods are over
 * twice the interval, so the periodic images of the input stay below 1e-11. */
void
test_ncm_fftlog_odd_full_size (void)
{
  const gdouble lnk_min = log (1.0e-3);
  const gdouble lnk_max = log (1.0e2);
  const gdouble lnk0    = 0.5 * (lnk_min + lnk_max);
  NcmFftlog *odd        = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (-lnk0, lnk0, lnk_max - lnk_min, 1000));
  NcmFftlog *even       = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (-lnk0, lnk0, lnk_max - lnk_min, 1000));
  gint i;

  ncm_fftlog_set_padding (odd, 1.187);
  ncm_fftlog_set_padding (even, 1.25);
  g_assert_cmpint (ncm_fftlog_get_full_size (odd) % 2, ==, 1);
  g_assert_cmpint (ncm_fftlog_get_full_size (even) % 2, ==, 0);

  ncm_fftlog_eval_by_function (odd, &_test_k3_gauss, NULL);
  ncm_fftlog_eval_by_function (even, &_test_k3_gauss, NULL);
  ncm_fftlog_prepare_splines (odd);
  ncm_fftlog_prepare_splines (even);

  for (i = 0; i <= 20; i++)
  {
    const gdouble lnr    = log (0.1) + i * log (1.0e2) / 20.0;
    const gdouble G_odd  = ncm_fftlog_eval_output (odd, 0, lnr);
    const gdouble G_even = ncm_fftlog_eval_output (even, 0, lnr);

    ncm_assert_cmpdouble_e (G_odd, ==, G_even, 1.0e-6, 0.0);
  }

  ncm_fftlog_free (odd);
  ncm_fftlog_free (even);
}

static gdouble
_test_power (const gdouble k, gpointer user_data)
{
  return pow (k, *(gdouble *) user_data);
}

/*
 * Gaussian window with F = k^(-1/2) on [1e-8, 1e3]: G(r) = Gamma(1/4) / (2 sqrt(r)). F grows
 * below k_min while F k decays, so the low-k continuation must keep that tail: cutting it on
 * the slope of F, not of F k, loses 1.2e-3 at r = 1e3 and 1.2e-2 at 1e5. Measured with the
 * slope of F k: 8.5e-5 and 8.5e-4, the tail beyond the taper at 0.4-0.8 of the padding (the
 * former blend of the two ends reached further, 6.8e-6 and 6.8e-5, but moved with N).
 */
void
test_ncm_fftlog_smooth_padding_negative_slope (void)
{
  gdouble s             = -0.5;
  NcmFftlog *fftlog     = NCM_FFTLOG (ncm_fftlog_gausswin2_new (-0.5 * log (1.0e-5), 0.5 * log (1.0e-5), log (1.0e11), 2000));
  const gdouble r_a[]   = {1.0e3, 1.0e5};
  const gdouble tol_a[] = {2.0e-4, 2.0e-3};
  guint i;

  ncm_fftlog_set_noring (fftlog, FALSE);
  ncm_fftlog_use_smooth_padding (fftlog, TRUE);
  ncm_fftlog_eval_by_function (fftlog, &_test_power, &s);
  ncm_fftlog_prepare_splines (fftlog);

  for (i = 0; i < G_N_ELEMENTS (r_a); i++)
  {
    const gdouble truth = 0.5 * tgamma (0.25) / sqrt (r_a[i]);

    g_assert_cmpfloat (fabs (ncm_fftlog_eval_output (fftlog, 0, log (r_a[i])) / truth - 1.0), <, tol_a[i]);
  }

  ncm_fftlog_free (fftlog);
}

/*
 * A constant input continued into the padding: the result must not move as the grid is
 * refined. With the continuation placed by slot index it moved by 2.3e-7, 1.2e-7, 5.8e-8 per
 * doubling of N (first order); placed in ln k it is stable to 1e-17. A fractional padding
 * rounds to a period L_T that changes with N (12.9973 at 3000 knots, 12.9983 at 4524 for
 * p = 0.3, L = 10): the former blend of the two ends, which mixes values whose k^b differ by
 * e^(b L_T), moved by 5.7e-2 of the peak between them with the bias of
 * ncm_fftlog_get_best_bias(); the tapered ends move by 2.5e-8, the roundoff of that bias at
 * the smallest r.
 */
void
test_ncm_fftlog_smooth_padding_refine_stable (void)
{
  gdouble s    = 0.0;
  gdouble prev = 0.0;
  guint n;

  for (n = 1000; n <= 4000; n *= 2)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 10.0, n));
    gdouble y;

    ncm_fftlog_set_noring (fftlog, FALSE);
    ncm_fftlog_use_smooth_padding (fftlog, TRUE);
    ncm_fftlog_eval_by_function (fftlog, &_test_power, &s);
    ncm_fftlog_prepare_splines (fftlog);

    y = ncm_fftlog_eval_output (fftlog, 0, 4.5);

    if (n > 1000)
      g_assert_cmpfloat (fabs (y - prev), <, 1.0e-13 * fabs (y));

    prev = y;
    ncm_fftlog_free (fftlog);
  }

  {
    const guint N_a[] = {3000, 4524};
    NcmVector *Gr[2];
    NcmFftlog *fftlog[2];
    NcmVector *lnr;
    gdouble bias = 0.0, peak;
    guint j, i;

    s = 0.0;

    for (j = 0; j < 2; j++)
    {
      fftlog[j] = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 10.0, N_a[j]));
      ncm_fftlog_set_padding (fftlog[j], 0.3);
      ncm_fftlog_set_noring (fftlog[j], FALSE);
      ncm_fftlog_use_smooth_padding (fftlog[j], TRUE);

      if (j == 0)
      {
        ncm_fftlog_eval_by_function (fftlog[j], &_test_power, &s);
        bias = ncm_fftlog_get_best_bias (fftlog[j]);
      }

      ncm_fftlog_set_bias (fftlog[j], bias);
      ncm_fftlog_eval_by_function (fftlog[j], &_test_power, &s);
      ncm_fftlog_prepare_splines (fftlog[j]);
      Gr[j] = ncm_fftlog_get_vector_Gr (fftlog[j], 0);
    }

    g_assert_cmpfloat (ncm_fftlog_get_full_length (fftlog[0]), !=, ncm_fftlog_get_full_length (fftlog[1]));
    peak = ncm_vector_get_max (Gr[1]);
    lnr  = ncm_fftlog_get_vector_lnr (fftlog[0]);

    for (i = 0; i < ncm_vector_len (Gr[0]); i++)
    {
      const gdouble G = ncm_fftlog_eval_output (fftlog[1], 0, ncm_vector_get (lnr, i));

      g_assert_cmpfloat (fabs (ncm_vector_get (Gr[0], i) - G), <, 1.0e-7 * peak);
    }

    ncm_vector_free (lnr);

    for (j = 0; j < 2; j++)
    {
      ncm_vector_free (Gr[j]);
      ncm_fftlog_free (fftlog[j]);
    }
  }
}

static gdouble
_test_zero (const gdouble k, gpointer user_data)
{
  return 0.0;
}

/* An identically zero transform is converged at once: 0/0 used to read as not converged. */
void
test_ncm_fftlog_calibrate_zero (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 100));
  NcmVector *Gr;
  guint i;

  ncm_fftlog_set_max_size (fftlog, 1000);
  ncm_fftlog_calibrate_size (fftlog, &_test_zero, NULL, 1.0e-10);

  Gr = ncm_fftlog_peek_output_vector (fftlog, 0);

  for (i = 0; i < ncm_vector_len (Gr); i++)
    g_assert_cmpfloat (ncm_vector_get (Gr, i), ==, 0.0);

  g_assert_cmpuint (ncm_fftlog_get_size (fftlog), <=, 1000);

  ncm_fftlog_free (fftlog);
}

void
test_ncm_fftlog_calibrate_max_n_traps (void)
{
  g_test_trap_subprocess ("/ncm/fftlog/calibrate/max_n/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exceeds the maximum (100)*");

  /* F = 1 with a padding of 0.3: the periodic image, e^(-L_T) = 2e-6 of the peak without a
   * bias, moves with the rounded period, and the calibration once stopped at 3770 knots on
   * a chance agreement. It must reach neither 1e-10 nor a false stop by 20000 knots. */
  g_test_trap_subprocess ("/ncm/fftlog/calibrate/short_padding/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exceeds the maximum (20000)*");
}

void
test_ncm_fftlog_calibrate_short_padding_subprocess (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 10.0, 100));
  gdouble s         = 0.0;

  ncm_fftlog_set_padding (fftlog, 0.3);
  ncm_fftlog_set_noring (fftlog, FALSE);
  ncm_fftlog_use_smooth_padding (fftlog, TRUE);
  ncm_fftlog_set_max_size (fftlog, 20000);
  ncm_fftlog_calibrate_size (fftlog, &_test_power, &s, 1.0e-10);
}

/* Starting at the maximum, the first growth step already passes it: the calibration must
 * stop before transforming at 120 knots, not return there. */
void
test_ncm_fftlog_calibrate_max_n_subprocess (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 100));

  ncm_fftlog_set_max_size (fftlog, 100);
  ncm_fftlog_calibrate_size (fftlog, &_test_gauss_k3, NULL, 1.0e-1);
}

/*
 * The low-ringing adjustment sets the phase of the Nyquist mode, so an odd full size, which
 * has none, must leave the output grid where it is; an even one moves it by less than one
 * knot. Either way the requested ln r0 stays as set: shifting it in place made each new
 * size shift from the last one.
 */
void
test_ncm_fftlog_noring_needs_nyquist (void)
{
  const gdouble pad_a[] = {1.187, 1.0}; /* full sizes 2187 (odd) and 2000 (even) */
  guint j;

  for (j = 0; j < G_N_ELEMENTS (pad_a); j++)
  {
    NcmFftlog *on      = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.3, -0.3, 20.0, 1000));
    NcmFftlog *off     = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.3, -0.3, 20.0, 1000));
    const gboolean odd = (j == 0);
    NcmVector *lnr_on, *lnr_off;
    gdouble shift;

    ncm_fftlog_set_padding (on, pad_a[j]);
    ncm_fftlog_set_padding (off, pad_a[j]);
    ncm_fftlog_set_noring (on, TRUE);
    ncm_fftlog_set_noring (off, FALSE);
    g_assert_cmpint (ncm_fftlog_get_full_size (on) % 2, ==, odd ? 1 : 0);

    ncm_fftlog_eval_by_function (on, &_test_gauss_k3, NULL);
    ncm_fftlog_eval_by_function (off, &_test_gauss_k3, NULL);

    lnr_on  = ncm_fftlog_get_vector_lnr (on);
    lnr_off = ncm_fftlog_get_vector_lnr (off);
    shift   = ncm_vector_get (lnr_on, 0) - ncm_vector_get (lnr_off, 0);

    if (odd)
      g_assert_cmpfloat (shift, ==, 0.0);
    else
      g_assert_cmpfloat (fabs (shift), <, 20.0 / ncm_fftlog_get_size (on));

    g_assert_cmpfloat (ncm_fftlog_get_lnr0 (on), ==, 0.3);

    ncm_vector_free (lnr_on);
    ncm_vector_free (lnr_off);

    ncm_fftlog_free (on);
    ncm_fftlog_free (off);
  }
}

/*
 * The bias changes what the transform represents, not G(r): the tophat truth of
 * test_ncm_fftlog_tophatwin2_truth() with b = 0.5 and 1. Measured 1.5e-14 and 1.6e-12 of
 * the peak; a positive bias multiplies the roundoff at small r by r^-(1 + b).
 */
void
test_ncm_fftlog_bias_tophatwin2_truth (void)
{
  const gint pos[]      = {-40, -10, 0, 10, 30};
  const gdouble truth[] = {
    0.49966783048136725739, 0.46163609991915580706, 0.34271556221491577223,
    0.11109960959342109508, 0.00071932026187923236137
  };
  const gdouble bias_a[] = {0.5, 1.0};
  const gdouble tol_a[]  = {5.0e-14, 5.0e-12};
  guint j;

  for (j = 0; j < G_N_ELEMENTS (bias_a); j++)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));
    NcmVector *Gr;
    gdouble peak;
    guint i, N_2;

    ncm_fftlog_set_noring (fftlog, FALSE);
    ncm_fftlog_set_bias (fftlog, bias_a[j]);
    ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);

    N_2  = ncm_fftlog_get_size (fftlog) / 2;
    Gr   = ncm_fftlog_get_vector_Gr (fftlog, 0);
    peak = ncm_vector_get_max (Gr);

    for (i = 0; i < G_N_ELEMENTS (pos); i++)
      g_assert_cmpfloat (fabs (ncm_vector_get (Gr, N_2 + pos[i]) - truth[i]), <, tol_a[j] * peak);

    ncm_vector_free (Gr);
    ncm_fftlog_free (fftlog);
  }
}

/*
 * Gaussian window with F = k^3 exp(-k^2): G(r) = 1 / (2 (1 + r^2)^2), and its first two
 * ln r derivatives, which carry the bias through the factor -(1 + b + a). Checked over the
 * central knots with b = 0.5 and 1; measured at most 4.9e-14 of each component's peak.
 */
static gdouble
_test_gauss_k3_gausswin2 (guint nd, const gdouble lnr)
{
  const gdouble r2   = exp (2.0 * lnr);
  const gdouble opr2 = 1.0 + r2;

  switch (nd)
  {
    case 0:
      return 0.5 / gsl_pow_2 (opr2);

    case 1:
      return -2.0 * r2 / gsl_pow_3 (opr2);

    default:
      return -4.0 * r2 * (1.0 - 2.0 * r2) / gsl_pow_4 (opr2);
  }
}

void
test_ncm_fftlog_bias_gausswin2_truth (void)
{
  const gdouble bias_a[] = {0.5, 1.0};
  guint j;

  for (j = 0; j < G_N_ELEMENTS (bias_a); j++)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 20.0, 250));
    NcmVector *lnr;
    gint i, N_2;
    guint nd;

    ncm_fftlog_set_noring (fftlog, FALSE);
    ncm_fftlog_set_nderivs (fftlog, 2);
    ncm_fftlog_set_bias (fftlog, bias_a[j]);
    ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);

    lnr = ncm_fftlog_get_vector_lnr (fftlog);
    N_2 = ncm_fftlog_get_size (fftlog) / 2;

    for (nd = 0; nd <= 2; nd++)
    {
      NcmVector *Gr = ncm_fftlog_get_vector_Gr (fftlog, nd);
      gdouble peak  = 0.0;

      for (i = 0; i < (gint) ncm_vector_len (lnr); i++)
        peak = GSL_MAX (peak, fabs (_test_gauss_k3_gausswin2 (nd, ncm_vector_get (lnr, i))));

      for (i = N_2 - 40; i <= N_2 + 30; i++)
      {
        const gdouble truth = _test_gauss_k3_gausswin2 (nd, ncm_vector_get (lnr, i));

        /* Measured: 5.5e-14 peak on x86-64 Linux, 2.4e-13 peak on macOS arm64 */
        g_assert_cmpfloat (fabs (ncm_vector_get (Gr, i) - truth), <, 5.0e-13 * peak);
      }

      ncm_vector_free (Gr);
    }

    ncm_vector_free (lnr);
    ncm_fftlog_free (fftlog);
  }
}

static gdouble
_test_gauss_kl2 (const gdouble k, gpointer user_data)
{
  const guint ell = GPOINTER_TO_UINT (user_data);

  return gsl_pow_uint (k, ell + 2) * exp (-0.5 * k * k);
}

/*
 * The Gaussian is its own Hankel transform: F = k^(l+2) exp(-k^2/2) gives
 * G = sqrt(pi/2) r^l exp(-r^2/2) and dG/dln r = G (l - r^2). Biases inside the margin of
 * ncm_fftlog_get_best_bias() (0.9 here), down to ell + 1 below zero; measured at most
 * 2.3e-12 (G) and 2.3e-10 (dG) of the peak over the central knots.
 */
void
test_ncm_fftlog_bias_sbessel_j_truth (void)
{
  const guint ell_a[]    = {0, 2, 5};
  const gdouble bias_a[] = {0.0, -1.5, -3.0};
  guint j;

  for (j = 0; j < G_N_ELEMENTS (ell_a); j++)
  {
    const guint ell   = ell_a[j];
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (ell, 0.0, 0.0, 20.0, 250));
    NcmVector *lnr, *G, *dG;
    gdouble peak = 0.0, dpeak = 0.0;
    gint i, N_2;

    ncm_fftlog_set_noring (fftlog, FALSE);
    ncm_fftlog_set_nderivs (fftlog, 1);
    ncm_fftlog_set_bias (fftlog, bias_a[j]);
    ncm_fftlog_eval_by_function (fftlog, &_test_gauss_kl2, GUINT_TO_POINTER (ell));

    lnr = ncm_fftlog_get_vector_lnr (fftlog);
    G   = ncm_fftlog_get_vector_Gr (fftlog, 0);
    dG  = ncm_fftlog_get_vector_Gr (fftlog, 1);
    N_2 = ncm_fftlog_get_size (fftlog) / 2;

    for (i = 0; i < (gint) ncm_vector_len (lnr); i++)
    {
      const gdouble r  = exp (ncm_vector_get (lnr, i));
      const gdouble Gt = sqrt (M_PI_2) * gsl_pow_uint (r, ell) * exp (-0.5 * r * r);

      peak  = GSL_MAX (peak, fabs (Gt));
      dpeak = GSL_MAX (dpeak, fabs (Gt * (ell - r * r)));
    }

    for (i = N_2 - 40; i <= N_2 + 30; i++)
    {
      const gdouble r  = exp (ncm_vector_get (lnr, i));
      const gdouble Gt = sqrt (M_PI_2) * gsl_pow_uint (r, ell) * exp (-0.5 * r * r);

      g_assert_cmpfloat (fabs (ncm_vector_get (G, i) - Gt), <, 1.0e-11 * peak);
      g_assert_cmpfloat (fabs (ncm_vector_get (dG, i) - Gt * (ell - r * r)), <, 1.0e-9 * dpeak);
    }

    ncm_vector_free (lnr);
    ncm_vector_free (G);
    ncm_vector_free (dG);
    ncm_fftlog_free (fftlog);
  }
}

static NcmFftlog *
_test_bias_power_law_fftlog (gboolean tophat, gdouble *s)
{
  const gdouble lnk_min = log (1.0e-10);
  const gdouble lnk_max = log (1.0e3);
  const gdouble lnk0    = 0.5 * (lnk_min + lnk_max);
  NcmFftlog *fftlog     = tophat ?
                          NCM_FFTLOG (ncm_fftlog_tophatwin2_new (-lnk0, lnk0, lnk_max - lnk_min, 100)) :
                          NCM_FFTLOG (ncm_fftlog_gausswin2_new (-lnk0, lnk0, lnk_max - lnk_min, 100));

  ncm_fftlog_set_padding (fftlog, 1.0);
  ncm_fftlog_set_nderivs (fftlog, 2);
  ncm_fftlog_use_smooth_padding (fftlog, TRUE);
  ncm_fftlog_eval_by_function (fftlog, &_test_power, s);

  return fftlog;
}

/*
 * The case of test_bias_castro.py: F = k^-0.6 over 13 decades. Unbiased, F grows about
 * e^25 across the table and its padding, and the change from N to 1.2N stalls at 6e-5 of
 * the peak (roundoff). With the bias of ncm_fftlog_get_best_bias() the padded function is
 * nearly flat; measured calibrations to 1e-9 at 2160 knots and to 1e-11 at 8000.
 */
void
test_ncm_fftlog_bias_power_law_converges (void)
{
  gdouble s         = -0.6;
  NcmFftlog *fftlog = _test_bias_power_law_fftlog (TRUE, &s);

  ncm_fftlog_set_bias (fftlog, ncm_fftlog_get_best_bias (fftlog));
  ncm_fftlog_set_max_size (fftlog, 10000);
  ncm_fftlog_calibrate_size (fftlog, &_test_power, &s, 1.0e-9);

  g_assert_cmpuint (ncm_fftlog_get_size (fftlog), <=, 10000);

  ncm_fftlog_free (fftlog);
}

static gdouble
_test_k3_over_1pk2_2 (const gdouble k, gpointer user_data)
{
  return gsl_pow_3 (k) / gsl_pow_2 (1.0 + k * k);
}

static gdouble
_test_k_plus_k_m02 (const gdouble k, gpointer user_data)
{
  return k + pow (k, -0.2);
}

/*
 * The rules of ncm_fftlog_get_best_bias(): the bias closest to zero between the end
 * slopes, their midpoint when F grows at both ends, kept ln(1/eps)/L_T inside the kernel's
 * range, and zero for a kernel without one.
 */
void
test_ncm_fftlog_bias_best (void)
{
  gdouble s_min, s_max, s;

  /* A power law: b = s, inside the range for both windows. */
  {
    NcmFftlog *fftlog;

    s      = 0.6;
    fftlog = _test_bias_power_law_fftlog (TRUE, &s);
    ncm_fftlog_get_end_slopes (fftlog, &s_min, &s_max);
    ncm_assert_cmpdouble_e (s_min, ==, s, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (s_max, ==, s, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (ncm_fftlog_get_best_bias (fftlog), ==, s, 1.0e-12, 0.0);
    ncm_fftlog_free (fftlog);
  }

  /* Near either end of the tophat range the bias stops at the margin; the Gaussian
   * range has no upper end. */
  {
    NcmFftlog *fftlog;

    s      = -0.6;
    fftlog = _test_bias_power_law_fftlog (TRUE, &s);
    g_assert_cmpfloat (ncm_fftlog_get_best_bias (fftlog), ==, -1.0 - log (GSL_DBL_EPSILON) / ncm_fftlog_get_full_length (fftlog));
    ncm_fftlog_free (fftlog);

    s      = 2.9;
    fftlog = _test_bias_power_law_fftlog (TRUE, &s);
    g_assert_cmpfloat (ncm_fftlog_get_best_bias (fftlog), ==, 3.0 + log (GSL_DBL_EPSILON) / ncm_fftlog_get_full_length (fftlog));
    ncm_fftlog_free (fftlog);

    fftlog = _test_bias_power_law_fftlog (FALSE, &s);
    ncm_assert_cmpdouble_e (ncm_fftlog_get_best_bias (fftlog), ==, s, 1.0e-12, 0.0);
    ncm_fftlog_free (fftlog);
  }

  /* Decaying at both ends (slopes 3 and -1): zero, as without the bias. */
  {
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));

    ncm_fftlog_use_smooth_padding (fftlog, TRUE);
    ncm_fftlog_eval_by_function (fftlog, &_test_k3_over_1pk2_2, NULL);
    ncm_fftlog_get_end_slopes (fftlog, &s_min, &s_max);
    ncm_assert_cmpdouble_e (s_min, ==, 3.0, 1.0e-6, 0.0);
    ncm_assert_cmpdouble_e (s_max, ==, -1.0, 1.0e-6, 0.0);
    g_assert_cmpfloat (ncm_fftlog_get_best_bias (fftlog), ==, 0.0);
    ncm_fftlog_free (fftlog);
  }

  /* Growing at both ends (slopes -0.2 and 1): the midpoint. */
  {
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));

    ncm_fftlog_use_smooth_padding (fftlog, TRUE);
    ncm_fftlog_eval_by_function (fftlog, &_test_k_plus_k_m02, NULL);
    ncm_fftlog_get_end_slopes (fftlog, &s_min, &s_max);
    ncm_assert_cmpdouble_e (s_min, ==, -0.2, 1.0e-4, 0.0);
    ncm_assert_cmpdouble_e (s_max, ==, 1.0, 1.0e-4, 0.0);
    g_assert_cmpfloat (ncm_fftlog_get_best_bias (fftlog), ==, 0.5 * (s_min + s_max));
    ncm_fftlog_free (fftlog);
  }

  /* The spherical Bessel range depends on ell. */
  {
    NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (3, 0.0, 0.0, 20.0, 250));
    gdouble bias_min, bias_max;

    ncm_fftlog_get_bias_range (fftlog, &bias_min, &bias_max);
    g_assert_cmpfloat (bias_min, ==, -4.0);
    g_assert_cmpfloat (bias_max, ==, 1.0);
    ncm_fftlog_free (fftlog);
  }

  /* No range: zero. */
  {
    NcmFftlog *fftlog = test_ncm_fftlog_no_bias_new ();
    gdouble bias_min, bias_max;

    ncm_fftlog_get_bias_range (fftlog, &bias_min, &bias_max);
    g_assert_cmpfloat (bias_min, ==, 0.0);
    g_assert_cmpfloat (bias_max, ==, 0.0);

    s = 0.6;
    ncm_fftlog_use_smooth_padding (fftlog, TRUE);
    ncm_fftlog_eval_by_function (fftlog, &_test_power, &s);
    g_assert_cmpfloat (ncm_fftlog_get_best_bias (fftlog), ==, 0.0);
    ncm_fftlog_free (fftlog);
  }
}

void
test_ncm_fftlog_bias_traps (void)
{
  g_test_trap_subprocess ("/ncm/fftlog/bias/invalid/range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*coefficients only for a bias in (-1, 3), got 3.5*");

  g_test_trap_subprocess ("/ncm/fftlog/bias/invalid/unsupported/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*coefficients only for a bias in (0, 0), got 0.5*");

  g_test_trap_subprocess ("/ncm/fftlog/bias/invalid/end_slopes/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*last evaluation with the smooth padding*");
}

/* The tophat coefficients diverge at b = 3, where W^2 t^b stops being integrable. */
void
test_ncm_fftlog_bias_invalid_range (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));

  ncm_fftlog_set_bias (fftlog, 3.5);
  ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);
}

/* A kernel whose coefficients ignore the bias must refuse one, not return a wrong G. */
void
test_ncm_fftlog_bias_invalid_unsupported (void)
{
  NcmFftlog *fftlog = test_ncm_fftlog_no_bias_new ();

  ncm_fftlog_set_bias (fftlog, 0.5);
  ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);
}

void
test_ncm_fftlog_bias_invalid_end_slopes (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (0.0, 0.0, 20.0, 250));
  gdouble s_min, s_max;

  ncm_fftlog_eval_by_function (fftlog, &_test_gauss_k3, NULL);
  ncm_fftlog_get_end_slopes (fftlog, &s_min, &s_max);
}

static gdouble
_test_power_m05 (const gdouble k, gpointer user_data)
{
  return 1.0 / sqrt (k);
}

/*
 * set_best_lnr0() puts k0 r0 at the first maximum x* of j_l (mpmath, 20 digits), and so the
 * output grid on the r the input determines: for F = k^(-1/2) on [1e-4, 1e4] and l = 10,
 * 91.9% of the grid is within 1e-6 of r^(-1/2) sqrt(pi) 2^(-3/2) Gamma(21/4) / Gamma(25/4)
 * (the earlier rule, balancing the kernel's size at the ends of the t range, kept 77.8%).
 */
void
test_ncm_fftlog_sbessel_j_best_lnr0 (void)
{
  const guint ell_a[]  = {1, 2, 5, 10, 50, 200};
  const gdouble xs_a[] = {
    2.0815759778181006105, 3.3420936573656941588, 6.7564563302041293231,
    12.143204100943153341, 53.420758806667608602, 205.19128656443803093
  };
  guint j;

  for (j = 0; j < G_N_ELEMENTS (ell_a); j++)
  {
    NcmFftlogSBesselJ *fftlog_jl = ncm_fftlog_sbessel_j_new (ell_a[j], 0.0, 0.3, 20.0, 100);
    NcmFftlog *fftlog            = NCM_FFTLOG (fftlog_jl);

    ncm_fftlog_sbessel_j_set_best_lnr0 (fftlog_jl);
    ncm_assert_cmpdouble_e (exp (ncm_fftlog_get_lnr0 (fftlog) + 0.3), ==, xs_a[j], 1.0e-13, 0.0);

    ncm_fftlog_set_lnr0 (fftlog, -0.7);
    ncm_fftlog_sbessel_j_set_best_lnk0 (fftlog_jl);
    ncm_assert_cmpdouble_e (exp (ncm_fftlog_get_lnk0 (fftlog) - 0.7), ==, xs_a[j], 1.0e-13, 0.0);

    ncm_fftlog_free (fftlog);
  }

  {
    NcmFftlogSBesselJ *fftlog_jl = ncm_fftlog_sbessel_j_new (0, 0.0, 0.3, 20.0, 100);

    ncm_fftlog_sbessel_j_set_best_lnr0 (fftlog_jl);
    g_assert_cmpfloat (ncm_fftlog_get_lnr0 (NCM_FFTLOG (fftlog_jl)), ==, -0.3);
    ncm_fftlog_free (NCM_FFTLOG (fftlog_jl));
  }

  {
    const guint ell              = 10;
    const gdouble lnk0           = 0.0;
    const gdouble Lk             = log (1.0e8);
    const gdouble Y              = sqrt (M_PI) * pow (2.0, -1.5) * exp (lgamma (5.25) - lgamma (6.25));
    NcmFftlogSBesselJ *fftlog_jl = ncm_fftlog_sbessel_j_new (ell, 0.0, lnk0, Lk, 4000);
    NcmFftlog *fftlog            = NCM_FFTLOG (fftlog_jl);
    NcmVector *lnr, *Gr;
    guint i, good = 0;

    ncm_fftlog_sbessel_j_set_best_lnr0 (fftlog_jl);
    ncm_fftlog_set_padding (fftlog, 1.0);
    ncm_fftlog_use_smooth_padding (fftlog, TRUE);
    ncm_fftlog_eval_by_function (fftlog, &_test_power_m05, NULL);
    ncm_fftlog_set_bias (fftlog, ncm_fftlog_get_best_bias (fftlog));
    ncm_fftlog_eval_by_function (fftlog, &_test_power_m05, NULL);

    lnr = ncm_fftlog_get_vector_lnr (fftlog);
    Gr  = ncm_fftlog_get_vector_Gr (fftlog, 0);

    for (i = 0; i < ncm_vector_len (Gr); i++)
    {
      const gdouble truth = Y * exp (-0.5 * ncm_vector_get (lnr, i));

      if (fabs (ncm_vector_get (Gr, i) / truth - 1.0) < 1.0e-6)
        good++;
    }

    g_assert_cmpfloat (good, >, 0.9 * ncm_vector_len (Gr));

    ncm_vector_free (lnr);
    ncm_vector_free (Gr);
    ncm_fftlog_free (fftlog);
  }
}

