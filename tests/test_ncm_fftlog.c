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

typedef struct _TestNcmFftlogK
{
  gsl_function Fk;
  gdouble lnr;
  gdouble w;
  guint ell;
  guint m;
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
void test_ncm_fftlog_sbessel_j_q0_5_new (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_sbessel_jljm_new (TestNcmFftlog *test, gconstpointer pdata);
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
void test_ncm_fftlog_sbessel_jljm_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_invalid_st (TestNcmFftlog *test, gconstpointer pdata);

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
  {"sbessel_j_q0_5", &test_ncm_fftlog_sbessel_j_q0_5_new},
#if defined (HAVE_FFTW3) && defined (HAVE_ACB_H)
  {"sbessel_jljm", &test_ncm_fftlog_sbessel_jljm_new},
#endif /* defined (HAVE_FFTW3) && defined (HAVE_ACB_H) */
};

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

/*
 *  g_test_add ("/ncm/fftlog/tophatwin2/traps", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_tophatwin2_new,
 *             &test_ncm_fftlog_tophatwin2_traps,
 *             &test_ncm_fftlog_free);
 *
 *  g_test_add ("/ncm/fftlog/gausswin2/traps", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_gausswin2_new,
 *             &test_ncm_fftlog_gausswin2_traps,
 *             &test_ncm_fftlog_free);
 *
 *  g_test_add ("/ncm/fftlog/sbessel_j/traps", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_sbessel_j_new,
 *             &test_ncm_fftlog_sbessel_j_traps,
 *             &test_ncm_fftlog_free);
 */
#if defined (HAVE_FFTW3) && defined (HAVE_ACB_H)

/*
 *  g_test_add ("/ncm/fftlog/sbessel_jljm/traps", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_sbessel_jljm_new,
 *             &test_ncm_fftlog_sbessel_jljm_traps,
 *             &test_ncm_fftlog_free);
 */
#endif /* defined (HAVE_FFTW3) && defined (HAVE_ACB_H) */

/*
 *  g_test_add ("/ncm/fftlog/tophatwin2/invalid/st/subprocess", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_tophatwin2_new,
 *             &test_ncm_fftlog_invalid_st,
 *             &test_ncm_fftlog_free);
 *  g_test_add ("/ncm/fftlog/gausswin2/invalid/st/subprocess", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_gausswin2_new,
 *             &test_ncm_fftlog_invalid_st,
 *             &test_ncm_fftlog_free);
 *  g_test_add ("/ncm/fftlog/sbessel_j/invalid/st/subprocess", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_sbessel_j_new,
 *             &test_ncm_fftlog_invalid_st,
 *             &test_ncm_fftlog_free);
 *  g_test_add ("/ncm/fftlog/sbessel_jljm/invalid/st/subprocess", TestNcmFftlog, NULL,
 *             &test_ncm_fftlog_sbessel_jljm_new,
 *             &test_ncm_fftlog_invalid_st,
 *             &test_ncm_fftlog_free);
 */
#ifdef HAVE_FFTW3
  g_test_run ();
#endif
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

static gdouble
_test_ncm_fftlog_sbessel_j_q0_5 (gdouble lnk, gpointer user_data)
{
  TestNcmFftlogK *args = (TestNcmFftlogK *) user_data;
  const gdouble kr     = exp (lnk + args->lnr);
  const gdouble k      = exp (lnk);

  return GSL_FN_EVAL (&args->Fk, k) * k * sqrt (kr) * ncm_sf_sbessel (args->ell, kr);
}

static gdouble
_test_ncm_fftlog_sbessel_jljm (gdouble lnk, gpointer user_data)
{
  TestNcmFftlogK *args = (TestNcmFftlogK *) user_data;
  const gdouble kr     = exp (lnk + args->lnr);
  const gdouble k      = exp (lnk);

/*
 *  printf ("% 22.15g % 22.15g\n", k,
 *       GSL_FN_EVAL (&args->Fk, k) * k * ncm_sf_sbessel (args->ell, kr * args->w) * ncm_sf_sbessel (args->m, kr / args->w));
 */
  return GSL_FN_EVAL (&args->Fk, k) * k * ncm_sf_sbessel (args->ell, kr * args->w) * ncm_sf_sbessel (args->m, kr / args->w);
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
  const guint N          = g_test_rand_int_range  (1000, 2000);
  const guint ell        = g_test_rand_int_range  (0, 10);
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (ell, 0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (1.0e+3), log (1.0e+4));

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
  arg->ns  = g_test_rand_double_range (0.5, 1.5);

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

void
test_ncm_fftlog_sbessel_j_q0_5_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (2000, 3000);
  const guint ell        = g_test_rand_int_range  (0, 10);
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (ell, 0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (1.0e+3), log (1.0e+4));

  test->fftlog       = fftlog;
  test->Fk.function  = &_test_ncm_fftlog_plaw;
  test->Fk.params    = arg;
  test->KFk.function = &_test_ncm_fftlog_sbessel_j_q0_5;
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
  ncm_fftlog_sbessel_j_set_q (NCM_FFTLOG_SBESSEL_J (fftlog), 0.5);

  ncm_fftlog_sbessel_j_set_best_lnr0 (NCM_FFTLOG_SBESSEL_J (fftlog));
  ncm_fftlog_sbessel_j_set_best_lnk0 (NCM_FFTLOG_SBESSEL_J (fftlog));

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));
  g_assert_true (NCM_IS_FFTLOG_SBESSEL_J (fftlog));
}

void
test_ncm_fftlog_sbessel_jljm_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = 1 * g_test_rand_int_range  (1000, 2000);
  const guint ell        = g_test_rand_int_range  (0, 10);
  const gint dell        = ell > 1 ? g_test_rand_int_range  (-2, 2) : g_test_rand_int_range  (-ell, ell + 2);
  const gdouble lnw      = 1.0 / 4.0 * log (g_test_rand_double_range (0.4, 1.0));
  NcmFftlog *fftlog      = NCM_FFTLOG (ncm_fftlog_sbessel_jljm_new (ell, dell, lnw, 0.0, 0.0, 20.0, N));
  TestNcmFftlogK *argK   = g_new (TestNcmFftlogK, 1);
  TestNcmFftlogPlaw *arg = g_new (TestNcmFftlogPlaw, 1);
  gdouble Lk             = g_test_rand_double_range (log (1.0e+6), log (1.0e+7));

  test->fftlog       = fftlog;
  test->Fk.function  = &_test_ncm_fftlog_plaw;
  test->Fk.params    = arg;
  test->KFk.function = &_test_ncm_fftlog_sbessel_jljm;
  test->KFk.params   = argK;
  test->argK         = argK;

  test->lnk_i = g_test_rand_double_range (log (1.0e-6), log (1.0e-4));
  test->lnk_f = test->lnk_i + Lk;

  test->ntests = NTESTS;

  arg->lnA = g_test_rand_double_range (log (1.0e-10), log (1.0e-9));
  arg->ns  = g_test_rand_double_range (0.92, 0.98);

  argK->lnr = 0.0;
  argK->Fk  = test->Fk;
  argK->ell = ell;
  argK->m   = ell + dell;
  argK->w   = exp (lnw);

  /*printf ("# %u %u % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", ell, N, Lk, argK->w, gsl_pow_4 (argK->w), exp (test->lnk_i), exp (test->lnk_f));*/

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);
  /*ncm_fftlog_set_lnr0 (fftlog, -0.5 * (test->lnk_i + test->lnk_f) + 3.0); */
  ncm_fftlog_sbessel_jljm_set_best_lnr0 (NCM_FFTLOG_SBESSEL_JLJM (fftlog));
  ncm_fftlog_sbessel_jljm_set_best_lnk0 (NCM_FFTLOG_SBESSEL_JLJM (fftlog));

  {
    const gdouble q = ncm_fftlog_sbessel_jljm_get_q (NCM_FFTLOG_SBESSEL_JLJM (fftlog));

    ncm_fftlog_sbessel_jljm_set_q (NCM_FFTLOG_SBESSEL_JLJM (fftlog), q);
    ncm_assert_cmpdouble_e (ncm_fftlog_sbessel_jljm_get_q (NCM_FFTLOG_SBESSEL_JLJM (fftlog)), ==, q, 1.0e-15, 0.0);
  }

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));
  g_assert_true (NCM_IS_FFTLOG_SBESSEL_JLJM (fftlog));
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

  ncm_fftlog_set_smooth_padding_scale (fftlog, 0.1);
  ncm_assert_cmpdouble_e (ncm_fftlog_get_smooth_padding_scale (fftlog), ==, 0.1, 1.0e-15, 0.0);

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
  guint nerr        = 0;
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-3, &res, &err);

    if ((fabs (fftlog_res / res - 1.0) > reltol) && (nerr < 15))
      nerr++;
    else
      ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);

    /* printf ("%u % 22.15g % 22.15g % 22.15g % 22.15e\n", l, test->argK->lnr, res, fftlog_res, fabs (res / fftlog_res - 1.0)); */
  }
}

void
test_ncm_fftlog_eval_vector (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  guint nerr        = 0;
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-3, &res, &err);

    if ((fabs (fftlog_res / res - 1.0) > reltol) && (nerr < 15))
      nerr++;
    else
      ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);

    /* printf ("%u % 22.15g % 22.15g % 22.15g % 22.15e\n", l, test->argK->lnr, res, fftlog_res, fabs (res / fftlog_res - 1.0)); */
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-3, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);
  }
}

void
test_ncm_fftlog_eval_calibrate_fail (TestNcmFftlog *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    NcmFftlog *fftlog = test->fftlog;

    ncm_fftlog_set_max_size (fftlog, 100);

    ncm_fftlog_calibrate_size (fftlog, test->Fk.function, test->Fk.params, 1.0e-1);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*maximum number of knots reached. Requested precision*");
}

void
test_ncm_fftlog_eval_serialized (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_serialize_dup_obj (ser, G_OBJECT (test->fftlog)));
  guint nerr        = 0;
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-3, &res, &err);

    if ((fabs (fftlog_res / res - 1.0) > reltol) && (nerr < 15))
      nerr++;
    else
      ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);

    /* printf ("%u % 22.15g % 22.15g % 22.15g % 22.15e\n", l, test->argK->lnr, res, fftlog_res, fabs (res / fftlog_res - 1.0)); */
  }

  ncm_serialize_free (ser);
  ncm_fftlog_free (fftlog);
}

void
test_ncm_fftlog_eval_deriv (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  guint nerr        = 0;
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-3, &res, &err);

    if ((fabs (fftlog_res / res - 1.0) > reltol) && (nerr < 15))
      nerr++;
    else
      ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);

    /* printf ("%u % 22.15g % 22.15g % 22.15g % 22.15e\n", l, test->argK->lnr, res, fftlog_res, fabs (res / fftlog_res - 1.0)); */
  }
}

void
test_ncm_fftlog_eval_use_eval_int (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  guint nerr        = 0;
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-3, &res, &err);

    if ((fabs (fftlog_res / res - 1.0) > reltol) && (nerr < 15))
      nerr++;
    else
      ncm_assert_cmpdouble_e (res, ==, fftlog_res, reltol, 0.0);

    /* printf ("%u % 22.15g % 22.15g % 22.15g % 22.15e\n", l, test->argK->lnr, res, fftlog_res, fabs (res / fftlog_res - 1.0)); */
  }
}

void
test_ncm_fftlog_eval_smooth_padding (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
  guint nerr        = 0;
  gdouble reltol    = 1.0e-1;
  NcmVector *lnr;
  guint i, len;

  ncm_fftlog_set_smooth_padding_scale (fftlog, 1.0e-4);
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
test_ncm_fftlog_sbessel_jljm_traps (TestNcmFftlog *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/fftlog/sbessel_jljm/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fftlog_invalid_st (TestNcmFftlog *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

