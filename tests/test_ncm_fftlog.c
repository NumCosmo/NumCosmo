/***************************************************************************
 *            test_ncm_fftlog.c
 *
 *  Sun September 03 11:45:13 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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
void test_ncm_fftlog_free (TestNcmFftlog *test, gconstpointer pdata);

void test_ncm_fftlog_eval (TestNcmFftlog *test, gconstpointer pdata);

void test_ncm_fftlog_tophatwin2_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_gausswin2_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_sbessel_j_traps (TestNcmFftlog *test, gconstpointer pdata);
void test_ncm_fftlog_invalid_st (TestNcmFftlog *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/fftlog/tophatwin2/eval", TestNcmFftlog, NULL,
              &test_ncm_fftlog_tophatwin2_new,
              &test_ncm_fftlog_eval,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/tophatwin2/traps", TestNcmFftlog, NULL,
              &test_ncm_fftlog_tophatwin2_new,
              &test_ncm_fftlog_tophatwin2_traps,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/gausswin2/eval", TestNcmFftlog, NULL,
              &test_ncm_fftlog_gausswin2_new,
              &test_ncm_fftlog_eval,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/gausswin2/traps", TestNcmFftlog, NULL,
              &test_ncm_fftlog_gausswin2_new,
              &test_ncm_fftlog_gausswin2_traps,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/sbessel_j/eval", TestNcmFftlog, NULL,
              &test_ncm_fftlog_sbessel_j_new,
              &test_ncm_fftlog_eval,
              &test_ncm_fftlog_free);

  g_test_add ("/ncm/fftlog/sbessel_j/traps", TestNcmFftlog, NULL,
              &test_ncm_fftlog_sbessel_j_new,
              &test_ncm_fftlog_sbessel_j_traps,
              &test_ncm_fftlog_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
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
#endif

#ifdef NUMCOSMO_HAVE_FFTW3  
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
  
  return GSL_FN_EVAL (&args->Fk, k) * k * exp (- kr * kr);
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
  const guint N          = g_test_rand_int_range  (1000, 2000);
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

  test->lnk_i        = g_test_rand_double_range (log (1.0e-4), log (1.0e0));
  test->lnk_f        = test->lnk_i + Lk;

  test->ntests       = NTESTS;

  arg->lnA           = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns            = g_test_rand_double_range (0.5, 1.5);

  argK->lnr          = 0.0;
  argK->Fk           = test->Fk;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_lnr0 (fftlog, -0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);
  
  g_assert (fftlog != NULL);
  g_assert (NCM_IS_FFTLOG (fftlog));
  g_assert (NCM_IS_FFTLOG_TOPHATWIN2 (fftlog));
}

void
test_ncm_fftlog_gausswin2_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (10000, 20000);
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

  test->lnk_i        = g_test_rand_double_range (log (1.0e-4), log (1.0e0));
  test->lnk_f        = test->lnk_i + Lk;

  test->ntests       = NTESTS;

  arg->lnA           = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns            = g_test_rand_double_range (0.5, 1.5);

  argK->lnr          = 0.0;
  argK->Fk           = test->Fk;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_lnr0 (fftlog, -0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);
  
  g_assert (fftlog != NULL);
  g_assert (NCM_IS_FFTLOG (fftlog));
  g_assert (NCM_IS_FFTLOG_GAUSSWIN2 (fftlog));
}

void
test_ncm_fftlog_sbessel_j_new (TestNcmFftlog *test, gconstpointer pdata)
{
  const guint N          = g_test_rand_int_range  (100000, 200000);
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

  test->lnk_i        = g_test_rand_double_range (log (1.0e-6), log (1.0e-4));
  test->lnk_f        = test->lnk_i + Lk;

  test->ntests       = NTESTS;

  arg->lnA           = g_test_rand_double_range (log (1.0e-10), log (1.0e+10));
  arg->ns            = g_test_rand_double_range (0.5, 1.5);

  argK->lnr          = 0.0;
  argK->Fk           = test->Fk;
  argK->ell          = ell;

  ncm_fftlog_set_lnk0 (fftlog, +0.5 * (test->lnk_i + test->lnk_f));
  ncm_fftlog_set_length (fftlog, Lk);

  ncm_fftlog_sbessel_j_set_best_lnr0 (NCM_FFTLOG_SBESSEL_J (fftlog));
  
  g_assert (fftlog != NULL);
  g_assert (NCM_IS_FFTLOG (fftlog));
  g_assert (NCM_IS_FFTLOG_SBESSEL_J (fftlog));
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
test_ncm_fftlog_eval (TestNcmFftlog *test, gconstpointer pdata)
{
  NcmFftlog *fftlog = test->fftlog;
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
    ncm_integral_locked_a_b (&test->KFk, test->lnk_i, test->lnk_f, 0.0, 1.0e-5, &res, &err);

    ncm_assert_cmpdouble_e (res, ==, fftlog_res, 1.0e-2, 0.0);
    /*printf ("%u % 22.15g % 22.15g % 22.15g % 22.15e\n", l, test->argK->lnr, res, fftlog_res, fabs (res / fftlog_res - 1.0));*/
  }
}

void
test_ncm_fftlog_tophatwin2_traps (TestNcmFftlog *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/ncm/fftlog/tophatwin2/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_fftlog_gausswin2_traps (TestNcmFftlog *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/ncm/fftlog/gausswin2/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_fftlog_sbessel_j_traps (TestNcmFftlog *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/ncm/fftlog/sbessel_j/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_fftlog_invalid_st (TestNcmFftlog *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}
