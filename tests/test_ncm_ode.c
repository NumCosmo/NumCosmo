/***************************************************************************
 *            test_ncm_ode.c
 *
 *  Mon January 07 20:00:57 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * test_ncm_ode.c
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
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

typedef struct _NcmODEEvalTestData
{
  gdouble w;
} NcmODEEvalTestData;

NCM_ODE_EVAL_DECLARE_IMPL (NcmODEEvalTest, ncm_ode_eval_test, NCM, ODE_EVAL_TEST, NcmODEEvalTestData)

gint
_ncm_ode_eval_test_df (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble * restrict df)
{
  NcmODEEvalTest *etest    = NCM_ODE_EVAL_TEST (ode_eval);
  NcmODEEvalTestData *data = _ncm_ode_eval_test_peek_ls (etest);
  
  df[0] = sin (data->w * t);
  
  return NCM_ODE_EVAL_RETURN_SUCCESS;
}

NCM_ODE_EVAL_DEFINE_IMPL (NcmODEEvalTest, ncm_ode_eval_test, NCM, ODE_EVAL_TEST, NcmODEEvalTestData, &_ncm_ode_eval_test_df, NULL, NULL)

typedef struct _TestNcmODE
{
  NcmODE *ode;
  NcmODEEval *ode_eval;
  guint ntests;
} TestNcmODE;

void test_ncm_ode_eval_test_new (TestNcmODE *test, gconstpointer pdata);
void test_ncm_ode_eval_free (TestNcmODE *test, gconstpointer pdata);
void test_ncm_ode_eval_test_df (TestNcmODE *test, gconstpointer pdata);

void test_ncm_ode_eval_test_traps (TestNcmODE *test, gconstpointer pdata);
void test_ncm_ode_eval_test_J (TestNcmODE *test, gconstpointer pdata);

void test_ncm_ode_new (TestNcmODE *test, gconstpointer pdata);
void test_ncm_ode_free (TestNcmODE *test, gconstpointer pdata);

void test_ncm_ode_traps (TestNcmODE *test, gconstpointer pdata);
void test_ncm_ode_invalid_st (TestNcmODE *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_set_nonfatal_assertions ();
  
  g_test_add ("/ncm/ode/eval/test/df", TestNcmODE, NULL,
              &test_ncm_ode_eval_test_new,
              &test_ncm_ode_eval_test_df,
              &test_ncm_ode_eval_free);
  
  g_test_add ("/ncm/ode/traps", TestNcmODE, NULL,
              &test_ncm_ode_new,
              &test_ncm_ode_traps,
              &test_ncm_ode_free);
  
  g_test_add ("/ncm/ode/eval/test/traps", TestNcmODE, NULL,
              &test_ncm_ode_eval_test_new,
              &test_ncm_ode_eval_test_traps,
              &test_ncm_ode_eval_free);
  
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/ode/invalid/st/subprocess", TestNcmODE, NULL,
              &test_ncm_ode_new,
              &test_ncm_ode_invalid_st,
              &test_ncm_ode_free);
  g_test_add ("/ncm/ode/eval/test/invalid/J/subprocess", TestNcmODE, NULL,
              &test_ncm_ode_eval_test_new,
              &test_ncm_ode_eval_test_J,
              &test_ncm_ode_eval_free);
#endif
  g_test_run ();
}

void
test_ncm_ode_new (TestNcmODE *test, gconstpointer pdata)
{
  /*NcmODE *ode = ncm_ode_new ();*/
  
  /*test->ode = ode;*/
  
  /*g_assert_true (ode != NULL);*/
  /*g_assert_true (NCM_IS_DIFF (ode));*/
}

void
test_ncm_ode_free (TestNcmODE *test, gconstpointer pdata)
{
  /*NcmODE *ode = test->ode;*/
  /*NCM_TEST_FREE (ncm_ode_free, ode);*/
}

void
test_ncm_ode_eval_test_new (TestNcmODE *test, gconstpointer pdata)
{
  NcmODEEvalTest *etest    = ncm_ode_eval_test_new ();
  NcmODEEvalTestData *data = _ncm_ode_eval_test_peek_ls (etest);
  
  test->ode_eval = NCM_ODE_EVAL (etest);
  data->w        = g_test_rand_double_range (-M_PI, +M_PI);
  
  g_assert_true (etest != NULL);
  g_assert_true (NCM_IS_ODE_EVAL_TEST (etest));
}

void
test_ncm_ode_eval_free (TestNcmODE *test, gconstpointer pdata)
{
  NcmODEEval *ode_eval = test->ode_eval;
  
  NCM_TEST_FREE (ncm_ode_eval_free, ode_eval);
}

void
test_ncm_ode_eval_test_df (TestNcmODE *test, gconstpointer pdata)
{
  NcmODEEvalTestData *data = _ncm_ode_eval_test_peek_ls (NCM_ODE_EVAL_TEST (test->ode_eval));
  gdouble t;
  
  for (t = 0.0; t <= 1.0; t += 0.01)
  {
    gdouble f  = 0.0;
    gdouble df = 0.0;
    
    ncm_ode_eval_df (test->ode_eval, 1, t, &f, &df);
    
    g_assert_cmpfloat (df, ==, sin (data->w * t));
  }
}

void
test_ncm_ode_traps (TestNcmODE *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/ode/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_ode_invalid_st (TestNcmODE *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

void
test_ncm_ode_eval_test_traps (TestNcmODE *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/ode/eval/test/invalid/J/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_ode_eval_test_J (TestNcmODE *test, gconstpointer pdata)
{
  const gdouble t = g_test_rand_double_range (0.0, 1.0);
  gdouble f = 0.0, *J, J_col;
  
  J = &J_col;
  
  ncm_ode_eval_J_dense (test->ode_eval, 1, t, &f, &J);
}

