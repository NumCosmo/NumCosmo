/***************************************************************************
 *            test_ncm_diff.c
 *
 *  Wed July 26 12:04:44 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * test_ncm_diff.c
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

typedef struct _TestNcmDiff
{
  NcmDiff *diff;
  guint ntests;
} TestNcmDiff;

void test_ncm_diff_new (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_free (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_sin (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_sin (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_sin (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_asin (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_asin (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_asin (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_tan (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_tan (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_tan (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_exp (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_exp (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_exp (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_log (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_log (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_log (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_poly3 (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_poly3 (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_poly3 (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_1_plaw (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_1_plaw (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_1_plaw (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_1_to_M_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_1_to_M_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_1_to_M_all (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_N_to_1_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_N_to_1_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_N_to_1_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rf_Hessian_N_to_1_all (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_rf_d1_N_to_M_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d1_N_to_M_all (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rc_d2_N_to_M_all (TestNcmDiff *test, gconstpointer pdata);

void test_ncm_diff_traps (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_invalid_st (TestNcmDiff *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_set_nonfatal_assertions ();
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/sin", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_sin,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/sin", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_sin,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/sin", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_sin,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/asin", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_asin,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/asin", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_asin,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/asin", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_asin,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/tan", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_tan,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/tan", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_tan,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/tan", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_tan,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/exp", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_exp,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/exp", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_exp,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/exp", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_exp,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/log", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_log,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/log", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_log,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/log", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_log,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/poly3", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_poly3,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/poly3", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_poly3,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/poly3", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_poly3,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_1/plaw", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_1_plaw,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_1/plaw", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_1_plaw,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_1/plaw", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_1_plaw,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/1_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_1_to_M_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/1_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_1_to_M_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/1_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_1_to_M_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_N_to_1_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_N_to_1_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_N_to_1_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/Hessian/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_Hessian_N_to_1_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rf/d1/N_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_N_to_M_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d1/N_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_N_to_M_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/rc/d2/N_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_N_to_M_all,
              &test_ncm_diff_free);
  
  g_test_add ("/ncm/diff/traps", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_traps,
              &test_ncm_diff_free);
  
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/diff/invalid/st/subprocess", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_invalid_st,
              &test_ncm_diff_free);
#endif
  g_test_run ();
}

void
test_ncm_diff_new (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = ncm_diff_new ();
  
  ncm_diff_set_round_off_pad (diff, 2.0e0);
  ncm_diff_set_trunc_error_pad (diff, 2.0e0);
  test->diff = diff;
  
  g_assert_true (diff != NULL);
  g_assert_true (NCM_IS_DIFF (diff));
}

void
test_ncm_diff_free (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  
  NCM_TEST_FREE (ncm_diff_free, diff);
}

/*
 * SIN
 */

static gdouble
_test_ncm_diff_sin (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return sin (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_dsin (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] * cos (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_d2sin (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return -gsl_pow_2 (w_ptr[0]) * sin (x * w_ptr[0]);
}

/*
 * ASIN
 */

static gdouble
_test_ncm_diff_asin (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;

  return asin (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_dasin (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] / sqrt (1.0 - gsl_pow_2 (x * w_ptr[0]));
}

static gdouble
_test_ncm_diff_d2asin (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return gsl_pow_3 (w_ptr[0]) * x / gsl_pow_3 (sqrt (1.0 - gsl_pow_2 (x * w_ptr[0])));
}

/*
 * TAN
 */
static gdouble
_test_ncm_diff_tan (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return tan (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_dtan (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] / gsl_pow_2 (cos (x * w_ptr[0]));
}

static gdouble
_test_ncm_diff_d2tan (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return 2.0 * gsl_pow_2 (w_ptr[0]) * tan (x * w_ptr[0]) / gsl_pow_2 (cos (x * w_ptr[0]));
}

/*
 * EXP
 */

static gdouble
_test_ncm_diff_exp (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return exp (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_dexp (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] * exp (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_d2exp (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return gsl_pow_2 (w_ptr[0]) * exp (x * w_ptr[0]);
}

/*
 * LOG
 */

static gdouble
_test_ncm_diff_log (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return log (x * w_ptr[0]);
}

static gdouble
_test_ncm_diff_dlog (const gdouble x, gpointer userdata)
{
  return 1.0 / x;
}

static gdouble
_test_ncm_diff_d2log (const gdouble x, gpointer userdata)
{
  return -1.0 / gsl_pow_2 (x);
}

/*
 * POLY3
 */
static gdouble
_test_ncm_diff_poly3 (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] + x * w_ptr[1] + x * x * w_ptr[2] + x * x * x * w_ptr[3];
}

static gdouble
_test_ncm_diff_dpoly3 (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[1] + 2.0 * x * w_ptr[2] + 3.0 * x * x * w_ptr[3];
}

static gdouble
_test_ncm_diff_d2poly3 (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return 2.0 * w_ptr[2] + 6.0 * x * w_ptr[3];
}

/*
 * PLAW
 */
static gdouble
_test_ncm_diff_plaw (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return pow (x, w_ptr[0]);
}

static gdouble
_test_ncm_diff_dplaw (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] * pow (x, w_ptr[0] - 1.0);
}

static gdouble
_test_ncm_diff_d2plaw (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;
  
  return w_ptr[0] * (w_ptr[0] - 1.0) * pow (x, w_ptr[0] - 2.0);
}

/*
 * 1 to M
 * ALL
 */

typedef struct _TestNcmDiffAll
{
  gdouble w_sin;
  gdouble w_asin;
  gdouble w_tan;
  gdouble w_exp;
  gdouble w_log;
  gdouble w_poly3[4];
  gdouble w_plaw;
} TestNcmDiffAll;

static void
_test_ncm_diff_all (const gdouble x, NcmVector *y, gpointer userdata)
{
  TestNcmDiffAll *arg = (TestNcmDiffAll *) userdata;
  
  g_assert_cmpuint (ncm_vector_len (y), ==, 7);
  
  ncm_vector_set (y, 0, _test_ncm_diff_sin   (x, &arg->w_sin));
  ncm_vector_set (y, 1, _test_ncm_diff_asin  (x, &arg->w_asin));
  ncm_vector_set (y, 2, _test_ncm_diff_tan   (x, &arg->w_tan));
  ncm_vector_set (y, 3, _test_ncm_diff_exp   (x, &arg->w_exp));
  ncm_vector_set (y, 4, _test_ncm_diff_log   (x, &arg->w_log));
  ncm_vector_set (y, 5, _test_ncm_diff_poly3 (x,  arg->w_poly3));
  ncm_vector_set (y, 6, _test_ncm_diff_plaw  (x, &arg->w_plaw));
}

static GArray *
_test_ncm_diff_dall (const gdouble x, gpointer userdata)
{
  TestNcmDiffAll *arg = (TestNcmDiffAll *) userdata;
  GArray *y_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *y        = NULL;
  
  g_array_set_size (y_a, 7);
  y = ncm_vector_new_array (y_a);
  
  ncm_vector_set (y, 0, _test_ncm_diff_dsin   (x, &arg->w_sin));
  ncm_vector_set (y, 1, _test_ncm_diff_dasin  (x, &arg->w_asin));
  ncm_vector_set (y, 2, _test_ncm_diff_dtan   (x, &arg->w_tan));
  ncm_vector_set (y, 3, _test_ncm_diff_dexp   (x, &arg->w_exp));
  ncm_vector_set (y, 4, _test_ncm_diff_dlog   (x, &arg->w_log));
  ncm_vector_set (y, 5, _test_ncm_diff_dpoly3 (x,  arg->w_poly3));
  ncm_vector_set (y, 6, _test_ncm_diff_dplaw  (x, &arg->w_plaw));
  
  ncm_vector_free (y);
  
  return y_a;
}

static GArray *
_test_ncm_diff_d2all (const gdouble x, gpointer userdata)
{
  TestNcmDiffAll *arg = (TestNcmDiffAll *) userdata;
  GArray *y_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *y        = NULL;
  
  g_array_set_size (y_a, 7);
  y = ncm_vector_new_array (y_a);
  
  ncm_vector_set (y, 0, _test_ncm_diff_d2sin   (x, &arg->w_sin));
  ncm_vector_set (y, 1, _test_ncm_diff_d2asin  (x, &arg->w_asin));
  ncm_vector_set (y, 2, _test_ncm_diff_d2tan   (x, &arg->w_tan));
  ncm_vector_set (y, 3, _test_ncm_diff_d2exp   (x, &arg->w_exp));
  ncm_vector_set (y, 4, _test_ncm_diff_d2log   (x, &arg->w_log));
  ncm_vector_set (y, 5, _test_ncm_diff_d2poly3 (x,  arg->w_poly3));
  ncm_vector_set (y, 6, _test_ncm_diff_d2plaw  (x, &arg->w_plaw));
  
  ncm_vector_free (y);
  
  return y_a;
}

/*
 * N to 1
 * ALL
 */

static gdouble
_test_ncm_diff_N_to_1_all (NcmVector *x, gpointer userdata)
{
  gdouble *w = (gdouble *) userdata;
  
  g_assert_cmpuint (ncm_vector_len (x), ==, 3);
  
  {
    const gdouble v1 = ncm_vector_get (x, 0);
    const gdouble v2 = ncm_vector_get (x, 1);
    const gdouble v3 = ncm_vector_get (x, 2);
    
    return sin (v1 * v2 * w[0]) * exp (v3 * w[1]);
  }
}

static GArray *
_test_ncm_diff_N_to_1_dall (GArray *x_a, gpointer userdata)
{
  gdouble *w   = (gdouble *) userdata;
  GArray *y_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *y = NULL;
  
  g_array_set_size (y_a, 3);
  y = ncm_vector_new_array (y_a);
  
  g_assert_cmpuint (x_a->len, ==, 3);
  
  {
    const gdouble v1 = g_array_index (x_a, gdouble, 0);
    const gdouble v2 = g_array_index (x_a, gdouble, 1);
    const gdouble v3 = g_array_index (x_a, gdouble, 2);
    
    ncm_vector_set (y, 0, v2 * w[0] * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 1, v1 * w[0] * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 2,      w[1] * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
  }
  
  ncm_vector_free (y);
  
  return y_a;
}

static GArray *
_test_ncm_diff_N_to_1_d2all (GArray *x_a, gpointer userdata)
{
  gdouble *w   = (gdouble *) userdata;
  GArray *y_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *y = NULL;
  
  g_array_set_size (y_a, 3);
  y = ncm_vector_new_array (y_a);
  
  g_assert_cmpuint (x_a->len, ==, 3);
  
  {
    const gdouble v1 = g_array_index (x_a, gdouble, 0);
    const gdouble v2 = g_array_index (x_a, gdouble, 1);
    const gdouble v3 = g_array_index (x_a, gdouble, 2);
    
    ncm_vector_set (y, 0, -gsl_pow_2 (v2 * w[0]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 1, -gsl_pow_2 (v1 * w[0]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 2,       gsl_pow_2 (w[1]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
  }
  
  ncm_vector_free (y);
  
  return y_a;
}

static GArray *
_test_ncm_diff_N_to_1_Hessian_all (GArray *x_a, gpointer userdata)
{
  gdouble *w   = (gdouble *) userdata;
  GArray *y_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmMatrix *y = NULL;
  
  g_array_set_size (y_a, 3 * 3);
  y = ncm_matrix_new_array (y_a, 3);
  
  g_assert_cmpuint (x_a->len, ==, 3);
  
  {
    const gdouble v1 = g_array_index (x_a, gdouble, 0);
    const gdouble v2 = g_array_index (x_a, gdouble, 1);
    const gdouble v3 = g_array_index (x_a, gdouble, 2);
    
    ncm_matrix_set (y, 0, 0, -gsl_pow_2 (v2 * w[0]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_matrix_set (y, 0, 1, -w[0] * exp (v3 * w[1]) * (w[0] * v1 * v2 * sin (v1 * v2 * w[0]) - cos (v1 * v2 * w[0])));
    ncm_matrix_set (y, 0, 2, w[0] * w[1] * v2 * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    
    ncm_matrix_set (y, 1, 0, -w[0] * exp (v3 * w[1]) * (w[0] * v1 * v2 * sin (v1 * v2 * w[0]) - cos (v1 * v2 * w[0])));
    ncm_matrix_set (y, 1, 1, -gsl_pow_2 (v1 * w[0]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_matrix_set (y, 1, 2, w[0] * w[1] * v1 * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    
    ncm_matrix_set (y, 2, 0, w[0] * w[1] * v2 * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_matrix_set (y, 2, 1, w[0] * w[1] * v1 * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_matrix_set (y, 2, 2,       gsl_pow_2 (w[1]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
  }
  
  ncm_matrix_free (y);
  
  return y_a;
}

/*
 * N to M
 * ALL
 */

static void
_test_ncm_diff_N_to_M_all (NcmVector *x, NcmVector *y, gpointer userdata)
{
  gdouble *w = (gdouble *) userdata;
  
  g_assert_cmpuint (ncm_vector_len (x), ==, 3);
  g_assert_cmpuint (ncm_vector_len (y), ==, 3);
  
  {
    const gdouble v1 = ncm_vector_get (x, 0);
    const gdouble v2 = ncm_vector_get (x, 1);
    const gdouble v3 = ncm_vector_get (x, 2);
    
    ncm_vector_set (y, 0, sin (v1 * v2 * w[0]) * exp (+v3 * w[1]));
    ncm_vector_set (y, 1, cos (v1 * v2 * w[0]) * exp (+v3 * w[1]));
    ncm_vector_set (y, 2, cos (v1 * v2 * w[0]) * exp (-v3 * w[1]));
  }
}

static GArray *
_test_ncm_diff_N_to_M_dall (GArray *x_a, gpointer userdata)
{
  gdouble *w   = (gdouble *) userdata;
  GArray *y_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *y = NULL;
  
  g_array_set_size (y_a, 3 * 3);
  y = ncm_vector_new_array (y_a);
  
  g_assert_cmpuint (x_a->len, ==, 3);
  
  {
    const gdouble v1 = g_array_index (x_a, gdouble, 0);
    const gdouble v2 = g_array_index (x_a, gdouble, 1);
    const gdouble v3 = g_array_index (x_a, gdouble, 2);
    
    ncm_vector_set (y, 0, v2 * w[0] * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 3, v1 * w[0] * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 6,      w[1] * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    
    ncm_vector_set (y, 1, -v2 * w[0] * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 4, -v1 * w[0] * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 7,       w[1] * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    
    ncm_vector_set (y, 2, -v2 * w[0] * sin (v1 * v2 * w[0]) * exp (-v3 * w[1]));
    ncm_vector_set (y, 5, -v1 * w[0] * sin (v1 * v2 * w[0]) * exp (-v3 * w[1]));
    ncm_vector_set (y, 8,     -w[1] * cos (v1 * v2 * w[0]) * exp (-v3 * w[1]));
  }
  
  ncm_vector_free (y);
  
  return y_a;
}

static GArray *
_test_ncm_diff_N_to_M_d2all (GArray *x_a, gpointer userdata)
{
  gdouble *w   = (gdouble *) userdata;
  GArray *y_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *y = NULL;
  
  g_array_set_size (y_a, 3 * 3);
  y = ncm_vector_new_array (y_a);
  
  g_assert_cmpuint (x_a->len, ==, 3);
  
  {
    const gdouble v1 = g_array_index (x_a, gdouble, 0);
    const gdouble v2 = g_array_index (x_a, gdouble, 1);
    const gdouble v3 = g_array_index (x_a, gdouble, 2);
    
    ncm_vector_set (y, 0, -gsl_pow_2 (v2 * w[0]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 3, -gsl_pow_2 (v1 * w[0]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 6,       gsl_pow_2 (w[1]) * sin (v1 * v2 * w[0]) * exp (v3 * w[1]));
    
    ncm_vector_set (y, 1, -gsl_pow_2 (v2 * w[0]) * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 4, -gsl_pow_2 (v1 * w[0]) * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    ncm_vector_set (y, 7,       gsl_pow_2 (w[1]) * cos (v1 * v2 * w[0]) * exp (v3 * w[1]));
    
    ncm_vector_set (y, 2, -gsl_pow_2 (v2 * w[0]) * cos (v1 * v2 * w[0]) * exp (-v3 * w[1]));
    ncm_vector_set (y, 5, -gsl_pow_2 (v1 * w[0]) * cos (v1 * v2 * w[0]) * exp (-v3 * w[1]));
    ncm_vector_set (y, 8,       gsl_pow_2 (w[1]) * cos (v1 * v2 * w[0]) * exp (-v3 * w[1]));
  }
  
  ncm_vector_free (y);
  
  return y_a;
}

/*
 * END FUNCS
 */

/*
 * SIN
 */

void
test_ncm_diff_rf_d1_1_to_1_sin (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0, 5.0);
    const gdouble x   = g_test_rand_double_range (-100.0, 100.0);
    const gdouble df  = ncm_diff_rf_d1_1_to_1 (diff, x, &_test_ncm_diff_sin, &w, &err);
    const gdouble Adf = _test_ncm_diff_dsin (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_sin (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0, 5.0);
    const gdouble x   = g_test_rand_double_range (-100.0, 100.0);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_sin, &w, &err);
    const gdouble Adf = _test_ncm_diff_dsin (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_sin (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0, 5.0);
    const gdouble x   = g_test_rand_double_range (-100.0, 100.0);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_sin, &w, &err);
    const gdouble Adf = _test_ncm_diff_d2sin (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * ASIN
 */

void
test_ncm_diff_rf_d1_1_to_1_asin (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-2, 1.0);
    const gdouble x   = g_test_rand_double_range (-0.95, 0.95);
    const gdouble df  = ncm_diff_rf_d1_1_to_1 (diff, x, &_test_ncm_diff_asin, &w, &err);
    const gdouble Adf = _test_ncm_diff_dasin (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_asin (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  gint nerr     = 1500;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-2, 1.0);
    const gdouble x   = g_test_rand_double_range (-0.95, 0.95);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_asin, &w, &err);
    const gdouble Adf = _test_ncm_diff_dasin (x, &w);

    if (((err == 0.0) || gsl_isnan (err)) && nerr)
    {
      nerr--;
      g_test_skip ("Unable to estimate error.");
      continue;
    }

    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_asin (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  gint nerr     = 5;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-2, 1.0);
    const gdouble x   = g_test_rand_double_range (-0.95, 0.95);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_asin, &w, &err);
    const gdouble Adf = _test_ncm_diff_d2asin (x, &w);
    
    /*printf ("%d %d % 22.15g % 22.15g % 22.15g % 22.15g\n", i, nerr, err, x, df, Adf);*/

    if (((err == 0.0) || gsl_isnan (err)) && nerr)
    {
      nerr--;
      g_test_skip ("Unable to estimate error.");
      continue;
    }

    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * TAN
 */

void
test_ncm_diff_rf_d1_1_to_1_tan (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-1, 0.5 * M_PI);
    const gdouble x   = g_test_rand_double_range (-1.0, 1.0);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_tan, &w, &err);
    const gdouble Adf = _test_ncm_diff_dtan (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_tan (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-1, 0.5 * M_PI);
    const gdouble x   = g_test_rand_double_range (-1.0, 1.0);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_tan, &w, &err);
    const gdouble Adf = _test_ncm_diff_dtan (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_tan (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  gint nerr     = 5;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-1, 0.5 * M_PI);
    const gdouble x   = g_test_rand_double_range (-1.0, 1.0);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_tan, &w, &err);
    const gdouble Adf = _test_ncm_diff_d2tan (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    if (((err == 0.0) || gsl_isnan (err)) && nerr)
    {
      nerr--;
      g_test_skip ("Unable to estimate error.");
      continue;
    }

    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * EXP
 */

void
test_ncm_diff_rf_d1_1_to_1_exp (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e2);
    const gdouble x   = g_test_rand_double_range (-1.0, 1.0);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_exp, &w, &err);
    const gdouble Adf = _test_ncm_diff_dexp (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_exp (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e2);
    const gdouble x   = g_test_rand_double_range (-1.0, 1.0);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_exp, &w, &err);
    const gdouble Adf = _test_ncm_diff_dexp (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_exp (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e2);
    const gdouble x   = g_test_rand_double_range (-1.0, 1.0);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_exp, &w, &err);
    const gdouble Adf = _test_ncm_diff_d2exp (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * LOG
 */

void
test_ncm_diff_rf_d1_1_to_1_log (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e2);
    const gdouble x   = g_test_rand_double_range (1.0e-5, 1.0e5);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_log, &w, &err);
    const gdouble Adf = _test_ncm_diff_dlog (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_log (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e2);
    const gdouble x   = g_test_rand_double_range (1.0e-5, 1.0e5);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_log, &w, &err);
    const gdouble Adf = _test_ncm_diff_dlog (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_log (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e2);
    const gdouble x   = g_test_rand_double_range (1.0e-5, 1.0e5);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_log, &w, &err);
    const gdouble Adf = _test_ncm_diff_d2log (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * POLY3
 */

void
test_ncm_diff_rf_d1_1_to_1_poly3 (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[4]      = {g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2)};
    const gdouble x   = g_test_rand_double_range (-1.0e3, 1.0e3);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_poly3, w, &err);
    const gdouble Adf = _test_ncm_diff_dpoly3 (x, w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_poly3 (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[4]      = {g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2)};
    const gdouble x   = g_test_rand_double_range (-1.0e3, 1.0e3);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_poly3, w, &err);
    const gdouble Adf = _test_ncm_diff_dpoly3 (x, w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_poly3 (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[4]      = {g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2), g_test_rand_double_range (-1.0e2, 1.0e2)};
    const gdouble x   = g_test_rand_double_range (-1.0e3, 1.0e3);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_poly3, w, &err);
    const gdouble Adf = _test_ncm_diff_d2poly3 (x, w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * PLAW
 */

void
test_ncm_diff_rf_d1_1_to_1_plaw (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e1);
    const gdouble x   = g_test_rand_double_range (1.0e-3, 1.0e3);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_plaw, &w, &err);
    const gdouble Adf = _test_ncm_diff_dplaw (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d1_1_to_1_plaw (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e1);
    const gdouble x   = g_test_rand_double_range (1.0e-3, 1.0e3);
    const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, x, &_test_ncm_diff_plaw, &w, &err);
    const gdouble Adf = _test_ncm_diff_dplaw (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

void
test_ncm_diff_rc_d2_1_to_1_plaw (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble err   = 0.0;
  guint ntests  = 1000;
  guint i;
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w         = g_test_rand_double_range (1.0e-3, 1.0e1);
    const gdouble x   = g_test_rand_double_range (1.0e-3, 1.0e3);
    const gdouble df  = ncm_diff_rc_d2_1_to_1 (diff, x, &_test_ncm_diff_plaw, &w, &err);
    const gdouble Adf = _test_ncm_diff_d2plaw (x, &w);
    
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, Adf, df, df / Adf - 1.0, err);*/
    ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
  }
}

/*
 * 1 to M
 * ALL
 */

void
test_ncm_diff_rf_d1_1_to_M_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 5;
  guint i, j;
  
  for (i = 0; i < ntests; i++)
  {
    TestNcmDiffAll arg =
    {
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (-0.95,         0.95),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI),
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (1.0e-3,    1.0e3),
      {
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2)
      },
      g_test_rand_double_range (1.0e-3,    1.0e2),
    };
    const gdouble x = g_test_rand_double_range (1.0e-3, 1.0);
    GArray *df_a    = ncm_diff_rf_d1_1_to_M (diff, x, 7, &_test_ncm_diff_all, &arg, &err_a);
    GArray *Adf_a   = _test_ncm_diff_dall (x, &arg);
    
    for (j = 0; j < 7; j++)
    {
      const gdouble df  = g_array_index (df_a,  gdouble, j);
      const gdouble Adf = g_array_index (Adf_a, gdouble, j);
      const gdouble err = g_array_index (err_a, gdouble, j);
      
      /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g\n", x, df, Adf, err);*/

      if (((err == 0.0) || gsl_isnan (err)) && nerr)
      {
        nerr--;
        g_test_skip ("Unable to estimate error.");
        continue;
      }

      ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
    }
    
    g_array_unref (df_a);
    g_array_unref (Adf_a);
    g_array_unref (err_a);
  }
}

void
test_ncm_diff_rc_d1_1_to_M_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 5;
  guint i, j;
  
  for (i = 0; i < ntests; i++)
  {
    TestNcmDiffAll arg =
    {
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (-0.95,         0.95),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI),
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (1.0e-3,    1.0e3),
      {
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2)
      },
      g_test_rand_double_range (1.0e-3,    1.0e2),
    };
    const gdouble x = g_test_rand_double_range (1.0e-3, 1.0);
    GArray *df_a    = ncm_diff_rc_d1_1_to_M (diff, x, 7, &_test_ncm_diff_all, &arg, &err_a);
    GArray *Adf_a   = _test_ncm_diff_dall (x, &arg);
    
    for (j = 0; j < 7; j++)
    {
      const gdouble df  = g_array_index (df_a,  gdouble, j);
      const gdouble Adf = g_array_index (Adf_a, gdouble, j);
      const gdouble err = g_array_index (err_a, gdouble, j);
      
      if (((err == 0.0) || gsl_isnan (err)) && nerr)
      {
        nerr--;
        g_test_skip ("Unable to estimate error.");
        continue;
      }

      ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
    }
    
    g_array_unref (df_a);
    g_array_unref (Adf_a);
    g_array_unref (err_a);
  }
}

void
test_ncm_diff_rc_d2_1_to_M_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 5;
  guint i, j;
  
  for (i = 0; i < ntests; i++)
  {
    TestNcmDiffAll arg =
    {
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (-0.95,         0.95),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI),
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (1.0e-3,    1.0e3),
      {
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2),
        g_test_rand_double_range (-1.0e2, 1.0e2)
      },
      g_test_rand_double_range (1.0e-3,    1.0e2),
    };
    const gdouble x = g_test_rand_double_range (1.0e-3, 1.0);
    GArray *df_a    = ncm_diff_rc_d2_1_to_M (diff, x, 7, &_test_ncm_diff_all, &arg, &err_a);
    GArray *Adf_a   = _test_ncm_diff_d2all (x, &arg);
    
    for (j = 0; j < 7; j++)
    {
      const gdouble df  = g_array_index (df_a,  gdouble, j);
      const gdouble Adf = g_array_index (Adf_a, gdouble, j);
      const gdouble err = g_array_index (err_a, gdouble, j);

      /*printf ("[%2u %2d] % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", j, nerr, x, Adf, df, df / Adf - 1.0, err);*/
      if (((err == 0.0) || gsl_isnan (err)) && nerr)
      {
        nerr--;
        g_test_skip ("Unable to estimate error.");
        continue;
      }

      ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
    }
    
    g_array_unref (df_a);
    g_array_unref (Adf_a);
    g_array_unref (err_a);
  }
}

/*
 * N to 1
 * ALL
 */

void
test_ncm_diff_rf_d1_N_to_1_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  guint i, j;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range ( -0.99,         0.99),
      g_test_rand_double_range ( -0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      GArray *df_a  = ncm_diff_rf_d1_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_dall (x_a, w);
      
      for (j = 0; j < x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

void
test_ncm_diff_rc_d1_N_to_1_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  guint i, j;
  gint nerr = 5;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range ( -0.99,         0.99),
      g_test_rand_double_range ( -0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      GArray *df_a  = ncm_diff_rc_d1_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_dall (x_a, w);
      
      for (j = 0; j < x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);
        
        if (((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

void
test_ncm_diff_rc_d2_N_to_1_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 5;
  guint i, j;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range ( -0.99,         0.99),
      g_test_rand_double_range ( -0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      GArray *df_a  = ncm_diff_rc_d2_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_d2all (x_a, w);
      
      for (j = 0; j < x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);
        
        if (((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

void
test_ncm_diff_rf_Hessian_N_to_1_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  guint nerr = 0;
  guint i, j;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range ( -0.99,         0.99),
      g_test_rand_double_range ( -0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      GArray *df_a  = ncm_diff_rf_Hessian_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_Hessian_all (x_a, w);
      
      for (j = 0; j < x_a->len * x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);
        
/*
 *       printf ("[%u, %u] (% 22.15g % 22.15g % 22.15g) % 22.15g % 22.15g % 22.15g % 22.15g [% 22.15g % 22.15g % 22.15g]\n",
 *               j / 3, j % 3, v1, v2, v3, Adf, df, df - Adf, err,
 *               w[0], w[1], w[2]);
 */
        if ((fabs (df - Adf) > err) && (nerr < 10))
          nerr++;
        else
          ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

/*
 * N to M
 * ALL
 */

void
test_ncm_diff_rf_d1_N_to_M_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 5;
  guint i, j;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-100.0,       100.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-10.0, 10.0);
    const gdouble v2 = g_test_rand_double_range (-10.0, 10.0);
    const gdouble v3 = g_test_rand_double_range (-10.0, 10.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      const guint dim = 3;
      GArray *df_a    = ncm_diff_rf_d1_N_to_M (diff, x_a, dim, &_test_ncm_diff_N_to_M_all, w, &err_a);
      GArray *Adf_a   = _test_ncm_diff_N_to_M_dall (x_a, w);
      
      for (j = 0; j < x_a->len * dim; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        if (((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

void
test_ncm_diff_rc_d1_N_to_M_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 5;
  guint i, j;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range ( -0.99,         0.99),
      g_test_rand_double_range ( -0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      const guint dim = 3;
      GArray *df_a    = ncm_diff_rc_d1_N_to_M (diff, x_a, dim, &_test_ncm_diff_N_to_M_all, w, &err_a);
      GArray *Adf_a   = _test_ncm_diff_N_to_M_dall (x_a, w);
      
      for (j = 0; j < x_a->len * dim; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);
        
        if (((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

void
test_ncm_diff_rc_d2_N_to_M_all (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  GArray *x_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a = NULL;
  guint ntests = 1000;
  gint nerr = 15;
  guint i, j;
  
  g_array_set_size (x_a, 3);
  
  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-1.0,          1.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };
    
    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);
    
    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;
    
    {
      const guint dim = 3;
      GArray *df_a    = ncm_diff_rc_d2_N_to_M (diff, x_a, dim, &_test_ncm_diff_N_to_M_all, w, &err_a);
      GArray *Adf_a   = _test_ncm_diff_N_to_M_d2all (x_a, w);
      
      for (j = 0; j < x_a->len * dim; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        /*printf ("(% 22.15g % 22.15g % 22.15g) % 22.15g % 22.15g % 22.15g\n", v1, v2, v3, df, Adf, err);*/
        if (((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
      }
      
      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }
  
  g_array_unref (x_a);
}

void
test_ncm_diff_traps (TestNcmDiff *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/diff/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_diff_invalid_st (TestNcmDiff *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

