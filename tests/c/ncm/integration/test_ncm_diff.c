/***************************************************************************
 *            test_ncm_diff.c
 *
 *  Wed July 26 12:04:44 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
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

void test_ncm_diff_misc (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_property_minima (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_log_tables (TestNcmDiff *test, gconstpointer pdata);

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
void test_ncm_diff_rf_Hessian_N_to_1_rosenbrock (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_1_to_1_tiny_x (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_1_to_1_tiny_x_domain (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_1_to_1_func_abs_precision (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_1_to_1_extreme_x (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_rf_Hessian_N_to_1_tiny_x (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_half_line (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_log (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_interval (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_second_step (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_spectral_fallback (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_Hessian (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_narrow (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_Hessian_upper (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_dual (TestNcmDiff *test, gconstpointer pdata);
void test_ncm_diff_domain_narrow_subprocess (TestNcmDiff *test, gconstpointer pdata);

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

  g_test_add ("/ncm/diff/property_minima", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_property_minima,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/misc", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_misc,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/log_tables", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_log_tables,
              &test_ncm_diff_free);

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

  g_test_add ("/ncm/diff/rf/d1/N_to_1/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d1/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d1/N_to_1/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d2/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d2/N_to_1/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rf/Hessian/N_to_1/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_Hessian_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rf/Hessian/N_to_1/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rf_Hessian_N_to_1_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rf/Hessian/N_to_1/rosenbrock", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_Hessian_N_to_1_rosenbrock,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rf/d1/N_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_N_to_M_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rf/d1/N_to_M/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rf_d1_N_to_M_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d1/N_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_N_to_M_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d1/N_to_M/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d1_N_to_M_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d2/N_to_M/all", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_N_to_M_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rc/d2/N_to_M/zero", TestNcmDiff, GINT_TO_POINTER (TRUE),
              &test_ncm_diff_new,
              &test_ncm_diff_rc_d2_N_to_M_all,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/1_to_1/tiny_x", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_1_to_1_tiny_x,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/1_to_1/extreme_x", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_1_to_1_extreme_x,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/1_to_1/tiny_x/domain", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_1_to_1_tiny_x_domain,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/1_to_1/func_abs_precision", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_1_to_1_func_abs_precision,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/rf/Hessian/N_to_1/tiny_x", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_rf_Hessian_N_to_1_tiny_x,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/half_line", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_half_line,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/log", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_log,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/interval", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_interval,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/second_step", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_second_step,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/spectral_fallback", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_spectral_fallback,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/Hessian", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_Hessian,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/Hessian/upper", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_Hessian_upper,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/dual", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_dual,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/narrow", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_narrow,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/domain/narrow/subprocess", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_domain_narrow_subprocess,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/traps", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_traps,
              &test_ncm_diff_free);

  g_test_add ("/ncm/diff/invalid/st/subprocess", TestNcmDiff, NULL,
              &test_ncm_diff_new,
              &test_ncm_diff_invalid_st,
              &test_ncm_diff_free);

  g_test_run ();
}

void
test_ncm_diff_new (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = ncm_diff_new ();

  test->diff = diff;

  g_assert_true (test->diff != NULL);
  g_assert_true (NCM_IS_DIFF (test->diff));
}

void
test_ncm_diff_free (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;

  NCM_TEST_FREE (ncm_diff_free, diff);
}

/* Every property accepts the minimum its specification documents. */
void
test_ncm_diff_property_minima (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = g_object_new (NCM_TYPE_DIFF,
                                "max-order", 1,
                                "richardson-step", 1.1,
                                "round-off-pad", 1.01,
                                "terr-pad", 1.1,
                                "ini-h", GSL_DBL_EPSILON,
                                "func-abs-precision", 0.0,
                                NULL);

  g_assert_cmpuint (ncm_diff_get_max_order (diff), ==, 1);
  g_assert_cmpfloat (ncm_diff_get_richardson_step (diff), ==, 1.1);
  g_assert_cmpfloat (ncm_diff_get_round_off_pad (diff), ==, 1.01);
  g_assert_cmpfloat (ncm_diff_get_trunc_error_pad (diff), ==, 1.1);
  g_assert_cmpfloat (ncm_diff_get_ini_h (diff), ==, GSL_DBL_EPSILON);
  g_assert_cmpfloat (ncm_diff_get_func_abs_precision (diff), ==, 0.0);

  ncm_diff_free (diff);
}

void
test_ncm_diff_misc (TestNcmDiff *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmDiff *diff     = NULL;

  g_assert_true (test->diff != NULL);
  g_assert_true (NCM_IS_DIFF (test->diff));

  diff = NCM_DIFF (ncm_serialize_dup_obj (ser, G_OBJECT (test->diff)));
  ncm_diff_ref (diff);
  ncm_diff_free (diff);
  ncm_diff_clear (&diff);

  ncm_serialize_free (ser);

  {
    guint max_order = ncm_diff_get_max_order (test->diff);

    ncm_diff_set_max_order (test->diff, max_order);
    ncm_diff_set_max_order (test->diff, max_order + 2);
    ncm_diff_set_max_order (test->diff, max_order - 2);
  }

  {
    const gdouble rs = ncm_diff_get_richardson_step (test->diff);
    guint max_order  = ncm_diff_get_max_order (test->diff);

    ncm_diff_set_richardson_step (test->diff, rs * 2.0);
    ncm_diff_set_richardson_step (test->diff, rs * 1.1);
    ncm_diff_set_max_order (test->diff, max_order);
    ncm_diff_set_max_order (test->diff, max_order + 2);
    ncm_diff_set_max_order (test->diff, max_order - 2);
  }

  test_ncm_diff_rf_d1_1_to_1_sin (test, pdata);
}

void
test_ncm_diff_log_tables (TestNcmDiff *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    ncm_diff_log_central_tables (test->diff);
    ncm_diff_log_forward_tables (test->diff);
    ncm_diff_log_backward_tables (test->diff);

    return;
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_passed ();
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

/* sin (w x + 1): of order one at a tiny x, unlike sin (w x). */
static gdouble
_test_ncm_diff_sin1 (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;

  return sin (x * w_ptr[0] + 1.0);
}

static gdouble
_test_ncm_diff_dsin1 (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;

  return w_ptr[0] * cos (x * w_ptr[0] + 1.0);
}

static gdouble
_test_ncm_diff_d2sin1 (const gdouble x, gpointer userdata)
{
  gdouble *w_ptr = (gdouble *) userdata;

  return -gsl_pow_2 (w_ptr[0]) * sin (x * w_ptr[0] + 1.0);
}

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

static gdouble
_test_ncm_diff_rosenbrock (NcmVector *x, gpointer userdata)
{
  const gdouble x1 = ncm_vector_get (x, 0);
  const gdouble x2 = ncm_vector_get (x, 1);

  return 0.1 * (100.0 * gsl_pow_2 (x2 - x1 * x1) + gsl_pow_2 (1.0 - x1));
}

static gdouble
_test_ncm_diff_exp12 (NcmVector *x, gpointer userdata)
{
  return exp (ncm_vector_get (x, 0) + 2.0 * ncm_vector_get (x, 1));
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

    /* printf ("%d %d % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", i, nerr, err, x, w, df, Adf); */

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
  guint ntests  = 1000;
  gint nerr     = 5;
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
  guint ntests  = 1000;
  gint nerr     = 5;
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
  guint ntests  = 1000;
  gint nerr     = 5;
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  guint i, j;

  g_array_set_size (x_a, 3);

  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };

    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);

    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      GArray *df_a  = ncm_diff_rf_d1_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_dall (x_a, w);
      gdouble scale = 0.0;

      for (j = 0; j < x_a->len; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

      for (j = 0; j < x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);

        /*
         * The error estimate of the zero coordinate must stay informative; the
         * other coordinates have steps of ini_h times their size.
         */
        if (zero && (j == i % 3))
          g_assert_cmpfloat (err, <=, 0.1 * scale);
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  guint i, j;
  gint nerr = 5;

  g_array_set_size (x_a, 3);

  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };

    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);

    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      GArray *df_a  = ncm_diff_rc_d1_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_dall (x_a, w);
      gdouble scale = 0.0;

      for (j = 0; j < x_a->len; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

      for (j = 0; j < x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        if (!zero && ((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);

        /*
         * The error estimate of the zero coordinate must stay informative; the
         * other coordinates have steps of ini_h times their size.
         */
        if (zero && (j == i % 3))
          g_assert_cmpfloat (err, <=, 0.1 * scale);
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  gint nerr           = 5;
  guint i, j;

  g_array_set_size (x_a, 3);

  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };

    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);

    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      GArray *df_a  = ncm_diff_rc_d2_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_d2all (x_a, w);
      gdouble scale = 0.0;

      /*
       * At a zero coordinate the second derivatives vanish by symmetry; the
       * estimate is measured against the curvature amplitude, the derivatives
       * without their sine factor. A curvature below 1.0e-4 of f is within the
       * padded cancellation error of the initial step and cannot be known to 10 %.
       */
      {
        const gdouble x1  = g_array_index (x_a, gdouble, 0);
        const gdouble x2  = g_array_index (x_a, gdouble, 1);
        const gdouble x3  = g_array_index (x_a, gdouble, 2);
        const gdouble amp = GSL_MAX (gsl_pow_2 (x2 * w[0]), GSL_MAX (gsl_pow_2 (x1 * w[0]), gsl_pow_2 (w[1])));

        scale = GSL_MAX (amp, 1.0e-4) * exp (x3 * w[1]);
      }

      for (j = 0; j < x_a->len; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

      for (j = 0; j < x_a->len; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        if (!zero && ((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);

        /*
         * The error estimate of the zero coordinate must stay informative; the
         * other coordinates have steps of ini_h times their size.
         */
        if (zero && (j == i % 3))
          g_assert_cmpfloat (err, <=, 0.1 * scale);
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  guint nerr          = 0;
  guint i, j;

  g_array_set_size (x_a, 3);

  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };

    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);

    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      GArray *df_a  = ncm_diff_rf_Hessian_N_to_1 (diff, x_a, &_test_ncm_diff_N_to_1_all, w, &err_a);
      GArray *Adf_a = _test_ncm_diff_N_to_1_Hessian_all (x_a, w);
      gdouble scale = 0.0;

      for (j = 0; j < x_a->len * x_a->len; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

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

        /* The error estimate must stay informative at a zero coordinate. */
        if (zero)
          g_assert_cmpfloat (err, <=, 0.1 * scale);
      }

      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }

  g_array_unref (x_a);
}

void
test_ncm_diff_rf_Hessian_N_to_1_rosenbrock (TestNcmDiff *test, gconstpointer pdata)
{
  /*
   * The Hessian of 0.1 [100 (x2 - x1^2)^2 + (1 - x1)^2] is
   * [[120 x1^2 - 40 x2 + 0.2, -40 x1], [-40 x1, 20]]; at points with a zero
   * coordinate the curvature along x2 dominates the mixed term.
   */
  const gdouble pts[3][2] = {
    {
      0.0, 1.0
    }, {
      0.0, 0.0
    }, {
      1.0, 0.0
    }
  };
  NcmDiff *diff = test->diff;
  GArray *x_a   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint i, j;

  g_array_set_size (x_a, 2);

  for (i = 0; i < 3; i++)
  {
    const gdouble x1   = pts[i][0];
    const gdouble x2   = pts[i][1];
    const gdouble H[4] = {
      120.0 * x1 * x1 - 40.0 * x2 + 0.2, -40.0 * x1, -40.0 * x1, 20.0
    };
    GArray *err_a = NULL;
    GArray *df_a;
    gdouble scale = 0.0;

    g_array_index (x_a, gdouble, 0) = x1;
    g_array_index (x_a, gdouble, 1) = x2;

    df_a = ncm_diff_rf_Hessian_N_to_1 (diff, x_a, &_test_ncm_diff_rosenbrock, NULL, &err_a);

    for (j = 0; j < 4; j++)
      scale = GSL_MAX (scale, fabs (H[j]));

    for (j = 0; j < 4; j++)
    {
      const gdouble err = g_array_index (err_a, gdouble, j);

      ncm_assert_cmpdouble_e (g_array_index (df_a, gdouble, j), ==, H[j], 0.0, err);
      g_assert_cmpfloat (err, <=, 0.1 * scale);
    }

    g_array_unref (df_a);
    g_array_unref (err_a);
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  gint nerr           = 5;
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

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      const guint dim = 3;
      GArray *df_a    = ncm_diff_rf_d1_N_to_M (diff, x_a, dim, &_test_ncm_diff_N_to_M_all, w, &err_a);
      GArray *Adf_a   = _test_ncm_diff_N_to_M_dall (x_a, w);
      gdouble scale   = 0.0;

      for (j = 0; j < x_a->len * dim; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

      for (j = 0; j < x_a->len * dim; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        if (!zero && ((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);

        /*
         * The error estimate of the zero coordinate must stay informative; the
         * other coordinates have steps of ini_h times their size.
         */
        if (zero && (j / dim == i % 3))
          g_assert_cmpfloat (err, <=, 0.1 * scale);
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  gint nerr           = 5;
  guint i, j;

  g_array_set_size (x_a, 3);

  for (i = 0; i < ntests; i++)
  {
    gdouble w[3] =
    {
      g_test_rand_double_range (-10.0,         10.0),
      g_test_rand_double_range (-0.99,         0.99),
      g_test_rand_double_range (-0.5 * M_PI,   0.5 * M_PI)
    };

    const gdouble v1 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v2 = g_test_rand_double_range (-1.0, 1.0);
    const gdouble v3 = g_test_rand_double_range (-1.0, 1.0);

    g_array_index (x_a, gdouble, 0) = v1;
    g_array_index (x_a, gdouble, 1) = v2;
    g_array_index (x_a, gdouble, 2) = v3;

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      const guint dim = 3;
      GArray *df_a    = ncm_diff_rc_d1_N_to_M (diff, x_a, dim, &_test_ncm_diff_N_to_M_all, w, &err_a);
      GArray *Adf_a   = _test_ncm_diff_N_to_M_dall (x_a, w);
      gdouble scale   = 0.0;

      for (j = 0; j < x_a->len * dim; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

      for (j = 0; j < x_a->len * dim; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        if (!zero && ((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);

        /*
         * The error estimate of the zero coordinate must stay informative; the
         * other coordinates have steps of ini_h times their size.
         */
        if (zero && (j / dim == i % 3))
          g_assert_cmpfloat (err, <=, 0.1 * scale);
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
  /* With pdata set, coordinate i % 3 of the i-th point is zero. */
  const gboolean zero = GPOINTER_TO_INT (pdata);
  NcmDiff *diff       = test->diff;
  GArray *x_a         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *err_a       = NULL;
  guint ntests        = 1000;
  gint nerr           = 15;
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

    if (zero)
      g_array_index (x_a, gdouble, i % 3) = 0.0;

    {
      const guint dim = 3;
      GArray *df_a    = ncm_diff_rc_d2_N_to_M (diff, x_a, dim, &_test_ncm_diff_N_to_M_all, w, &err_a);
      GArray *Adf_a   = _test_ncm_diff_N_to_M_d2all (x_a, w);
      gdouble scale   = 0.0;

      /*
       * At a zero coordinate the second derivatives vanish by symmetry; the
       * estimate is measured against the curvature amplitude, the derivatives
       * without their sine and cosine factors. A curvature below 1.0e-4 of f
       * is within the padded cancellation error of the initial step and cannot be
       * known to 10 %.
       */
      {
        const gdouble x1  = g_array_index (x_a, gdouble, 0);
        const gdouble x2  = g_array_index (x_a, gdouble, 1);
        const gdouble x3  = g_array_index (x_a, gdouble, 2);
        const gdouble amp = GSL_MAX (gsl_pow_2 (x2 * w[0]), GSL_MAX (gsl_pow_2 (x1 * w[0]), gsl_pow_2 (w[1])));

        scale = GSL_MAX (amp, 1.0e-4) * exp (fabs (x3 * w[1]));
      }

      for (j = 0; j < x_a->len * dim; j++)
        scale = GSL_MAX (scale, fabs (g_array_index (Adf_a, gdouble, j)));

      for (j = 0; j < x_a->len * dim; j++)
      {
        const gdouble df  = g_array_index (df_a,  gdouble, j);
        const gdouble Adf = g_array_index (Adf_a, gdouble, j);
        const gdouble err = g_array_index (err_a, gdouble, j);

        /*printf ("(% 22.15g % 22.15g % 22.15g) % 22.15g % 22.15g % 22.15g\n", v1, v2, v3, df, Adf, err);*/
        if (!zero && ((err == 0.0) || gsl_isnan (err)) && nerr)
        {
          nerr--;
          g_test_skip ("Unable to estimate error.");
          continue;
        }

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);

        /*
         * The error estimate of the zero coordinate must stay informative; the
         * other coordinates have steps of ini_h times their size.
         */
        if (zero && (j / dim == i % 3))
          g_assert_cmpfloat (err, <=, 0.1 * scale);
      }

      g_array_unref (df_a);
      g_array_unref (Adf_a);
      g_array_unref (err_a);
    }
  }

  g_array_unref (x_a);
}

/*
 * A tiny nonzero x gives a step ini_h |x| far below the scale of a function
 * of order one there: the quotients are cancellation noise. The derivative
 * must still come out inside an informative estimate.
 */
void
test_ncm_diff_1_to_1_tiny_x (TestNcmDiff *test, gconstpointer pdata)
{
  typedef gdouble (*Method) (NcmDiff *, const gdouble, NcmDiffFunc1to1, gpointer, gdouble *);

  const Method methods[5] = {
    &ncm_diff_rf_d1_1_to_1, &ncm_diff_rc_d1_1_to_1, &ncm_diff_rc_d2_1_to_1,
    &ncm_diff_sc_d1_1_to_1, &ncm_diff_sc_d2_1_to_1
  };
  const gchar *names[5]             = {"rf_d1", "rc_d1", "rc_d2", "sc_d1", "sc_d2"};
  const NcmDiffFunc1to1 funcs[2]    = {&_test_ncm_diff_exp, &_test_ncm_diff_sin1};
  const NcmDiffFunc1to1 exact[2][5] = {
    {&_test_ncm_diff_dexp, &_test_ncm_diff_dexp, &_test_ncm_diff_d2exp, &_test_ncm_diff_dexp, &_test_ncm_diff_d2exp},
    {&_test_ncm_diff_dsin1, &_test_ncm_diff_dsin1, &_test_ncm_diff_d2sin1, &_test_ncm_diff_dsin1, &_test_ncm_diff_d2sin1}
  };
  const gdouble x0s[3] = {1.0e-4, 1.0e-8, 1.0e-12};
  NcmDiff *diff        = test->diff;
  gdouble w            = 1.0;
  guint m, k, n;

  for (m = 0; m < G_N_ELEMENTS (methods); m++)
  {
    for (k = 0; k < 2; k++)
    {
      for (n = 0; n < 3; n++)
      {
        gdouble err       = 0.0;
        const gdouble df  = methods[m](diff, x0s[n], funcs[k], &w, &err);
        const gdouble Adf = exact[k][m](x0s[n], &w);

        g_test_message ("%s f%u x = %.0e: % .15e exact % .15e err %.2e", names[m], k, x0s[n], df, Adf, err);
        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
        g_assert_cmpfloat (err, <=, 1.0e-3 * fabs (Adf));
      }
    }
  }
}

/* Counts the evaluations of the function it wraps. */
typedef struct _TestNcmDiffCount
{
  NcmDiffFunc1to1 f;
  gdouble w;
  guint n;
} TestNcmDiffCount;

static gdouble
_test_ncm_diff_counted (const gdouble x, gpointer userdata)
{
  TestNcmDiffCount *cnt = (TestNcmDiffCount *) userdata;

  cnt->n++;

  return cnt->f (x, &cnt->w);
}

/*
 * The coordinates closest to zero the step formula meets: the step moved up
 * from ini_h |x| must stay finite, and the spectral search must reach the
 * scale of f in a bounded number of evaluations.
 */
void
test_ncm_diff_1_to_1_extreme_x (TestNcmDiff *test, gconstpointer pdata)
{
  typedef gdouble (*Method) (NcmDiff *, const gdouble, NcmDiffFunc1to1, gpointer, gdouble *);

  const Method methods[5] = {
    &ncm_diff_rf_d1_1_to_1, &ncm_diff_rc_d1_1_to_1, &ncm_diff_rc_d2_1_to_1,
    &ncm_diff_sc_d1_1_to_1, &ncm_diff_sc_d2_1_to_1
  };
  const gchar *names[5]             = {"rf_d1", "rc_d1", "rc_d2", "sc_d1", "sc_d2"};
  const gboolean spectral[5]        = {FALSE, FALSE, FALSE, TRUE, TRUE};
  const NcmDiffFunc1to1 funcs[2]    = {&_test_ncm_diff_exp, &_test_ncm_diff_sin1};
  const NcmDiffFunc1to1 exact[2][5] = {
    {&_test_ncm_diff_dexp, &_test_ncm_diff_dexp, &_test_ncm_diff_d2exp, &_test_ncm_diff_dexp, &_test_ncm_diff_d2exp},
    {&_test_ncm_diff_dsin1, &_test_ncm_diff_dsin1, &_test_ncm_diff_d2sin1, &_test_ncm_diff_dsin1, &_test_ncm_diff_d2sin1}
  };
  const gdouble x0s[2] = {1.0e-160, 1.0e-300};
  NcmDiff *diff        = test->diff;
  guint m, k, n;

  for (m = 0; m < G_N_ELEMENTS (methods); m++)
  {
    for (k = 0; k < 2; k++)
    {
      for (n = 0; n < 2; n++)
      {
        TestNcmDiffCount cnt = {funcs[k], 1.0, 0};
        gdouble err          = 0.0;
        const gdouble df     = methods[m](diff, x0s[n], &_test_ncm_diff_counted, &cnt, &err);
        const gdouble Adf    = exact[k][m](x0s[n], &cnt.w);

        g_test_message ("%s f%u x = %.0e: % .15e exact % .15e err %.2e evals %u", names[m], k, x0s[n], df, Adf, err, cnt.n);
        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
        g_assert_cmpfloat (err, <=, 1.0e-3 * fabs (Adf));

        if (spectral[m])
          g_assert_cmpuint (cnt.n, <, 300);
      }
    }
  }
}

/*
 * Functions whose scale is x itself are fine at a tiny x and have no
 * values on the other side of zero; any larger step must not spoil them.
 */
void
test_ncm_diff_1_to_1_tiny_x_domain (TestNcmDiff *test, gconstpointer pdata)
{
  typedef gdouble (*Method) (NcmDiff *, const gdouble, NcmDiffFunc1to1, gpointer, gdouble *);

  const Method methods[3]           = {&ncm_diff_rf_d1_1_to_1, &ncm_diff_rc_d1_1_to_1, &ncm_diff_rc_d2_1_to_1};
  const NcmDiffFunc1to1 funcs[2]    = {&_test_ncm_diff_log, &_test_ncm_diff_plaw};
  const NcmDiffFunc1to1 exact[2][3] = {
    {&_test_ncm_diff_dlog, &_test_ncm_diff_dlog, &_test_ncm_diff_d2log},
    {&_test_ncm_diff_dplaw, &_test_ncm_diff_dplaw, &_test_ncm_diff_d2plaw}
  };
  const gdouble x0s[2] = {1.0e-6, 1.0e-10};
  gdouble w[2]         = {1.0, 0.5};
  NcmDiff *diff        = test->diff;
  guint m, k, n;

  for (m = 0; m < 3; m++)
  {
    for (k = 0; k < 2; k++)
    {
      for (n = 0; n < 2; n++)
      {
        gdouble err       = 0.0;
        const gdouble df  = methods[m](diff, x0s[n], funcs[k], &w[k], &err);
        const gdouble Adf = exact[k][m](x0s[n], &w[k]);

        ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
        g_assert_cmpfloat (err, <=, 1.0e-3 * fabs (Adf));
      }
    }
  }
}

/* ln (1 + x) evaluated as written: 1 + x rounds, an absolute error of eps / 2. */
static gdouble
_test_ncm_diff_log1p_naive (const gdouble x, gpointer userdata)
{
  return log (1.0 + x);
}

/*
 * ln (1 + x) for x from 1e-10 to 1: its values carry an absolute error of
 * eps / 2, far above eps |f| at small x. With that precision stated every
 * estimate covers the error.
 */
void
test_ncm_diff_1_to_1_func_abs_precision (TestNcmDiff *test, gconstpointer pdata)
{
  typedef gdouble (*Method) (NcmDiff *, const gdouble, NcmDiffFunc1to1, gpointer, gdouble *);

  const Method methods[3] = {&ncm_diff_rf_d1_1_to_1, &ncm_diff_rc_d1_1_to_1, &ncm_diff_rc_d2_1_to_1};
  const guint orders[3]   = {1, 1, 2};
  NcmDiff *diff           = test->diff;
  gdouble prec            = 0.0;
  guint ntests            = 1000;
  guint m, i;

  g_assert_cmpfloat (ncm_diff_get_func_abs_precision (diff), ==, 0.0);
  g_object_set (diff, "func-abs-precision", 0.5 * GSL_DBL_EPSILON, NULL);
  g_object_get (diff, "func-abs-precision", &prec, NULL);
  g_assert_cmpfloat (prec, ==, 0.5 * GSL_DBL_EPSILON);

  for (m = 0; m < G_N_ELEMENTS (methods); m++)
  {
    for (i = 0; i < ntests; i++)
    {
      const gdouble x   = exp (g_test_rand_double_range (log (1.0e-10), 0.0));
      const gdouble Adf = (orders[m] == 1) ? 1.0 / (1.0 + x) : -1.0 / gsl_pow_2 (1.0 + x);
      gdouble err       = 0.0;
      const gdouble df  = methods[m](diff, x, &_test_ncm_diff_log1p_naive, NULL, &err);

      ncm_assert_cmpdouble_e (df, ==, Adf, 0.0, err);
    }
  }

  ncm_diff_set_func_abs_precision (diff, 0.0);
}

/* exp (x1 + 2 x2) at tiny coordinates: H = exp (x1 + 2 x2) [[1, 2], [2, 4]]. */
void
test_ncm_diff_rf_Hessian_N_to_1_tiny_x (TestNcmDiff *test, gconstpointer pdata)
{
  const gdouble pts[3][2] = {
    {
      1.0e-8, 1.0e-12
    }, {
      1.0e-4, 1.0
    }, {
      1.0e-12, 0.0
    }
  };
  NcmDiff *diff = test->diff;
  GArray *x_a   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint i, j;

  g_array_set_size (x_a, 2);

  for (i = 0; i < 3; i++)
  {
    const gdouble e    = exp (pts[i][0] + 2.0 * pts[i][1]);
    const gdouble H[4] = {e, 2.0 * e, 2.0 * e, 4.0 * e};
    GArray *err_a      = NULL;
    GArray *df_a;

    g_array_index (x_a, gdouble, 0) = pts[i][0];
    g_array_index (x_a, gdouble, 1) = pts[i][1];

    df_a = ncm_diff_rf_Hessian_N_to_1 (diff, x_a, &_test_ncm_diff_exp12, NULL, &err_a);

    for (j = 0; j < 4; j++)
    {
      const gdouble err = g_array_index (err_a, gdouble, j);

      g_test_message ("point %u entry %u: % .15e exact % .15e err %.2e", i, j, g_array_index (df_a, gdouble, j), H[j], err);
      ncm_assert_cmpdouble_e (g_array_index (df_a, gdouble, j), ==, H[j], 0.0, err);
      g_assert_cmpfloat (err, <=, 1.0e-3 * H[j]);
    }

    g_array_unref (df_a);
    g_array_unref (err_a);
  }

  g_array_unref (x_a);
}

/* exp (x) recording the extreme points it was evaluated at. */
typedef struct _TestNcmDiffDom
{
  gdouble xmin;
  gdouble xmax;
} TestNcmDiffDom;

static gdouble
_test_ncm_diff_exp_dom (const gdouble x, gpointer userdata)
{
  TestNcmDiffDom *dom = (TestNcmDiffDom *) userdata;

  dom->xmin = GSL_MIN (dom->xmin, x);
  dom->xmax = GSL_MAX (dom->xmax, x);

  return exp (x);
}

static gdouble
_test_ncm_diff_exp12_dom (NcmVector *x, gpointer userdata)
{
  TestNcmDiffDom *dom = (TestNcmDiffDom *) userdata;
  const gdouble x1    = ncm_vector_get (x, 0);

  dom->xmin = GSL_MIN (dom->xmin, x1);
  dom->xmax = GSL_MAX (dom->xmax, x1);

  return exp (x1 + 2.0 * ncm_vector_get (x, 1));
}

/* Sets a one-dimensional domain [lb, ub] on the differentiator. */
static void
_test_ncm_diff_set_domain_1d (NcmDiff *diff, const gdouble lb, const gdouble ub)
{
  NcmVector *lb_v = ncm_vector_new (1);
  NcmVector *ub_v = ncm_vector_new (1);

  ncm_vector_set (lb_v, 0, lb);
  ncm_vector_set (ub_v, 0, ub);
  ncm_diff_set_domain (diff, lb_v, ub_v);

  ncm_vector_free (lb_v);
  ncm_vector_free (ub_v);
}

/*
 * exp on [0, inf) at a tiny x: the central step cannot grow without leaving
 * the domain, so the scheme falls back to a forward difference, with one
 * warning per call unless they are off; the result is informative.
 */
void
test_ncm_diff_domain_half_line (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff        = test->diff;
  const gdouble x0s[2] = {1.0e-8, 1.0e-12};
  TestNcmDiffDom dom   = {GSL_POSINF, GSL_NEGINF};
  gdouble w            = 1.0;
  guint pass, n;

  _test_ncm_diff_set_domain_1d (diff, 0.0, GSL_POSINF);

  for (pass = 0; pass < 2; pass++)
  {
    const gboolean warn = (pass == 0);

    ncm_diff_set_domain_warnings (diff, warn);

    for (n = 0; n < 2; n++)
    {
      gdouble err1 = 0.0, err2 = 0.0, df1, df2;

      if (warn)
        g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*");

      df1 = ncm_diff_rc_d1_1_to_1 (diff, x0s[n], &_test_ncm_diff_exp_dom, &dom, &err1);

      if (warn)
        g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*");

      df2 = ncm_diff_rc_d2_1_to_1 (diff, x0s[n], &_test_ncm_diff_exp_dom, &dom, &err2);

      g_test_assert_expected_messages ();
      g_assert_cmpfloat (dom.xmin, >=, 0.0);

      ncm_assert_cmpdouble_e (df1, ==, _test_ncm_diff_dexp (x0s[n], &w), 0.0, err1);
      g_assert_cmpfloat (err1, <=, 1.0e-3 * _test_ncm_diff_dexp (x0s[n], &w));
      ncm_assert_cmpdouble_e (df2, ==, _test_ncm_diff_d2exp (x0s[n], &w), 0.0, err2);
      g_assert_cmpfloat (err2, <=, 1.0e-3 * _test_ncm_diff_d2exp (x0s[n], &w));
    }
  }

  ncm_diff_clear_domain (diff);
}

/* log on [0, inf): its scale is x itself, central stays inside, no fallback. */
void
test_ncm_diff_domain_log (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff        = test->diff;
  const gdouble x0s[2] = {1.0e-6, 1.0e-3};
  gdouble w            = 1.0;
  guint n;

  _test_ncm_diff_set_domain_1d (diff, 0.0, GSL_POSINF);

  for (n = 0; n < 2; n++)
  {
    gdouble err1 = 0.0, err2 = 0.0;
    const gdouble df1 = ncm_diff_rc_d1_1_to_1 (diff, x0s[n], &_test_ncm_diff_log, &w, &err1);
    const gdouble df2 = ncm_diff_rc_d2_1_to_1 (diff, x0s[n], &_test_ncm_diff_log, &w, &err2);

    ncm_assert_cmpdouble_e (df1, ==, _test_ncm_diff_dlog (x0s[n], &w), 0.0, err1);
    g_assert_cmpfloat (err1, <=, 1.0e-3 * fabs (_test_ncm_diff_dlog (x0s[n], &w)));
    ncm_assert_cmpdouble_e (df2, ==, _test_ncm_diff_d2log (x0s[n], &w), 0.0, err2);
    g_assert_cmpfloat (err2, <=, 1.0e-3 * fabs (_test_ncm_diff_d2log (x0s[n], &w)));
  }

  ncm_diff_clear_domain (diff);
}

/*
 * exp on [1, 2]: near the lower edge the fallback is forward, near the upper
 * one backward, in the middle central stays central and nothing is warned.
 */
void
test_ncm_diff_domain_interval (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff        = test->diff;
  const gdouble x0s[3] = {1.0 + 1.0e-8, 2.0 - 1.0e-8, 1.5};
  TestNcmDiffDom dom   = {GSL_POSINF, GSL_NEGINF};
  gdouble w            = 1.0;
  guint n;

  _test_ncm_diff_set_domain_1d (diff, 1.0, 2.0);

  for (n = 0; n < 3; n++)
  {
    const gboolean edge = (n < 2);
    gdouble err1        = 0.0, err2 = 0.0, df1, df2;

    if (edge)
      g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, (n == 0) ? "*NcmDiff*domain*forward*" : "*NcmDiff*domain*backward*");

    df1 = ncm_diff_rc_d1_1_to_1 (diff, x0s[n], &_test_ncm_diff_exp_dom, &dom, &err1);

    if (edge)
      g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, (n == 0) ? "*NcmDiff*domain*forward*" : "*NcmDiff*domain*backward*");

    df2 = ncm_diff_rc_d2_1_to_1 (diff, x0s[n], &_test_ncm_diff_exp_dom, &dom, &err2);

    g_test_assert_expected_messages ();
    g_assert_cmpfloat (dom.xmin, >=, 1.0);
    g_assert_cmpfloat (dom.xmax, <=, 2.0);

    ncm_assert_cmpdouble_e (df1, ==, _test_ncm_diff_dexp (x0s[n], &w), 0.0, err1);
    g_assert_cmpfloat (err1, <=, 1.0e-3 * _test_ncm_diff_dexp (x0s[n], &w));
    ncm_assert_cmpdouble_e (df2, ==, _test_ncm_diff_d2exp (x0s[n], &w), 0.0, err2);
    g_assert_cmpfloat (err2, <=, 1.0e-3 * _test_ncm_diff_d2exp (x0s[n], &w));
  }

  ncm_diff_clear_domain (diff);
}

/* Both points of a one-sided second difference must stay in the domain. */
void
test_ncm_diff_domain_second_step (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  guint dual, side;

  ncm_diff_set_richardson_step (diff, 3.0);
  ncm_diff_set_domain_warnings (diff, FALSE);

  for (dual = 0; dual < 2; dual++)
  {
    ncm_diff_set_dual_series (diff, dual);

    for (side = 0; side < 2; side++)
    {
      const gdouble lo   = side ? -0.01 : 0.0;
      const gdouble hi   = side ? 0.0 : 0.01;
      TestNcmDiffDom dom = {GSL_POSINF, GSL_NEGINF};
      gdouble err, value;

      _test_ncm_diff_set_domain_1d (diff, lo, hi);
      value = ncm_diff_rc_d2_1_to_1 (diff, 0.0, &_test_ncm_diff_exp_dom, &dom, &err);
      g_assert_cmpfloat (dom.xmin, >=, lo);
      g_assert_cmpfloat (dom.xmax, <=, hi);
      ncm_assert_cmpdouble_e (value, ==, 1.0, 0.0, err);
    }
  }

  ncm_diff_clear_domain (diff);
}

/*
 * exp on [0, 1] at x = 1e-12: no symmetric window fits, the spectral call
 * falls back to Richardson, which falls back to a forward difference; both
 * are warned about, nothing with domain-warnings off.
 */
void
test_ncm_diff_domain_spectral_fallback (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff      = test->diff;
  const gdouble x0   = 1.0e-12;
  TestNcmDiffDom dom = {GSL_POSINF, GSL_NEGINF};
  guint pass;

  _test_ncm_diff_set_domain_1d (diff, 0.0, 1.0);

  for (pass = 0; pass < 2; pass++)
  {
    const gboolean warn = (pass == 0);
    gdouble err         = 0.0;
    gdouble df;

    ncm_diff_set_domain_warnings (diff, warn);

    if (warn)
    {
      g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*spectral*Richardson*");
      g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*");
    }

    df = ncm_diff_sc_d1_1_to_1 (diff, x0, &_test_ncm_diff_exp_dom, &dom, &err);
    g_test_assert_expected_messages ();

    g_assert_cmpfloat (dom.xmin, >=, 0.0);
    g_assert_cmpfloat (dom.xmax, <=, 1.0);
    ncm_assert_cmpdouble_e (df, ==, exp (x0), 0.0, err);
    g_assert_cmpfloat (err, <=, 1.0e-3);
  }

  ncm_diff_clear_domain (diff);
}

/*
 * exp (x1 + 2 x2) with x1 on [0, inf) at x1 = 1e-8: only the diagonal along
 * x1 is central and falls back, one warning; the cross term is forward.
 */
void
test_ncm_diff_domain_Hessian (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff      = test->diff;
  NcmVector *lb_v    = ncm_vector_new (2);
  NcmVector *ub_v    = ncm_vector_new (2);
  GArray *x_a        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble e    = exp (1.0e-8 + 2.0);
  const gdouble H[4] = {e, 2.0 * e, 2.0 * e, 4.0 * e};
  TestNcmDiffDom dom = {GSL_POSINF, GSL_NEGINF};
  GArray *err_a      = NULL;
  GArray *df_a;
  guint j;

  ncm_vector_set (lb_v, 0, 0.0);
  ncm_vector_set (lb_v, 1, GSL_NEGINF);
  ncm_vector_set (ub_v, 0, GSL_POSINF);
  ncm_vector_set (ub_v, 1, GSL_POSINF);
  ncm_diff_set_domain (diff, lb_v, ub_v);

  g_array_set_size (x_a, 2);
  g_array_index (x_a, gdouble, 0) = 1.0e-8;
  g_array_index (x_a, gdouble, 1) = 1.0;

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*");
  df_a = ncm_diff_rf_Hessian_N_to_1 (diff, x_a, &_test_ncm_diff_exp12_dom, &dom, &err_a);
  g_test_assert_expected_messages ();
  g_assert_cmpfloat (dom.xmin, >=, 0.0);

  for (j = 0; j < 4; j++)
  {
    const gdouble err = g_array_index (err_a, gdouble, j);

    ncm_assert_cmpdouble_e (g_array_index (df_a, gdouble, j), ==, H[j], 0.0, err);
    g_assert_cmpfloat (err, <=, 1.0e-3 * H[j]);
  }

  ncm_diff_clear_domain (diff);
  g_array_unref (df_a);
  g_array_unref (err_a);
  g_array_unref (x_a);
  ncm_vector_free (lb_v);
  ncm_vector_free (ub_v);
}

/*
 * exp (x1 + 2 x2) with x1 on (-inf, 0] at x1 = -1e-8: the diagonal along x1
 * falls back from central to backward and the cross term from forward to
 * backward, two warnings; no point has x1 > 0.
 */
void
test_ncm_diff_domain_Hessian_upper (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff      = test->diff;
  NcmVector *lb_v    = ncm_vector_new (2);
  NcmVector *ub_v    = ncm_vector_new (2);
  GArray *x_a        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble e    = exp (-1.0e-8 + 2.0);
  const gdouble H[4] = {e, 2.0 * e, 2.0 * e, 4.0 * e};
  TestNcmDiffDom dom = {GSL_POSINF, GSL_NEGINF};
  GArray *err_a      = NULL;
  GArray *df_a;
  guint j;

  ncm_vector_set (lb_v, 0, GSL_NEGINF);
  ncm_vector_set (lb_v, 1, GSL_NEGINF);
  ncm_vector_set (ub_v, 0, 0.0);
  ncm_vector_set (ub_v, 1, GSL_POSINF);
  ncm_diff_set_domain (diff, lb_v, ub_v);

  g_array_set_size (x_a, 2);
  g_array_index (x_a, gdouble, 0) = -1.0e-8;
  g_array_index (x_a, gdouble, 1) = 1.0;

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*central*backward*");
  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*forward*backward*");
  df_a = ncm_diff_rf_Hessian_N_to_1 (diff, x_a, &_test_ncm_diff_exp12_dom, &dom, &err_a);
  g_test_assert_expected_messages ();
  g_assert_cmpfloat (dom.xmax, <=, 0.0);

  for (j = 0; j < 4; j++)
  {
    const gdouble err = g_array_index (err_a, gdouble, j);

    ncm_assert_cmpdouble_e (g_array_index (df_a, gdouble, j), ==, H[j], 0.0, err);
    g_assert_cmpfloat (err, <=, 1.0e-3 * H[j]);
  }

  ncm_diff_clear_domain (diff);
  g_array_unref (df_a);
  g_array_unref (err_a);
  g_array_unref (x_a);
  ncm_vector_free (lb_v);
  ncm_vector_free (ub_v);
}

/* The half-line case with the dual-series scheme: same fallback, same bounds. */
void
test_ncm_diff_domain_dual (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff      = test->diff;
  TestNcmDiffDom dom = {GSL_POSINF, GSL_NEGINF};
  const gdouble x0   = 1.0e-8;
  gdouble w          = 1.0;
  gdouble err1       = 0.0, err2 = 0.0, df1, df2;

  ncm_diff_set_dual_series (diff, TRUE);
  _test_ncm_diff_set_domain_1d (diff, 0.0, GSL_POSINF);

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*");
  df1 = ncm_diff_rc_d1_1_to_1 (diff, x0, &_test_ncm_diff_exp_dom, &dom, &err1);
  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*NcmDiff*domain*");
  df2 = ncm_diff_rc_d2_1_to_1 (diff, x0, &_test_ncm_diff_exp_dom, &dom, &err2);
  g_test_assert_expected_messages ();

  g_assert_cmpfloat (dom.xmin, >=, 0.0);
  ncm_assert_cmpdouble_e (df1, ==, _test_ncm_diff_dexp (x0, &w), 0.0, err1);
  g_assert_cmpfloat (err1, <=, 1.0e-3 * _test_ncm_diff_dexp (x0, &w));
  ncm_assert_cmpdouble_e (df2, ==, _test_ncm_diff_d2exp (x0, &w), 0.0, err2);
  g_assert_cmpfloat (err2, <=, 1.0e-3 * _test_ncm_diff_d2exp (x0, &w));

  ncm_diff_set_dual_series (diff, FALSE);
  ncm_diff_clear_domain (diff);
}

/* A domain too narrow for any useful step is an error, not a warning. */
void
test_ncm_diff_domain_narrow (TestNcmDiff *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/diff/domain/narrow/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcmDiff*domain*");
}

void
test_ncm_diff_domain_narrow_subprocess (TestNcmDiff *test, gconstpointer pdata)
{
  NcmDiff *diff = test->diff;
  gdouble w     = 1.0;
  gdouble err   = 0.0;

  _test_ncm_diff_set_domain_1d (diff, 1.0, 1.0 + 1.0e-13);
  ncm_diff_rc_d1_1_to_1 (diff, 1.0 + 5.0e-14, &_test_ncm_diff_exp, &w, &err);
}

void
test_ncm_diff_traps (TestNcmDiff *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/diff/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_diff_invalid_st (TestNcmDiff *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

