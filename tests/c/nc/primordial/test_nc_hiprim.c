/***************************************************************************
 *            test_nc_hiprim.c
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

/*
 * #NcHIPrimPowerLaw is defined by two closed forms,
 *
 *   P_SA(k) = A_s (k/k*)^(n_s - 1),   P_T(k) = r A_s (k/k*)^(n_T - 1),
 *
 * so the spectra are checked against those rather than against themselves. The tensor
 * side had no test at all: it is reached in production only when CLASS is handed a
 * primordial model through nc_cbe_set_prim(), which is an acceptance-tier run, and the
 * evaluation itself needs none of that.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_LN10E10ASA 3.1
#define TEST_T_SA_RATIO 0.15
#define TEST_N_SA 0.96
#define TEST_N_T (-0.02)

typedef struct _TestNcHIPrim
{
  NcHIPrim *prim;
} TestNcHIPrim;

static void
test_nc_hiprim_new (TestNcHIPrim *test, gconstpointer pdata)
{
  NcHIPrimPowerLaw *pl = nc_hiprim_power_law_new ();
  NcHIPrim *prim       = NC_HIPRIM (pl);

  ncm_model_orig_param_set (NCM_MODEL (prim), NC_HIPRIM_POWER_LAW_LN10E10ASA, TEST_LN10E10ASA);
  ncm_model_orig_param_set (NCM_MODEL (prim), NC_HIPRIM_POWER_LAW_T_SA_RATIO, TEST_T_SA_RATIO);
  ncm_model_orig_param_set (NCM_MODEL (prim), NC_HIPRIM_POWER_LAW_N_SA,       TEST_N_SA);
  ncm_model_orig_param_set (NCM_MODEL (prim), NC_HIPRIM_POWER_LAW_N_T,        TEST_N_T);

  test->prim = prim;

  g_assert_true (NC_IS_HIPRIM (prim));
  g_assert_true (NC_IS_HIPRIM_POWER_LAW (prim));
}

static void
test_nc_hiprim_free (TestNcHIPrim *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hiprim_free, test->prim);
}

/* Both spectra against their closed forms, over a range of k around the pivot. */
static void
test_nc_hiprim_spectra (TestNcHIPrim *test, gconstpointer pdata)
{
  const gdouble lnk_pivot = nc_hiprim_get_lnk_pivot (test->prim);
  guint i;

  ncm_assert_cmpdouble_e (nc_hiprim_get_k_pivot (test->prim), ==, exp (lnk_pivot), 1.0e-14, 0.0);

  for (i = 0; i <= 16; i++)
  {
    const gdouble lnk   = lnk_pivot - 4.0 + 0.5 * i;
    const gdouble ln_ka = lnk - lnk_pivot;
    const gdouble lnSA  = (TEST_N_SA - 1.0) * ln_ka + TEST_LN10E10ASA - 10.0 * M_LN10;
    const gdouble lnT   = TEST_N_T * ln_ka + TEST_LN10E10ASA - 10.0 * M_LN10 + log (TEST_T_SA_RATIO);

    ncm_assert_cmpdouble_e (nc_hiprim_lnSA_powspec_lnk (test->prim, lnk), ==, lnSA, 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (nc_hiprim_lnT_powspec_lnk (test->prim, lnk), ==, lnT, 1.0e-13, 0.0);

    /* The linear forms are the exponentials of the logarithmic ones. */
    ncm_assert_cmpdouble_e (nc_hiprim_SA_powspec_k (test->prim, exp (lnk)), ==, exp (lnSA), 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (nc_hiprim_T_powspec_k (test->prim, exp (lnk)), ==, exp (lnT), 1.0e-12, 0.0);
  }
}

/* The amplitudes are the spectra at the pivot, and their ratio is the parameter that
 * sets it -- which is the definition of r. */
static void
test_nc_hiprim_amplitudes (TestNcHIPrim *test, gconstpointer pdata)
{
  const gdouble k_pivot = nc_hiprim_get_k_pivot (test->prim);
  const gdouble SA      = nc_hiprim_SA_Ampl (test->prim);
  const gdouble T       = nc_hiprim_T_Ampl (test->prim);

  ncm_assert_cmpdouble_e (SA, ==, nc_hiprim_SA_powspec_k (test->prim, k_pivot), 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (T, ==, nc_hiprim_T_powspec_k (test->prim, k_pivot), 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (nc_hiprim_T_SA_ratio (test->prim), ==, TEST_T_SA_RATIO, 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (T / SA, ==, TEST_T_SA_RATIO, 1.0e-12, 0.0);
}

/* CLASS asks the model whether it implements the tensor spectrum before installing the
 * callback that reads it; a model that answered wrongly would be accepted and then
 * evaluated through a vfunc it does not have. */
static void
test_nc_hiprim_impl (TestNcHIPrim *test, gconstpointer pdata)
{
  g_assert_true (ncm_model_check_impl_opt (NCM_MODEL (test->prim), NC_HIPRIM_IMPL_lnSA_powspec_lnk));
  g_assert_true (ncm_model_check_impl_opt (NCM_MODEL (test->prim), NC_HIPRIM_IMPL_lnT_powspec_lnk));
}

static void
test_nc_hiprim_serialize (TestNcHIPrim *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcHIPrim *dup     = NC_HIPRIM (ncm_serialize_dup_obj (ser, G_OBJECT (test->prim)));
  const gdouble lnk = nc_hiprim_get_lnk_pivot (test->prim) + 1.0;

  g_assert_true (G_OBJECT_TYPE (dup) == G_OBJECT_TYPE (test->prim));

  ncm_assert_cmpdouble_e (nc_hiprim_lnSA_powspec_lnk (dup, lnk), ==,
                          nc_hiprim_lnSA_powspec_lnk (test->prim, lnk), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_hiprim_lnT_powspec_lnk (dup, lnk), ==,
                          nc_hiprim_lnT_powspec_lnk (test->prim, lnk), 1.0e-15, 0.0);

  nc_hiprim_free (dup);
  ncm_serialize_free (ser);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/hiprim/power_law/spectra", TestNcHIPrim, NULL,
              &test_nc_hiprim_new, &test_nc_hiprim_spectra, &test_nc_hiprim_free);
  g_test_add ("/nc/hiprim/power_law/amplitudes", TestNcHIPrim, NULL,
              &test_nc_hiprim_new, &test_nc_hiprim_amplitudes, &test_nc_hiprim_free);
  g_test_add ("/nc/hiprim/power_law/impl", TestNcHIPrim, NULL,
              &test_nc_hiprim_new, &test_nc_hiprim_impl, &test_nc_hiprim_free);
  g_test_add ("/nc/hiprim/power_law/serialize", TestNcHIPrim, NULL,
              &test_nc_hiprim_new, &test_nc_hiprim_serialize, &test_nc_hiprim_free);

  g_test_run ();
}

