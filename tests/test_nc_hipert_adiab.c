/***************************************************************************
 *            test_hipert_adiab.c
 *
 *  Thu August 11 14:34:27 2022
 *  Copyright  2022  Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Eduardo José Barroso <eduardo.jsbarroso@uel.br>
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

typedef struct _TestNcHIPertAdiab
{
  NcmCSQ1D *pa;
  NcmModel *model;
  NcDataClusterNCount *ncdata;
  gdouble area;
  guint ntests;
} TestNcHIPertAdiab;

void test_nc_hipert_adiab_new (TestNcHIPertAdiab *test, gconstpointer pdata);
void test_nc_hipert_adiab_free (TestNcHIPertAdiab *test, gconstpointer pdata);
void test_nc_hipert_adiab_sanity (TestNcHIPertAdiab *test, gconstpointer pdata);
static gdouble _test_nc_hipert_iadiab_eval_N (const gdouble x, gpointer userdata);
static gdouble _test_nc_hipert_iadiab_eval_xi (const gdouble x, gpointer userdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/hipert_adiab/sanity", TestNcHIPertAdiab, NULL,
              &test_nc_hipert_adiab_new,
              &test_nc_hipert_adiab_sanity,
              &test_nc_hipert_adiab_free);
 
 
  g_test_run ();
}

void
test_nc_hipert_adiab_new (TestNcHIPertAdiab *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoSFB");
  test->pa 		      = NCM_CSQ1D(nc_hipert_adiab_new ());
  test->model 		      = NCM_MODEL (cosmo);
  g_assert_true (test->pa != NULL);
  g_assert_true (NCM_IS_CSQ1D (test->pa));
}

void
test_nc_hipert_adiab_free (TestNcHIPertAdiab *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_csq1d_free, test->pa);
  NCM_TEST_FREE (ncm_model_free, test->model);
}

void
test_nc_hipert_adiab_sanity (TestNcHIPertAdiab *test, gconstpointer pdata)
{
 NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
 gdouble tau = ncm_rng_uniform_gen(rng, -50.0, -10.0);
 gdouble k, F1, F2, testF1, testF2, N, Ntest;
 gboolean found;
 gint i;
 gdouble *tauf = malloc(sizeof *tauf);
 gdouble *taua = malloc(sizeof *taua);
 gdouble *taui = malloc(sizeof *taui);
 gdouble *J11 = malloc(sizeof *J11);
 gdouble *J12 = malloc(sizeof *J12);
 gdouble *J22 = malloc(sizeof *J22);  
 gdouble w    = g_test_rand_double_range (1.0e-3, 1.0e2);  
 gdouble err       = 0.0;
 gdouble *nut = malloc(sizeof *nut);
 gdouble *xit = malloc(sizeof *xit);
 gdouble *F1t = malloc(sizeof *F1t);
 
 k    = 1.0e5;
 ncm_csq1d_set_k (test->pa, k);
 ncm_csq1d_find_adiab_time_limit (test->pa, test->model, -100.0, 20.0, 1.0e-2, taui);
 ncm_csq1d_set_ti (test->pa, *taui);
 ncm_csq1d_set_init_cond_adiab (test->pa, test->model, *taui);
 ncm_csq1d_prepare(test->pa, test->model);
  
 GArray *t_array = ncm_csq1d_get_time_array (test->pa, taua);
 ncm_csq1d_get_J_at(test->pa, test->model, tau, J11, J12, J22);
 F1 = nc_hipert_iadiab_eval_F1(NC_HIPERT_IADIAB(test->model), tau, k);
 F2 = nc_hipert_iadiab_eval_F2(NC_HIPERT_IADIAB(test->model), tau, k);
 Ntest = _test_nc_hipert_iadiab_eval_N (tau, &w);
 N = 1.0 / nc_hipert_iadiab_eval_H (NC_HIPERT_IADIAB(test->model), tau, k); 
 nc_hipert_iadiab_eval_system(NC_HIPERT_IADIAB(test->model), tau, k, nut, xit, F1t);
 ncm_assert_cmpdouble_e (N, ==, Ntest, 0.0, err);
 {
  NcmDiff *diff = ncm_diff_new();
  const gdouble df  = ncm_diff_rc_d1_1_to_1 (diff, tau, &_test_nc_hipert_iadiab_eval_N, &w, &err);
  const gdouble dl  = ncm_diff_rc_d1_1_to_1 (diff, tau, &_test_nc_hipert_iadiab_eval_xi, &w, &err);
  const gdouble d2l  = ncm_diff_rc_d2_1_to_1 (diff, tau, &_test_nc_hipert_iadiab_eval_N, &w, &err);
  const gdouble nu = nc_hipert_iadiab_eval_nu (NC_HIPERT_IADIAB(test->model), tau, k);
  const gdouble F1Test = 1.0 * df/ (2.0 * nu * N ); 
  const gdouble F1lTest = 1.0 * dl/ (2.0 * nu);
  const gdouble F2Test = 1.0  / (4.0 * nu * nu * N) * (d2l - 2.0 * pow(dl, 2) / N - dl);
  /*ncm_assert_cmpdouble_e (d2f, ==, F2, 0.0, 1.0e-7);*/
  /*ncm_assert_cmpdouble_e (F1Test, ==, -F1, 0.0, 1.0e0);
  ncm_assert_cmpdouble_e (F1lTest, ==, -F1, 0.0, 1.0e0);*/
  ncm_assert_cmpdouble_e (F2Test, ==, F2, 0.0, 1.0e-2);
 }

  

 printf("% 22.15g", N);
 g_clear_pointer (&t_array, g_array_unref);
 ncm_rng_free (rng);
}

static gdouble
_test_nc_hipert_iadiab_eval_xi (const gdouble x, gpointer userdata)
{
  const gdouble k = 1.0e5;
  const gdouble X_B = 1.0e25;
  const gdouble w = 0.33;
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (x);
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = 1.0;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble _x = exp(lnX_B - tabs);
  const gdouble mult = 3.0 * (1.0 + w) / (2.0 * w);
  const gdouble x3   = _x * _x * _x;
  const gdouble z2 = mult * 1.0 / (x3 * N);
  const gdouble m = 2.0 * z2;
  const gdouble nu = sqrt (w) * k * N * _x;
  const gdouble mnu =  m * nu;
  
  return log (mnu);
}

static gdouble
_test_nc_hipert_iadiab_eval_N (const gdouble x, gpointer userdata)
{
  const gdouble k = 1.0e5;
  const gdouble X_B = 1.0e25;
  const gdouble w = 0.33;
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (x);
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = 1.0;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  
  return 1.0 / E;
}