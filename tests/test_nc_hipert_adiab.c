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
void test_nc_hipert_adiab_derivatives (TestNcHIPertAdiab *test, gconstpointer pdata);
static gdouble _test_nc_hipert_iadiab_eval_nu (const gdouble x, gpointer userdata);
static gdouble _test_nc_hipert_iadiab_eval_xi (const gdouble x, gpointer userdata);
static gdouble _test_nc_hipert_iadiab_eval_E2 (const gdouble x, gpointer userdata);

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
 
  g_test_add ("/nc/hipert_adiab/derivatives", TestNcHIPertAdiab, NULL,
              &test_nc_hipert_adiab_new,
              &test_nc_hipert_adiab_derivatives,
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
 gdouble k;
 NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
 gdouble tau = ncm_rng_uniform_gen(rng, -50.0, -10.0);
 gdouble *taua = malloc(sizeof *taua);
 gdouble *taui = malloc(sizeof *taui);
 gdouble *J11 = malloc(sizeof *J11);
 gdouble *J12 = malloc(sizeof *J12);
 gdouble *J22 = malloc(sizeof *J22);  
 gdouble *nut = malloc(sizeof *nut);
 gdouble *xit = malloc(sizeof *xit);
 gdouble *F1t = malloc(sizeof *F1t);

 k    = 1.0e9;
 ncm_csq1d_set_k (test->pa, k);
 ncm_csq1d_set_reltol (test->pa, 1.0e-5);
 ncm_csq1d_find_adiab_time_limit (test->pa, test->model, -1.0e2, 20.0, 1.0e-5, taui);
 ncm_csq1d_set_ti (test->pa, *taui);
 ncm_csq1d_set_tf (test->pa, 0.0);
 ncm_csq1d_set_init_cond_adiab (test->pa, test->model, *taui);
 ncm_csq1d_prepare(test->pa, test->model);
 GArray *t_array = ncm_csq1d_get_time_array (test->pa, taua);
 ncm_csq1d_get_J_at(test->pa, test->model, tau, J11, J12, J22);
 nc_hipert_iadiab_eval_system(NC_HIPERT_IADIAB(test->model), tau, k, nut, xit, F1t);
 
 {
  const gdouble nur = nc_hipert_iadiab_eval_nu(NC_HIPERT_IADIAB(test->model), tau, k); 
  const gdouble xir = nc_hipert_iadiab_eval_xi(NC_HIPERT_IADIAB(test->model), tau, k);
  const gdouble F1r = nc_hipert_iadiab_eval_F1(NC_HIPERT_IADIAB(test->model), tau, k);
  const gdouble m   = nc_hipert_iadiab_eval_m(NC_HIPERT_IADIAB(test->model), tau, k);
  
  ncm_assert_cmpdouble_e (*xit, ==, xir, 0.0, 0.0);
  ncm_assert_cmpdouble_e (*nut, ==, nur, 0.0, 0.0);
  ncm_assert_cmpdouble_e (*F1t, ==, F1r, 0.0, 0.0);
  ncm_assert_cmpdouble_e (xir, ==, log (m * nur), 0.0, 0.0);
 }
 g_clear_pointer (&t_array, g_array_unref);
 ncm_rng_free (rng);
}

void
test_nc_hipert_adiab_derivatives (TestNcHIPertAdiab *test, gconstpointer pdata)
{
 gdouble k, F1, F2, dE2dz, d2E2dz2;
 NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
 gdouble tau = ncm_rng_uniform_gen(rng, -50.0, -10.0);
 gdouble z   = ncm_rng_uniform_gen(rng, 11.0, 1.0);
 gdouble w   = g_test_rand_double_range (1.0e-3, 1.0e2);
 gdouble err = 0.0;
 
 k = 1.0e9;
 F1 = nc_hipert_iadiab_eval_F1(NC_HIPERT_IADIAB(test->model), tau, k);
 F2 = nc_hipert_iadiab_eval_F2(NC_HIPERT_IADIAB(test->model), tau, k);
 dE2dz = nc_hicosmo_dE2_dz (NC_HICOSMO(test->model), z);
 d2E2dz2 = nc_hicosmo_d2E2_dz2 (NC_HICOSMO(test->model), z);
 ncm_assert_cmpdouble_e (nc_hicosmo_E2 (NC_HICOSMO(test->model), z), ==, _test_nc_hipert_iadiab_eval_E2(z, &w), 0.0, 1.0e-3); 
 
 {
  NcmDiff *diff = ncm_diff_new();
  const gdouble dnu  = ncm_diff_rc_d1_1_to_1 (diff, tau, &_test_nc_hipert_iadiab_eval_nu, &w, &err);
  const gdouble dxi  = ncm_diff_rc_d1_1_to_1 (diff, tau, &_test_nc_hipert_iadiab_eval_xi, &w, &err);
  const gdouble d2xi  = ncm_diff_rc_d2_1_to_1 (diff, tau, &_test_nc_hipert_iadiab_eval_xi, &w, &err);
  const gdouble dEz  = ncm_diff_rc_d1_1_to_1 (diff, z, &_test_nc_hipert_iadiab_eval_E2, &w, &err);
  const gdouble d2Ez  = ncm_diff_rc_d2_1_to_1 (diff, z, &_test_nc_hipert_iadiab_eval_E2, &w, &err);
  const gdouble nu = nc_hipert_iadiab_eval_nu (NC_HIPERT_IADIAB(test->model), tau, k);
  const gdouble F1lTest = 1.0 * dxi/ (2.0 * nu);
  const gdouble F2Test = 1.0  / (4.0 * nu * nu) * (d2xi - dxi * dnu / nu);
  
  ncm_assert_cmpdouble_e ((F1lTest - F1) / F1, ==, 0.0, 0.0, 1.0e-3);
  ncm_assert_cmpdouble_e ((F2Test - F2) / F2, ==, 0.0, 0.0, 1.0e-3);
  ncm_assert_cmpdouble_e ((dEz - dE2dz) / dE2dz, ==, 0.0, 0.0, 1.0e-3);  
  ncm_assert_cmpdouble_e ((d2Ez - d2E2dz2) / d2E2dz2, ==, 0.0, 0.0, 1.0e-3);
  ncm_diff_free (diff);
 }
 
ncm_rng_free (rng);
}

static gdouble
_test_nc_hipert_iadiab_eval_xi (const gdouble x, gpointer userdata)
{
  const gdouble k = 1.0e9;
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
_test_nc_hipert_iadiab_eval_nu (const gdouble x, gpointer userdata)
{
  const gdouble k = 1.0e9;
  const gdouble X_B = 1.0e25;
  const gdouble w = 0.33;
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (x);
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = 1.0;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble _x = exp(lnX_B - tabs);
  return sqrt(w) * _x * k * 1.0 / E;
}

static gdouble
_test_nc_hipert_iadiab_eval_E2 (const gdouble x, gpointer userdata)
{
  const gdouble X_B = 1.0e25;
  const gdouble w = 0.33;
  const gdouble lnX_B = log (X_B);
  const gdouble _x = x + 1.0;
  const gdouble lnx = log(_x);
  const gdouble x_3_1pw = pow (_x, 3.0 * (1.0 + w));
  const gdouble Omega_w = 1.0;
  const gdouble e_3t_1mw = expm1 (3.0 *(1.0 - w) * (lnx - lnX_B));

  return -Omega_w *x_3_1pw * e_3t_1mw;
}