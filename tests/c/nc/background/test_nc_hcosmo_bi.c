/***************************************************************************
 *            test_nc_hcosmo_bi.c
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_hcosmo_bi.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>

typedef struct _TestNcHCosmoBI
{
  NcHCosmoBI *cosmo;
} TestNcHCosmoBI;

typedef struct _TestNcHCosmoBICase
{
  gdouble Omega_m0;
  gdouble H[3];
} TestNcHCosmoBICase;

/* Two parameter sets of cluster_sampler_v2.py and one with B < 0 (B of the NcHCosmoBI description) */
static const TestNcHCosmoBICase cases[] = {
  {
    0.334, {
      0.85 * 73.04, 73.04, 1.15 * 73.04
    }
  },
  {
    0.334, {
      0.90 * 73.04, 0.95 * 73.04, 1.15 * 73.04
    }
  },
  {
    0.1,   {
      0.6 * 70.0, 70.0, 1.4 * 70.0
    }
  },
};

void test_nc_hcosmo_bi_new (TestNcHCosmoBI *test, gconstpointer pdata);
void test_nc_hcosmo_bi_free (TestNcHCosmoBI *test, gconstpointer pdata);

void test_nc_hcosmo_bi_today (TestNcHCosmoBI *test, gconstpointer pdata);
void test_nc_hcosmo_bi_lcdm (TestNcHCosmoBI *test, gconstpointer pdata);
void test_nc_hcosmo_bi_ode (TestNcHCosmoBI *test, gconstpointer pdata);
void test_nc_hcosmo_bi_valid (TestNcHCosmoBI *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/hcosmo_bi/today", TestNcHCosmoBI, NULL,
              &test_nc_hcosmo_bi_new,
              &test_nc_hcosmo_bi_today,
              &test_nc_hcosmo_bi_free);

  g_test_add ("/nc/hcosmo_bi/lcdm", TestNcHCosmoBI, NULL,
              &test_nc_hcosmo_bi_new,
              &test_nc_hcosmo_bi_lcdm,
              &test_nc_hcosmo_bi_free);

  g_test_add ("/nc/hcosmo_bi/ode", TestNcHCosmoBI, NULL,
              &test_nc_hcosmo_bi_new,
              &test_nc_hcosmo_bi_ode,
              &test_nc_hcosmo_bi_free);

  g_test_add ("/nc/hcosmo_bi/valid", TestNcHCosmoBI, NULL,
              &test_nc_hcosmo_bi_new,
              &test_nc_hcosmo_bi_valid,
              &test_nc_hcosmo_bi_free);

  g_test_run ();
}

static void
_test_nc_hcosmo_bi_set (NcHCosmoBI *cosmo, const TestNcHCosmoBICase *c)
{
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_H1, c->H[0]);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_H2, c->H[1]);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_H3, c->H[2]);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_OMEGA_M, c->Omega_m0);
}

/*
 * Reference for test_nc_hcosmo_bi_ode(), integrated with CVODE (Adams, reltol 1e-15).
 * Bianchi I field equations in t_lb = H_0 (t_0 - t), 8 pi G = 1, units of H_0:
 * y = (ln a_i, V H_i / H_0), d ln a_i / dt_lb = -y_{3+i} / V and
 * d y_{3+i} / dt_lb = -(3 Omega_m0 / 2 + 3 Omega_L0 V).
 */
static gint
_test_nc_hcosmo_bi_ode_f (sunrealtype t_lb, N_Vector y, N_Vector ydot, gpointer user_data)
{
  const gdouble *Omega = user_data;
  const gdouble V      = exp (NV_Ith_S (y, 0) + NV_Ith_S (y, 1) + NV_Ith_S (y, 2));
  const gdouble src    = -(1.5 * Omega[0] + 3.0 * Omega[1] * V);
  guint i;

  for (i = 0; i < 3; i++)
  {
    NV_Ith_S (ydot, i)     = -NV_Ith_S (y, 3 + i) / V;
    NV_Ith_S (ydot, 3 + i) = src;
  }

  return 0;
}

void
test_nc_hcosmo_bi_new (TestNcHCosmoBI *test, gconstpointer pdata)
{
  test->cosmo = nc_hcosmo_bi_new ();

  g_assert_true (NC_IS_HCOSMO_BI (test->cosmo));
}

void
test_nc_hcosmo_bi_free (TestNcHCosmoBI *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hcosmo_bi_free, test->cosmo);
}

void
test_nc_hcosmo_bi_today (TestNcHCosmoBI *test, gconstpointer pdata)
{
  guint k;

  for (k = 0; k < G_N_ELEMENTS (cases); k++)
  {
    const TestNcHCosmoBICase *c = &cases[k];
    gdouble ln_a[3], E[3], H0;

    _test_nc_hcosmo_bi_set (test->cosmo, c);
    H0 = nc_hcosmo_bi_H0 (test->cosmo);

    nc_hcosmo_bi_ln_a (test->cosmo, 0.0, &ln_a[0], &ln_a[1], &ln_a[2]);
    nc_hcosmo_bi_E (test->cosmo, 0.0, &E[0], &E[1], &E[2]);

    g_assert_cmpfloat (fabs (ln_a[0]), <, 1.0e-15);
    g_assert_cmpfloat (fabs (ln_a[1]), <, 1.0e-15);
    g_assert_cmpfloat (fabs (ln_a[2]), <, 1.0e-15);
    ncm_assert_cmpdouble_e (E[0], ==, c->H[0] / H0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (E[1], ==, c->H[1] / H0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (E[2], ==, c->H[2] / H0, 1.0e-15, 0.0);
  }
}

void
test_nc_hcosmo_bi_lcdm (TestNcHCosmoBI *test, gconstpointer pdata)
{
  const TestNcHCosmoBICase c = {
    0.3, {
      70.0, 70.0, 70.0
    }
  };
  const gdouble Omega_L0 = 0.7;
  const gdouble r        = sqrt (Omega_L0);
  const gdouble age      = asinh (sqrt (Omega_L0 / c.Omega_m0)) / (1.5 * r);
  guint j;

  _test_nc_hcosmo_bi_set (test->cosmo, &c);

  ncm_assert_cmpdouble_e (nc_hcosmo_bi_Omega_s0 (test->cosmo), ==, 0.0, 0.0, 0.0);
  ncm_assert_cmpdouble_e (nc_hcosmo_bi_age (test->cosmo), ==, age, 1.0e-15, 0.0);

  for (j = 0; j < 20; j++)
  {
    const gdouble t_lb   = 0.95 * age * j / 19.0;
    const gdouble a      = cbrt (c.Omega_m0 / Omega_L0) * pow (sinh (1.5 * r * (age - t_lb)), 2.0 / 3.0);
    const gdouble E_lcdm = sqrt (c.Omega_m0 / gsl_pow_3 (a) + Omega_L0);
    gdouble ln_a[3], E[3];

    nc_hcosmo_bi_ln_a (test->cosmo, t_lb, &ln_a[0], &ln_a[1], &ln_a[2]);
    nc_hcosmo_bi_E (test->cosmo, t_lb, &E[0], &E[1], &E[2]);

    ncm_assert_cmpdouble_e (exp (ln_a[0]), ==, a, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (exp (ln_a[1]), ==, a, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (exp (ln_a[2]), ==, a, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (E[0], ==, E_lcdm, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (E[1], ==, E_lcdm, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (E[2], ==, E_lcdm, 1.0e-14, 0.0);
  }
}

void
test_nc_hcosmo_bi_ode (TestNcHCosmoBI *test, gconstpointer pdata)
{
  SUNContext sunctx;
  guint k;

  g_assert_cmpint (SUNContext_Create (SUN_COMM_NULL, &sunctx), ==, 0);

  for (k = 0; k < G_N_ELEMENTS (cases); k++)
  {
    const TestNcHCosmoBICase *c = &cases[k];
    N_Vector y                  = N_VNew_Serial (6, sunctx);
    SUNMatrix A                 = SUNDenseMatrix (6, 6, sunctx);
    SUNLinearSolver LS          = SUNLinSol_Dense (y, A, sunctx);
    gpointer cvode              = CVodeCreate (CV_ADAMS, sunctx);
    gdouble Omega[2];
    gdouble age, H0;
    guint i, j;

    _test_nc_hcosmo_bi_set (test->cosmo, c);
    H0       = nc_hcosmo_bi_H0 (test->cosmo);
    age      = nc_hcosmo_bi_age (test->cosmo);
    Omega[0] = c->Omega_m0;
    Omega[1] = nc_hcosmo_bi_Omega_L0 (test->cosmo);

    for (i = 0; i < 3; i++)
    {
      NV_Ith_S (y, i)     = 0.0;
      NV_Ith_S (y, 3 + i) = c->H[i] / H0;
    }

    g_assert_cmpint (CVodeInit (cvode, &_test_nc_hcosmo_bi_ode_f, 0.0, y), ==, CV_SUCCESS);
    g_assert_cmpint (CVodeSStolerances (cvode, 1.0e-15, 1.0e-30), ==, CV_SUCCESS);
    g_assert_cmpint (CVodeSetUserData (cvode, Omega), ==, CV_SUCCESS);
    g_assert_cmpint (CVodeSetMaxNumSteps (cvode, 100000), ==, CV_SUCCESS);
    g_assert_cmpint (CVodeSetLinearSolver (cvode, LS, A), ==, CV_SUCCESS);

    for (j = 1; j < 20; j++)
    {
      const gdouble t_lb_j = 0.95 * age * j / 19.0;
      gdouble ln_a[3], E[3], V, t_lb;

      g_assert_cmpint (CVode (cvode, t_lb_j, y, &t_lb, CV_NORMAL), ==, CV_SUCCESS);

      nc_hcosmo_bi_ln_a (test->cosmo, t_lb_j, &ln_a[0], &ln_a[1], &ln_a[2]);
      nc_hcosmo_bi_E (test->cosmo, t_lb_j, &E[0], &E[1], &E[2]);
      V = exp (NV_Ith_S (y, 0) + NV_Ith_S (y, 1) + NV_Ith_S (y, 2));

      for (i = 0; i < 3; i++)
      {
        ncm_assert_cmpdouble_e (exp (ln_a[i]), ==, exp (NV_Ith_S (y, i)), 1.0e-11, 0.0);
        ncm_assert_cmpdouble_e (E[i], ==, NV_Ith_S (y, 3 + i) / V, 1.0e-10, 0.0);
      }
    }

    CVodeFree (&cvode);
    SUNLinSolFree (LS);
    SUNMatDestroy (A);
    N_VDestroy (y);
  }

  SUNContext_Free (&sunctx);
}

void
test_nc_hcosmo_bi_valid (TestNcHCosmoBI *test, gconstpointer pdata)
{
  const TestNcHCosmoBICase valid = {
    0.3, {
      50.0, 70.0, 90.0
    }
  };
  const TestNcHCosmoBICase invalid = {
    0.99, {
      50.0, 70.0, 90.0
    }
  };

  _test_nc_hcosmo_bi_set (test->cosmo, &valid);
  g_assert_true (ncm_model_params_valid (NCM_MODEL (test->cosmo)));

  _test_nc_hcosmo_bi_set (test->cosmo, &invalid);
  g_assert_cmpfloat (nc_hcosmo_bi_Omega_L0 (test->cosmo), <, 0.0);
  g_assert_false (ncm_model_params_valid (NCM_MODEL (test->cosmo)));
}

