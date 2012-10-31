/***************************************************************************
 *            quantum-gravity.c
 *
 *  Wed Jan 28 11:55:21 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_hyperg.h>

static gint nsigma = -1;
static gint max_iter = 10000;
static gint snia_id = -1;
static gint bao_id = -1;
static gint cmb_id = -1;
static gint H_id = -1;
static gint cr_params[2] = {-1, -1};
static gint err_param = -1;
static gboolean fit_model = FALSE;
static gboolean calc_grid = FALSE;
static gboolean with_H = FALSE;
static gboolean with_BBN = FALSE;
static gboolean resample = FALSE;
static gboolean verbose = FALSE;

static GOptionEntry entries[] =
{
  { "fit",           'f', 0, G_OPTION_ARG_NONE,   &fit_model,     "Fit model", NULL },
  { "n-sigma",       'n', 0, G_OPTION_ARG_INT,    &nsigma,        "Confidence region probability", NULL },
  { "cr-x",          'x', 0, G_OPTION_ARG_INT,    &cr_params[0],  "Confidence region x parameter", NULL },
  { "cr-y",          'y', 0, G_OPTION_ARG_INT,    &cr_params[1],  "Confidence region y parameter", NULL },
  { "max-iter",      'm', 0, G_OPTION_ARG_INT,    &max_iter,      "Max number of iterations used by the minimization algorithms", NULL },
  { "sn-id",         'S', 0, G_OPTION_ARG_INT,    &snia_id,       "ID of the SNe Ia sample to use", NULL },
  { "bao-id",        'B', 0, G_OPTION_ARG_INT,    &bao_id,        "ID of the BAO sample to use", NULL },
  { "cmb-id",        'C', 0, G_OPTION_ARG_INT,    &cmb_id,        "ID of the CMB sample to use", NULL },
  { "H-id",          'E', 0, G_OPTION_ARG_INT,    &H_id,          "Use the H(z_i) data sample", NULL },
  { "with-H",        'H', 0, G_OPTION_ARG_NONE,   &with_H,        "Use the HST H0", NULL },
  { "with-BBN",      'N', 0, G_OPTION_ARG_NONE,   &with_BBN,      "Use BBN", NULL },
  { "err-param",       0, 0, G_OPTION_ARG_INT,    &err_param,     "Calculate the error on param - i", "i" },
  { "resample",      'R', 0, G_OPTION_ARG_NONE,   &resample,      "Resample model using default params", NULL },
  { "grid",          'G', 0, G_OPTION_ARG_NONE,   &calc_grid,     "Calculate grid", NULL },
  { "verbose",       'v', 0, G_OPTION_ARG_NONE,   &verbose,       "Be verbose", NULL },
  { NULL }
};


gdouble
teste_func_filon (gdouble x, gpointer data)
{
  return gsl_pow_2(cos(5.0 * x)) * exp(M_PI * x);
}

gint
main(gint argc, gchar *argv[])
{
  NcHICosmo *model;
  GError *error = NULL;
  GOptionContext *context;

  ncm_cfg_init ();

  context = g_option_context_new ("- test the quantum gravity model");
  g_option_context_add_main_entries (context, entries, NULL);
  g_option_context_parse (context, &argc, &argv, &error);

  model = nc_hicosmo_qg_new ();

  if (TRUE)
  {
	gint i;
	GTimer *bench = g_timer_new ();
	gdouble nu = M_PI;
	gdouble xtot = -1e3;
	gdouble elapsed1 = 0.0, elapsed2 = 0.0;
	gint n = 10000;
	gdouble rest1[n];
	gdouble rest2[n];
	g_timer_start (bench);
	for (i = 0; i < n; i++)
	{
	  gdouble x = xtot / (n + 1.0) * (i + 1);
	  gdouble res1 = ncm_sf_0F1 (nu, x);
	  rest1[i] = res1;
	  //gdouble res2 = gsl_sf_hyperg_0F1 (nu, x);
	  //printf ("% 20.15g % 20.15g % 20.15g | %e\n", x, res1, res2, fabs((res1 - res2)/res1));
	}
	elapsed1 = g_timer_elapsed (bench, NULL);
	g_timer_start (bench);
	for (i = 0; i < n; i++)
	{
	  gdouble x = xtot / (n + 1.0) * (i + 1);
	  gdouble res2 = gsl_sf_hyperg_0F1 (nu, x);
	  rest2[i] = res2;
	}
	elapsed2 = g_timer_elapsed (bench, NULL);
	for (i = 0; i < n; i++)
	{
	  gdouble x = xtot / (n + 1.0) * (i + 1);
	  printf ("%20.15g % 20.15g % 20.15g %e\n", x, rest1[i], rest2[i], fabs((rest2[i]-rest2[i])/rest1[i]));
	}
	printf ("MY => %f | GSL => %f\n", elapsed1/n, elapsed2/n);
	exit (0);
  }


  if (FALSE)
  {
	gsl_matrix *B = gsl_matrix_alloc (2, 2);
	gsl_matrix *exp_B = gsl_matrix_alloc (2, 2);

	gsl_matrix_set (B, 0, 0, -5.0);
	gsl_matrix_set (B, 0, 1, sqrt(2.0));
	gsl_matrix_set (B, 1, 0, -100.0);
	gsl_matrix_set (B, 1, 1, 2.0);

	ncm_matrix_exp_2x2 (B, exp_B);

	printf ("%.15g %.15g\n", gsl_matrix_get (exp_B, 0, 0), gsl_matrix_get (exp_B, 0, 1));
	printf ("%.15g %.15g\n", gsl_matrix_get (exp_B, 1, 0), gsl_matrix_get (exp_B, 1, 1));

	exit(0);
  }

  {
	gdouble max, trans;
	NcHICosmoLCDM *lcdm = nc_hicosmo_lcdm_new ();
	gint i;

	ncm_model_param_set (NCM_MODEL (lcdm), 0, nc_hicosmo_H0 (model));
	ncm_model_param_set (NCM_MODEL (lcdm), 1, ncm_model_param_get (NCM_MODEL (model), 1));
	ncm_model_param_set (NCM_MODEL (lcdm), 2, 1.0 - ncm_model_param_get (NCM_MODEL (model), 1) - nc_hicosmo_Omega_r (model));

	ncm_model_param_set (NCM_MODEL (model), 1, 0.25 - 1e-6);
	ncm_model_param_set (NCM_MODEL (model), 2, log (1e-30));
	ncm_model_param_set (NCM_MODEL (model), 3, 1e-6);
	ncm_model_param_set (NCM_MODEL (model), 4, 1e-18);
	printf("#"); ncm_model_params_print_all (NCM_MODEL (model), stdout);
	nc_hicosmo_qg_max_z (NCM_MODEL (model), &max, &trans);
	NcScaleFactor *a = nc_scale_factor_new (NC_TIME_CONFORMAL, 1e20);
	nc_scale_factor_prepare (a, NC_HICOSMO (lcdm));

	printf ("# %g %g %g\n", nc_hicosmo_Omega_r (model), nc_hicosmo_Omega_m (model), 1.0 - nc_hicosmo_Omega_r(model) - nc_hicosmo_Omega_m (model));
	nc_hicosmo_E2 (model, max);
	printf ("# max = %g, trans = %g, E2(max) = %g, E2(trans) = %g, ETA_B = %g, T_0 = %g | %"G_GSIZE_FORMAT" %"G_GSIZE_FORMAT"\n",
	        max, trans, nc_hicosmo_E2 (model, max), nc_hicosmo_E2 (model, trans),
	        nc_hicosmo_qg_get_eta_b (NCM_MODEL (model), NULL), nc_hicosmo_qg_get_eta_b (NCM_MODEL (model), NULL) / ncm_c_c (),
	        sizeof(long double), sizeof(double));

	if (FALSE)
	{
	  for (i = 0; i <= 100000; i++)
	  {
		gdouble alpha = -20.0 + 40.0 / (100000.0) * i;
		gdouble x = 1.0e26 * exp (-alpha * alpha / 2.0);
		gdouble alphap2 = nc_hicosmo_qg_alphaprime2 (NCM_MODEL (model), alpha, NULL);
		gdouble dalphap2dalpha = nc_hicosmo_qg_dalphaprime2_dalpha (NCM_MODEL (model), alpha, NULL);
		printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", alpha, x, alphap2, dalphap2dalpha, dalphap2dalpha / alphap2);
	  }
	}

	if (FALSE)
	{
	  gdouble xi = 1e-7;
	  NcHICosmoQGMode *mode1 = nc_hicosmo_qg_modefunc (NCM_MODEL (model), 1.0e-4L, xi, 1e20L);
	  NcHICosmoQGMode *mode2 = nc_hicosmo_qg_pert_new (NCM_MODEL (model), 1e-20, xi);

	  mode2->type = NC_HICOSMO_QG_PERT_CURVATURE;
	  nc_hicosmo_qg_modefunc_init (mode1);
	  nc_hicosmo_qg_pert_init (mode2, 1.0e-4L);
	  nc_hicosmo_qg_pert_R_to_h (mode2, xi, N_VGetArrayPointer (mode2->y));

	  printf ("% 20.15e % 20.15e % 20.15e % 20.15e | % 20.15e % 20.15e\n", NV_Ith_S(mode1->y, 0), NV_Ith_S(mode1->y, 1), NV_Ith_S(mode1->y, 2), NV_Ith_S(mode1->y, 3),
	          hypot (NV_Ith_S(mode1->y, 0), NV_Ith_S(mode1->y, 2)), hypot (NV_Ith_S(mode1->y, 1), NV_Ith_S(mode1->y, 3)));
	  printf ("% 20.15e % 20.15e % 20.15e % 20.15e | % 20.15e % 20.15e\n", NV_Ith_S(mode2->y, 0), NV_Ith_S(mode2->y, 1), NV_Ith_S(mode2->y, 2), NV_Ith_S(mode2->y, 3),
	          hypot (NV_Ith_S(mode2->y, 0), NV_Ith_S(mode2->y, 2)), hypot (NV_Ith_S(mode2->y, 1), NV_Ith_S(mode2->y, 3)));
	}


	if (TRUE)
	{
	  NcHICosmoQGMode *mode = nc_hicosmo_qg_modefunc (NCM_MODEL (model), 1.0e6L, 1e-14L, 1e32L);
	  GTimer *bench = g_timer_new ();
	  gdouble part;
	  printf ("# x_eq = %.15g\n", pow((1.0 - 7.83178468907193e-05) / 7.83178468907193e-05, 1.0 / (1.0 - 3.0 * 1e-4)));
	  //nc_hicosmo_qg_modefunc_prepare_pw_spline (mode, TRUE);

	  for (i = 0; FALSE && i < 1000; i++)
	  {
		long double x = powl (10.0L, -8.0L + 20.0L / 999.0L * i);
		long double Vm1_4, x2cs2_E2;
		nc_hicosmo_qg_evolfunc (NCM_MODEL (model), x, &Vm1_4, &x2cs2_E2);
		printf ("%.20Lg %.20Lg %.20Lg %.15g %.15g\n", x, Vm1_4 + 1.0L / 4.0L,
		        x2cs2_E2, nc_hicosmo_qg_V (NCM_MODEL (model), x, NULL), nc_hicosmo_qg_cs2 (NCM_MODEL (model), x, NULL));
	  }

	  for (i = 0; FALSE && i < 1000; i++)
	  {
		gdouble logk = -4.0L + 8.0L / 999.0L * i;
		printf ("%.15g %.15g\n", logk, ncm_spline_eval (mode->pw_spline, logk));
	  }

	  for (i = 0; TRUE && i <= 1000; i++)
	  {
		mode->k = powl (10.0L, -3.0 + 8.0L / 1000.0L * i);
		printf ("# Calculating mode: %.15Lg\n", mode->k);
		nc_hicosmo_qg_modefunc_init (mode);
		part = g_timer_elapsed (bench, NULL);
		nc_hicosmo_qg_modefunc_evolve (mode);
		printf ("# Time to evolve: %f, total time: %f\n", g_timer_elapsed (bench, NULL)-part, g_timer_elapsed (bench, NULL));
		fflush (stdout);
	  }
	}

	if (FALSE)
	{
	  gdouble lambda_f = nc_hicosmo_qg_get_lambda_f (NCM_MODEL (model), NULL);
	  gdouble lambda_i = nc_hicosmo_qg_get_lambda_i (NCM_MODEL (model), NULL);
	  gdouble log_lambda_f_lambda_i = log (lambda_f / lambda_i);
	  for (i = 0; i < 10000; i++)
	  {
		gdouble lambda = lambda_i * exp (log_lambda_f_lambda_i / (10000.0 - 1.0) * i);
		gdouble eta = nc_hicosmo_qg_eta_lambda (NCM_MODEL (model), lambda, FALSE);
		gdouble alpha = nc_hicosmo_qg_int_1_zeta2_lambda (NCM_MODEL (model), lambda, FALSE);
		gdouble x = nc_hicosmo_qg_x_lambda (NCM_MODEL (model), lambda, FALSE);

		printf ("%.15g %.15g %.15g %.15g\n", lambda, x, eta, alpha);
	  }
	}

	if (FALSE)
	{
	  NcHICosmoQGMode *qgmode = nc_hicosmo_qg_pert_new (NCM_MODEL (model), 1e-30, 1e-15);
	  for (i = 0; i < 10; i++)
	  {
		nc_hicosmo_qg_pert_init (qgmode, 1e-2 * exp (5.0 * M_LN10 / 99.0 * i));
		nc_hicosmo_qg_pert_evolve (qgmode);
		printf ("\n\n");
		fflush (stdout);
	  }
	}
  }

  ncm_model_free (NCM_MODEL (model));

  return 0;
}
