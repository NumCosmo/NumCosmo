/***************************************************************************
 *            recomb.c
 *
 *  Fri Oct 10 00:53:41 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_statistics_ulong.h>

static gint nsigma = -1;
static gint max_iter = 10000;
static gint snia_id = -1;
static gint bao_id = -1;
static gint cmb_id = -1;
static gint H_id = -1;
static gint cr_params[2] = {-1, -1};
static gint r_inj = 10;
static gchar *wmodel = "NcHICosmoDEXcdm";
static gint err_param = -1;
static gdouble alpha = 50.0f;
static gdouble sigma_alpha = 6.0f;
static gboolean calc_grid = FALSE;
static gboolean with_top = FALSE;
static gboolean with_H = FALSE;
static gboolean with_BBN = FALSE;
static gboolean test_TOP = FALSE;
static gboolean spherical = FALSE;
static gboolean resample = FALSE;
static gboolean change_params = FALSE;
static gboolean verbose = FALSE;
static gboolean flat = FALSE;

static GOptionEntry entries[] =
{
  { "n-sigma",       'n', 0, G_OPTION_ARG_INT,    &nsigma,        "Confidence region probability", NULL },
  { "cr-x",          'x', 0, G_OPTION_ARG_INT,    &cr_params[0],  "Confidence region x parameter", NULL },
  { "cr-y",          'y', 0, G_OPTION_ARG_INT,    &cr_params[1],  "Confidence region y parameter", NULL },
  { "max-iter",      'm', 0, G_OPTION_ARG_INT,    &max_iter,      "Max number of iterations used by the minimization algorithms", NULL },
  { "sn-id",         'S', 0, G_OPTION_ARG_INT,    &snia_id,       "ID of the SNe Ia sample to use", NULL },
  { "bao-id",        'B', 0, G_OPTION_ARG_INT,    &bao_id,        "ID of the BAO sample to use", NULL },
  { "cmb-id",        'C', 0, G_OPTION_ARG_INT,    &cmb_id,        "ID of the CMB sample to use", NULL },
  { "r-inj-n",       'r', 0, G_OPTION_ARG_INT,    &r_inj,         "r_inj n", NULL },
  { "model",         'w', 0, G_OPTION_ARG_STRING, &wmodel,        "Which model", NULL },
  { "alpha",         'a', 0, G_OPTION_ARG_DOUBLE, &alpha,         "alpha", NULL },
  { "sigma_alpha",   'A', 0, G_OPTION_ARG_DOUBLE, &sigma_alpha,   "Sigma_alpha", NULL },
  { "with-H",        'H', 0, G_OPTION_ARG_NONE,   &with_H,        "Use the HST H0 prior", NULL },
  { "H-id",          'E', 0, G_OPTION_ARG_INT,    &H_id,          "Use the H(z_i) data sample", NULL },
  { "with-BBN",      'N', 0, G_OPTION_ARG_NONE,   &with_BBN,      "Use BBN Prior", NULL },
  { "err-param",       0, 0, G_OPTION_ARG_INT,    &err_param,     "Calculate the error on param - i", "i" },
  { "test-TOP",        0, 0, G_OPTION_ARG_NONE,   &test_TOP,      "Test several topologies* and print the chi2", NULL },
  { "resample",      'R', 0, G_OPTION_ARG_NONE,   &resample,      "Resample model using default params", NULL },
  { "grid",          'G', 0, G_OPTION_ARG_NONE,   &calc_grid,     "Calculate grid", NULL },
  { "flat",          'p', 0, G_OPTION_ARG_NONE,   &flat,          "Plane model", NULL },
  { "with-TOP",      'T', 0, G_OPTION_ARG_NONE,   &with_top,      "Use Topology prior", NULL },
  { "spherical",       0, 0, G_OPTION_ARG_NONE,   &spherical,     "Impose Omega_k < 0", NULL },
  { "change-params", 'c', 0, G_OPTION_ARG_NONE,   &change_params, "Change the parametrization: \\Omega_L -> \\Omega_k", NULL },
  { "verbose",       'v', 0, G_OPTION_ARG_NONE,   &verbose,       "Be verbose", NULL },
  { NULL }
};

gint
main (gint argc, gchar *argv[])
{
  NcHICosmo *model;
  NcDataSet *ds;
  NcLikelihood *lh;
  NcmFit *fit = NULL;
  NcDistance *dist;
  NcmMSet *mset;
  GError *error = NULL;
  GOptionContext *context;
  gboolean tofit = TRUE;

  ncm_cfg_init ();

  context = g_option_context_new ("- test the recombination algorothims");
  g_option_context_add_main_entries (context, entries, NULL);
  g_option_context_parse (context, &argc, &argv, &error);

  dist = nc_distance_new (2.0);
  ds = nc_dataset_new ();
  model = nc_hicosmo_new_from_name (NC_TYPE_MODEL_DE, wmodel);
  mset = ncm_mset_new (NCM_MODEL (model), NULL);

  lh = nc_likelihood_new (ds);

  if (snia_id != -1)
  {
	NcData *snia = nc_data_distance_modulus_snia (dist, snia_id);
	nc_dataset_append_data (ds, snia);
  }

  if (cmb_id != -1)
  {
	NcData *cmb_data = nc_data_cmb (dist, cmb_id);
	nc_dataset_append_data (ds, cmb_data);
  }

  if (bao_id != -1)
  {
	NcData *bao_data = nc_data_bao (dist, bao_id);
	nc_dataset_append_data (ds, bao_data);
  }

  if (H_id != -1)
  {
	NcData *H_data = nc_data_hubble_function (H_id);
	nc_dataset_append_data (ds, H_data);
	ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_H0, NCM_PARAM_TYPE_FREE);
  }

  if (with_BBN)
	nc_hicosmo_de_new_add_bbn (lh);
  if (with_H)
  {
	nc_prior_add_gaussian_data (lh, NC_HICOSMO_ID, NC_HICOSMO_DE_H0, 72.0, 8.0);
	ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_H0, NCM_PARAM_TYPE_FREE);
  }

  if (snia_id == -1 && cmb_id == -1 && bao_id == -1 && H_id == -1)
  {
	printf ("# Will not fit.\n");
	tofit = FALSE;
  }

  if (spherical)
	nc_prior_add_oneside_a_inf_const_func (lh, nc_hicosmo_func0_new (&nc_hicosmo_Omega_k),
	                                       0.0, 1.0e-6);

  if (flat)
  {
	nc_hicosmo_de_omega_x2omega_k (model);
	ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FIXED);
  }
  else if (change_params)
	nc_hicosmo_de_omega_x2omega_k (model);

  if (tofit)
  {
	fit = ncm_fit_new (lh, mset, NCM_FIT_TYPE_LEAST_SQUARES, NCM_FIT_GRAD_ANALYTICAL);
	ncm_fit_run (fit, NC_BF_MAX_ITER, verbose);
	ncm_fit_log_info (fit);
	ncm_fit_numdiff_m2lnL_covar (fit);
	ncm_fit_log_covar (fit);
	ncm_mset_free (mset);
  }

  printf ("# Omega_r %.15f | %.15f | %.15f\n", nc_hicosmo_Omega_r (model), NC_C_RADIATION_TEMP_TO_h2OMEGA_R (2.725), NC_C_RADIATION_h2OMEGA_R_TO_TEMP(0.49 * nc_hicosmo_Omega_r (model)));
  printf ("#");
  ncm_model_params_print_all (NCM_MODEL (model), stdout);
  printf ("# z_r = %g\n", nc_distance_decoupling_redshift (dist, model));

  {
	gdouble x_step = 1.0;
	gdouble g_step = 0.0;
	FILE *nul = fopen ("/dev/null", "a");
	NcThermodynRecomb *recomb = nc_thermodyn_recomb_new (model, NC_THERMODYN_RECOMB_FULL);

	nul = stdout;
	nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
	printf ("# Ob %g h %g Obh2 %g => Or %g, zstar = %g\n",
	        nc_hicosmo_Omega_b (model),
	        nc_hicosmo_h (model),
	        nc_hicosmo_Omega_b (model) *
	        nc_hicosmo_h (model) *
	        nc_hicosmo_h (model),
	        nc_hicosmo_Omega_r (model),
	        recomb->x_rec - 1.0);
	nc_thermodyn_recomb_taubar (recomb, -log(NC_PERTURBATION_START_X));
	//    printf ("\n\n");
	//for (z_step = 3.0/T_PHOTON_0-1.0; z_step <= 1800.1; z_step += 1.0 / T_PHOTON_0)
	for (g_step = -log(NC_PERTURBATION_START_X); g_step < 0.0009 && FALSE; g_step += 1.0e-2)
	  //for (z_step = 1800.0; z_step >= -0.1 ; z_step -= 1.0e0)
	  //while (1)
	{
	  x_step = exp(-g_step);
	  fprintf (nul, "% 10.4f %.15g % .5g % .5g % .5g % 10.4f % 1.15e % 1.15g % 1.15g % 1.15g % 1.15g % 1.15g % 1.15g\n",
	           exp(-g_step),
	           nc_thermodyn_recomb_get_Xe_at (recomb, x_step),
	           nc_thermodyn_H_ionization_saha (model, x_step),
	           nc_thermodyn_HeI_ionization_saha (model, x_step),
	           nc_thermodyn_HeII_ionization_saha (model, x_step),
	           nc_hicosmo_T_gamma0 (model) * x_step,
	           nc_thermodyn_recomb_taubar (recomb, x_step) * sqrt(nc_hicosmo_E2 (model, x_step - 1.0)) / x_step,
	           nc_thermodyn_recomb_optical_depth (recomb, x_step),
	           nc_thermodyn_recomb_optical_depth_x0_x1 (recomb, x_step,   1.000001 * NC_PERTURBATION_START_X),
	           nc_thermodyn_recomb_optical_depth_R_x0_x1 (recomb, x_step, 1.000001 * NC_PERTURBATION_START_X),
	           nc_thermodyn_recomb_g (recomb, x_step),
	           nc_thermodyn_recomb_gbar (recomb, x_step),
	           nc_thermodyn_recomb_gbbar (recomb, x_step)
	           );
	  fflush (stdout);
	}

	if (TRUE)
	{
	  NcLinearPert *pert = nc_pert_linear_new (NCM_MODEL (model), 1 << 4, 1e-7, 1e-7, 1e-100, 1e-100);
	  NcLinearPertTF *pert_tf = nc_pert_transfer_function_new (pert, 1.0e-2, 1.0e3, 200.0);
	  gint i;
	  if (FALSE)
	  {
		nc_pert_transfer_function_prepare (pert_tf);
		for (i = 0; i < 1000; i++)
		{
		  gdouble logk = log (1.0e-2) + log (1.0e5) / (1000.0 - 1.0) * i;
		  gdouble k = exp (logk);
		  gdouble phi = nc_pert_transfer_function_get (pert_tf, k);
		  printf ("% 20.15e % 20.15e\n", k, phi);
		}
	  }
	  pert->pws->k = 1.0e0;
	  pert->solver->init (pert);
	  pert->solver->evol (pert, pert->g0);
	  pert->solver->init (pert);
	  pert->solver->evol (pert, pert->g0);
	  pert->solver->init (pert);
	  pert->solver->evol (pert, pert->g0);
	  pert->solver->init (pert);
	  pert->solver->evol (pert, pert->g0);
	  pert->solver->init (pert);
	  pert->solver->evol (pert, pert->g0);
	  pert->solver->init (pert);
	  pert->solver->evol (pert, pert->g0);
	}


	if (FALSE)
	{
	  NcLinearPert *pert = nc_pert_linear_new (NCM_MODEL (model), 1 << 3, 1e-7, 1e-7, 1e-10, 1e-10);
	  NcLinearPertSplines *pspline = nc_pert_linear_splines_new (pert, NC_LINEAR_PERTURBATIONS_SPLINE_ALL, 60, 220, 1.0e-2, 2000.0);
	  GArray *los_table = nc_pert_get_default_los_table (1000);
	  gint i = 0;
	  GTimer *total_time = g_timer_new ();
	  NcHICosmoLCDM *lcdm = nc_hicosmo_lcdm_new ();
	  NcHICosmoQGMode *qgmode;
	  gdouble mean_time = 0.0, elapsed, elapsed_total = 0.0;

	  nc_pert_linear_prepare_splines (pspline);
	  for (i = 0; i < 1000 && FALSE; i++)
	  {
		g_timer_start (total_time);
		pert->pws->k = 1e-1 * exp ( log (1e4) / 999.0 * i );
		pert->solver->init (pert);
		pert->solver->evol (pert, pert->g0);
		elapsed = g_timer_elapsed (total_time, NULL);
		elapsed_total += elapsed;
		mean_time = (i * mean_time + elapsed) / (i + 1.0);
		printf ("# elapsed %6.5f %6.5f %6.5f %6.5f\n", elapsed, mean_time, mean_time * 1000, elapsed_total);
	  }

	  printf("#[POS]");ncm_model_params_print_all (NCM_MODEL (model), stdout);
	  ncm_model_param_set (NCM_MODEL (lcdm), 1, 1.0 - 2.83178468907193e-07);
	  ncm_model_param_set (NCM_MODEL (lcdm), 2, log(1e-26));
	  ncm_model_param_set (NCM_MODEL (lcdm), 3, 2.83178468907193e-07);
	  ncm_model_param_set (NCM_MODEL (lcdm), 4, 1e-4);
	  printf("#[PRE]");ncm_model_params_print_all (NCM_MODEL (lcdm), stdout);

	  qgmode = nc_hicosmo_qg_modefunc (NCM_MODEL (lcdm), 1.0e6L, 1e-8L, 1e20L);
	  //nc_hicosmo_qg_modefunc_prepare_pw_spline (qgmode, TRUE);
	  //nc_pert_cov_direct (pert);
	  //nc_pert_linear_calculate_Ncs (pert);

	  for (i = 0; i < 100 && FALSE; i++)
	  {
		gint n;
		gdouble k = pspline->k0 + (pspline->k1 - pspline->k0) / 99.0 * i;
		for (n = 0; n < ncm_vector_len (pspline->ga); n++)
		{
		  gdouble eta = ncm_vector_get (pspline->ga, n);
		  gdouble xi = nc_scale_factor_z_t (pert->a, pert->eta0 - eta) + 1.0;
		  gdouble gi = log (xi);
		  gdouble taui = nc_thermodyn_recomb_optical_depth (pert->recomb, xi);
		  gdouble E2 = nc_hicosmo_E2 (model, xi - 1.0);
		  gdouble x_E = xi / sqrt (E2);
		  printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", k, eta, xi, gi, taui,
		          ncm_spline_eval (pspline->Sk[0][n], k),
		          ncm_spline_eval (pspline->Sk[1][n], k),
		          ncm_spline_eval (pspline->Sk[2][n], k),
		          x_E
		          );
		}
		printf ("\n\n");
	  }

	  if (FALSE)
	  {
		gint n;
		for (n = 0; n < ncm_vector_len (pspline->ga); n++)
		{
		  gdouble eta = ncm_vector_get (pspline->ga, n);
		  gdouble xi = nc_scale_factor_z_t (pert->a, pert->eta0 - eta) + 1.0;
		  gdouble gi = log (xi);
		  gdouble taui = nc_thermodyn_recomb_optical_depth (pert->recomb, xi);
		  gdouble E2 = nc_hicosmo_E2 (model, xi - 1.0);
		  gdouble x_E = xi / sqrt (E2);
		  for (i = 0; i < 300; i++)
		  {
			gdouble k = pspline->k0 + (pspline->k1 - pspline->k0) / 299.0 * i;
			printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", k, eta, xi, gi, taui,
			        ncm_spline_eval (pspline->Sk[0][n], k),
			        ncm_spline_eval (pspline->Sk[1][n], k),
			        ncm_spline_eval (pspline->Sk[2][n], k),
			        x_E
			        );
		  }
		  printf ("\n\n");
		}
	  }


	  if (FALSE)
	  {
		nc_pert_linear_calc_Nc_spline (pspline, qgmode->pw_spline, los_table, 500);
		for (i = 2; i <= 2000 && TRUE; i++)
		  printf ("%d %.15g %.15g\n", i, ncm_spline_eval (pspline->Nc, i), i * (i + 1.0) * ncm_spline_eval (pspline->Nc, i));
		printf ("\n\n");
	  }

	  if (TRUE)
	  {
		for (i = 0; i < 2000 && TRUE; i++)
		{
		  gdouble k = pspline->k0 * exp (log(pspline->k1 / pspline->k0 )/ (2000.0 - 1.0) * i);
		  gdouble phi = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_PHI), k);
		  gdouble theta0 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_THETA0), k);
		  gdouble c0 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_C0), k);
		  gdouble b0 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_B0), k);
		  gdouble theta1 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_THETA1), k);
		  gdouble theta2 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_THETA2), k);
		  gdouble c1 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_C1), k);
		  gdouble b1 = ncm_spline_eval (NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline, NC_PERT_B1), k);

		  printf ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
		          k, k / NC_C_HUBBLE_RADIUS, phi, theta0, c0, b0, theta1, c1, b1, theta2);
		}
	  }

	  return 0;
	  for (i = 0; i <= 100 && FALSE; i++)
		//      while (1)
	  {
		GTimer *bench = g_timer_new ();
		gdouble last_g;
		//gdouble timee;
		gint niter = 0;
		//pert->pws->k = 1.0e3 + (1000.0 - 0.1) / 100.0 * i;
		//pert->pws->k = 1.0e2; //0.10943257094418121 * NC_C / (100.0e3);
		//pert->solver->init (pert);

		g_timer_start (bench);
		//timee = g_timer_elapsed (bench, NULL);
		//        for (step = pert->a->ti*(1.0 + STEP); step <= pert->tf; step += (pert->tf-pert->a->ti) * STEP)
		//        for (step = pert->tf; step <= pert->tf; step += (pert->tf-pert->a->ti) * 1.0e-4)
		//        nc_pert_linear_evolve_to_end (pert);
		//nc_pert_linear_prepare_source_splines (pert, 0.1 + (1000.0 - 0.1) / 100.0 * i, 100);
		//nc_pert_linear_los_theta (pert, 2000);
		//pert->solver->evol (pert, pert->gf);
		while (FALSE && !pert->solver->evol_step (pert, pert->gf))
		{
		  gdouble z, x;
		  //gdouble timehere;
		  gdouble E;
		  gdouble taubar;
		  niter++;
		  //timee = g_timer_elapsed (bench, NULL);
		  //timehere = g_timer_elapsed (bench, NULL) - timee;
		  z = pert->solver->get_z (pert);
		  E = sqrt(nc_hicosmo_E2 (model, z));
		  x = 1.0 + z;
		  taubar = nc_thermodyn_recomb_taubar (recomb, x);

		  //          CVodeGetQuad (pert->cvode, &ttt, pert->yQ);
		  //          printf ("%.15g %.15g", pert->g, x);
		  //          for (ii = 50; ii <= 50; ii++)
		  //            printf (" %.15g", NV_Ith_S  (pert->yQ, ii));
		  //          printf ("\n");
		  //          printf ("%.15g\n", pert->g);
		  if (FALSE && (pert->pws->g > last_g * 0.9995))
		  {
			gint j;
			printf ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g",
			        pert->pws->g,                      //1
			        taubar,                       //2
			        pert->pws->k,                      //3
			        E,                            //4
			        x,                            //5
			        pert->solver->get_phi (pert), //6
			        pert->solver->get_c0 (pert),  //7
			        pert->solver->get_b0 (pert),  //8
			        pert->solver->get_c1 (pert),  //9
			        pert->solver->get_b1 (pert)   //10
			        );
			for (j = 0; j <= pert->lmax && j <= 5; j++)
			  printf (" %.15g", pert->solver->get_theta (pert, j));
			for (j = 0; j <= pert->lmax && j <= 5; j++)
			  printf (" %.15g", pert->solver->get_theta_p (pert, j));
			printf ("\n");
			last_g = pert->pws->g;
			fflush (stdout);
		  }
		  if (FALSE)
		  {
			gint j;
			for (j = 3; j <= pert->lmax; j++)
			  printf ("%.15g %d %.15g\n", pert->pws->g, j, pert->solver->get_theta (pert, j));
			printf ("\n");
		  }
		}
		if (FALSE)
		{
		  gint j;
		  printf ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g",
		          pert->pws->g,                             //1
		          pert->pws->k,                             //2
		          pert->solver->get_phi (pert),        //3
		          pert->solver->get_c0 (pert),         //4
		          pert->solver->get_b0 (pert),         //5
		          pert->solver->get_c1 (pert),         //6
		          pert->solver->get_b1 (pert)          //7
		          );
		  for (j = 0; j <= pert->lmax; j++)
			printf (" %.15g", pert->solver->get_theta (pert, j));
		  printf ("\n");
		}
		if (TRUE)
		{
		  printf ("# mode (%f) gi (%f) time %fs | total time %fs | niter %d\n", pert->pws->k, pert->pws->g, g_timer_elapsed (bench, NULL), g_timer_elapsed (total_time, NULL), niter);
		  pert->solver->print_stats (pert);
		  fflush (stdout);
		}
	  }
	  ncm_model_free (NCM_MODEL (lcdm));
	}
  }

  ncm_model_free (NCM_MODEL (model));
  nc_likelihood_free (lh);
  nc_dataset_free0 (ds, TRUE);

  return 0;
}

/* FIXME LIMBO: Code to calculate the set multiplication
 *
 g_timer_start (bench);
 if (TRUE)
 {
   glong index[300];
   glong last = 1000;
   glong lasta = 300;
   glong mina = 0;
   glong min = 1000000;
   glong ma[300];
   glong count = 0;
   glong tcount = 0;
   glong first = 0;
   for (i = 0; i < 300; i++)
   index[i] = 1;

   while (1)
   {
	 if (index[lasta-1] > last)
	 break;
	 for (i = first; i < lasta; i++)
	 {
	   glong val = (i + 1) * index[i];
	   if (index[i] > last) {first = i + 1; continue;}
	   if (val < min) { min = val; mina = 1; ma[0] = i;}
	   else if (val == min) { ma[mina] = i; mina++; }
	   if ((i != lasta - 1) && (index[i] == index[lasta - 1])) break;
	   }
	   for (i = 0; i < mina; i++)
	   index[ma[i]]++;
	   count++;
	   //printf ("# min here %ld count %ld | first %ld\n", min, mina, first);
	   min = 1000000;
	   tcount += mina;
	   }

	   printf ("# %ld*%ld = %ld | %ld %ld\n", last, lasta, last * lasta, count, tcount);
	   }
	   elapsed = g_timer_elapsed (bench, &times[0]);
	   printf ("# time: %lu | %.15g\n", times[0], elapsed);

	   */
