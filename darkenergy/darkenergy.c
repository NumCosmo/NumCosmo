/***************************************************************************
 *            darkenergy.c
 *
 *  Wed Apr 21 11:49:05 2010 (pulled from cosmolib/tools)
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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

#include "de_options.h"
#include "data_cluster.h"
#include "savedata.h"

#include <math.h>
#include <glib.h>
#include <gsl/gsl_blas.h>

gint
main (gint argc, gchar *argv[])
{
  NcHICosmo *model;
  NcDataSet *ds;
  NcLikelihood *lh;
  NcmFit *fit;
  NcDEModelEntries de_model = NC_DE_MODEL_ENTRIES;
  NcDEDataSimpleEntries de_data_simple = NC_DE_DATA_SIMPLE_ENTRIES;
  NcDEDataClusterEntries de_data_cluster = NC_DE_DATA_CLUSTER_ENTRIES;
  NcDEFitEntries de_fit = NC_DE_FIT_ENTRIES;
  NcDistance *dist;
  NcmMSet *mset, *fiduc;
  GError *error = NULL;
  GOptionContext *context;
  GPtrArray *ca_array = NULL;
  gchar *full_cmd_line = ncm_cfg_command_line (argv, argc);

  ncm_cfg_init ();

  context = g_option_context_new ("- test the dark energy models");
  g_option_context_set_summary (context, "DE Summary <FIXME>");
  g_option_context_set_description (context, "DE Description <FIXME>");

  g_option_context_set_main_group (context, nc_de_opt_get_model_group (&de_model));
  g_option_context_add_group (context, nc_de_opt_get_data_simple_group (&de_data_simple));
  g_option_context_add_group (context, nc_de_opt_get_data_cluster_group (&de_data_cluster));
  g_option_context_add_group (context, nc_de_opt_get_fit_group (&de_fit));

  if(!g_option_context_parse (context, &argc, &argv, &error))
  {
	g_warning ("Invalid options. (%s)", full_cmd_line);
	printf (g_option_context_get_help (context, TRUE, NULL));
	g_option_context_free (context);
	return 0;
  }

  if (de_fit.file_out != NULL)
	ncm_cfg_set_logfile (de_fit.file_out);

  if (de_data_simple.snia_id == -1 &&
      de_data_simple.cmb_id == -1 &&
      de_data_simple.bao_id == -1 &&
      de_data_simple.H_id == -1 &&
      de_data_simple.cluster_id == -1)
  {
	g_warning ("No dataset was chosen.");
	printf (g_option_context_get_help (context, TRUE, NULL));
	return 0;
  }

  g_option_context_free (context);

  nc_message ("# NumCosmo Version -- "NUMCOSMO_VERSION"\n");
  nc_message ("# Command Line: %s\n", full_cmd_line);

  ds = nc_dataset_new ();
  model = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, de_model.model_name);
  mset = ncm_mset_new (NCM_MODEL (model), NULL);
  dist = nc_distance_new (2.0);

  if (de_model.help_names)
  {
	gint i;
	nc_message ("# Model name -- %s\n", ncm_model_name (NCM_MODEL (model)));
	for (i = 0; i < ncm_model_len (NCM_MODEL (model)); i++)
	  nc_message ("# Model parameter [%02d] = %s\n", i, ncm_model_param_name (NCM_MODEL (model), i));
	return 0;
  }

  lh = nc_likelihood_new (ds);

  {
	ncm_model_params_set_all (NCM_MODEL (model),
	                           de_model.H0, de_model.Omega_c, de_model.Omega_x,
	                           de_model.T_gamma, de_model.Omega_b, de_model.n_s,
	                           de_model.sigma_8, de_model.w[0], de_model.w[1],
	                           de_model.w[2]);
  }

  if (de_fit.fit_params != NULL)
  {
	gchar **cmds, **cmds_orig;
	gboolean have_free = TRUE;

	cmds = g_strsplit (de_fit.fit_params, ",", 0);
	cmds_orig = cmds;

	//ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FIXED);

	for (; cmds[0] != NULL && FALSE; cmds = &cmds[1])
	{
	  gchar **name_cmd = g_strsplit (cmds[0], "=", 2);
	  gint index;

	  if (name_cmd[1] == NULL)
		g_error ("Usage --fit-params H0=fit, with the equal symbol =");

	  index = ncm_model_param_index_from_name (NCM_MODEL (model), name_cmd[0]);
	  if (index < 0)
		g_error ("Parameter (%s) not found, use --help-names to obtain a list of parameters names", name_cmd[0]);

	  if ((strncasecmp (name_cmd[1], "fix", 3) == 0) || (strncasecmp (name_cmd[1], "0", 1) == 0))
		ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, index, NCM_PARAM_TYPE_FIXED);
	  else if ((strncasecmp (name_cmd[1], "fit", 3) == 0) || (strncasecmp (name_cmd[1], "1", 1) == 0))
	  {
		ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, index, NCM_PARAM_TYPE_FREE);
		have_free = TRUE;
	  }
	  else
		g_error ("Command (%s) not recognised, use fix or 0 to keep the parameter fixed and fit or 1 to fit this parameter", name_cmd[1]);

	  g_strfreev (name_cmd);
	}
	g_strfreev (cmds_orig);

	if (!have_free)
	  g_error ("No parameter was set free");
  }

  if (de_model.flat)
  {
	nc_hicosmo_de_omega_x2omega_k (model);
	ncm_model_param_set (NCM_MODEL (model), NC_HICOSMO_DE_OMEGA_X, 0.0);
	ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FIXED);
  }
  else if (de_model.Omega_k)
  {
	nc_hicosmo_de_omega_x2omega_k (model);
  }

  if (de_data_simple.snia_id != -1)
  {
	NcData *snia = nc_data_distance_modulus_snia (dist, de_data_simple.snia_id);
	nc_dataset_append_data (ds, snia);
  }

  if (de_data_simple.cmb_id != -1)
  {
	NcData *cmb_data = nc_data_cmb (dist, de_data_simple.cmb_id);
	nc_dataset_append_data (ds, cmb_data);
  }

  if (de_data_simple.bao_id != -1)
  {
	NcData *bao_data = nc_data_bao (dist, de_data_simple.bao_id);
	nc_dataset_append_data (ds, bao_data);
  }

  if (de_data_simple.H_id != -1)
  {
	NcData *H_data = nc_data_hubble_function (de_data_simple.H_id);
	nc_dataset_append_data (ds, H_data);
  }

  if (de_data_simple.cluster_id != -1)
	ca_array = nc_de_data_cluster_new (dist, mset, &de_data_cluster, ds, de_data_simple.cluster_id);

  if (de_data_simple.BBN)
	nc_hicosmo_de_new_add_bbn (lh);

  if (de_data_simple.H0_Hst)
  {
	if (ncm_mset_param_get_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_H0) == NCM_PARAM_TYPE_FIXED)
	  g_warning ("A prior on H0 was required but H0 was kept fixed. This prior will be useless.");

	nc_prior_add_gaussian_data (lh, NC_HICOSMO_ID, NC_HICOSMO_DE_H0, 73.8, 2.4);
  }

  fiduc = ncm_mset_copy_all (mset);
  fit = ncm_fit_new (lh, mset, de_fit.min_algo, de_fit.diff_algo);

  de_fit.fisher = (de_fit.fisher || (de_fit.nsigma_fisher != -1) || (de_fit.nsigma != -1) || (de_fit.onedim_cr != NULL));
  de_fit.fit = (de_fit.fit || de_fit.fisher);
  de_fit.save_best_fit = (de_fit.save_best_fit || de_fit.save_fisher || (de_fit.nsigma_fisher != -1) || (de_fit.nsigma != -1));

  if (de_fit.fit)
  {
	if (de_fit.fit_params == NULL)
	  g_error ("No parameter specified, use --fit-params, or --help-names to obtain a list of parameters avaliable.");
	fit->params_prec_target = 1e-5;
	ncm_fit_run (fit, de_fit.max_iter, de_fit.msg_level);
	ncm_fit_log_info (fit);
  }

  if (de_fit.save_best_fit)
  {
	FILE *f_bf;
	gchar *bfile = NULL;

	f_bf = nc_de_open_dataout_file (model, "best_fit", &bfile);
	ncm_mset_params_pretty_print (fit->mset, f_bf, full_cmd_line);
	fclose (f_bf);

	nc_message ("# Params file: %s \n", bfile);

	g_free (bfile);
  }

  if (de_fit.fisher)
  {
	ncm_fit_numdiff_m2lnL_covar (fit);
	ncm_fit_log_covar (fit);
	if (de_fit.save_fisher)
	{
	  FILE *f_MF;
	  gchar *mfile = NULL;

	  f_MF = nc_de_open_dataout_file (model, "MF", &mfile);
	  ncm_fit_fishermatrix_print (fit, f_MF, full_cmd_line);
	  fclose (f_MF);

	  nc_message ("# FM file: %s \n", mfile);

	  g_free (mfile);

	}
  }

  if (de_fit.montecarlo > 0)
  {
	NcmMSet *resample_mset;
	NcmMatrix *param_matrix;

	if (de_fit.fiducial)
	  resample_mset = fiduc;
	else
	  resample_mset = fit->mset;

	param_matrix = ncm_fit_montecarlo_matrix (fit, resample_mset, de_fit.max_iter, de_fit.mc_ni, de_fit.montecarlo, de_fit.msg_level);
	ncm_fit_montecarlo_matrix_mean_covar (fit, param_matrix);
	ncm_fit_log_covar (fit);
	if (de_fit.mc_data)
	  ncm_fit_montecarlo_matrix_print (fit, param_matrix);

	ncm_matrix_free (param_matrix);
  }

  if (de_fit.onedim_cr != NULL)
  {
	while (de_fit.onedim_cr[0] != NULL)
	{
	  gint p_n = atoi (de_fit.onedim_cr[0]);
	  gdouble prob_sigma, err_inf, err_sup;
	  de_fit.onedim_cr = &de_fit.onedim_cr[1];
	  switch (de_fit.nsigma)
	  {
		case 1:
		  prob_sigma = NC_C_STATS_1SIGMA;
		  break;
		case 2:
		  prob_sigma = NC_C_STATS_2SIGMA;
		  break;
		case 3:
		  prob_sigma = NC_C_STATS_3SIGMA;
		  break;
		default:
		  prob_sigma = NC_C_STATS_1SIGMA;
		  break;
	  }
	  ncm_fit_cr_1dim (fit, NC_HICOSMO_ID, p_n, prob_sigma, 1, &err_inf, &err_sup);
	  nc_message ("#  One dimension confidence region for %s[%02d] = % .5g (% .5g, % .5g)\n",
	              ncm_model_param_name (NCM_MODEL (model), p_n), p_n,
	              ncm_model_param_get (NCM_MODEL (model), p_n), err_inf, err_sup);
	}
  }

  if (de_fit.resample)
  {
	nc_message ("# Resample from bestfit values\n");
	nc_dataset_resample (ds, mset, TRUE);
	ncm_fit_run (fit, NC_BF_MAX_ITER, de_fit.msg_level);
	ncm_fit_log_info (fit);
	ncm_fit_numdiff_m2lnL_covar (fit);
	ncm_fit_log_covar (fit);
  }

  if (de_fit.nsigma >= 0 && (de_fit.bidim_cr[0] != -1) && (de_fit.bidim_cr[1] != -1))
  {
	GList *points_1sigma;
	GList *points_2sigma;
	GList *points_3sigma;
	points_1sigma = NULL;
	points_2sigma = NULL;
	points_3sigma = NULL;
	switch (de_fit.nsigma)
	{
	  case 1:
		points_1sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_1SIGMA);
		break;
	  case 2:
		points_2sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_2SIGMA);
		break;
	  case 3:
		points_3sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_3SIGMA);
		break;
	  default:
		points_1sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_1SIGMA);
		points_2sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_2SIGMA);
		points_3sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_3SIGMA);
		break;
	}

	{
	  FILE *f_PL;
	  gchar *pfile = NULL;

	  f_PL = nc_de_open_dataout_file (model, "PL", &pfile);

	  fprintf (f_PL, "# %s\n", full_cmd_line);
	  if (points_1sigma != NULL)
	  {
		ncm_fit_cr_points_print (points_1sigma, f_PL);
		fprintf (f_PL, "\n\n");
	  }
	  if (points_2sigma != NULL)
	  {
		ncm_fit_cr_points_print (points_2sigma, f_PL);
		fprintf (f_PL, "\n\n");
	  }
	  if (points_3sigma != NULL)
	  {
		ncm_fit_cr_points_print (points_3sigma, f_PL);
		fprintf (f_PL, "\n\n");
	  }

	  fclose (f_PL);

	  nc_message ("# PL file: %s \n", pfile);

	  g_free (pfile);
	}
	ncm_fit_cr_points_free (points_1sigma);
	ncm_fit_cr_points_free (points_2sigma);
	ncm_fit_cr_points_free (points_3sigma);
  }

  if (de_fit.nsigma_fisher >= 0 && (de_fit.bidim_cr[0] != -1) && (de_fit.bidim_cr[1] != -1))
  {
	GList *points_1sigma;
	GList *points_2sigma;
	GList *points_3sigma;
	points_1sigma = NULL;
	points_2sigma = NULL;
	points_3sigma = NULL;
	switch (de_fit.nsigma_fisher)
	{
	  case 1:
		points_1sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_1SIGMA);
		break;
	  case 2:
		points_2sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_2SIGMA);
		break;
	  case 3:
		points_3sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_3SIGMA);
		break;
	  default:
		points_1sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_1SIGMA);
		points_2sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_2SIGMA);
		points_3sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], NC_C_STATS_3SIGMA);
		break;
	}

	{
	  FILE *f_MFcr;
	  gchar *mcrfile = NULL;

	  f_MFcr = nc_de_open_dataout_file (model, "MF_cr", &mcrfile);

	  fprintf (f_MFcr, "# %s\n", full_cmd_line);
	  if (points_1sigma != NULL)
	  {
		ncm_fit_cr_points_print (points_1sigma, f_MFcr);
		fprintf (f_MFcr, "\n\n");
	  }
	  if (points_2sigma != NULL)
	  {
		ncm_fit_cr_points_print (points_2sigma, f_MFcr);
		fprintf (f_MFcr, "\n\n");
	  }
	  if (points_3sigma != NULL)
	  {
		ncm_fit_cr_points_print (points_3sigma, f_MFcr);
		fprintf (f_MFcr, "\n\n");
	  }

	  fclose (f_MFcr);

	  nc_message ("# MF confidence regions file: %s \n", mcrfile);

	  g_free (mcrfile);
	}

	ncm_fit_cr_points_free (points_1sigma);
	ncm_fit_cr_points_free (points_2sigma);
	ncm_fit_cr_points_free (points_3sigma);
  }

  if (de_data_cluster.print_mass_function == TRUE && ca_array != NULL)
  {
	gint i;
	FILE *f_mf;
	gchar *mfile = NULL;

	f_mf = nc_de_open_dataout_file (model, "MassFunction", &mfile);
	for (i = 0; i < ca_array->len; i++)
	{
	  NcData *dca_unbinned = g_ptr_array_index (ca_array, i);
	  nc_mass_function_print (dca_unbinned, model, f_mf, full_cmd_line);
	  fprintf (f_mf, "\n\n");
	}
	fclose (f_mf);
	g_ptr_array_free (ca_array, TRUE);
	ca_array = NULL;

	nc_message ("# MassFunction file: %s \n", mfile);

	g_free (mfile);
  }

  if (ca_array != NULL)
  {
	g_ptr_array_free (ca_array, TRUE);
  }

  ncm_model_free (NCM_MODEL (model));
  ncm_mset_free (fiduc);
  ncm_mset_free (mset);
  ncm_fit_free (fit);
  nc_likelihood_free (lh);
  nc_dataset_free0 (ds, TRUE);
  g_free (full_cmd_line);

  return 0;
}
