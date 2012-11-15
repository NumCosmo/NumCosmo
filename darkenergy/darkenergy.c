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
  NcDERunEntries de_run = NC_DE_RUN_ENTRIES;
  NcDEModelEntries de_model = NC_DE_MODEL_ENTRIES;
  NcDEDataSimpleEntries de_data_simple = NC_DE_DATA_SIMPLE_ENTRIES;
  NcDEDataClusterEntries de_data_cluster = NC_DE_DATA_CLUSTER_ENTRIES;
  NcDEFitEntries de_fit = NC_DE_FIT_ENTRIES;
  NcDistance *dist;
  NcmMSet *mset, *fiduc;
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry *de_model_entries = NULL;
  GOptionEntry *de_data_simple_entries = NULL;
  GOptionEntry *de_data_cluster_entries = NULL;
  GOptionEntry *de_fit_entries = NULL;
  GPtrArray *ca_array = NULL;
  gchar *full_cmd_line = NULL;
  gchar *runconf_cmd_line = NULL;

  ncm_cfg_init ();

  context = g_option_context_new ("- test the dark energy models");
  g_option_context_set_summary (context, "DE Summary <FIXME>");
  g_option_context_set_description (context, "DE Description <FIXME>");

  g_option_context_set_main_group (context, nc_de_opt_get_run_group (&de_run));
  g_option_context_add_group (context, nc_de_opt_get_model_group (&de_model, &de_model_entries));
  g_option_context_add_group (context, nc_de_opt_get_data_simple_group (&de_data_simple, &de_data_simple_entries));
  g_option_context_add_group (context, nc_de_opt_get_data_cluster_group (&de_data_cluster, &de_data_cluster_entries));
  g_option_context_add_group (context, nc_de_opt_get_fit_group (&de_fit, &de_fit_entries));


  {
	gint i;
	for (i = 0; i < argc; i++)
	{
	  if (strcmp (argv[i], "--runconf") == 0 || strcmp (argv[i], "-c") == 0)
	  {
		if (i + 1 == argc || argv[i+1] == NULL)
		{
		  fprintf (stderr, "Invalid run options:\n  %s.\n", error->message);
		  printf ("%s", g_option_context_get_help (context, TRUE, NULL));
		  g_option_context_free (context);
		  return 0;
		}
		else
		  de_run.runconf = argv[i + 1];
	  }
	}
  }

  if (de_run.runconf != NULL)
  {
	gchar **runconf_argv = g_new0 (gchar *, 1000);
	gint runconf_argc = 0;
	GKeyFile *runconf = g_key_file_new ();

	runconf_argv[0] = g_strdup (argv[0]);
	runconf_argc++;
	
	if (!g_key_file_load_from_file (runconf, de_run.runconf, G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS, &error))
	{
	  fprintf (stderr, "Invalid run configuration file: %s\n  %s\n", de_run.runconf, error->message);
	  printf ("%s", g_option_context_get_help (context, TRUE, NULL));
	  return 0;
	}

	ncm_cfg_keyfile_to_arg (runconf, "DarkEnergy Model", de_model_entries,        runconf_argv, &runconf_argc);
	ncm_cfg_keyfile_to_arg (runconf, "Data Simple",      de_data_simple_entries,  runconf_argv, &runconf_argc);
	ncm_cfg_keyfile_to_arg (runconf, "Data Cluster",     de_data_cluster_entries, runconf_argv, &runconf_argc);
	ncm_cfg_keyfile_to_arg (runconf, "Fit",              de_fit_entries,          runconf_argv, &runconf_argc);
	
	g_key_file_free (runconf);

	runconf_cmd_line = ncm_cfg_command_line (&runconf_argv[1], runconf_argc - 1);
	if (!g_option_context_parse (context, &runconf_argc, &runconf_argv, &error))
	{
	  fprintf (stderr, "Invalid configuration file options:\n  %s.\n", error->message);
	  printf ("%s", g_option_context_get_help (context, TRUE, NULL));
	  g_option_context_free (context);
	  return 0;
	}
  }

  full_cmd_line = ncm_cfg_command_line (argv, argc);
  if (!g_option_context_parse (context, &argc, &argv, &error))
  {
	fprintf (stderr, "Invalid configuration file options:\n  %s.\n", error->message);
	printf ("%s", g_option_context_get_help (context, TRUE, NULL));
	g_option_context_free (context);
	return 0;
  }

  if (runconf_cmd_line != NULL)
  {
	gchar *tmp = g_strconcat (full_cmd_line, " ", runconf_cmd_line, NULL);
	g_free (full_cmd_line);
	g_free (runconf_cmd_line);
	runconf_cmd_line = NULL;
	full_cmd_line = tmp;
  }
  
  if (de_run.saverun != NULL)
  {
	GKeyFile *runconf = g_key_file_new ();
	gchar *runconf_data;
	gsize len;

    ncm_cfg_entries_to_keyfile (runconf, "DarkEnergy Model", de_model_entries);
    ncm_cfg_entries_to_keyfile (runconf, "Data Simple",      de_data_simple_entries);
    ncm_cfg_entries_to_keyfile (runconf, "Data Cluster",     de_data_cluster_entries);
    ncm_cfg_entries_to_keyfile (runconf, "Fit",              de_fit_entries);

	runconf_data = g_key_file_to_data (runconf, &len, &error);
	if (error != NULL)
	  fprintf (stderr, "Error converting options to configuration file:\n  %s\n", error->message);
	if (!g_file_set_contents (de_run.saverun, runconf_data, len, &error))
	  fprintf (stderr, "Error saving configuration file to disk:\n  %s\n", error->message);
	g_free (runconf_data);
	g_key_file_free (runconf);
  }

  if (de_fit.file_out != NULL)
	ncm_cfg_set_logfile (de_fit.file_out);

  if (de_data_simple.snia_id == NULL &&
      de_data_simple.cmb_id == NULL &&
      de_data_simple.bao_id == NULL &&
      de_data_simple.H_id == NULL &&
      de_data_simple.cluster_id == NULL)
  {
	printf ("No action or data was chosen.\n");
	printf ("%s", g_option_context_get_help (context, TRUE, NULL));
	g_option_context_free (context);
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

  if (de_model.flat)
  {
	nc_hicosmo_de_omega_x2omega_k (model);
	ncm_model_param_set (NCM_MODEL (model), NC_HICOSMO_DE_OMEGA_X, 0.0);
	ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FIXED);
  }
  else if (de_model.Omega_k)
  {
	nc_hicosmo_de_omega_x2omega_k (model);
	ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FREE);
  }

  if (de_model.pos_Omega_x)
  {
	nc_prior_add_positive (lh, NC_HICOSMO_ID, NC_HICOSMO_DE_OMEGA_X);
  }

  if (de_data_simple.snia_id != NULL)
  {
	const GEnumValue *snia_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_DISTANCE_MU_SN_IA_ID,
	                                                              de_data_simple.snia_id);
	if (snia_id != NULL)
	{
	  NcData *snia = nc_data_distance_mu_snia (dist, snia_id->value);
	  nc_dataset_append_data (ds, snia);
	}
	else
	  g_error ("Supernovae sample '%s' not found run --snia-list to list the available options", de_data_simple.snia_id);
  }

  if (de_data_simple.cmb_id != NULL)
  {
	const GEnumValue *cmb_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_CMB_ID,
	                                                             de_data_simple.cmb_id);
	if (cmb_id != NULL)
	{
	  NcData *cmb_data = nc_data_cmb (dist, cmb_id->value);
	  nc_dataset_append_data (ds, cmb_data);
	}
	else
	  g_error ("CMB sample '%s' not found run --snia-list to list the available options", de_data_simple.cmb_id);
  }

  if (de_data_simple.bao_id != NULL)
  {
	const GEnumValue *bao_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_BAO_ID,
	                                                             de_data_simple.bao_id);
	if (bao_id != NULL)
	{
	  NcData *bao_data = nc_data_bao (dist, bao_id->value);
	  nc_dataset_append_data (ds, bao_data);
	}
	else
	  g_error ("BAO sample '%s' not found run --snia-list to list the available options", de_data_simple.bao_id);
  }

  if (de_data_simple.H_id != NULL)
  {
	const GEnumValue *H_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_HUBBLE_ID,
	                                                           de_data_simple.H_id);
	if (H_id != NULL)
	{
	  NcData *H_data = nc_data_hubble (H_id->value);
	  nc_dataset_append_data (ds, H_data);
	}
	else
	  g_error ("Hubble sample '%s' not found run --snia-list to list the available options", de_data_simple.H_id);
  }

  if (de_data_simple.cluster_id != NULL)
  {
	const GEnumValue *cluster_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_CLUSTER_ABUNDANCE_ID,
	                                                                 de_data_simple.cluster_id);
	if (cluster_id != NULL)
	{
	ca_array = nc_de_data_cluster_new (dist, mset, &de_data_cluster, ds, cluster_id->value);
	}
	else
	  g_error ("Cluster sample '%s' not found run --snia-list to list the available options", de_data_simple.cluster_id);
  }

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
  de_fit.save_best_fit = (de_fit.save_best_fit || de_fit.save_fisher);
  
  if (de_fit.fit)
  {
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

	nc_message ("#---------------------------------------------------------------------------------- \n", bfile);
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

	  nc_message ("#---------------------------------------------------------------------------------- \n", mfile);
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

	nc_dataset_log_info (fit->lh->ds);
	ncm_mset_pretty_log (resample_mset);
		
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
		  prob_sigma = ncm_c_stats_1sigma ();
		  break;
		case 2:
		  prob_sigma = ncm_c_stats_2sigma ();
		  break;
		case 3:
		  prob_sigma = ncm_c_stats_3sigma ();
		  break;
		default:
		  prob_sigma = ncm_c_stats_1sigma ();
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
		points_1sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_1sigma ());
		break;
	  case 2:
		points_2sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_2sigma ());
		break;
	  case 3:
		points_3sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_3sigma ());
		break;
	  default:
		points_1sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_1sigma ());
		points_2sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_2sigma ());
		points_3sigma = ncm_fit_cr2 (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_3sigma ());
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
		points_1sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_1sigma ());
		break;
	  case 2:
		points_2sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_2sigma ());
		break;
	  case 3:
		points_3sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_3sigma ());
		break;
	  default:
		points_1sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_1sigma ());
		points_2sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_2sigma ());
		points_3sigma = ncm_fit_cr2_fisher (fit, NC_HICOSMO_ID, de_fit.bidim_cr[0], NC_HICOSMO_ID, de_fit.bidim_cr[1], ncm_c_stats_3sigma ());
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
