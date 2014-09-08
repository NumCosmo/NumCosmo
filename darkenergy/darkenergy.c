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

#include <unistd.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_blas.h>

gint
main (gint argc, gchar *argv[])
{
  NcHICosmo *model;
  NcmDataset *dset;
  NcmLikelihood *lh;
  NcmFit *fit;
  NcDERunEntries de_run = NC_DE_RUN_ENTRIES;
  NcDEModelEntries de_model = NC_DE_MODEL_ENTRIES;
  NcDEDataSimpleEntries de_data_simple = NC_DE_DATA_SIMPLE_ENTRIES;
  NcDEDataClusterEntries de_data_cluster = NC_DE_DATA_CLUSTER_ENTRIES;
  NcDEFitEntries de_fit = NC_DE_FIT_ENTRIES;
  NcDistance *dist;
  NcmMSet *mset, *fiduc = NULL;
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry *de_model_entries = NULL;
  GOptionEntry *de_data_simple_entries = NULL;
  GOptionEntry *de_data_cluster_entries = NULL;
  GOptionEntry *de_fit_entries = NULL;
  GPtrArray *ca_array = NULL;
  gchar *full_cmd_line = NULL;
  gchar *runconf_cmd_line = NULL;
  gboolean is_de = FALSE;
  NcmRNG *rng = ncm_rng_pool_get ("darkenergy");

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
        if (i + 1 == argc || argv[i + 1] == NULL)
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
    gchar **runconf_argv_m = g_new0 (gchar *, 1000);
    gint runconf_argc = 0;
    GKeyFile *runconf = g_key_file_new ();
    guint i;

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

    for (i = 0; i < runconf_argc; i++)
    {
      runconf_argv_m[i] = runconf_argv[i]; 
    }
    runconf_argv_m[i] = NULL;
    
    runconf_cmd_line = ncm_cfg_command_line (&runconf_argv[1], runconf_argc - 1);
    if (!g_option_context_parse (context, &runconf_argc, &runconf_argv, &error))
    {
      fprintf (stderr, "Invalid configuration file options:\n  %s.\n", error->message);
      printf ("%s", g_option_context_get_help (context, TRUE, NULL));
      g_option_context_free (context);
      return 0;
    }

    g_strfreev (runconf_argv_m);
    g_free (runconf_argv);
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
    printf ("No action or data were chosen.\n");
    printf ("%s", g_option_context_get_help (context, TRUE, NULL));
    g_option_context_free (context);
    return 0;
  }

  g_option_context_free (context);

  ncm_message ("# NumCosmo Version -- "NUMCOSMO_VERSION"\n");
  ncm_message ("# Command Line: %s\n", full_cmd_line);

  dset = ncm_dataset_new ();
  model = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, de_model.model_name);
  mset = ncm_mset_new (NCM_MODEL (model), NULL);
  dist = nc_distance_new (2.0);
  
  ncm_serialize_global_set (dist, "MainDist", FALSE);

  if (g_type_is_a (G_OBJECT_TYPE (model), NC_TYPE_HICOSMO_DE))
    is_de = TRUE;

  if (de_model.help_names)
  {
    gint i;
    ncm_message ("# Model name -- %s\n", ncm_model_name (NCM_MODEL (model)));
    for (i = 0; i < ncm_model_len (NCM_MODEL (model)); i++)
      ncm_message ("# Model parameter [%02d] = %s\n", i, ncm_model_param_name (NCM_MODEL (model), i));
    return 0;
  }

  lh = ncm_likelihood_new (dset);

  if (de_model.flat)
  {
    if (!is_de)
      g_error ("flat option is valid only for darkenergy models");
    nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (model));
    ncm_model_param_set (NCM_MODEL (model), NC_HICOSMO_DE_OMEGA_X, 0.0);
    ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FIXED);
  }
  else if (de_model.Omega_k)
  {
    if (!is_de)
      g_error ("omegak option is valid only for darkenergy models");
    nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (model));
    ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FREE);
  }

  if (de_model.pos_Omega_x)
  {
    if (!is_de)
      g_error ("omegak > 0 option is valid only for darkenergy models");
    ncm_prior_add_positive (lh, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_X);
  }
  
  if (de_data_simple.snia_id != NULL)
  {
    const GEnumValue *snia_id = 
      ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_SNIA_ID,
                                        de_data_simple.snia_id);
    if (snia_id == NULL)
      g_error ("Supernovae sample '%s' not found run --snia-list to list"
               " the available options", de_data_simple.snia_id);

    if (snia_id->value >= NC_DATA_SNIA_SIMPLE_START && snia_id->value <= NC_DATA_SNIA_SIMPLE_END)
    {
      NcmData *snia = nc_data_dist_mu_new (dist, snia_id->value);
      ncm_dataset_append_data (dset, snia);
      ncm_data_free (snia);
    }
    else if (snia_id->value >= NC_DATA_SNIA_COV_START && snia_id->value <= NC_DATA_SNIA_COV_END)
    {
      NcSNIADistCov *dcov;
      NcmData *data;

      if (de_data_simple.snia_objser == NULL)
        dcov = nc_snia_dist_cov_new (dist);
      else
        dcov = NC_SNIA_DIST_COV (ncm_serialize_global_from_string (de_data_simple.snia_objser));

      g_assert (NC_IS_SNIA_DIST_COV (dcov));

      data = nc_data_snia_cov_new (de_data_simple.snia_use_det);
      nc_data_snia_load_cat (NC_DATA_SNIA_COV (data), snia_id->value);
      
      ncm_mset_set (mset, NCM_MODEL (dcov));
      ncm_dataset_append_data (dset, data);
      ncm_data_free (data);
      nc_snia_dist_cov_free (dcov);
    }
    else
      g_assert_not_reached ();
  }

 
  if (de_data_simple.cmb_id != NULL)
  {
    const GEnumValue *cmb_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_CMB_ID,
                                                                 de_data_simple.cmb_id);
    if (cmb_id != NULL)
    {
      NcmData *cmb_data = nc_data_cmb_create (dist, cmb_id->value);
      ncm_dataset_append_data (dset, cmb_data);
      ncm_data_free (cmb_data);
    }
    else
      g_error ("CMB sample '%s' not found run --cmb-list to list the available options", de_data_simple.cmb_id);
  }

  if (de_data_simple.bao_id != NULL)
  {
    guint i;
    guint nbao = g_strv_length (de_data_simple.bao_id);
    for (i = 0; i < nbao; i++)
    {
      gchar *bao_id_i = de_data_simple.bao_id[i];
      const GEnumValue *bao_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_BAO_ID, bao_id_i);
      if (bao_id != NULL)
      {
        NcmData *bao_data = nc_data_bao_create (dist, bao_id->value);
        ncm_dataset_append_data (dset, bao_data);
        ncm_data_free (bao_data);
      }
      else
        g_error ("BAO sample '%s' not found run --bao-list to list the available options", bao_id_i);
    }
  }

  if (de_data_simple.H_id != NULL)
  {
    guint i;
    guint nHz = g_strv_length (de_data_simple.H_id);
    
    for (i = 0; i < nHz; i++)
    {
      gchar *H_id_i = de_data_simple.H_id[i];
      const GEnumValue *H_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_HUBBLE_ID, H_id_i);
      if (H_id != NULL)
      {
        NcmData *H_data = nc_data_hubble_new (H_id->value);
        ncm_dataset_append_data (dset, H_data);
        ncm_data_free (H_data);
      }
      else
        g_error ("Hubble sample '%s' not found run --H-list to list the available options", H_id_i);
    }
  }

  if (de_data_simple.H_BAO_id != NULL)
  {
    guint i;
    guint nHrs = g_strv_length (de_data_simple.H_BAO_id);
    
    for (i = 0; i < nHrs; i++)
    {
      gchar *Hrs_id_i = de_data_simple.H_BAO_id[i];
      const GEnumValue *Hrs_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_HUBBLE_BAO_ID, Hrs_id_i);
      if (Hrs_id != NULL)
      {
        NcmData *Hrs_data = nc_data_hubble_bao_new (dist, Hrs_id->value);
        ncm_dataset_append_data (dset, Hrs_data);
        ncm_data_free (Hrs_data);
      }
      else
        g_error ("Hubble BAO sample '%s' not found run --H-BAO-list to list the available options", Hrs_id_i);
    }
  }

  if (de_data_simple.cluster_id != NULL)
  {
    const GEnumValue *cluster_id = ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_CLUSTER_ABUNDANCE_ID,
                                                                     de_data_simple.cluster_id);
    if (cluster_id != NULL)
    {
      ca_array = nc_de_data_cluster_new (dist, mset, &de_data_cluster, dset, cluster_id->value, rng);
    }
    else
      g_error ("Cluster sample '%s' not found run --cluster-list to list the available options", de_data_simple.cluster_id);

  }

  if (de_data_simple.BBN)
    nc_hicosmo_de_new_add_bbn (lh);

  if (de_data_simple.BBN_Ob)
  {
    NcmMSetFunc *Omega_bh2 = nc_hicosmo_create_mset_func0 (&nc_hicosmo_Omega_bh2);
    ncm_prior_add_gaussian_const_func (lh, Omega_bh2, 0.022, 0.002);
    ncm_mset_func_free (Omega_bh2);
  }

  if (de_fit.qspline_cp && TRUE)
  {
    if (!NC_IS_HICOSMO_QSPLINE (model))
      g_error ("Continuity priors are only valid for NcHICosmoQSPline model");
    else
    {
      NcHICosmoQSplineContPrior *qspline_cp = 
        nc_hicosmo_qspline_add_continuity_priors (NC_HICOSMO_QSPLINE (model), lh, 1.0e-10, de_fit.qspline_cp_sigma);
      nc_hicosmo_qspline_cont_prior_free (qspline_cp);
    }
  }

  if (de_fit.fiducial != NULL)
    fiduc = ncm_mset_load (de_fit.fiducial);
  else
    fiduc = ncm_mset_ref (mset);

  {
    const GEnumValue *fit_type_id = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_TYPE, de_fit.fit_type);
    const GEnumValue *fit_diff_id = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_GRAD_TYPE, de_fit.fit_diff);
    if (fit_type_id == NULL)
      g_error ("Fit type '%s' not found run --fit-list to list the available options", de_fit.fit_type);
    if (fit_diff_id == NULL)
      g_error ("Fit type '%s' not found run --fit-list to list the available options", de_fit.fit_diff);

    fit = ncm_fit_new (fit_type_id->value, de_fit.fit_algo, lh, mset, fit_diff_id->value);
    ncm_fit_set_m2lnL_reltol (fit, de_fit.fit_reltol);
    ncm_fit_set_params_reltol (fit, de_fit.fit_params_reltol);
  }

  if (de_fit.qspline_cp && FALSE)
  {
    if (!NC_IS_HICOSMO_QSPLINE (model))
      g_error ("Continuity priors are only valid for NcHICosmoQSPline model");
    else
    {
      NcHICosmoQSplineContPrior *qspline_cp = 
        nc_hicosmo_qspline_add_continuity_constraints (NC_HICOSMO_QSPLINE (model), fit, de_fit.qspline_cp_sigma);
      nc_hicosmo_qspline_cont_prior_free (qspline_cp);
    }
  }

  if (de_data_simple.priors_gauss != NULL)
  {
    guint i;
    guint npriors = g_strv_length (de_data_simple.priors_gauss);
    
    for (i = 0; i < npriors; i++)
    {
      GError *error = NULL;
      gchar *priors_str = de_data_simple.priors_gauss[i];
      GVariant *prior_hash = g_variant_parse (G_VARIANT_TYPE ("a{sv}"), priors_str, NULL, NULL, &error);
      if (prior_hash == NULL)
        g_error ("Invalid prior string: %s.", error->message);

      {
        NcmMSetPIndex p_i;
        gchar *model_ns = NULL;
        gchar *p_name = NULL;
        gdouble mu, sigma;
        
        if (!g_variant_lookup (prior_hash, "model", "s", &model_ns))
          g_error ("Prior must contain `model' key.");
        
        p_i.mid = ncm_mset_get_id_by_ns (model_ns);
        if (p_i.mid < 0)
          g_error ("Model %s not found.", model_ns);

        if (ncm_mset_peek (mset, p_i.mid) == NULL)
          g_error ("Model `%s' not present in NcmMSet.", model_ns);

        if (!g_variant_lookup (prior_hash, "param", "s", &p_name))
          g_error ("Prior must contain `param' key.");

        if (!ncm_model_param_index_from_name (ncm_mset_peek (mset, p_i.mid), p_name, &p_i.pid))
          g_error ("Parameter `%s' not found in model `%s'.", p_name, model_ns);
        
        if (!g_variant_lookup (prior_hash, "mean", "d", &mu))
          g_error ("Prior must contain `mean' key.");
        
        if (!g_variant_lookup (prior_hash, "sigma", "d", &sigma))
          g_error ("Prior must contain `sigma' key.");

        g_assert_cmpfloat (sigma, >, 0.0);

        ncm_prior_add_gaussian_data (lh, p_i.mid, p_i.pid, mu, sigma);
        g_variant_unref (prior_hash);
        g_free (model_ns);
        g_free (p_name);
      }
    }
  }


  
  /* All initializations done! */

  if (de_fit.resample)
  {
    NcmRNG *resample_rng;
    if (de_fit.mc_seed > -1)
      resample_rng = ncm_rng_seeded_new (NULL, de_fit.mc_seed);
    else
      resample_rng = ncm_rng_new (NULL);

    ncm_cfg_msg_sepa ();
    ncm_message ("# Resampling from fiducial model.\n");
    ncm_dataset_resample (dset, fiduc, resample_rng);
    ncm_rng_free (resample_rng);
  }

  de_fit.fisher = (de_fit.fisher || (de_fit.nsigma_fisher != -1) || (de_fit.nsigma != -1) || (de_fit.onedim_cr != NULL));
  de_fit.fit = (de_fit.fit || de_fit.fisher);
  de_fit.save_best_fit = (de_fit.save_best_fit || de_fit.save_fisher);

  if (de_fit.fit)
  {
    ncm_fit_set_maxiter (fit, de_fit.max_iter);
    ncm_fit_run (fit, de_fit.msg_level);
    
    ncm_fit_log_info (fit);
  }

  if (de_fit.save_best_fit)
  {
    FILE *f_bf;
    gchar *bfile = NULL;

    f_bf = nc_de_open_dataout_file (model, "best_fit", &bfile);
    ncm_mset_params_pretty_print (fit->mset, f_bf, full_cmd_line);
    fclose (f_bf);

    ncm_cfg_msg_sepa ();
    ncm_message ("# Params file: %s \n", bfile);

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

      ncm_message ("#---------------------------------------------------------------------------------- \n", mfile);
      ncm_message ("# FM file: %s \n", mfile);

      g_free (mfile);

    }
  }

  if (de_fit.montecarlo > 0)
  {
    if (de_fit.mcbs <= 0)
    {
      NcmFitMC *mc = ncm_fit_mc_new (fit, de_fit.mc_rtype, de_fit.msg_level);
      gdouble m2lnL = 0.0;
      if (fit->fstate->is_best_fit)
        m2lnL = ncm_fit_state_get_m2lnL_curval (fit->fstate);

      if (de_fit.fiducial != NULL)
        ncm_fit_mc_set_fiducial (mc, fiduc);

      if (de_fit.mc_nthreads > 1)
        ncm_fit_mc_set_nthreads (mc, de_fit.mc_nthreads);

      if (de_fit.mc_seed > -1)
      {
        NcmRNG *mc_rng = ncm_rng_seeded_new (NULL, de_fit.mc_seed);
        ncm_fit_mc_set_rng (mc, mc_rng);
        ncm_rng_free (mc_rng);
      }
      
      if (de_fit.mc_data != NULL)
        ncm_fit_mc_set_data_file (mc, de_fit.mc_data);
      
      ncm_fit_mc_start_run (mc);
      
      if (de_fit.mc_ni >= 0)
        ncm_fit_mc_set_first_sample_id (mc, de_fit.mc_ni);

      ncm_fit_mc_run_lre (mc, de_fit.montecarlo, de_fit.mc_lre);
      ncm_fit_mc_end_run (mc);
      
      ncm_fit_mc_mean_covar (mc);
      ncm_fit_catalog_param_pdf (mc->fcat, 0);
      ncm_fit_log_covar (fit);
      {
        gdouble p_value = ncm_fit_catalog_param_pdf_pvalue (mc->fcat, m2lnL, FALSE);
        ncm_message ("#   - pvalue for fitted model [% 20.15g] %04.2f%%.\n#\n", m2lnL, 100.0 * p_value);
      }
      ncm_fit_mc_clear (&mc);
    }
    else
    {
      NcmMSet *resample_mset;
      NcmFitMCBS *mcbs = ncm_fit_mcbs_new (fit);
      gdouble m2lnL = ncm_fit_state_get_m2lnL_curval (fit->fstate);

      if (de_fit.fiducial != NULL)
        resample_mset = fiduc;
      else
        resample_mset = fit->mset;

      if (de_fit.mc_seed > -1)
      {
        NcmRNG *mc_rng = ncm_rng_seeded_new (NULL, de_fit.mc_seed);
        ncm_fit_mcbs_set_rng (mcbs, mc_rng);
        ncm_rng_free (mc_rng);
      }
      
      if (de_fit.mc_data != NULL)
        ncm_fit_mcbs_set_filename (mcbs, de_fit.mc_data);

      ncm_fit_mcbs_run (mcbs, resample_mset, de_fit.mc_ni, de_fit.montecarlo, de_fit.mcbs, de_fit.mc_rtype, de_fit.msg_level, de_fit.mc_nthreads);

      ncm_fit_catalog_param_pdf (mcbs->fcat, 0);
      ncm_fit_log_covar (fit);
      {
        gdouble p_value = ncm_fit_catalog_param_pdf_pvalue (mcbs->fcat, m2lnL, FALSE);
        ncm_message ("#   - pvalue for fitted model [% 20.15g] %04.2f%%.\n#\n", m2lnL, 100.0 * p_value);
      }
      ncm_fit_mcbs_clear (&mcbs);
    }
  }

  if (de_fit.mcmc)
  {
    NcmMCSamplerGauss *mcsg = ncm_mc_sampler_gauss_new (0);
    NcmFitMCMC *mcmc = ncm_fit_mcmc_new (fit, NCM_MC_SAMPLER (mcsg), de_fit.msg_level);

    if (de_fit.fisher)
    {
      ncm_mc_sampler_gauss_set_cov (mcsg, fit->fstate->covar);
    }
    else
    {
      ncm_mc_sampler_gauss_set_cov_from_scale (mcsg);
    }
    
    if (de_fit.mc_seed > -1)
    {
      NcmRNG *mcmc_rng = ncm_rng_seeded_new (NULL, de_fit.mc_seed);
      ncm_fit_mcmc_set_rng (mcmc, mcmc_rng);
      ncm_rng_free (mcmc_rng);
    }

    if (de_fit.mc_data != NULL)
      ncm_fit_mcmc_set_data_file (mcmc, de_fit.mc_data);

    ncm_fit_mcmc_start_run (mcmc);
    ncm_fit_mcmc_run_lre (mcmc, 0, de_fit.mc_lre);
    ncm_fit_mcmc_end_run (mcmc);

    ncm_fit_mcmc_mean_covar (mcmc);
    ncm_fit_catalog_param_pdf (mcmc->fcat, 0);
    ncm_fit_log_covar (fit);
    ncm_fit_mcmc_clear (&mcmc);    
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
      {
        NcmMSetPIndex pi = {nc_hicosmo_id (), p_n};
        NcmLHRatio1d *lhr1d = ncm_lh_ratio1d_new (fit, &pi);
        ncm_lh_ratio1d_find_bounds (lhr1d, fit->mtype, prob_sigma, &err_inf, &err_sup);
        ncm_lh_ratio1d_free (lhr1d);
      }
      //ncm_fit_cr_1dim (fit, nc_hicosmo_id (), p_n, prob_sigma, 1, &err_inf, &err_sup);
      ncm_message ("#  One dimension confidence region for %s[%02d] = % .5g (% .5g, % .5g)\n",
                   ncm_model_param_name (NCM_MODEL (model), p_n), p_n,
                   ncm_model_param_get (NCM_MODEL (model), p_n), err_inf, err_sup);
    }
  }

  if (de_fit.nsigma >= 0 && (de_fit.bidim_cr[0] != -1) && (de_fit.bidim_cr[1] != -1))
  {
    NcmLHRatio2d *lhr2d;
    NcmMSetPIndex pi1;
    NcmMSetPIndex pi2;
    NcmLHRatio2dRegion *rg_1sigma = NULL;
    NcmLHRatio2dRegion *rg_2sigma = NULL;
    NcmLHRatio2dRegion *rg_3sigma = NULL;

    pi1.mid = nc_hicosmo_id ();
    pi2.mid = nc_hicosmo_id ();
    pi1.pid = de_fit.bidim_cr[0];
    pi2.pid = de_fit.bidim_cr[1];

    lhr2d = ncm_lh_ratio2d_new (fit, &pi1, &pi2);
    
    switch (de_fit.nsigma)
    {
      case 1:
      {
        rg_1sigma = ncm_lh_ratio2d_conf_region (lhr2d, ncm_c_stats_1sigma (), 0.0, de_fit.msg_level);
        break;
      }
      case 2:
        rg_2sigma = ncm_lh_ratio2d_conf_region (lhr2d, ncm_c_stats_2sigma (), 0.0, de_fit.msg_level);
        break;
      case 3:
        rg_3sigma = ncm_lh_ratio2d_conf_region (lhr2d, ncm_c_stats_3sigma (), 0.0, de_fit.msg_level);
        break;
      default:
        rg_1sigma = ncm_lh_ratio2d_conf_region (lhr2d, ncm_c_stats_1sigma (), 0.0, de_fit.msg_level);
        rg_2sigma = ncm_lh_ratio2d_conf_region (lhr2d, ncm_c_stats_2sigma (), 0.0, de_fit.msg_level);
        rg_3sigma = ncm_lh_ratio2d_conf_region (lhr2d, ncm_c_stats_3sigma (), 0.0, de_fit.msg_level);
        break;
    }

    {
      FILE *f_PL;
      gchar *pfile = NULL;

      f_PL = nc_de_open_dataout_file (model, "PL", &pfile);

      fprintf (f_PL, "# %s\n", full_cmd_line);
      if (rg_1sigma != NULL)
      {
       ncm_lh_ratio2d_region_print (rg_1sigma, f_PL);
        fprintf (f_PL, "\n\n");
      }
      if (rg_2sigma != NULL)
      {
        ncm_lh_ratio2d_region_print (rg_2sigma, f_PL);
        fprintf (f_PL, "\n\n");
      }
      if (rg_3sigma != NULL)
      {
        ncm_lh_ratio2d_region_print (rg_3sigma, f_PL);
        fprintf (f_PL, "\n\n");
      }

      fclose (f_PL);

      ncm_message ("# PL file: %s \n", pfile);

      g_free (pfile);
    }
    ncm_lh_ratio2d_region_clear (&rg_1sigma);
    ncm_lh_ratio2d_region_clear (&rg_2sigma);
    ncm_lh_ratio2d_region_clear (&rg_3sigma);
    ncm_lh_ratio2d_free (lhr2d);
  }

  if (de_fit.nsigma_fisher >= 0 && (de_fit.bidim_cr[0] != -1) && (de_fit.bidim_cr[1] != -1))
  {
    NcmLHRatio2d *lhr2d;
    NcmMSetPIndex pi1;
    NcmMSetPIndex pi2;
    NcmLHRatio2dRegion *rg_1sigma = NULL;
    NcmLHRatio2dRegion *rg_2sigma = NULL;
    NcmLHRatio2dRegion *rg_3sigma = NULL;

    pi1.mid = nc_hicosmo_id ();
    pi2.mid = nc_hicosmo_id ();
    pi1.pid = de_fit.bidim_cr[0];
    pi2.pid = de_fit.bidim_cr[1];

    lhr2d = ncm_lh_ratio2d_new (fit, &pi1, &pi2);

    switch (de_fit.nsigma_fisher)
    {
      case 1:
        rg_1sigma = ncm_lh_ratio2d_fisher_border (lhr2d, ncm_c_stats_1sigma (), 0.0, de_fit.msg_level);
        break;
      case 2:
        rg_2sigma = ncm_lh_ratio2d_fisher_border (lhr2d, ncm_c_stats_2sigma (), 0.0, de_fit.msg_level);
        break;
      case 3:
        rg_3sigma = ncm_lh_ratio2d_fisher_border (lhr2d, ncm_c_stats_3sigma (), 0.0, de_fit.msg_level);
        break;
      default:
        rg_1sigma = ncm_lh_ratio2d_fisher_border (lhr2d, ncm_c_stats_1sigma (), 0.0, de_fit.msg_level);
        rg_2sigma = ncm_lh_ratio2d_fisher_border (lhr2d, ncm_c_stats_2sigma (), 0.0, de_fit.msg_level);
        rg_3sigma = ncm_lh_ratio2d_fisher_border (lhr2d, ncm_c_stats_3sigma (), 0.0, de_fit.msg_level);
        break;
    }

    {
      FILE *f_MFcr;
      gchar *mcrfile = NULL;

      f_MFcr = nc_de_open_dataout_file (model, "MF_cr", &mcrfile);

      fprintf (f_MFcr, "# %s\n", full_cmd_line);
      if (rg_1sigma != NULL)
      {
        ncm_lh_ratio2d_region_print (rg_1sigma, f_MFcr);
        fprintf (f_MFcr, "\n\n");
      }
      if (rg_2sigma != NULL)
      {
        ncm_lh_ratio2d_region_print (rg_2sigma, f_MFcr);
        fprintf (f_MFcr, "\n\n");
      }
      if (rg_3sigma != NULL)
      {
        ncm_lh_ratio2d_region_print (rg_3sigma, f_MFcr);
        fprintf (f_MFcr, "\n\n");
      }

      fclose (f_MFcr);

      ncm_message ("# MF confidence regions file: %s \n", mcrfile);

      g_free (mcrfile);
    }

    ncm_lh_ratio2d_region_clear (&rg_1sigma);
    ncm_lh_ratio2d_region_clear (&rg_2sigma);
    ncm_lh_ratio2d_region_clear (&rg_3sigma);
  }

  if (de_data_cluster.print_mass_function == TRUE && ca_array != NULL)
  {
    gint i;
    FILE *f_mf;
    gchar *mfile = NULL;

    f_mf = nc_de_open_dataout_file (model, "MassFunction", &mfile);
    for (i = 0; i < ca_array->len; i++)
    {
      NcDataClusterNCount *dca_unbinned = g_ptr_array_index (ca_array, i);
      nc_data_cluster_ncount_print (dca_unbinned, model, f_mf, full_cmd_line);
      fprintf (f_mf, "\n\n");
    }
    fclose (f_mf);
    g_ptr_array_free (ca_array, TRUE);
    ca_array = NULL;

    ncm_message ("# MassFunction file: %s \n", mfile);

    g_free (mfile);
  }

  if (de_fit.kinematics_sigma)
  {
    NcmMSetFunc *dec_param_func = nc_hicosmo_create_mset_func1 (nc_hicosmo_q);
    NcmMSetFunc *E2_func = nc_hicosmo_create_mset_func1 (nc_hicosmo_E2);
    NcmMSetFunc *Em2_func = nc_hicosmo_create_mset_func1 (nc_hicosmo_Em2);
    NcmMSetFunc *dec_func = nc_hicosmo_create_mset_func1 (nc_hicosmo_dec);
    NcmMSetFunc *wec_func = nc_hicosmo_create_mset_func1 (nc_hicosmo_wec);
    gint i;

    ncm_message ("# Kinematics data:\n");
    for (i = 0; i < de_fit.kinematics_n; i++)
    {
      gdouble z = de_fit.kinematics_z / (de_fit.kinematics_n - 1.0) * i;
      gdouble q_z, E2_z, Em2_z, dec, wec;
      gdouble sigma_q_z, sigma_E2_z, sigma_Em2_z, sigma_dec, sigma_wec;
      ncm_fit_function_error (fit, dec_param_func, &z, FALSE, &q_z, &sigma_q_z);
      ncm_fit_function_error (fit, E2_func, &z, FALSE, &E2_z, &sigma_E2_z);
      ncm_fit_function_error (fit, Em2_func, &z, FALSE, &Em2_z, &sigma_Em2_z);
      ncm_fit_function_error (fit, dec_func, &z, FALSE, &dec, &sigma_dec);
      ncm_fit_function_error (fit, wec_func, &z, FALSE, &wec, &sigma_wec);
      ncm_message ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", 
                   z, 
                   q_z, sigma_q_z, 
                   E2_z, sigma_E2_z, 
                   Em2_z, sigma_Em2_z,
                   dec, sigma_dec,
                   wec, sigma_wec);
    }
    ncm_mset_func_free (dec_param_func);
    ncm_mset_func_free (E2_func);
    ncm_mset_func_free (Em2_func);
  }

  if (ca_array != NULL)
  {
    g_ptr_array_free (ca_array, TRUE);
  }

  if (de_fit.save_mset != NULL)
    ncm_mset_save (mset, de_fit.save_mset, TRUE);

  g_free (de_model_entries); de_model_entries = NULL;
  g_free (de_data_simple_entries); de_data_simple_entries = NULL;
  g_free (de_data_cluster_entries); de_data_cluster_entries = NULL;
  g_free (de_fit_entries); de_fit_entries = NULL;

  if (fiduc != NULL)
    ncm_mset_free (fiduc);
  
  ncm_serialize_global_reset ();
  ncm_model_free (NCM_MODEL (model));
  ncm_mset_free (mset);
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_rng_free (rng);
  nc_distance_free (dist);
  g_free (full_cmd_line);

  return 0;
}
