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
#include <gsl/gsl_cdf.h>

gint
main (gint argc, gchar *argv[])
{
  NcHICosmo *cosmo;
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
  NcmObjArray *funcs_oa = NULL;
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
  gboolean is_gcg = FALSE;
  gboolean is_idem2 = FALSE;  
  NcmRNG *rng = ncm_rng_pool_get ("darkenergy");
  NcmMSetCatalog *mcat = NULL;
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_cfg_init_full_ptr (&argc, &argv);
  
  context = g_option_context_new ("- test the dark energy models");
  g_option_context_set_summary (context, "general darkenergy and kinematic models analyzer");
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
      if ((strncmp (argv[i], "--runconf", 9) == 0) || (strcmp (argv[i], "-c") == 0))
      {
        guint rc_size = strlen (argv[i]);
        if ((rc_size == 9) || (rc_size == 2))
        {
          if (i + 1 == argc || argv[i + 1] == NULL)
          {
            gchar *msg = g_option_context_get_help (context, TRUE, NULL);
            fprintf (stderr, "Invalid run options: missing argument for `%s'\n", argv[i]);
            printf ("%s", msg);
            g_free (msg);
            g_option_context_free (context);
            return 0;
          }
          else
            de_run.runconf = argv[i + 1];
        }
        else if ((rc_size > 10) && (argv[i][9] == '='))
        {
          de_run.runconf = &argv[i][10];
        }
        else
        {
          gchar *msg = g_option_context_get_help (context, TRUE, NULL);
          fprintf (stderr, "Invalid run options:\n");
          printf ("%s", msg);
          g_free (msg);
          g_option_context_free (context);
          return 0;
        }
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

  if (de_run.main_seed >= 0)
    ncm_rng_set_seed (rng, de_run.main_seed);

  if (de_run.nthreads > 0)
    ncm_func_eval_set_max_threads (de_run.nthreads);
  
  if (de_data_simple.snia_id    == NULL &&
      de_data_simple.cmb_id     == NULL &&
      de_data_simple.bao_id     == NULL &&
      de_data_simple.H_id       == NULL &&
      de_data_simple.cluster_id == NULL &&
      de_data_simple.Planck     == NULL &&
      de_data_simple.data_files == NULL)
  {
    printf ("No action or data were chosen.\n");
    printf ("%s", g_option_context_get_help (context, TRUE, NULL));
    g_option_context_free (context);
    return 0;
  }

  g_option_context_free (context);

  ncm_message ("# NumCosmo Version -- "NUMCOSMO_VERSION"\n");
  ncm_message ("# Command Line: %s\n", full_cmd_line);

  dist = nc_distance_new (2.0);
  ncm_serialize_global_set (dist, "MainDist", FALSE);

  if (de_model.mset_file != NULL)
    mset = ncm_mset_load (de_model.mset_file, ser);
  else
    mset = ncm_mset_empty_new ();

  dset  = ncm_dataset_new ();
  if ((cosmo = NC_HICOSMO (ncm_mset_get (mset, nc_hicosmo_id ()))) == NULL)
  {
    if (de_model.model_name != NULL)
      cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, de_model.model_name);
    else
      cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  }

  if (ncm_mset_peek (mset, nc_hireion_id ()) == NULL)
  {
    NcHIReion *reion;

    if (de_model.model_reion != NULL)
      reion = nc_hireion_new_from_name (NC_TYPE_HIREION, de_model.model_reion);
    else
      reion = NC_HIREION (nc_hireion_camb_new ());
    
    ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
    nc_hireion_free (reion);
  }

  if (ncm_mset_peek (mset, nc_hiprim_id ()) == NULL)
  {
    NcHIPrim *prim;

    if (de_model.model_prim != NULL)
      prim = nc_hiprim_new_from_name (NC_TYPE_HIPRIM, de_model.model_prim);
    else
      prim = NC_HIPRIM (nc_hiprim_power_law_new ());
    
    ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
    nc_hiprim_free (prim);
  }
  
  if (ncm_mset_peek (mset, nc_hicosmo_id ()) == NULL)
    ncm_mset_push (mset, NCM_MODEL (cosmo));

  if (g_type_is_a (G_OBJECT_TYPE (cosmo), NC_TYPE_HICOSMO_DE))
    is_de = TRUE;
  else if (g_type_is_a (G_OBJECT_TYPE (cosmo), NC_TYPE_HICOSMO_GCG))
    is_gcg = TRUE;
  else if (g_type_is_a (G_OBJECT_TYPE (cosmo), NC_TYPE_HICOSMO_IDEM2))
    is_idem2 = TRUE;

  if (de_model.help_names)
  {
    gint i;
    ncm_message ("# Cosmological model name -- %s\n", ncm_model_name (NCM_MODEL (cosmo)));
    for (i = 0; i < ncm_model_len (NCM_MODEL (cosmo)); i++)
      ncm_message ("# Model parameter [%02d] = %s\n", i, ncm_model_param_name (NCM_MODEL (cosmo), i));
    return 0;
  }

  lh = ncm_likelihood_new (dset);

  if (de_model.flat)
  {
    if (is_de)
    {
      nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));
      ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
      ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FIXED);
    }
    else if (is_gcg)
    {
      nc_hicosmo_gcg_omega_x2omega_k (NC_HICOSMO_GCG (cosmo));
      ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_GCG_OMEGA_X, 0.0);
      ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_GCG_OMEGA_X, NCM_PARAM_TYPE_FIXED);
    }
    else if (is_idem2)
    {
      nc_hicosmo_idem2_omega_x2omega_k (NC_HICOSMO_IDEM2 (cosmo));
      ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_IDEM2_OMEGA_X, 0.0);
      ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_IDEM2_OMEGA_X, NCM_PARAM_TYPE_FIXED);
    }
    else
      g_error ("flat option is valid only for darkenergy models");
  }
  else if (de_model.Omega_k)
  {
    if (is_de)
    {
      nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));
      ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_X, NCM_PARAM_TYPE_FREE);
    }
    else if (is_gcg)
    {
      nc_hicosmo_gcg_omega_x2omega_k (NC_HICOSMO_GCG (cosmo));
      ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_GCG_OMEGA_X, NCM_PARAM_TYPE_FREE);
    }
    else if (is_idem2)
    {
      nc_hicosmo_idem2_omega_x2omega_k (NC_HICOSMO_IDEM2 (cosmo));
      ncm_mset_param_set_ftype (mset, nc_hicosmo_id (), NC_HICOSMO_IDEM2_OMEGA_X, NCM_PARAM_TYPE_FREE);
    }    
    else
      g_error ("omegak option is valid only for darkenergy models");
  }

  if (de_model.pos_Omega_x)
  {
    if (is_de)
    {
      ncm_likelihood_priors_add_flat_param (lh, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_X, 0.0, HUGE_VAL, 1.0);
    }
    else if (is_gcg)
    {
      ncm_likelihood_priors_add_flat_param (lh, nc_hicosmo_id (), NC_HICOSMO_GCG_OMEGA_X, 0.0, HUGE_VAL, 1.0);
    }
    else if (is_idem2)
    {
      ncm_likelihood_priors_add_flat_param (lh, nc_hicosmo_id (), NC_HICOSMO_IDEM2_OMEGA_X, 0.0, HUGE_VAL, 1.0);
    }
    else
      g_error ("omegak > 0 option is valid only for darkenergy models");
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
      NcDataDistMu *dist_mu = nc_data_dist_mu_new_from_id (dist, snia_id->value);
      ncm_dataset_append_data (dset, NCM_DATA (dist_mu));
      ncm_data_free (NCM_DATA (dist_mu));
    }
    else if (snia_id->value >= NC_DATA_SNIA_COV_START && snia_id->value <= NC_DATA_SNIA_COV_END)
    {
      NcmData *data = nc_data_snia_cov_new (de_data_simple.snia_use_det);
      guint sigma_int_len;
      
      nc_data_snia_cov_load_cat (NC_DATA_SNIA_COV (data), snia_id->value);

      sigma_int_len = nc_data_snia_cov_sigma_int_len (NC_DATA_SNIA_COV (data));

      if (ncm_mset_peek (mset, nc_snia_dist_cov_id ()) == NULL)
      {
        NcSNIADistCov *dcov;
        if (de_data_simple.snia_objser == NULL)
        {
          dcov = nc_snia_dist_cov_new (dist, sigma_int_len);
          nc_snia_dist_cov_set_default_params_by_id (dcov, snia_id->value);
        }
        else
          dcov = NC_SNIA_DIST_COV (ncm_serialize_global_from_string (de_data_simple.snia_objser));

        g_assert (NC_IS_SNIA_DIST_COV (dcov));
        ncm_mset_set (mset, NCM_MODEL (dcov));
        nc_snia_dist_cov_free (dcov);
      }
      ncm_dataset_append_data (dset, data);
      ncm_data_free (data);
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
        NcDataHubble *H_data = nc_data_hubble_new_from_id (H_id->value);
        ncm_dataset_append_data (dset, NCM_DATA (H_data));
        ncm_data_free (NCM_DATA (H_data));
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

  if (de_data_simple.Planck != NULL)
  {
    guint i;
    guint nPlanck = g_strv_length (de_data_simple.Planck);
    NcHIPertBoltzmannCBE *boltz = nc_hipert_boltzmann_cbe_new ();

    if (ncm_mset_peek (mset, nc_planck_fi_id ()) == NULL)
    {
      NcPlanckFI *planck_fi = nc_planck_fi_new_from_name (de_data_simple.PlanckFI == NULL ? "NcPlanckFICorTT" : de_data_simple.PlanckFI);
      ncm_mset_push (mset, NCM_MODEL (planck_fi));
      nc_planck_fi_free (planck_fi);
    }
    
    for (i = 0; i < nPlanck; i++)
    {
      gchar *Planck_i       = de_data_simple.Planck[i];
      NcDataPlanckLKL *plik = nc_data_planck_lkl_full_new (Planck_i, NC_HIPERT_BOLTZMANN (boltz));
      ncm_dataset_append_data (dset, NCM_DATA (plik));

      ncm_data_free (NCM_DATA (plik));
    }

    nc_hipert_boltzmann_cbe_free (boltz);
  }

  if (de_data_simple.data_files != NULL)
  {
    guint ndata_files = g_strv_length (de_data_simple.data_files);
    guint i;
    
    for (i = 0; i < ndata_files; i++)
    {
      if (!g_file_test (de_data_simple.data_files[i], G_FILE_TEST_EXISTS))
      {
        g_warning ("data file: `%s' not found, skipping.", de_data_simple.data_files[i]);
        continue;
      }
      else
      {
        NcmData *data = ncm_data_new_from_file (de_data_simple.data_files[i]);
        ncm_dataset_append_data (dset, data);
        ncm_data_free (data);        
      }
    }
    
  }
  
  if (de_data_simple.BBN)
    nc_hicosmo_de_new_add_bbn (lh);

  if (de_data_simple.BBN_Ob)
  {
    NcmMSetFunc *Omega_b0h2 = NCM_MSET_FUNC (ncm_mset_func_list_new ("NcHICosmo:Omega_b0h2", NULL));
    ncm_likelihood_priors_add_gauss_func (lh, Omega_b0h2, 0.022, 0.002, 0.0);
    ncm_mset_func_free (Omega_b0h2);
  }

  if (de_fit.qspline_cp)
  {
    if (!NC_IS_HICOSMO_QSPLINE (cosmo))
      g_error ("Continuity priors are only valid for NcHICosmoQSPline model");
    else
    {
      NcHICosmoQSplineContPrior *qspline_cp = 
        nc_hicosmo_qspline_add_continuity_priors (NC_HICOSMO_QSPLINE (cosmo), lh, 1.0e-10, de_fit.qspline_cp_sigma);

      ncm_mset_set (mset, NCM_MODEL (qspline_cp));
      nc_hicosmo_qspline_cont_prior_free (qspline_cp);
    }
  }

  if (de_data_simple.PlanckPriors)
  {
    NcPlanckFI *planck_fi = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
    if (planck_fi == NULL)
      g_warning ("Planck Priors enabled but not NcPlanckFI model set.");
    else
    {
      if (NC_IS_PLANCK_FI_COR_TT (planck_fi))
      {
        nc_planck_fi_cor_tt_add_all_default_priors (lh);
      }
    }
  }

  if (de_fit.fiducial != NULL)
    fiduc = ncm_mset_load (de_fit.fiducial, ser);
  else
    fiduc = ncm_mset_ref (mset);

  if (de_fit.fit_type == NULL)
  {
#ifdef NUMCOSMO_HAVE_NLOPT
    de_fit.fit_type = g_strdup ("nlopt");
#else
    de_fit.fit_type = g_strdup ("gsl-mms");
#endif
  }
  if (de_fit.fit_diff == NULL)
    de_fit.fit_diff = g_strdup ("numdiff-forward");
  
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

  if (de_data_simple.priors_gauss != NULL)
  {
    guint i;
    guint npriors = g_strv_length (de_data_simple.priors_gauss);
    
    for (i = 0; i < npriors; i++)
    {
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

        ncm_likelihood_priors_add_gauss_param (lh, p_i.mid, p_i.pid, mu, sigma);
        g_variant_unref (prior_hash);
        g_free (model_ns);
        g_free (p_name);
      }
    }
  }

  /* All initializations done! */

  if (de_fit.resample)
  {
    ncm_cfg_msg_sepa ();
    ncm_message ("# Resampling from fiducial model.\n");
    ncm_dataset_resample (dset, fiduc, rng);
  }

  de_fit.fisher        = !de_fit.fisher && ((de_fit.nsigma_fisher != -1) || (de_fit.nsigma != -1) || (de_fit.onedim_cr != NULL)) ? 1 : de_fit.fisher;
  de_fit.fit           = (de_fit.fit || de_fit.fisher);
  de_fit.save_best_fit = (de_fit.save_best_fit || de_fit.save_fisher);

  if (de_fit.fit)
  {
    ncm_fit_set_maxiter (fit, de_fit.max_iter);
    if (de_fit.restart)
      ncm_fit_run_restart (fit, de_fit.msg_level, de_fit.restart_abstol, de_fit.restart_reltol, NULL, de_fit.save_mset);
    else
      ncm_fit_run (fit, de_fit.msg_level);
    
    ncm_fit_log_info (fit);
  }

  if (de_fit.save_best_fit)
  {
    FILE *f_bf;
    gchar *bfile = NULL;

    f_bf = nc_de_open_dataout_file (cosmo, "best_fit", &bfile);
    ncm_mset_params_pretty_print (fit->mset, f_bf, full_cmd_line);
    fclose (f_bf);

    ncm_cfg_msg_sepa ();
    ncm_message ("# Params file: %s \n", bfile);

    g_free (bfile);
  }

  if (de_fit.fisher)
  {
    switch (de_fit.fisher)
    {
      case 1:
        ncm_fit_obs_fisher (fit);
        break;
      case 2:
        ncm_fit_fisher (fit);
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    ncm_fit_log_covar (fit);
    
    if (de_fit.save_fisher)
    {
      FILE *f_MF;
      gchar *mfile = NULL;

      f_MF = nc_de_open_dataout_file (cosmo, "MF", &mfile);
      ncm_fit_fishermatrix_print (fit, f_MF, full_cmd_line);
      fclose (f_MF);

      ncm_message ("#---------------------------------------------------------------------------------- \n");
      ncm_message ("# FM file: %s\n", mfile);

      g_free (mfile);

    }
  }

  if (de_fit.funcs != NULL)
  {
    const guint len = g_strv_length (de_fit.funcs);
    guint i;

    funcs_oa = ncm_obj_array_new ();

    for (i = 0; i < len; i++)
    {
      NcmMSetFunc *func = NULL;
      gdouble *x = NULL;
      gchar *func_name;
      guint len;

      func_name = ncm_util_function_params (de_fit.funcs[i], &x, &len);
      if (func_name == NULL)
        g_error ("darkenergy: invalid function name: `%s'.", de_fit.funcs[i]);

      if (ncm_mset_func_list_has_ns_name ("NcHICosmo", func_name))
      {
        func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcHICosmo", func_name, NULL));
      }
      else if (ncm_mset_func_list_has_ns_name ("NcDistance", func_name))
      {
        func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcDistance", func_name, G_OBJECT (dist)));
      }
      else
      {
        g_warning ("darkenergy: function `%s' not found, skipping!", de_fit.funcs[i]);
        func = NULL;
        continue;
      }

      if (len > 0)
        ncm_mset_func_set_eval_x (func, x, len);

      g_clear_pointer (&func_name, g_free);
      g_clear_pointer (&x, g_free);

      g_assert (ncm_mset_func_is_scalar (func));
      g_assert (ncm_mset_func_is_const (func));

      ncm_obj_array_add (funcs_oa, G_OBJECT (func));
    }
  }
  
  if (de_fit.mc)
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

    if (de_fit.mc_unordered)
      ncm_fit_mc_keep_order (mc, FALSE);

    ncm_fit_mc_start_run (mc);

    if (de_fit.mc_ni >= 0)
      ncm_fit_mc_set_first_sample_id (mc, de_fit.mc_ni);

    ncm_fit_mc_run_lre (mc, de_fit.mc_prerun, de_fit.mc_lre);
    ncm_fit_mc_end_run (mc);

    ncm_fit_mc_mean_covar (mc);
    ncm_mset_catalog_param_pdf (mc->mcat, 0);
    ncm_fit_log_covar (fit);
    {
      gdouble p_value = ncm_mset_catalog_param_pdf_pvalue (mc->mcat, m2lnL, FALSE);
      ncm_message ("#   - pvalue for fitted model [% 20.15g] %04.2f%%.\n#\n", m2lnL, 100.0 * p_value);
    }
    ncm_mset_catalog_clear (&mcat);
    mcat = ncm_fit_mc_get_catalog (mc);
    ncm_fit_mc_clear (&mc);
  }

  if (de_fit.mcbs)
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

    ncm_fit_mcbs_run (mcbs, resample_mset, de_fit.mc_ni, de_fit.mc_prerun, de_fit.mcbs_nbootstraps, de_fit.mc_rtype, de_fit.msg_level, de_fit.mc_nthreads);

    ncm_mset_catalog_param_pdf (mcbs->mcat, 0);
    ncm_fit_log_covar (fit);
    {
      gdouble p_value = ncm_mset_catalog_param_pdf_pvalue (mcbs->mcat, m2lnL, FALSE);
      ncm_message ("#   - pvalue for fitted model [% 20.15g] %04.2f%%.\n#\n", m2lnL, 100.0 * p_value);
    }

    ncm_mset_catalog_clear (&mcat);
    mcat = ncm_fit_mcbs_get_catalog (mcbs);
    ncm_fit_mcbs_clear (&mcbs);
  }
  
  if (de_fit.mcmc)
  {
    NcmMSetTransKernGauss *mcsg = ncm_mset_trans_kern_gauss_new (0);
    NcmFitMCMC *mcmc = ncm_fit_mcmc_new (fit, NCM_MSET_TRANS_KERN (mcsg), de_fit.msg_level);

    if (de_fit.fisher)
    {
      NcmMatrix *covar = ncm_matrix_dup (fit->fstate->covar);
      ncm_matrix_scale (covar, 2.0);
      ncm_mset_trans_kern_gauss_set_cov (mcsg, covar);
      ncm_matrix_free (covar);
    }
    else
    {
      ncm_mset_trans_kern_gauss_set_cov_from_rescale (mcsg, 0.1);
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
    ncm_fit_mcmc_run_lre (mcmc, de_fit.mc_prerun, de_fit.mc_lre);
    ncm_fit_mcmc_end_run (mcmc);

    ncm_fit_mcmc_mean_covar (mcmc);
    ncm_mset_catalog_param_pdf (mcmc->mcat, 0);
    ncm_fit_log_covar (fit);

    ncm_mset_catalog_clear (&mcat);
    mcat = ncm_fit_mcmc_get_catalog (mcmc);
    ncm_fit_mcmc_clear (&mcmc);    
  }
  
  if (de_fit.esmcmc)
  {
    NcmMSetTransKernGauss *init_sampler = ncm_mset_trans_kern_gauss_new (0);
    NcmFitESMCMC *esmcmc;

    if (de_fit.esmcmc_walk)
    {
      NcmFitESMCMCWalkerWalk *walk = ncm_fit_esmcmc_walker_walk_new (de_fit.mc_nwalkers);
      
      esmcmc = ncm_fit_esmcmc_new_funcs_array (fit, 
                                               de_fit.mc_nwalkers, 
                                               NCM_MSET_TRANS_KERN (init_sampler), 
                                               NCM_FIT_ESMCMC_WALKER (walk),
                                               de_fit.msg_level,
                                               funcs_oa);
      
      ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (walk));
    }
    else if (de_fit.esmcmc_aps)
    {
      NcmFitESMCMCWalkerAPS *aps = ncm_fit_esmcmc_walker_aps_new (de_fit.mc_nwalkers, ncm_mset_fparams_len (mset));
      
      esmcmc = ncm_fit_esmcmc_new_funcs_array (fit, 
                                               de_fit.mc_nwalkers, 
                                               NCM_MSET_TRANS_KERN (init_sampler), 
                                               NCM_FIT_ESMCMC_WALKER (aps),
                                               de_fit.msg_level,
                                               funcs_oa);
      
      ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (aps));
    }
    else
    {
      NcmFitESMCMCWalkerStretch *stretch = ncm_fit_esmcmc_walker_stretch_new (de_fit.mc_nwalkers, ncm_mset_fparams_len (mset));
      esmcmc = ncm_fit_esmcmc_new_funcs_array (fit, 
                                               de_fit.mc_nwalkers, 
                                               NCM_MSET_TRANS_KERN (init_sampler), 
                                               NCM_FIT_ESMCMC_WALKER (stretch),
                                               de_fit.msg_level,
                                               funcs_oa);
      
      if (de_fit.esmcmc_sbox)
      {
        ncm_fit_esmcmc_walker_stretch_set_box_mset (stretch, mset);
      }
      ncm_fit_esmcmc_walker_stretch_multi (stretch, de_fit.esmcmc_ms);

      ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (stretch));
    }

    ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
    ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));

    if (de_fit.mc_nthreads > 1)
      ncm_fit_esmcmc_set_nthreads (esmcmc, de_fit.mc_nthreads);
    
    if (de_fit.fisher)
    {
      NcmMatrix *covar = ncm_matrix_dup (fit->fstate->covar);
      ncm_matrix_scale (covar, 2.0);
      ncm_mset_trans_kern_gauss_set_cov (init_sampler, covar);
      ncm_matrix_free (covar);
    }
    else
    {
      ncm_mset_trans_kern_gauss_set_cov_from_rescale (init_sampler, 0.01);
    }

    if (de_fit.mc_seed > -1)
    {
      NcmRNG *esmcmc_rng = ncm_rng_seeded_new (NULL, de_fit.mc_seed);
      ncm_fit_esmcmc_set_rng (esmcmc, esmcmc_rng);
      ncm_rng_free (esmcmc_rng);
    }

    if (de_fit.mc_data != NULL)
      ncm_fit_esmcmc_set_data_file (esmcmc, de_fit.mc_data);

    ncm_fit_esmcmc_start_run (esmcmc);
    ncm_fit_esmcmc_run_lre (esmcmc, de_fit.mc_prerun, de_fit.mc_lre);
    ncm_fit_esmcmc_end_run (esmcmc);

    ncm_mset_catalog_clear (&mcat);
    mcat = ncm_fit_esmcmc_get_catalog (esmcmc);

		ncm_fit_esmcmc_mean_covar (esmcmc);
    ncm_mset_catalog_param_pdf (mcat, 0);
    ncm_fit_log_covar (fit);

    ncm_fit_esmcmc_clear (&esmcmc);
    ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (init_sampler));
  }

  if (de_fit.onedim_cr != NULL)
  {
    gchar **onedim_cr = de_fit.onedim_cr;
    while (onedim_cr[0] != NULL)
    {
      const gchar *pname = onedim_cr[0];
      const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, onedim_cr[0]);
      onedim_cr = &onedim_cr[1];

      if (pi == NULL)
      {
        g_warning ("darkenergy: parameter `%s' not found.", pname);
        continue;
      }
      else
      {
        NcmLHRatio1d *lhr1d = ncm_lh_ratio1d_new (fit, pi);
        gdouble prob_sigma, err_inf, err_sup;

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

        ncm_lh_ratio1d_find_bounds (lhr1d, prob_sigma, fit->mtype, &err_inf, &err_sup);
        ncm_lh_ratio1d_free (lhr1d);

        ncm_message ("#  One dimension confidence region for %s[%05d:%02d] = % .5g (% .5g, % .5g)\n",
                     ncm_mset_param_name (mset, pi->mid, pi->pid), pi->mid, pi->pid,
                     ncm_mset_param_get (mset, pi->mid, pi->pid), err_inf, err_sup);
      }
    }
  }

  if (de_fit.nsigma >= 0 && (de_fit.bidim_cr[0] != NULL) && (de_fit.bidim_cr[1] != NULL))
  {
    NcmLHRatio2d *lhr2d;
    const NcmMSetPIndex *pi1 = ncm_mset_fparam_get_pi_by_name (mset, de_fit.bidim_cr[0]);
    const NcmMSetPIndex *pi2 = ncm_mset_fparam_get_pi_by_name (mset, de_fit.bidim_cr[1]);
    NcmLHRatio2dRegion *rg_1sigma = NULL;
    NcmLHRatio2dRegion *rg_2sigma = NULL;
    NcmLHRatio2dRegion *rg_3sigma = NULL;

    if (pi1 == NULL)
      g_error ("darkenergy: cannot find parameter named `%s'", de_fit.bidim_cr[0]);
    if (pi2 == NULL)
      g_error ("darkenergy: cannot find parameter named `%s'", de_fit.bidim_cr[1]);
    
    lhr2d = ncm_lh_ratio2d_new (fit, pi1, pi2, de_fit.lhr_prec);
    
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

      f_PL = nc_de_open_dataout_file (cosmo, "PL", &pfile);

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

  if (de_fit.nsigma_fisher >= 0 && (de_fit.bidim_cr[0] != NULL) && (de_fit.bidim_cr[1] != NULL))
  {
    NcmLHRatio2d *lhr2d;
    const NcmMSetPIndex *pi1 = ncm_mset_fparam_get_pi_by_name (mset, de_fit.bidim_cr[0]);
    const NcmMSetPIndex *pi2 = ncm_mset_fparam_get_pi_by_name (mset, de_fit.bidim_cr[1]);
    NcmLHRatio2dRegion *rg_1sigma = NULL;
    NcmLHRatio2dRegion *rg_2sigma = NULL;
    NcmLHRatio2dRegion *rg_3sigma = NULL;

    if (pi1 == NULL)
      g_error ("darkenergy: cannot find parameter named `%s'", de_fit.bidim_cr[0]);
    if (pi2 == NULL)
      g_error ("darkenergy: cannot find parameter named `%s'", de_fit.bidim_cr[1]);

    lhr2d = ncm_lh_ratio2d_new (fit, pi1, pi2, de_fit.lhr_prec);

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

      f_MFcr = nc_de_open_dataout_file (cosmo, "MF_cr", &mcrfile);

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

    f_mf = nc_de_open_dataout_file (cosmo, "MassFunction", &mfile);
    for (i = 0; i < ca_array->len; i++)
    {
      NcDataClusterNCount *dca_unbinned = g_ptr_array_index (ca_array, i);
      nc_data_cluster_ncount_print (dca_unbinned, cosmo, f_mf, full_cmd_line);
      fprintf (f_mf, "\n\n");
    }
    fclose (f_mf);
    g_ptr_array_free (ca_array, TRUE);
    ca_array = NULL;

    ncm_message ("# MassFunction file: %s \n", mfile);

    g_free (mfile);
  }

  if (ca_array != NULL)
  {
    g_ptr_array_free (ca_array, TRUE);
  }

  if (de_fit.save_mset != NULL)
    ncm_mset_save (mset, ser, de_fit.save_mset, TRUE);

  g_free (de_model_entries); de_model_entries = NULL;
  g_free (de_data_simple_entries); de_data_simple_entries = NULL;
  g_free (de_data_cluster_entries); de_data_cluster_entries = NULL;
  g_free (de_fit_entries); de_fit_entries = NULL;

  if (fiduc != NULL)
    ncm_mset_free (fiduc);

  ncm_serialize_clear (&ser);
  ncm_mset_catalog_clear (&mcat);
  ncm_serialize_global_reset (FALSE);
  ncm_model_free (NCM_MODEL (cosmo));
  ncm_mset_free (mset);
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_rng_free (rng);
  nc_distance_free (dist);
  g_free (full_cmd_line);

  if (TRUE)
  {
    g_clear_pointer (&de_model.model_name, g_free);

    g_clear_pointer (&de_data_simple.snia_id,     g_free);
    g_clear_pointer (&de_data_simple.snia_objser, g_free);
    g_clear_pointer (&de_data_simple.cmb_id,      g_free);
    g_clear_pointer (&de_data_simple.cluster_id,  g_free);

    g_clear_pointer (&de_data_simple.bao_id,       g_strfreev);
    g_clear_pointer (&de_data_simple.H_id,         g_strfreev);
    g_clear_pointer (&de_data_simple.H_BAO_id,     g_strfreev);
    g_clear_pointer (&de_data_simple.priors_gauss, g_strfreev);
    g_clear_pointer (&de_data_simple.data_files,   g_strfreev);

    g_clear_pointer (&de_data_cluster.filter_type,       g_free);
    g_clear_pointer (&de_data_cluster.ps_type,           g_free);
    g_clear_pointer (&de_data_cluster.multiplicity_name, g_free);
    g_clear_pointer (&de_data_cluster.clusterm_ser,      g_free);
    g_clear_pointer (&de_data_cluster.clusterz_ser,      g_free);
    g_clear_pointer (&de_data_cluster.save_cata,         g_free);

    g_clear_pointer (&de_data_cluster.cata_file, g_strfreev);

    g_clear_pointer (&de_fit.file_out,    g_free);
    g_clear_pointer (&de_fit.fit_type,    g_free);
    g_clear_pointer (&de_fit.fit_diff,    g_free);
    g_clear_pointer (&de_fit.fit_algo,    g_free);
    g_clear_pointer (&de_fit.bidim_cr[0], g_free);
    g_clear_pointer (&de_fit.bidim_cr[1], g_free);
    g_clear_pointer (&de_fit.fiducial,    g_free);
    g_clear_pointer (&de_fit.mc_data,     g_free);
    g_clear_pointer (&de_fit.save_mset,   g_free);
    
    g_clear_pointer (&de_fit.onedim_cr, g_strfreev);
    g_clear_pointer (&de_fit.funcs, g_strfreev);
  }
  
  return 0;
}
