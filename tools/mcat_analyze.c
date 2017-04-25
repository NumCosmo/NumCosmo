/***************************************************************************
 *            mcat_analyze.c
 *
 *  Mon March 16 11:04:23 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti & Mariana Penna Lima (January 4th 2016)
 *  <sandro@iaoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * mcat_analyze.c
 *
 * Copyright (C) 2015 - Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <gsl/gsl_cdf.h>

typedef enum _NcMCatFunc {
  NC_MCAT_FUNC_HICOSMO = 0,
  NC_MCAT_FUNC_DIST,
} NcMCatFunc;

void
_nc_bestfit_error (NcmStatsDist1d *sd1, gdouble Pa, gdouble center, gchar *center_desc, gdouble *lb, gdouble *ub)
{
  gdouble P_bf, Pa_2;
  P_bf = ncm_stats_dist1d_eval_pdf (sd1, center);
  Pa_2 = Pa * 0.5;
  
  if (P_bf < Pa_2)
  {
    ncm_message ("# Right-skewed distribution - %s close to the minimum value: maximum left tail correspond to %6.2f%% confidence interval (CI), where the required probability (right and left) was %6.2f%% CI.\n", 
                 center_desc, P_bf * 100.0, Pa_2 * 100.0);
    ncm_message ("# The lower and upper error bounds correspond to %6.2f%% and %6.2f%% CI, respectively. Total = %6.2f%% (%5.2f-sigma).\n", 
                 P_bf * 100.0, (Pa - P_bf) * 100.0, Pa * 100.0, sqrt (gsl_cdf_chisq_Pinv (Pa, 1.0)));
    
    lb[0] = sd1->xi;
    ub[0] = ncm_stats_dist1d_eval_inv_pdf (sd1, Pa);
  }
  else if (P_bf > 1.0 - Pa_2)
  {
    ncm_message ("# Left-skewed distribution - %s close to the maximum value: maximum right tail correspond to %6.2f%% confidence interval (CI), where the required probability (right and left) was %6.2f%% CI.\n", 
                 center_desc, (1.0 - P_bf) * 100.0, Pa_2 * 100.0);
    ncm_message ("# The lower and upper error bounds correspond to %6.2f%% and %6.2f%% CI, respectively. Total = %6.2f%% (%5.2f-sigma).\n", 
                 (Pa - 1.0 + P_bf) * 100.0, (1.0 - P_bf) * 100.0, Pa * 100.0, sqrt (gsl_cdf_chisq_Pinv (Pa, 1.0)));
    lb[0] = ncm_stats_dist1d_eval_inv_pdf_tail (sd1, Pa);
    ub[0] = sd1->xf;
  }
  else
  {
    lb[0] = ncm_stats_dist1d_eval_inv_pdf (sd1, P_bf - Pa_2);
    ub[0] = ncm_stats_dist1d_eval_inv_pdf (sd1, P_bf + Pa_2);
  }

  lb[0] = center - lb[0];
  ub[0] = ub[0] - center;
  
}

gint
main (gint argc, gchar *argv[])
{
  gchar *cat_filename     = NULL;
  gboolean auto_trim      = FALSE;
  gboolean info           = FALSE;
  gboolean info_scf       = FALSE;
  gboolean diag_chains    = FALSE;
  gboolean chain_evol     = FALSE;
  gboolean list_all       = FALSE;
  gboolean list_hicosmo   = FALSE;
  gboolean list_hicosmo_z = FALSE;
  gboolean list_dist      = FALSE;
  gboolean list_dist_z    = FALSE;
  gboolean use_direct     = FALSE;
  gboolean dump           = FALSE;
  gint dump_chain         = -1;
  gint trim               = -1;
  gdouble zi = 0.0;
  gdouble zf = 1.0;
  gint nsteps = 100;
  gint burnin = 0;
  gint ntests = 100;
  gchar **funcs          = NULL;
  gchar **distribs       = NULL;
  gchar **params         = NULL;
  gchar **params_evol    = NULL;
  gchar **mode_errors    = NULL;
  gchar **median_errors  = NULL;
  gchar **bestfit_errors = NULL;
  gchar **funcs_pvalue   = NULL;
  gchar **visual_hw      = NULL;
  guint i;
  
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry entries[] =
  {
    { "catalog",        'c', 0, G_OPTION_ARG_FILENAME,     &cat_filename,   "Catalog filename.", NULL },
    { "auto-trim",      'a', 0, G_OPTION_ARG_NONE,         &auto_trim,      "Auto trim the catalog.", NULL },
    { "info",           'i', 0, G_OPTION_ARG_NONE,         &info,           "Print catalog information.", NULL },
    { "info-scor-full", 'C', 0, G_OPTION_ARG_NONE,         &info_scf,       "Calculates the selfcorrelation time using the full sample not the ensemble averages.", NULL },
    { "diag-chains",      0, 0, G_OPTION_ARG_NONE,         &diag_chains,    "Applies diagnostics to all chains.", NULL },
    { "ntests",         'n', 0, G_OPTION_ARG_INT,          &ntests,         "Number of tests to use in the diagnostics.", NULL },
    { "chain-evol",     'I', 0, G_OPTION_ARG_NONE,         &chain_evol,     "Print chain evolution.", NULL },
    { "list",           'l', 0, G_OPTION_ARG_NONE,         &list_all,       "Print all available functions.", NULL },
    { "list-hicosmo",     0, 0, G_OPTION_ARG_NONE,         &list_hicosmo,   "Print available constant functions from NcHICosmo.", NULL },
    { "list-hicosmo-z",   0, 0, G_OPTION_ARG_NONE,         &list_hicosmo_z, "Print available redshift functions from NcHICosmo.", NULL },
    { "list-dist",        0, 0, G_OPTION_ARG_NONE,         &list_dist,      "Print available constant functions from NcDistance.", NULL },
    { "list-dist-z",      0, 0, G_OPTION_ARG_NONE,         &list_dist_z,    "Print available redshift functions from NcDistance.", NULL },
    { "use-direct",       0, 0, G_OPTION_ARG_NONE,         &use_direct,     "Whether to use the direct quantile algorithm (much more memory intensive but faster for small samples).", NULL },
    { "burn-in",        'b', 0, G_OPTION_ARG_INT,          &burnin,         "Burn-in size (default 0).", NULL },
    { "zi",               0, 0, G_OPTION_ARG_DOUBLE,       &zi,             "Initial redshift (default 0).", NULL },
    { "zf",               0, 0, G_OPTION_ARG_DOUBLE,       &zf,             "Final redshift (default 1).", NULL },
    { "nsteps",           0, 0, G_OPTION_ARG_INT,          &nsteps,         "Number of points in the functions grid (default 100).", NULL },
    { "function",       'f', 0, G_OPTION_ARG_STRING_ARRAY, &funcs,          "Redshift functions to be analyzed.", NULL},
    { "distribution",   'd', 0, G_OPTION_ARG_STRING_ARRAY, &distribs,       "Function distributions to be analyzed.", NULL},
    { "parameter",      'p', 0, G_OPTION_ARG_STRING_ARRAY, &params,         "Model parameters' to be analyzed.", NULL},
    { "parameter-evol", 'P', 0, G_OPTION_ARG_STRING_ARRAY, &params_evol,    "Calculate the time evolution of the parameter.", NULL},
    { "mode-error",     'o', 0, G_OPTION_ARG_STRING_ARRAY, &mode_errors,    "Print mode and 1-3 sigma asymmetric error bars of the model parameters' to be analyzed.", NULL},
    { "median-error",   'e', 0, G_OPTION_ARG_STRING_ARRAY, &median_errors,  "Print median and 1-3 sigma asymmetric error bars of the model parameters' to be analyzed.", NULL},
    { "bestfit-error",  's', 0, G_OPTION_ARG_STRING_ARRAY, &bestfit_errors, "Print best fit and 1-3 sigma asymmetric error bars of the model parameters' to be analyzed.", NULL},
    { "funcs-pvalue",   'F', 0, G_OPTION_ARG_STRING_ARRAY, &funcs_pvalue,   "Print the p-value of the function at each redshift, giving the upper integration limits.", NULL },
    { "visual-hw",      'V', 0, G_OPTION_ARG_STRING_ARRAY, &visual_hw,      "Print the points to the visual HW test.", NULL },
    { "dump",           'D', 0, G_OPTION_ARG_NONE,         &dump,           "Print all chains interweaved.", NULL },
    { "dump-chain",       0, 0, G_OPTION_ARG_INT,          &dump_chain,     "Print all points from the N-th chain.", "N"},
    { "trim",           't', 0, G_OPTION_ARG_INT,          &trim,           "Trim the catalog at T.", "T" },
    { NULL }
  };

  ncm_cfg_init ();
  
  context = g_option_context_new ("- analyze catalogs from Monte Carlo (MC, MCMC, ESMCMC, bootstrap MC).");
  g_option_context_set_summary (context, "catalog analyzer");
  g_option_context_set_description (context, "CA Description <FIXME>");

  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
  {
    g_print ("option parsing failed: %s\n", error->message);
    exit (1);
  }

  if (list_hicosmo || list_all)
  {
    GArray *func_table = ncm_mset_func_list_select ("NcHICosmo", 0, 1); 
    ncm_cfg_msg_sepa ();
    g_message ("# Available constant functions from NcHICosmo models:\n");    
    for (i = 0; i < func_table->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (func_table, NcmMSetFuncListStruct, i);
      g_message ("# - %32s: %s\n", fdata->name, fdata->desc);
    }
    g_array_unref (func_table);
  }
  if (list_hicosmo_z || list_all)
  {
    GArray *func_z_table = ncm_mset_func_list_select ("NcHICosmo", 1, 1); 
    ncm_cfg_msg_sepa ();
    g_message ("# Available redshift functions from NcHICosmo models:\n");    
    for (i = 0; i < func_z_table->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (func_z_table, NcmMSetFuncListStruct, i);
      g_message ("# - %32s: %s\n", fdata->name, fdata->desc);
    }
    g_array_unref (func_z_table);
  }
  if (list_dist || list_all)
  {
    GArray *func_table = ncm_mset_func_list_select ("NcDistance", 0, 1); 
    ncm_cfg_msg_sepa ();
    g_message ("# Available constant functions from NcDistance:\n");    
    for (i = 0; i < func_table->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (func_table, NcmMSetFuncListStruct, i);
      g_message ("# - %32s: %s\n", fdata->name, fdata->desc);
    }
    g_array_unref (func_table);
  }
  if (list_dist_z || list_all)
  {
    GArray *func_z_table = ncm_mset_func_list_select ("NcDistance", 1, 1);
    ncm_cfg_msg_sepa ();
    g_message ("# Available redshift functions from NcDistance:\n");    
    for (i = 0; i < func_z_table->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (func_z_table, NcmMSetFuncListStruct, i);
      g_message ("# - %32s: %s\n", fdata->name, fdata->desc);
    }
    g_array_unref (func_z_table);
  }
  
  if (cat_filename == NULL)
  {
    if (list_all || list_hicosmo || list_hicosmo_z || list_dist || list_dist_z)
    {
      exit (0);
    }
    else
    {
      g_print ("No catalog filename, use --catalog/-c.\n");
      exit (1);
    }
  }
  else
  {
    NcmMSetCatalog *mcat = ncm_mset_catalog_new_from_file_ro (cat_filename, burnin);
    NcmMSet *mset = ncm_mset_catalog_get_mset (mcat);

    if (auto_trim)
    {
      ncm_mset_catalog_trim_by_type (mcat, ntests, NCM_MSET_CATALOG_TRIM_TYPE_ESS, NCM_FIT_RUN_MSGS_FULL);
      if (trim > 0)
        g_warning ("mcat_analize: --auto-trim enabled ignoring --trim");
    }
    else if (trim > 0)
    {
      ncm_mset_catalog_trim (mcat, trim);
    }

    ncm_mset_catalog_estimate_autocorrelation_tau (mcat, info_scf);
    if (info)
    {
      /*NcmMatrix *cov = NULL;*/
      const gchar *rtype_str = ncm_mset_catalog_get_run_type (mcat);
      ncm_cfg_msg_sepa ();
      g_message ("# Catalog run type: `%s'.\n", rtype_str);
      g_message ("# Catalog size:      %u.\n", ncm_mset_catalog_len (mcat));
      g_message ("# Catalog n-chains:  %u.\n", mcat->nchains);
      g_message ("# Catalog nadd-vals: %u.\n", mcat->nadd_vals);
      g_message ("# Catalog weighted:  %s.\n", mcat->weighted ? "yes" : "no");
      ncm_mset_catalog_log_current_chain_stats (mcat);

      {
        gdouble max_ess = 0.0;
        ncm_mset_catalog_calc_max_ess_time (mcat, ntests, &max_ess, NCM_FIT_RUN_MSGS_FULL);
        ncm_mset_catalog_calc_heidel_diag (mcat, ntests, 0.0, NCM_FIT_RUN_MSGS_FULL);
      }
      
      g_message ("#\n");
      
      ncm_mset_pretty_log (mset);
      /*ncm_mset_catalog_get_covar (mcat, &cov);*/
      /*ncm_mset_fparams_log_covar (mset, cov);*/
      ncm_mset_catalog_log_full_covar (mcat);
      /*ncm_matrix_clear (&cov);*/
      ncm_mset_catalog_log_current_stats (mcat);
    }

    if (diag_chains)
    {
      gdouble max_ess   = 0.0;
      gdouble wp_pvalue = 0.0;
      
      ncm_mset_catalog_max_ess_time_by_chain (mcat, ntests, &max_ess, NCM_FIT_RUN_MSGS_FULL);
      ncm_mset_catalog_heidel_diag_by_chain (mcat, ntests, 0.0, &wp_pvalue, NCM_FIT_RUN_MSGS_FULL);
    }

    if (chain_evol)
    {
      guint last_t = ncm_mset_catalog_max_time (mcat);
      const gdouble nchains = mcat->nchains;
      guint i;
      
      ncm_message ("# Chain evolution from 0 to %u\n", last_t);
      for (i = 0; i < last_t; i++)
      {
        NcmVector *e_mean = ncm_mset_catalog_peek_e_mean_t (mcat, i);
        NcmVector *e_var  = ncm_mset_catalog_peek_e_var_t (mcat, i);
        const guint len = ncm_vector_len (e_mean);
        guint j;
        
        ncm_message ("%10u", i);
        for (j = 0; j < len; j++)
        {
          const gdouble mean_j = ncm_stats_vec_get_mean (mcat->pstats, j);
          const gdouble var_j  = ncm_stats_vec_get_var (mcat->pstats, j);
          const gdouble sd_j   = sqrt (var_j / nchains);
          
          ncm_message (" % 20.15g % 20.15g", (ncm_vector_get (e_mean, j) - mean_j) / sd_j, sqrt (ncm_vector_get (e_var, j) / nchains) / sd_j);
        }
        ncm_message ("\n");
      }
    }
    
    /*********************************************************************************************************
     * 
     * Functions
     * 
     *********************************************************************************************************/
    if (funcs != NULL)
    {
      guint nfuncs = g_strv_length (funcs);
      NcmVector *z_vec = ncm_vector_new (nsteps);
      GArray *p_val = g_array_new (FALSE, FALSE, sizeof (gdouble));
      NcDistance *dist = nc_distance_new (zf + 0.2);
      
      g_array_set_size (p_val, 3);
      g_array_index (p_val, gdouble, 0) = gsl_cdf_chisq_P (1.0, 1.0);
      g_array_index (p_val, gdouble, 1) = gsl_cdf_chisq_P (4.0, 1.0);
      g_array_index (p_val, gdouble, 2) = gsl_cdf_chisq_P (9.0, 1.0);
      
      for (i = 0; i < nsteps; i++)
      {
        const gdouble z = zi + (zf - zi) / (nsteps - 1.0) * i;
        ncm_vector_set (z_vec, i, z);
      }

      for (i = 0; i < nfuncs; i++)
      {
        NcmMSetFunc *mset_func = NULL;
        NcmMatrix *res = NULL;
        guint k;

        if (mset_func == NULL)
        {
          if (ncm_mset_func_list_has_ns_name ("NcHICosmo", funcs[i]))
          {
            mset_func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcHICosmo", funcs[i], NULL));
            if (ncm_mset_func_get_dim (mset_func) != 1 || ncm_mset_func_get_nvar (mset_func) != 1)
            {
              g_warning ("# Function `%s' is not R => R, skipping.", ncm_mset_func_peek_name (mset_func));
              ncm_mset_func_clear (&mset_func);
              continue;
            }
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcHICosmo z function: `%s' in [% 20.15g % 20.15g].\n", ncm_mset_func_peek_desc (mset_func), zi, zf);
          }
        }
        if (mset_func == NULL)
        {
          if (ncm_mset_func_list_has_ns_name ("NcDistance", funcs[i]))
          {
            mset_func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcDistance", funcs[i], G_OBJECT (dist)));
            if (ncm_mset_func_get_dim (mset_func) != 1 || ncm_mset_func_get_nvar (mset_func) != 1)
            {
              g_warning ("# Function `%s' is not R => R, skipping.", ncm_mset_func_peek_name (mset_func));
              ncm_mset_func_clear (&mset_func);
              continue;
            }
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcDistance z function: `%s' in [% 20.15g % 20.15g].\n", ncm_mset_func_peek_desc (mset_func), zi, zf);
          }
        }
        if (mset_func == NULL)
        {
          g_warning ("# Function `%s' not found, skipping...\n", funcs[i]);
          continue;
        }

        if (use_direct)
          res = ncm_mset_catalog_calc_ci_direct (mcat, mset_func, z_vec, p_val);
        else
          res = ncm_mset_catalog_calc_ci_interp (mcat, mset_func, z_vec, p_val, 100, NCM_FIT_RUN_MSGS_SIMPLE);

        for (k = 0; k < nsteps; k++)
        {
          ncm_message ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", 
                       ncm_vector_get (z_vec, k), 
                       ncm_matrix_get (res, k, 0), 
                       ncm_matrix_get (res, k, 1), ncm_matrix_get (res, k, 2), 
                       ncm_matrix_get (res, k, 3), ncm_matrix_get (res, k, 4), 
                       ncm_matrix_get (res, k, 5), ncm_matrix_get (res, k, 6));
        }
        ncm_message ("\n\n");

        ncm_mset_func_free (mset_func);
        ncm_matrix_free (res);
      }
      g_array_unref (p_val);
      ncm_vector_clear (&z_vec);
      nc_distance_clear (&dist);
    }
    /*********************************************************************************************************
     * 
     * Distributions
     * 
     *********************************************************************************************************/
    if (distribs != NULL)
    {
      guint ndistribs = g_strv_length (distribs);
      NcDistance *dist = nc_distance_new (zf + 0.2);

      for (i = 0; i < ndistribs; i++)
      {
        NcmMSetFunc *mset_func = NULL;
        gdouble *x = NULL;
        guint len;
        guint k;
        gchar *func_name = ncm_util_function_params (distribs[i], &x, &len);
        g_assert (func_name != NULL);

        if (mset_func == NULL)
        {
          if (ncm_mset_func_list_has_ns_name ("NcHICosmo", func_name))
          {
            mset_func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcHICosmo", func_name, NULL));
            if (len > 0)
              ncm_mset_func_set_eval_x (mset_func, x, len);
            
            if (ncm_mset_func_get_dim (mset_func) != 1 || !ncm_mset_func_is_const (mset_func))
            {
              g_warning ("# Function `%s' is not constant, skipping.", ncm_mset_func_peek_name (mset_func));
              ncm_mset_func_clear (&mset_func);
              continue;
            }
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcHICosmo function: `%s'.\n", ncm_mset_func_peek_desc (mset_func));
          }
        }
        if (mset_func == NULL)
        {
          if (ncm_mset_func_list_has_ns_name ("NcDistance", func_name))
          {
            mset_func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcDistance", func_name, G_OBJECT (dist)));
            if (len > 0)
              ncm_mset_func_set_eval_x (mset_func, x, len);

            if (ncm_mset_func_get_dim (mset_func) != 1 || !ncm_mset_func_is_const (mset_func))
            {
              g_warning ("# Function `%s' is not constant, skipping.", ncm_mset_func_peek_name (mset_func));
              ncm_mset_func_clear (&mset_func);
              continue;
            }
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcDistance function: `%s'.\n", ncm_mset_func_peek_desc (mset_func));
          }
        }
        if (mset_func == NULL)
        {
          g_warning ("# Function `%s' not found, skipping...\n", func_name);
          continue;
        }

        {
          NcmStatsDist1d *sd1 = ncm_mset_catalog_calc_distrib (mcat, mset_func, NCM_FIT_RUN_MSGS_SIMPLE);
          for (k = 0; k < nsteps; k++)
          {
            gdouble x = sd1->xi + (sd1->xf - sd1->xi) / (nsteps - 1.0) * k;
            gdouble u = 1.0 / (nsteps - 1.0) * k;
            ncm_message ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", 
                         x, 
                         ncm_stats_dist1d_eval_p (sd1, x),
                         ncm_stats_dist1d_eval_pdf (sd1, x),
                         u, 
                         ncm_stats_dist1d_eval_inv_pdf (sd1, u),
                         ncm_stats_dist1d_eval_inv_pdf_tail (sd1, u));
          }
          ncm_stats_dist1d_free (sd1);
        }
        ncm_message ("\n\n");

        g_free (func_name);
        ncm_mset_func_free (mset_func);
      }
      nc_distance_clear (&dist);
    }

    /*********************************************************************************************************
     * 
     * Parameters
     * 
     *********************************************************************************************************/
    if (params != NULL)
    {
      guint nparams = g_strv_length (params);

      for (i = 0; i < nparams; i++)
      {
        const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, params[i]);
        gchar *end_ptr = NULL;
        glong add_param = strtol (params[i], &end_ptr, 10);
        guint k;

        if (pi == NULL && (params[i] == end_ptr))
        {
          g_warning ("# Parameter `%s' not found, skipping...\n", params[i]);
          continue;
        }

        {
          NcmStatsDist1d *sd1;
          if (pi != NULL)
            sd1 = ncm_mset_catalog_calc_param_distrib (mcat, pi, NCM_FIT_RUN_MSGS_SIMPLE);
          else
            sd1 = ncm_mset_catalog_calc_add_param_distrib (mcat, add_param, NCM_FIT_RUN_MSGS_SIMPLE);
          
          for (k = 0; k < nsteps; k++)
          {
            gdouble x = sd1->xi + (sd1->xf - sd1->xi) / (nsteps - 1.0) * k;
            gdouble u = 1.0 / (nsteps - 1.0) * k;
            ncm_message ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", 
                         x, 
                         ncm_stats_dist1d_eval_p (sd1, x),
                         ncm_stats_dist1d_eval_pdf (sd1, x),
                         u, 
                         ncm_stats_dist1d_eval_inv_pdf (sd1, u),
                         ncm_stats_dist1d_eval_inv_pdf_tail (sd1, u));
          }
          ncm_stats_dist1d_free (sd1);
        }
        ncm_message ("\n\n");
      }
    }
    
    /*********************************************************************************************************
     * 
     * Parameters evolution
     * 
     *********************************************************************************************************/
    if (params_evol != NULL)
    {
      guint nparams = g_strv_length (params_evol);

      for (i = 0; i < nparams; i++)
      {
        const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, params_evol[i]);
        gchar *end_ptr = NULL;
        glong add_param = strtol (params_evol[i], &end_ptr, 10);
        guint t, k;

        if (pi == NULL && (params_evol[i] == end_ptr))
        {
          g_warning ("# Parameter `%s' not found, skipping...\n", params_evol[i]);
          continue;
        }

        {
          NcmVector *pv;
          NcmMatrix *evol_t;
          gchar *out_file = g_strdup_printf ("param_evol_%s.dat", params_evol[i]);
          FILE *out = fopen (out_file, "w");
          
          if (pi != NULL)
            ncm_mset_catalog_calc_param_ensemble_evol (mcat, pi, nsteps, NCM_FIT_RUN_MSGS_SIMPLE, &pv, &evol_t);
          else
            ncm_mset_catalog_calc_add_param_ensemble_evol (mcat, add_param, nsteps, NCM_FIT_RUN_MSGS_SIMPLE, &pv, &evol_t); 
          
          for (t = 0; t < ncm_matrix_nrows (evol_t); t++)
          {
            for (k = 0; k < ncm_matrix_ncols (evol_t); k++)
            {
              fprintf (out, "%u % 20.15g % 20.15g\n", t, ncm_vector_get (pv, k), ncm_matrix_get (evol_t, t, k));
            }
            fprintf (out, "\n\n");
          }

          fclose (out);
          ncm_vector_free (pv);
          ncm_matrix_free (evol_t);
        }
      }
    }
    
   /*********************************************************************************************************
     * 
     * Parameters - Mode and error bars
     * 
     *********************************************************************************************************/
    if (mode_errors != NULL)
    {
      guint nparams = g_strv_length (mode_errors);
      gdouble Pa1 = gsl_cdf_chisq_P (1.0, 1.0);
      gdouble Pa2 = gsl_cdf_chisq_P (4.0, 1.0);
      gdouble Pa3 = gsl_cdf_chisq_P (9.0, 1.0);
      ncm_message ("# Computing mode and 1-3 sigma asymmetric error bars - lower (l) and upper (u) bounds, respectively.\n");

      for (i = 0; i < nparams; i++)
      {
        const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, mode_errors[i]);
        gchar *end_ptr = NULL;
        glong add_param = strtol (mode_errors[i], &end_ptr, 10);
        ncm_message ("# Parameter `%s'| mode | 1l 1u | 2l 2u | 3l 3u\n", mode_errors[i]);

        if (pi == NULL && (mode_errors[i] == end_ptr))
        {
          g_warning ("# Parameter `%s' not found, skipping...\n", mode_errors[i]);
          continue;
        }

        {
          NcmStatsDist1d *sd1;
          if (pi != NULL)
            sd1 = ncm_mset_catalog_calc_param_distrib (mcat, pi, NCM_FIT_RUN_MSGS_SIMPLE);
          else
            sd1 = ncm_mset_catalog_calc_add_param_distrib (mcat, add_param, NCM_FIT_RUN_MSGS_SIMPLE);

          {
            gdouble mode = ncm_stats_dist1d_eval_mode (sd1);
            gdouble x_l1 = 0.0, x_u1 = 0.0;
            gdouble x_l2 = 0.0, x_u2 = 0.0;
            gdouble x_l3 = 0.0, x_u3 = 0.0;
                            
            _nc_bestfit_error (sd1, Pa1, mode, "mode", &x_l1, &x_u1);
            _nc_bestfit_error (sd1, Pa2, mode, "mode", &x_l2, &x_u2);
            _nc_bestfit_error (sd1, Pa3, mode, "mode", &x_l3, &x_u3);
            
              ncm_message (" % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n",
                           mode, x_l1, x_u1, x_l2, x_u2, x_l3, x_u3);
        }

          ncm_stats_dist1d_free (sd1);
        }
        ncm_message ("\n\n");
      }
    }

    /*********************************************************************************************************
     * 
     * Parameters - Median and error bars
     * 
     *********************************************************************************************************/
    if (median_errors != NULL)
    {
      guint nparams = g_strv_length (median_errors);
      gdouble v1 = (1.0 - gsl_cdf_chisq_P (1.0, 1.0)) * 0.5;
      gdouble v2 = (1.0 - gsl_cdf_chisq_P (4.0, 1.0)) * 0.5;
      gdouble v3 = (1.0 - gsl_cdf_chisq_P (9.0, 1.0)) * 0.5;

      ncm_message ("# Computing median and 1-3 sigma asymmetric error bars - lower (l) and upper (u) bounds, respectively.\n");

      for (i = 0; i < nparams; i++)
      {
        const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, median_errors[i]);
        gchar *end_ptr = NULL;
        glong add_param = strtol (median_errors[i], &end_ptr, 10);
        ncm_message ("# Parameter `%s' | 1l 1u | 2l 2u | 3l 3u\n", median_errors[i]);

        if (pi == NULL && (median_errors[i] == end_ptr))
        {
          g_warning ("# Parameter `%s' not found, skipping...\n", median_errors[i]);
          continue;
        }

        {
          NcmStatsDist1d *sd1;
          if (pi != NULL)
            sd1 = ncm_mset_catalog_calc_param_distrib (mcat, pi, NCM_FIT_RUN_MSGS_SIMPLE);
          else
            sd1 = ncm_mset_catalog_calc_add_param_distrib (mcat, add_param, NCM_FIT_RUN_MSGS_SIMPLE);

          {
              const gdouble median   = ncm_stats_dist1d_eval_inv_pdf (sd1, 0.5);
                            
              const gdouble l1 = ncm_stats_dist1d_eval_inv_pdf (sd1, v1);
              const gdouble u1 = ncm_stats_dist1d_eval_inv_pdf_tail (sd1, v1);
              const gdouble l2 = ncm_stats_dist1d_eval_inv_pdf (sd1, v2);
              const gdouble u2 = ncm_stats_dist1d_eval_inv_pdf_tail (sd1, v2);
              const gdouble l3 = ncm_stats_dist1d_eval_inv_pdf (sd1, v3);
              const gdouble u3 = ncm_stats_dist1d_eval_inv_pdf_tail (sd1, v3);
              
              const gdouble x_l1   = median - l1;
              const gdouble x_u1   = - median + u1;
              const gdouble x_l2   = median - l2;
              const gdouble x_u2   = - median + u2;
              const gdouble x_l3   = median - l3;
              const gdouble x_u3   = - median + u3;

              ncm_message (" % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n",
                           median, x_l1, x_u1, x_l2, x_u2, x_l3, x_u3);
            }

          ncm_stats_dist1d_free (sd1);
        }
        ncm_message ("\n\n");
      }
    }

    /*********************************************************************************************************
     * 
     * Parameters - best fit and error bars
     * 
     *********************************************************************************************************/
    if (bestfit_errors != NULL)
    {
      guint nparams = g_strv_length (bestfit_errors);
      gdouble Pa1 = gsl_cdf_chisq_P (1.0, 1.0);
      gdouble Pa2 = gsl_cdf_chisq_P (4.0, 1.0);
      gdouble Pa3 = gsl_cdf_chisq_P (9.0, 1.0);
      ncm_message ("# Best-fit and 1-3 sigma asymmetric error bars - lower (l) and upper (u) bounds, respectively.\n");

      for (i = 0; i < nparams; i++)
      {
        gchar **name_val = g_strsplit (bestfit_errors[i], "=", 2);
        guint n_name_val = g_strv_length (name_val);
        if (n_name_val != 2)
        {
          g_error ("bestfit_errors: String not valid");
        }
        else
        {
          gchar *name = g_strstrip (name_val[0]);
          gchar *val  = g_strstrip (name_val[1]);
          gchar *endptr = NULL;
          gdouble bestfit = g_ascii_strtod (val, &endptr);
          if (endptr == val)
          {
            g_error ("Convertion to double failed for string `%s'", val);
          }
          else
          {
            const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, name);
            gchar *end_ptr = NULL;
            glong add_param = strtol (name, &end_ptr, 10);
            ncm_message ("# Parameter `%s'| best fit | 1l 1u | 2l 2u | 3l 3u\n", name);

            if (pi == NULL && (bestfit_errors[i] == end_ptr))
            {
              g_warning ("# Parameter `%s' not found, skipping...\n", name);
              continue;
            }

            {
              NcmStatsDist1d *sd1;
              if (pi != NULL)
                sd1 = ncm_mset_catalog_calc_param_distrib (mcat, pi, NCM_FIT_RUN_MSGS_SIMPLE);
              else
                sd1 = ncm_mset_catalog_calc_add_param_distrib (mcat, add_param, NCM_FIT_RUN_MSGS_SIMPLE);

              {
                  gdouble x_l1 = 0.0, x_u1 = 0.0;
                  gdouble x_l2 = 0.0, x_u2 = 0.0;
                  gdouble x_l3 = 0.0, x_u3 = 0.0;

                  _nc_bestfit_error (sd1, Pa1, bestfit, "best fit", &x_l1, &x_u1);
                  _nc_bestfit_error (sd1, Pa2, bestfit, "best fit", &x_l2, &x_u2);
                  _nc_bestfit_error (sd1, Pa3, bestfit, "best fit", &x_l3, &x_u3);

                  ncm_message (" % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n",
                               bestfit, x_l1, x_u1, x_l2, x_u2, x_l3, x_u3);
                }

              ncm_stats_dist1d_free (sd1);
            }
            ncm_message ("\n\n");
          }
        }
      }
    }

    /*********************************************************************************************************
    * 
    * P-values of the Functions at each readshift 
    * 
    *********************************************************************************************************/
    if (funcs_pvalue != NULL)
    {
      guint nfuncs_pv = g_strv_length (funcs_pvalue);
      NcmVector *z_vec = ncm_vector_new (nsteps);
      NcDistance *dist = nc_distance_new (zf + 0.2);

      for (i = 0; i < nsteps; i++)
      {
        const gdouble z = zi + (zf - zi) / (nsteps - 1.0) * i;
        ncm_vector_set (z_vec, i, z);
      }
      
      for (i = 0; i < nfuncs_pv; i++)
      {
        NcmMSetFunc *mset_func = NULL;
        NcmMatrix *res = NULL;
        gdouble *x = NULL;
        guint len;
        guint j,k;
        gchar *func_name = ncm_util_function_params (funcs_pvalue[i], &x, &len);
        GArray *lims = g_array_new (FALSE, FALSE, sizeof (gdouble));

        g_assert (func_name != NULL);
        g_assert_cmpuint (len, !=, 0);

        g_array_append_vals (lims, x, len);
        
        if (mset_func == NULL)
        {
          if (ncm_mset_func_list_has_ns_name ("NcHICosmo", func_name))
          {
            mset_func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcHICosmo", func_name, NULL));
            if (ncm_mset_func_get_dim (mset_func) != 1 || ncm_mset_func_get_nvar (mset_func) != 1)
            {
              g_warning ("# Function `%s' is not R => R, skipping.", ncm_mset_func_peek_name (mset_func));
              ncm_mset_func_clear (&mset_func);
              continue;
            }
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcHICosmo z function: `%s' in [% 20.15g % 20.15g].\n", ncm_mset_func_peek_desc (mset_func), zi, zf);
          }
        }
        if (mset_func == NULL)
        {
          if (ncm_mset_func_list_has_ns_name ("NcDistance", func_name))
          {
            mset_func = NCM_MSET_FUNC (ncm_mset_func_list_new_ns_name ("NcDistance", func_name, G_OBJECT (dist)));
            if (ncm_mset_func_get_dim (mset_func) != 1 || ncm_mset_func_get_nvar (mset_func) != 1)
            {
              g_warning ("# Function `%s' is not R => R, skipping.", ncm_mset_func_peek_name (mset_func));
              ncm_mset_func_clear (&mset_func);
              continue;
            }
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcDistance z function: `%s' in [% 20.15g % 20.15g].\n", ncm_mset_func_peek_desc (mset_func), zi, zf);
          }
        }       
        if (mset_func == NULL)
        {
          g_warning ("# Function `%s' not found, skipping...\n", func_name);
          continue;
        }        

        res = ncm_mset_catalog_calc_pvalue (mcat, mset_func, z_vec, lims, 100, NCM_FIT_RUN_MSGS_SIMPLE);

        for (k = 0; k < nsteps; k++)
        {
          ncm_message ("% 20.15g", ncm_vector_get (z_vec, k)); 
          for (j = 0; j < len; j++)
          {
            ncm_message (" % 20.15g", ncm_matrix_get (res, k, j));
          }
          ncm_message ("\n");
        }
        ncm_message ("\n\n");

        ncm_mset_func_free (mset_func);
        ncm_matrix_free (res);
        g_free (func_name);
        g_free (x);
        g_array_unref (lims);
      }
      
      ncm_vector_clear (&z_vec);
      nc_distance_clear (&dist);
    }

    /*********************************************************************************************************
     * 
     * Visual HW test
     * 
     *********************************************************************************************************/
    if (visual_hw != NULL)
    {
      guint nparams = g_strv_length (visual_hw);

      for (i = 0; i < nparams; i++)
      {
        const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, visual_hw[i]);
        gchar *end_ptr          = NULL;
        glong add_param         = strtol (visual_hw[i], &end_ptr, 10);
        guint k;

        if (pi == NULL && (visual_hw[i] == end_ptr))
        {
          g_warning ("# Parameter `%s' not found for the HW test, skipping...\n", visual_hw[i]);
          continue;
        }

        {
          gint fpi          = (pi != NULL) ? (ncm_mset_fparam_get_fpi (mset, pi->mid, pi->pid) + mcat->nadd_vals) : add_param;
          gdouble mean      = 0.0;
          gdouble var       = 0.0;
          NcmVector *cumsum = ncm_stats_vec_visual_heidel_diag (mcat->pstats, fpi, 0, &mean, &var);
          const gdouble sd  = sqrt (var);

          for (k = 0; k < ncm_vector_len (cumsum); k++)
          {
            const gdouble mean_n   = mean * (k + 1.0);
            const gdouble cumsum_k = ncm_vector_get (cumsum, k);

            ncm_message ("%f % 20.15g % 20.15g % 20.15g\n", 
                         (gdouble)k + 1.0, 
                         mean_n,
                         cumsum_k,
                         sd);
          }
          ncm_vector_free (cumsum);
        }
        ncm_message ("\n\n");
      }
    }

    if (dump)
    {
      const guint cat_size = ncm_mset_catalog_len (mcat);
      guint i;

      for (i = 0; i < cat_size; i++)
      {
        NcmVector *vec = ncm_mset_catalog_peek_row (mcat, i);
        ncm_vector_log_vals (vec, "", "% 22.15g", TRUE);
      }
    }

    if (dump_chain >= 0)
    {
      if (dump_chain >= mcat->nchains)
      {
        g_error ("# Chain number %d does not exists, the catalog contains %u chains!", dump_chain, mcat->nchains);
      }
      else
      {
        NcmStatsVec *chain = g_ptr_array_index (mcat->chain_pstats, dump_chain);
        guint nitens = ncm_stats_vec_nitens (chain);
        guint i;

        for (i = 0; i < nitens; i++)
        {
          NcmVector *vec = ncm_stats_vec_peek_row (chain, i);
          ncm_vector_log_vals (vec, "", "% 22.15g", TRUE);
        }
      }
    }
    
    ncm_mset_clear (&mset);
    ncm_mset_catalog_clear (&mcat);
  } 

  return 0;
}