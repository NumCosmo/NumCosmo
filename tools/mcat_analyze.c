/***************************************************************************
 *            mcat_analyze.c
 *
 *  Mon March 16 11:04:23 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@iaoftware.com.br>
 ****************************************************************************/
/*
 * mcat_analyze.c
 *
 * Copyright (C) 2015 - Sandro Dias Pinto Vitenti
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

gint
main (gint argc, gchar *argv[])
{
  gchar *cat_filename = NULL;
  gboolean info           = FALSE;
  gboolean list_all       = FALSE;
  gboolean list_hicosmo   = FALSE;
  gboolean list_hicosmo_z = FALSE;
  gboolean list_dist      = FALSE;
  gboolean list_dist_z    = FALSE;
  gboolean use_direct     = FALSE;
  gdouble zi = 0.0;
  gdouble zf = 1.0;
  gint nsteps = 100;
  gint burnin = 0;
  gchar **funcs = NULL;
  gchar **distribs = NULL;
  gchar **params = NULL;
  guint i;
  
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry entries[] =
  {
    { "catalog",        'c', 0, G_OPTION_ARG_FILENAME,     &cat_filename,   "Catalog filename.", NULL },
    { "info",           'i', 0, G_OPTION_ARG_NONE,         &info,           "Print catalog informations.", NULL },
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
    { "distribution",   'd', 0, G_OPTION_ARG_STRING_ARRAY, &distribs,       "Functions distributions to be analyzed.", NULL},
    { "parameter",      'p', 0, G_OPTION_ARG_STRING_ARRAY, &params,         "Model parameters' to be analyzed.", NULL},
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
    GArray *func_table = nc_hicosmo_class_func_table ();
    ncm_cfg_msg_sepa ();
    g_message ("# Available constant functions from NcHICosmo models:\n");    
    for (i = 0; i < func_table->len; i++)
    {
      NcHICosmoFunc *f = &g_array_index (func_table, NcHICosmoFunc, i);
      g_message ("# - %32s: %s\n", f->name, f->desc);
    }
    g_array_unref (func_table);
  }
  if (list_hicosmo_z || list_all)
  {
    GArray *func_z_table = nc_hicosmo_class_func_z_table ();
    ncm_cfg_msg_sepa ();
    g_message ("# Available redshift functions from NcHICosmo models:\n");    
    for (i = 0; i < func_z_table->len; i++)
    {
      NcHICosmoFuncZ *fz = &g_array_index (func_z_table, NcHICosmoFuncZ, i);
      g_message ("# - %32s: %s\n", fz->name, fz->desc);
    }
    g_array_unref (func_z_table);
  }
  if (list_dist || list_all)
  {
    GArray *func_table = nc_distance_class_func_table ();
    ncm_cfg_msg_sepa ();
    g_message ("# Available constant functions from NcDistance:\n");    
    for (i = 0; i < func_table->len; i++)
    {
      NcDistanceFunc *f = &g_array_index (func_table, NcDistanceFunc, i);
      g_message ("# - %32s: %s\n", f->name, f->desc);
    }
    g_array_unref (func_table);
  }
  if (list_dist_z || list_all)
  {
    GArray *func_z_table = nc_distance_class_func_z_table ();
    ncm_cfg_msg_sepa ();
    g_message ("# Available redshift functions from NcDistance:\n");    
    for (i = 0; i < func_z_table->len; i++)
    {
      NcDistanceFuncZ *fz = &g_array_index (func_z_table, NcDistanceFuncZ, i);
      g_message ("# - %32s: %s\n", fz->name, fz->desc);
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
    NcmMSetCatalog *mcat = ncm_mset_catalog_new_from_file (cat_filename);
    NcmMSet *mset = ncm_mset_catalog_get_mset (mcat);
    
    if (info)
    {
      NcmMatrix *cov = NULL;
      const gchar *rtype_str = ncm_mset_catalog_get_run_type (mcat);
      ncm_cfg_msg_sepa ();
      g_message ("# Catalog run type: `%s'.\n", rtype_str);
      ncm_mset_pretty_log (mset);
      ncm_mset_catalog_get_covar (mcat, &cov);
      ncm_mset_fparams_log_covar (mset, cov);
      ncm_matrix_clear (&cov);
      ncm_mset_catalog_log_current_stats (mcat);
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
          NcHICosmoFuncZ *fz = nc_hicosmo_class_get_func_z (funcs[i]);
          if (fz != NULL)
          {
            if (!ncm_model_check_impl (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ())), fz->impl))
            {
              g_warning ("# function `%s' not supported by NcHICosmo `%s', skipping...", funcs[i],
                         ncm_model_name (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ()))));
              continue;
            }
            mset_func = nc_hicosmo_create_mset_arrayfunc1 (fz->f, nsteps);
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcHICosmo z function: `%s' in [% 20.15g % 20.15g].\n", fz->desc, zi, zf);
          }
        }
        if (mset_func == NULL)
        {
          NcDistanceFuncZ *fz = nc_distance_class_get_func_z (funcs[i]);
          if (fz != NULL)
          {
            if (!ncm_model_check_impl (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ())), fz->impl))
            {
              g_warning ("# function `%s' not supported by NcHICosmo `%s', skipping...", funcs[i],
                         ncm_model_name (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ()))));
              continue;
            }
            mset_func = nc_distance_create_mset_arrayfunc1 (dist, fz->f, nsteps);
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcDistance z function: `%s' in [% 20.15g % 20.15g].\n", fz->desc, zi, zf);
          }
        }
        if (mset_func == NULL)
        {
          g_warning ("# Function `%s' not found, skipping...\n", funcs[i]);
          continue;
        }

        if (use_direct)
          res = ncm_mset_catalog_calc_ci_direct (mcat, burnin, mset_func, ncm_vector_ptr (z_vec, 0), p_val);
        else
          res = ncm_mset_catalog_calc_ci_interp (mcat, burnin, mset_func, ncm_vector_ptr (z_vec, 0), p_val, 100, NCM_FIT_RUN_MSGS_SIMPLE);

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
        guint k;

        if (mset_func == NULL)
        {
          NcHICosmoFunc *f = nc_hicosmo_class_get_func (distribs[i]);
          if (f != NULL)
          {
            if (!ncm_model_check_impl (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ())), f->impl))
            {
              g_warning ("# function `%s' not supported by NcHICosmo `%s', skipping...", distribs[i],
                         ncm_model_name (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ()))));
              continue;
            }
            mset_func = nc_hicosmo_create_mset_func0 (f->f);
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcHICosmo function: `%s'.\n", f->desc);
          }
        }
        if (mset_func == NULL)
        {
          NcDistanceFunc *f = nc_distance_class_get_func (distribs[i]);
          if (f != NULL)
          {
            if (!ncm_model_check_impl (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ())), f->impl))
            {
              g_warning ("# function `%s' not supported by NcHICosmo `%s', skipping...", distribs[i],
                         ncm_model_name (NCM_MODEL (ncm_mset_peek (mset, nc_hicosmo_id ()))));
              continue;
            }
            mset_func = nc_distance_func0_new (dist, f->f);
            ncm_cfg_msg_sepa ();
            g_message ("# Printing NcDistance function: `%s'.\n", f->desc);
          }
        }
        if (mset_func == NULL)
        {
          g_warning ("# Function `%s' not found, skipping...\n", distribs[i]);
          continue;
        }

        {
          NcmStatsDist1d *sd1 = ncm_mset_catalog_calc_distrib (mcat, burnin, mset_func, NCM_FIT_RUN_MSGS_SIMPLE);
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
        NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (mset, params[i]);
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
            sd1 = ncm_mset_catalog_calc_param_distrib (mcat, burnin, pi, NCM_FIT_RUN_MSGS_SIMPLE);
          else
            sd1 = ncm_mset_catalog_calc_add_param_distrib (mcat, burnin, add_param, NCM_FIT_RUN_MSGS_SIMPLE);
          
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
    
    ncm_mset_clear (&mset);
    ncm_mset_catalog_clear (&mcat);
  }
  
  return 0;
}
