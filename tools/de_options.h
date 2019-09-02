/***************************************************************************
 *            de_options.h
 *
 *  Sat Apr 24 14:40:43 2010
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

typedef struct _NcDERunEntries NcDERunEntries;

/**
 * NcDERunEntries:
 *
 * FIXME
 */
struct _NcDERunEntries
{
  gchar *runconf;
  gchar *saverun;
  gint64 main_seed;
  guint nthreads;
};

#define NC_DE_RUN_ENTRIES {NULL, NULL, -1, 0}

typedef struct _NcDEModelEntries NcDEModelEntries;

/**
 * NcDEModelEntries:
 *
 * FIXME
 */
struct _NcDEModelEntries
{
  gchar *mset_file;
  gchar *model_name;
  gchar *model_reion;
  gchar *model_prim;
  gboolean flat;
  gboolean pos_Omega_x;
  gboolean Omega_k;
  gboolean help_names;
};

#define NC_DE_MODEL_ENTRIES {NULL, NULL, NULL, NULL, FALSE, FALSE, FALSE, FALSE}

typedef struct _NcDEDataSimpleEntries NcDEDataSimpleEntries;

/**
 * NcDEDataSimpleEntries:
 *
 * FIXME
 */
struct _NcDEDataSimpleEntries
{
  gchar *snia_id;
  gchar *snia_objser;
  gchar **bao_id;
  gchar *cmb_id;
  gchar **H_id;
  gchar **H_BAO_id;
  gchar *cluster_id;
  gchar **priors_gauss;
  gchar **Planck;
  gchar *PlanckFI;
  gboolean PlanckPriors;
  gboolean BBN;
  gboolean BBN_Ob;
  gboolean snia_use_det;
  gchar **data_files;
};

#define NC_DE_DATA_SIMPLE_ENTRIES {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, FALSE, FALSE, FALSE, FALSE, NULL}

typedef struct _NcDEDataClusterEntries NcDEDataClusterEntries;

/**
 * NcDEDataNcusterEntries:
 *
 * FIXME
 * 
 */
struct _NcDEDataClusterEntries
{
  gchar *filter_type;
  gchar *ps_type;
  gchar *multiplicity_name;
  gchar *clusterm_ser;
  gchar *clusterz_ser;
  gint mf_ds_index;
  gboolean use_true_data;
  gboolean binned;
  gboolean binmass;
  gboolean use_Mobs_local;
  gboolean use_selection;
  gdouble area_survey;
  gint n_bins;
  gchar **cata_file;
  gchar *save_cata;
  gboolean print_mass_function;
};

#define NC_DE_DATA_CLUSTER_ENTRIES {NULL, NULL, NULL, NULL, NULL, 0, FALSE, FALSE, FALSE, FALSE, FALSE, 5000.0, 10, NULL, NULL, FALSE}

typedef struct _NcDEFitEntries NcDEFitEntries;

/**
 * NcDEFitEntries:
 *
 * FIXME
 * 
 */
struct _NcDEFitEntries
{
  gchar *file_out;
  gchar *fit_type;
  gchar *fit_diff;
  gchar *fit_algo;
  gdouble fit_reltol;
  gdouble fit_params_reltol;
  gint nsigma;
  gint nsigma_fisher;
  gchar *bidim_cr[2];
  gchar **onedim_cr;
  gchar **funcs;
  gdouble lhr_prec;
  gint max_iter;
  gboolean resample;
  gint msg_level;
  gint mc_rtype;
  gint mc_ni;
  gint mc_nthreads;
  glong mc_seed;
  gint mc_nwalkers;
  gint mc_prerun;
  gint mcbs_nbootstraps;
  gdouble mc_lre;
  gchar *fiducial;
  gchar *mc_data;
  gboolean mc_unordered;
  gboolean fit;
  gboolean restart;
  gdouble  restart_abstol;
  gdouble  restart_reltol;
  gboolean mc;
  gboolean mcbs;
  gboolean mcmc;
  gboolean esmcmc;
  gboolean esmcmc_walk;
  gboolean esmcmc_aps;
  gboolean esmcmc_sbox;
  gboolean esmcmc_ms;
  gint fisher;
  gboolean qspline_cp;
  gdouble qspline_cp_sigma;
  gboolean save_fisher;
  gboolean save_best_fit;
  gchar *save_mset;
};

#define NC_DE_FIT_ENTRIES { NULL, NULL, NULL, NULL, 1e-8, 1e-5, -1, -1, {NULL, NULL}, NULL, NULL, 1.0e-5, NCM_FIT_DEFAULT_MAXITER, FALSE, NCM_FIT_RUN_MSGS_SIMPLE, NCM_FIT_MC_RESAMPLE_FROM_MODEL, 0, 0, -1, 100, 0, 100, 1.0e3, NULL, NULL, FALSE, FALSE, FALSE, 0.0, 1.0e-4, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 0, FALSE, 1.0, FALSE, FALSE, NULL}

GOptionGroup *nc_de_opt_get_run_group (NcDERunEntries *de_run);
GOptionGroup *nc_de_opt_get_model_group (NcDEModelEntries *de_model, GOptionEntry **de_model_entries);
GOptionGroup *nc_de_opt_get_data_simple_group (NcDEDataSimpleEntries *de_data_simple, GOptionEntry **de_data_simple_entries);
GOptionGroup *nc_de_opt_get_data_cluster_group (NcDEDataClusterEntries *de_data_cluster, GOptionEntry **de_data_cluster_entries);
GOptionGroup *nc_de_opt_get_fit_group (NcDEFitEntries *de_fit_group, GOptionEntry **de_fit_entries);
