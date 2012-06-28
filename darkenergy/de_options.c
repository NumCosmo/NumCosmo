/***************************************************************************
 *            de_options.c
 *
 *  Sat Apr 24 14:40:36 2010
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

#include <stdio.h>
#include <glib.h>
#include "de_options.h"

/************************************************************************************************************
 *  Main group entries
 ************************************************************************************************************/

GOptionGroup *
nc_de_opt_get_model_group (NcDEModelEntries *de_model)
{
  GOptionEntry model_entries[] =
  {
	{ "model",    'm', 0, G_OPTION_ARG_STRING, &de_model->model_name, "Name of the darkenergy model to be used.",                      NULL},
	{ "H0",       'h', 0, G_OPTION_ARG_DOUBLE, &de_model->H0,         "Hubble parameter (Km/s/Mpc)",                                 "H0" },
	{ "Omega_c",  'c', 0, G_OPTION_ARG_DOUBLE, &de_model->Omega_c,    "Dark matter density",                                         NULL },
	{ "Omega_x",  'x', 0, G_OPTION_ARG_DOUBLE, &de_model->Omega_x,    "Dark energy density (when changing variable Omega_k)",        NULL },
	{ "Omega_b",  'b', 0, G_OPTION_ARG_DOUBLE, &de_model->Omega_b,    "Baryonic density",                                            NULL },
	{ "T_gamma",  't', 0, G_OPTION_ARG_DOUBLE, &de_model->T_gamma,    "CMB radiation temperature in the present epoch",              NULL },
	{ "n_s",      'n', 0, G_OPTION_ARG_DOUBLE, &de_model->n_s,        "Spectral index",                                              NULL },
	{ "sigma_8",  's', 0, G_OPTION_ARG_DOUBLE, &de_model->sigma_8,    "Sigma_8",                                                     NULL },
	{ "omega_0",  'w', 0, G_OPTION_ARG_DOUBLE, &de_model->w[0],       "Equation of state parameter w0",                              "w0" },
	{ "omega_1",    0, 0, G_OPTION_ARG_DOUBLE, &de_model->w[1],       "Equation of state parameter w1",                              "w1" },
	{ "omega_2",    0, 0, G_OPTION_ARG_DOUBLE, &de_model->w[2],       "Equation of state parameter w2",                              "w2" },
	{ "Omega_k",    0, 0, G_OPTION_ARG_NONE,   &de_model->Omega_k,    "Change variable Omega_x -> Omega_k.",                         NULL },
	{ "flat",       0, 0, G_OPTION_ARG_NONE,   &de_model->flat,       "Change variable Omega_x -> Omega_k and set Omega_k to zero.", NULL },
	{ "help-names", 0, 0, G_OPTION_ARG_NONE,   &de_model->help_names, "Print the parameters names of the chosen model", NULL },
	{ NULL }
  };
  GOptionGroup *model_group = g_option_group_new ("model", " - Dark energy model options", "Show help options related to dark energy model", NULL, NULL);
  g_option_group_add_entries (model_group, model_entries);

  return model_group;
}

GOptionGroup *
nc_de_opt_get_data_simple_group (NcDEDataSimpleEntries *de_data_simple)
{
  GOptionEntry data_simple_entries[] =
  {
	{ "snia-id",    'S', 0, G_OPTION_ARG_INT,    &de_data_simple->snia_id,    "ID of the SNe Ia sample to use",                      NULL },
	{ "bao-id",     'B', 0, G_OPTION_ARG_INT,    &de_data_simple->bao_id,     "ID of the BAO sample to use",                         NULL },
	{ "cmb-id",     'C', 0, G_OPTION_ARG_INT,    &de_data_simple->cmb_id,     "ID of the CMB sample to use",                         NULL },
	{ "H-id",       'E', 0, G_OPTION_ARG_INT,    &de_data_simple->H_id,       "Use the H(z_i) data sample",                          NULL },
	{ "cluster-id", 'U', 0, G_OPTION_ARG_INT,    &de_data_simple->cluster_id, "Use cluster abundance data",                          NULL },
	{ "H0_Hst",     'H', 0, G_OPTION_ARG_NONE,   &de_data_simple->H0_Hst,     "Use the HST H0 data (single gaussian (H0 - 73.8) / 2.4)", NULL },
	{ "BBN",        'N', 0, G_OPTION_ARG_NONE,   &de_data_simple->BBN,        "Use BBN Prior",                                       NULL },
	{ NULL }
  };
  GOptionGroup *data_simple_group = g_option_group_new ("data", " - Cosmological data options", "Show help options related to the dataset to be used.", NULL, NULL);
  g_option_group_add_entries (data_simple_group, data_simple_entries);

  return data_simple_group;
}

GOptionGroup *
nc_de_opt_get_data_cluster_group (NcDEDataClusterEntries *de_data_cluster)
{
  GOptionEntry entries_cluster[] =
  {
	{ "window",        0, 0, G_OPTION_ARG_STRING,         &de_data_cluster->window_name,           "Which window to use, NcWindowTophat: Top Hat (default), NcWindowGaussian: Gaussian", NULL },
	{ "transfer",      0, 0, G_OPTION_ARG_STRING,         &de_data_cluster->transfer_name,         "Which transfer function to use, NcTransferFuncBBKS: BBKS, NcTransferFuncEH: EH (default)", NULL },
	{ "multiplicity",  0, 0, G_OPTION_ARG_STRING,         &de_data_cluster->multiplicity_name,     "Which multiplicity function to use, NcMultiplicityFuncPS: Press-Schechter, NcMultiplicityFuncST: Sheth-Tormen, NcMultiplicityFuncJenkins: Jenkins et al., NcMultiplicityFuncWarren: Warren et al., NcMultiplicityFuncTinker: Tinker (default), NcMultiplicityFuncTinkerMean: Tinker et al. (mean), NcMultiplicityFuncTinkerCrit: Tinker et al. (critical)", NULL },
	{ "cluster-m",     0, 0, G_OPTION_ARG_STRING,         &de_data_cluster->clusterm_ser,          "Which NcClusterMass object to use", NULL },
	{ "cluster-z",     0, 0, G_OPTION_ARG_STRING,         &de_data_cluster->clusterz_ser,          "Which NcClusterRedshift object to use", NULL },
	{ "mf_ds_index",   0, 0, G_OPTION_ARG_INT,            &de_data_cluster->mf_ds_index,           "Determines the coefficients of the multiplicity function.", NULL },
	{ "use-true-data", 0, 0, G_OPTION_ARG_NONE,           &de_data_cluster->use_true_data,         "Use true mass and redshift, must be avaliable.", NULL },
	{ "binned",        0, 0, G_OPTION_ARG_NONE,           &de_data_cluster->binned,                "Binned analyses.", NULL },
	{ "photoz",        0, 0, G_OPTION_ARG_NONE,           &de_data_cluster->use_photoz,            "Include photometric redshift uncertainty.", NULL },
	{ "Mobs",          0, 0, G_OPTION_ARG_NONE,           &de_data_cluster->use_Mobs,              "Include mass observable relation.", NULL },
	{ "binmass",       0, 0, G_OPTION_ARG_NONE,           &de_data_cluster->binmass,               "Use dNdz integrating the full mass function in each z.", NULL },
	{ "Mi",            0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->Mi,                    "Ncuster mass minimum", "M_min" },
	{ "Mf",            0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->Mf,                    "Ncuster mass maximum", "M_max" },
	{ "area",          0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->area_survey,           "User must provide the area in square degree. The conversion to steradian is done internally.", NULL },
	{ "z_initial",     0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->z_initial,             "Initial redshift", NULL },
	{ "z_final",       0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->z_final,               "Final redshift", NULL },
	{ "photoz_sigma0", 0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->photoz_sigma0,         "Standard deviation of photometric redshift sigma = sigma0 (1 + z).", NULL },
	{ "photoz_bias",   0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->photoz_bias,           "Bias of the photometric redshift relation.", NULL },
	{ "lnM_sigma0",    0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->lnM_sigma0,            "Standard deviation of mass observable relation sigma0.", NULL },
	{ "lnM_bias",      0, 0, G_OPTION_ARG_DOUBLE,         &de_data_cluster->lnM_bias,              "Bias of the mass-observable relation.", NULL },
	{ "n_bins",        0, 0, G_OPTION_ARG_INT,            &de_data_cluster->n_bins,                "Number of bins", NULL },
	{ "catalog",       0, 0, G_OPTION_ARG_FILENAME_ARRAY, &de_data_cluster->cata_file,             "Use the folowing catalog as the observational data. It can be used multiple times.", "catalog.dat"},
	{ "save-cat",      0, 0, G_OPTION_ARG_FILENAME,       &de_data_cluster->save_cata,             "Use this option to save the catalog used. (will overwrite)", NULL },
	{ "print_mf",      0, 0, G_OPTION_ARG_NONE,           &de_data_cluster->print_mass_function,   "Create a file and print the mass function from the used FITS catalog and the theoretical one.", NULL },
	{ NULL }
  };
  GOptionGroup *data_cluster_group = g_option_group_new ("cluster", " - Include cluster number counts\n\t - Use --cluster-id 0 and --catalog to use a fit file catalog\n\t - Use --cluster-id 1 and --catalog to use a text file catalog (plain two columns redshift and ln mass)\n\t - Use --cluster-id 2 to make a mock catalog from theory\n", "Show help options related to cluster", NULL, NULL);
  g_option_group_add_entries (data_cluster_group, entries_cluster);
  return data_cluster_group;
}

static gboolean
_nc_de_print_fit_list (const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
  printf ("<FIXME> PRINT HERE ALL ALGORITHMS\n");
  return FALSE;
}

GOptionGroup *
nc_de_opt_get_fit_group (NcDEFitEntries *de_fit)
{
  GOptionEntry fit_entries[] =
  {
	{ "out",            0, 0, G_OPTION_ARG_FILENAME,     &de_fit->file_out,      "Output filename.", "output.dat" },
	{ "minalgo",        0, 0, G_OPTION_ARG_INT,          &de_fit->min_algo,      "Minimization algorithim to be used.", NULL },
	{ "diffalgo",       0, 0, G_OPTION_ARG_INT,          &de_fit->diff_algo,     "Differentiation algorithim to be used.", NULL },
	{ "fit-list",       0, G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, &_nc_de_print_fit_list,  "Print all the minimization/differentiation algorithims avaliable.", NULL },
	{ "n-sigma",        0, 0, G_OPTION_ARG_INT,          &de_fit->nsigma,        "Confidence region probability 1, 2 or 3 sigmas. A zero value calculate all three confidence regions.", NULL },
	{ "n-sigma-fisher", 0, 0, G_OPTION_ARG_INT,          &de_fit->nsigma_fisher, "Confidence region (Fisher matrix) probability 1, 2 or 3 sigmas. A zero value calculate all three confidence regions.", NULL },
	{ "cr-x",           0, 0, G_OPTION_ARG_INT,          &de_fit->bidim_cr[0],   "Confidence region x parameter", NULL },
	{ "cr-y",           0, 0, G_OPTION_ARG_INT,          &de_fit->bidim_cr[1],   "Confidence region y parameter", NULL },
	{ "max-iter",       0, 0, G_OPTION_ARG_INT,          &de_fit->max_iter,      "Max number of iterations used by the minimization algorithms", NULL },
	{ "err-param",      0, 0, G_OPTION_ARG_STRING_ARRAY, &de_fit->onedim_cr,     "Calculate the one dimensional confidence region", NULL },
	{ "fit-params",     0, 0, G_OPTION_ARG_STRING,       &de_fit->fit_params,    "Parameters to be fitted, use fit-params H0=fit,Omega_m=fix,...", NULL },
	{ "resample",       0, 0, G_OPTION_ARG_NONE,         &de_fit->resample,      "Resample model using default params", NULL },
	{ "msg-level",      0, 0, G_OPTION_ARG_INT,          &de_fit->msg_level,     "Fit message level (0: no msg, 1: simple, 2: full)", NULL },
	{ "montecarlo",     0, 0, G_OPTION_ARG_INT,          &de_fit->montecarlo,    "Resample the original data 'montecarlo' times.", NULL},
	{ "mc-ni",          0, 0, G_OPTION_ARG_INT,          &de_fit->mc_ni,         "Start the 'montecarlo' at the ni realization.", NULL},
	{ "fiducial",       0, 0, G_OPTION_ARG_NONE,         &de_fit->fiducial,      "Use the fiducial model to resample.", NULL},
	{ "mc-data",        0, 0, G_OPTION_ARG_NONE,         &de_fit->mc_data,       "Print all data from monte carlo.", NULL},
	{ "fit",            0, 0, G_OPTION_ARG_NONE,         &de_fit->fit,           "Fit model using the selected data.", NULL},
	{ "fisher",         0, 0, G_OPTION_ARG_NONE,         &de_fit->fisher,        "Calculated the Fisher matrix.", NULL},
	{ "save-fisher",    0, 0, G_OPTION_ARG_NONE,         &de_fit->save_fisher,   "Create a file and print the Fisher matrix.", NULL},
	{ "save-best-fit",  0, 0, G_OPTION_ARG_NONE,         &de_fit->save_best_fit, "Create a file and print the cosmological parameters (both best-fit and fixed ones).", NULL},
	{ NULL }
  };
  GOptionGroup *fit_group = g_option_group_new ("fit", " - Choices of statistical analysis", "Show help options related to model fitting", NULL, NULL);
  g_option_group_add_entries (fit_group, fit_entries);
  return fit_group;
}
