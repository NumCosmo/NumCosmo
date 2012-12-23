/***************************************************************************
 *            q-piecewise.c
 *
 *  Mon Jun 18 19:46:26 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
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
#include <gmp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>

#define RESAMPLE_OMEGA_M 0.30f
#define RESAMPLE_OMEGA_LAMBDA ((1.0f - RESAMPLE_OMEGA_M) + 0.00f)
#define RESAMPLE_OMEGA -1.0f

static gdouble last_z = 0.80;
static gdouble omega_k = 0.001;
static gdouble z0 = 0.10;
static gdouble interval = 0.30;
static gdouble j_sigma = 0.0;
static gint nsigma = 0;
static gint max_iter = 10000;
static gint snia_id = -1;
static gint H_id = -1;
static gint msg_level = NCM_FIT_RUN_MSGS_SIMPLE;
static gboolean curved = FALSE;
static gboolean with_E = FALSE;
static gboolean with_E_CABRE = FALSE;
static gboolean with_BAO = FALSE;
static gboolean resample = FALSE;
static gboolean print_data = FALSE;
static gboolean print_E = FALSE;
static gboolean least_squares = FALSE;
static gboolean change_params = FALSE;
static gchar *file_out = NULL;
static gint min_algo = NCM_FIT_TYPE_NLOPT;
static gint diff_algo = NCM_FIT_GRAD_NUMDIFF_FORWARD;

static GOptionEntry entries[] =
{
  { "last-z",        'l', 0, G_OPTION_ARG_DOUBLE,   &last_z,        "Last redshift of the partitioned interval (default 0.80)", NULL },
  { "z0",            'z', 0, G_OPTION_ARG_DOUBLE,   &z0,            "Redshift of the variable change", NULL },
  { "omega_k",       'O', 0, G_OPTION_ARG_DOUBLE,   &omega_k,       "Omega_k", NULL },
  { "interval",      'i', 0, G_OPTION_ARG_DOUBLE,   &interval,      "Interval size (default 0.30)", NULL },
  { "j-sigma",       'j', 0, G_OPTION_ARG_DOUBLE,   &j_sigma,       "Standard deviation on the j continuity prior (default 0.0 off)", NULL },
  { "n-sigma",       'n', 0, G_OPTION_ARG_INT,      &nsigma,        "Confidence region probability", NULL },
  { "max-iter",      'k', 0, G_OPTION_ARG_INT,      &max_iter,      "Max number of iterations used by the minimization algorithms", NULL },
  { "snia-id",       's', 0, G_OPTION_ARG_INT,      &snia_id,       "ID of the sample to use", NULL },
  { "H-type",        'E', 0, G_OPTION_ARG_INT,      &H_id,          "Use the H(z_i) data sample", NULL },
  { "curved",        'C', 0, G_OPTION_ARG_NONE,     &curved,        "Fit also Omega_k", NULL },
  { "model-id",      'm', 0, G_OPTION_ARG_NONE,     &curved,        "Fit also Omega_k", NULL },
  { "resample",      'r', 0, G_OPTION_ARG_NONE,     &resample,      "Resample using LCDM (0.30, 0.70)", NULL },
  { "with-E",        'E', 0, G_OPTION_ARG_NONE,     &with_E,        "Use the H(z_i) data sample", NULL },
  { "with-E-CABRE",    0, 0, G_OPTION_ARG_NONE,     &with_E_CABRE,  "Use the H(z_i) data sample from 0807.3551 [astro-ph]", NULL },
  { "with-BAO",      'B', 0, G_OPTION_ARG_NONE,     &with_BAO,      "Use the BAO prior", NULL },
  { "out",            0,  0, G_OPTION_ARG_FILENAME, &file_out,      "Output filename.", "output.dat" },
  { "minalgo",        0,  0, G_OPTION_ARG_INT,      &min_algo,      "Minimization algorithim to be used.", NULL },
  { "diffalgo",       0,  0, G_OPTION_ARG_INT,      &diff_algo,     "Differentiation algorithim to be used.", NULL },
  { "change-params", 'c', 0, G_OPTION_ARG_NONE,     &change_params, "Change the parametrization to q_0 -> lnE(z0) and q_i -> q(z0)", NULL },
  { "print-data",    'd', 0, G_OPTION_ARG_NONE,     &print_data,    "Print the fitted data", NULL },
  { "print-E",       'P', 0, G_OPTION_ARG_NONE,     &print_E,       "Print values of E in the interval", NULL },
  { "msg-level",     'v', 0, G_OPTION_ARG_INT,      &msg_level,     "Be verbose", NULL },
  { NULL }
};

gint
main(gint argc, gchar *argv[])
{
  GError *error = NULL;
  GOptionContext *context;
  NcHICosmoQPW *qpw;
  NcHICosmoLCDM *lcdm;
  NcmDataset *dset;
  NcmLikelihood *lh;
  NcmFit *fit = NULL;
  NcDistance *dist = nc_distance_new (2.0);
  NcmMSet *mset, *mset_lcdm;
  gint i, j;

  ncm_cfg_init ();

  context = g_option_context_new ("- test the q piecewise/spline model");
  g_option_context_add_main_entries (context, entries, NULL);
  g_option_context_parse (context, &argc, &argv, &error);

  dset = ncm_dataset_new ();
  qpw = nc_hicosmo_qpw_new (interval, last_z, !curved);
  lcdm = nc_hicosmo_lcdm_new ();

  mset = ncm_mset_new (NCM_MODEL (qpw), NULL);
  mset_lcdm = ncm_mset_new (NCM_MODEL (lcdm), NULL);

  if (snia_id != -1)
  {
    NcmData *snia = nc_data_dist_mu_new (dist, snia_id);
    ncm_dataset_append_data (dset, snia);
    ncm_data_free (snia);
  }

  lh = ncm_likelihood_new (dset);

  if (H_id != -1)
  {
    NcmData *H_data = nc_data_hubble_new (H_id);
    ncm_dataset_append_data (dset, H_data);
    ncm_mset_param_set_ftype (mset, NC_HICOSMO_ID, NC_HICOSMO_QPW_H0, NCM_PARAM_TYPE_FREE);
    ncm_data_free (H_data);
  }

  if (resample)
    ncm_dataset_resample (dset, mset_lcdm);

  if (change_params)
  {
    gint wpiece = nc_hicosmo_qpw_index (qpw, z0) + 3;
    nc_hicosmo_qpw_change_params (qpw, z0);
    ncm_model_param_set (NCM_MODEL (qpw), 2, 0.01);
    ncm_model_param_set (NCM_MODEL (qpw), wpiece, -0.2);
  }

  if (with_BAO)
  {
    NcmData *bao_data = nc_data_bao_create (dist, NC_DATA_BAO_DVDV_PERCIVAL2007);
    ncm_dataset_append_data (dset, bao_data);
    ncm_data_free (bao_data);
  }

  if (j_sigma != 0.0f)
    nc_hicosmo_qpw_add_continuity_priors (qpw, lh, j_sigma);

#define SIZE 10000
  if (least_squares)
  {
    gint wpiece = nc_hicosmo_qpw_index (qpw, z0) + 3;
    //    nc_hicosmo_qpw_add_asymptotic_cdm_prior (lh, 500.0f, 1.0f/2.0f, 0.0001f);
    //    nc_hicosmo_qpw_add_asymptotic_cdm_prior (lh, 600.0f, 1.0f/2.0f, 0.0001f);
    fit = ncm_fit_new (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_ANALYTICAL);
    ncm_fit_run (fit, msg_level);
    ncm_fit_log_info (fit);
    ncm_fit_numdiff_m2lnL_covar (fit);
    ncm_fit_log_covar (fit);
    printf ("# jerk: % 12.4g\n", nc_hicosmo_j (NC_HICOSMO (qpw), 0.0));

    //ncm_fit_error (fit, 2, ncm_c_stats_1sigma (), 1, ERR(0));printf("\n");

    while (FALSE)
    {
      gdouble old_s = 0.0;
      gdouble s = 0.0;
      guint n;
      for (n = 0; n < ncm_likelihood_priors_length (lh); n++)
      {
        NcmMSetFunc *pdata = NCM_MSET_FUNC (ncm_likelihood_priors_peek (lh, n));
        NcHICosmoQPWContPrior *cprior = (NcHICosmoQPWContPrior *)pdata->obj;
        old_s = cprior->sigma;
        s += gsl_pow_2 (ncm_mset_func_eval0 (pdata, mset));
        printf ("Bunga[%d] = %g\n", n, gsl_pow_2 (ncm_mset_func_eval0 (pdata, mset)));
      }
      s = ((old_s) + s/n) / 2.0;
      printf ("sigma_qp = %g\n", s);

      for (n = 0; n < ncm_likelihood_priors_length (lh); n++)
      {
        NcmMSetFunc *pdata = NCM_MSET_FUNC (ncm_likelihood_priors_peek (lh, n));
        NcHICosmoQPWContPrior *cprior = (NcHICosmoQPWContPrior *)pdata->obj;
        cprior->sigma = s;
      }
      ncm_fit_run (fit, FALSE);
    }

    if (FALSE)
    {
      gdouble chi2_data, chi2_priors;
      ncm_fit_data_m2lnL_val (fit, &chi2_data);
      ncm_fit_priors_m2lnL_val (fit, &chi2_priors);
      printf ("%6f %6f %6f %6f %6f\n", j_sigma, fit->fstate->m2lnL, chi2_data, chi2_priors, chi2_data + chi2_priors);
    }
    //printf ("%g %g\n",1000.0f, nc_hicosmo_q (cp, 1000.0f));
    //printf("#===>%g<===\n",ncm_fit_GoF_wmean (fit));

    if (FALSE)
    {
      gdouble sigma_t = 0.0;
      gdouble chi2d, step;
      for (i = 0; i < ncm_model_len (NCM_MODEL (qpw)); i++)
        sigma_t += ncm_fit_covar_var (fit, NC_HICOSMO_ID, i);
      sigma_t = sqrt(sigma_t);
      ncm_fit_m2lnL_val (fit, &chi2d);
      printf ("#%g %g %g\n", j_sigma, sigma_t, chi2d);
      for (step = 0.0f; step <= 1.01; step += 1e-4)
      {
        gdouble mu_i = nc_distance_comoving (dist, NC_HICOSMO (qpw), step);
        gdouble mu_i1 = nc_distance_comoving (dist, NC_HICOSMO (qpw), step + 1e-4);
        printf ("%g %g\n", step, 1.0/((mu_i1 - mu_i) / 1e-4));
        if ((mu_i1 - mu_i) / 1e-4 < 0)
          break;
      }
    }

    if (FALSE)
    {
      gdouble pprob = 1.0f-ncm_fit_lr_test (fit, NC_HICOSMO_ID, wpiece, 0.0, 2);
      printf ("%g %g %g\n", z0, pprob, sqrt(gsl_cdf_chisq_Qinv (1.0 - pprob, 1)));
    }

    if (FALSE)
    {
      printf ("%g %g %g %g\n", z0,
              ncm_fit_lr_test (fit, NC_HICOSMO_ID, wpiece, -1.0, 1),
              ncm_fit_lr_test (fit, NC_HICOSMO_ID, wpiece,  0.0, 1),
              ncm_fit_lr_test (fit, NC_HICOSMO_ID, wpiece,  2.0, 1));
      //      printf ("%g %g %g %g\n", z0,
      //              ncm_fit_lr_test (fit, 0, -1.0f),
      //              ncm_fit_lr_test (fit, 0,  0.0f),
      //              ncm_fit_lr_test (fit, 0,  2.0f));
    }

    //ncm_fit_error (fit, wpiece, ncm_c_stats_1sigma (), ERR(wpiece));
    //ncm_fit_error (fit, 0, ncm_c_stats_1sigma (), ERR(0));
    //ncm_fit_error (fit, 1, ncm_c_stats_1sigma (), ERR(1));
    //ncm_fit_error (fit, 2, ncm_c_stats_1sigma (), ERR(2));
    //ncm_fit_jackknife (fit, stdout, verbose);
    //nc_likelihood_jackknife_print (fit, stdout);

    if (FALSE)
    {
      GTimer *prob_time = g_timer_new();
      gdouble range1 = ncm_fit_prob (fit, NC_HICOSMO_ID, wpiece, -5.0, -1.0);
      gdouble range2 = ncm_fit_prob (fit, NC_HICOSMO_ID, wpiece, -1.0, 0.0);
      gdouble range3 = ncm_fit_prob (fit, NC_HICOSMO_ID, wpiece, 0.0, 2.0);
      gdouble range4 = ncm_fit_prob (fit, NC_HICOSMO_ID, wpiece, 2.0, 5.0);
      gdouble norm = range1 + range2 + range3 + range4;
      range1 /= norm;range2 /= norm;range3 /= norm;range4 /= norm;
      //printf ("%g %g %g %g %g\n", z0, range1, range2, range3, range4);
      ncm_fit_dprob (fit, NC_HICOSMO_ID, wpiece, -5.0, 5.0, 0.01, norm);
      printf ("# time %g\n", g_timer_elapsed (prob_time, NULL));
      fflush (stdout);
    }

    if (FALSE)
    {
      gdouble step;
      for (step = 0.10; step <= 0.81; step += 0.05)
      {
        wpiece = nc_hicosmo_qpw_index (qpw, step) + 1;
        nc_hicosmo_qpw_change_params (qpw, step);
        ncm_model_param_set (NCM_MODEL (qpw), 0, 0.0);
        ncm_model_param_set (NCM_MODEL (qpw), wpiece, -0.2);
        fit = ncm_fit_new (NCM_FIT_TYPE_GSL_LS, NULL, lh, mset, NCM_FIT_GRAD_ANALYTICAL);
        ncm_fit_run (fit, msg_level);
        printf ("estimate E(%-8.6f) = % -8.6f, q(%-8.6f) = % -8.6f | err E = % -8.6f, q = % -8.6f\n",
                step, exp(ncm_model_param_get (NCM_MODEL (qpw), 0)),
                step, ncm_model_param_get (NCM_MODEL (qpw), wpiece),
                ncm_fit_covar_sd (fit, NC_HICOSMO_ID, 0), ncm_fit_covar_sd (fit, NC_HICOSMO_ID, wpiece)
                );
        printf ("modelval E(%-8.6f) = % -8.6f, q(%-8.6f) = % -8.6f\n",
                step, sqrt(0.3f*(1.0f+step)*(1.0f+step)*(1.0f+step) + 0.7f),
                step, 1.0f/2.0f * (0.3*(1.0f+step)*(1.0f+step)*(1.0f+step) - 1.4f)/(0.3f*(1.0f+step)*(1.0f+step)*(1.0f+step) + 0.7f)
                );
        printf ("#################################################\n");
      }
    }

    if (FALSE)
    {
      GTimer *timer = g_timer_new ();
      //NcParams *orig_cp = lcdm_cp;
      NcmMSet *orig_mset = ncm_mset_dup (mset);

      for (i = 0 ; i < SIZE ; i++)
      {
        //        gint j;
        if (FALSE)
          ncm_dataset_resample (lh->dset, orig_mset);
 

        ncm_fit_run (fit, FALSE);
        //printf ("# %g\n",fit->m2lnL);
        //nc_params_print_all (cp, stdout);
        printf ("%.16g\n", fit->fstate->m2lnL);
        fflush (stdout);
        if (i % 100 == 99)
        {
          gdouble elapsed = g_timer_elapsed(timer, NULL);
          printf ("# sample size = %d, time = %f, time/sample = %f\n", i+1, elapsed, elapsed/(i+1.0f));
          fflush(stdout);
        }
      }
    }

    if (FALSE)
    {
      gdouble range = 1.01f;
      gsl_matrix *q_matrix = gsl_matrix_alloc (SIZE, 1000);
      gsl_vector *q_bias = gsl_vector_alloc (1000);
      gsl_vector *q_var = gsl_vector_alloc (1000);
      gsl_matrix *jerk_matrix = gsl_matrix_alloc (SIZE, 1000);
      gsl_vector *jerk_bias = gsl_vector_alloc (1000);
      gsl_vector *jerk_var = gsl_vector_alloc (1000);
      GTimer *timer = g_timer_new ();
      for (i = 0 ; i < SIZE ; i++)
      {
        ncm_dataset_resample (lh->dset, mset_lcdm);
        ncm_fit_run (fit, FALSE);
        for (j = 0; j < 1000; j++)
        {
          gdouble step = range / 1000.0f * j;
          gdouble q_pw = nc_hicosmo_q (NC_HICOSMO (qpw), step);
          gdouble q_lcdm = nc_hicosmo_q (NC_HICOSMO (lcdm), step);
          gdouble jerk_pw = nc_hicosmo_j (NC_HICOSMO (qpw), step);
          gdouble jerk_lcdm = nc_hicosmo_j (NC_HICOSMO (lcdm), step);
          gsl_matrix_set (q_matrix, i, j, (q_pw - q_lcdm));
          gsl_matrix_set (jerk_matrix, i, j, (jerk_pw - jerk_lcdm));
        }

        if ( i%100 == 99 )
        {
          gdouble elapsed;
          for (j = 0; j < 1000; j++)
          {
            gdouble step = range / 1000.0f * j;
            gsl_vector_view q_i = gsl_matrix_column (q_matrix, j);
            gdouble q_bias_i = gsl_stats_mean (q_i.vector.data, q_i.vector.stride, i+1);
            gdouble q_var_i = gsl_stats_sd (q_i.vector.data, q_i.vector.stride, i+1);
            gsl_vector_view jerk_i = gsl_matrix_column (jerk_matrix, j);
            gdouble jerk_bias_i = gsl_stats_mean (jerk_i.vector.data, jerk_i.vector.stride, i+1);
            gdouble jerk_var_i = gsl_stats_sd (jerk_i.vector.data, jerk_i.vector.stride, i+1);
            gsl_vector_set(q_bias, j, fabs(q_bias_i));
            gsl_vector_set(q_var, j, q_var_i);
            gsl_vector_set(jerk_bias, j, fabs(jerk_bias_i));
            gsl_vector_set(jerk_var, j, jerk_var_i);
            printf ("%d %g %g %g %g %g\n", i+1, step, q_bias_i, q_var_i, jerk_bias_i, jerk_var_i);
          }
          printf ("\n\n");
          elapsed = g_timer_elapsed(timer, NULL);
          printf ("# sample size = %d, time = %f, time/sample = %f\n", i+1, elapsed, elapsed/(i+1.0f));
          printf ("# q bias (0, r/2): meam = % 8.6f, sigma = % 8.6f; (0, r): meam = % 8.6f, sigma = % 8.6f;\n",
                  gsl_stats_mean (q_bias->data, q_bias->stride, q_bias->size/2),
                  gsl_stats_sd (q_bias->data, q_bias->stride, q_bias->size/2),
                  gsl_stats_mean (q_bias->data, q_bias->stride, q_bias->size),
                  gsl_stats_sd (q_bias->data, q_bias->stride, q_bias->size)
                  );
          printf ("# q var  (0, r/2): meam = % 8.6f, sigma = % 8.6f; (0, r): meam = % 8.6f, sigma = % 8.6f;\n",
                  gsl_stats_mean (q_var->data, q_var->stride, q_var->size/2),
                  gsl_stats_sd (q_var->data, q_var->stride, q_var->size/2),
                  gsl_stats_mean (q_var->data, q_var->stride, q_var->size),
                  gsl_stats_sd (q_var->data, q_var->stride, q_var->size)
                  );
          printf ("# jerk bias (0, r/2): meam = % 8.6f, sigma = % 8.6f; (0, r): meam = % 8.6f, sigma = % 8.6f;\n",
                  gsl_stats_mean (jerk_bias->data, jerk_bias->stride, jerk_bias->size/2),
                  gsl_stats_sd (jerk_bias->data, jerk_bias->stride, jerk_bias->size/2),
                  gsl_stats_mean (jerk_bias->data, jerk_bias->stride, jerk_bias->size),
                  gsl_stats_sd (jerk_bias->data, jerk_bias->stride, jerk_bias->size)
                  );
          printf ("# jerk var  (0, r/2): meam = % 8.6f, sigma = % 8.6f; (0, r): meam = % 8.6f, sigma = % 8.6f;\n",
                  gsl_stats_mean (jerk_var->data, jerk_var->stride, jerk_var->size/2),
                  gsl_stats_sd (jerk_var->data, jerk_var->stride, jerk_var->size/2),
                  gsl_stats_mean (jerk_var->data, jerk_var->stride, jerk_var->size),
                  gsl_stats_sd (jerk_var->data, jerk_var->stride, jerk_var->size)
                  );
          fflush (stdout);
        }
      }
    }
  }

  if (TRUE)
  {
    printf ("# Number of knots in the model %u\n", ncm_model_len (NCM_MODEL (qpw)) - 2);
    fit = ncm_fit_new (NCM_FIT_TYPE_GSL_MMS, "NCM_FIT_GSL_MMS_NMSIMPLEX2", lh, mset, diff_algo);

    ncm_fit_run (fit, msg_level);
    ncm_fit_log_info (fit);
    ncm_fit_numdiff_m2lnL_covar (fit);
    ncm_fit_log_covar (fit);
  }

  //  ncm_fit_cr (lh, 0, 1, 0.954f);

  if (print_data && fit != NULL)
  {
    gdouble q = ncm_mset_param_get (fit->mset, NC_HICOSMO_ID, 0);
    gdouble pos = 0.0f;
    gdouble sigma_q = ncm_fit_covar_sd (fit, NC_HICOSMO_ID, 0);
    //    gdouble E = 1.0f;

    for (i = 1; i < ncm_model_len (NCM_MODEL (qpw)); i++)
    {
      gdouble step;
      gdouble J = ncm_mset_param_get (fit->mset, NC_HICOSMO_ID, i);
      gdouble sigma_j = ncm_fit_covar_sd (fit, NC_HICOSMO_ID, i);
      gdouble cov_0_i = ncm_fit_covar_cov (fit, NC_HICOSMO_ID, 0, NC_HICOSMO_ID, i);
      gboolean end = (i + 1 == ncm_model_len (NCM_MODEL (qpw)) ? TRUE : FALSE);
      for (j = 1; j < i; j++)
        cov_0_i += ncm_fit_covar_cov (fit, NC_HICOSMO_ID, i, NC_HICOSMO_ID, j) * interval;

      for (step = 0.0f; (step < interval) || (end && (step + pos) <= 1.8); step += 0.01)
      {
        gdouble sigma_q_z = sqrt(sigma_q*sigma_q + step*step*sigma_j*sigma_j + 2.0f * step * cov_0_i);
        printf ("%g %g %g\n", pos + step,
                nc_hicosmo_q (NC_HICOSMO (qpw), pos + step), sigma_q_z
                );
      }

      q += interval * J;
      sigma_q = sqrt(sigma_q*sigma_q + interval*interval*sigma_j*sigma_j + 2.0f * interval * cov_0_i);
      pos += interval;
    }
  }

  ncm_fit_clear (&fit);
  ncm_model_free (NCM_MODEL (qpw));
  ncm_model_free (NCM_MODEL (lcdm));
  ncm_mset_free (mset);
  ncm_mset_free (mset_lcdm);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  return 0;
}
