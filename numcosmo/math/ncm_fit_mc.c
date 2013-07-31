/***************************************************************************
 *            ncm_fit_mc.c
 *
 *  Sat December 01 17:19:10 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * 
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

/**
 * SECTION:ncm_fit_mc
 * @title: Monte Carlo Analysis
 * @short_description: Object implementing Monte Carlo analysis
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_mc.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_statistics_double.h>

enum
{
  PROP_0,
  PROP_FIT
};

G_DEFINE_TYPE (NcmFitMC, ncm_fit_mc, G_TYPE_OBJECT);

static void
ncm_fit_mc_init (NcmFitMC *mc)
{
  mc->fit    = NULL;
  mc->fparam = NULL;
  mc->m2lnL  = NULL;
  mc->nt     = ncm_timer_new ();
  mc->n      = 0;
  mc->h      = NULL;
  mc->h_pdf  = NULL;
}

static void
ncm_fit_mc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMC *mc = NCM_FIT_MC (object);
  g_return_if_fail (NCM_IS_FIT_MC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      mc->fit = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_mc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitMC *mc = NCM_FIT_MC (object);
  g_return_if_fail (NCM_IS_FIT_MC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, mc->fit);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_mc_dispose (GObject *object)
{
  NcmFitMC *mc = NCM_FIT_MC (object);

  ncm_fit_clear (&mc->fit);
  ncm_stats_vec_clear (&mc->fparam);
  ncm_vector_clear (&mc->m2lnL);
  ncm_timer_clear (&mc->nt);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->dispose (object);
}

static void
ncm_fit_mc_finalize (GObject *object)
{
  NcmFitMC *mc = NCM_FIT_MC (object);

  if (mc->h != NULL)
    gsl_histogram_free (mc->h);
  if (mc->h_pdf != NULL)
    gsl_histogram_pdf_free (mc->h_pdf);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->finalize (object);
}

static void
ncm_fit_mc_class_init (NcmFitMCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_fit_mc_set_property;
  object_class->get_property = &ncm_fit_mc_get_property;
  object_class->dispose      = &ncm_fit_mc_dispose;
  object_class->finalize     = &ncm_fit_mc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_fit_mc_new:
 * @fit: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitMC *
ncm_fit_mc_new (NcmFit *fit)
{
  return g_object_new (NCM_TYPE_FIT_MC, 
                       "fit", fit,
                       NULL);
}

/**
 * ncm_fit_mc_free:
 * @mc: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_free (NcmFitMC *mc)
{
  g_object_unref (mc);
}

/**
 * ncm_fit_mc_clear:
 * @mc: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_clear (NcmFitMC **mc)
{
  g_clear_object (mc);
}

/**
 * ncm_fit_mc_run:
 * @mc: FIXME
 * @fiduc: FIXME
 * @ni: FIXME
 * @nf: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_run (NcmFitMC *mc, NcmMSet *fiduc, guint ni, guint nf, NcmFitRunMsgs mtype)
{
  const guint n = nf - ni;
  guint free_params_len = ncm_mset_fparams_len (mc->fit->mset);
  guint param_len = ncm_mset_param_len (mc->fit->mset);
  NcmVector *bf = ncm_vector_new (param_len);
  gint i;

  g_assert (nf > ni);
  g_assert (ncm_mset_cmp (mc->fit->mset, fiduc, FALSE));

  mc->n = n;
  mc->m2lnL_min = GSL_POSINF;
  mc->m2lnL_max = GSL_NEGINF;

  ncm_stats_vec_clear (&mc->fparam);
  mc->fparam = ncm_stats_vec_new (free_params_len, NCM_STATS_VEC_COV, TRUE);
  ncm_vector_clear (&mc->m2lnL);
  mc->m2lnL = ncm_vector_new (n);

  ncm_mset_param_get_vector (mc->fit->mset, bf);
  
  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# Starting Montecarlo...\n");
    ncm_dataset_log_info (mc->fit->lh->dset);
    ncm_cfg_msg_sepa ();
    g_message ("# Fiducial model set:\n");
    ncm_mset_pretty_log (fiduc);
    ncm_cfg_msg_sepa ();
    g_message ("# Fitting model set:\n");
    ncm_mset_pretty_log (mc->fit->mset);

    ncm_cfg_msg_sepa ();
    g_message ("# Calculating [%06d] Montecarlo fits\n", n);
    if (ni > 0)
      g_message ("# Resampling %u-times\n#", ni);
  }
  
  for (i = 0; i < ni; i++)
  {
    ncm_dataset_resample (mc->fit->lh->dset, fiduc);
    if (mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message (".");
  }
  if (ni > 0 && mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("\n");

  ncm_timer_task_start (mc->nt, nf - ni);
  ncm_timer_set_name (mc->nt, "Montecarlo");
  if (mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (mc->nt);
  
  for (; i < nf; i++)
  {
    NcmVector *bf_i = ncm_stats_vec_peek_x (mc->fparam);
    NcmFitRunMsgs mcrun_msg = (mtype == NCM_FIT_RUN_MSGS_FULL) ? NCM_FIT_RUN_MSGS_SIMPLE : NCM_FIT_RUN_MSGS_NONE;
    gdouble m2lnL;

    ncm_mset_param_set_vector (mc->fit->mset, bf);
    ncm_dataset_resample (mc->fit->lh->dset, fiduc);
    
    ncm_fit_run (mc->fit, mcrun_msg);
    
    ncm_mset_fparams_get_vector (mc->fit->mset, bf_i);
    m2lnL = ncm_fit_state_get_m2lnL_curval (mc->fit->fstate);
    ncm_vector_set (mc->m2lnL, i - ni, m2lnL);

    //ncm_vector_log_vals (bf_i, "# Last bestfit: ", "% 12.5g");
    ncm_stats_vec_update (mc->fparam);
    g_message ("# Improvements: mean = % 8.7f, var = % 8.7f, cov = % 8.7f\n", mc->fparam->mean_inc, mc->fparam->var_inc, mc->fparam->cov_inc);
    ncm_vector_log_vals (mc->fparam->mean, "# Current mean: ", "% 12.5g");
    ncm_vector_log_vals (mc->fparam->var, "# Current var:  ", "% 12.5g");

    mc->m2lnL_min = GSL_MIN (mc->m2lnL_min, m2lnL);
    mc->m2lnL_max = GSL_MAX (mc->m2lnL_max, m2lnL);

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_timer_task_increment (mc->nt);
      ncm_timer_task_log_elapsed (mc->nt);
      ncm_timer_task_log_mean_time (mc->nt);
      ncm_timer_task_log_time_left (mc->nt);
      ncm_timer_task_log_end_datetime (mc->nt);
    }
  }
  ncm_timer_task_end (mc->nt);
  ncm_vector_free (bf);
}

/**
 * ncm_fit_mc_print:
 * @mc: FIXME
 *
 * FIXME
 */
void 
ncm_fit_mc_print (NcmFitMC *mc)
{
  guint i;
  guint len = mc->fparam->saved_x->len;

  for (i = 0; i < len; i++)
  {
    gint j;
    NcmVector *x_i = g_ptr_array_index (mc->fparam->saved_x, i);
    for (j = 0; j < ncm_vector_len (x_i); j++)
      g_message ("% -20.15g ", ncm_vector_get (x_i, j));
    g_message ("\n");
  }
}

/**
 * ncm_fit_mc_mean_covar:
 * @mc: FIXME
 *
 * FIXME
 */
void
ncm_fit_mc_mean_covar (NcmFitMC *mc)
{
  ncm_stats_vec_get_cov_matrix (mc->fparam, mc->fit->fstate->covar);
  ncm_stats_vec_get_mean_vector (mc->fparam, mc->fit->fstate->fparams);
  ncm_mset_fparams_set_vector (mc->fit->mset, mc->fit->fstate->fparams);
  mc->fit->fstate->has_covar = TRUE;
}

/**
 * ncm_fit_mc_gof_pdf:
 * @mc: FIXME
 *
 * FIXME
 * 
 */
void
ncm_fit_mc_gof_pdf (NcmFitMC *mc)
{
  const guint nbins = mc->n / 10 >= 10 ? mc->n / 10 : 10;
  gint i;
	
  if (mc->h != NULL && mc->h->n != nbins)
  {
    gsl_histogram_free (mc->h);
    mc->h = NULL;
  }
  if (mc->h == NULL)
    mc->h = gsl_histogram_alloc (nbins);

  if (mc->h_pdf != NULL && mc->h_pdf->n != nbins)
  {
    gsl_histogram_pdf_free (mc->h_pdf);
    mc->h_pdf = NULL;
  }
  if (mc->h_pdf == NULL)
    mc->h_pdf = gsl_histogram_pdf_alloc (nbins);

  gsl_histogram_set_ranges_uniform (mc->h, mc->m2lnL_min, mc->m2lnL_max);

  for (i = 0; i < mc->n; i++)
    gsl_histogram_increment (mc->h, ncm_vector_get (mc->m2lnL, i));

  gsl_histogram_pdf_init (mc->h_pdf, mc->h);  
}

/**
 * ncm_fit_mc_gof_pdf_print:
 * @mc: FIXME
 *
 * FIXME
 * 
 */
void
ncm_fit_mc_gof_pdf_print (NcmFitMC *mc)
{
  g_assert (mc->h != NULL);

  ncm_cfg_msg_sepa ();
  g_message ("# Monte Carlo Goodness of fit.\n");
  g_message ("#  m2lnL distribution:\n#\n");
  g_message ("#   - interval [% 20.15g % 20.15g],\n", mc->m2lnL_min, mc->m2lnL_max);
  g_message ("#   - mean               % 20.15g,\n", gsl_histogram_mean (mc->h));
  g_message ("#   - standard deviation % 20.15g.\n#\n", gsl_histogram_sigma (mc->h));
  
}

/**
 * ncm_fit_mc_gof_pdf_pvalue:
 * @mc: FIXME
 * @m2lnL: FIXME
 * @both: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble
ncm_fit_mc_gof_pdf_pvalue (NcmFitMC *mc, gdouble m2lnL, gboolean both)
{
  gsize i;
  g_assert (mc->h_pdf != NULL);

  if (m2lnL < mc->m2lnL_min || m2lnL > mc->m2lnL_max)
  {
    g_warning ("ncm_fit_mc_gof_pdf_pvalue: value % 20.15g outside mc obtained interval [% 20.15g % 20.15g]. Assuming 0 pvalue.",
               m2lnL, mc->m2lnL_min, mc->m2lnL_max);
    return 0.0;
  }

  gsl_histogram_find (mc->h, m2lnL, &i);
  g_assert_cmpint (i, <=, mc->h_pdf->n);
  
  if (i == 0)
    return 1.0;
  else
    return (1.0 - mc->h_pdf->sum[i - 1]);
}
