/***************************************************************************
 *            ncm_fit_gsl_mm.c
 *
 *  Mon Jun 11 12:08:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_fit_gsl_mm
 * @title: NcmFitGSLMM
 * @short_description: Best-fit finder -- GSL non-linear minimization algorithms.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_gsl_mm.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_ALGO,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmFitGSLMM, ncm_fit_gsl_mm, NCM_TYPE_FIT);

static void
ncm_fit_gsl_mm_init (NcmFitGSLMM *fit_gsl_mm)
{
  fit_gsl_mm->mm = NULL;
  fit_gsl_mm->algo = 0;
  fit_gsl_mm->desc = NULL;
}

static gdouble nc_residual_multimin_f (const gsl_vector *x, gpointer p);
static void nc_residual_multimin_df (const gsl_vector *x, gpointer p, gsl_vector *df);
static void nc_residual_multimin_fdf (const gsl_vector *x, gpointer p, gdouble *f, gsl_vector *df);

static void
_ncm_fit_gsl_mm_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_gsl_mm_parent_class)->constructed (object);
  {
    NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (object);
    NcmFit *fit = NCM_FIT (fit_gsl_mm);
    guint i;

    fit_gsl_mm->err_a = GSL_POSINF;
    fit_gsl_mm->err_b = 1.0e-1;

    for (i = 0; i < fit->fstate->fparam_len; i++)
    {
      gdouble pscale = ncm_mset_fparam_get_scale (fit->mset, i);
      fit_gsl_mm->err_a = GSL_MIN (fit_gsl_mm->err_a, pscale);
    }

    fit_gsl_mm->f.f      = &nc_residual_multimin_f;
    fit_gsl_mm->f.df     = &nc_residual_multimin_df;
    fit_gsl_mm->f.fdf    = &nc_residual_multimin_fdf;
    fit_gsl_mm->f.n      = fit->fstate->fparam_len;
    fit_gsl_mm->f.params = fit_gsl_mm;

    ncm_fit_gsl_mm_set_algo (fit_gsl_mm, fit_gsl_mm->algo);
  }
}

static void
_ncm_fit_gsl_mm_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (object);
  g_return_if_fail (NCM_IS_FIT_GSL_MM (object));

  switch (prop_id)
  {
    case PROP_ALGO:
    {
      if (fit_gsl_mm->mm == NULL)
        fit_gsl_mm->algo = g_value_get_enum (value);
      else
        ncm_fit_gsl_mm_set_algo (fit_gsl_mm, g_value_get_enum (value));
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_gsl_mm_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (object);
  g_return_if_fail (NCM_IS_FIT_GSL_MM (object));

  switch (prop_id)
  {
    case PROP_ALGO:
    g_value_set_enum (value, fit_gsl_mm->algo);
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_gsl_mm_finalize (GObject *object)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (object);

  gsl_multimin_fdfminimizer_free (fit_gsl_mm->mm);
  fit_gsl_mm->mm = NULL;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_gsl_mm_parent_class)->finalize (object);
}

static NcmFit *_ncm_fit_gsl_mm_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
static void _ncm_fit_gsl_mm_reset (NcmFit *fit);
static gboolean _ncm_fit_gsl_mm_run (NcmFit *fit, NcmFitRunMsgs mtype);
static const gchar *_ncm_fit_gsl_mm_get_desc (NcmFit *fit);

static void
ncm_fit_gsl_mm_class_init (NcmFitGSLMMClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitClass* fit_class     = NCM_FIT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_gsl_mm_constructed;
  object_class->set_property = &_ncm_fit_gsl_mm_set_property;
  object_class->get_property = &_ncm_fit_gsl_mm_get_property;
  object_class->finalize     = &ncm_fit_gsl_mm_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ALGO,
                                   g_param_spec_enum ("algorithm",
                                                      NULL,
                                                      "GSL multidimensional minimization algorithm",
                                                      NCM_TYPE_FIT_GSLMM_ALGOS, NCM_FIT_GSL_MM_VECTOR_BFGS2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  fit_class->copy_new = &_ncm_fit_gsl_mm_copy_new;
  fit_class->reset    = &_ncm_fit_gsl_mm_reset;
  fit_class->run      = &_ncm_fit_gsl_mm_run;
  fit_class->get_desc = &_ncm_fit_gsl_mm_get_desc;
}

static NcmFit *
_ncm_fit_gsl_mm_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (fit);
  return ncm_fit_gsl_mm_new (lh, mset, gtype, fit_gsl_mm->algo);
}

static void 
_ncm_fit_gsl_mm_reset (NcmFit *fit)
{
  /* Chain up : start */
  NCM_FIT_CLASS (ncm_fit_gsl_mm_parent_class)->reset (fit);
  {
    NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (fit);
    if (fit_gsl_mm->f.n != fit->fstate->fparam_len)
    {
      gsl_multimin_fdfminimizer_free (fit_gsl_mm->mm);
      fit_gsl_mm->mm = NULL;
      fit_gsl_mm->f.n = fit->fstate->fparam_len;
      ncm_fit_gsl_mm_set_algo (fit_gsl_mm, fit_gsl_mm->algo);
    }
  }
}

static gboolean
_ncm_fit_gsl_mm_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (fit);
  gint status;
  const gdouble prec = ncm_fit_get_params_reltol (fit);
  gdouble last_min   = GSL_POSINF;
  guint restart      = 10;
  const gdouble rfac = 0.99;

  if (ncm_fit_has_equality_constraints (fit) || ncm_fit_has_inequality_constraints (fit))
    g_error ("_ncm_fit_gsl_mm_run: GSL algorithms do not support constraints.");

  g_assert (fit->fstate->fparam_len != 0);
  
  ncm_mset_fparams_get_vector (fit->mset, fit->fstate->fparams);
  gsl_multimin_fdfminimizer_set (fit_gsl_mm->mm, &fit_gsl_mm->f, ncm_vector_gsl (fit->fstate->fparams), fit_gsl_mm->err_a, fit_gsl_mm->err_b);

  do
  {
    gdouble pscale;
    fit->fstate->niter++;
    status = gsl_multimin_fdfminimizer_iterate (fit_gsl_mm->mm);
    pscale = prec * fabs (fit_gsl_mm->mm->f != 0.0 ? fit_gsl_mm->mm->f : 1.0);

    if ((fit->fstate->niter == 1) && !gsl_finite(fit_gsl_mm->mm->f))
    {
      ncm_fit_params_set_vector (fit, fit->fstate->fparams);
      return FALSE;
    }

    if (status == GSL_ENOPROG)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_fit_log_step_error (fit, gsl_strerror (status));
      status = GSL_SUCCESS;
    }
    else
      status = gsl_multimin_test_gradient (fit_gsl_mm->mm->gradient, pscale);

    if ((restart > 0) && (status == GSL_SUCCESS))
    {
      if (fit_gsl_mm->mm->f < (last_min * rfac))
      {
        gsl_multimin_fdfminimizer_restart (fit_gsl_mm->mm);
        status = GSL_CONTINUE;
        restart--;
        last_min = fit_gsl_mm->mm->f;
      }
    }

    ncm_fit_state_set_m2lnL_curval (fit->fstate, fit_gsl_mm->mm->f);
    ncm_fit_log_step (fit);
  }
  while ( (status == GSL_CONTINUE) && (fit->fstate->niter < fit->maxiter) );

  ncm_mset_fparams_get_vector (fit->mset, fit->fstate->fparams);
  ncm_fit_state_set_m2lnL_curval (fit->fstate, fit_gsl_mm->mm->f);
  ncm_fit_state_set_m2lnL_prec (fit->fstate, gsl_blas_dnrm2 (fit_gsl_mm->mm->gradient) / fit_gsl_mm->mm->f);

  ncm_fit_params_set_gsl_vector (fit, fit_gsl_mm->mm->x);

  return TRUE;
}

static gdouble
nc_residual_multimin_f (const gsl_vector *x, gpointer p)
{
  NcmFit *fit = NCM_FIT (p);
  gdouble result;

  ncm_fit_params_set_gsl_vector (fit, x);
  if (!ncm_mset_params_valid (fit->mset))
    return GSL_EDOM;
  
  ncm_fit_m2lnL_val (fit, &result);
  return result;
}

static void
nc_residual_multimin_df (const gsl_vector *x, gpointer p, gsl_vector *df)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *dfv = ncm_vector_new_gsl_static (df);
  
  ncm_fit_params_set_gsl_vector (fit, x);
  if (!ncm_mset_params_valid (fit->mset))
    g_warning ("nc_residual_multimin_df: stepping in a invalid parameter point, continuing anyway.");

  ncm_fit_m2lnL_grad (fit, dfv);
  ncm_vector_free (dfv);
}

static void
nc_residual_multimin_fdf (const gsl_vector *x, gpointer p, gdouble *f, gsl_vector *df)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *dfv = ncm_vector_new_gsl_static (df);
  
  ncm_fit_params_set_gsl_vector (fit, x);
  if (!ncm_mset_params_valid (fit->mset))
    g_warning ("nc_residual_multimin_fdf: stepping in a invalid parameter point, continuing anyway.");

  ncm_fit_m2lnL_val_grad (fit, f, dfv);

  ncm_vector_free (dfv);
}

static const gchar *
_ncm_fit_gsl_mm_get_desc (NcmFit *fit)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (fit);
  if (fit_gsl_mm->desc == NULL)
  {
    fit_gsl_mm->desc = g_strdup_printf ("GSL Multidimensional Minimization:%s", 
                                        fit_gsl_mm->mm != NULL ? gsl_multimin_fdfminimizer_name (fit_gsl_mm->mm) : "not-set");
  }
  return fit_gsl_mm->desc;
}

/**
 * ncm_fit_gsl_mm_new:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 * @algo: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_gsl_mm_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitGSLMMAlgos algo)
{
  return g_object_new (NCM_TYPE_FIT_GSL_MM, 
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       "algorithm", algo,
                       NULL
                       );
}

/**
 * ncm_fit_gsl_mm_new_default:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_gsl_mm_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return g_object_new (NCM_TYPE_FIT_GSL_MM, 
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       NULL
                       );
}

/**
 * ncm_fit_gsl_mm_new_by_name:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 * @algo_name: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_gsl_mm_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name)
{
  if (algo_name != NULL)
  {
    const GEnumValue *algo = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_GSLMM_ALGOS,
                                                               algo_name);
    if (algo == NULL)
      g_error ("ncm_fit_gsl_mm_new_by_name: algorithm %s not found.", algo_name);
    return ncm_fit_gsl_mm_new (lh, mset, gtype, algo->value);
  }
  else
    return ncm_fit_gsl_mm_new_default (lh, mset, gtype);
}

/**
 * ncm_fit_gsl_mm_set_algo:
 * @fit_gsl_mm: a #NcmFitGSLMM.
 * @algo: a #gsl_mm_algorithm.
 *
 * FIXME
 *
 */
void
ncm_fit_gsl_mm_set_algo (NcmFitGSLMM *fit_gsl_mm, NcmFitGSLMMAlgos algo)
{
  const gsl_multimin_fdfminimizer_type *ncm_fit_gsl_mm_algos[] = {
    gsl_multimin_fdfminimizer_conjugate_fr,
    gsl_multimin_fdfminimizer_conjugate_pr,
    gsl_multimin_fdfminimizer_vector_bfgs,
    gsl_multimin_fdfminimizer_vector_bfgs2,
    gsl_multimin_fdfminimizer_steepest_descent,
  };
  NcmFit *fit = NCM_FIT (fit_gsl_mm);

  g_assert (fit_gsl_mm->algo < NCM_FIT_GSL_MM_NUM_ALGOS);

  if (fit_gsl_mm->algo != algo)
  {
    fit_gsl_mm->algo = algo;
    if (fit_gsl_mm->mm != NULL)
      gsl_multimin_fdfminimizer_free (fit_gsl_mm->mm);
    fit_gsl_mm->mm = NULL;
    if (fit_gsl_mm->desc != NULL)
      g_free (fit_gsl_mm->desc);
  }

  if (fit_gsl_mm->mm == NULL)
    fit_gsl_mm->mm = gsl_multimin_fdfminimizer_alloc (ncm_fit_gsl_mm_algos[fit_gsl_mm->algo], fit->fstate->fparam_len);
}
