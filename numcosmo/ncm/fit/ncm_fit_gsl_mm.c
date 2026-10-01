/***************************************************************************
 *            ncm_fit_gsl_mm.c
 *
 *  Mon Jun 11 12:08:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmFitGSLMM:
 *
 * Best-fit finder using the GSL gradient-based minimizers (gsl_multimin_fdfminimizer).
 *
 * It minimizes $-2\ln L$ with the gradient from ncm_fit_m2lnL_grad(). The first trial
 * step is the smallest free-parameter scale. At a point the models report invalid
 * (ncm_mset_params_valid()) $-2\ln L$ is $+\infty$, so the line minimization
 * backtracks; parameter bounds are not enforced. Constraints are not supported and
 * abort the run.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_fit_gsl_mm.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_ALGO,
  PROP_SIZE,
};

struct _NcmFitGSLMM
{
  /*< private >*/
  NcmFit parent_instance;
  gsl_multimin_fdfminimizer *mm;
  gsl_multimin_function_fdf f;
  NcmFitGSLMMAlgos algo;
  gchar *desc;
  gdouble err_a;
  gdouble err_b;
};

G_DEFINE_TYPE (NcmFitGSLMM, ncm_fit_gsl_mm, NCM_TYPE_FIT)

static void
ncm_fit_gsl_mm_init (NcmFitGSLMM *fit_gsl_mm)
{
  fit_gsl_mm->mm    = NULL;
  fit_gsl_mm->algo  = 0;
  fit_gsl_mm->desc  = NULL;
  fit_gsl_mm->err_a = GSL_POSINF;
  fit_gsl_mm->err_b = GSL_POSINF;
}

static gdouble nc_residual_multimin_f (const gsl_vector *x, gpointer p);
static void nc_residual_multimin_df (const gsl_vector *x, gpointer p, gsl_vector *df);
static void nc_residual_multimin_fdf (const gsl_vector *x, gpointer p, gdouble *f, gsl_vector *df);
static void _ncm_fit_gsl_mm_prepare (NcmFitGSLMM *fit_gsl_mm);

static void
_ncm_fit_gsl_mm_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_gsl_mm_parent_class)->constructed (object);
  {
    NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (object);

    fit_gsl_mm->err_b = 1.0e-1;

    fit_gsl_mm->f.f      = &nc_residual_multimin_f;
    fit_gsl_mm->f.df     = &nc_residual_multimin_df;
    fit_gsl_mm->f.fdf    = &nc_residual_multimin_fdf;
    fit_gsl_mm->f.n      = 0;
    fit_gsl_mm->f.params = fit_gsl_mm;

    _ncm_fit_gsl_mm_prepare (fit_gsl_mm);
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
      ncm_fit_gsl_mm_set_algo (fit_gsl_mm, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_fit_gsl_mm_finalize (GObject *object)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (object);

  g_clear_pointer (&fit_gsl_mm->mm, gsl_multimin_fdfminimizer_free);
  g_clear_pointer (&fit_gsl_mm->desc, g_free);

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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmFitClass *fit_class     = NCM_FIT_CLASS (klass);

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

  _ncm_fit_gsl_mm_prepare (NCM_FIT_GSL_MM (fit));
}

static const gsl_multimin_fdfminimizer_type *_ncm_fit_gsl_mm_type (NcmFitGSLMMAlgos algo);

/* Sets the first trial step to the smallest free-parameter scale and allocates the
 * minimizer for the current number of free parameters; without free parameters there
 * is none. */
static void
_ncm_fit_gsl_mm_prepare (NcmFitGSLMM *fit_gsl_mm)
{
  NcmFit *fit            = NCM_FIT (fit_gsl_mm);
  NcmMSet *mset          = ncm_fit_peek_mset (fit);
  const guint fparam_len = ncm_fit_state_get_fparam_len (ncm_fit_peek_state (fit));
  guint i;

  fit_gsl_mm->err_a = GSL_POSINF;

  for (i = 0; i < fparam_len; i++)
    fit_gsl_mm->err_a = GSL_MIN (fit_gsl_mm->err_a, ncm_mset_fparam_get_scale (mset, i));

  if (fit_gsl_mm->f.n != fparam_len)
  {
    g_clear_pointer (&fit_gsl_mm->mm, gsl_multimin_fdfminimizer_free);
    fit_gsl_mm->f.n = fparam_len;
  }

  if ((fit_gsl_mm->mm == NULL) && (fparam_len > 0))
    fit_gsl_mm->mm = gsl_multimin_fdfminimizer_alloc (_ncm_fit_gsl_mm_type (fit_gsl_mm->algo), fparam_len);
}

static const gsl_multimin_fdfminimizer_type *
_ncm_fit_gsl_mm_type (NcmFitGSLMMAlgos algo)
{
  switch (algo)
  {
    case NCM_FIT_GSL_MM_CONJUGATE_FR:
      return gsl_multimin_fdfminimizer_conjugate_fr;

    case NCM_FIT_GSL_MM_CONJUGATE_PR:
      return gsl_multimin_fdfminimizer_conjugate_pr;

    case NCM_FIT_GSL_MM_VECTOR_BFGS:
      return gsl_multimin_fdfminimizer_vector_bfgs;

    case NCM_FIT_GSL_MM_VECTOR_BFGS2:
      return gsl_multimin_fdfminimizer_vector_bfgs2;

    case NCM_FIT_GSL_MM_STEEPEST_DESCENT:
      return gsl_multimin_fdfminimizer_steepest_descent;

    default:                                                         /* LCOV_EXCL_LINE */
      g_error ("_ncm_fit_gsl_mm_type: unknown algorithm %d.", algo); /* LCOV_EXCL_LINE */

      return NULL; /* LCOV_EXCL_LINE */
  }
}

static gboolean
_ncm_fit_gsl_mm_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (fit);
  NcmFitState *fstate     = ncm_fit_peek_state (fit);
  NcmMSet *mset           = ncm_fit_peek_mset (fit);
  const gdouble prec      = ncm_fit_get_params_reltol (fit);
  gdouble last_min        = GSL_POSINF;
  guint restart           = 10;
  const gdouble rfac      = 0.99;
  gint status;

  if (ncm_fit_equality_constraints_len (fit) || ncm_fit_inequality_constraints_len (fit))
    g_error ("_ncm_fit_gsl_mm_run: GSL algorithms do not support constraints.");

  g_assert (ncm_fit_state_get_fparam_len (fstate) != 0);

  ncm_mset_fparams_get_vector (mset, ncm_fit_state_peek_fparams (fstate));
  gsl_multimin_fdfminimizer_set (fit_gsl_mm->mm, &fit_gsl_mm->f, ncm_vector_gsl (ncm_fit_state_peek_fparams (fstate)), fit_gsl_mm->err_a, fit_gsl_mm->err_b);

  do {
    gdouble pscale;

    ncm_fit_state_add_iter (fstate, 1);

    status = gsl_multimin_fdfminimizer_iterate (fit_gsl_mm->mm);
    pscale = prec * fabs (fit_gsl_mm->mm->f != 0.0 ? fit_gsl_mm->mm->f : 1.0);

    if ((ncm_fit_state_get_niter (fstate) == 1) && !gsl_finite (fit_gsl_mm->mm->f))
    {
      ncm_fit_params_set_vector (fit, ncm_fit_state_peek_fparams (fstate));

      return FALSE;
    }

    /* No progress means the line minimization cannot lower -2 ln L: converged. */
    if (status == GSL_ENOPROG)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_fit_log_step_error (fit, gsl_strerror (status));

      status = GSL_SUCCESS;
    }
    else if (status != GSL_SUCCESS)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_fit_log_step_error (fit, gsl_strerror (status));

      break;
    }
    else
    {
      status = gsl_multimin_test_gradient (fit_gsl_mm->mm->gradient, pscale);
    }

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

    ncm_fit_state_set_m2lnL_curval (fstate, fit_gsl_mm->mm->f);
    ncm_fit_log_step (fit);
  } while ((status == GSL_CONTINUE) && (ncm_fit_state_get_niter (fstate) < ncm_fit_get_maxiter (fit)));

  ncm_fit_params_set_gsl_vector (fit, fit_gsl_mm->mm->x);
  ncm_mset_fparams_get_vector (mset, ncm_fit_state_peek_fparams (fstate));
  ncm_fit_state_set_m2lnL_curval (fstate, fit_gsl_mm->mm->f);
  ncm_fit_state_set_m2lnL_prec (fstate, fabs (gsl_blas_dnrm2 (fit_gsl_mm->mm->gradient) / fit_gsl_mm->mm->f));

  return (status == GSL_SUCCESS);
}

static gdouble
nc_residual_multimin_f (const gsl_vector *x, gpointer p)
{
  NcmFit *fit   = NCM_FIT (p);
  NcmMSet *mset = ncm_fit_peek_mset (fit);
  gdouble result;

  ncm_fit_params_set_gsl_vector (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_POSINF;

  ncm_fit_m2lnL_val (fit, &result);

  return result;
}

static void
nc_residual_multimin_df (const gsl_vector *x, gpointer p, gsl_vector *df)
{
  NcmFit *fit    = NCM_FIT (p);
  NcmMSet *mset  = ncm_fit_peek_mset (fit);
  NcmVector *dfv = ncm_vector_new_gsl_static (df);

  ncm_fit_params_set_gsl_vector (fit, x);

  if (!ncm_mset_params_valid (mset))
    g_warning ("nc_residual_multimin_df: stepping in a invalid parameter point, continuing anyway.");

  ncm_fit_m2lnL_grad (fit, dfv);
  ncm_vector_free (dfv);
}

static void
nc_residual_multimin_fdf (const gsl_vector *x, gpointer p, gdouble *f, gsl_vector *df)
{
  NcmFit *fit    = NCM_FIT (p);
  NcmMSet *mset  = ncm_fit_peek_mset (fit);
  NcmVector *dfv = ncm_vector_new_gsl_static (df);

  ncm_fit_params_set_gsl_vector (fit, x);

  /* The line minimization backtracks from an infinite value. */
  if (!ncm_mset_params_valid (mset))
  {
    *f = GSL_POSINF;
    gsl_vector_set_zero (df);
    ncm_vector_free (dfv);

    return;
  }

  ncm_fit_m2lnL_val_grad (fit, f, dfv);

  ncm_vector_free (dfv);
}

static const gchar *
_ncm_fit_gsl_mm_get_desc (NcmFit *fit)
{
  NcmFitGSLMM *fit_gsl_mm = NCM_FIT_GSL_MM (fit);

  if (fit_gsl_mm->desc == NULL)
    fit_gsl_mm->desc = g_strdup_printf ("GSL Multidimensional Minimization:%s",
                                        _ncm_fit_gsl_mm_type (fit_gsl_mm->algo)->name);

  return fit_gsl_mm->desc;
}

/**
 * ncm_fit_gsl_mm_new:
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 * @algo: a #NcmFitGSLMMAlgos
 *
 * Creates a #NcmFitGSLMM for @lh and @mset using @algo, with gradients computed as
 * @gtype says.
 *
 * Returns: (transfer full): a new #NcmFitGSLMM.
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
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 *
 * Creates a #NcmFitGSLMM as ncm_fit_gsl_mm_new() with #NCM_FIT_GSL_MM_VECTOR_BFGS2.
 *
 * Returns: (transfer full): a new #NcmFitGSLMM.
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
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 * @algo_name: (nullable): name or nick of a #NcmFitGSLMMAlgos
 *
 * Creates a #NcmFitGSLMM as ncm_fit_gsl_mm_new() with the algorithm named @algo_name,
 * or #NCM_FIT_GSL_MM_VECTOR_BFGS2 when @algo_name is %NULL. An unknown name aborts.
 *
 * Returns: (transfer full): a new #NcmFitGSLMM.
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
  {
    return ncm_fit_gsl_mm_new_default (lh, mset, gtype);
  }
}

/**
 * ncm_fit_gsl_mm_set_algo:
 * @fit_gsl_mm: a #NcmFitGSLMM
 * @algo: a #NcmFitGSLMMAlgos
 *
 * Sets the minimization algorithm of @fit_gsl_mm to @algo.
 *
 */
void
ncm_fit_gsl_mm_set_algo (NcmFitGSLMM *fit_gsl_mm, NcmFitGSLMMAlgos algo)
{
  g_assert_cmpint (algo, <, NCM_FIT_GSL_MM_NUM_ALGOS);

  if (fit_gsl_mm->algo != algo)
  {
    fit_gsl_mm->algo = algo;

    g_clear_pointer (&fit_gsl_mm->mm, gsl_multimin_fdfminimizer_free);
    g_clear_pointer (&fit_gsl_mm->desc, g_free);
  }

  if ((fit_gsl_mm->mm == NULL) && (fit_gsl_mm->f.n > 0))
    fit_gsl_mm->mm = gsl_multimin_fdfminimizer_alloc (_ncm_fit_gsl_mm_type (fit_gsl_mm->algo), fit_gsl_mm->f.n);
}

