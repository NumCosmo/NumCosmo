/***************************************************************************
 *            ncm_fit_state.c
 *
 *  Thu November 29 15:27:11 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_fit_state
 * @title: NcmFitState
 * @short_description: State of a NcmFit object.
 *
 * Object that stores the state of a NcmFit object.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit.h"
#include "math/ncm_cfg.h"
#include "ncm_fit_state.h"

enum
{
  PROP_0,
  PROP_DATA_LEN,
  PROP_FPARAM_LEN,
  PROP_DOF,
  PROP_IS_LEAST_SQUARES,
  PROP_NITERS,
  PROP_FUNC_EVAL,
  PROP_GRAD_EVAL,
  PROP_IS_BEST_FIT,
};

struct _NcmFitState
{
  /*< private >*/
  GObject parent_instance;
  guint data_len;
  guint fparam_len;
  guint alloc_data_len;
  guint alloc_fparam_len;
  gint dof;
  guint niter;
  guint func_eval;
  guint grad_eval;
  gdouble m2lnL_prec;
  gdouble params_prec;
  gdouble elapsed_time;
  gdouble m2lnL_curval;
  NcmVector *dm2lnL;
  NcmVector *fparams;
  NcmVector *ls_f;
  NcmMatrix *ls_J;
  NcmMatrix *covar;
  NcmMatrix *hessian;
  gboolean is_best_fit;
  gboolean is_least_squares;
  gboolean has_covar;
};


G_DEFINE_TYPE (NcmFitState, ncm_fit_state, G_TYPE_OBJECT);

static void
ncm_fit_state_init (NcmFitState *fstate)
{
  fstate->data_len     = 0;
  fstate->fparam_len   = 0;
  fstate->covar        = NULL;
  fstate->hessian      = NULL;
  fstate->fparams      = NULL;
  fstate->dm2lnL       = NULL;
  fstate->ls_J         = NULL;
  fstate->m2lnL_curval = 0.0;
  fstate->niter        = 0;
  fstate->func_eval    = 0;
  fstate->grad_eval    = 0;
  fstate->m2lnL_prec   = 0.0;
  fstate->params_prec  = 0.0;

  fstate->alloc_data_len   = 0;
  fstate->alloc_fparam_len = 0;

  fstate->is_best_fit      = FALSE;
  fstate->is_least_squares = FALSE;
  fstate->has_covar        = FALSE;
}

static void
ncm_fit_state_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_state_parent_class)->constructed (object);
  {
    NcmFitState *fstate = NCM_FIT_STATE (object);

    ncm_fit_state_realloc (fstate);
  }
}

static void
ncm_fit_state_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitState *fstate = NCM_FIT_STATE (object);

  g_return_if_fail (NCM_IS_FIT_STATE (object));

  switch (prop_id)
  {
    case PROP_DATA_LEN:
      fstate->data_len = g_value_get_uint (value);
      ncm_fit_state_realloc (fstate);
      break;
    case PROP_FPARAM_LEN:
      fstate->fparam_len = g_value_get_uint (value);
      ncm_fit_state_realloc (fstate);
      break;
    case PROP_DOF:
      fstate->dof = g_value_get_int (value);
      break;
    case PROP_IS_LEAST_SQUARES:
      fstate->is_least_squares = g_value_get_boolean (value);
      ncm_fit_state_realloc (fstate);
      break;
    case PROP_NITERS:
      fstate->niter = g_value_get_uint (value);
      break;
    case PROP_FUNC_EVAL:
      fstate->func_eval = g_value_get_uint (value);
      break;
    case PROP_GRAD_EVAL:
      fstate->grad_eval = g_value_get_uint (value);
      break;
    case PROP_IS_BEST_FIT:
      fstate->is_best_fit = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_state_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitState *fstate = NCM_FIT_STATE (object);

  g_return_if_fail (NCM_IS_FIT_STATE (object));

  switch (prop_id)
  {
    case PROP_DATA_LEN:
      g_value_set_uint (value, fstate->data_len);
      break;
    case PROP_FPARAM_LEN:
      g_value_set_uint (value, fstate->fparam_len);
      break;
    case PROP_DOF:
      g_value_set_int (value, fstate->dof);
      break;
    case PROP_IS_LEAST_SQUARES:
      g_value_set_boolean (value, fstate->is_least_squares);
      break;
    case PROP_NITERS:
      g_value_set_uint (value, fstate->niter);
      break;
    case PROP_FUNC_EVAL:
      g_value_set_uint (value, fstate->func_eval);
      break;
    case PROP_GRAD_EVAL:
      g_value_set_uint (value, fstate->grad_eval);
      break;
    case PROP_IS_BEST_FIT:
      g_value_set_boolean (value, fstate->is_best_fit);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_state_dispose (GObject *object)
{
  NcmFitState *fstate = NCM_FIT_STATE (object);

  ncm_vector_clear (&fstate->fparams);
  ncm_vector_clear (&fstate->dm2lnL);
  ncm_vector_clear (&fstate->ls_f);
  ncm_matrix_clear (&fstate->covar);
  ncm_matrix_clear (&fstate->hessian);
  ncm_matrix_clear (&fstate->ls_J);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_state_parent_class)->dispose (object);
}

static void
ncm_fit_state_class_init (NcmFitStateClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &ncm_fit_state_constructed;
  object_class->set_property = &ncm_fit_state_set_property;
  object_class->get_property = &ncm_fit_state_get_property;
  object_class->dispose      = &ncm_fit_state_dispose;

  g_object_class_install_property (object_class,
                                   PROP_DATA_LEN,
                                   g_param_spec_uint ("data-len",
                                                      NULL,
                                                      "Data length",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_FPARAM_LEN,
                                   g_param_spec_uint ("fparam-len",
                                                      NULL,
                                                      "Free parameters length",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DOF,
                                   g_param_spec_int ("dof",
                                                     NULL,
                                                     "Degrees of freedom",
                                                     G_MININT, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_IS_LEAST_SQUARES,
                                   g_param_spec_boolean ("is-least-squares",
                                                         NULL,
                                                         "Is a least squares fit state",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NITERS,
                                   g_param_spec_uint ("niters",
                                                      NULL,
                                                      "Number of interations",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_FUNC_EVAL,
                                   g_param_spec_uint ("func-eval",
                                                      NULL,
                                                      "Number of function evaluations",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_GRAD_EVAL,
                                   g_param_spec_uint ("grad-eval",
                                                      NULL,
                                                      "Number of gradient evaluations",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_IS_BEST_FIT,
                                   g_param_spec_boolean ("is-best-fit",
                                                         NULL,
                                                         "Is a best fit state",
                                                         FALSE,
                                                         G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_fit_state_new:
 * @data_len: FIXME
 * @fparam_len: FIXME
 * @dof: FIXME
 * @is_least_squares: FIXME
 *
 * Instantiates a new #NcmFitState.
 *
 * Returns: (transfer full): a newly allocated #NcmFitState.
 */
NcmFitState *
ncm_fit_state_new (guint data_len, guint fparam_len, gint dof, gboolean is_least_squares)
{
  return g_object_new (NCM_TYPE_FIT_STATE,
                       "data-len", data_len,
                       "fparam-len", fparam_len,
                       "dof", dof,
                       "is-least-squares", is_least_squares,
                       NULL);
}

/**
 * ncm_fit_state_ref:
 * @fstate: a #NcmFitState
 *
 * Increases the reference count of @fstate by one.
 *
 * Returns: (transfer full): @fstate.
 */
NcmFitState *
ncm_fit_state_ref (NcmFitState *fstate)
{
  return g_object_ref (fstate);
}

/**
 * ncm_fit_state_free:
 * @fstate: a #NcmFitState
 *
 * Decreases the reference count of @fstate by one.
 *
 */
void
ncm_fit_state_free (NcmFitState *fstate)
{
  g_object_unref (fstate);
}

/**
 * ncm_fit_state_clear:
 * @fstate: a #NcmFitState
 *
 * Decreases the reference count of *@fstate by one, and sets the pointer *@fstate to NULL.
 *
 */
void
ncm_fit_state_clear (NcmFitState **fstate)
{
  g_clear_object (fstate);
}

/**
 * ncm_fit_state_set_all:
 * @fstate: a #NcmFitState
 * @data_len: Number of data points
 * @fparam_len: Number of free parameters
 * @dof: Degrees of freedom
 * @is_least_squares: whether it is a least squares fit
 *
 * Sets all the properties of @fstate.
 *
 */
void
ncm_fit_state_set_all (NcmFitState *fstate, guint data_len, guint fparam_len, gint dof, gboolean is_least_squares)
{
  fstate->data_len         = data_len;
  fstate->fparam_len       = fparam_len;
  fstate->dof              = dof;
  fstate->is_least_squares = is_least_squares;

  ncm_fit_state_realloc (fstate);
}

/**
 * ncm_fit_state_reset:
 * @fstate: a #NcmFitState
 *
 * Resets the #NcmFitState to its initial state.
 *
 */
void
ncm_fit_state_reset (NcmFitState *fstate)
{
  fstate->niter        = 0;
  fstate->func_eval    = 0;
  fstate->grad_eval    = 0;
  fstate->elapsed_time = 0.0;
  fstate->m2lnL_prec   = 0.0;
  fstate->params_prec  = 0.0;
  fstate->m2lnL_curval = 0.0;

  fstate->is_best_fit = FALSE;
}

/**
 * ncm_fit_state_realloc:
 * @fstate: a #NcmFitState
 *
 * Reallocates the #NcmFitState to its current state.
 *
 */
void
ncm_fit_state_realloc (NcmFitState *fstate)
{
  gboolean fparam_diff_len = FALSE;

  if (fstate->data_len == 0)
    return;

  if (fstate->alloc_fparam_len != fstate->fparam_len)
  {
    ncm_vector_clear (&fstate->fparams);
    ncm_vector_clear (&fstate->dm2lnL);
    ncm_matrix_clear (&fstate->covar);
    ncm_matrix_clear (&fstate->hessian);
    fparam_diff_len = TRUE;

    if (fstate->fparam_len > 0)
    {
      fstate->covar   = ncm_matrix_new (fstate->fparam_len, fstate->fparam_len);
      fstate->hessian = ncm_matrix_new (fstate->fparam_len, fstate->fparam_len);
      fstate->fparams = ncm_vector_new (fstate->fparam_len);
      fstate->dm2lnL  = ncm_vector_new (fstate->fparam_len);
    }

    fstate->alloc_fparam_len = fstate->fparam_len;
  }

  if (fstate->is_least_squares)
  {
    gboolean data_diff_len = fstate->alloc_data_len != fstate->data_len;

    if ((fparam_diff_len || data_diff_len || (fstate->ls_J == NULL)) && (fstate->fparam_len > 0))
    {
      ncm_matrix_clear (&fstate->ls_J);
      fstate->ls_J = ncm_matrix_new (fstate->data_len, fstate->fparam_len);
    }

    if (data_diff_len || (fstate->ls_f == NULL))
    {
      ncm_vector_clear (&fstate->ls_f);
      fstate->ls_f = ncm_vector_new (fstate->data_len);
    }
  }
  else
  {
    ncm_matrix_clear (&fstate->ls_J);
    ncm_vector_clear (&fstate->ls_f);
  }

  fstate->alloc_data_len = fstate->data_len;
}

/**
 * ncm_fit_state_set_ls:
 * @fstate: a #NcmFitState
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * Sets the least squares data of @fstate.
 *
 * This method is used by #NcmFit implementations to set the precision of the parameters.
 * It should not be used by the user.
 *
 */
void
ncm_fit_state_set_ls (NcmFitState *fstate, NcmVector *f, NcmMatrix *J)
{
  g_assert (fstate->is_least_squares);

  fstate->m2lnL_curval = ncm_vector_dnrm2 (f);

  ncm_vector_memcpy (fstate->ls_f, f);

  gsl_blas_dgemv (CblasTrans, 2.0, ncm_matrix_gsl (fstate->ls_J),
                  ncm_vector_gsl (fstate->ls_f), 0.0,
                  ncm_vector_gsl (fstate->dm2lnL));

  fstate->m2lnL_prec = sqrt (ncm_vector_dnrm2 (fstate->dm2lnL));

  if (fabs (fstate->m2lnL_curval) > 1.0e-3)
    fstate->m2lnL_prec = fabs (fstate->m2lnL_prec / fstate->m2lnL_curval);

  ncm_matrix_memcpy (fstate->ls_J, J);
}

/**
 * ncm_fit_state_set_fparam_len:
 * @fstate: a #NcmFitState
 * @fparam_len: Number of free parameters
 *
 * Sets the number of free parameters of @fstate.
 *
 * This method is used by #NcmFit implementations to set @fstate state.
 * It should not be used by the user.
 */
void
ncm_fit_state_set_fparam_len (NcmFitState *fstate, guint fparam_len)
{
  fstate->fparam_len = fparam_len;
}

/**
 * ncm_fit_state_get_fparam_len:
 * @fstate: a #NcmFitState
 *
 * Gets the number of free parameters of @fstate.
 *
 * Returns: Number of free parameters
 */
guint
ncm_fit_state_get_fparam_len (NcmFitState *fstate)
{
  return fstate->fparam_len;
}

/**
 * ncm_fit_state_set_niter:
 * @fstate: a #NcmFitState
 * @niter: Number of iterations
 *
 * Sets the number of iterations of @fstate.
 *
 * This method is used by #NcmFit implementations to set @fstate state.
 * It should not be used by the user.
 *
 */
void
ncm_fit_state_set_niter (NcmFitState *fstate, guint niter)
{
  fstate->niter = niter;
}

/**
 * ncm_fit_state_get_niter:
 * @fstate: a #NcmFitState
 *
 * Gets the number of iterations of @fstate.
 *
 * Returns: Number of iterations
 */
guint
ncm_fit_state_get_niter (NcmFitState *fstate)
{
  return fstate->niter;
}

/**
 * ncm_fit_state_set_m2lnL_prec:
 * @fstate: a #NcmFitState
 * @prec: Precision of the m2lnL
 *
 * Sets the precision of the m2lnL of @fstate.
 *
 * This method is used by #NcmFit implementations to set @fstate state.
 * It should not be used by the user.
 *
 */
void
ncm_fit_state_set_m2lnL_prec (NcmFitState *fstate, gdouble prec)
{
  fstate->m2lnL_prec = prec;
}

/**
 * ncm_fit_state_get_m2lnL_prec:
 * @fstate: a #NcmFitState
 *
 * Gets the precision of the m2lnL of @fstate.
 *
 * Returns: Precision of the m2lnL
 */
gdouble
ncm_fit_state_get_m2lnL_prec (NcmFitState *fstate)
{
  return fstate->m2lnL_prec;
}

/**
 * ncm_fit_state_set_m2lnL_curval:
 * @fstate: a #NcmFitState
 * @m2lnL: Current value of the m2lnL
 *
 * Sets the current value of the m2lnL of @fstate.
 *
 * This method is used by #NcmFit implementations to set @fstate state.
 * It should not be used by the user.
 *
 */
void
ncm_fit_state_set_m2lnL_curval (NcmFitState *fstate, gdouble m2lnL)
{
  fstate->m2lnL_curval = m2lnL;
}

/**
 * ncm_fit_state_get_m2lnL_curval:
 * @fstate: a #NcmFitState
 *
 * Gets the current value of the m2lnL of @fstate.
 *
 * Returns: Current value of the m2lnL
 */
gdouble
ncm_fit_state_get_m2lnL_curval (NcmFitState *fstate)
{
  return fstate->m2lnL_curval;
}

/**
 * ncm_fit_state_set_params_prec:
 * @fstate: a #NcmFitState
 * @prec: Precision of the parameters
 *
 * Sets the precision of the parameters of @fstate.
 *
 * This method is used by #NcmFit implementations to set the precision of the parameters.
 * It should not be used by the user.
 *
 */
void
ncm_fit_state_set_params_prec (NcmFitState *fstate, gdouble prec)
{
  fstate->params_prec = prec;
}

/**
 * ncm_fit_state_get_params_prec:
 * @fstate: a #NcmFitState
 *
 * Gets the precision of the parameters of @fstate.
 *
 * Returns: Precision of the parameters
 */
gdouble
ncm_fit_state_get_params_prec (NcmFitState *fstate)
{
  return fstate->params_prec;
}

/**
 * ncm_fit_state_get_data_len:
 * @fstate: a #NcmFitState
 *
 * Gets the number of data points used to compute @fstate.
 *
 * Returns: Number of data points
 */
guint
ncm_fit_state_get_data_len (NcmFitState *fstate)
{
  return fstate->data_len;
}

/**
 * ncm_fit_state_peek_fparams:
 * @fstate: a #NcmFitState
 *
 * Gets the free parameters #NcmVector of @fstate.
 *
 * Returns: (transfer none): Free parameters #NcmVector
 */
NcmVector *
ncm_fit_state_peek_fparams (NcmFitState *fstate)
{
  return fstate->fparams;
}

/**
 * ncm_fit_state_peek_covar:
 * @fstate: a #NcmFitState
 *
 * Gets the covariance #NcmMatrix of @fstate.
 *
 * Returns: (transfer none): Covariance #NcmMatrix
 */
NcmMatrix *
ncm_fit_state_peek_covar (NcmFitState *fstate)
{
  return fstate->covar;
}

