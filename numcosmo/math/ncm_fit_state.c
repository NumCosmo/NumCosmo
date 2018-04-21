/***************************************************************************
 *            ncm_fit_state.c
 *
 *  Thu November 29 15:27:11 2012
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
 * SECTION:ncm_fit_state
 * @title: NcmFitState
 * @short_description: State of a NcmFit object.
 *
 * FIXME
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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
 * FIXME
 * 
 * Returns: FIXME
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
 * @fstate: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmFitState *
ncm_fit_state_ref (NcmFitState *fstate)
{
  return g_object_ref (fstate);
}

/**
 * ncm_fit_state_free:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_fit_state_free (NcmFitState *fstate)
{
  g_object_unref (fstate);
}

/**
 * ncm_fit_state_clear:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_fit_state_clear (NcmFitState **fstate)
{
  g_clear_object (fstate);
}

/**
 * ncm_fit_state_set_all:
 * @fstate: FIXME
 * @data_len: FIXME
 * @fparam_len: FIXME
 * @dof: FIXME
 * @is_least_squares: FIXME
 * 
 * FIXME
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
 * @fstate: FIXME
 * 
 * FIXME
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
 * @fstate: FIXME
 * 
 * FIXME
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
    
    if ((fparam_diff_len || data_diff_len || fstate->ls_J == NULL) && fstate->fparam_len > 0)
    {
      ncm_matrix_clear (&fstate->ls_J);
      fstate->ls_J = ncm_matrix_new (fstate->data_len, fstate->fparam_len);
    }
    
    if (data_diff_len || fstate->ls_f == NULL)
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
 * @fstate: FIXME
 * @f: FIXME
 * @J: FIXME
 * 
 * FIXME
 * 
 */
/**
 * ncm_fit_state_set_niter:
 * @fstate: FIXME
 * @niter: FIXME
 * 
 * FIXME
 * 
 */
/**
 * ncm_fit_state_get_niter:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fit_state_set_m2lnL_prec:
 * @fstate: FIXME
 * @prec: FIXME
 * 
 * FIXME
 * 
 */
/**
 * ncm_fit_state_get_m2lnL_prec:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fit_state_set_m2lnL_curval:
 * @fstate: FIXME
 * @m2lnL: FIXME
 * 
 * FIXME
 * 
 */
/**
 * ncm_fit_state_get_m2lnL_curval:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fit_state_set_params_prec:
 * @fstate: FIXME
 * @prec: FIXME
 * 
 * FIXME
 * 
 */
/**
 * ncm_fit_state_get_params_prec:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fit_state_get_data_len:
 * @fstate: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */

