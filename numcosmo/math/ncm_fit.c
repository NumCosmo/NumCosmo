/***************************************************************************
 *            ncm_fit.c
 *
 *  Sat Aug 16 16:22:13 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit
 * @title: Model Fitting Abstract Class
 * @short_description: Class for implementing fitting methods
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
#include "math/util.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_fit_gsl_ls.h"
#include "math/ncm_fit_gsl_mm.h"
#include "math/ncm_fit_gsl_mms.h"
#include "math/ncm_fit_levmar.h"
#include "math/ncm_fit_nlopt.h"
#include "ncm_enum_types.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>

enum
{
  PROP_0,
  PROP_MAXITER,
  PROP_LIKELIHOOD,
  PROP_MSET,
  PROP_STATE,
  PROP_GRAD_TYPE,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmFit, ncm_fit, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcmFitConstraint, ncm_fit_constraint, (GBoxedCopyFunc)&ncm_fit_constraint_dup, (GBoxedFreeFunc)&ncm_fit_constraint_free);

/**
 * ncm_fit_constraint_new:
 * @func: FIXME
 * @tot: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmFitConstraint *
ncm_fit_constraint_new (NcmFit *fit, NcmMSetFunc *func, gdouble tot)
{
  NcmFitConstraint *fitc = g_new (NcmFitConstraint, 1);
  g_assert (ncm_mset_func_is_scalar (func) && ncm_mset_func_is_const (func));
  fitc->fit = ncm_fit_ref (fit);
  fitc->func = ncm_mset_func_ref (func);
  fitc->tot = tot;
  return fitc;
}

/**
 * ncm_fit_constraint_dup:
 * @fitc: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmFitConstraint *
ncm_fit_constraint_dup (NcmFitConstraint *fitc)
{
  return ncm_fit_constraint_new (fitc->fit, fitc->func, fitc->tot);
}

/**
 * ncm_fit_constraint_free:
 * @fitc: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fit_constraint_free (NcmFitConstraint *fitc)
{
  ncm_mset_func_free (fitc->func);
  g_free (fitc);
}

static void
ncm_fit_init (NcmFit *fit)
{
  fit->maxiter       = NCM_FIT_MAXITER;
  fit->m2lnL_reltol  = NCM_FIT_DEFAULT_M2LNL_RELTOL * 0.0 + 1e-8;
  fit->params_reltol = NC_HICOSMO_DEFAULT_PARAMS_RELTOL * 0.0 + 1e-8;
  fit->timer = g_timer_new ();
  fit->mtype = NCM_FIT_RUN_MSGS_NONE;

  fit->equality_constraints = g_ptr_array_sized_new (10);
  g_ptr_array_set_free_func (fit->equality_constraints, (GDestroyNotify) &ncm_fit_constraint_free);

  fit->inequality_constraints = g_ptr_array_sized_new (10);
  g_ptr_array_set_free_func (fit->inequality_constraints, (GDestroyNotify) &ncm_fit_constraint_free);
}

static void
_ncm_fit_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_parent_class)->constructed (object);
  {
    NcmFit *fit = NCM_FIT (object);
    gint n = ncm_dataset_get_n (fit->lh->dset);
    gint n_priors = ncm_likelihood_priors_length (fit->lh);
    gint data_dof = ncm_dataset_get_dof (fit->lh->dset);

    g_assert (ncm_dataset_all_init (fit->lh->dset));

    ncm_mset_prepare_fparam_map (fit->mset);

    /* 
     * It is no longer an error to fit 0 parameters, it just sets the value
     * of m2lnL in the fit object.
     * 
     * if (ncm_mset_fparam_len (fit->mset) == 0)
     * g_warning ("ncm_fit_new: mset object has 0 free parameters");
     * 
     */
    {
      guint data_len   = n + n_priors;
      guint fparam_len = ncm_mset_fparam_len (fit->mset);
      gint dof         = data_dof + n_priors - fparam_len;

      fit->fstate = ncm_fit_state_new (data_len, fparam_len, dof,
                                       NCM_FIT_GET_CLASS (fit)->is_least_squares);

      g_assert (data_len > 0);
    }
  }
}

static void
_ncm_fit_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFit *fit = NCM_FIT (object);
  g_return_if_fail (NCM_IS_FIT (object));

  switch (prop_id)
  {
    case PROP_MAXITER:
      ncm_fit_set_maxiter (fit, g_value_get_uint (value));
      break;
    case PROP_LIKELIHOOD:
      ncm_likelihood_clear (&fit->lh);
      fit->lh = g_value_dup_object (value);
      break;
    case PROP_MSET:
      ncm_mset_clear (&fit->mset);
      fit->mset = g_value_dup_object (value);
      break;
    case PROP_GRAD_TYPE:
      ncm_fit_set_grad_type (fit, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFit *fit = NCM_FIT (object);
  g_return_if_fail (NCM_IS_FIT (object));

  switch (prop_id)
  {
    case PROP_MAXITER:
      g_value_set_uint (value, ncm_fit_get_maxiter (fit));
      break;
    case PROP_LIKELIHOOD:
      g_value_set_object (value, fit->lh);
      break;
    case PROP_MSET:
      g_value_set_object (value, fit->mset);
      break;
    case PROP_GRAD_TYPE:
      g_value_set_enum (value, fit->grad.gtype);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_dispose (GObject *object)
{
  NcmFit *fit = NCM_FIT (object);

  ncm_likelihood_clear (&fit->lh);
  ncm_mset_clear (&fit->mset);
  ncm_fit_state_clear (&fit->fstate);

  if (fit->equality_constraints != NULL)
  {
    g_ptr_array_unref (fit->equality_constraints);
    fit->equality_constraints = NULL;
  }
  
  if (fit->inequality_constraints != NULL)
  {
    g_ptr_array_unref (fit->inequality_constraints);
    fit->inequality_constraints = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_parent_class)->dispose (object);
}

static void
ncm_fit_finalize (GObject *object)
{
  NcmFit *fit = NCM_FIT (object);

  g_timer_destroy (fit->timer);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_parent_class)->finalize (object);
}

static void _ncm_fit_reset (NcmFit *fit);

static void
ncm_fit_class_init (NcmFitClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_constructed;
  object_class->set_property = &_ncm_fit_set_property;
  object_class->get_property = &_ncm_fit_get_property;
  object_class->dispose      = &ncm_fit_dispose;
  object_class->finalize     = &ncm_fit_finalize;

  klass->is_least_squares = FALSE;
  klass->reset            = &_ncm_fit_reset;

  g_object_class_install_property (object_class,
                                   PROP_MAXITER,
                                   g_param_spec_uint ("maxiter",
                                                      NULL,
                                                      "Maximum number of interations",
                                                      0, G_MAXUINT32, NCM_FIT_MAXITER,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LIKELIHOOD,
                                   g_param_spec_object ("likelihood",
                                                        NULL,
                                                        "Likelihood object",
                                                        NCM_TYPE_LIKELIHOOD,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "Model set object",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_GRAD_TYPE,
                                   g_param_spec_enum ("grad-type",
                                                      NULL,
                                                      "Differentiation method",
                                                      NCM_TYPE_FIT_GRAD_TYPE, NCM_FIT_GRAD_NUMDIFF_FORWARD,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  

}

static void
_ncm_fit_reset (NcmFit *fit)
{
  ncm_mset_prepare_fparam_map (fit->mset);
  {
    guint n          = ncm_dataset_get_n (fit->lh->dset);
    gint n_priors    = ncm_likelihood_priors_length (fit->lh);
    guint data_len   = n + n_priors;
    guint fparam_len = ncm_mset_fparam_len (fit->mset);
    gint data_dof    = ncm_dataset_get_dof (fit->lh->dset);
    gint dof         = data_dof - fparam_len;

    g_assert (ncm_dataset_all_init (fit->lh->dset));
    g_assert (data_len > 0);

    ncm_fit_state_set_all (fit->fstate, data_len, fparam_len, dof,
                           NCM_FIT_GET_CLASS (fit)->is_least_squares);
    ncm_fit_state_reset (fit->fstate);
  }
}

/**
 * ncm_fit_new:
 * @ftype: a #NcmFitType
 * @algo_name: name of the algorithm to be used
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFit *
ncm_fit_new (NcmFitType ftype, gchar *algo_name, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  switch (ftype)
  {
    case NCM_FIT_TYPE_GSL_LS:
      return ncm_fit_gsl_ls_new (lh, mset, gtype);
      break;
    case NCM_FIT_TYPE_GSL_MM:
      return ncm_fit_gsl_mm_new_by_name (lh, mset, gtype, algo_name);
      break;
    case NCM_FIT_TYPE_GSL_MMS:
      return ncm_fit_gsl_mms_new_by_name (lh, mset, gtype, algo_name);
      break;
#ifdef NUMCOSMO_HAVE_LEVMAR
    case NCM_FIT_TYPE_LEVMAR:
      return ncm_fit_levmar_new_by_name (lh, mset, gtype, algo_name);
      break;
#endif /* NUMCOSMO_HAVE_LEVMAR */
#ifdef NUMCOSMO_HAVE_NLOPT
    case NCM_FIT_TYPE_NLOPT:
      return ncm_fit_nlopt_new_by_name (lh, mset, gtype, algo_name);
      break;
#endif /* NUMCOSMO_HAVE_NLOPT */
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_fit_ref:
 * @fit: a #NcmFit.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmFit *
ncm_fit_ref (NcmFit *fit)
{
  return g_object_ref (fit);
}

/**
 * ncm_fit_copy_new:
 * @fit: a #NcmFit
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 *
 * Duplicates the #NcmFit object with new references for its contents.
 *
 * Returns: (transfer full): FIXME
 */
NcmFit *
ncm_fit_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return NCM_FIT_GET_CLASS (fit)->copy_new (fit, lh, mset, gtype);
}

/**
 * ncm_fit_free:
 * @fit: a #NcmFit.
 *
 * FIXME
 *
 */
void
ncm_fit_free (NcmFit *fit)
{
  g_object_unref (fit);
}

/**
 * ncm_fit_clear:
 * @fit: a #NcmFit.
 *
 * FIXME
 *
 */
void
ncm_fit_clear (NcmFit **fit)
{
  g_clear_object (fit);
}

static NcmFitGrad _ncm_fit_grad_analitical = {
  NCM_FIT_GRAD_ANALYTICAL,
  "Analytical gradient",
  &ncm_fit_ls_J_an,
  &ncm_fit_ls_f_J_an,
  &ncm_fit_m2lnL_grad_an,
  &ncm_fit_m2lnL_val_grad_an,
};

static NcmFitGrad _ncm_fit_grad_numdiff_forward = {
  NCM_FIT_GRAD_NUMDIFF_FORWARD,
  "Numerical differentiantion (forward)",
  &ncm_fit_ls_J_nd_fo,
  &ncm_fit_ls_f_J_nd_fo,
  &ncm_fit_m2lnL_grad_nd_fo,
  &ncm_fit_m2lnL_val_grad_nd_fo,
};

static NcmFitGrad _ncm_fit_grad_numdiff_central = {
  NCM_FIT_GRAD_NUMDIFF_CENTRAL,
  "Numerical differentiantion (central)",
  &ncm_fit_ls_J_nd_ce,
  &ncm_fit_ls_f_J_nd_ce,
  &ncm_fit_m2lnL_grad_nd_ce,
  &ncm_fit_m2lnL_val_grad_nd_ce,
};

static NcmFitGrad _ncm_fit_grad_numdiff_accurate = {
  NCM_FIT_GRAD_NUMDIFF_ACCURATE,
  "Numerical differentiantion (Richardson extrapolation)",
  &ncm_fit_ls_J_nd_ce,
  &ncm_fit_ls_f_J_nd_ce,
  &ncm_fit_m2lnL_grad_nd_ac,
  &ncm_fit_m2lnL_val_grad_nd_ac,
};

/**
 * ncm_fit_set_grad_type:
 * @fit: a #NcmLikelihood.
 * @gtype: a #NcmFitGradType.
 *
 * FIXME
 *
 */
void
ncm_fit_set_grad_type (NcmFit *fit, NcmFitGradType gtype)
{
  fit->grad.gtype = gtype;
  
  if ((gtype == NCM_FIT_GRAD_ANALYTICAL) && !ncm_likelihood_has_m2lnL_grad (fit->lh))
    g_error ("Likelihood do not support analytical gradient, try to use a numerical algorithm.");

  switch (gtype)
  {
    case NCM_FIT_GRAD_ANALYTICAL:
      fit->grad = _ncm_fit_grad_analitical;
      break;
    case NCM_FIT_GRAD_NUMDIFF_FORWARD:
      fit->grad = _ncm_fit_grad_numdiff_forward;
      break;
    case NCM_FIT_GRAD_NUMDIFF_CENTRAL:
      fit->grad = _ncm_fit_grad_numdiff_central;
      break;
    case NCM_FIT_GRAD_NUMDIFF_ACCURATE:
      fit->grad = _ncm_fit_grad_numdiff_accurate;
      break;
    default:
      g_error ("Invalid gtype %d", gtype);
      break;
  }
}

/**
 * ncm_fit_ls_set_state:
 * @fit: a #NcmFit.
 * @prec: FIXME
 * @x: a #NcmVector.
 * @f: a #NcmVector.
 * @J: a #NcmMatrix.
 *
 * FIXME
 */
void
ncm_fit_ls_set_state (NcmFit *fit, gdouble prec, NcmVector *x, NcmVector *f, NcmMatrix *J)
{
  ncm_mset_fparams_set_vector (fit->mset, x);
  g_assert (fit->fstate->is_least_squares);

  fit->fstate->m2lnL_prec = prec;
  fit->fstate->m2lnL = gsl_blas_dnrm2 (ncm_vector_gsl (f));

  ncm_vector_memcpy (fit->fstate->ls_f, f);
  
  gsl_blas_dgemv (CblasTrans, 2.0, NCM_MATRIX_GSL(fit->fstate->ls_J), 
                  ncm_vector_gsl (fit->fstate->ls_f), 0.0, 
                  ncm_vector_gsl (fit->fstate->dm2lnL));

  ncm_matrix_memcpy (fit->fstate->ls_J, J);
}

/**
 * ncm_fit_set_maxiter:
 * @fit: a #NcmFit.
 * @maxiter: FIXME.
 *
 * FIXME
 */
void 
ncm_fit_set_maxiter (NcmFit *fit, guint maxiter)
{
  fit->maxiter = maxiter;
}

/**
 * ncm_fit_get_maxiter:
 * @fit: a #NcmFit.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
guint 
ncm_fit_get_maxiter (NcmFit *fit)
{
  return fit->maxiter;
}

/**
 * ncm_fit_add_equality_constraint:
 * @fit: a #NcmFit.
 * @func: FIXME
 * @tot: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_fit_add_equality_constraint (NcmFit *fit, NcmMSetFunc *func, gdouble tot)
{
  NcmFitConstraint *fitc = ncm_fit_constraint_new (fit, func, tot);
  g_ptr_array_add (fit->equality_constraints, fitc);
}

/**
 * ncm_fit_add_inequality_constraint:
 * @fit: a #NcmFit.
 * @func: FIXME
 * @tot: FIXME
 *
 * FIXME
 * 
 */
void
ncm_fit_add_inequality_constraint (NcmFit *fit, NcmMSetFunc *func, gdouble tot)
{
  NcmFitConstraint *fitc = ncm_fit_constraint_new (fit, func, tot);
  g_ptr_array_add (fit->inequality_constraints, fitc);
}

/**
 * ncm_fit_remove_equality_constraints:
 * @fit: a #NcmFit.
 *
 * FIXME
 * 
 */
void 
ncm_fit_remove_equality_constraints (NcmFit *fit)
{
  g_ptr_array_set_size (fit->equality_constraints, 0);
}

/**
 * ncm_fit_remove_inequality_constraints:
 * @fit: a #NcmFit.
 *
 * FIXME
 * 
 */
void 
ncm_fit_remove_inequality_constraints (NcmFit *fit)
{
  g_ptr_array_set_size (fit->inequality_constraints, 0);
}

/**
 * ncm_fit_has_equality_constraints:
 * @fit: a #NcmFit.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
guint 
ncm_fit_has_equality_constraints (NcmFit *fit)
{
  return fit->equality_constraints->len;
}

/**
 * ncm_fit_has_inequality_constraints:
 * @fit: a #NcmFit.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
guint
ncm_fit_has_inequality_constraints (NcmFit *fit)
{
  return fit->inequality_constraints->len;
}

/**
 * ncm_fit_ls_covar:
 * @fit: a #NcmFit.
 *
 * FIXME
 * 
 */
void
ncm_fit_ls_covar (NcmFit *fit)
{
  g_assert (fit->fstate->is_least_squares);
  ncm_fit_ls_J (fit, fit->fstate->ls_J);
  gsl_multifit_covar (NCM_MATRIX_GSL (fit->fstate->ls_J), 0.0, 
                      NCM_MATRIX_GSL (fit->fstate->covar));
  fit->fstate->has_covar = TRUE;
}

/**
 * ncm_fit_covar_fparam_var:
 * @fit: a #NcmFit.
 * @fpi: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_var (NcmFit *fit, guint fpi)
{
  g_assert (fit->fstate->has_covar);
  return ncm_matrix_get (fit->fstate->covar, fpi, fpi);
}

/**
 * ncm_fit_covar_fparam_sd:
 * @fit: a #NcmFit.
 * @fpi: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_sd (NcmFit *fit, guint fpi)
{
  return sqrt (ncm_fit_covar_fparam_var (fit, fpi));
}

/**
 * ncm_fit_covar_fparam_cov:
 * @fit: a #NcmFit.
 * @fpi1: FIXME
 * @fpi2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_cov (NcmFit *fit, guint fpi1, guint fpi2)
{
  g_assert (fit->fstate->has_covar);
  return ncm_matrix_get (fit->fstate->covar, fpi1, fpi2);
}

/**
 * ncm_fit_covar_fparam_cor:
 * @fit: a #NcmFit.
 * @fpi1: FIXME
 * @fpi2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_cor (NcmFit *fit, guint fpi1, guint fpi2)
{
  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2) / (ncm_fit_covar_fparam_sd (fit, fpi1) * ncm_fit_covar_fparam_sd (fit, fpi2));
}

/**
 * ncm_fit_covar_var:
 * @fit: a #NcmFit.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_var (NcmFit *fit, NcmModelID mid, guint pid)
{
  gint fpi = ncm_mset_fparam_get_fpi (fit->mset, mid, pid);
  if (fpi < 0)
    g_error ("Parameter (%d:%u) was not fit.", mid, pid);
  return ncm_fit_covar_fparam_var (fit, fpi);
}

/**
 * ncm_fit_covar_sd:
 * @fit: a #NcmFit.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_sd (NcmFit *fit, NcmModelID mid, guint pid)
{
  return sqrt (ncm_fit_covar_var (fit, mid, pid));
}

/**
 * ncm_fit_covar_cov:
 * @fit: a #NcmFit
 * @mid1: a #NcmModelID.
 * @pid1: FIXME
 * @mid2: a #NcmModelID.
 * @pid2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_cov (NcmFit *fit, NcmModelID mid1, guint pid1, NcmModelID mid2, guint pid2)
{
  gint fpi1 = ncm_mset_fparam_get_fpi (fit->mset, mid1, pid1);
  gint fpi2 = ncm_mset_fparam_get_fpi (fit->mset, mid2, pid2);

  if (fpi1 < 0 || fpi1 < 0)
    g_error ("Parameters (%d:%u, %d:%u) were not fit.", mid1, pid1, mid2, pid2);

  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2);
}

/**
 * ncm_fit_covar_cor:
 * @fit: a #NcmFit
 * @mid1: a #NcmModelID.
 * @pid1: FIXME
 * @mid2: a #NcmModelID.
 * @pid2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_cor (NcmFit *fit, NcmModelID mid1, guint pid1, NcmModelID mid2, guint pid2)
{
  gint fpi1 = ncm_mset_fparam_get_fpi (fit->mset, mid1, pid1);
  gint fpi2 = ncm_mset_fparam_get_fpi (fit->mset, mid2, pid2);

  if (fpi1 < 0 || fpi1 < 0)
    g_error ("Parameters (%d:%u, %d:%u) were not fit.", mid1, pid1, mid2, pid2);

  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2) / (ncm_fit_covar_fparam_sd (fit, fpi1) * ncm_fit_covar_fparam_sd (fit, fpi2));
}

gboolean
_ncm_fit_run_empty (NcmFit *fit, NcmFitRunMsgs mtype)
{
  fit->mtype = mtype;

  ncm_fit_m2lnL_val (fit, &fit->fstate->m2lnL);
  ncm_fit_log_step (fit, fit->fstate->m2lnL);

  fit->fstate->m2lnL_prec   = 0.0;
  fit->fstate->has_covar    = FALSE;
  fit->fstate->is_best_fit  = TRUE;
  fit->fstate->elapsed_time = 0.0;
  fit->fstate->niter        = 0;
  fit->fstate->func_eval    = 0;
  fit->fstate->grad_eval    = 0;
  
  return TRUE;
}

/**
 * ncm_fit_reset:
 * @fit: a #NcmFit
 *
 * FIXME
 *
 */
void 
ncm_fit_reset (NcmFit *fit)
{
  NCM_FIT_GET_CLASS (fit)->reset (fit);
}

/**
 * ncm_fit_run:
 * @fit: a #NcmFit
 * @mtype: a #NcmFitRunMsgs
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  gboolean run;

  fit->mtype = mtype;

  ncm_fit_reset (fit);
  ncm_fit_log_start (fit);
  g_timer_start (fit->timer);

  if (ncm_mset_fparam_len (fit->mset) == 0)
    run = _ncm_fit_run_empty (fit, mtype);
  else
    run = NCM_FIT_GET_CLASS (fit)->run (fit, mtype);

  fit->fstate->elapsed_time = g_timer_elapsed (fit->timer, NULL);
  fit->fstate->is_best_fit = run;
  
  ncm_fit_log_end (fit);

  return run;
}

/**
 * ncm_fit_is_least_squares:
 * @fit: a #NcmFit
 *
 * FIXME
 * 
 * Returns: whenever the fit object use a least squares method.
 */
gboolean 
ncm_fit_is_least_squares (NcmFit *fit)
{
  return NCM_FIT_GET_CLASS (fit)->is_least_squares;
}

/**
 * ncm_fit_get_desc:
 * @fit: a #NcmFit
 *
 * FIXME
 * 
 * Returns: (transfer none): fit object description.
 */
const gchar *
ncm_fit_get_desc (NcmFit *fit)
{
  return NCM_FIT_GET_CLASS (fit)->get_desc (fit);
}

/**
 * ncm_fit_log_start:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_start (NcmFit *fit)
{
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# Model fitting. Interating using:\n");
    g_message ("#  - solver:            %s\n", ncm_fit_get_desc (fit));
    g_message ("#  - differentiation:   %s\n", fit->grad.diff_name);

    if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("#");
  }
}

/**
 * ncm_fit_log_step_error:
 * @fit: a #NcmFit
 * @strerror: FIXME
 * @...: FIXME
 *
 * FIXME
 */
void
ncm_fit_log_step_error (NcmFit *fit, const gchar *strerror, ...)
{
  va_list ap;

  va_start (ap, strerror);
  if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
  {
    gchar *errmsg = g_strdup_vprintf (strerror, ap);
    g_message ("\n#  [%s] error = %s\n#", ncm_fit_get_desc (fit), errmsg);
    g_free (errmsg);
  }
  else if (fit->mtype == NCM_FIT_RUN_MSGS_FULL)
  {
    gchar *errmsg = g_strdup_vprintf (strerror, ap);
    g_message ("#  [%s] error = %s\n", ncm_fit_get_desc (fit), errmsg);
    g_free (errmsg);
  }
  va_end (ap);
  return;
}

/**
 * ncm_fit_log_end:
 * @fit: a #NcmFit
 *
 * This function prints in the log the precision with which the best-fit was found.
 */
void
ncm_fit_log_end (NcmFit *fit)
{
  g_timer_stop (fit->timer);
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("\n");
    g_message ("#  Minimum found with precision: % 8.5e (|df|/f or |dx|)\n", fit->fstate->m2lnL_prec);
  }
  ncm_fit_log_state (fit, fit->fstate->m2lnL);
  return;
}

/**
 * ncm_fit_log_state:
 * @fit: a #NcmFit
 * @m2lnL: minus two times the logarithm base e of the likelihood.
 *
 * This function prints in the log the current state.
 * 
 */
void
ncm_fit_log_state (NcmFit *fit, gdouble m2lnL)
{
  gint i;
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    gdouble elap_sec = g_timer_elapsed (fit->timer, NULL);
    gulong elap_min = elap_sec / 60;
    gulong elap_hour = elap_min / 60;
    gulong elap_day = elap_hour / 24;

    elap_sec = fmod (elap_sec, 60);
    elap_min = elap_min % 60;
    elap_hour =  elap_hour % 24;
    g_message ("#  Elapsed time: %02lu days, %02lu:%02lu:%010.7f\n", elap_day, elap_hour, elap_min, elap_sec);
    g_message ("#  iteration            [%06d]\n", fit->fstate->niter);
    g_message ("#  function evaluations [%06d]\n", fit->fstate->func_eval);
    g_message ("#  gradient evaluations [%06d]\n", fit->fstate->grad_eval);
    g_message ("#  degrees of freedom   [%06d]\n", fit->fstate->dof);
    g_message ("#  m2lnL     = %20.15g\n", m2lnL);
    g_message ("#  Fit parameters:\n#    ");
    for (i = 0; i < ncm_mset_fparam_len (fit->mset); i++)
      g_message ("[% -20.15g] ", ncm_mset_fparam_get (fit->mset, i));
    g_message ("\n");
  }
  return;
}

/**
 * ncm_fit_log_step:
 * @fit: a #NcmFit
 * @m2lnL: minus two times the logarithm base e of the likelihood.
 *
 * FIXME
 */
void
ncm_fit_log_step (NcmFit *fit, gdouble m2lnL)
{
  if (fit->mtype == NCM_FIT_RUN_MSGS_FULL)
    ncm_fit_log_state (fit, m2lnL);
  else if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
  {
    if (fit->fstate->niter < 10 || ((fit->fstate->niter % 10) == 0))
      ncm_message (".");
  }
  return;
}

/**
 * ncm_fit_log_finish:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_finish (NcmFit *fit)
{
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    const gdouble m2lnL = fit->fstate->m2lnL;
    const gint dof = fit->fstate->dof;
    const gdouble m2lnL_dof = m2lnL / dof;
    
    g_message ("#  m2lnL/dof = %20.15g\n", m2lnL_dof);
    g_message ("#  |m2lnL-dof|/sqrt(2*dof) = %20.15g,\n", fabs(m2lnL_dof) / sqrt(2.0 * dof));
    g_message ("#  GoF_tt = %4.2f%% = (%4.2f + %4.2f)%%; GoF = %4.2f%%\n",
               gsl_cdf_chisq_Q (dof + fabs(dof - m2lnL), dof) * 100.0 +
               gsl_cdf_chisq_P (dof - fabs(dof - m2lnL), dof) * 100.0,
               gsl_cdf_chisq_P (dof - fabs(dof - m2lnL), dof) * 100.0,
               gsl_cdf_chisq_Q (dof + fabs(dof - m2lnL), dof) * 100.0,
               gsl_cdf_chisq_Q (m2lnL, dof) * 100.0);
  }
  return;
}

/**
 * ncm_fit_log_info:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_info (NcmFit *fit)
{
  ncm_dataset_log_info (fit->lh->dset);
  ncm_mset_pretty_log (fit->mset);
  return;
}

/**
 * ncm_fit_log_covar:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_covar (NcmFit *fit)
{
  gint i, j;
  guint name_size = ncm_mset_max_fparam_name (fit->mset);
  guint free_params_len = ncm_mset_fparam_len (fit->mset);
  gchar *box = "---------------";
  g_assert (fit->fstate->has_covar);

  ncm_cfg_msg_sepa ();
  g_message ("# Fitted parameters covariance matrix\n");
  g_message ("#                                        ");
  for (i = 0; i < name_size; i++) g_message (" ");

  for (i = 0; i < free_params_len; i++)
    i ? g_message ("%s",box) : g_message ("-%s",box);
  if (i)
    g_message ("\n");

  for (i = 0; i < free_params_len; i++)
  {
    NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (fit->mset, i);
    const gchar *pname = ncm_mset_fparam_name (fit->mset, i);
    g_message ("# %*s[%02d%02d] = % -12.4g +/- % -12.4g |",
               name_size, pname, pi->mid, pi->pid,
               ncm_mset_fparam_get (fit->mset, i),
               ncm_fit_covar_fparam_sd (fit, i));
    for (j = 0; j < free_params_len; j++)
    {
      g_message (" % -12.4g |", ncm_fit_covar_fparam_cor (fit, i, j));
    }
    g_message ("\n");
  }
  g_message ("#                                        ");
  for (i = 0; i < name_size; i++) g_message (" ");
  for (i = 0; i < free_params_len; i++)
    i ? g_message ("%s",box) : g_message ("-%s",box);
  if (i)
    g_message ("\n");

  return;
}

/**
 * ncm_fit_fisher_matrix_print:
 * @fit: a #NcmFit
 * @out: name of the file
 * @header: pointer to the command line
 *
 * This function print the command line (first line, commented), the cosmological
 * parameters' names which were fitted (second line, commented) and the Fisher Matrix.
 *
 */
void
ncm_fit_fishermatrix_print (NcmFit *fit, FILE *out, gchar *header)
{
  gint i, j;
  guint name_size = ncm_mset_max_param_name (fit->mset);
  guint free_params_len = ncm_mset_fparam_len (fit->mset);

  if (header != NULL)
    fprintf (out, "# %s\n# ", header);
  else
    fprintf (out, "# ");

  for (i = 0; i < free_params_len; i++)
  {
    const gchar *pname = ncm_mset_fparam_name (fit->mset, i);
    fprintf (out, "%*s[%02d] ", name_size, pname, i);
  }
  fprintf (out, "\n");

  for (i = 0; i < free_params_len; i++)
  {
    for (j = 0; j < free_params_len; j++)
    {
      fprintf (out, " % -20.15g", ncm_fit_covar_fparam_cov (fit, i, j));
    }
    fprintf (out, "\n");
  }
}

/**
 * ncm_fit_data_m2lnL_val:
 * @fit: a #NcmFit
 * @data_m2lnL: (out): minus two times the logarithm base e of the likelihood.
 *
 * This function computes minus two times the logarithm base e of the likelihood 
 * using only the data set and not considering any prior. The result is set 
 * on @data_m2lnL.
 * 
 */
void
ncm_fit_data_m2lnL_val (NcmFit *fit, gdouble *data_m2lnL)
{
  ncm_dataset_m2lnL_val (fit->lh->dset, fit->mset, data_m2lnL);
}

/**
 * ncm_fit_priors_m2lnL_val:
 * @fit: a #NcmFit
 * @priors_m2lnL: (out): minus two times the logarithm base e of the likelihood.
 *
 * This function computes minus two times the logarithm base e of the likelihood 
 * using the data set and taking into account the assumed priors. The result is 
 * set on @priors_m2lnL.
 * 
 */
void
ncm_fit_priors_m2lnL_val (NcmFit *fit, gdouble *priors_m2lnL)
{
  ncm_likelihood_priors_m2lnL_val (fit->lh, fit->mset, priors_m2lnL);
}

/**
 * ncm_fit_m2lnL_val:
 * @fit: a #NcmFit
 * @m2lnL: (out): minus two times the logarithm base e of the likelihood. 
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL)
{
  ncm_likelihood_m2lnL_val (fit->lh, fit->mset, m2lnL);
  
  fit->fstate->func_eval++;
}

/**
 * ncm_fit_ls_f:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_ls_f (NcmFit *fit, NcmVector *f)
{
  ncm_likelihood_leastsquares_f (fit->lh, fit->mset, f);
  fit->fstate->func_eval++;
}

/**
 * ncm_fit_m2lnL_grad:
 * @fit: a #NcmFit
 * @df: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df)
{
  fit->grad.m2lnL_grad (fit, df);
}

/**
 * ncm_fit_m2lnL_grad_an:
 * @fit: a #NcmFit
 * @df: a #NcmVector
 *
 * Analytical gradient.
 * 
 */
void
ncm_fit_m2lnL_grad_an (NcmFit *fit, NcmVector *df)
{
  ncm_likelihood_m2lnL_grad (fit->lh, fit->mset, df);
  fit->fstate->grad_eval++;  
}


/**
 * ncm_fit_m2lnL_grad_nd_fo:
 * @fit: a #NcmFit
 * @grad: a #NcmVector
 *
 * Numerical differentiation (forward).
 * 
 */
void
ncm_fit_m2lnL_grad_nd_fo (NcmFit *fit, NcmVector *grad)
{
  gdouble m2lnL;
  ncm_fit_m2lnL_val_grad_nd_fo (fit, &m2lnL, grad);
}

/**
 * ncm_fit_m2lnL_grad_nd_ce:
 * @fit: a #NcmFit
 * @grad: a #NcmVector
 *
 * Numerical differentiation (central).
 * 
 */
void
ncm_fit_m2lnL_grad_nd_ce (NcmFit *fit, NcmVector *grad)
{
  gint i;
  guint fparam_len = ncm_mset_fparam_len (fit->mset);

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p = ncm_mset_fparam_get (fit->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (fit->mset, i));
    const gdouble h = p_scale * GSL_ROOT3_DBL_EPSILON;
    const gdouble pph = p + h;
    const gdouble pmh = p - h;
    const gdouble twoh = pph - pmh;
    const gdouble one_2h = 1.0 / twoh;
    gdouble m2lnL_pph, m2lnL_pmh;

    ncm_mset_fparam_set (fit->mset, i, pph);
    ncm_fit_m2lnL_val (fit, &m2lnL_pph);

    ncm_mset_fparam_set (fit->mset, i, pmh);
    ncm_fit_m2lnL_val (fit, &m2lnL_pmh);

    ncm_vector_set (grad, i, (m2lnL_pph - m2lnL_pmh) * one_2h);
    ncm_mset_fparam_set (fit->mset, i, p);
  }

  fit->fstate->grad_eval++;
}

typedef struct __ncm_fit_numdiff_1
{
  NcmFit *fit;
  guint n;
} __ncm_fit_numdiff_1;

static gdouble
_ncm_fit_numdiff_1_m2lnL (gdouble x, gpointer userdata)
{
  gdouble res;
  __ncm_fit_numdiff_1 *nd = (__ncm_fit_numdiff_1 *)userdata;
  ncm_mset_fparam_set (nd->fit->mset, nd->n, x);
  ncm_fit_m2lnL_val (nd->fit, &res);
  return res;
}

/**
 * ncm_fit_m2lnL_grad_nd_ac:
 * @fit: a #NcmFit
 * @grad: a #NcmVector
 *
 * Numerical differentiation (accurate).
 * 
 */
void
ncm_fit_m2lnL_grad_nd_ac (NcmFit *fit, NcmVector *grad)
{
  gsl_function F;
  __ncm_fit_numdiff_1 nd;
  guint fparam_len = ncm_mset_fparam_len (fit->mset);
  gint i;

  nd.fit = fit;
  F.params = &nd;
  F.function = &_ncm_fit_numdiff_1_m2lnL;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p = ncm_mset_fparam_get (fit->mset, i);
    const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
    gdouble err, diff;
    nd.n = i;
    diff = ncm_numdiff_1 (&F, p, p_scale, &err);
    ncm_vector_set (grad, i, diff);
    ncm_mset_fparam_set (fit->mset, i, p);
  }

  fit->fstate->grad_eval++;
}

/**
 * ncm_fit_m2lnL_val_grad:
 * @fit: a #NcmFit
 * @result: (out): FIXME
 * @df: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df)
{
  fit->grad.m2lnL_val_grad (fit, result, df);
}

/**
 * ncm_fit_m2lnL_val_grad_an:
 * @fit: a #NcmFit
 * @result: (out): FIXME
 * @df: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val_grad_an (NcmFit *fit, gdouble *result, NcmVector *df)
{
  ncm_likelihood_m2lnL_val_grad (fit->lh, fit->mset, result, df);

  fit->fstate->func_eval++;
  fit->fstate->grad_eval++;
}

/**
 * ncm_fit_m2lnL_val_grad_nd_fo:
 * @fit: a #NcmFit
 * @m2lnL: (out): minus two times the logarithm base e of the likelihood.
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val_grad_nd_fo (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  gint i;
  guint fparam_len = ncm_mset_fparam_len (fit->mset);

  ncm_fit_m2lnL_val (fit, m2lnL);

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p = ncm_mset_fparam_get (fit->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (fit->mset, i));
    const gdouble htilde = GSL_SQRT_DBL_EPSILON * p_scale;
    const gdouble pph = p + htilde;
    const gdouble h = pph - p;
    const gdouble one_h = 1.0 / h;
    gdouble m2lnL_pph;
    ncm_mset_fparam_set (fit->mset, i, pph);
    ncm_fit_m2lnL_val (fit, &m2lnL_pph);
    ncm_vector_set (grad, i, (m2lnL_pph - *m2lnL) * one_h);
    ncm_mset_fparam_set (fit->mset, i, p);
  }
  fit->fstate->func_eval++;
  fit->fstate->grad_eval++;
}

/**
 * ncm_fit_m2lnL_val_grad_nd_ce:
 * @fit: a #NcmFit
 * @m2lnL: minus two times the logarithm base e of the likelihood.
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val_grad_nd_ce (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  ncm_fit_m2lnL_val (fit, m2lnL);
  ncm_fit_m2lnL_grad_nd_ce (fit, grad);
}

/**
 * ncm_fit_m2lnL_val_grad_nd_ac:
 * @fit: a #NcmFit
 * @m2lnL: (out): minus two times the logarithm base e of the likelihood.
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val_grad_nd_ac (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  ncm_fit_m2lnL_val (fit, m2lnL);
  ncm_fit_m2lnL_grad_nd_ac (fit, grad);
}

/**
 * ncm_fit_ls_J:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_ls_J (NcmFit *fit, NcmMatrix *J)
{
  fit->grad.ls_J (fit, J);
}

/**
 * ncm_fit_ls_J_an:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 * 
 */
void
ncm_fit_ls_J_an (NcmFit *fit, NcmMatrix *J)
{
  ncm_likelihood_leastsquares_J (fit->lh, fit->mset, J);
  fit->fstate->grad_eval++;

  return;
}

/**
 * ncm_fit_ls_J_nd_fo:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_ls_J_nd_fo (NcmFit *fit, NcmMatrix *J)
{
  gint i;
  guint fparam_len = ncm_mset_fparam_len (fit->mset);

  ncm_fit_ls_f (fit, fit->fstate->ls_f);

  for (i = 0; i < fparam_len; i++)
  {
    NcmVector *J_col_i = ncm_matrix_get_col (J, i);
    const gdouble p = ncm_mset_fparam_get (fit->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (fit->mset, i));
    const gdouble htilde = p_scale * GSL_SQRT_DBL_EPSILON;
    const gdouble pph = p + htilde;
    const gdouble h = pph - p;
    const gdouble one_h = 1.0 / h;
    gint k;
    ncm_mset_fparam_set (fit->mset, i, pph);
    ncm_fit_ls_f (fit, J_col_i);
    for (k = 0; k < ncm_vector_len (fit->fstate->ls_f); k++)
      ncm_vector_set (J_col_i, k,
                      (ncm_vector_get (J_col_i, k) - ncm_vector_get (fit->fstate->ls_f, k)) * one_h
                      );
    ncm_mset_fparam_set (fit->mset, i, p);
    ncm_vector_free (J_col_i);
  }
  fit->fstate->grad_eval++;
}

/**
 * ncm_fit_ls_J_nd_ce:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 * 
 */
void
ncm_fit_ls_J_nd_ce (NcmFit *fit, NcmMatrix *J)
{
  gint i;
  guint fparam_len = ncm_mset_fparam_len (fit->mset);

  for (i = 0; i < fparam_len; i++)
  {
    NcmVector *J_col_i = ncm_matrix_get_col (J, i);
    const gdouble p = ncm_mset_fparam_get (fit->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (fit->mset, i));
    const gdouble h = p_scale * GSL_ROOT3_DBL_EPSILON;
    const gdouble pph = p + h;
    const gdouble pmh = p - h;
    const gdouble twoh = pph - pmh;
    const gdouble one_2h = 1.0 / twoh;
    gint k;

    ncm_mset_fparam_set (fit->mset, i, pph);
    ncm_fit_ls_f (fit, J_col_i);

    ncm_mset_fparam_set (fit->mset, i, pmh);
    ncm_fit_ls_f (fit, fit->fstate->ls_f);

    for (k = 0; k < ncm_vector_len (fit->fstate->ls_f); k++)
      ncm_vector_set (J_col_i, k,
                      (ncm_vector_get (J_col_i, k) - ncm_vector_get (fit->fstate->ls_f, k)) * one_2h
                      );
    ncm_mset_fparam_set (fit->mset, i, p);
    ncm_vector_free (J_col_i);
  }

  fit->fstate->grad_eval++;
}


/**
 * ncm_fit_ls_f_J:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_ls_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  fit->grad.ls_f_J (fit, f, J);
}

/**
 * ncm_fit_ls_f_J_an:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_ls_f_J_an (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_likelihood_leastsquares_f_J (fit->lh, fit->mset, f, J);

  fit->fstate->func_eval++;
  fit->fstate->grad_eval++;

  return;
}

/**
 * ncm_fit_ls_f_J_nd_fo:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_ls_f_J_nd_fo (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_fit_ls_J_nd_fo (fit, J);
  ncm_vector_memcpy (f, fit->fstate->ls_f);
}

/**
 * ncm_fit_ls_f_J_nd_ce:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_ls_f_J_nd_ce (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_fit_ls_f (fit, f);
  ncm_fit_ls_J_nd_ce (fit, J);
}

typedef struct __ncm_fit_numdiff_2
{
  NcmFit *fit;
  guint n1;
  guint n2;
  gdouble v;
  gdouble p1_scale;
  gdouble p2_scale;
} _ncm_fit_numdiff_2;

static gdouble
_ncm_fit_numdiff_2_m2lnL (gdouble u, gpointer userdata)
{
  gdouble res;
  _ncm_fit_numdiff_2 *nd = (_ncm_fit_numdiff_2 *)userdata;
  if (nd->n1 == nd->n2)
    ncm_mset_fparam_set (nd->fit->mset, nd->n1, u);
  else
  {
    ncm_mset_fparam_set (nd->fit->mset, nd->n1, (u + nd->v) / nd->p1_scale);
    ncm_mset_fparam_set (nd->fit->mset, nd->n2, (u - nd->v) / nd->p2_scale);
  }
  ncm_fit_m2lnL_val (nd->fit, &res);
  return res;
}

/**
 * ncm_fit_numdiff_m2lnL_hessian:
 * @fit: a #NcmFit
 * @H: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_numdiff_m2lnL_hessian (NcmFit *fit, NcmMatrix *H)
{
  gsl_function F;
  _ncm_fit_numdiff_2 nd;
  gint i, j;
  gdouble fx;
  const gdouble target_err = 1e-5;
  guint free_params_len = ncm_mset_fparams_len (fit->mset);

  ncm_fit_m2lnL_val (fit, &fx);

  nd.fit = fit;
  F.params = &nd;
  F.function = &_ncm_fit_numdiff_2_m2lnL;

  /* Diagonal */
  for (i = 0; i < free_params_len; i++)
  {
    const gdouble p = ncm_mset_fparam_get (fit->mset, i);
    const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
    gdouble err, diff;
    nd.n1 = i;
    nd.n2 = i;
    diff = ncm_numdiff_2_err (&F, &fx, p, p_scale, target_err, &err);
    if (fabs(err / diff) > target_err)
      g_warning ("ncm_fit_numdiff_m2lnL_hessian: effective error (% 20.15e) bigger then (% 20.15e)", fabs(err / diff), target_err);
    ncm_matrix_set (H, i, i, diff);
    ncm_mset_fparam_set (fit->mset, i, p);
  }

  for (i = 0; i < free_params_len; i++)
  {
    for (j = i + 1; j < free_params_len; j++)
    {
      const gdouble p1_scale = 1.0 / ncm_mset_fparam_get_scale (fit->mset, i);
      const gdouble p2_scale = 1.0 / ncm_mset_fparam_get_scale (fit->mset, j);
      const gdouble p1 = ncm_mset_fparam_get (fit->mset, i);
      const gdouble p2 = ncm_mset_fparam_get (fit->mset, j);
      const gdouble u = (p1_scale * p1 + p2_scale * p2) / 2.0;
      const gdouble v = (p1_scale * p1 - p2_scale * p2) / 2.0;
      const gdouble u_scale = 1.0;
      gdouble err, diff;
      nd.n1 = i;
      nd.n2 = j;
      nd.v = v;
      nd.p1_scale = p1_scale;
      nd.p2_scale = p2_scale;
      diff = ncm_numdiff_2_err (&F, &fx, u, u_scale, target_err, &err);
      if (fabs(err / diff) > target_err)
        g_warning ("ncm_fit_numdiff_m2lnL_hessian: effective error (% 20.15e) bigger then (% 20.15e)", fabs(err / diff), target_err);
      ncm_matrix_set (H, i, j,
                      0.5 * ( p1_scale * p2_scale * diff -
                             (p2_scale / p1_scale) * ncm_matrix_get (H, i, i) -
                             (p1_scale / p2_scale) * ncm_matrix_get (H, j, j)
                             )
                      );
      ncm_matrix_set (H, j, i, ncm_matrix_get (H, i, j));
      //printf ("d2[%d %d] = % 20.15g\n", i, j, (diff - gsl_matrix_get (H, k, k) - gsl_matrix_get (H, l, l)) / 2.0);
      ncm_mset_fparam_set (fit->mset, i, p1);
      ncm_mset_fparam_set (fit->mset, j, p2);
    }
  }
}

/**
 * ncm_fit_numdiff_m2lnL_covar:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_numdiff_m2lnL_covar (NcmFit *fit)
{
  gint ret;

  if (ncm_mset_fparam_len (fit->mset) == 0)
    g_error ("ncm_fit_numdiff_m2lnL_covar: mset object has 0 free parameters");

  ncm_fit_numdiff_m2lnL_hessian (fit, fit->fstate->hessian);
  ncm_matrix_memcpy (fit->fstate->covar, fit->fstate->hessian);
  ncm_matrix_scale (fit->fstate->covar, 0.5);

  ret = gsl_linalg_cholesky_decomp (NCM_MATRIX_GSL (fit->fstate->covar));
  if (ret == GSL_SUCCESS)
  {
    ret = gsl_linalg_cholesky_invert (NCM_MATRIX_GSL (fit->fstate->covar));
    NCM_TEST_GSL_RESULT ("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_cholesky_invert]", ret);
  }
  else if (ret == GSL_EDOM)
  {
    NcmMatrix *LU = ncm_matrix_dup (fit->fstate->hessian);
    gsl_permutation *p = gsl_permutation_alloc (ncm_mset_fparam_len (fit->mset));
    gint signum;
    gint ret1;

    ncm_matrix_scale (LU, 0.5);
    
    g_warning ("ncm_fit_numdiff_m2lnL_covar: covariance matrix not positive definite, errors are not trustworthy.");
    
    ret1 = gsl_linalg_LU_decomp (NCM_MATRIX_GSL (LU), p, &signum);
    NCM_TEST_GSL_RESULT ("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_LU_decomp]", ret1);
    
    ret1 = gsl_linalg_LU_invert (NCM_MATRIX_GSL (LU), p, NCM_MATRIX_GSL (fit->fstate->covar));
    NCM_TEST_GSL_RESULT ("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_LU_invert]", ret1);
    
    gsl_permutation_free (p);
    ncm_matrix_free (LU);
  }
  else
    NCM_TEST_GSL_RESULT ("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_cholesky_decomp]", ret);
  fit->fstate->has_covar = TRUE;
}

typedef struct _FitDProb
{
  NcmMSetPIndex pi;
  NcmFit *fit_val;
  NcmFit *fit;
} FitDProb;

static gdouble
fit_dprob(gdouble val, gpointer p)
{
  FitDProb *dprob_arg = (FitDProb *)p;
  ncm_mset_param_set (dprob_arg->fit_val->mset, dprob_arg->pi.mid, dprob_arg->pi.pid, val);
  ncm_fit_run (dprob_arg->fit_val, NCM_FIT_RUN_MSGS_NONE);
  if (dprob_arg->fit->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message (".");
  return exp (-(dprob_arg->fit_val->fstate->m2lnL - dprob_arg->fit->fstate->m2lnL) / 2.0);
}

/**
 * ncm_fit_prob:
 * @fit: a #NcmFit.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 * @a: FIXME
 * @b: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_prob (NcmFit *fit, NcmModelID mid, guint pid, gdouble a, gdouble b)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_dup (fit->mset);
  gsl_integration_workspace **w;
  FitDProb dprob_arg;
  gdouble result, error;
  gint error_code;
  gsl_function F;

  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, fit->lh, mset_val, fit->grad.gtype);

  dprob_arg.pi.mid = mid;
  dprob_arg.pi.pid  = pid;
  dprob_arg.fit     = fit;
  dprob_arg.fit_val = fit_val;

  F.function = &fit_dprob;
  F.params = &dprob_arg;

  w = ncm_integral_get_workspace();
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("#");
  error_code = gsl_integration_qags (&F, a, b, 1e-10, NCM_INTEGRAL_ERROR, NCM_INTEGRAL_PARTITION, *w, &result, &error);
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("\n");
  NCM_TEST_GSL_RESULT ("ncm_fit_prob", error_code);
  ncm_memory_pool_return (w);

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);
  return result;
}

/**
 * ncm_fit_dprob:
 * @fit: a #NcmFit.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 * @a: FIXME
 * @b: FIXME
 * @step: FIXME
 * @norm: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_fit_dprob (NcmFit *fit, NcmModelID mid, guint pid, gdouble a, gdouble b, gdouble step, gdouble norm)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_dup (fit->mset);
  FitDProb dprob_arg;
  gdouble point;

  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, fit->lh, mset_val, fit->grad.gtype);

  dprob_arg.pi.mid = mid;
  dprob_arg.pi.pid  = pid;
  dprob_arg.fit     = fit;
  dprob_arg.fit_val = fit_val;

  for (point = a; point <= b; point += step)
  {
    g_message ("%g %g\n", point, fit_dprob (point, &dprob_arg) / norm);
  }

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);

  return;
}

/**
 * ncm_fit_lr_test_range:
 * @fit: a #NcmFit
 * @mid: a #NcmModelID.
 * @pid: FIXME
 * @start: FIXME
 * @stop: FIXME
 * @step: FIXME
 *
 * FIXME
 */
void
ncm_fit_lr_test_range (NcmFit *fit, NcmModelID mid, guint pid, gdouble start, gdouble stop, gdouble step)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_dup (fit->mset);
  gdouble walk;

  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, fit->lh, mset_val, fit->grad.gtype);

  for (walk = start; walk <= stop; walk += step)
  {
    ncm_mset_param_set (fit_val->mset, mid, pid, walk);
    ncm_fit_run (fit_val, NCM_FIT_RUN_MSGS_NONE);
    //g_message ("%g %g %g\n", walk, (fit_val->fstate->m2lnL - fit->fstate->m2lnL), gsl_cdf_chisq_Q (fit_val->fstate->m2lnL - fit->fstate->m2lnL, 1));
    g_message ("%g %g %g %g %g\n", walk,
               (fit_val->fstate->m2lnL - fit->fstate->m2lnL),
               gsl_ran_chisq_pdf (fit_val->fstate->m2lnL - fit->fstate->m2lnL, 1),
               gsl_cdf_chisq_Q (fit_val->fstate->m2lnL - fit->fstate->m2lnL, 1),
               gsl_cdf_ugaussian_Q (sqrt(fit_val->fstate->m2lnL - fit->fstate->m2lnL))
               );
  }

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);
  return;
}

/**
 * ncm_fit_lr_test:
 * @fit: a #NcmFit.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 * @val: FIXME
 * @dof: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_lr_test (NcmFit *fit, NcmModelID mid, guint pid, gdouble val, gint dof)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_dup (fit->mset);
  gdouble result;

  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, fit->lh, mset_val, fit->grad.gtype);

  ncm_mset_param_set (fit_val->mset, mid, pid, val);
  ncm_fit_run (fit_val, NCM_FIT_RUN_MSGS_NONE);
  result = gsl_cdf_chisq_Q (fit_val->fstate->m2lnL - fit->fstate->m2lnL, dof);

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);
  return result;
}

/*
 * ncm_fit_chisq_test:
 * @fit: a #NcmFit
 * @bins: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 *
 gdouble
 ncm_fit_chisq_test (NcmFit *fit, size_t bins)
 {
   NcmLikelihood *lh = fit->lh;
   NcmDataset *ds = lh->ds;
   NcmData *data = nc_dataset_get_data (ds, 0);
   NcFunction distance = nc_function_get (data->res_type);
   gsl_histogram *obs = gsl_histogram_alloc (bins);
   gsl_histogram *exp = gsl_histogram_alloc (bins);
   gsl_vector_view y_exp;
   gdouble y_min, y_max;
   gint i;
   if (data == NULL)
   return GSL_NAN;

   ncm_distance_thread_pool (distance.f, fit->cp, data->x, fit->f, NULL, NULL);
   y_exp = gsl_vector_subvector (fit->f, 0, data->x->size);
   y_min = GSL_MIN (gsl_vector_min (data->y), gsl_vector_min (&y_exp.vector))*0.999;
   y_max = GSL_MAX (gsl_vector_max (data->y), gsl_vector_max (&y_exp.vector))*1.001;

   gsl_histogram_set_ranges_uniform (obs, y_min, y_max);
   gsl_histogram_set_ranges_uniform (exp, y_min, y_max);

   for (i = 0; i < data->x->size; i++)
   {
     gsl_histogram_increment (obs, gsl_vector_get (data->y, i));
     gsl_histogram_increment (exp, gsl_vector_get (&y_exp.vector, i));
     }

     gsl_histogram_sub (obs, exp);
     //  if (gsl_histogram_min_val (exp) == 0.0)
     //    g_error ("Cannot chisq test");
     gsl_histogram_mul (obs, obs);
     gsl_histogram_div (obs, exp);
     //gsl_histogram_fprintf (fit->log, obs, "%f", "%f");
     g_message ("\n");
     //gsl_histogram_fprintf (fit->log, exp, "%f", "%f");
     g_message ("\n");
     g_message ("BLA [%g, %g] %zu %g [%g]\n", y_min, y_max, bins, gsl_histogram_sum (obs), gsl_histogram_sum (exp));
     return gsl_histogram_sum (obs);
     }
     */

/**
 * ncm_fit_function_error:
 * @fit: a #NcmFit
 * @func: a #NcmMSetFunc
 * @x: FIXME
 * @pretty_print: FIXME
 * @f: (out): FIXME
 * @sigma_f: (out): FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_function_error (NcmFit *fit, NcmMSetFunc *func, gdouble *x, gboolean pretty_print, gdouble *f, gdouble *sigma_f)
{
  NcmVector *v = ncm_mset_func_numdiff_fparams (func, fit->mset, x, NULL);
  NcmVector *tmp1 = ncm_vector_dup (v);
  gdouble result;
  gint ret;

  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->fstate->covar), ncm_vector_gsl (v), 0.0, ncm_vector_gsl (tmp1));
  NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
  ret = gsl_blas_ddot (ncm_vector_gsl (v), ncm_vector_gsl (tmp1), &result);
  NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);

  *f = ncm_mset_func_eval (func, fit->mset, x);
  *sigma_f = sqrt(result);
  if (pretty_print)
    g_message ("# % -12.4g +/- % -12.4g\n", *f, *sigma_f);

  ncm_vector_free (v);
  ncm_vector_free (tmp1);

  return;
}

/**
 * ncm_fit_function_cov:
 * @fit: a #NcmFit
 * @func1: a #NcmMSetFunc
 * @z1: FIXME
 * @func2: a #NcmMSetFunc
 * @z2: FIXME
 * @pretty_print: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_function_cov (NcmFit *fit, NcmMSetFunc *func1, gdouble z1, NcmMSetFunc *func2, gdouble z2, gboolean pretty_print)
{
  g_assert_not_reached ();
  /*
   gdouble result, cor, s1, s2;
   gint ret;
   NcmVector *tmp1 = ncm_fit_params_get_tmp_vector (fit->pt, fit->mset);
   NcmVector *tmp2 = ncm_fit_params_get_tmp_vector (fit->pt, fit->mset);

   NCM_FUNC_DF (func1, fit->mset, fit->pt, z1, tmp1);
   NCM_FUNC_DF (func2, fit->mset, fit->pt, z2, tmp2);

   ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->fstate->covar), ncm_vector_gsl (tmp1), 0.0, ncm_vector_gsl (fit->df));
   NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
   ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp1), &result);
   NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
   s1 = sqrt(result);

   ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->fstate->covar), ncm_vector_gsl (tmp2), 0.0, ncm_vector_gsl (fit->df));
   NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
   ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp2), &result);
   NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
   s2 = sqrt(result);

   ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->fstate->covar), ncm_vector_gsl (tmp1), 0.0, ncm_vector_gsl (fit->df));
   NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
   ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp2), &result);
   NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
   cor = result / (s1*s2);

   if (pretty_print)
   {
     g_message ("# % -12.4g\t| % -12.4g\n", NCM_FUNC_F(func1, fit->mset, z1), NCM_FUNC_F(func2, fit->mset, z2));
     g_message ("#---------------------------------------------\n");
     g_message ("# % -12.4g\t | % -12.4g\t| % -12.4g\n", s1, 1.0, cor);
     g_message ("# % -12.4g\t | % -12.4g\t| % -12.4g\n", s2, cor, 1.0);
     g_message ("#---------------------------------------------\n");
}

ncm_fit_params_return_tmp_vector (fit->pt, tmp1);
ncm_fit_params_return_tmp_vector (fit->pt, tmp2);

return sqrt(result);
*/
}
