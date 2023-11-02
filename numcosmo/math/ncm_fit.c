/***************************************************************************
 *            ncm_fit.c
 *
 *  Sat Aug 16 16:22:13 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * @title: NcmFit
 * @short_description: Abstract class for implementing fitting methods.
 *
 * This object implements a abstract class for implementing fitting methods.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_fit_gsl_ls.h"
#include "math/ncm_fit_gsl_mm.h"
#include "math/ncm_fit_gsl_mms.h"
#include "math/ncm_fit_levmar.h"
#include "math/ncm_fit_nlopt.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_vector.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_LIKELIHOOD,
  PROP_MSET,
  PROP_STATE,
  PROP_GRAD_TYPE,
  PROP_MAXITER,
  PROP_M2LNL_RELTOL,
  PROP_M2LNL_ABSTOL,
  PROP_PARAMS_RELTOL,
  PROP_EQC,
  PROP_EQC_TOT,
  PROP_INEQC,
  PROP_INEQC_TOT,
  PROP_SUBFIT,
  PROP_SIZE,
};

typedef struct _NcmFitPrivate
{
  /*< private >*/
  GObject parent_instance;
  NcmLikelihood *lh;
  NcmMSet *mset;
  NcmFitState *fstate;
  NcmFitRunMsgs mtype;
  NcmFitGrad grad;
  guint maxiter;
  gdouble m2lnL_reltol;
  gdouble m2lnL_abstol;
  gdouble params_reltol;
  GTimer *timer;
  NcmObjArray *equality_constraints;
  GArray *equality_constraints_tot;
  NcmObjArray *inequality_constraints;
  GArray *inequality_constraints_tot;
  NcmFit *sub_fit;
  NcmDiff *diff;
} NcmFitPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmFit, ncm_fit, G_TYPE_OBJECT);

static void
ncm_fit_init (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->maxiter       = 0;
  self->m2lnL_reltol  = 0.0;
  self->m2lnL_abstol  = 0.0;
  self->params_reltol = 0.0;
  self->timer         = g_timer_new ();
  self->mtype         = NCM_FIT_RUN_MSGS_NONE;

  self->equality_constraints       = ncm_obj_array_new ();
  self->inequality_constraints     = ncm_obj_array_new ();
  self->equality_constraints_tot   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->inequality_constraints_tot = g_array_new (FALSE, FALSE, sizeof (gdouble));

  self->sub_fit = NULL;
  self->diff    = ncm_diff_new ();
}

static void
_ncm_fit_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFit *fit         = NCM_FIT (object);
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_return_if_fail (NCM_IS_FIT (object));

  switch (prop_id)
  {
    case PROP_LIKELIHOOD:
      ncm_likelihood_clear (&self->lh);
      self->lh = g_value_dup_object (value);
      g_assert (self->lh != NULL);
      break;
    case PROP_MSET:
      ncm_mset_clear (&self->mset);
      self->mset = g_value_dup_object (value);
      g_assert (self->mset != NULL);
      break;
    case PROP_STATE:
      ncm_fit_state_clear (&self->fstate);
      self->fstate = g_value_dup_object (value);
      break;
    case PROP_GRAD_TYPE:
      ncm_fit_set_grad_type (fit, g_value_get_enum (value));
      break;
    case PROP_MAXITER:
      ncm_fit_set_maxiter (fit, g_value_get_uint (value));
      break;
    case PROP_M2LNL_RELTOL:
      ncm_fit_set_m2lnL_reltol (fit, g_value_get_double (value));
      break;
    case PROP_M2LNL_ABSTOL:
      ncm_fit_set_m2lnL_abstol (fit, g_value_get_double (value));
      break;
    case PROP_PARAMS_RELTOL:
      ncm_fit_set_params_reltol (fit, g_value_get_double (value));
      break;
    case PROP_EQC:
      g_clear_pointer (&self->equality_constraints, g_ptr_array_unref);
      self->equality_constraints = g_value_dup_boxed (value);
      break;
    case PROP_EQC_TOT:
    {
      NcmVector *v = g_value_get_object (value);

      g_clear_pointer (&self->equality_constraints_tot, g_array_unref);

      if (v != NULL)
        self->equality_constraints_tot = ncm_vector_dup_array (v);
      else
        self->equality_constraints_tot = g_array_new (FALSE, FALSE, sizeof (gdouble));

      break;
    }
    case PROP_INEQC:
      g_clear_pointer (&self->inequality_constraints, g_ptr_array_unref);
      self->inequality_constraints = g_value_dup_boxed (value);
      break;
    case PROP_INEQC_TOT:
    {
      NcmVector *v = g_value_get_object (value);

      g_clear_pointer (&self->inequality_constraints_tot, g_array_unref);

      if (v != NULL)
        self->inequality_constraints_tot = ncm_vector_dup_array (v);
      else
        self->inequality_constraints_tot = g_array_new (FALSE, FALSE, sizeof (gdouble));

      break;
    }
    case PROP_SUBFIT:
      ncm_fit_set_sub_fit (fit, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fit_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFit *fit         = NCM_FIT (object);
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_return_if_fail (NCM_IS_FIT (object));

  switch (prop_id)
  {
    case PROP_LIKELIHOOD:
      g_value_set_object (value, self->lh);
      break;
    case PROP_MSET:
      g_value_set_object (value, self->mset);
      break;
    case PROP_STATE:
      g_value_set_object (value, self->fstate);
      break;
    case PROP_GRAD_TYPE:
      g_value_set_enum (value, self->grad.gtype);
      break;
    case PROP_MAXITER:
      g_value_set_uint (value, ncm_fit_get_maxiter (fit));
      break;
    case PROP_M2LNL_RELTOL:
      g_value_set_double (value, ncm_fit_get_m2lnL_reltol (fit));
      break;
    case PROP_M2LNL_ABSTOL:
      g_value_set_double (value, ncm_fit_get_m2lnL_abstol (fit));
      break;
    case PROP_PARAMS_RELTOL:
      g_value_set_double (value, ncm_fit_get_params_reltol (fit));
      break;
    case PROP_EQC:
      g_value_set_boxed (value, self->equality_constraints);
      break;
    case PROP_EQC_TOT:

      if (self->equality_constraints_tot->len > 0)
        g_value_take_object (value, ncm_vector_new_array (self->equality_constraints_tot));
      else
        g_value_set_object (value, NULL);

      break;
    case PROP_INEQC:
      g_value_set_boxed (value, self->inequality_constraints);
      break;
    case PROP_INEQC_TOT:

      if (self->inequality_constraints_tot->len > 0)
        g_value_take_object (value, ncm_vector_new_array (self->inequality_constraints_tot));
      else
        g_value_set_object (value, NULL);

      break;
    case PROP_SUBFIT:
      g_value_set_object (value, self->sub_fit);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fit_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_parent_class)->constructed (object);
  {
    NcmFit *fit         = NCM_FIT (object);
    NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
    gint n              = ncm_dataset_get_n (self->lh->dset);
    gint n_priors       = ncm_likelihood_priors_length_f (self->lh) +
                          ncm_likelihood_priors_length_m2lnL (self->lh);
    gint data_dof = ncm_dataset_get_dof (self->lh->dset);

    g_assert (ncm_dataset_all_init (self->lh->dset));

    if (!self->mset->valid_map)
      ncm_mset_prepare_fparam_map (self->mset);

    /*
     * It is no longer an error to fit 0 parameters, it just sets the value
     * of m2lnL in the fit object.
     *
     * if (ncm_mset_fparam_len (self->mset) == 0)
     * g_warning ("ncm_fit_new: mset object has 0 free parameters");
     *
     */

    {
      guint data_len   = n + n_priors;
      guint fparam_len = ncm_mset_fparam_len (self->mset);
      gint dof         = data_dof + n_priors - fparam_len;

      if (self->fstate == NULL)
        self->fstate = ncm_fit_state_new (data_len, fparam_len, dof,
                                          NCM_FIT_GET_CLASS (fit)->is_least_squares);
      else
        ncm_fit_state_set_all (self->fstate, data_len, fparam_len, dof,
                               NCM_FIT_GET_CLASS (fit)->is_least_squares);

      g_assert (data_len > 0);
    }
  }
}

static void
ncm_fit_dispose (GObject *object)
{
  NcmFit *fit         = NCM_FIT (object);
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_likelihood_clear (&self->lh);
  ncm_mset_clear (&self->mset);
  ncm_fit_state_clear (&self->fstate);

  g_clear_pointer (&self->equality_constraints, ncm_obj_array_unref);
  g_clear_pointer (&self->inequality_constraints, ncm_obj_array_unref);
  g_clear_pointer (&self->equality_constraints_tot, g_array_unref);
  g_clear_pointer (&self->inequality_constraints_tot, g_array_unref);

  ncm_fit_clear (&self->sub_fit);

  ncm_diff_clear (&self->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_parent_class)->dispose (object);
}

static void
ncm_fit_finalize (GObject *object)
{
  NcmFit *fit         = NCM_FIT (object);
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_timer_destroy (self->timer);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_parent_class)->finalize (object);
}

static void _ncm_fit_reset (NcmFit *fit);

static void
ncm_fit_class_init (NcmFitClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_constructed;
  object_class->set_property = &_ncm_fit_set_property;
  object_class->get_property = &_ncm_fit_get_property;
  object_class->dispose      = &ncm_fit_dispose;
  object_class->finalize     = &ncm_fit_finalize;

  klass->is_least_squares = FALSE;
  klass->reset            = &_ncm_fit_reset;

  g_object_class_install_property (object_class,
                                   PROP_LIKELIHOOD,
                                   g_param_spec_object ("likelihood",
                                                        NULL,
                                                        "Likelihood object",
                                                        NCM_TYPE_LIKELIHOOD,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "Model set object",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_STATE,
                                   g_param_spec_object ("state",
                                                        NULL,
                                                        "Fit state object",
                                                        NCM_TYPE_FIT_STATE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_GRAD_TYPE,
                                   g_param_spec_enum ("grad-type",
                                                      NULL,
                                                      "Differentiation method",
                                                      NCM_TYPE_FIT_GRAD_TYPE, NCM_FIT_GRAD_NUMDIFF_FORWARD,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MAXITER,
                                   g_param_spec_uint ("maxiter",
                                                      NULL,
                                                      "Maximum number of interations",
                                                      0, G_MAXUINT32, NCM_FIT_DEFAULT_MAXITER,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_M2LNL_RELTOL,
                                   g_param_spec_double ("m2lnL-reltol",
                                                        NULL,
                                                        "Relative tolarence in m2lnL",
                                                        0.0, G_MAXDOUBLE, NCM_FIT_DEFAULT_M2LNL_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M2LNL_ABSTOL,
                                   g_param_spec_double ("m2lnL-abstol",
                                                        NULL,
                                                        "Absolute tolarence in m2lnL",
                                                        0.0, G_MAXDOUBLE, NCM_FIT_DEFAULT_M2LNL_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PARAMS_RELTOL,
                                   g_param_spec_double ("params-reltol",
                                                        NULL,
                                                        "Relative tolarence in fitted parameters",
                                                        0.0, G_MAXDOUBLE, NCM_FIT_DEFAULT_PARAMS_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EQC,
                                   g_param_spec_boxed ("equality-constraints",
                                                       NULL,
                                                       "Equality constraints array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EQC_TOT,
                                   g_param_spec_object ("equality-constraints-tot",
                                                        NULL,
                                                        "Equality constraints tolerance",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INEQC,
                                   g_param_spec_boxed ("inequality-constraints",
                                                       NULL,
                                                       "Inequality constraints array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INEQC_TOT,
                                   g_param_spec_object ("inequality-constraints-tot",
                                                        NULL,
                                                        "Inequality constraints tolerance",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SUBFIT,
                                   g_param_spec_object ("sub-fit",
                                                        NULL,
                                                        "Subsidiary fit",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
_ncm_fit_reset (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  /*ncm_mset_prepare_fparam_map (self->mset);*/
  {
    guint n          = ncm_dataset_get_n (self->lh->dset);
    gint n_priors    = ncm_likelihood_priors_length_f (self->lh) + ncm_likelihood_priors_length_m2lnL (self->lh);
    guint data_len   = n + n_priors;
    guint fparam_len = ncm_mset_fparam_len (self->mset);
    gint data_dof    = ncm_dataset_get_dof (self->lh->dset);
    gint dof         = data_dof - fparam_len;

    g_assert (ncm_dataset_all_init (self->lh->dset));
    g_assert (data_len > 0);

    ncm_fit_state_set_all (self->fstate, data_len, fparam_len, dof,
                           NCM_FIT_GET_CLASS (fit)->is_least_squares);
    ncm_fit_state_reset (self->fstate);
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
 * Creates a new #NcmFit object.
 *
 * Returns: (transfer full): a new #NcmFit object.
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
    case NCM_FIT_TYPE_LEVMAR:
      return ncm_fit_levmar_new_by_name (lh, mset, gtype, algo_name);

      break;
#ifdef NUMCOSMO_HAVE_NLOPT
    case NCM_FIT_TYPE_NLOPT:
      return ncm_fit_nlopt_new_by_name (lh, mset, gtype, algo_name);

      break;
#endif /* NUMCOSMO_HAVE_NLOPT */
    default:
      g_error ("ncm_fit_new: fit-type not found %d, try to compile the library with the optional package NLOpt.", ftype);
      break;
  }
}

/**
 * ncm_fit_ref:
 * @fit: a #NcmFit
 *
 * Increases the reference count of @fit.
 *
 * Returns: (transfer full): @fit.
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
 * Returns: (transfer full): Copy of @fit.
 */
NcmFit *
ncm_fit_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return NCM_FIT_GET_CLASS (fit)->copy_new (fit, lh, mset, gtype);
}

/**
 * ncm_fit_dup:
 * @fit: a #NcmFit
 * @ser: a #NcmSerialize
 *
 * Duplicates the #NcmFit object duplicating all its contents.
 *
 * Returns: (transfer full): Duplicate of @fit.
 */
NcmFit *
ncm_fit_dup (NcmFit *fit, NcmSerialize *ser)
{
  return NCM_FIT (ncm_serialize_dup_obj (ser, G_OBJECT (fit)));
}

/**
 * ncm_fit_free:
 * @fit: a #NcmFit
 *
 * Atomically decrements the reference count of @fit by one. If the reference count drops to 0,
 * all memory allocated by @fit is released.
 *
 */
void
ncm_fit_free (NcmFit *fit)
{
  g_object_unref (fit);
}

/**
 * ncm_fit_clear:
 * @fit: a #NcmFit
 *
 * The reference count of @fit is decreased and the pointer is set to NULL.
 *
 */
void
ncm_fit_clear (NcmFit **fit)
{
  g_clear_object (fit);
}

/**
 * ncm_fit_set_sub_fit:
 * @fit: a #NcmFit
 * @sub_fit: a #NcmFit
 *
 * Sets a #NcmFit object to be used as subsidiary fit.
 *
 */
void
ncm_fit_set_sub_fit (NcmFit *fit, NcmFit *sub_fit)
{
  NcmFitPrivate *self     = ncm_fit_get_instance_private (fit);
  NcmFitPrivate *sub_self = ncm_fit_get_instance_private (sub_fit);

  ncm_fit_clear (&self->sub_fit);

  if (self->mset == sub_self->mset)
    g_error ("ncm_fit_set_sub_fit: cannot use the same mset in both fit and sub_fit.");

  if (!ncm_mset_is_subset (self->mset, sub_self->mset))
    g_error ("ncm_fit_set_sub_fit: sub_fit must contain a NcmMSet which is a subset of the NcmMSet of fit.");

  {
    guint fparams_len = ncm_mset_fparams_len (sub_self->mset);
    guint i;

    g_assert_cmpuint (fparams_len, >, 0);

    for (i = 0; i < fparams_len * 0; i++)
    {
      const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (sub_self->mset, i);

      if (ncm_mset_param_get_ftype (self->mset, pi->mid, pi->pid) != NCM_PARAM_TYPE_FIXED)
        g_error ("ncm_fit_set_sub_fit: parameter [%d %u] (%s) is free in both fit and sub_fit.",
                 pi->mid, pi->pid, ncm_mset_param_name (self->mset, pi->mid, pi->pid));
    }
  }

  self->sub_fit = ncm_fit_ref (sub_fit);
}

/**
 * ncm_fit_get_sub_fit:
 * @fit: a #NcmFit
 *
 * Gets a #NcmFit object to be used as subsidiary fit.
 *
 * Returns: (transfer full): a #NcmFit object.
 */
NcmFit *
ncm_fit_get_sub_fit (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return ncm_fit_ref (self->sub_fit);
}

static void _ncm_fit_m2lnL_grad_nd_fo (NcmFit *fit, NcmVector *grad);
static void _ncm_fit_m2lnL_grad_nd_ce (NcmFit *fit, NcmVector *grad);
static void _ncm_fit_m2lnL_grad_nd_ac (NcmFit *fit, NcmVector *grad);
static void _ncm_fit_m2lnL_val_grad_nd_fo (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
static void _ncm_fit_m2lnL_val_grad_nd_ce (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
static void _ncm_fit_m2lnL_val_grad_nd_ac (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
static void _ncm_fit_ls_J_nd_fo (NcmFit *fit, NcmMatrix *J);
static void _ncm_fit_ls_J_nd_ce (NcmFit *fit, NcmMatrix *J);
static void _ncm_fit_ls_J_nd_ac (NcmFit *fit, NcmMatrix *J);
static void _ncm_fit_ls_f_J_nd_fo (NcmFit *fit, NcmVector *f, NcmMatrix *J);
static void _ncm_fit_ls_f_J_nd_ce (NcmFit *fit, NcmVector *f, NcmMatrix *J);
static void _ncm_fit_ls_f_J_nd_ac (NcmFit *fit, NcmVector *f, NcmMatrix *J);

static NcmFitGrad _ncm_fit_grad_numdiff_forward = {
  NCM_FIT_GRAD_NUMDIFF_FORWARD,
  "Numerical differentiantion (forward)",
  &_ncm_fit_ls_J_nd_fo,
  &_ncm_fit_ls_f_J_nd_fo,
  &_ncm_fit_m2lnL_grad_nd_fo,
  &_ncm_fit_m2lnL_val_grad_nd_fo,
};

static NcmFitGrad _ncm_fit_grad_numdiff_central = {
  NCM_FIT_GRAD_NUMDIFF_CENTRAL,
  "Numerical differentiantion (central)",
  &_ncm_fit_ls_J_nd_ce,
  &_ncm_fit_ls_f_J_nd_ce,
  &_ncm_fit_m2lnL_grad_nd_ce,
  &_ncm_fit_m2lnL_val_grad_nd_ce,
};

static NcmFitGrad _ncm_fit_grad_numdiff_accurate = {
  NCM_FIT_GRAD_NUMDIFF_ACCURATE,
  "Numerical differentiantion (Richardson extrapolation)",
  &_ncm_fit_ls_J_nd_ac,
  &_ncm_fit_ls_f_J_nd_ac,
  &_ncm_fit_m2lnL_grad_nd_ac,
  &_ncm_fit_m2lnL_val_grad_nd_ac,
};

/**
 * ncm_fit_set_grad_type:
 * @fit: a #NcmLikelihood
 * @gtype: a #NcmFitGradType
 *
 * Sets the differentiation method to be used.
 *
 */
void
ncm_fit_set_grad_type (NcmFit *fit, NcmFitGradType gtype)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->grad.gtype = gtype;

  switch (gtype)
  {
    case NCM_FIT_GRAD_NUMDIFF_FORWARD:
      self->grad = _ncm_fit_grad_numdiff_forward;
      break;
    case NCM_FIT_GRAD_NUMDIFF_CENTRAL:
      self->grad = _ncm_fit_grad_numdiff_central;
      break;
    case NCM_FIT_GRAD_NUMDIFF_ACCURATE:
      self->grad = _ncm_fit_grad_numdiff_accurate;
      break;
    default:
      g_error ("Invalid gtype %d", gtype);
      break;
  }
}

/**
 * ncm_fit_set_maxiter:
 * @fit: a #NcmFit
 * @maxiter: maximum number of interations
 *
 * Sets the maximum number of iterations.
 *
 */
void
ncm_fit_set_maxiter (NcmFit *fit, guint maxiter)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->maxiter = maxiter;
}

/**
 * ncm_fit_get_grad_type:
 * @fit: a #NcmFit
 *
 * Gets the differentiation method to be used.
 *
 */
NcmFitGradType
ncm_fit_get_grad_type (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->grad.gtype;
}

/**
 * ncm_fit_get_maxiter:
 * @fit: a #NcmFit.
 *
 * Gets the maximum number of iterations.
 *
 * Returns: a integer (maxiter) that corresponds to the maximum number of iterations.
 */
guint
ncm_fit_get_maxiter (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->maxiter;
}

/**
 * ncm_fit_set_m2lnL_reltol:
 * @fit: a #NcmFit
 * @tol: relative tolarance
 *
 * Sets the relative tolerance for the m2lnL.
 *
 */
void
ncm_fit_set_m2lnL_reltol (NcmFit *fit, gdouble tol)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->m2lnL_reltol = tol;
}

/**
 * ncm_fit_get_m2lnL_reltol:
 * @fit: a #NcmFit
 *
 * Gets the relative tolerance for the m2lnL.
 *
 * Returns: the relative tolerance (double).
 */
gdouble
ncm_fit_get_m2lnL_reltol (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->m2lnL_reltol;
}

/**
 * ncm_fit_set_m2lnL_abstol:
 * @fit: a #NcmFit
 * @tol: absolute tolerance
 *
 * Sets the absolute tolerance for the m2lnL.
 *
 */
void
ncm_fit_set_m2lnL_abstol (NcmFit *fit, gdouble tol)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->m2lnL_abstol = tol;
}

/**
 * ncm_fit_get_m2lnL_abstol:
 * @fit: a #NcmFit
 *
 * Gets the absolute tolerance for the m2lnL.
 *
 * Returns: the absolute tolerance (double).
 */
gdouble
ncm_fit_get_m2lnL_abstol (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->m2lnL_abstol;
}

/**
 * ncm_fit_set_params_reltol:
 * @fit: a #NcmFit
 * @tol: relative tolerance
 *
 * Sets the relative tolerance for the fitted parameters.
 *
 */
void
ncm_fit_set_params_reltol (NcmFit *fit, gdouble tol)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->params_reltol = tol;
}

/**
 * ncm_fit_set_messages:
 * @fit: a #NcmFit
 * @mtype: a #NcmFitRunMsgs
 *
 * Sets the log level for the messages to be printed during the fit.
 *
 */
void
ncm_fit_set_messages (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->mtype = mtype;
}

/**
 * ncm_fit_get_params_reltol:
 * @fit: a #NcmFit
 *
 * Gets the relative tolerance for the fitted parameters.
 *
 * Returns: the relative tolerance (double).
 */
gdouble
ncm_fit_get_params_reltol (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->params_reltol;
}

/**
 * ncm_fit_get_messages:
 * @fit: a #NcmFit
 *
 * Gets the log level for the messages to be printed during the fit.
 *
 * Returns: a #NcmFitRunMsgs.
 */
NcmFitRunMsgs
ncm_fit_get_messages (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->mtype;
}

/**
 * ncm_fit_peek_mset:
 * @fit: a #NcmFit
 *
 * Peeks the #NcmMSet object.
 *
 * Returns: (transfer none): a #NcmMSet object.
 */
NcmMSet *
ncm_fit_peek_mset (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->mset;
}

/**
 * ncm_fit_peek_state:
 * @fit: a #NcmFit
 *
 * Peeks the #NcmFitState object.
 *
 * Returns: (transfer none): a #NcmFitState object.
 */
NcmFitState *
ncm_fit_peek_state (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->fstate;
}

/**
 * ncm_fit_peek_likelihood:
 * @fit: a #NcmFit
 *
 * Peeks the #NcmLikelihood object.
 *
 * Returns: (transfer none): a #NcmLikelihood object.
 */
NcmLikelihood *
ncm_fit_peek_likelihood (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->lh;
}

/**
 * ncm_fit_peek_diff:
 * @fit: a #NcmFit
 *
 * Peeks the #NcmDiff object.
 *
 * Returns: (transfer none): a #NcmDiff object.
 */
NcmDiff *
ncm_fit_peek_diff (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->diff;
}

/**
 * ncm_fit_params_set:
 * @fit: a #NcmFit
 * @i: the parameter index
 * @x: a double
 *
 * Sets the parameters vector.
 *
 */
void
ncm_fit_params_set (NcmFit *fit, guint i, const gdouble x)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_mset_fparam_set (self->mset, i, x);
  ncm_fit_params_update (fit);
}

/**
 * ncm_fit_params_set_vector:
 * @fit: a #NcmFit
 * @x: a #NcmVector
 *
 * Sets the parameters vector.
 *
 */
void
ncm_fit_params_set_vector (NcmFit *fit, NcmVector *x)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_mset_fparams_set_vector (self->mset, x);
  ncm_fit_params_update (fit);
}

/**
 * ncm_fit_params_set_vector_offset:
 * @fit: a #NcmFit
 * @x: a #NcmVector
 * @offset: FIXME
 *
 * Sets the parameters from vector @x starting at @offset.
 *
 */
void
ncm_fit_params_set_vector_offset (NcmFit *fit, NcmVector *x, guint offset)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_mset_fparams_set_vector_offset (self->mset, x, offset);
  ncm_fit_params_update (fit);
}

/**
 * ncm_fit_params_set_array:
 * @fit: a #NcmFit
 * @x: (in) (array) (element-type gdouble): an array of gdoubles
 *
 * Sets the parameters from array @x.
 *
 */
void
ncm_fit_params_set_array (NcmFit *fit, const gdouble *x)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_mset_fparams_set_array (self->mset, x);
  ncm_fit_params_update (fit);
}

/**
 * ncm_fit_params_set_gsl_vector: (skip)
 * @fit: a #NcmFit
 * @x: a gsl_vector
 *
 * Sets the parameters from a gsl_vector @x.
 *
 */
void
ncm_fit_params_set_gsl_vector (NcmFit *fit, const gsl_vector *x)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_mset_fparams_set_gsl_vector (self->mset, x);
  ncm_fit_params_update (fit);
}

/**
 * ncm_fit_params_update:
 * @fit: a #NcmFit
 *
 * Updates the parameters vector.
 *
 */
void
ncm_fit_params_update (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  if (self->sub_fit != NULL)
  {
    NcmFitRunMsgs mlevel = self->mtype == NCM_FIT_RUN_MSGS_FULL ? NCM_FIT_RUN_MSGS_SIMPLE : NCM_FIT_RUN_MSGS_NONE;

    ncm_fit_run (self->sub_fit, mlevel);
  }
}

/**
 * ncm_fit_add_equality_constraint:
 * @fit: a #NcmFit
 * @func: a #NcmMSetFunc
 * @tot: tolerance
 *
 * Adds an equality constraint with the function @func and the tolerance @tot.
 *
 */
void
ncm_fit_add_equality_constraint (NcmFit *fit, NcmMSetFunc *func, const gdouble tot)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_obj_array_add (self->equality_constraints, G_OBJECT (func));
  g_array_append_val (self->equality_constraints_tot, tot);
}

/**
 * ncm_fit_add_inequality_constraint:
 * @fit: a #NcmFit
 * @func: a #NcmMSetFunc
 * @tot: tolerance
 *
 * Adds an inequality constraint with the function @func and the tolerance @tot.
 *
 */
void
ncm_fit_add_inequality_constraint (NcmFit *fit, NcmMSetFunc *func, const gdouble tot)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_obj_array_add (self->inequality_constraints, G_OBJECT (func));
  g_array_append_val (self->inequality_constraints_tot, tot);
}

/**
 * ncm_fit_remove_equality_constraints:
 * @fit: a #NcmFit
 *
 * Removes all equality constraints.
 *
 */
void
ncm_fit_remove_equality_constraints (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_ptr_array_set_size (self->equality_constraints, 0);
  g_array_set_size (self->equality_constraints_tot, 0);
}

/**
 * ncm_fit_remove_inequality_constraints:
 * @fit: a #NcmFit
 *
 * Removes all inequality constraints.
 *
 */
void
ncm_fit_remove_inequality_constraints (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_ptr_array_set_size (self->inequality_constraints, 0);
  g_array_set_size (self->inequality_constraints_tot, 0);
}

/**
 * ncm_fit_equality_constraints_len:
 * @fit: a #NcmFit
 *
 * Gets the number of equality constraints.
 *
 * Returns: the number of equality constraints.
 */
guint
ncm_fit_equality_constraints_len (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->equality_constraints->len;
}

/**
 * ncm_fit_inequality_constraints_len:
 * @fit: a #NcmFit
 *
 * Gets the number of inequality constraints.
 *
 * Returns: the number of inequality constraints.
 */
guint
ncm_fit_inequality_constraints_len (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  return self->inequality_constraints->len;
}

/**
 * ncm_fit_get_equality_constraint:
 * @fit: a #NcmFit
 * @i: index of the equality constraint
 * @func: (out) (transfer none): a #NcmMSetFunc
 * @tot: (out): a double
 *
 * Gets the equality constraint at index @i.
 *
 */
void
ncm_fit_get_equality_constraint (NcmFit *fit, guint i, NcmMSetFunc **func, gdouble *tot)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_assert_cmpuint (i, <, self->equality_constraints->len);
  g_assert (func != NULL);
  g_assert (tot != NULL);

  *func = NCM_MSET_FUNC (g_ptr_array_index (self->equality_constraints, i));
  *tot  = g_array_index (self->equality_constraints_tot, gdouble, i);
}

/**
 * ncm_fit_get_inequality_constraint:
 * @fit: a #NcmFit
 * @i: index of the inequality constraint
 * @func: (out) (transfer none): a #NcmMSetFunc
 * @tot: (out): a double
 *
 * Gets the inequality constraint at index @i.
 *
 */
void
ncm_fit_get_inequality_constraint (NcmFit *fit, guint i, NcmMSetFunc **func, gdouble *tot)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_assert_cmpuint (i, <, self->inequality_constraints->len);
  g_assert (func != NULL);
  g_assert (tot != NULL);

  *func = NCM_MSET_FUNC (g_ptr_array_index (self->inequality_constraints, i));
  *tot  = g_array_index (self->inequality_constraints_tot, gdouble, i);
}

/**
 * ncm_fit_ls_covar:
 * @fit: a #NcmFit
 *
 * Computes the covariance matrix using the least squares method,
 * and fills up the internal structure matrix.
 *
 */
void
ncm_fit_ls_covar (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmMatrix *J        = ncm_fit_state_peek_J (self->fstate);
  NcmMatrix *covar    = ncm_fit_state_peek_covar (self->fstate);

  ncm_fit_ls_J (fit, J);
  gsl_multifit_covar (ncm_matrix_gsl (J), 0.0,
                      ncm_matrix_gsl (covar));

  ncm_fit_state_set_has_covar (self->fstate, TRUE);
}

/**
 * ncm_fit_covar_fparam_var:
 * @fit: a #NcmFit
 * @fpi: index of a free parameter
 *
 * Computes the variance of the fitted parameter @fpi.
 * This index refers to the list of all FREE parameters set in the MSet.
 *
 * See also the similar function ncm_fit_covar_var() to which one has to provide
 * the respective model of the parameter.
 *
 * Returns: the variance of the fitted parameter @fpi
 */
gdouble
ncm_fit_covar_fparam_var (NcmFit *fit, guint fpi)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_assert (ncm_fit_state_has_covar (self->fstate));

  {
    NcmMatrix *covar = ncm_fit_state_peek_covar (self->fstate);

    return ncm_matrix_get (covar, fpi, fpi);
  }
}

/**
 * ncm_fit_covar_fparam_sd:
 * @fit: a #NcmFit
 * @fpi: index of a free parameter
 *
 * Computes the standard deviation of the fitted parameter @fpi.
 * This index refers to the list of all FREE parameters set in the MSet.
 *
 * See also the similar function ncm_fit_covar_sd() to which one has to provide
 * the respective model of the parameter.
 *
 * Returns: the standard deviation of the fitted parameter @fpi
 */
gdouble
ncm_fit_covar_fparam_sd (NcmFit *fit, guint fpi)
{
  return sqrt (ncm_fit_covar_fparam_var (fit, fpi));
}

/**
 * ncm_fit_covar_fparam_cov:
 * @fit: a #NcmFit
 * @fpi1: index of a free parameter
 * @fpi2: index of a free parameter
 *
 * Computes the covariance between the fitted parameters @fpi1 and @fpi2.
 * These indices refers to the list of all FREE parameters set in the MSet.
 *
 * See also the similar function ncm_fit_covar_cov() to which one has to provide
 * the respective models of the parameters.
 *
 * Returns: the covariance between the fitted parameters @pdi1 and @fpdi2
 */
gdouble
ncm_fit_covar_fparam_cov (NcmFit *fit, guint fpi1, guint fpi2)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_assert (ncm_fit_state_has_covar (self->fstate));

  {
    NcmMatrix *covar = ncm_fit_state_peek_covar (self->fstate);

    return ncm_matrix_get (covar, fpi1, fpi2);
  }
}

/**
 * ncm_fit_covar_fparam_cor:
 * @fit: a #NcmFit
 * @fpi1: index of a free parameter
 * @fpi2: index of a free parameter
 *
 * Computes the correlation between the fitted parameters @fpi1 and @fpi2.
 * These indices refers to the list of all FREE parameters set in the MSet.
 *
 * See also the similar function ncm_fit_covar_cor() to which one has to provide
 * the respective models of the parameters.
 *
 * Returns: the correlation between the fitted parameters @pdi1 and @fpdi2
 */
gdouble
ncm_fit_covar_fparam_cor (NcmFit *fit, guint fpi1, guint fpi2)
{
  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2) / (ncm_fit_covar_fparam_sd (fit, fpi1) * ncm_fit_covar_fparam_sd (fit, fpi2));
}

/**
 * ncm_fit_covar_var:
 * @fit: a #NcmFit
 * @mid: a #NcmModelID
 * @pid: the parameter's index of the model @mid (integer)
 *
 * Computes the variance of the fitted parameter @pid of the model @mid.
 *
 * Returns: the variance of the fitted parameter @pid
 */
gdouble
ncm_fit_covar_var (NcmFit *fit, NcmModelID mid, guint pid)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  gint fpi            = ncm_mset_fparam_get_fpi (self->mset, mid, pid);

  if (fpi < 0)
    g_error ("Parameter (%d:%u) was not fit.", mid, pid);

  return ncm_fit_covar_fparam_var (fit, fpi);
}

/**
 * ncm_fit_covar_sd:
 * @fit: a #NcmFit
 * @mid: a #NcmModelID
 * @pid: the parameter's index of the model @mid (integer)
 *
 * Computes the standard deviation of the fitted parameter @pid of the model @mid.
 *
 * Returns: the standard deviation of the fitted parameter @pid
 */
gdouble
ncm_fit_covar_sd (NcmFit *fit, NcmModelID mid, guint pid)
{
  return sqrt (ncm_fit_covar_var (fit, mid, pid));
}

/**
 * ncm_fit_covar_cov:
 * @fit: a #NcmFit
 * @mid1: a #NcmModelID
 * @pid1: the parameter's index of the model @mid1 (integer)
 * @mid2: a #NcmModelID
 * @pid2: the parameter's index of the model @mid1 (integer)
 *
 * Computes the covariance between the parameters @pid1 and @pid2 of the models
 * @mid1 and @mid2, respectively.
 *
 * Returns: the covariance between @pid1 and @pid2
 */
gdouble
ncm_fit_covar_cov (NcmFit *fit, NcmModelID mid1, guint pid1, NcmModelID mid2, guint pid2)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  gint fpi1           = ncm_mset_fparam_get_fpi (self->mset, mid1, pid1);
  gint fpi2           = ncm_mset_fparam_get_fpi (self->mset, mid2, pid2);

  if ((fpi1 < 0) || (fpi1 < 0))
    g_error ("Parameters (%d:%u, %d:%u) were not fit.", mid1, pid1, mid2, pid2);

  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2);
}

/**
 * ncm_fit_covar_cor:
 * @fit: a #NcmFit
 * @mid1: a #NcmModelID
 * @pid1: the parameter's index of the model @mid1 (integer)
 * @mid2: a #NcmModelID
 * @pid2: the parameter's index of the model @mid1 (integer)
 *
 * Computes the correlation between the parameters @pid1 and @pid2 of the models
 * @mid1 and @mid2, respectively.
 *
 * Returns: the correlation between @pid1 and @pid2
 */
gdouble
ncm_fit_covar_cor (NcmFit *fit, NcmModelID mid1, guint pid1, NcmModelID mid2, guint pid2)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  gint fpi1           = ncm_mset_fparam_get_fpi (self->mset, mid1, pid1);
  gint fpi2           = ncm_mset_fparam_get_fpi (self->mset, mid2, pid2);

  if ((fpi1 < 0) || (fpi1 < 0))
    g_error ("Parameters (%d:%u, %d:%u) were not fit.", mid1, pid1, mid2, pid2);

  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2) / (ncm_fit_covar_fparam_sd (fit, fpi1) * ncm_fit_covar_fparam_sd (fit, fpi2));
}

gboolean
_ncm_fit_run_empty (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  gdouble m2lnL_curval;

  self->mtype = mtype;

  ncm_fit_m2lnL_val (fit, &m2lnL_curval);
  ncm_fit_state_set_m2lnL_curval (self->fstate, m2lnL_curval);
  ncm_fit_log_step (fit);

  ncm_fit_state_set_m2lnL_prec (self->fstate, 0.0);
  ncm_fit_state_set_params_prec (self->fstate, 0.0);

  ncm_fit_state_set_has_covar (self->fstate, FALSE);
  ncm_fit_state_set_is_best_fit (self->fstate, TRUE);
  ncm_fit_state_set_elapsed_time (self->fstate, 0.0);
  ncm_fit_state_set_niter (self->fstate, 0);
  ncm_fit_state_set_func_eval (self->fstate, 0);
  ncm_fit_state_set_grad_eval (self->fstate, 0);

  return TRUE;
}

/**
 * ncm_fit_reset:
 * @fit: a #NcmFit
 *
 * Resets the fit.
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
 * Computes the minimization.
 *
 * Returns: TRUE if the minimization went through.
 */
gboolean
ncm_fit_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  gboolean run;
  gdouble m2lnL_i;

  self->mtype = mtype;

  ncm_fit_reset (fit);
  ncm_fit_log_start (fit);
  g_timer_start (self->timer);

  ncm_fit_m2lnL_val (fit, &m2lnL_i);

  if (gsl_finite (m2lnL_i))
  {
    if (ncm_mset_fparam_len (self->mset) == 0)
      run = _ncm_fit_run_empty (fit, mtype);
    else
      run = NCM_FIT_GET_CLASS (fit)->run (fit, mtype);
  }
  else
  {
    g_warning ("ncm_fit_run: initial point provides m2lnL = % 22.15g, giving up.", m2lnL_i);
    run = FALSE;
  }

  ncm_fit_state_set_elapsed_time (self->fstate, g_timer_elapsed (self->timer, NULL));
  ncm_fit_state_set_is_best_fit (self->fstate, run);

  ncm_fit_log_end (fit);

  return run;
}

/**
 * ncm_fit_run_restart:
 * @fit: a #NcmFit
 * @mtype: a #NcmFitRunMsgs
 * @abstol: restart absolute tolerance
 * @reltol: restart relative tolarence
 * @save_mset: (nullable): the #NcmMSet used to save progress
 * @mset_file: (nullable): the file name to save progress
 *
 * Re-runs the fit until the difference between fits are less
 * than the required tolerance, i.e.,
 * $$ m2lnL_{i-1} - m2lnL_i < \mathrm{abstol} + \mathrm{reltol}\vert m2lnL_{i-1}\vert. $$
 *
 */
void
ncm_fit_run_restart (NcmFit *fit, NcmFitRunMsgs mtype, const gdouble abstol, const gdouble reltol, NcmMSet *save_mset, const gchar *mset_file)
{
  NcmFitPrivate *self    = ncm_fit_get_instance_private (fit);
  gboolean save_progress = (mset_file == NULL) ? FALSE : TRUE;
  NcmMSet *mset_out      = (save_mset == NULL) ? self->mset : save_mset;
  NcmSerialize *ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  self->mtype = mtype;

  if (ncm_mset_fparam_len (self->mset) == 0)
  {
    ncm_fit_reset (fit);
    ncm_fit_log_start (fit);

    g_timer_start (self->timer);
    _ncm_fit_run_empty (fit, mtype);

    ncm_fit_state_set_elapsed_time (self->fstate, g_timer_elapsed (self->timer, NULL));
    ncm_fit_state_set_is_best_fit (self->fstate, TRUE);

    ncm_fit_log_end (fit);

    return;
  }
  else
  {
    gboolean restart;
    gdouble last_m2lnL = 0.0, m2lnL = 0.0;
    gint n = 0;

    ncm_fit_m2lnL_val (fit, &last_m2lnL);

    do {
      ncm_fit_run (fit, mtype);

      m2lnL   = ncm_fit_state_get_m2lnL_curval (self->fstate);
      restart = (last_m2lnL - m2lnL) >= (abstol + reltol * fabs (last_m2lnL));

      if (self->mset != mset_out)
        ncm_mset_param_set_mset (mset_out, self->mset);

      if (save_progress)
        ncm_mset_save (mset_out, ser, mset_file, TRUE);

      if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
      {
        ncm_cfg_msg_sepa ();
        g_message ("# Restarting:              %s\n", restart ? "yes" : "no");
        g_message ("#  - absolute improvement: %-22.15g\n", (last_m2lnL - m2lnL));
        g_message ("#  - relative improvement: %-22.15g\n", (last_m2lnL - m2lnL) / last_m2lnL);
        g_message ("#  - m2lnL_%-3d:            %-22.15g\n", n,     last_m2lnL);
        g_message ("#  - m2lnL_%-3d:            %-22.15g\n", n + 1, m2lnL);
      }

      last_m2lnL = m2lnL;

      n++;
    } while (restart);
  }

  ncm_serialize_free (ser);
}

/**
 * ncm_fit_is_least_squares:
 * @fit: a #NcmFit
 *
 * Indicates if the least squares fitting is being used (TRUE) or not (FALSE).
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
 * Gets the fit object description.
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
 * This function prints in the log the initial state.
 *
 */
void
ncm_fit_log_start (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# Model fitting. Interating using:\n");
    g_message ("#  - solver:            %s\n", ncm_fit_get_desc (fit));
    g_message ("#  - differentiation:   %s\n", self->grad.diff_name);

    if (self->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("#");
  }
}

/**
 * ncm_fit_log_step_error:
 * @fit: a #NcmFit
 * @strerror: error message
 * @...: arguments
 *
 * This function prints in the log the error message.
 *
 */
void
ncm_fit_log_step_error (NcmFit *fit, const gchar *strerror, ...)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  va_list ap;

  va_start (ap, strerror);

  if (self->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
  {
    gchar *errmsg = g_strdup_vprintf (strerror, ap);

    g_message ("\n#  [%s] error = %s\n#", ncm_fit_get_desc (fit), errmsg);
    g_free (errmsg);
  }
  else if (self->mtype == NCM_FIT_RUN_MSGS_FULL)
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
 *
 */
void
ncm_fit_log_end (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_timer_stop (self->timer);

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (self->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("\n");

    g_message ("#  Minimum found with precision: |df|/f = % 8.5e and |dx| = % 8.5e\n",
               ncm_fit_state_get_m2lnL_prec (self->fstate),
               ncm_fit_state_get_params_prec (self->fstate));
  }

  ncm_fit_log_state (fit);

  return;
}

/**
 * ncm_fit_log_state:
 * @fit: a #NcmFit
 *
 * This function prints in the log the current state.
 *
 */
void
ncm_fit_log_state (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  guint i;

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    gdouble elap_sec = g_timer_elapsed (self->timer, NULL);
    gulong elap_min  = elap_sec / 60;
    gulong elap_hour = elap_min / 60;
    gulong elap_day  = elap_hour / 24;

    elap_sec  = fmod (elap_sec, 60);
    elap_min  = elap_min % 60;
    elap_hour =  elap_hour % 24;
    g_message ("#  Elapsed time: %02lu days, %02lu:%02lu:%010.7f\n", elap_day, elap_hour, elap_min, elap_sec);
    g_message ("#  iteration            [%06d]\n", ncm_fit_state_get_niter (self->fstate));
    g_message ("#  function evaluations [%06d]\n", ncm_fit_state_get_func_eval (self->fstate));
    g_message ("#  gradient evaluations [%06d]\n", ncm_fit_state_get_grad_eval (self->fstate));
    g_message ("#  degrees of freedom   [%06d]\n", ncm_fit_state_get_dof (self->fstate));

    if (self->lh->m2lnL_v != NULL)
    {
      g_message ("#  m2lnL     = % 20.15g ( ", ncm_fit_state_get_m2lnL_curval (self->fstate));

      for (i = 0; i < ncm_vector_len (self->lh->m2lnL_v); i++)
      {
        g_message ("% 13.8g ", ncm_vector_get (self->lh->m2lnL_v, i));
      }

      g_message (")\n");
    }
    else
    {
      g_message ("#  m2lnL     = % 20.15g\n", ncm_fit_state_get_m2lnL_curval (self->fstate));
    }

    g_message ("#  Fit parameters:\n#    ");

    for (i = 0; i < ncm_mset_fparam_len (self->mset); i++)
      g_message ("% -20.15g ", ncm_mset_fparam_get (self->mset, i));

    g_message ("\n");
  }

  return;
}

/**
 * ncm_fit_log_step:
 * @fit: a #NcmFit
 *
 * This function prints in the log one step of the minimization.
 *
 */
void
ncm_fit_log_step (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  if (self->mtype == NCM_FIT_RUN_MSGS_FULL)
  {
    ncm_fit_log_state (fit);
  }
  else if (self->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
  {
    const guint niter = ncm_fit_state_get_niter (self->fstate);

    if ((niter < 10) || ((niter % 10) == 0))
      ncm_message (".");
  }

  return;
}

/**
 * ncm_fit_log_finish:
 * @fit: a #NcmFit
 *
 * This function prints in the log the final state.
 *
 */
void
ncm_fit_log_finish (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    const gdouble m2lnL     = ncm_fit_state_get_m2lnL_curval (self->fstate);
    const gint dof          = ncm_fit_state_get_dof (self->fstate);
    const gdouble m2lnL_dof = m2lnL / dof;

    g_message ("#  m2lnL/dof = %20.15g\n", m2lnL_dof);
    g_message ("#  |m2lnL-dof|/sqrt(2*dof) = %20.15g,\n", fabs (m2lnL_dof) / sqrt (2.0 * dof));
    g_message ("#  GoF_tt = %4.2f%% = (%4.2f + %4.2f)%%; GoF = %4.2f%%\n",
               gsl_cdf_chisq_Q (dof + fabs (dof - m2lnL), dof) * 100.0 +
               gsl_cdf_chisq_P (dof - fabs (dof - m2lnL), dof) * 100.0,
               gsl_cdf_chisq_P (dof - fabs (dof - m2lnL), dof) * 100.0,
               gsl_cdf_chisq_Q (dof + fabs (dof - m2lnL), dof) * 100.0,
               gsl_cdf_chisq_Q (m2lnL, dof) * 100.0);
  }

  return;
}

/**
 * ncm_fit_log_info:
 * @fit: a #NcmFit
 *
 * Prints to the log the data set and the model set.
 *
 */
void
ncm_fit_log_info (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_dataset_log_info (self->lh->dset);
  ncm_mset_pretty_log (self->mset);

  return;
}

/**
 * ncm_fit_log_covar:
 * @fit: a #NcmFit
 *
 * Prints to the log file the names and indices of the fitted parameters, their best-fit
 * values, standard deviations and correlation matrix.
 */
void
ncm_fit_log_covar (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_assert (ncm_fit_state_has_covar (self->fstate));
  ncm_mset_fparams_log_covar (self->mset, ncm_fit_state_peek_covar (self->fstate));

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
  NcmFitPrivate *self   = ncm_fit_get_instance_private (fit);
  guint name_size       = ncm_mset_max_param_name (self->mset);
  guint free_params_len = ncm_mset_fparam_len (self->mset);
  guint i, j;

  if (header != NULL)
    fprintf (out, "# %s\n# ", header);
  else
    fprintf (out, "# ");

  for (i = 0; i < free_params_len; i++)
  {
    const gchar *pname = ncm_mset_fparam_name (self->mset, i);

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
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_dataset_m2lnL_val (self->lh->dset, self->mset, data_m2lnL);
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
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_likelihood_priors_m2lnL_val (self->lh, self->mset, priors_m2lnL);
}

/**
 * ncm_fit_m2lnL_val:
 * @fit: a #NcmFit
 * @m2lnL: (out): minus two times the logarithm base e of the likelihood.
 *
 * Computes minus two times the logarithm base e of the likelihood.
 *
 */
void
ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_likelihood_m2lnL_val (self->lh, self->mset, m2lnL);
  ncm_fit_state_add_func_eval (self->fstate, 1);
}

/**
 * ncm_fit_ls_f:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 *
 * Computes the residuals vector.
 *
 */
void
ncm_fit_ls_f (NcmFit *fit, NcmVector *f)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_likelihood_leastsquares_f (self->lh, self->mset, f);
  ncm_fit_state_add_func_eval (self->fstate, 1);
}

/**
 * ncm_fit_m2lnL_grad:
 * @fit: a #NcmFit
 * @df: a #NcmVector
 *
 * Computes the gradient of the minus two times the logarithm base e of the likelihood.
 *
 */
void
ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->grad.m2lnL_grad (fit, df);
}

static gdouble _ncm_fit_numdiff_m2lnL_val (NcmVector *x, gpointer user_data);
static void _ncm_fit_numdiff_ls_f (NcmVector *x, NcmVector *y, gpointer user_data);

static void
_ncm_fit_m2lnL_grad_nd_fo (NcmFit *fit, NcmVector *grad)
{
  gdouble m2lnL;

  _ncm_fit_m2lnL_val_grad_nd_fo (fit, &m2lnL, grad);
}

static void
_ncm_fit_m2lnL_grad_nd_ce (NcmFit *fit, NcmVector *grad)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  guint fparam_len    = ncm_mset_fparam_len (self->mset);
  guint i;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p       = ncm_mset_fparam_get (self->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (self->mset, i));
    const gdouble h       = p_scale * GSL_ROOT3_DBL_EPSILON;
    const gdouble pph     = p + h;
    const gdouble pmh     = p - h;
    const gdouble twoh    = pph - pmh;
    const gdouble one_2h  = 1.0 / twoh;
    gdouble m2lnL_pph, m2lnL_pmh;

    ncm_fit_params_set (fit, i, pph);
    ncm_likelihood_m2lnL_val (self->lh, self->mset, &m2lnL_pph);
    /*ncm_fit_m2lnL_val (fit, &m2lnL_pph);*/

    ncm_fit_params_set (fit, i, pmh);
    ncm_likelihood_m2lnL_val (self->lh, self->mset, &m2lnL_pmh);
    /*ncm_fit_m2lnL_val (fit, &m2lnL_pmh);*/

    ncm_vector_set (grad, i, (m2lnL_pph - m2lnL_pmh) * one_2h);
    ncm_fit_params_set (fit, i, p);
  }

  ncm_fit_state_add_grad_eval (self->fstate, 1);
}

static void
_ncm_fit_m2lnL_grad_nd_ac (NcmFit *fit, NcmVector *grad)
{
  NcmFitPrivate *self         = ncm_fit_get_instance_private (fit);
  const guint free_params_len = ncm_mset_fparams_len (self->mset);
  GArray *x_a                 = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x                = NULL;
  GArray *grad_a              = NULL;

  g_array_set_size (x_a, free_params_len);
  x = ncm_vector_new_array (x_a);

  ncm_mset_fparams_get_vector (self->mset, x);
  grad_a = ncm_diff_rf_d1_N_to_1 (self->diff, x_a, _ncm_fit_numdiff_m2lnL_val, fit, NULL);

  ncm_vector_set_array (grad, grad_a);
  ncm_mset_fparams_set_vector (self->mset, x);

  g_array_unref (x_a);
  g_array_unref (grad_a);
  ncm_vector_free (x);

  ncm_fit_state_add_grad_eval (self->fstate, 1);
}

/**
 * ncm_fit_m2lnL_val_grad:
 * @fit: a #NcmFit
 * @result: (out): the minus two times the logarithm base e of the likelihood
 * @df: a #NcmVector
 *
 * Computes the minus two times the logarithm base e of the likelihood and its
 * gradient.
 *
 */
void
ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->grad.m2lnL_val_grad (fit, result, df);
}

static void
_ncm_fit_m2lnL_val_grad_nd_fo (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  guint fparam_len    = ncm_mset_fparam_len (self->mset);
  guint i;

  ncm_fit_m2lnL_val (fit, m2lnL);

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p       = ncm_mset_fparam_get (self->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (self->mset, i));
    const gdouble htilde  = GSL_SQRT_DBL_EPSILON * p_scale;
    const gdouble pph     = p + htilde;
    const gdouble h       = pph - p;
    const gdouble one_h   = 1.0 / h;
    gdouble m2lnL_pph;

    ncm_fit_params_set (fit, i, pph);
    ncm_fit_m2lnL_val (fit, &m2lnL_pph);
    ncm_vector_set (grad, i, (m2lnL_pph - *m2lnL) * one_h);
    ncm_fit_params_set (fit, i, p);
  }

  ncm_fit_state_add_func_eval (self->fstate, 1);
  ncm_fit_state_add_grad_eval (self->fstate, 1);
}

static void
_ncm_fit_m2lnL_val_grad_nd_ce (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  ncm_fit_m2lnL_val (fit, m2lnL);
  _ncm_fit_m2lnL_grad_nd_ce (fit, grad);
}

static void
_ncm_fit_m2lnL_val_grad_nd_ac (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  ncm_fit_m2lnL_val (fit, m2lnL);
  _ncm_fit_m2lnL_grad_nd_ac (fit, grad);
}

/**
 * ncm_fit_ls_J:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * Computes the Jacobian matrix for the least squares problem.
 *
 */
void
ncm_fit_ls_J (NcmFit *fit, NcmMatrix *J)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->grad.ls_J (fit, J);
}

static void
_ncm_fit_ls_J_nd_fo (NcmFit *fit, NcmMatrix *J)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  guint fparam_len    = ncm_mset_fparam_len (self->mset);
  NcmVector *f        = ncm_fit_state_peek_f (self->fstate);
  guint i;

  ncm_fit_ls_f (fit, f);

  for (i = 0; i < fparam_len; i++)
  {
    NcmVector *J_col_i    = ncm_matrix_get_col (J, i);
    const gdouble p       = ncm_mset_fparam_get (self->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (self->mset, i));
    const gdouble htilde  = p_scale * GSL_SQRT_DBL_EPSILON;
    const gdouble pph     = p + htilde;
    const gdouble h       = pph - p;
    const gdouble one_h   = 1.0 / h;
    guint k;

    ncm_fit_params_set (fit, i, pph);
    ncm_fit_ls_f (fit, J_col_i);

    for (k = 0; k < ncm_vector_len (f); k++)
      ncm_vector_set (J_col_i, k,
                      (ncm_vector_get (J_col_i, k) - ncm_vector_get (f, k)) * one_h
                     );

    ncm_fit_params_set (fit, i, p);
    ncm_vector_free (J_col_i);
  }

  ncm_fit_state_add_grad_eval (self->fstate, 1);
}

static void
_ncm_fit_ls_J_nd_ce (NcmFit *fit, NcmMatrix *J)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  guint fparam_len    = ncm_mset_fparam_len (self->mset);
  NcmVector *f        = ncm_fit_state_peek_f (self->fstate);
  guint i;

  for (i = 0; i < fparam_len; i++)
  {
    NcmVector *J_col_i    = ncm_matrix_get_col (J, i);
    const gdouble p       = ncm_mset_fparam_get (self->mset, i);
    const gdouble p_scale = GSL_MAX (fabs (p), ncm_mset_fparam_get_scale (self->mset, i));
    const gdouble h       = p_scale * GSL_ROOT3_DBL_EPSILON;
    const gdouble pph     = p + h;
    const gdouble pmh     = p - h;
    const gdouble twoh    = pph - pmh;
    const gdouble one_2h  = 1.0 / twoh;
    guint k;

    ncm_fit_params_set (fit, i, pph);
    ncm_fit_ls_f (fit, J_col_i);

    ncm_fit_params_set (fit, i, pmh);
    ncm_fit_ls_f (fit, f);

    for (k = 0; k < ncm_vector_len (f); k++)
      ncm_vector_set (J_col_i, k,
                      (ncm_vector_get (J_col_i, k) - ncm_vector_get (f, k)) * one_2h
                     );

    ncm_fit_params_set (fit, i, p);
    ncm_vector_free (J_col_i);
  }

  ncm_fit_state_add_grad_eval (self->fstate, 1);
}

static void
_ncm_fit_ls_J_nd_ac (NcmFit *fit, NcmMatrix *J)
{
  NcmFitPrivate *self    = ncm_fit_get_instance_private (fit);
  const guint fparam_len = ncm_mset_fparam_len (self->mset);
  const guint data_len   = ncm_fit_state_get_data_len (self->fstate);
  GArray *x_a            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x           = NULL;
  GArray *J_a            = NULL;

  g_array_set_size (x_a, fparam_len);
  x = ncm_vector_new_array (x_a);

  ncm_mset_fparams_get_vector (self->mset, x);
  J_a = ncm_diff_rf_d1_N_to_M (self->diff, x_a, data_len, _ncm_fit_numdiff_ls_f, fit, NULL);

  ncm_matrix_set_from_array (J, J_a);
  ncm_mset_fparams_set_vector (self->mset, x);

  g_array_unref (x_a);
  g_array_unref (J_a);
  ncm_vector_free (x);

  ncm_fit_state_add_grad_eval (self->fstate, 1);
}

/**
 * ncm_fit_ls_f_J:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * Computes the residuals vector and the Jacobian matrix for the least squares problem.
 *
 */
void
ncm_fit_ls_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  self->grad.ls_f_J (fit, f, J);
}

static void
_ncm_fit_ls_f_J_nd_fo (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmVector *fit_f    = ncm_fit_state_peek_f (self->fstate);

  _ncm_fit_ls_J_nd_fo (fit, J);
  ncm_vector_memcpy (f, fit_f);
}

static void
_ncm_fit_ls_f_J_nd_ce (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_fit_ls_f (fit, f);
  _ncm_fit_ls_J_nd_ce (fit, J);
}

static void
_ncm_fit_ls_f_J_nd_ac (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_fit_ls_f (fit, f);
  _ncm_fit_ls_J_nd_ac (fit, J);
}

static gdouble
_ncm_fit_numdiff_m2lnL_val (NcmVector *x, gpointer user_data)
{
  NcmFit *fit         = NCM_FIT (user_data);
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  gdouble res         = 0.0;

  ncm_mset_fparams_set_vector (self->mset, x);

  ncm_likelihood_m2lnL_val (self->lh, self->mset, &res);

  return res;
}

static void
_ncm_fit_numdiff_ls_f (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmFit *fit         = NCM_FIT (user_data);
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  ncm_mset_fparams_set_vector (self->mset, x);

  ncm_fit_ls_f (fit, y);
}

/**
 * ncm_fit_numdiff_m2lnL_hessian:
 * @fit: a #NcmFit
 * @H: a #NcmMatrix
 * @reltol: relative tolerance.
 *
 * Calculates the Hessian of $-2\ln(L)$ with respect to the free parameters.
 *
 */
void
ncm_fit_numdiff_m2lnL_hessian (NcmFit *fit, NcmMatrix *H, gdouble reltol)
{
  NcmFitPrivate *self         = ncm_fit_get_instance_private (fit);
  const guint free_params_len = ncm_mset_fparams_len (self->mset);
  GArray *x_a                 = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x                = NULL;
  GArray *H_a                 = NULL;

  g_array_set_size (x_a, free_params_len);
  x = ncm_vector_new_array (x_a);

  ncm_mset_fparams_get_vector (self->mset, x);
  H_a = ncm_diff_rf_Hessian_N_to_1 (self->diff, x_a, _ncm_fit_numdiff_m2lnL_val, fit, NULL);

  ncm_matrix_set_from_array (H, H_a);
  ncm_mset_fparams_set_vector (self->mset, x);

  g_array_unref (x_a);
  g_array_unref (H_a);
  ncm_vector_free (x);
}

/**
 * ncm_fit_fisher_to_covar:
 * @fit: a #NcmFit
 * @fisher: a #NcmMatrix
 *
 * Inverts the matrix @fisher and sets as the covariance matrix
 * of @fit. The Fisher matrix used can be both the Fisher or the
 * observed Fisher matrices.
 *
 */
void
ncm_fit_fisher_to_covar (NcmFit *fit, NcmMatrix *fisher)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmMatrix *covar    = ncm_fit_state_peek_covar (self->fstate);
  gint ret;

  if (ncm_mset_fparam_len (self->mset) == 0)
    g_error ("ncm_fit_fisher_to_covar: mset object has 0 free parameters");

  ncm_matrix_memcpy (covar, fisher);

  ret = ncm_matrix_cholesky_decomp (covar, 'U');

  if (ret == 0)
  {
    ret = ncm_matrix_cholesky_inverse (covar, 'U');

    if (ret != 0)
      g_error ("ncm_fit_fisher_to_covar[ncm_matrix_cholesky_decomp]: %d.", ret);

    ncm_matrix_copy_triangle (covar, 'U');
  }
  else if (ret > 0)
  {
    NcmMatrix *LU      = ncm_matrix_dup (fisher);
    gsl_permutation *p = gsl_permutation_alloc (ncm_mset_fparam_len (self->mset));
    gint signum;
    gint ret1;

    ncm_matrix_scale (LU, 0.5);

    g_warning ("ncm_fit_fisher_to_covar: covariance matrix not positive definite, errors are not trustworthy.");

    ret1 = gsl_linalg_LU_decomp (ncm_matrix_gsl (LU), p, &signum);
    NCM_TEST_GSL_RESULT ("ncm_fit_fisher_to_covar[gsl_linalg_LU_decomp]", ret1);

    ret1 = gsl_linalg_LU_invert (ncm_matrix_gsl (LU), p, ncm_matrix_gsl (covar));
    NCM_TEST_GSL_RESULT ("ncm_fit_fisher_to_covar[gsl_linalg_LU_invert]", ret1);

    gsl_permutation_free (p);
    ncm_matrix_free (LU);
  }
  else
  {
    g_error ("ncm_fit_fisher_to_covar[ncm_matrix_cholesky_decomp]: %d.", ret);
  }

  ncm_fit_state_set_has_covar (self->fstate, TRUE);
}

/**
 * ncm_fit_obs_fisher:
 * @fit: a #NcmFit
 *
 * Calculates the covariance from the observed Fisher
 * matrix, see ncm_fit_numdiff_m2lnL_covar().
 *
 */
void
ncm_fit_obs_fisher (NcmFit *fit)
{
  ncm_fit_numdiff_m2lnL_covar (fit);
}

/**
 * ncm_fit_fisher:
 * @fit: a #NcmFit
 *
 * Calculates the covariance from the Fisher matrix, see
 * ncm_dataset_fisher_matrix().
 *
 */
void
ncm_fit_fisher (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmMatrix *IM       = NULL;

  if ((ncm_likelihood_priors_length_f (self->lh) > 0) || (ncm_likelihood_priors_length_m2lnL (self->lh) > 0))
    g_warning ("ncm_fit_fisher: the analysis contains priors which are ignored in the Fisher matrix calculation.");

  ncm_dataset_fisher_matrix (self->lh->dset, self->mset, &IM);
  ncm_fit_fisher_to_covar (fit, IM);

  ncm_matrix_clear (&IM);
}

/**
 * ncm_fit_numdiff_m2lnL_covar:
 * @fit: a #NcmFit
 *
 * Calcualtes the covariance matrix using the inverse of the
 * Hessian matrix $\partial_i\partial_j -\ln(L)$, where
 * the derivatives are taken with respect to the free parameters.
 *
 */
void
ncm_fit_numdiff_m2lnL_covar (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmMatrix *hessian  = ncm_fit_state_peek_hessian (self->fstate);

  if (ncm_mset_fparam_len (self->mset) == 0)
    g_error ("ncm_fit_numdiff_m2lnL_covar: mset object has 0 free parameters");

  ncm_fit_numdiff_m2lnL_hessian (fit, hessian, self->params_reltol);
  ncm_matrix_scale (hessian, 0.5);

  ncm_fit_fisher_to_covar (fit, hessian);
}

/**
 * ncm_fit_numdiff_m2lnL_lndet_covar:
 * @fit: a #NcmFit
 *
 * Calculates the logarithm of the determinant of the covariance matrix
 * using the inverse of the Hessian matrix $\partial_i\partial_j -\ln(L)$,
 * where the derivatives are taken with respect to the free parameters.
 *
 * Returns: the logarithm of the determinant of the covariance matrix.
 */
gdouble
ncm_fit_numdiff_m2lnL_lndet_covar (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmMatrix *covar    = ncm_fit_state_peek_covar (self->fstate);
  NcmMatrix *hessian  = ncm_fit_state_peek_hessian (self->fstate);
  const guint len     = ncm_matrix_nrows (covar);
  gdouble lndetC      = 0.0;
  guint i;
  gint ret;

  if (ncm_mset_fparam_len (self->mset) == 0)
    g_error ("ncm_fit_numdiff_m2lnL_covar: mset object has 0 free parameters");

  ncm_fit_numdiff_m2lnL_hessian (fit, hessian, 1.0e-2);
  ncm_matrix_memcpy (covar, hessian);
  ncm_matrix_scale (covar, 0.5);

  ret = ncm_matrix_cholesky_decomp (covar, 'U');

  if (ret == 0)
  {
    for (i = 0; i < len; i++)
      lndetC += log (ncm_matrix_get (covar, i, i));

    lndetC = -2.0 * lndetC;
  }
  else if (ret > 0)
  {
    NcmMatrix *LU      = ncm_matrix_dup (hessian);
    gsl_permutation *p = gsl_permutation_alloc (ncm_mset_fparam_len (self->mset));
    gint signum;
    gint ret1;

    ncm_matrix_scale (LU, 0.5);

    g_warning ("ncm_fit_numdiff_m2lnL_covar: covariance matrix not positive definite, errors are not trustworthy.");

    ret1 = gsl_linalg_LU_decomp (ncm_matrix_gsl (LU), p, &signum);
    NCM_TEST_GSL_RESULT ("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_LU_decomp]", ret1);

    for (i = 0; i < len; i++)
      lndetC += log (fabs (ncm_matrix_get (LU, i, i)));

    lndetC = -lndetC;

    gsl_permutation_free (p);
    ncm_matrix_free (LU);
  }
  else
  {
    g_error ("ncm_fit_numdiff_m2lnL_lndet_covar[ncm_matrix_cholesky_decomp]: %d.", ret);
  }

  return lndetC;
}

/**
 * ncm_fit_get_covar:
 * @fit: a #NcmFit
 *
 * Returns a copy of the covariance matrix (pre-calculated by, e.g, ncm_fit_numdiff_m2lnL_covar()).
 *
 * Returns: (transfer full): the covariance matrix
 */
NcmMatrix *
ncm_fit_get_covar (NcmFit *fit)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  g_assert (ncm_fit_state_has_covar (self->fstate));

  return ncm_matrix_dup (ncm_fit_state_peek_covar (self->fstate));
}

typedef struct _FitDProb
{
  NcmMSetPIndex pi;
  NcmFit *fit_val;
  NcmFit *fit;
} FitDProb;

static gdouble
fit_dprob (gdouble val, gpointer p)
{
  FitDProb *dprob_arg     = (FitDProb *) p;
  NcmFitPrivate *self     = ncm_fit_get_instance_private (dprob_arg->fit);
  NcmFitPrivate *self_val = ncm_fit_get_instance_private (dprob_arg->fit_val);

  ncm_mset_param_set (self_val->mset, dprob_arg->pi.mid, dprob_arg->pi.pid, val);
  ncm_fit_run (dprob_arg->fit_val, NCM_FIT_RUN_MSGS_NONE);

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message (".");

  {
    const gdouble m2lnL_fv = ncm_fit_state_get_m2lnL_curval (self_val->fstate);
    const gdouble m2lnL    = ncm_fit_state_get_m2lnL_curval (self->fstate);

    return exp (-(m2lnL_fv - m2lnL) / 2.0);
  }
}

/**
 * ncm_fit_prob:
 * @fit: a #NcmFit
 * @mid: a #NcmModelID
 * @pid: a parameter id
 * @a: probability lower limit
 * @b: probability upper limit
 *
 * Computes the probability of the parameter @pid of the model @mid
 * to be in the interval [@a, @b].
 *
 * Returns: the probability
 */
gdouble
ncm_fit_prob (NcmFit *fit, NcmModelID mid, guint pid, gdouble a, gdouble b)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmSerialize *ser   = ncm_serialize_global ();
  NcmMSet *mset_val   = ncm_mset_dup (self->mset, ser);
  gsl_integration_workspace **w;
  NcmFit *fit_val;
  FitDProb dprob_arg;
  gdouble result, error;
  gint error_code;
  gsl_function F;

  ncm_serialize_free (ser);
  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, self->lh, mset_val, self->grad.gtype);

  dprob_arg.pi.mid  = mid;
  dprob_arg.pi.pid  = pid;
  dprob_arg.fit     = fit;
  dprob_arg.fit_val = fit_val;

  F.function = &fit_dprob;
  F.params   = &dprob_arg;

  w = ncm_integral_get_workspace ();

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("#");

  error_code = gsl_integration_qags (&F, a, b, 1e-10, NCM_INTEGRAL_ERROR, NCM_INTEGRAL_PARTITION, *w, &result, &error);

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    g_message ("\n");

  NCM_TEST_GSL_RESULT ("ncm_fit_prob", error_code);
  ncm_memory_pool_return (w);

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);

  return result;
}

/**
 * ncm_fit_dprob:
 * @fit: a #NcmFit
 * @mid: a #NcmModelID
 * @pid: a parameter id
 * @a: probability lower limit
 * @b: probability upper limit
 * @step: probability density step
 * @norm: normalization factor
 *
 * Computes the probability density of the parameter @pid of the model @mid
 * in the interval [@a, @b] with a step @step.
 *
 */
void
ncm_fit_dprob (NcmFit *fit, NcmModelID mid, guint pid, gdouble a, gdouble b, gdouble step, gdouble norm)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmSerialize *ser   = ncm_serialize_global ();
  NcmMSet *mset_val   = ncm_mset_dup (self->mset, ser);
  NcmFit *fit_val;
  FitDProb dprob_arg;
  gdouble point;

  ncm_serialize_free (ser);
  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, self->lh, mset_val, self->grad.gtype);

  dprob_arg.pi.mid  = mid;
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
 * @pid: the parameter id
 * @start: starting value
 * @stop: ending value
 * @step: step size
 *
 * Computes the likelihood ratio test for the parameter @pid of the model @mid
 * in the interval [@start, @stop] with a step @step. The result is printed
 * to the log.
 *
 */
void
ncm_fit_lr_test_range (NcmFit *fit, NcmModelID mid, guint pid, gdouble start, gdouble stop, gdouble step)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmSerialize *ser   = ncm_serialize_global ();
  NcmMSet *mset_val   = ncm_mset_dup (self->mset, ser);
  NcmFit *fit_val;
  gdouble walk;

  ncm_serialize_free (ser);
  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, self->lh, mset_val, self->grad.gtype);

  {
    NcmFitPrivate *self_val = ncm_fit_get_instance_private (fit_val);

    for (walk = start; walk <= stop; walk += step)
    {
      ncm_mset_param_set (self_val->mset, mid, pid, walk);
      ncm_fit_run (fit_val, NCM_FIT_RUN_MSGS_NONE);

      {
        const gdouble m2lnL_fv = ncm_fit_state_get_m2lnL_curval (self_val->fstate);
        const gdouble m2lnL    = ncm_fit_state_get_m2lnL_curval (self->fstate);

        g_message ("%g %g %g %g %g\n", walk,
                   (m2lnL_fv - m2lnL),
                   gsl_ran_chisq_pdf (m2lnL_fv - m2lnL, 1),
                   gsl_cdf_chisq_Q (m2lnL_fv - m2lnL, 1),
                   gsl_cdf_ugaussian_Q (sqrt (m2lnL_fv - m2lnL))
                  );
      }
    }

    ncm_fit_free (fit_val);
    ncm_mset_free (mset_val);
  }

  return;
}

/**
 * ncm_fit_lr_test:
 * @fit: a #NcmFit.
 * @mid: a #NcmModelID.
 * @pid: the parameter id
 * @val: the parameter value
 * @dof: degrees of freedom
 *
 * Computes the likelihood ratio test for the parameter @pid of the model @mid
 * with the value @val.
 *
 * Returns: the probability of the null hypothesis.
 */
gdouble
ncm_fit_lr_test (NcmFit *fit, NcmModelID mid, guint pid, gdouble val, gint dof)
{
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);
  NcmSerialize *ser   = ncm_serialize_global ();
  NcmMSet *mset_val   = ncm_mset_dup (self->mset, ser);
  NcmFit *fit_val;
  gdouble result;

  ncm_serialize_free (ser);
  ncm_mset_param_set_ftype (mset_val, mid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_copy_new (fit, self->lh, mset_val, self->grad.gtype);

  {
    NcmFitPrivate *self_val = ncm_fit_get_instance_private (fit_val);

    ncm_mset_param_set (self_val->mset, mid, pid, val);
    ncm_fit_run (fit_val, NCM_FIT_RUN_MSGS_NONE);
    {
      const gdouble m2lnL_fv = ncm_fit_state_get_m2lnL_curval (self_val->fstate);
      const gdouble m2lnL    = ncm_fit_state_get_m2lnL_curval (self->fstate);

      result = gsl_cdf_chisq_Q (m2lnL_fv - m2lnL, dof);
    }
    ncm_fit_free (fit_val);
    ncm_mset_free (mset_val);
  }

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
 *  gdouble
 *  ncm_fit_chisq_test (NcmFit *fit, size_t bins)
 *  {
 *  NcmLikelihood *lh = self->lh;
 *  NcmDataset *ds = lh->ds;
 *  NcmData *data = nc_dataset_get_data (ds, 0);
 *  NcFunction distance = nc_function_get (data->res_type);
 *  gsl_histogram *obs = gsl_histogram_alloc (bins);
 *  gsl_histogram *exp = gsl_histogram_alloc (bins);
 *  gsl_vector_view y_exp;
 *  gdouble y_min, y_max;
 *  gint i;
 *  if (data == NULL)
 *  return GSL_NAN;
 *
 *  ncm_distance_thread_pool (distance.f, self->cp, data->x, self->f, NULL, NULL);
 *  y_exp = gsl_vector_subvector (self->f, 0, data->x->size);
 *  y_min = GSL_MIN (gsl_vector_min (data->y), gsl_vector_min (&y_exp.vector))*0.999;
 *  y_max = GSL_MAX (gsl_vector_max (data->y), gsl_vector_max (&y_exp.vector))*1.001;
 *
 *  gsl_histogram_set_ranges_uniform (obs, y_min, y_max);
 *  gsl_histogram_set_ranges_uniform (exp, y_min, y_max);
 *
 *  for (i = 0; i < data->x->size; i++)
 *  {
 *    gsl_histogram_increment (obs, gsl_vector_get (data->y, i));
 *    gsl_histogram_increment (exp, gsl_vector_get (&y_exp.vector, i));
 *    }
 *
 *    gsl_histogram_sub (obs, exp);
 *    //  if (gsl_histogram_min_val (exp) == 0.0)
 *    //    g_error ("Cannot chisq test");
 *    gsl_histogram_mul (obs, obs);
 *    gsl_histogram_div (obs, exp);
 *    //gsl_histogram_fprintf (self->log, obs, "%f", "%f");
 *    g_message ("\n");
 *    //gsl_histogram_fprintf (self->log, exp, "%f", "%f");
 *    g_message ("\n");
 *    g_message ("BLA [%g, %g] %zu %g [%g]\n", y_min, y_max, bins, gsl_histogram_sum (obs), gsl_histogram_sum (exp));
 *    return gsl_histogram_sum (obs);
 *    }
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
  NcmFitPrivate *self = ncm_fit_get_instance_private (fit);

  if (ncm_mset_fparam_len (self->mset) < 1)
  {
    g_warning ("ncm_fit_function_error: called but no free parameters were set in #NcmMSet.");
    *f       = ncm_mset_func_eval_nvar (func, self->mset, x);
    *sigma_f = 0.0;
  }
  else
  {
    NcmMatrix *covar = ncm_fit_state_peek_covar (self->fstate);
    NcmVector *v     = NULL;
    NcmVector *tmp1  = NULL;
    gdouble result;
    gint ret;

    ncm_mset_func_numdiff_fparams (func, self->mset, x, &v);
    tmp1 = ncm_vector_dup (v);

    if (!ncm_fit_state_has_covar (self->fstate))
      g_error ("ncm_fit_function_error: called without any covariance matrix calculated.");

    ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (covar), ncm_vector_gsl (v), 0.0, ncm_vector_gsl (tmp1));
    NCM_TEST_GSL_RESULT ("ncm_fit_function_error[covar.v]", ret);
    ret = gsl_blas_ddot (ncm_vector_gsl (v), ncm_vector_gsl (tmp1), &result);
    NCM_TEST_GSL_RESULT ("ncm_fit_function_error[v.covar.v]", ret);

    *f       = ncm_mset_func_eval_nvar (func, self->mset, x);
    *sigma_f = sqrt (result);

    if (pretty_print)
      g_message ("# % -12.4g +/- % -12.4g\n", *f, *sigma_f);

    ncm_vector_free (v);
    ncm_vector_free (tmp1);

    return;
  }
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
  NCM_UNUSED (fit);
  NCM_UNUSED (func1);
  NCM_UNUSED (z1);
  NCM_UNUSED (func2);
  NCM_UNUSED (z2);
  NCM_UNUSED (pretty_print);
  g_assert_not_reached ();

  /*
   *  gdouble result, cor, s1, s2;
   *  gint ret;
   *  NcmVector *tmp1 = ncm_fit_params_get_tmp_vector (self->pt, self->mset);
   *  NcmVector *tmp2 = ncm_fit_params_get_tmp_vector (self->pt, self->mset);
   *
   *  NCM_FUNC_DF (func1, self->mset, self->pt, z1, tmp1);
   *  NCM_FUNC_DF (func2, self->mset, self->pt, z2, tmp2);
   *
   *  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (self->fstate->covar), ncm_vector_gsl (tmp1), 0.0, ncm_vector_gsl (self->df));
   *  NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
   *  ret = gsl_blas_ddot (ncm_vector_gsl (self->df), ncm_vector_gsl (tmp1), &result);
   *  NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
   *  s1 = sqrt(result);
   *
   *  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (self->fstate->covar), ncm_vector_gsl (tmp2), 0.0, ncm_vector_gsl (self->df));
   *  NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
   *  ret = gsl_blas_ddot (ncm_vector_gsl (self->df), ncm_vector_gsl (tmp2), &result);
   *  NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
   *  s2 = sqrt(result);
   *
   *  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (self->fstate->covar), ncm_vector_gsl (tmp1), 0.0, ncm_vector_gsl (self->df));
   *  NCM_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
   *  ret = gsl_blas_ddot (ncm_vector_gsl (self->df), ncm_vector_gsl (tmp2), &result);
   *  NCM_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
   *  cor = result / (s1*s2);
   *
   *  if (pretty_print)
   *  {
   *  g_message ("# % -12.4g\t| % -12.4g\n", NCM_FUNC_F(func1, self->mset, z1), NCM_FUNC_F(func2, self->mset, z2));
   *  g_message ("#---------------------------------------------\n");
   *  g_message ("# % -12.4g\t | % -12.4g\t| % -12.4g\n", s1, 1.0, cor);
   *  g_message ("# % -12.4g\t | % -12.4g\t| % -12.4g\n", s2, cor, 1.0);
   *  g_message ("#---------------------------------------------\n");
   *  }
   *
   *  ncm_fit_params_return_tmp_vector (self->pt, tmp1);
   *  ncm_fit_params_return_tmp_vector (self->pt, tmp2);
   *
   *  return sqrt(result);
   */
}

