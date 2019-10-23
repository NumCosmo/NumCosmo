/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_csq1d.c
 *
 *  Mon September 09 13:56:19 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_csq1d.c
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_csq1d
 * @title: NcmCSQ1D
 * @short_description: Abstract class for Harmonic Oscillator calculation through complex structure quantization.
 *
 * 
 * The system:
 * \begin{align}
 * q^\prime &= \frac{\Pi_q}{m},
 * \Pi_q^\prime &= m\nu^2q.
 * \end{align}
 * 
 * \begin{equation}
 * \xi = \ln (m\nu)
 * \end{equation}
 * 
 * \begin{equation}
 * F^n = \left(\frac{1}{2\nu}\frac{\partial}{\partial t}\right)^n \xi. 
 * \end{equation}
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include <complex.h>

#include "math/ncm_csq1d.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_diff.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_multifit.h>
#include <nlopt.h>

#include <nvector/nvector_serial.h>

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <sundials/sundials_types.h> 
#endif /* NUMCOSMO_GIR_SCAN */

#define PRINT_EVOL TRUE

typedef enum _NcmCSQ1DEvolStop
{
  NCM_CSQ1D_EVOL_STOP_ERROR = -1,
  NCM_CSQ1D_EVOL_STOP_FINISHED = 0,
  NCM_CSQ1D_EVOL_STOP_ADIABATIC_START,
  NCM_CSQ1D_EVOL_STOP_UP_START,
  NCM_CSQ1D_EVOL_STOP_UM_START,
} NcmCSQ1DEvolStop;

typedef enum _NcmCSQ1DEvolState
{
  NCM_CSQ1D_EVOL_STATE_INVALID = 0,
  NCM_CSQ1D_EVOL_STATE_ADIABATIC,
  NCM_CSQ1D_EVOL_STATE_UP,
  NCM_CSQ1D_EVOL_STATE_UM,
} NcmCSQ1DEvolState;

struct _NcmCSQ1DPrivate
{
  gdouble reltol;
  gdouble abstol;
  gdouble k;
  gdouble ti;
  gdouble tf;
  gdouble t;
  NcmCSQ1DEvolState state;
  gdouble adiab_threshold;
  gboolean save_evol;
  gboolean sing_detect;
  NcmModelCtrl *ctrl;
  gpointer cvode;
  gpointer cvode_Up;
  gpointer cvode_Um;
  gpointer arkode;
  gboolean cvode_init;
  gboolean cvode_Up_init;
  gboolean cvode_Um_init;
  gboolean arkode_init;
  N_Vector y;
  N_Vector y_Up;
  N_Vector y_Um;
  SUNMatrix A;
  SUNMatrix A_Up;
  SUNMatrix A_Um;
  SUNLinearSolver LS;
  SUNLinearSolver LS_Up;
  SUNLinearSolver LS_Um;
  NcmSpline *alpha_s;
  NcmSpline *dgamma_s;
  NcmDiff *diff;
};

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_K,
  PROP_TI,
  PROP_TF,
  PROP_ADIAB_THRESHOLD,
  PROP_SAVE_EVOL,
  PROP_SING_DETECT,
};

G_DEFINE_BOXED_TYPE (NcmCSQ1DSingFitUp, ncm_csq1d_sing_fit_up, ncm_csq1d_sing_fit_up_dup, ncm_csq1d_sing_fit_up_free);
G_DEFINE_TYPE_WITH_PRIVATE (NcmCSQ1D, ncm_csq1d, G_TYPE_OBJECT);

static void
ncm_csq1d_init (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv = ncm_csq1d_get_instance_private (csq1d);

  self->reltol          = 0.0;
  self->abstol          = 0.0;
  self->k               = 0.0;
  self->ti              = 0.0;
  self->tf              = 0.0;
  self->state           = NCM_CSQ1D_EVOL_STATE_INVALID;
  self->adiab_threshold = 0.0;
  self->save_evol       = FALSE;
  self->ctrl            = ncm_model_ctrl_new (NULL);

  self->cvode           = NULL;
  self->cvode_init      = FALSE;
  self->cvode_Up        = NULL;
  self->cvode_Up_init   = FALSE;
  self->cvode_Um        = NULL;
  self->cvode_Um_init   = FALSE;
  self->arkode          = NULL;
  self->arkode_init     = FALSE;

  self->y               = N_VNew_Serial (2);
  self->y_Up            = N_VNew_Serial (2);
  self->y_Um            = N_VNew_Serial (2);
  
  self->A               = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );  

  self->A_Up           = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer)self->A_Up, "SUNDenseMatrix", 0, );  

  self->A_Um           = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer)self->A_Um, "SUNDenseMatrix", 0, );  

  self->LS              = SUNDenseLinearSolver (self->y, self->A);
  NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );

  self->LS_Up          = SUNDenseLinearSolver (self->y_Up, self->A_Up);
  NCM_CVODE_CHECK ((gpointer)self->LS_Up, "SUNDenseLinearSolver", 0, );

  self->LS_Um          = SUNDenseLinearSolver (self->y_Um, self->A_Um);
  NCM_CVODE_CHECK ((gpointer)self->LS_Um, "SUNDenseLinearSolver", 0, );

  self->alpha_s         = ncm_spline_cubic_notaknot_new ();
  self->dgamma_s        = ncm_spline_cubic_notaknot_new ();

  self->diff            = ncm_diff_new ();
}

static void
_ncm_csq1d_dispose (GObject *object)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = csq1d->priv;
  
  ncm_model_ctrl_clear (&self->ctrl);
  ncm_spline_clear (&self->alpha_s);
  ncm_spline_clear (&self->dgamma_s);

  ncm_diff_clear (&self->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_csq1d_parent_class)->dispose (object);
}

static void
_ncm_csq1d_finalize (GObject *object)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = csq1d->priv;

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode      = NULL;
    self->cvode_init = FALSE;
  }
  if (self->cvode_Up != NULL)
  {
    CVodeFree (&self->cvode_Up);
    self->cvode_Up      = NULL;
    self->cvode_Up_init = FALSE;
  }
  if (self->cvode_Um != NULL)
  {
    CVodeFree (&self->cvode_Um);
    self->cvode_Um      = NULL;
    self->cvode_Um_init = FALSE;
  }
  if (self->arkode != NULL)
  {
    ARKStepFree (&self->arkode);
    self->arkode      = NULL;
    self->arkode_init = FALSE;
  }

  g_clear_pointer (&self->y,    N_VDestroy);
  g_clear_pointer (&self->y_Up, N_VDestroy);
  g_clear_pointer (&self->y_Um, N_VDestroy);

  if (self->A != NULL)
    SUNMatDestroy (self->A);

  if (self->A_Up != NULL)
    SUNMatDestroy (self->A_Up);

  if (self->A_Um != NULL)
    SUNMatDestroy (self->A_Um);

  if (self->LS != NULL)
  {
    gint flag = SUNLinSolFree (self->LS);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  if (self->LS_Up != NULL)
  {
    gint flag = SUNLinSolFree (self->LS_Up);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  if (self->LS_Um != NULL)
  {
    gint flag = SUNLinSolFree (self->LS_Um);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_csq1d_parent_class)->finalize (object);
}

static void
_ncm_csq1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  g_return_if_fail (NCM_IS_CSQ1D (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      ncm_csq1d_set_reltol (csq1d, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_csq1d_set_abstol (csq1d, g_value_get_double (value));
      break;
    case PROP_K:
      ncm_csq1d_set_k (csq1d, g_value_get_double (value));
      break;
    case PROP_TI:
      ncm_csq1d_set_ti (csq1d, g_value_get_double (value));
      break;
    case PROP_TF:
      ncm_csq1d_set_tf (csq1d, g_value_get_double (value));
      break;
    case PROP_ADIAB_THRESHOLD:
      ncm_csq1d_set_adiab_threshold (csq1d, g_value_get_double (value));
      break;
    case PROP_SAVE_EVOL:
      ncm_csq1d_set_save_evol (csq1d, g_value_get_boolean (value));
      break;
    case PROP_SING_DETECT:
      ncm_csq1d_set_sing_detect (csq1d, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_csq1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  g_return_if_fail (NCM_IS_CSQ1D (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, ncm_csq1d_get_reltol (csq1d));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_csq1d_get_abstol (csq1d));
      break;
    case PROP_K:
      g_value_set_double (value, ncm_csq1d_get_k (csq1d));
      break;
    case PROP_TI:
      g_value_set_double (value, ncm_csq1d_get_ti (csq1d));
      break;
    case PROP_TF:
      g_value_set_double (value, ncm_csq1d_get_tf (csq1d));
      break;
    case PROP_ADIAB_THRESHOLD:
      g_value_set_double (value, ncm_csq1d_get_adiab_threshold (csq1d));
      break;
    case PROP_SAVE_EVOL:
      g_value_set_boolean (value, ncm_csq1d_get_save_evol (csq1d));
      break;
    case PROP_SING_DETECT:
      g_value_set_boolean (value, ncm_csq1d_get_sing_detect (csq1d));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble _ncm_csq1d_eval_xi  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_dxi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_nu  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_m   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_dm  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_F1  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_F2  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_FN  (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);

static void
ncm_csq1d_class_init (NcmCSQ1DClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_csq1d_set_property;
  object_class->get_property = &_ncm_csq1d_get_property;
  object_class->dispose      = &_ncm_csq1d_dispose;
  object_class->finalize     = &_ncm_csq1d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NCM_DEFAULT_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance tolerance",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_K,
                                   g_param_spec_double ("k", 
                                                        NULL, 
                                                        "Mode k",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TI,
                                   g_param_spec_double ("ti",
                                                        NULL,
                                                        "The initial time t_i",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TF,
                                   g_param_spec_double ("tf",
                                                        NULL,
                                                        "The final time t_f",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ADIAB_THRESHOLD,
                                   g_param_spec_double ("adiab-threshold",
                                                        NULL,
                                                        "The adiabatic threshold",
                                                        0.0, G_MAXDOUBLE, 1.0e-1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SAVE_EVOL,
                                   g_param_spec_boolean ("save-evol",
                                                         NULL,
                                                         "Save the system evolution",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SING_DETECT,
                                   g_param_spec_boolean ("sing-detect",
                                                         NULL,
                                                         "Singularity detection",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->eval_xi             = &_ncm_csq1d_eval_xi;
  klass->eval_dxi            = &_ncm_csq1d_eval_dxi;
  klass->eval_nu             = &_ncm_csq1d_eval_nu;
  klass->eval_nu2            = &_ncm_csq1d_eval_nu2;
  klass->eval_m              = &_ncm_csq1d_eval_m;
  klass->eval_dm             = &_ncm_csq1d_eval_dm;
  klass->eval_F1             = &_ncm_csq1d_eval_F1;
  klass->eval_F2             = &_ncm_csq1d_eval_F2;
  klass->eval_FN             = &_ncm_csq1d_eval_FN;
  klass->eval_powspec_factor = &_ncm_csq1d_eval_powspec_factor;
}

static gdouble 
_ncm_csq1d_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_xi: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_dxi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_F1 (csq1d, model, t, k) * 2.0 * NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t, k);
}

static gdouble 
_ncm_csq1d_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_nu: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return gsl_pow_2 (NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t, k));
}

static gdouble 
_ncm_csq1d_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return exp (NCM_CSQ1D_GET_CLASS (csq1d)->eval_xi (csq1d, model, t, k)) / NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t, k);
}

static gdouble 
_ncm_csq1d_eval_dm (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_dm: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_F1: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_F2: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_FN (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_FN: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_powspec_factor: not implemented."); 
  return 0.0;
}

/***********************************************************************************
 * Singular Up model
 ***********************************************************************************/
/**
 * ncm_csq1d_sing_fit_up_new:
 * @chi_dim: model dimension for $\chi$ fitting
 * @Up_dim: model dimension for $\chi$ fitting
 * 
 * New singular model for fitting $\chi$ and $\Upsilon_+$
 * across a singularity.
 * 
 * Returns: (transfer full): a new NcmCSQ1DSingFit. 
 */
NcmCSQ1DSingFitUp *
ncm_csq1d_sing_fit_up_new (const gint chi_dim, const gint Up_dim)
{
  NcmCSQ1DSingFitUp *sing_up = g_new0 (NcmCSQ1DSingFitUp, 1);

  g_assert_cmpint (chi_dim,     >=, 2);
  g_assert_cmpint (chi_dim % 2, ==, 0);

  g_assert_cmpint (Up_dim,      >=, 1);
  g_assert_cmpint (Up_dim  % 2, ==, 1);

  sing_up->chi_dim = chi_dim;
  sing_up->Up_dim  = Up_dim;
  sing_up->chi_c   = ncm_vector_new (chi_dim);
  sing_up->Up_c    = ncm_vector_new (Up_dim);

  return sing_up;
}

/**
 * ncm_csq1d_sing_fit_up_dup:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * 
 * Duplicates @sing_up, without duplicating the contents. 
 * 
 * Returns: (transfer full): a shallow copy of @sing_up. 
 */
NcmCSQ1DSingFitUp *
ncm_csq1d_sing_fit_up_dup (NcmCSQ1DSingFitUp *sing_up)
{
  NcmCSQ1DSingFitUp *sing_up_dup = g_new0 (NcmCSQ1DSingFitUp, 0);

  sing_up_dup->chi_dim = sing_up->chi_dim;
  sing_up_dup->Up_dim  = sing_up->Up_dim;
  sing_up_dup->chi_c   = ncm_vector_ref (sing_up->chi_c);
  sing_up_dup->Up_c    = ncm_vector_ref (sing_up->Up_c);

  return sing_up_dup;
}

/**
 * ncm_csq1d_sing_fit_up_free:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * 
 * Frees @sing_up. 
 */
void
ncm_csq1d_sing_fit_up_free (NcmCSQ1DSingFitUp *sing_up)
{
  ncm_vector_clear (&sing_up->chi_c);
  ncm_vector_clear (&sing_up->Up_c);

  g_free (sing_up);
}

/**
 * ncm_csq1d_sing_fit_up_fit:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time vector
 * @chim_t: $m\chi/t$ vector
 * @exp_Up: $\exp(\Upsilon_+)$ vector
 * 
 * Fits the models using the data in @t, @chim and @exp_Up.
 *  
 */
void
ncm_csq1d_sing_fit_up_fit (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, NcmVector *t, NcmVector *chim_t, NcmVector *exp_Up)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;

  g_assert_cmpint (ncm_vector_len (t), ==, ncm_vector_len (chim_t));
  g_assert_cmpint (ncm_vector_len (t), ==, ncm_vector_len (exp_Up));
  
  {
    const gint len                      = ncm_vector_len (chim_t);
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (len, sing_up->chi_dim);
    NcmMatrix *cov                      = ncm_matrix_new (sing_up->chi_dim, sing_up->chi_dim);
    NcmMatrix *X                        = ncm_matrix_new (len, sing_up->chi_dim);
    gdouble chisq;
    gint ret;
    gint i;

    printf ("Fitting chi * m / t.\n");
    for (i = 0; i < len; i++)
    {
      const gdouble t_i = ncm_vector_get (t, i);
      const gdouble m_i = ncm_csq1d_eval_m (csq1d, model, t_i, self->k);
      gdouble Ta        = 1.0;
      gdouble Tb        = m_i / t_i;
      gdouble lfac      = 1.0;
      gint j;

      for (j = 0; j < sing_up->chi_dim; j += 2)
      {
        ncm_matrix_set (X, i, j + 0, Ta);
        ncm_matrix_set (X, i, j + 1, Tb);

        Ta = Ta * t_i / lfac;
        Tb = Tb * t_i / lfac;
        lfac += 1.0;
      }
    }

    ret = gsl_multifit_linear (ncm_matrix_gsl (X), ncm_vector_gsl (chim_t), ncm_vector_gsl (sing_up->chi_c), ncm_matrix_gsl (cov), &chisq, work);
    g_assert_cmpint (ret, ==, 0);

    printf ("MULTIFIT RET: %d, chisq % 22.15g, sqrt (chisq / len) % 22.15g\n", ret, chisq, sqrt (chisq / len));

    ncm_vector_log_vals (sing_up->chi_c, "COEFS: ", "% 22.15g", TRUE);
    ncm_matrix_log_vals (cov,            "COV:  ", "% 22.15g");

    ncm_matrix_free (cov);
    ncm_matrix_free (X);
    gsl_multifit_linear_free (work);
  }

  {
    const gint len                      = ncm_vector_len (exp_Up);
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (len, sing_up->Up_dim);
    NcmMatrix *cov                      = ncm_matrix_new (sing_up->Up_dim, sing_up->Up_dim);
    NcmMatrix *X                        = ncm_matrix_new (len, sing_up->Up_dim);
    gdouble chisq;
    gint ret, i;

    for (i = 0; i < len; i++)
    {
      const gdouble t_i = ncm_vector_get (t, i);
      const gdouble m_i = ncm_csq1d_eval_m (csq1d, model, t_i, self->k);
      gdouble Ta        = t_i * t_i;
      gdouble Tb        = m_i * t_i;
      gdouble lfac      = 1.0;
      gint j;

      ncm_matrix_set (X, i, 0, 1.0);

      for (j = 1; j < sing_up->Up_dim; j += 2)
      {
        ncm_matrix_set (X, i, j + 0, Ta);
        ncm_matrix_set (X, i, j + 1, Tb);

        Ta = Ta * t_i / lfac;
        Tb = Tb * t_i / lfac;
        lfac += 1.0;
      }
    }

    ret = gsl_multifit_linear (ncm_matrix_gsl (X), ncm_vector_gsl (exp_Up), ncm_vector_gsl (sing_up->Up_c), ncm_matrix_gsl (cov), &chisq, work);
    g_assert_cmpint (ret, ==, 0);

    printf ("MULTIFIT RET: %d, chisq % 22.15g, sqrt (chisq / len) % 22.15g\n", ret, chisq, sqrt (chisq / len));

    ncm_vector_log_vals (sing_up->Up_c, "COEFS: ", "% 22.15g", TRUE);
    ncm_matrix_log_vals (cov,           "COV:  ",  "% 22.15g");

    ncm_matrix_free (cov);
    ncm_matrix_free (X);
    gsl_multifit_linear_free (work);
  }  
}

/**
 * ncm_csq1d_sing_fit_up_eval_chi:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * 
 * Evaluates the model for $\chi$ at $t$.
 *  
 */
gdouble
ncm_csq1d_sing_fit_up_eval_chi (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble m = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  gdouble Ta      = t / m;
  gdouble Tb      = 1.0;
  gdouble lfac    = 1.0;
  gdouble res     = 0.0;
  gint j;

  for (j = 0; j < sing_up->chi_dim; j += 2)
  {
    res += Ta * ncm_vector_get (sing_up->chi_c, j + 0);
    res += Tb * ncm_vector_get (sing_up->chi_c, j + 1);
    Ta = Ta * t / lfac;
    Tb = Tb * t / lfac;
    lfac += 1.0;
  }

  return res;
}

/**
 * ncm_csq1d_sing_fit_up_eval_exp_Up:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * 
 * Evaluates the model for $\exp(\Upsilon_+)$ at $t$.
 *  
 */
gdouble
ncm_csq1d_sing_fit_up_eval_exp_Up (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble m = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  gdouble Ta      = t * t;
  gdouble Tb      = m * t;
  gdouble lfac    = 1.0;
  gdouble res     = ncm_vector_get (sing_up->Up_c, 0);
  gint j;

  for (j = 1; j < sing_up->Up_dim; j += 2)
  {
    res += Ta * ncm_vector_get (sing_up->Up_c, j + 0);
    res += Tb * ncm_vector_get (sing_up->Up_c, j + 1);
    Ta = Ta * t / lfac;
    Tb = Tb * t / lfac;
    lfac += 1.0;
  }

  return res;
}

/**
 * ncm_csq1d_sing_fit_up_eval_dchi:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * 
 * Evaluates the model for $\chi^\prime$ at $t$.
 *  
 */
gdouble
ncm_csq1d_sing_fit_up_eval_dchi (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble m  = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  const gdouble dm = ncm_csq1d_eval_dm (csq1d, model, t, self->k);
  gdouble Ta       = 1.0 / m;
  gdouble Ta1      = -t * dm / (m * m);
  gdouble Tb       = 1.0 / t;
  gdouble lfac     = 1.0;
  gdouble res      = 0.0;
  gdouble n        = 0.0;
  gint j;

  for (j = 0; j < sing_up->chi_dim; j += 2)
  {
    res += Ta  * ncm_vector_get (sing_up->chi_c, j + 0) * (n + 1.0);
    res += Ta1 * ncm_vector_get (sing_up->chi_c, j + 0);
    res += Tb  * ncm_vector_get (sing_up->chi_c, j + 1) * (n + 0.0);

    Ta  = Ta  * t / lfac;
    Ta1 = Ta1 * t / lfac;
    Tb  = Tb  * t / lfac;

    lfac += 1.0;
    n    += 1.0;
  }

  return res;
}

/**
 * ncm_csq1d_sing_fit_up_eval_dexp_Up:
 * @sing_up: a #NcmCSQ1DSingFitUp
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * 
 * Evaluates the model for $\exp(\Upsilon_+)^\prime$ at $t$.
 *  
 */
gdouble
ncm_csq1d_sing_fit_up_eval_dexp_Up (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble m  = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  const gdouble dm = ncm_csq1d_eval_dm (csq1d, model, t, self->k);
  gdouble Ta       = t;
  gdouble Tb       = m;
  gdouble Tb1      = dm * t;
  gdouble lfac     = 1.0;
  gdouble res      = 0.0;
  gdouble n        = 0.0;
  gint j;

  for (j = 1; j < sing_up->Up_dim; j += 2)
  {
    res += Ta  * ncm_vector_get (sing_up->Up_c, j + 0) * (n + 2.0);
    res += Tb  * ncm_vector_get (sing_up->Up_c, j + 1) * (n + 1.0);
    res += Tb1 * ncm_vector_get (sing_up->Up_c, j + 1);

    Ta  = Ta  * t / lfac;
    Tb  = Tb  * t / lfac;
    Tb1 = Tb1 * t / lfac;

    lfac += 1.0;
    n    += 1.0;
  }

  return res;
}


/**
 * ncm_csq1d_ref:
 * @csq1d: a #NcmCSQ1D
 *
 * Increases the reference count of @csq1d.
 *
 * Returns: (transfer full): @csq1d.
 */
NcmCSQ1D *
ncm_csq1d_ref (NcmCSQ1D *csq1d)
{
  return g_object_ref (csq1d);
}

/**
 * ncm_csq1d_free:
 * @csq1d: a #NcmCSQ1D
 *
 * Decreases the reference count of @csq1d.
 *
 */
void
ncm_csq1d_free (NcmCSQ1D *csq1d)
{
  g_object_unref (csq1d);
}

/**
 * ncm_csq1d_clear:
 * @csq1d: a #NcmCSQ1D
 *
 * Decreases the reference count of *@csq1d and sets the pointer *@csq1d to NULL.
 *
 */
void
ncm_csq1d_clear (NcmCSQ1D **csq1d)
{
  g_clear_object (csq1d);
}

/**
 * ncm_csq1d_set_reltol:
 * @csq1d: a #NcmCSQ1D
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance to @reltol.
 *
 */
void 
ncm_csq1d_set_reltol (NcmCSQ1D *csq1d, const gdouble reltol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->reltol = reltol;
}

/**
 * ncm_csq1d_set_abstol:
 * @csq1d: a #NcmCSQ1D
 * @abstol: absolute tolerance
 *
 * Sets the absolute tolerance to @abstol.
 *
 */
void 
ncm_csq1d_set_abstol (NcmCSQ1D *csq1d, const gdouble abstol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->abstol = abstol;
}

/**
 * ncm_csq1d_set_k:
 * @csq1d: a #NcmCSQ1D
 * @k: mode $k$
 *
 * Sets the mode $k$ to @k.
 *
 */
void 
ncm_csq1d_set_k (NcmCSQ1D *csq1d, const gdouble k)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->k = k;
}

/**
 * ncm_csq1d_set_ti:
 * @csq1d: a #NcmCSQ1D
 * @ti: mode $t_i$
 *
 * Sets the initial time $t_i$ to @ti.
 *
 */
void 
ncm_csq1d_set_ti (NcmCSQ1D *csq1d, const gdouble ti)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->ti = ti;
}

/**
 * ncm_csq1d_set_tf:
 * @csq1d: a #NcmCSQ1D
 * @tf: mode $t_f$
 *
 * Sets the initial time $t_f$ to @tf.
 *
 */
void 
ncm_csq1d_set_tf (NcmCSQ1D *csq1d, const gdouble tf)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->tf = tf;
}

/**
 * ncm_csq1d_set_adiab_threshold:
 * @csq1d: a #NcmCSQ1D
 * @adiab_threshold: mode $A_t$
 *
 * Sets the adiabatic threshold $A_t$.
 *
 */
void 
ncm_csq1d_set_adiab_threshold (NcmCSQ1D *csq1d, const gdouble adiab_threshold)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->adiab_threshold = adiab_threshold;
}

/**
 * ncm_csq1d_set_save_evol:
 * @csq1d: a #NcmCSQ1D
 * @save: whether to save all evolution
 *
 * If true saves all evolution to be evaluted later through FIXME
 *
 */
void 
ncm_csq1d_set_save_evol (NcmCSQ1D *csq1d, gboolean save_evol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  if (self->save_evol != save_evol)
  {
    ncm_model_ctrl_force_update (self->ctrl);    
    self->save_evol = save_evol;
  }
}

/**
 * ncm_csq1d_set_sing_detect:
 * @csq1d: a #NcmCSQ1D
 * @enable: whether to try to detect all mass singularities
 *
 * If true it tries to detect all singularities caused
 * by the mass crossing $0$ or diverging to $\pm\infty$.
 *
 */
void 
ncm_csq1d_set_sing_detect (NcmCSQ1D *csq1d, gboolean enable)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  if (self->sing_detect != enable)
  {
    ncm_model_ctrl_force_update (self->ctrl);    
    self->sing_detect = enable;
  }
}

/**
 * ncm_csq1d_get_reltol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the relative tolerance.
 */
gdouble 
ncm_csq1d_get_reltol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->reltol;
}

/**
 * ncm_csq1d_get_abstol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the absolute tolerance to @abstol.
 */
gdouble
ncm_csq1d_get_abstol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->abstol;
}

/**
 * ncm_csq1d_get_k:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the mode $k$.
 */
gdouble
ncm_csq1d_get_k (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->k;
}

/**
 * ncm_csq1d_get_ti:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the initial time $t_i$.
 */
gdouble
ncm_csq1d_get_ti (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->ti;
}

/**
 * ncm_csq1d_get_tf:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the initial time $t_f$.
 */
gdouble
ncm_csq1d_get_tf (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->tf;
}

/**
 * ncm_csq1d_get_adiab_threshold:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the adiabatic threshold $A_t$.
 */
gdouble
ncm_csq1d_get_adiab_threshold (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->adiab_threshold;
}

/**
 * ncm_csq1d_get_save_evol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: whether the evolution will be saved.
 */
gboolean
ncm_csq1d_get_save_evol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->save_evol;
}

/**
 * ncm_csq1d_get_sing_detect:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: whether the evolution will be saved.
 */
gboolean
ncm_csq1d_get_sing_detect (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->sing_detect;
}

typedef struct _NcmCSQ1DWS
{
  NcmCSQ1D *csq1d;
  NcmModel *model;
  gdouble reltol;
} NcmCSQ1DWS;

static gint _ncm_csq1d_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint _ncm_csq1d_f_Up (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Up (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint _ncm_csq1d_f_Um (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Um (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void
_ncm_csq1d_set_init_cond (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t0)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;

  gdouble alpha, dgamma;

  ncm_csq1d_eval_adiab_at (csq1d, model, t0, &alpha, &dgamma, NULL, NULL);

  NV_Ith_S (self->y, 0) = alpha;
  NV_Ith_S (self->y, 1) = dgamma;

  self->t     = t0;
  self->state = NCM_CSQ1D_EVOL_STATE_ADIABATIC;
}

static void
_ncm_csq1d_prepare_integrator (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  gint flag;

  if (!self->cvode_init)
  {
    self->cvode = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode, &_ncm_csq1d_f, self->t, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode, self->t, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode, 0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode, &_ncm_csq1d_J);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode, fabs (self->t) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

  flag = CVodeSetUserData (self->cvode, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode, self->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetMaxStep (self->cvode, (self->tf - self->t) / 100.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
}

static void
_ncm_csq1d_prepare_integrator_Up (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  gint flag;

  if (!self->cvode_Up_init)
  {
    self->cvode_Up = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode_Up, &_ncm_csq1d_f_Up, self->t, self->y_Up);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_Up_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode_Up, self->t, self->y_Up);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode_Up, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode_Up, 100000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode_Up, self->LS_Up, self->A_Up);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode_Up, &_ncm_csq1d_J_Up);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode_Up, fabs (self->t) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );  

  flag = CVodeSetUserData (self->cvode_Up, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode_Up, self->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetMaxStep (self->cvode_Up, (self->tf - self->t) / 100.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
/*
  flag = CVodeSetMaxOrd (self->cvode_Up, 2);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxOrd", 1, );

  flag = CVodeSetMinStep (self->cvode_Up, 1.0e-5);
  NCM_CVODE_CHECK (&flag, "CVodeSetMinStep", 1, );
*/
}

static void
_ncm_csq1d_prepare_integrator_Um (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  gint flag;

  if (!self->cvode_Um_init)
  {
    self->cvode_Um = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode_Um, &_ncm_csq1d_f_Um, self->t, self->y_Um);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_Um_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode_Um, self->t, self->y_Um);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode_Um, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode_Um, 100000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode_Um, self->LS_Um, self->A_Um);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode_Um, &_ncm_csq1d_J_Um);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode_Um, fabs (self->t) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );  

  flag = CVodeSetUserData (self->cvode_Um, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode_Um, self->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetMaxStep (self->cvode_Um, (self->tf - self->t) / 100.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
}

static gint
_ncm_csq1d_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble alpha  = NV_Ith_S (y, 0);
  const gdouble dgamma = NV_Ith_S (y, 1);

  const gdouble nu     = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble F1     = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu  = 2.0 * nu;
  
  NV_Ith_S (ydot, 0) = - twonu * sinh (dgamma);
  NV_Ith_S (ydot, 1) = + twonu * (-F1 + cosh (dgamma) * tanh (alpha));

  return 0;
}

static gint
_ncm_csq1d_f_Up (realtype t, N_Vector y_Up, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble alpha    = NV_Ith_S (y_Up, 0);
  const gdouble Up       = NV_Ith_S (y_Up, 1);

  const gdouble m        = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2      = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Up   = exp (Up);
  const gdouble chi      = sinh (alpha);
  const gdouble ch_alpha = cosh (alpha);

  /*printf ("### % 22.15g % 22.15g % 22.15g % 22.15g\n", t, + m * nu2 * ch_alpha / exp_Up, - exp_Up / (m * ch_alpha), + 2.0 * m * nu2 * chi / exp_Up);*/
  
  NV_Ith_S (ydot, 0) = + m * nu2 * ch_alpha / exp_Up - exp_Up / (m * ch_alpha);
  NV_Ith_S (ydot, 1) = + 2.0 * m * nu2 * chi / exp_Up;
  
  return 0;
}

static gint
_ncm_csq1d_f_Um (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble alpha = NV_Ith_S (y, 0);
  const gdouble Um    = NV_Ith_S (y, 1);

  const gdouble m        = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2      = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Um   = exp (Um);
  const gdouble chi      = sinh (alpha);
  const gdouble ch_alpha = cosh (alpha);
/*
  printf ("### % 22.15g % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g\n", 
          t, + m * nu2 * exp_Um / ch_alpha, - ch_alpha / (m * exp_Um), - 2.0 * chi / (m * exp_Um),
          alpha, Um, m);
*/  
  NV_Ith_S (ydot, 0) = + m * nu2 * exp_Um / ch_alpha - ch_alpha / (m * exp_Um);
  NV_Ith_S (ydot, 1) = - 2.0 * chi / (m * exp_Um);

  return 0;
}

static gint
_ncm_csq1d_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble alpha  = NV_Ith_S (y, 0);
  const gdouble dgamma = NV_Ith_S (y, 1);
  
  const gdouble nu     = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu  = 2.0 * nu;

  /* - twonu * sinh (dgamma); */
  SM_ELEMENT_D (J, 0, 0) = 0.0;
  SM_ELEMENT_D (J, 0, 1) = - twonu * cosh (dgamma);

  /* + twonu * (-F1 + cosh (dgamma) * tanh (alpha)); */
  SM_ELEMENT_D (J, 1, 0) = + twonu * cosh (dgamma) / gsl_pow_2 (cosh (alpha));
  SM_ELEMENT_D (J, 1, 1) = + twonu * sinh (dgamma) * tanh (alpha);
  
  return 0;
}

static gint
_ncm_csq1d_J_Up (realtype t, N_Vector y_Up, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble alpha  = NV_Ith_S (y_Up, 0);
  const gdouble Up     = NV_Ith_S (y_Up, 1);
  
  const gdouble nu2      = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble m        = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Up   = exp (Up);
  const gdouble chi      = sinh (alpha);
  const gdouble ch_alpha = cosh (alpha);

  /* + m * nu2 * ch_alpha / exp_Up - exp_Up / (m * ch_alpha); */
  SM_ELEMENT_D (J, 0, 0) = + m * nu2 * chi / exp_Up + exp_Up * tanh (alpha) / (m * ch_alpha);
  SM_ELEMENT_D (J, 0, 1) = - m * nu2 * ch_alpha / exp_Up - exp_Up / (m * ch_alpha);

  /* + 2.0 * m * nu2 * chi / exp_Up; */
  SM_ELEMENT_D (J, 1, 0) = + 2.0 * m * nu2 * ch_alpha / exp_Up;
  SM_ELEMENT_D (J, 1, 1) = - 2.0 * m * nu2 * chi   / exp_Up;

  return 0;
}

static gint
_ncm_csq1d_J_Um (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble alpha   = NV_Ith_S (y, 0);
  const gdouble Um      = NV_Ith_S (y, 1);

  const gdouble m        = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2      = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Um   = exp (Um);
  const gdouble chi      = sinh (alpha);
  const gdouble ch_alpha = cosh (alpha);

  /* + m * nu2 * exp_Um / ch_alpha - ch_alpha / (m * exp_Um); */
  SM_ELEMENT_D (J, 0, 0) = + m * nu2 * exp_Um * tanh (alpha) / ch_alpha - chi / (m * exp_Um);
  SM_ELEMENT_D (J, 0, 1) = + m * nu2 * exp_Um / ch_alpha + ch_alpha / (m * exp_Um);

  /* - 2.0 * chi / (m * exp_Um); */
  SM_ELEMENT_D (J, 1, 0) = - 2.0 * ch_alpha / (m * exp_Um);
  SM_ELEMENT_D (J, 1, 1) = + 2.0 * chi      / (m * exp_Um);

  return 0;
}

/**
 * ncm_csq1d_eval_xi: (virtual eval_xi)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\xi$ 
 */
/**
 * ncm_csq1d_eval_dxi: (virtual eval_dxi)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\mathrm{d}\xi/\mathrm{d}t$ 
 */
/**
 * ncm_csq1d_eval_nu: (virtual eval_nu)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\nu$ 
 */
/**
 * ncm_csq1d_eval_nu2: (virtual eval_nu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\nu^2$ 
 */
/**
 * ncm_csq1d_eval_m: (virtual eval_m)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $m$ 
 */
/**
 * ncm_csq1d_eval_dm: (virtual eval_dm)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\mathrm{d}m/\mathrm{d}t$ 
 */
/**
 * ncm_csq1d_eval_F1: (virtual eval_F1)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $F_1$ 
 */
/**
 * ncm_csq1d_eval_F2: (virtual eval_F2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $F_2$ 
 */
/**
 * ncm_csq1d_eval_FN: (virtual eval_FN)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @n: order $n$
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $F_n$ 
 */
/**
 * ncm_csq1d_eval_powspec_factor: (virtual eval_powspec_factor)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @k: mode $k$
 *
 * Returns: $F_n$ 
 */

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_adiabatic (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop reason = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t    = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_ADIABATIC);

  if (PRINT_EVOL)
    printf ("# ENTER EVOL ADIABATIC: TIME % 22.15g\n", self->t);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble alpha, dgamma;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode, self->tf, self->y, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_adiabatic]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    alpha   = NV_Ith_S (self->y, 0);
    dgamma  = NV_Ith_S (self->y, 1);

    if (PRINT_EVOL)
    {
      const gdouble xi         = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
      const gdouble gamma      = dgamma + xi;
      const gdouble chi        = sinh (alpha);
      const gdouble lnch_alpha = gsl_sf_lncosh (alpha);
      const gdouble Up         = + gamma + lnch_alpha;
      const gdouble Um         = - gamma + lnch_alpha;
      printf ("# E[AD] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g\n", self->t, alpha, dgamma, chi, Up, Um, xi);
    }

    is_finished = (self->t == self->tf);
    
    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }
    
    if ((fabs (dgamma) > self->adiab_threshold) || (fabs (alpha) > self->adiab_threshold))
    {
      if (dgamma > 0.0)
        reason = NCM_CSQ1D_EVOL_STOP_UP_START;
      else
        reason = NCM_CSQ1D_EVOL_STOP_UM_START;
      break;
    }

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }
  }
  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Up (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop reason    = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t       = asinh (self->t);
  GArray *t_a                = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *exp_Up_a           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *chim_t_a           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_UP);
  if (PRINT_EVOL)
    printf ("# ENTER EVOL UP: TIME % 22.15g\n", self->t);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble alpha, chim, chi, dgamma, gamma, Up, xi, m, exp_Up, chim_t;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode_Up, self->tf, self->y_Up, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_Up]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
    m       = ncm_csq1d_eval_m (csq1d, model, self->t, self->k);
    alpha   = NV_Ith_S (self->y_Up, 0);
    Up      = NV_Ith_S (self->y_Up, 1);
    chi     = sinh (alpha);
    chim    = m * chi;
    gamma   = Up - gsl_sf_lncosh (alpha);
    dgamma  = gamma - xi;
    exp_Up  = exp (Up);
    chim_t  = chim / self->t;

    if (PRINT_EVOL)
    {
      const gdouble Um = - Up + 2.0 * gsl_sf_lncosh (alpha);
      printf ("# E[UP] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g\n", self->t, alpha, dgamma, chi, Up, Um, xi);
    }

    is_finished = (self->t == self->tf);
    
    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }
    
    g_array_append_val (t_a,      self->t);
    g_array_append_val (exp_Up_a, exp_Up);
    g_array_append_val (chim_t_a, chim_t);
    
    if ((fabs (dgamma) < self->adiab_threshold) && (fabs (chi) < self->adiab_threshold))
    {
      reason = NCM_CSQ1D_EVOL_STOP_ADIABATIC_START;
      break;
    }

    if (dgamma < 0.0)
    {
      reason = NCM_CSQ1D_EVOL_STOP_UM_START;
      break;
    }

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }
  }

  g_array_unref (t_a);
  g_array_unref (exp_Up_a);
  g_array_unref (chim_t_a);

  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Um (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop reason = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t    = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_UM);
  if (PRINT_EVOL)
    printf ("# ENTER EVOL UM: TIME % 22.15g\n", self->t);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble alpha, chi, dgamma, gamma, Um, xi;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode_Um, self->tf, self->y_Um, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_Um]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
    alpha   = NV_Ith_S (self->y_Um, 0);
    Um      = NV_Ith_S (self->y_Um, 1);
    chi     = sinh (alpha);
    gamma   = -Um + gsl_sf_lncosh (alpha);
    dgamma  = gamma - xi;

    if (PRINT_EVOL)
    {
      const gdouble Up        = - Um + 2.0 * gsl_sf_lncosh (alpha);
      printf ("# E[UM] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g\n", self->t, alpha, dgamma, chi, Up, Um, xi);
    }

    is_finished = (self->t == self->tf);
    
    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }
    
    if ((fabs (dgamma) < self->adiab_threshold) && (fabs (chi) < self->adiab_threshold))
    {
      reason = NCM_CSQ1D_EVOL_STOP_ADIABATIC_START;
      break;
    }
    if (dgamma > 0.0)
    {
      reason = NCM_CSQ1D_EVOL_STOP_UP_START;
      break;
    }
    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }
  }
  return reason;
}


static void
_ncm_csq1d_evol_save (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DWS *ws, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop stop;

  if (PRINT_EVOL)
    printf ("# ENTERING EVOL SAVE: STATE %d, TIME % 22.15g\n", self->state, self->t);
  
  switch (self->state)
  {
    case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
    {
      _ncm_csq1d_prepare_integrator (csq1d, ws);
      stop = _ncm_csq1d_evol_adiabatic (csq1d, model, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STATE_UP:
    {
      _ncm_csq1d_prepare_integrator_Up (csq1d, ws);
      stop = _ncm_csq1d_evol_Up (csq1d, ws, model, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STATE_UM:
    {
      _ncm_csq1d_prepare_integrator_Um (csq1d, ws);
      stop = _ncm_csq1d_evol_Um (csq1d, model, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  switch (stop)
  {
    case NCM_CSQ1D_EVOL_STOP_UP_START:
    {
      switch (self->state)
      {
        case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
        {
          const gdouble alpha  = NV_Ith_S (self->y, 0);
          const gdouble dgamma = NV_Ith_S (self->y, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = xi + dgamma;
          const gdouble Up     = + gamma + gsl_sf_lncosh (alpha);

          NV_Ith_S (self->y_Up, 0) = alpha;
          NV_Ith_S (self->y_Up, 1) = Up;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UM:
        {
          const gdouble alpha = NV_Ith_S (self->y_Um, 0);
          const gdouble Um    = NV_Ith_S (self->y_Um, 1);
          const gdouble Up    = -Um + 2.0 * gsl_sf_lncosh (alpha);

          NV_Ith_S (self->y_Up, 0) = alpha;
          NV_Ith_S (self->y_Up, 1) = Up;
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }
      self->state = NCM_CSQ1D_EVOL_STATE_UP;
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STOP_UM_START:
    {
      switch (self->state)
      {
        case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
        {
          const gdouble alpha  = NV_Ith_S (self->y, 0);
          const gdouble dgamma = NV_Ith_S (self->y, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = xi + dgamma;
          const gdouble Um     = - gamma + gsl_sf_lncosh (alpha);

          NV_Ith_S (self->y_Um, 0) = alpha;
          NV_Ith_S (self->y_Um, 1) = Um;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UP:
        {
          const gdouble alpha = NV_Ith_S (self->y_Up, 0);
          const gdouble Up    = NV_Ith_S (self->y_Up, 1);
          const gdouble Um    = -Up + 2.0 * gsl_sf_lncosh (alpha);

          NV_Ith_S (self->y_Um, 0) = alpha;
          NV_Ith_S (self->y_Um, 1) = Um;
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }

      self->state = NCM_CSQ1D_EVOL_STATE_UM;
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STOP_ADIABATIC_START:
    {
      switch (self->state)
      {
        case NCM_CSQ1D_EVOL_STATE_UP:
        {
          const gdouble alpha  = NV_Ith_S (self->y_Up, 0);
          const gdouble Up     = NV_Ith_S (self->y_Up, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = Up - gsl_sf_lncosh (alpha);
          const gdouble dgamma = gamma - xi;

          NV_Ith_S (self->y, 0) = alpha;
          NV_Ith_S (self->y, 1) = dgamma;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UM:
        {
          const gdouble alpha  = NV_Ith_S (self->y_Um, 0);
          const gdouble Um     = NV_Ith_S (self->y_Um, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = - Um + gsl_sf_lncosh (alpha);
          const gdouble dgamma = gamma - xi;

          NV_Ith_S (self->y, 0) = alpha;
          NV_Ith_S (self->y, 1) = dgamma;
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }

      self->state = NCM_CSQ1D_EVOL_STATE_ADIABATIC;      
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STOP_FINISHED:
      return;
      break;
    default:
      g_error ("_ncm_csq1d_evol_save: error evolving the ode system.");
      break;
  }  
}

static gdouble
_ncm_csq1d_mass_root (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;
  return ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
}

static gdouble
_ncm_csq1d_sing_detect (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws, NcmModel *model, gdouble ti, gdouble tf)
{
  /*NcmCSQ1DPrivate * const self = csq1d->priv;*/
  gdouble tsing = 0.5 * (tf + ti);
  gint iter = 0, max_iter = 1000;
  const gdouble root_reltol = 1.0e-12;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gint status;

  F.function = &_ncm_csq1d_mass_root;
  F.params   = ws;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, ti, tf);

  if (PRINT_EVOL)
    printf ("# Searching for singularities:\n");
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    tsing  = gsl_root_fsolver_root (s);
    ti     = gsl_root_fsolver_x_lower (s);
    tf     = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (ti, tf, 0.0, root_reltol);

    if (PRINT_EVOL)
      printf ("#   %5d [%.7e, %.7e] %.7e %+.7e\n", iter, ti, tf, tsing, tf - ti);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return tsing;
}

/**
 * ncm_csq1d_prepare: (virtual prepare)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 *
 * Prepares the object using @model.
 *
 */
void 
ncm_csq1d_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DWS ws = {csq1d, model, 0.0};
  
  if (NCM_CSQ1D_GET_CLASS (csq1d)->prepare != NULL)
  {
    NCM_CSQ1D_GET_CLASS (csq1d)->prepare (csq1d, model);
  }

  g_assert_cmpfloat (self->tf, >, self->ti);

  if (self->sing_detect)
    _ncm_csq1d_sing_detect (csq1d, &ws, model, self->ti, self->tf);
  _ncm_csq1d_set_init_cond (csq1d, model, self->ti);

  if (self->save_evol)
  {
    GArray *asinh_t_a    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *alpha_a      = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *dgamma_a     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);

    _ncm_csq1d_evol_save (csq1d, model, &ws, asinh_t_a, alpha_a, dgamma_a);

    ncm_spline_set_array (self->alpha_s,  asinh_t_a, alpha_a,  TRUE);
    ncm_spline_set_array (self->dgamma_s, asinh_t_a, dgamma_a, TRUE);

    g_array_unref (asinh_t_a);
    g_array_unref (alpha_a);
    g_array_unref (dgamma_a);
  }
  else
  {
    g_assert_not_reached ();
  }
}

/**
 * ncm_csq1d_get_time_array:
 * @csq1d: a #NcmCSQ1D
 * @smallest_t: (out) (allow-none): the smallest absolute value of $t$ in the array
 *
 * Returns: (transfer full) (element-type gdouble): the time array of the computed steps.
 */
GArray *
ncm_csq1d_get_time_array (NcmCSQ1D *csq1d, gdouble *smallest_t)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmVector *asinh_t_v = ncm_spline_get_xv (self->alpha_s);
  const guint len      = ncm_vector_len (asinh_t_v);
  GArray *t_a          = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);
  gdouble s_t          = 1.0e300;
  gint i;
  
  for (i = 0; i < len; i++)
  {
    const gdouble t_i = sinh (ncm_vector_fast_get (asinh_t_v, i));
    const gdouble a_t = fabs (t_i);
    
    g_array_append_val (t_a, t_i);
    s_t = MIN (a_t, s_t);
  }

  if (smallest_t != NULL)
    smallest_t[0] = s_t;

  return t_a;
}

gdouble
_ncm_csq1d_find_adiab_time_limit_f (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) params;
  /*NcmCSQ1DPrivate * const self = ws->csq1d->priv;*/

  gdouble alpha, dgamma, cmp, alpha_reltol, dgamma_reltol;

  ncm_csq1d_eval_adiab_at (ws->csq1d, ws->model, t, &alpha, &dgamma, &alpha_reltol, &dgamma_reltol);

  cmp = MAX (alpha_reltol, dgamma_reltol);

  return log (cmp / ws->reltol);
}

/**
 * ncm_csq1d_find_adiab_time_limit:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t0: time lower bound $t_0$
 * @t1: time upper bound $t_1$
 * @reltol: relative tolerance
 * @ti: (out): adiabatic time limit $t_i$
 *
 * Computes the time upper limit $t_i \in [t_0, t_1]$ where the adiabatic
 * approximation is satisfied up to @reltol.
 * 
 * Returns: whether the time limit was found.
 */
gboolean 
ncm_csq1d_find_adiab_time_limit (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble reltol, gdouble *ti)
{
  /*NcmCSQ1DPrivate * const self = csq1d->priv;*/
  gdouble alpha0, dgamma0, alpha_reltol0, dgamma_reltol0;
  gdouble alpha1, dgamma1, alpha_reltol1, dgamma_reltol1;
  
  gboolean adiab0, adiab1;

  g_assert_cmpfloat (reltol, >, 0.0);

  ncm_csq1d_eval_adiab_at (csq1d, model, t0, &alpha0, &dgamma0, &alpha_reltol0, &dgamma_reltol0);
  ncm_csq1d_eval_adiab_at (csq1d, model, t1, &alpha1, &dgamma1, &alpha_reltol1, &dgamma_reltol1);

  adiab0 = ((fabs (alpha_reltol0) < reltol) && (fabs (dgamma_reltol0) < reltol));
  adiab1 = ((fabs (alpha_reltol1) < reltol) && (fabs (dgamma_reltol1) < reltol));
  
  if ((adiab0 && adiab1) || (!adiab0 && !adiab1))
  {
    if (PRINT_EVOL)
      g_warning ("# Impossible to find the adiabatic limit: t0 % 22.15g % 22.15g % 22.15g t1 % 22.15g % 22.15g % 22.15g\n", t0, alpha0, dgamma0, t1, alpha1, dgamma1);
    return FALSE;
  }
  else
  {
    NcmCSQ1DWS ws = {csq1d, model, reltol};
    gint iter = 0, max_iter = 1000;
    const gdouble root_reltol = 1.0e-2;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    gint status;

    ti[0] = 0.5 * (t0 + t1);

    F.function = &_ncm_csq1d_find_adiab_time_limit_f;
    F.params   = &ws;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, t0, t1);

    if (PRINT_EVOL)
      printf ("# Searching for the adiabatic time limit: \n");
    
    do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      ti[0]  = gsl_root_fsolver_root (s);
      t0     = gsl_root_fsolver_x_lower (s);
      t1     = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (t0, t1, 0.0, root_reltol);

      if (PRINT_EVOL)
        printf ("#  %5d [%.7e, %.7e] %.7e %+.7e\n", iter, t0, t1, ti[0], t1 - t0);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);    
  }

  return TRUE;
}

static gdouble 
_ncm_csq1d_F2_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;
  const gdouble F2 = ncm_csq1d_eval_F2 (ws->csq1d, ws->model, t, self->k);

  return F2;
}

static gdouble 
_ncm_csq1d_lnnu_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;
  const gdouble nu = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);

  return log (nu);
}


/**
 * ncm_csq1d_eval_adiab_at:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @alpha: (out): value of $\alpha(t)$
 * @dgamma: (out): value of $\Delta\gamma(t)$
 * @alpha_reltol: (out) (allow-none): estimated error on $\alpha(t)$
 * @dgamma_reltol: (out) (allow-none): estimated error on $\Delta\gamma(t)$
 *
 * Computes the value of the adiabatic approximation of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 *
 */
void 
ncm_csq1d_eval_adiab_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DWS ws       = {csq1d, model, 0.0};
  const gdouble F1    = ncm_csq1d_eval_F1 (csq1d, model, t, self->k);
  const gdouble F2    = ncm_csq1d_eval_F2 (csq1d, model, t, self->k);
  const gdouble F1_2  = F1 * F1;
  const gdouble F1_3  = F1_2 * F1;
  const gdouble nu    = ncm_csq1d_eval_nu (csq1d, model, t, self->k);
  const gdouble twonu = 2.0 * nu;
  gdouble err, F3, d2F2, dlnnu, F4, alpha_reltol0, dgamma_reltol0;

  F3             = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err) / twonu;
  d2F2           = ncm_diff_rc_d2_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err);
  dlnnu          = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_lnnu_func, &ws, &err);
  F4             = d2F2 / gsl_pow_2 (twonu) - dlnnu * F3 / twonu;
  alpha_reltol0  = gsl_pow_2 ((F1_3 / 3.0 - F3) / F1);
  dgamma_reltol0 = gsl_pow_2 ((F4 - F1_2 * F2) / F2);

  alpha[0]  = + F1 + F1_3 / 3.0 - F3;
  dgamma[0] = - (1.0 + F1_2) * F2 + F4;

  if (alpha_reltol != NULL)
    alpha_reltol[0] = alpha_reltol0;
  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = dgamma_reltol0;
}

/**
 * ncm_csq1d_eval_at:
 * @csq1d: a #NcmCSQ1D
 * @t: time $t$
 * @alpha: (out): value of $\alpha(t)$
 * @dgamma: (out): value of $\Delta\gamma(t)$
 *
 * Computes the value of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 *
 */
void 
ncm_csq1d_eval_at (NcmCSQ1D *csq1d, const gdouble t, gdouble *alpha, gdouble *dgamma)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble a_t = asinh (t);
  
  alpha[0]  = ncm_spline_eval (self->alpha_s, a_t);
  dgamma[0] = ncm_spline_eval (self->dgamma_s, a_t);
}

/**
 * ncm_csq1d_alpha_dgamma_to_phi_Pphi:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @alpha: value of $\alpha(t)$
 * @dgamma: value of $\Delta\gamma(t)$
 * @phi: (out caller-allocates) (array fixed-size=2): real and imaginary parts of $\phi$, i.e., $[\mathrm{Re}(\phi), \mathrm{Im}(\phi)]$
 * @Pphi: (out caller-allocates) (array fixed-size=2): real and imaginary parts of $\Pi_\phi$, i.e., $[\mathrm{Re}(\Pi_\phi), \mathrm{Im}(\Pi_\phi)]$
 *
 * Computes the value of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 *
 */
void 
ncm_csq1d_alpha_dgamma_to_phi_Pphi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble alpha, const gdouble dgamma, gdouble *phi, gdouble *Pphi)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble gamma               = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;
  const gdouble exp_gamma_p_alpha_2 = exp (0.5 * (gamma + alpha));
  const gdouble exp_gamma_m_alpha_2 = exp (0.5 * (gamma - alpha));
  
  phi[0]  = +0.5 / exp_gamma_m_alpha_2;
  phi[1]  = -0.5 / exp_gamma_p_alpha_2;

  Pphi[0] = -0.5 * exp_gamma_p_alpha_2;
  Pphi[1] = -0.5 * exp_gamma_m_alpha_2;
}

/**
 * ncm_csq1d_get_J_at:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @J11: (out): $J_{11}$
 * @J12: (out): $J_{12}$
 * @J22: (out): $J_{22}$
 * 
 * Computes the complex structure matrix.
 * 
 */
void 
ncm_csq1d_get_J_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *J11, gdouble *J12, gdouble *J22)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble a_t = asinh (t);
  
  const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
  const gdouble gamma  = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  J11[0] = cosh (alpha) * exp (-gamma);
  J22[0] = cosh (alpha) * exp (+gamma);
  J12[0] = -sinh (alpha);
}


/*
    if (fabs (prev_Up / Up - 1.0) <= self->reltol * 1.0e-4)
      stable_Up++;
    else
      stable_Up = 0;

    prev_Up = Up;
    
    if (!fit && (stable_Up >= 3) && TRUE)
    {
      const gdouble l_chim_t = g_array_index (chim_t_a, gdouble, t_a->len - 1);
      const gdouble l_exp_Up = g_array_index (exp_Up_a, gdouble, t_a->len - 1);
      gint i, n, len;

      for (i = t_a->len - 1; i >= 0; i--)
      {
        const gdouble t_i        = g_array_index (t_a,      gdouble, i);
        const gdouble chim_t_i   = g_array_index (chim_t_a, gdouble, i);
        const gdouble exp_Up_i   = g_array_index (exp_Up_a, gdouble, i);
        const gdouble cmp_chim_t = fabs (chim_t_i / l_chim_t - 1.0);
        const gdouble cmp_exp_Up = fabs (exp_Up_i / l_exp_Up - 1.0);
        
        printf ("% 22.15g % 22.15g % 22.15g % 22.15g %e %e\n", t_i, t_i * ncm_csq1d_eval_dxi (csq1d, model, t_i, self->k),  chim_t_i, exp_Up_i, cmp_exp_Up, cmp_chim_t);

        if ((cmp_chim_t > 0.15) || (cmp_exp_Up > 0.15))
        {
          printf ("Found! %d % 22.15g len %d\n", i, t_i, t_a->len - i);
          break;
        }
      }
      n   = i;
      len = t_a->len - n;

      {
        NcmVector *t_v      = ncm_vector_new_data_static (&g_array_index (t_a,      gdouble, n), len, 1);
        NcmVector *chim_t_v = ncm_vector_new_data_static (&g_array_index (chim_t_a, gdouble, n), len, 1);
        NcmVector *exp_Up_v = ncm_vector_new_data_static (&g_array_index (exp_Up_a, gdouble, n), len, 1);
        gdouble f_err_chi = 0.0;
        gdouble f_err_Up  = 0.0;

        ncm_csq1d_sing_fit_up_fit (sing_up, csq1d, model, t_v, chim_t_v, exp_Up_v);

        for (i = len - 1; i >= 0; i--)
        {
          const gdouble t_i      = -ncm_vector_get (t_v, i);
          const gdouble xi       = ncm_csq1d_eval_xi (csq1d, model, t_i, self->k);
          const gdouble m        = ncm_csq1d_eval_m (csq1d, model, t_i, self->k);
          const gdouble nu2      = ncm_csq1d_eval_nu2 (csq1d, model, t_i, self->k);
          const gdouble chi      = ncm_csq1d_sing_fit_up_eval_chi (sing_up, csq1d, model, t_i);
          const gdouble exp_Up   = ncm_csq1d_sing_fit_up_eval_exp_Up (sing_up, csq1d, model, t_i);
          const gdouble Up       = log (exp_Up);
          const gdouble dexp_Up  = ncm_csq1d_sing_fit_up_eval_dexp_Up (sing_up, csq1d, model, t_i);
          const gdouble dchi     = ncm_csq1d_sing_fit_up_eval_dchi (sing_up, csq1d, model, t_i);
          const gdouble f_exp_Up = 2.0 * m * nu2 * chi;
          const gdouble f_chi    = -exp_Up / m + m * nu2 * (1.0 + chi * chi) / exp_Up;
          const gdouble err_chi  = fabs (dchi / f_chi - 1.0);
          const gdouble err_Up   = fabs (dexp_Up / f_exp_Up - 1.0);
          const gdouble alpha    = asinh (chi);
          const gdouble asinh_t  = asinh (t_i);
          const gdouble gamma    = Up - gsl_sf_lncosh (alpha);
          const gdouble dgamma   = gamma - xi;

          if (i + 1 == len)
          {
            f_err_chi = err_chi;
            f_err_Up  = err_Up;
          }

          g_array_append_val (asinh_t_a, asinh_t);
          g_array_append_val (alpha_a,   alpha);
          g_array_append_val (dgamma_a,  dgamma);
          last_asinh_t = asinh_t;
          
          
          printf ("# CMP[%5d] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 12.7e | % 22.15g % 22.15g % 12.7e\n",
                  i, t_i, chi, Up, dexp_Up, f_exp_Up, err_Up, dchi, f_chi, err_chi);

          if ((err_chi / f_err_chi > 1.1) || (err_Up / f_err_Up > 1.1))
          {
            self->t                  = t_i;
            NV_Ith_S (self->y_Up, 0) = alpha;
            NV_Ith_S (self->y_Up, 1) = Up;
            _ncm_csq1d_prepare_integrator_Up (csq1d, ws);
            break;
          }
        }
      }

      fit = TRUE;
    }

 */