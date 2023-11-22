/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_csq1d.c
 *
 *  Mon September 09 13:56:19 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_csq1d.c
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>

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

#define PRINT_EVOL FALSE

typedef enum _NcmCSQ1DEvolStop
{
  NCM_CSQ1D_EVOL_STOP_ERROR = -1,
  NCM_CSQ1D_EVOL_STOP_FINISHED = 0,
  NCM_CSQ1D_EVOL_STOP_ADIABATIC_START,
  NCM_CSQ1D_EVOL_STOP_UP_START,
  NCM_CSQ1D_EVOL_STOP_UM_START,
} NcmCSQ1DEvolStop;

typedef struct _NcmCSQ1DPrivate
{
  gdouble reltol;
  gdouble abstol;
  gdouble k;
  gdouble ti;
  gdouble tf;
  gdouble t;
  gboolean init_cond_set;
  NcmCSQ1DEvolState state;
  gdouble adiab_threshold;
  gdouble prop_threshold;
  gboolean save_evol;
  gboolean sing_detect;
  NcmModelCtrl *ctrl;
  gpointer cvode;
  gpointer cvode_Up;
  gpointer cvode_Um;
  gpointer cvode_Prop;
  gpointer arkode;
  gboolean cvode_init;
  gboolean cvode_Up_init;
  gboolean cvode_Um_init;
  gboolean cvode_Prop_init;
  gboolean arkode_init;
  N_Vector y;
  N_Vector y_Up;
  N_Vector y_Um;
  N_Vector y_Prop;
  SUNMatrix A;
  SUNMatrix A_Up;
  SUNMatrix A_Um;
  SUNMatrix A_Prop;
  SUNLinearSolver LS;
  SUNLinearSolver LS_Up;
  SUNLinearSolver LS_Um;
  SUNLinearSolver LS_Prop;
  NcmSpline *alpha_s;
  NcmSpline *dgamma_s;
  NcmDiff *diff;
  NcmSpline *R[4];
  gdouble tf_Prop;
} NcmCSQ1DPrivate;

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_K,
  PROP_TI,
  PROP_TF,
  PROP_ADIAB_THRESHOLD,
  PROP_PROP_THRESHOLD,
  PROP_SAVE_EVOL,
  PROP_SING_DETECT,
};

G_DEFINE_BOXED_TYPE (NcmCSQ1DSingFitUp, ncm_csq1d_sing_fit_up, ncm_csq1d_sing_fit_up_dup, ncm_csq1d_sing_fit_up_free)
G_DEFINE_TYPE_WITH_PRIVATE (NcmCSQ1D, ncm_csq1d, G_TYPE_OBJECT)

static void
ncm_csq1d_init (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  self->reltol          = 0.0;
  self->abstol          = 0.0;
  self->k               = 0.0;
  self->ti              = 0.0;
  self->tf              = 0.0;
  self->init_cond_set   = FALSE;
  self->state           = NCM_CSQ1D_EVOL_STATE_INVALID;
  self->adiab_threshold = 0.0;
  self->prop_threshold  = 0.0;
  self->save_evol       = FALSE;
  self->ctrl            = ncm_model_ctrl_new (NULL);

  self->cvode           = NULL;
  self->cvode_init      = FALSE;
  self->cvode_Up        = NULL;
  self->cvode_Up_init   = FALSE;
  self->cvode_Um        = NULL;
  self->cvode_Um_init   = FALSE;
  self->cvode_Prop      = NULL;
  self->cvode_Prop_init = FALSE;
  self->arkode          = NULL;
  self->arkode_init     = FALSE;

  self->y      = N_VNew_Serial (2);
  self->y_Up   = N_VNew_Serial (2);
  self->y_Um   = N_VNew_Serial (2);
  self->y_Prop = N_VNew_Serial (4);

  self->A = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer) self->A, "SUNDenseMatrix", 0, );

  self->A_Up = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer) self->A_Up, "SUNDenseMatrix", 0, );

  self->A_Um = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer) self->A_Um, "SUNDenseMatrix", 0, );

  self->A_Prop = SUNDenseMatrix (4, 4);
  NCM_CVODE_CHECK ((gpointer) self->A_Prop, "SUNDenseMatrix", 0, );

  self->LS = SUNDenseLinearSolver (self->y, self->A);
  NCM_CVODE_CHECK ((gpointer) self->LS, "SUNDenseLinearSolver", 0, );

  self->LS_Up = SUNDenseLinearSolver (self->y_Up, self->A_Up);
  NCM_CVODE_CHECK ((gpointer) self->LS_Up, "SUNDenseLinearSolver", 0, );

  self->LS_Um = SUNDenseLinearSolver (self->y_Um, self->A_Um);
  NCM_CVODE_CHECK ((gpointer) self->LS_Um, "SUNDenseLinearSolver", 0, );

  self->LS_Prop = SUNDenseLinearSolver (self->y_Prop, self->A_Prop);
  NCM_CVODE_CHECK ((gpointer) self->LS_Um, "SUNDenseLinearSolver", 0, );

  self->alpha_s  = ncm_spline_cubic_notaknot_new ();
  self->dgamma_s = ncm_spline_cubic_notaknot_new ();

  {
    gint i;

    for (i = 0; i < 4; i++)
      self->R[i] = ncm_spline_cubic_notaknot_new ();
  }
  self->tf_Prop = 0.0;

  self->diff = ncm_diff_new ();
}

static void
_ncm_csq1d_dispose (GObject *object)
{
  NcmCSQ1D *csq1d              = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  ncm_model_ctrl_clear (&self->ctrl);
  ncm_spline_clear (&self->alpha_s);
  ncm_spline_clear (&self->dgamma_s);

  {
    gint i;

    for (i = 0; i < 4; i++)
      ncm_spline_clear (&self->R[i]);
  }

  ncm_diff_clear (&self->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_csq1d_parent_class)->dispose (object);
}

static void
_ncm_csq1d_finalize (GObject *object)
{
  NcmCSQ1D *csq1d              = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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

  if (self->cvode_Prop != NULL)
  {
    CVodeFree (&self->cvode_Prop);
    self->cvode_Prop      = NULL;
    self->cvode_Prop_init = FALSE;
  }

  if (self->arkode != NULL)
  {
    ARKStepFree (&self->arkode);
    self->arkode      = NULL;
    self->arkode_init = FALSE;
  }

  g_clear_pointer (&self->y,      N_VDestroy);
  g_clear_pointer (&self->y_Up,   N_VDestroy);
  g_clear_pointer (&self->y_Um,   N_VDestroy);
  g_clear_pointer (&self->y_Prop, N_VDestroy);

  if (self->A != NULL)
    SUNMatDestroy (self->A);

  if (self->A_Up != NULL)
    SUNMatDestroy (self->A_Up);

  if (self->A_Um != NULL)
    SUNMatDestroy (self->A_Um);

  if (self->A_Prop != NULL)
    SUNMatDestroy (self->A_Prop);

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

  if (self->LS_Prop != NULL)
  {
    gint flag = SUNLinSolFree (self->LS_Prop);

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
    case PROP_PROP_THRESHOLD:
      ncm_csq1d_set_prop_threshold (csq1d, g_value_get_double (value));
      break;
    case PROP_SAVE_EVOL:
      ncm_csq1d_set_save_evol (csq1d, g_value_get_boolean (value));
      break;
    case PROP_SING_DETECT:
      ncm_csq1d_set_sing_detect (csq1d, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    case PROP_PROP_THRESHOLD:
      g_value_set_double (value, ncm_csq1d_get_prop_threshold (csq1d));
      break;
    case PROP_SAVE_EVOL:
      g_value_set_boolean (value, ncm_csq1d_get_save_evol (csq1d));
      break;
    case PROP_SING_DETECT:
      g_value_set_boolean (value, ncm_csq1d_get_sing_detect (csq1d));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _ncm_csq1d_eval_xi         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_dxi        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_nu         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_nu2        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_m          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_int_1_m    (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_int_mnu2   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_int_qmnu2  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_dm         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_F1         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_F2         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_FN         (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);

static void
ncm_csq1d_class_init (NcmCSQ1DClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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
                                   PROP_PROP_THRESHOLD,
                                   g_param_spec_double ("prop-threshold",
                                                        NULL,
                                                        "The propagator threshold",
                                                        0.0, 1.0, 1.0e-1,
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
  klass->eval_int_1_m        = &_ncm_csq1d_eval_int_1_m;
  klass->eval_int_mnu2       = &_ncm_csq1d_eval_int_mnu2;
  klass->eval_int_qmnu2      = &_ncm_csq1d_eval_int_qmnu2;
  klass->eval_int_q2mnu2     = &_ncm_csq1d_eval_int_q2mnu2;
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
_ncm_csq1d_eval_int_1_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_int_1_m: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_int_mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_int_mnu2: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_int_qmnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_int_qmnu2: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_int_q2mnu2: not implemented.");

  return 0.0;
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

typedef struct _NcmCSQ1DWS
{
  NcmCSQ1D *csq1d;
  NcmModel *model;
  gdouble reltol;
} NcmCSQ1DWS;

static gdouble _ncm_csq1d_F1_func (const gdouble t, gpointer user_data);

static gdouble
_ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0};
  const gdouble nu             = ncm_csq1d_eval_nu (csq1d, model, t, self->k);
  const gdouble twonu          = 2.0 * nu;

  gdouble err;

  return ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F1_func, &ws, &err) / twonu;
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
  NcmCSQ1DSingFitUp *sing_up_dup = g_new0 (NcmCSQ1DSingFitUp, 1);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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

    /*printf ("Fitting chi * m / t.\n");*/
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

        Ta    = Ta * t_i / lfac;
        Tb    = Tb * t_i / lfac;
        lfac += 1.0;
      }
    }

    ret = gsl_multifit_linear (ncm_matrix_gsl (X), ncm_vector_gsl (chim_t), ncm_vector_gsl (sing_up->chi_c), ncm_matrix_gsl (cov), &chisq, work);
    g_assert_cmpint (ret, ==, 0);

    /*printf ("MULTIFIT RET: %d, chisq % 22.15g, sqrt (chisq / len) % 22.15g\n", ret, chisq, sqrt (chisq / len));*/

    /*ncm_vector_log_vals (sing_up->chi_c, "COEFS: ", "% 22.15g", TRUE);*/
    /*ncm_matrix_log_vals (cov,            "COV:  ", "% 22.15g");*/

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

        Ta    = Ta * t_i / lfac;
        Tb    = Tb * t_i / lfac;
        lfac += 1.0;
      }
    }

    ret = gsl_multifit_linear (ncm_matrix_gsl (X), ncm_vector_gsl (exp_Up), ncm_vector_gsl (sing_up->Up_c), ncm_matrix_gsl (cov), &chisq, work);
    g_assert_cmpint (ret, ==, 0);

    /*printf ("MULTIFIT RET: %d, chisq % 22.15g, sqrt (chisq / len) % 22.15g\n", ret, chisq, sqrt (chisq / len));*/

    /*ncm_vector_log_vals (sing_up->Up_c, "COEFS: ", "% 22.15g", TRUE);*/
    /*ncm_matrix_log_vals (cov,           "COV:  ",  "% 22.15g");*/

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble m              = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  gdouble Ta                   = t / m;
  gdouble Tb                   = 1.0;
  gdouble lfac                 = 1.0;
  gdouble res                  = 0.0;
  gint j;

  for (j = 0; j < sing_up->chi_dim; j += 2)
  {
    res  += Ta * ncm_vector_get (sing_up->chi_c, j + 0);
    res  += Tb * ncm_vector_get (sing_up->chi_c, j + 1);
    Ta    = Ta * t / lfac;
    Tb    = Tb * t / lfac;
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble m              = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  gdouble Ta                   = t * t;
  gdouble Tb                   = m * t;
  gdouble lfac                 = 1.0;
  gdouble res                  = ncm_vector_get (sing_up->Up_c, 0);
  gint j;

  for (j = 1; j < sing_up->Up_dim; j += 2)
  {
    res  += Ta * ncm_vector_get (sing_up->Up_c, j + 0);
    res  += Tb * ncm_vector_get (sing_up->Up_c, j + 1);
    Ta    = Ta * t / lfac;
    Tb    = Tb * t / lfac;
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble m              = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  const gdouble dm             = ncm_csq1d_eval_dm (csq1d, model, t, self->k);
  gdouble Ta                   = 1.0 / m;
  gdouble Ta1                  = -t * dm / (m * m);
  gdouble Tb                   = 1.0 / t;
  gdouble lfac                 = 1.0;
  gdouble res                  = 0.0;
  gdouble n                    = 0.0;
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble m              = ncm_csq1d_eval_m (csq1d, model, t, self->k);
  const gdouble dm             = ncm_csq1d_eval_dm (csq1d, model, t, self->k);
  gdouble Ta                   = t;
  gdouble Tb                   = m;
  gdouble Tb1                  = dm * t;
  gdouble lfac                 = 1.0;
  gdouble res                  = 0.0;
  gdouble n                    = 0.0;
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  if (self->ti != ti)
  {
    self->ti            = ti;
    self->init_cond_set = FALSE;
  }
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  self->adiab_threshold = adiab_threshold;
}

/**
 * ncm_csq1d_set_prop_threshold:
 * @csq1d: a #NcmCSQ1D
 * @prop_threshold: mode $P_t$
 *
 * Sets the propagator threshold $P_t$.
 *
 */
void
ncm_csq1d_set_prop_threshold (NcmCSQ1D *csq1d, const gdouble prop_threshold)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  self->prop_threshold = prop_threshold;
}

/**
 * ncm_csq1d_set_save_evol:
 * @csq1d: a #NcmCSQ1D
 * @save: whether to save all evolution
 *
 * If true saves all evolution to be evaluted later through
 * ncm_csq1d_eval_at() and related methods.
 *
 */
void
ncm_csq1d_set_save_evol (NcmCSQ1D *csq1d, gboolean save_evol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  if (self->sing_detect != enable)
  {
    ncm_model_ctrl_force_update (self->ctrl);
    self->sing_detect = enable;
  }
}

/**
 * ncm_csq1d_set_init_cond:
 * @csq1d: a #NcmCSQ1D
 * @state: a #NcmCSQ1DEvolState
 * @ti: initial time $t_i$
 * @x: $\alpha$ or $\chi$ depending on the @state
 * @y: $\delta\gamma$, $U_+$ or $U_-$ depending on the @state
 *
 * Sets the values of the initial conditions at $t_i$.
 * This method also updates the value of $t_i$.
 *
 */
void
ncm_csq1d_set_init_cond (NcmCSQ1D *csq1d, NcmCSQ1DEvolState state, const gdouble ti, const gdouble x, const gdouble y)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (state)
  {
    case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
      NV_Ith_S (self->y, 0) = x;
      NV_Ith_S (self->y, 1) = y;
      break;
    case NCM_CSQ1D_EVOL_STATE_UP:
      NV_Ith_S (self->y_Up, 0) = x;
      NV_Ith_S (self->y_Up, 1) = y;
      break;
    case NCM_CSQ1D_EVOL_STATE_UM:
      NV_Ith_S (self->y_Um, 0) = x;
      NV_Ith_S (self->y_Um, 1) = y;
      break;
    default:
      g_error ("ncm_csq1d_set_init_cond: state %d not supported", state);
      break;
  }

  self->t     = ti;
  self->state = state;

  ncm_csq1d_set_ti (csq1d, ti);
  self->init_cond_set = TRUE;
}

/**
 * ncm_csq1d_set_init_cond_adiab:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @ti: initial time $t_i$
 *
 * Sets the values of the initial conditions at $t_i$.
 * This method also updates the value of $t_i$.
 *
 */
void
ncm_csq1d_set_init_cond_adiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  gdouble alpha, dgamma;

  ncm_csq1d_eval_adiab_at (csq1d, model, ti, &alpha, &dgamma, NULL, NULL);

  if ((fabs (dgamma) > self->adiab_threshold) || (fabs (alpha) > self->adiab_threshold))
    g_error ("ncm_csq1d_set_init_cond_adiab: time ti == % 22.15g is not a valid adiabatic time alpha, dgamma == (% 22.15g, % 22.15g)",
             ti, alpha, dgamma);
  else
    ncm_csq1d_set_init_cond (csq1d, NCM_CSQ1D_EVOL_STATE_ADIABATIC, ti, alpha, dgamma);
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->adiab_threshold;
}

/**
 * ncm_csq1d_get_prop_threshold:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the propagator threshold $P_t$.
 */
gdouble
ncm_csq1d_get_prop_threshold (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->prop_threshold;
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->sing_detect;
}

static gint _ncm_csq1d_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint _ncm_csq1d_f_Up (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Up (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint _ncm_csq1d_f_Um (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Um (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void
_ncm_csq1d_prepare_integrator (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
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
 *  flag = CVodeSetMaxOrd (self->cvode_Up, 2);
 *  NCM_CVODE_CHECK (&flag, "CVodeSetMaxOrd", 1, );
 *
 *  flag = CVodeSetMinStep (self->cvode_Up, 1.0e-5);
 *  NCM_CVODE_CHECK (&flag, "CVodeSetMinStep", 1, );
 */
}

static void
_ncm_csq1d_prepare_integrator_Um (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
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
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble alpha  = NV_Ith_S (y, 0);
  const gdouble dgamma = NV_Ith_S (y, 1);

  const gdouble nu    = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble F1    = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu = 2.0 * nu;

  NV_Ith_S (ydot, 0) = -twonu *sinh (dgamma);

  NV_Ith_S (ydot, 1) = +twonu * (-F1 + cosh (dgamma) * tanh (alpha));

  return 0;
}

static gint
_ncm_csq1d_f_Up (realtype t, N_Vector y_Up, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble chi = NV_Ith_S (y_Up, 0);
  const gdouble Up  = NV_Ith_S (y_Up, 1);

  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Up    = exp (Up);
  const gdouble ch2_alpha = 1.0 + chi * chi;

  NV_Ith_S (ydot, 0) = +m * nu2 * ch2_alpha / exp_Up - exp_Up / m;
  NV_Ith_S (ydot, 1) = +2.0 * m * nu2 * chi / exp_Up;

  return 0;
}

static gint
_ncm_csq1d_f_Um (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble chi = NV_Ith_S (y, 0);
  const gdouble Um  = NV_Ith_S (y, 1);

  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Um    = exp (Um);
  const gdouble ch2_alpha = 1.0 + chi * chi;

  NV_Ith_S (ydot, 0) = +m * nu2 * exp_Um - ch2_alpha / (m * exp_Um);
  NV_Ith_S (ydot, 1) = -2.0 * chi / (m * exp_Um);

  return 0;
}

static gint
_ncm_csq1d_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble alpha  = NV_Ith_S (y, 0);
  const gdouble dgamma = NV_Ith_S (y, 1);

  const gdouble nu    = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu = 2.0 * nu;

  /* - twonu * sinh (dgamma); */
  SM_ELEMENT_D (J, 0, 0) = 0.0;
  SM_ELEMENT_D (J, 0, 1) = -twonu *cosh (dgamma);

  /* + twonu * (-F1 + cosh (dgamma) * tanh (alpha)); */
  SM_ELEMENT_D (J, 1, 0) = +twonu *cosh (dgamma) / gsl_pow_2 (cosh (alpha));

  SM_ELEMENT_D (J, 1, 1) = +twonu *sinh (dgamma) * tanh (alpha);

  return 0;
}

static gint
_ncm_csq1d_J_Up (realtype t, N_Vector y_Up, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble chi = NV_Ith_S (y_Up, 0);
  const gdouble Up  = NV_Ith_S (y_Up, 1);

  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Up    = exp (Up);
  const gdouble ch2_alpha = 1.0 + chi * chi;

  /* + m * nu2 * ch2_alpha / exp_Up - exp_Up / m; */
  SM_ELEMENT_D (J, 0, 0) = +2.0 * m * nu2 * chi / exp_Up;
  SM_ELEMENT_D (J, 0, 1) = -m * nu2 * ch2_alpha / exp_Up - exp_Up / m;

  /* + 2.0 * m * nu2 * chi / exp_Up; */
  SM_ELEMENT_D (J, 1, 0) = +2.0 * m * nu2       / exp_Up;
  SM_ELEMENT_D (J, 1, 1) = -2.0 * m * nu2 * chi / exp_Up;

  return 0;
}

static gint
_ncm_csq1d_J_Um (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble chi = NV_Ith_S (y, 0);
  const gdouble Um  = NV_Ith_S (y, 1);

  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble exp_Um    = exp (Um);
  const gdouble ch2_alpha = 1.0 + chi * chi;

  /* + m * nu2 * exp_Um - ch2_alpha / (m * exp_Um); */
  SM_ELEMENT_D (J, 0, 0) = -2.0 * chi / (m * exp_Um);
  SM_ELEMENT_D (J, 0, 1) = +m * nu2 * exp_Um + ch2_alpha / (m * exp_Um);

  /* - 2.0 * chi / (m * exp_Um); */
  SM_ELEMENT_D (J, 1, 0) = -2.0       / (m * exp_Um);
  SM_ELEMENT_D (J, 1, 1) = +2.0 * chi / (m * exp_Um);

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
 * ncm_csq1d_eval_int_1_m: (virtual eval_int_1_m)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\int 1/m \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_int_mnu2: (virtual eval_int_mnu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\int m\nu^2 \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_int_qmnu2: (virtual eval_int_qmnu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\int \left(\int 1/m \mathrm{d}t\right) m\nu^2 \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_int_q2mnu2: (virtual eval_int_q2mnu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\int \left(\int 1/m \mathrm{d}t\right)^2 m\nu^2 \mathrm{d}t$.
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DEvolStop reason      = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t         = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_ADIABATIC);

  if (PRINT_EVOL)
    printf ("# ENTER EVOL ADIABATIC: TIME % 22.15g % 22.15g\n", self->t, self->tf);

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
      const gdouble Up         = +gamma + lnch_alpha;
      const gdouble Um         = -gamma + lnch_alpha;

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

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }

    if ((fabs (dgamma) > self->adiab_threshold) || (fabs (alpha) > self->adiab_threshold))
    {
      if (dgamma > 0.0)
        reason = NCM_CSQ1D_EVOL_STOP_UP_START;
      else
        reason = NCM_CSQ1D_EVOL_STOP_UM_START;

      break;
    }
  }

  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Up (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DEvolStop reason      = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t         = asinh (self->t);
  GArray *t_a                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *exp_Up_a             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *chim_t_a             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_UP);

  if (PRINT_EVOL)
    printf ("# ENTER EVOL UP: TIME % 22.15g => % 22.15g\n", self->t, self->tf);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble chim, chi, dgamma, gamma, Up, xi, m, exp_Up, chim_t;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode_Up, self->tf, self->y_Up, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_Up]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
    m       = ncm_csq1d_eval_m (csq1d, model, self->t, self->k);
    chi     = NV_Ith_S (self->y_Up, 0);
    Up      = NV_Ith_S (self->y_Up, 1);
    chim    = m * chi;
    gamma   = Up - 0.5 * log1p (chi * chi);
    dgamma  = gamma - xi;
    exp_Up  = exp (Up);
    chim_t  = chim / self->t;

    if (PRINT_EVOL)
    {
      const gdouble Um = -Up + log1p (chi * chi);

      printf ("# E[UP] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g\n", self->t, asinh (chi), dgamma, chi, Up, Um, xi);
    }

    is_finished = (self->t == self->tf);

    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      const gdouble alpha = asinh (chi);

      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }

    g_array_append_val (t_a,      self->t);
    g_array_append_val (exp_Up_a, exp_Up);
    g_array_append_val (chim_t_a, chim_t);

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }

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
  }

  g_array_unref (t_a);
  g_array_unref (exp_Up_a);
  g_array_unref (chim_t_a);

  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Um (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DEvolStop reason      = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t         = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_UM);

  if (PRINT_EVOL)
    printf ("# ENTER EVOL UM: TIME % 22.15g => % 22.15g\n", self->t, self->tf);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble chi, dgamma, gamma, Um, xi;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode_Um, self->tf, self->y_Um, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_Um]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
    chi     = NV_Ith_S (self->y_Um, 0);
    Um      = NV_Ith_S (self->y_Um, 1);
    gamma   = -Um + 0.5 * log1p (chi * chi);
    dgamma  = gamma - xi;

    if (PRINT_EVOL)
    {
      const gdouble Up = -Um + log1p (chi * chi);

      printf ("# E[UM] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g\n", self->t, asinh (chi), dgamma, chi, Up, Um, xi);
    }

    is_finished = (self->t == self->tf);

    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      const gdouble alpha = asinh (chi);

      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
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
  }

  return reason;
}

static void
_ncm_csq1d_evol_save (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DWS *ws, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DEvolStop stop;

  if (PRINT_EVOL)
    printf ("# ENTERING EVOL SAVE: STATE %d, TIME % 22.15g => % 22.15g\n", self->state, self->t, self->tf);

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
          const gdouble Up     = +gamma + gsl_sf_lncosh (alpha);

          NV_Ith_S (self->y_Up, 0) = sinh (alpha);
          NV_Ith_S (self->y_Up, 1) = Up;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UM:
        {
          const gdouble chi = NV_Ith_S (self->y_Um, 0);
          const gdouble Um  = NV_Ith_S (self->y_Um, 1);
          const gdouble Up  = -Um + log1p (chi * chi);

          NV_Ith_S (self->y_Up, 0) = chi;
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
          const gdouble Um     = -gamma + gsl_sf_lncosh (alpha);

          NV_Ith_S (self->y_Um, 0) = sinh (alpha);
          NV_Ith_S (self->y_Um, 1) = Um;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UP:
        {
          const gdouble chi = NV_Ith_S (self->y_Up, 0);
          const gdouble Up  = NV_Ith_S (self->y_Up, 1);
          const gdouble Um  = -Up + log1p (chi * chi);

          NV_Ith_S (self->y_Um, 0) = chi;
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
          const gdouble chi    = NV_Ith_S (self->y_Up, 0);
          const gdouble Up     = NV_Ith_S (self->y_Up, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = Up - 0.5 * log1p (chi * chi);
          const gdouble dgamma = gamma - xi;

          NV_Ith_S (self->y, 0) = asinh (chi);
          NV_Ith_S (self->y, 1) = dgamma;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UM:
        {
          const gdouble chi    = NV_Ith_S (self->y_Um, 0);
          const gdouble Um     = NV_Ith_S (self->y_Um, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = -Um + 0.5 * log1p (chi * chi);
          const gdouble dgamma = gamma - xi;

          NV_Ith_S (self->y, 0) = asinh (chi);
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
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  return ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
}

static gdouble
_ncm_csq1d_sing_detect (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws, NcmModel *model, gdouble ti, gdouble tf)
{
  /*NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);*/
  gdouble tsing = 0.5 * (tf + ti);
  guint iter = 0, max_iter = 1000;
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

  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    tsing  = gsl_root_fsolver_root (s);
    ti     = gsl_root_fsolver_x_lower (s);
    tf     = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (ti, tf, 0.0, root_reltol);

    if (PRINT_EVOL)
      printf ("#   %5d [%.7e, %.7e] %.7e %+.7e\n", iter, ti, tf, tsing, tf - ti);
  } while (status == GSL_CONTINUE && iter < max_iter);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0};

  if (NCM_CSQ1D_GET_CLASS (csq1d)->prepare != NULL)
    NCM_CSQ1D_GET_CLASS (csq1d)->prepare (csq1d, model);

  g_assert (self->init_cond_set);
  g_assert_cmpfloat (self->tf, >, self->ti);

  if (self->sing_detect)
    _ncm_csq1d_sing_detect (csq1d, &ws, model, self->ti, self->tf);

  if (self->save_evol)
  {
    GArray *asinh_t_a = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *alpha_a   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *dgamma_a  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmVector *asinh_t_v         = ncm_spline_get_xv (self->alpha_s);
  const guint len              = ncm_vector_len (asinh_t_v);
  GArray *t_a                  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);
  gdouble s_t                  = 1.0e300;
  guint i;

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
  /*NcmCSQ1DPrivate * const self = ws->ncm_csq1d_get_instance_private (csq1d);*/

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
  /*NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);*/
  gdouble alpha0, dgamma0, alpha_reltol0, dgamma_reltol0;
  gdouble alpha1, dgamma1, alpha_reltol1, dgamma_reltol1;

  gboolean adiab0, adiab1;

  g_assert_cmpfloat (reltol, >, 0.0);

  ncm_csq1d_eval_adiab_at (csq1d, model, t0, &alpha0, &dgamma0, &alpha_reltol0, &dgamma_reltol0);
  ncm_csq1d_eval_adiab_at (csq1d, model, t1, &alpha1, &dgamma1, &alpha_reltol1, &dgamma_reltol1);

  g_assert (gsl_finite (alpha0));
  g_assert (gsl_finite (alpha1));
  g_assert (gsl_finite (dgamma0));
  g_assert (gsl_finite (dgamma1));
  g_assert (gsl_finite (alpha_reltol0));
  g_assert (gsl_finite (alpha_reltol1));
  g_assert (gsl_finite (dgamma_reltol0));
  g_assert (gsl_finite (dgamma_reltol1));

  adiab0 = ((fabs (alpha_reltol0) < reltol) && (fabs (dgamma_reltol0) < reltol));
  adiab1 = ((fabs (alpha_reltol1) < reltol) && (fabs (dgamma_reltol1) < reltol));

  if ((adiab0 && adiab1) || (!adiab0 && !adiab1))
  {
    if (PRINT_EVOL)
      g_warning ("# Impossible to find the adiabatic limit: \n\tt0 % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n\tt1 % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n",
                 t0, alpha0, alpha_reltol0, dgamma0, dgamma_reltol0,
                 t1, alpha1, alpha_reltol1, dgamma1, dgamma_reltol1);

    return FALSE;
  }
  else
  {
    NcmCSQ1DWS ws = {csq1d, model, reltol};
    guint iter = 0, max_iter = 1000;
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

    do {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      ti[0]  = gsl_root_fsolver_root (s);
      t0     = gsl_root_fsolver_x_lower (s);
      t1     = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (t0, t1, 0.0, root_reltol);

      if (PRINT_EVOL)
        printf ("#  %5d [%.7e, %.7e] %.7e %+.7e\n", iter, t0, t1, ti[0], t1 - t0);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
  }

  return TRUE;
}

static gdouble
_ncm_csq1d_F1_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, t, self->k);

  return F1;
}

static gdouble
_ncm_csq1d_F2_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F2             = ncm_csq1d_eval_F2 (ws->csq1d, ws->model, t, self->k);

  return F2;
}

static gdouble
_ncm_csq1d_lnnu_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble nu             = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);

  return log (nu);
}

static gdouble _ncm_csq1d_abs_F1_logt (gdouble at, gpointer user_data);
static gdouble _ncm_csq1d_ln_abs_F1_eps_logt (gdouble at, gpointer user_data);

/**
 * ncm_csq1d_find_adiab_max:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t0: time lower bound $t_0$
 * @t1: time upper bound $t_1$
 * @border_eps: border epsilon $\epsilon$
 * @F1_min: (out): the value of $F_1(t_\mathrm{min})$
 * @t_Bl: (out): the value of $t_{B,\mathrm{lower}}$
 * @t_Bu: (out): the value of $t_{B,\mathrm{upper}}$
 *
 * Computes the time $t_\mathrm{min}$ that minimizes $F_1(t)$.
 *
 * Returns: the time $t_\mathrm{min}$.
 */
gdouble
ncm_csq1d_find_adiab_max (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble border_eps, gdouble *F1_min, gdouble *t_Bl, gdouble *t_Bu)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, border_eps};

  gsl_min_fminimizer *fmin = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  const gdouble atl        = log (t0);
  const gdouble atu        = log (t1);
  gdouble at0              = atl;
  gdouble at1              = atu;
  gdouble atm              = (at0 + at1) * 0.5;
  guint iter               = 0;
  gint status;
  gsl_function F;
  gint ret;

  F.params   = &ws;
  F.function = &_ncm_csq1d_abs_F1_logt;

  ret = gsl_min_fminimizer_set (fmin, &F, atm, at0, at1);
  NCM_TEST_GSL_RESULT ("ncm_csq1d_find_adiab_max", ret);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate (fmin);

    if (status)
      g_error ("ncm_csq1d_find_adiab_max: Cannot find minimum (%s)", gsl_strerror (status));

    atm = gsl_min_fminimizer_x_minimum (fmin);
    at0 = gsl_min_fminimizer_x_lower (fmin);
    at1 = gsl_min_fminimizer_x_upper (fmin);

    status = gsl_min_test_interval (at0, at1, 0.0, self->reltol);

    /*ncm_message ("[%d] % 22.15e % 22.15e % 22.15e\n", status, exp (atm), exp (at0), exp (at1));*/
  } while (status == GSL_CONTINUE && iter < 1000);

  gsl_min_fminimizer_free (fmin);

  {
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    guint max_iter = 1000;

    iter = 0;

    F.function = &_ncm_csq1d_ln_abs_F1_eps_logt;
    F.params   = &ws;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);

    gsl_root_fsolver_set (s, &F, atl, atm);

    do {
      iter++;
      status  = gsl_root_fsolver_iterate (s);
      t_Bl[0] = gsl_root_fsolver_root (s);
      at0     = gsl_root_fsolver_x_lower (s);
      at1     = gsl_root_fsolver_x_upper (s);
      status  = gsl_root_test_interval (at0, at1, 0.0, 1.0e-7);

      /*ncm_message ("Bl: [%d] % 22.15e % 22.15e % 22.15e\n", status, exp (t_Bl[0]), exp (at0), exp (at1));*/
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_set (s, &F, atm, atu);

    do {
      iter++;
      status  = gsl_root_fsolver_iterate (s);
      t_Bu[0] = gsl_root_fsolver_root (s);
      at0     = gsl_root_fsolver_x_lower (s);
      at1     = gsl_root_fsolver_x_upper (s);
      status  = gsl_root_test_interval (at0, at1, 0.0, 1.0e-7);

      /*ncm_message ("Bu: [%d] % 22.15e % 22.15e % 22.15e\n", status, exp (t_Bu[0]), exp (at0), exp (at1));*/
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
  }

  {
    const gdouble tm = exp (atm);

    F1_min[0] = ncm_csq1d_eval_F1 (csq1d, model, tm, self->k);
    t_Bl[0]   = exp (t_Bl[0]);
    t_Bu[0]   = exp (t_Bu[0]);

    return tm;
  }
}

static gdouble
_ncm_csq1d_abs_F1_logt (gdouble at, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, exp (at), self->k);

  return fabs (F1);
}

static gdouble
_ncm_csq1d_ln_abs_F1_eps_logt (gdouble at, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, exp (at), self->k);

  return fabs (F1 / ws->reltol) - 1.0;
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0};
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t, self->k);
  const gdouble F2             = ncm_csq1d_eval_F2 (csq1d, model, t, self->k);
  const gdouble F1_2           = F1 * F1;
  const gdouble F1_3           = F1_2 * F1;
  const gdouble nu             = ncm_csq1d_eval_nu (csq1d, model, t, self->k);
  const gdouble twonu          = 2.0 * nu;
  gdouble err, F3, d2F2, dlnnu, F4, alpha_reltol0, dgamma_reltol0;

  F3             = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err) / twonu;
  d2F2           = ncm_diff_rc_d2_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err);
  dlnnu          = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_lnnu_func, &ws, &err);
  F4             = d2F2 / gsl_pow_2 (twonu) - dlnnu * F3 / twonu;
  alpha_reltol0  = gsl_pow_2 ((F1_3 / 3.0 - F3) / F1);
  dgamma_reltol0 = gsl_pow_2 ((F4 - F1_2 * F2) / F2);

  alpha[0]  = +F1 + F1_3 / 3.0 - F3;
  dgamma[0] = -(1.0 + F1_2) * F2 + F4;

  if (alpha_reltol != NULL)
    alpha_reltol[0] = alpha_reltol0;

  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = dgamma_reltol0;
}

/**
 * ncm_csq1d_eval_nonadiab_at:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @nonadiab_frame: frame
 * @t: time $t$
 * @chi: (out): value of $\chi(t)$
 * @Up: (out): value of $U_+$
 *
 * Computes the value of the non-adiabatic VDC order order two
 * in the variables $\chi$ and $U$ at $t$. The VDC can be computed
 * at different frames by choosing @nonadiab_frame.
 *
 */
void
ncm_csq1d_eval_nonadiab_at (NcmCSQ1D *csq1d, NcmModel *model, guint nonadiab_frame, const gdouble t, gdouble *chi, gdouble *Up)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble q0             = ncm_csq1d_eval_int_1_m (csq1d, model, t, self->k);
  const gdouble q1             = ncm_csq1d_eval_int_q2mnu2 (csq1d, model, t, self->k);
  const gdouble p1             = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t, self->k);
  const gdouble r1             = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t, self->k);

  switch (nonadiab_frame)
  {
    case 0:
      chi[0] = (2.0 * p1 - 1.0) * (q0 + q1) + 2.0 * r1;
      Up[0]  = -2.0 * p1;
      break;
    case 1:
      chi[0] = +2.0 * r1;
      Up[0]  = -2.0 * p1;
      break;
    case 2:
      chi[0] = 0.0;
      Up[0]  = 0.0;
      break;
  }
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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);

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
  NcmCSQ1DPrivate * const self      = ncm_csq1d_get_instance_private (csq1d);
  const gdouble gamma               = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;
  const gdouble exp_gamma_p_alpha_2 = exp (0.5 * (gamma + alpha));
  const gdouble exp_gamma_m_alpha_2 = exp (0.5 * (gamma - alpha));

  phi[0] = +0.5 / exp_gamma_m_alpha_2;
  phi[1] = -0.5 / exp_gamma_p_alpha_2;

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
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);

  const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
  const gdouble gamma  = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  J11[0] = cosh (alpha) * exp (-gamma);
  J22[0] = cosh (alpha) * exp (+gamma);
  J12[0] = -sinh (alpha);
}

/**
 * ncm_csq1d_get_H_poincare_hp:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_H_poincare_hp (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  x[0]   = 0.0;
  lny[0] = -ncm_csq1d_eval_xi (csq1d, model, t, self->k);
}

/**
 * ncm_csq1d_get_poincare_hp:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_poincare_hp (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);

  const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
  const gdouble gamma  = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  x[0]   = exp (-gamma) * tanh (alpha);
  lny[0] = -gamma - gsl_sf_lncosh (alpha);
}

#define _ALPHA_TO_THETA(alpha) (atan (sinh (alpha)))
#define _THETA_TO_ALPHA(theta) (asinh (tan (theta)))

/*
 *  static void
 *  _ncm_csq1d_ct_g0g2 (const gdouble alpha, const gdouble gamma, const gdouble p, gdouble *alphap, gdouble *gammap)
 *  {
 *  alphap[0] = alpha;
 *  gammap[0] = gamma - p;
 *  }
 */

static void
_ncm_csq1d_ct_g0g1 (const gdouble alpha, const gdouble gamma, const gdouble p, gdouble *alphap, gdouble *gammap)
{
  const gdouble theta = _ALPHA_TO_THETA (alpha);
  gdouble thetap      = 0.0;

  ncm_util_mln_1mIexpzA_1pIexpmzA (gamma, theta, tanh (0.5 * p), gammap, &thetap);
  alphap[0] = _THETA_TO_ALPHA (thetap);
}

/**
 * ncm_csq1d_get_poincare_hp_frame:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @adiab_frame: which adiabatic frame to use
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_poincare_hp_frame (NcmCSQ1D *csq1d, NcmModel *model, guint adiab_frame, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);

  switch (adiab_frame)
  {
    case 0:
    {
      const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
      const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
      const gdouble gamma  = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

      x[0]   = exp (-gamma) * tanh (alpha);
      lny[0] = -gamma - gsl_sf_lncosh (alpha);
      break;
    }
    case 1:
    {
      const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
      const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);

      x[0]   = exp (-dgamma) * tanh (alpha);
      lny[0] = -dgamma - gsl_sf_lncosh (alpha);
      break;
    }
    case 2:
    {
      const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
      const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
      const gdouble theta  = _ALPHA_TO_THETA (alpha);
      const gdouble F1     = ncm_csq1d_eval_F1 (csq1d, model, t, self->k);
      const gdouble xi1    = atanh (-F1);

      gdouble thetap = 0.0;
      gdouble alphap = 0.0;
      gdouble gammap = 0.0;

      ncm_util_mln_1mIexpzA_1pIexpmzA (dgamma, theta, tanh (0.5 * xi1), &gammap, &thetap);

      alphap = _THETA_TO_ALPHA (thetap);

      x[0]   = exp (-gammap) * tanh (alphap);
      lny[0] = -gammap - gsl_sf_lncosh (alphap);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_csq1d_alpha_dgamma_to_minkowski_frame:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @adiab_frame: which adiabatic frame to use
 * @t: time $t$
 * @alpha: $\alpha$
 * @dgamma: $\delta\gamma$
 * @x1: (out): $x_1$
 * @x2: (out): $x_2$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_alpha_dgamma_to_minkowski_frame (NcmCSQ1D *csq1d, NcmModel *model, guint adiab_frame, const gdouble t, const gdouble alpha, const gdouble dgamma, gdouble *x1, gdouble *x2)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (adiab_frame)
  {
    case 0:
    {
      const gdouble gamma = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

      x1[0] = +sinh (alpha);
      x2[0] = -cosh (alpha) * sinh (gamma);
      break;
    }
    case 1:
    {
      x1[0] = +sinh (alpha);
      x2[0] = -cosh (alpha) * sinh (dgamma);
      break;
    }
    case 2:
    {
      const gdouble theta = _ALPHA_TO_THETA (alpha);
      const gdouble F1    = ncm_csq1d_eval_F1 (csq1d, model, t, self->k);
      const gdouble xi1   = atanh (-F1);

      gdouble thetap = 0.0;
      gdouble alphap = 0.0;
      gdouble gammap = 0.0;

      ncm_util_mln_1mIexpzA_1pIexpmzA (dgamma, theta, tanh (0.5 * xi1), &gammap, &thetap);

      alphap = _THETA_TO_ALPHA (thetap);

      x1[0] = +sinh (alphap);
      x2[0] = -cosh (alphap) * sinh (gammap);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_csq1d_get_minkowski_frame:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @adiab_frame: which adiabatic frame to use
 * @t: time $t$
 * @x1: (out): $x_1$
 * @x2: (out): $x_2$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_minkowski_frame (NcmCSQ1D *csq1d, NcmModel *model, guint adiab_frame, const gdouble t, gdouble *x1, gdouble *x2)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);
  const gdouble alpha          = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble dgamma         = ncm_spline_eval (self->dgamma_s, a_t);

  ncm_csq1d_alpha_dgamma_to_minkowski_frame (csq1d, model, adiab_frame, t, alpha, dgamma, x1, x2);
}

/**
 * ncm_csq1d_get_Hadiab_poincare_hp:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_Hadiab_poincare_hp (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  gdouble alpha, dgamma, gamma;

  ncm_csq1d_eval_adiab_at (csq1d, model, t, &alpha, &dgamma, NULL, NULL);

  gamma = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  x[0]   = exp (-gamma) * tanh (alpha);
  lny[0] = -gamma - gsl_sf_lncosh (alpha);
}

/**
 * ncm_csq1d_get_H_poincare_disc:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_H_poincare_disc (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble xi             = ncm_csq1d_eval_xi (csq1d, model, t, self->k);

  x[0]   = 0.0;
  lny[0] = -tanh (xi / 2);
}

/**
 * ncm_csq1d_get_Hadiab_poincare_disc:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_Hadiab_poincare_disc (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  gdouble alpha, dgamma, gamma;

  ncm_csq1d_eval_adiab_at (csq1d, model, t, &alpha, &dgamma, NULL, NULL);

  gamma = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  x[0]   = +sinh (alpha) / (1.0 + cosh (alpha) * cosh (gamma));
  lny[0] = -sinh (gamma) * cosh (alpha) / (1.0 + cosh (alpha) * cosh (gamma));
}

/**
 * ncm_csq1d_get_poincare_disc:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @x: (out): $x$
 * @lny: (out): $\ln(y)$
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_get_poincare_disc (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);

  const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
  const gdouble gamma  = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  x[0]   = +sinh (alpha) / (1.0 + cosh (alpha) * cosh (gamma));
  lny[0] = -sinh (gamma) * cosh (alpha) / (1.0 + cosh (alpha) * cosh (gamma));
}

static gint _ncm_csq1d_f_Prop (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Prop (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static gint _ncm_csq1d_root_Prop (realtype lambda, N_Vector y, realtype *gout, gpointer user_data);

static void
_ncm_csq1d_prepare_integrator_Prop (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws, const gdouble ti, const gdouble tf)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  gint flag;

  if (!self->cvode_Prop_init)
  {
    self->cvode_Prop = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode_Prop, &_ncm_csq1d_f_Prop, ti, self->y_Prop);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_Prop_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode_Prop, ti, self->y_Prop);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode_Prop, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode_Prop, 100000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode_Prop, self->LS_Prop, self->A_Prop);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode_Prop, &_ncm_csq1d_J_Prop);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode_Prop, fabs (ti) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

  flag = CVodeSetUserData (self->cvode_Prop, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode_Prop, tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  /*flag = CVodeSetMaxStep (self->cvode_Prop, (self->tf - self->t) / 100.0);*/
  /*NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );*/

  flag = CVodeRootInit (self->cvode_Prop, 3, &_ncm_csq1d_root_Prop);
  NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );
}

static gint
_ncm_csq1d_root_Prop (realtype t, N_Vector y, realtype *gout, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  /*NcmCSQ1DPrivate * const self = ws->ncm_csq1d_get_instance_private (csq1d);*/
  const gdouble a_2  = NV_Ith_S (y, 0) * NV_Ith_S (y, 0);
  const gdouble u1_2 = -NV_Ith_S (y, 1) * NV_Ith_S (y, 2) - NV_Ith_S (y, 3) * NV_Ith_S (y, 3);

  gout[0] = u1_2 - ws->reltol;
  gout[1] = fabs ((a_2 + u1_2) - 1.0) - 1.0e-1;
  gout[2] = fabs ((a_2 + u1_2) - 1.0) - 5.0e-1;

  return 0;
}

static gint
_ncm_csq1d_f_Prop (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble a = NV_Ith_S (y, 0);
  const gdouble b = NV_Ith_S (y, 1);
  const gdouble c = NV_Ith_S (y, 2);
  const gdouble h = NV_Ith_S (y, 3);

  const gdouble m       = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2     = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble mnu2    = m * nu2;
  const gdouble mnu2_2  = 0.5 * mnu2;
  const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t, self->k);
  const gdouble q       = -int_1_m;
  const gdouble q2      = q * q;

  if (G_UNLIKELY (mnu2 == 0.0))
  {
    NV_Ith_S (ydot, 0) = 0.0;
    NV_Ith_S (ydot, 1) = 0.0;
    NV_Ith_S (ydot, 2) = 0.0;
    NV_Ith_S (ydot, 3) = 0.0;

    return 0;
  }

  NV_Ith_S (ydot, 0) = -mnu2_2 * (+b - q2 * c + 2.0 * q * h);
  NV_Ith_S (ydot, 1) = -mnu2   * (+q * b - q2 * (a - h));
  NV_Ith_S (ydot, 2) = -mnu2   * (-q * c +      (a + h));
  NV_Ith_S (ydot, 3) = -mnu2_2 * (-b - q2 * c + 2.0 * q * a);

  return 0;
}

static gint
_ncm_csq1d_J_Prop (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble m       = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2     = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble mnu2    = m * nu2;
  const gdouble mnu2_2  = 0.5 * mnu2;
  const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t, self->k);
  const gdouble q       = -int_1_m;
  const gdouble q2      = q * q;

  if (G_UNLIKELY (mnu2 == 0.0))
    return 0;

  /* -mnu2_2 * (+b - q2 * c + 2.0 * q * h); */
  SM_ELEMENT_D (J, 0, 0) = 0.0;
  SM_ELEMENT_D (J, 0, 1) = -mnu2_2;
  SM_ELEMENT_D (J, 0, 2) = -mnu2_2 * q2;
  SM_ELEMENT_D (J, 0, 3) = -mnu2 * q;

  /* -mnu2   * (+q * b - q2 * (a - h)); */
  SM_ELEMENT_D (J, 1, 0) = +mnu2 * q2;
  SM_ELEMENT_D (J, 1, 1) = -mnu2 * q;
  SM_ELEMENT_D (J, 1, 2) = 0.0;
  SM_ELEMENT_D (J, 1, 3) = -mnu2 * q2;

  /* -mnu2   * (-q * c +      (a + h)); */
  SM_ELEMENT_D (J, 2, 0) = -mnu2;
  SM_ELEMENT_D (J, 2, 1) = 0.0;
  SM_ELEMENT_D (J, 2, 2) = +mnu2 * q;
  SM_ELEMENT_D (J, 2, 3) = -mnu2;

  /* -mnu2_2 * (-b - q2 * c + 2.0 * q * a); */
  SM_ELEMENT_D (J, 3, 0) = -mnu2 * q;
  SM_ELEMENT_D (J, 3, 1) = +mnu2_2;
  SM_ELEMENT_D (J, 3, 2) = +mnu2_2 * q2;
  SM_ELEMENT_D (J, 3, 3) = 0.0;

  return 0;
}

static gdouble _ncm_csq1d_prepare_prop_mnu2 (gdouble t, gpointer params);
static gdouble _ncm_csq1d_prepare_prop_qmnu2 (gdouble t, gpointer params);
static gdouble _ncm_csq1d_prepare_prop_q2mnu2 (gdouble t, gpointer params);

static void
_ncm_csq1d_prepare_prop_eval_u1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti, const gdouble t, gdouble u1[3])
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0};

  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gsl_function F;
  gdouble int_mnu2, int_qmnu2, int_q2mnu2, err;

  F.function = &_ncm_csq1d_prepare_prop_mnu2;
  F.params   = &ws;

  gsl_integration_qag (&F, ti, t, 0.0, self->reltol, NCM_INTEGRAL_PARTITION, 6, *w, &int_mnu2, &err);

  F.function = &_ncm_csq1d_prepare_prop_qmnu2;
  F.params   = &ws;

  gsl_integration_qag (&F, ti, t, 0.0, self->reltol, NCM_INTEGRAL_PARTITION, 6, *w, &int_qmnu2, &err);

  F.function = &_ncm_csq1d_prepare_prop_q2mnu2;
  F.params   = &ws;

  gsl_integration_qag (&F, ti, t, 0.0, self->reltol, NCM_INTEGRAL_PARTITION, 6, *w, &int_q2mnu2, &err);

  u1[0] = -int_qmnu2;
  u1[1] = +int_q2mnu2;
  u1[2] = -int_mnu2;

  ncm_memory_pool_return (w);
}

/**
 * ncm_csq1d_prepare_prop:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @ti: initial time $t_i$
 * @tii: integral approximation time $t_{\mathrm{i}i}$
 * @tf: max time $t_f$
 *
 *
 * Computes the complex structure matrix.
 *
 */
void
ncm_csq1d_prepare_prop (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti, const gdouble tii, const gdouble tf)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, self->prop_threshold};
  GArray *t_a                  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  gdouble t                    = 0.0;
  gboolean tf_Prop_found       = FALSE;
  GArray *R_a[4];
  gdouble u1[3];
  gint i;

  _ncm_csq1d_prepare_prop_eval_u1 (csq1d, model, ti, tii, u1);

/*
 *  ncm_message ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", ti, tii,
 *     u1[0], +ncm_csq1d_eval_int_qmnu2 (csq1d, model, tii, self->k),
 *     u1[1], +ncm_csq1d_eval_int_q2mnu2 (csq1d, model, tii, self->k),
 *     u1[2], -ncm_csq1d_eval_int_mnu2 (csq1d, model, tii, self->k)
 *     );
 */

  NV_Ith_S (self->y_Prop, 0) = 1.0;
  NV_Ith_S (self->y_Prop, 1) = u1[1];
  NV_Ith_S (self->y_Prop, 2) = u1[2];
  NV_Ith_S (self->y_Prop, 3) = u1[0];

  g_array_append_val (t_a, ti);

  for (i = 0; i < 4; i++)
  {
    R_a[i] = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    g_array_append_val (R_a[i], NV_Ith_S (self->y_Prop, i));
  }

  _ncm_csq1d_prepare_integrator_Prop (csq1d, &ws, tii, tf);

  self->tf_Prop = GSL_NAN;

  while (TRUE)
  {
    gint flag = CVode (self->cvode_Prop, tf, self->y_Prop, &t, CV_ONE_STEP);

    NCM_CVODE_CHECK (&flag, "CVode[ncm_csq1d_prepare_prop]", 1, );

    g_array_append_val (t_a, t);

    for (i = 0; i < 4; i++)
      g_array_append_val (R_a[i], NV_Ith_S (self->y_Prop, i));

    if (flag == CV_ROOT_RETURN)
    {
      gint roots[4];

      flag = CVodeGetRootInfo (self->cvode_Prop, roots);
      NCM_CVODE_CHECK (&flag, "CVodeGetRootInfo[ncm_csq1d_prepare_prop]", 1, );

      /*ncm_message ("% 22.15g %d %d %d\n", t, roots[0], roots[1], roots[2]);*/
      if (roots[0] && !tf_Prop_found)
      {
        self->tf_Prop = t;
        tf_Prop_found = TRUE;
      }
      else if ((roots[1] != 0) || (roots[2] != 0))
      {
        break;
      }
    }

    if (t >= tf)
      break;
  }

  for (i = 0; i < 4; i++)
  {
    ncm_spline_set_array (self->R[i], t_a, R_a[i], TRUE);
    g_array_unref (R_a[i]);
  }

  g_array_unref (t_a);
}

static gdouble
_ncm_csq1d_prepare_prop_mnu2 (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble m   = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2 = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);

  return m * nu2;
}

static gdouble
_ncm_csq1d_prepare_prop_qmnu2 (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble m    = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2  = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble mnu2 = m * nu2;

  if (mnu2 == 0)
  {
    return 0.0;
  }
  else
  {
    const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t, self->k);
    const gdouble q       = -int_1_m;

    return q * mnu2;
  }
}

static gdouble
_ncm_csq1d_prepare_prop_q2mnu2 (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble m    = ncm_csq1d_eval_m (ws->csq1d, ws->model, t, self->k);
  const gdouble nu2  = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t, self->k);
  const gdouble mnu2 = m * nu2;

  if (mnu2 == 0)
  {
    return 0.0;
  }
  else
  {
    const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t, self->k);
    const gdouble q       = -int_1_m;

    return q * q * mnu2;
  }
}

/**
 * ncm_csq1d_get_tf_prop:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: current final time $t_f$ for the propagator.
 */
gdouble
ncm_csq1d_get_tf_prop (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->tf_Prop;
}

/**
 * ncm_csq1d_get_prop_vector_chi_Up:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @chi: (out): $\chi$
 * @Up: (out): $U_+$
 *
 */
void
ncm_csq1d_get_prop_vector_chi_Up (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *chi, gdouble *Up)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble b              = ncm_spline_eval (self->R[1], t);
  const gdouble c              = ncm_spline_eval (self->R[2], t);
  const gdouble h              = ncm_spline_eval (self->R[3], t);
  const gdouble n0             = sqrt (-(b * c + h * h));

  chi[0] = h / n0;
  Up[0]  = log (-c / n0);
}

/**
 * ncm_csq1d_evolve_prop_vector_chi_Up:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @nonadiab_frame: frame
 * @chi_i: $\chi_i$
 * @Up_i: $U_{+i}$
 * @chi: (out): $\chi$
 * @Up: (out): $U_+$
 *
 *
 */
void
ncm_csq1d_evolve_prop_vector_chi_Up (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const guint nonadiab_frame, gdouble chi_i, gdouble Up_i, gdouble *chi, gdouble *Up)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble int_1_m        = ncm_csq1d_eval_int_1_m (csq1d, model, t, self->k);
  const gdouble q0             = int_1_m;
  const gdouble q1             = ncm_csq1d_eval_int_q2mnu2 (csq1d, model, t, self->k);
  const gdouble a              = ncm_spline_eval (self->R[0], t);
  const gdouble b              = ncm_spline_eval (self->R[1], t);
  const gdouble c              = ncm_spline_eval (self->R[2], t);
  const gdouble h              = ncm_spline_eval (self->R[3], t);
  const gdouble a11            = a + h;
  const gdouble a22            = a - h;
  const gdouble a12            = b;
  const gdouble a21            = c;
  const gdouble onepchi2       = 1.0 + chi_i * chi_i;
  const gdouble exp_Up_i       = exp (Up_i);

  chi[0] = a12 * a21 * chi_i + a22 * (a11 * chi_i - a12 * exp_Up_i) - a11 * a21 * onepchi2 / exp_Up_i;
  Up[0]  = -(2.0 * a21 * a22 * chi_i - a22 * a22 * exp_Up_i - a21 * a21 * onepchi2 / exp_Up_i);

  switch (nonadiab_frame)
  {
    case 0:
      chi[0] = (chi[0] - Up[0] * q0);
      Up[0]  = log (Up[0]);
      break;
    case 1:
      /* q1 L2p */
      chi[0] = (chi[0] + Up[0] * q1);
      Up[0]  = log (Up[0]);
      break;
    case 2:
    {
      const gdouble p1 = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t, self->k);

      /* q1 L2p */
      chi[0] = (chi[0] + Up[0] * q1);
      Up[0]  = log (Up[0]);

      /* p1 g1 */
      Up[0] = Up[0] + 2.0 * p1;

      /* r1 g2 */
      {
        const gdouble alpha = asinh (chi[0]);
        const gdouble gamma = Up[0] - 0.5 * log1p (chi[0] * chi[0]);
        const gdouble r1    = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t, self->k);

        gdouble alphap, gammap;

        _ncm_csq1d_ct_g0g1 (alpha, gamma, -2.0 * r1, &alphap, &gammap);

        chi[0] = sinh (alphap);
        Up[0]  = gammap + 0.5 * log1p (chi[0] * chi[0]);
      }
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_csq1d_alpha_gamma_circle:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @alpha: $\alpha$ variable
 * @gamma: $\gamma$ variable
 * @r: radius
 * @theta: angle
 * @alphap: (out): output $\alpha$ variable
 * @gammap: (out): output $\gamma$ variable
 *
 * Computes the complex structure matrix parameters for a circle
 * around the point $(\alpha, \gamma)$ with radius $r$ and angle
 * $\theta$.
 *
 */
void
ncm_csq1d_alpha_gamma_circle (NcmCSQ1D *csq1d, NcmModel *model, const gdouble alpha, const gdouble gamma, const gdouble r, const gdouble theta, gdouble *alphap, gdouble *gammap)
{
  const gdouble ct = cos (theta);
  const gdouble st = sin (theta);
  const gdouble tr = tanh (r);
  const gdouble ca = cosh (alpha);
  const gdouble sa = sinh (alpha);
  const gdouble t1 = cosh (r) * sa + ct * ca * sinh (r);
  const gdouble t2 = -(2.0 * st * tr / (ca + tr * (st + ct * sa)));

  alphap[0] = asinh (t1);
  gammap[0] = gamma + 0.5 * log1p (t2);
}

