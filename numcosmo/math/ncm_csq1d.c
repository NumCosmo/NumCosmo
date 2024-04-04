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
#include "ncm_enum_types.h"

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
  gdouble ti;
  gdouble tf;
  gdouble t;
  gboolean init_cond_set;
  NcmCSQ1DEvolState state;
  gdouble adiab_threshold;
  gdouble prop_threshold;
  gboolean save_evol;
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
  NcmSpline *gamma_s;
  NcmDiff *diff;
  NcmSpline *R[4];
  gdouble tf_Prop;
  gdouble ti_Prop;
  NcmCSQ1DState *cur_state;
  NcmCSQ1DVacuumType initial_condition_type;
  gdouble vacuum_reltol;
  gdouble vacuum_max_time;
  gdouble vacuum_final_time;
} NcmCSQ1DPrivate;

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_TI,
  PROP_TF,
  PROP_ADIAB_THRESHOLD,
  PROP_PROP_THRESHOLD,
  PROP_SAVE_EVOL,
  PROP_VACUUM_TYPE,
  PROP_VACUUM_RELTOL,
  PROP_VACUUM_MAX_TIME,
};

G_DEFINE_BOXED_TYPE (NcmCSQ1DState, ncm_csq1d_state, ncm_csq1d_state_copy, ncm_csq1d_state_free)
G_DEFINE_TYPE_WITH_PRIVATE (NcmCSQ1D, ncm_csq1d, G_TYPE_OBJECT)

static void
ncm_csq1d_init (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  self->reltol          = 0.0;
  self->abstol          = 0.0;
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
  self->gamma_s  = ncm_spline_cubic_notaknot_new ();

  {
    gint i;

    for (i = 0; i < 4; i++)
      self->R[i] = ncm_spline_cubic_notaknot_new ();
  }
  self->tf_Prop = 0.0;
  self->ti_Prop = 0.0;

  self->diff      = ncm_diff_new ();
  self->cur_state = ncm_csq1d_state_new ();

  self->initial_condition_type = NCM_CSQ1D_VACUUM_TYPE_LENGTH;
  self->vacuum_reltol          = 0.0;
  self->vacuum_max_time        = 0.0;
  self->vacuum_final_time      = 0.0;
}

static void
_ncm_csq1d_dispose (GObject *object)
{
  NcmCSQ1D *csq1d              = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  ncm_model_ctrl_clear (&self->ctrl);
  ncm_spline_clear (&self->alpha_s);
  ncm_spline_clear (&self->dgamma_s);
  ncm_spline_clear (&self->gamma_s);

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

  ncm_csq1d_state_free (self->cur_state);

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
    case PROP_VACUUM_TYPE:
      ncm_csq1d_set_initial_condition_type (csq1d, g_value_get_enum (value));
      break;
    case PROP_VACUUM_RELTOL:
      ncm_csq1d_set_vacuum_reltol (csq1d, g_value_get_double (value));
      break;
    case PROP_VACUUM_MAX_TIME:
      ncm_csq1d_set_vacuum_max_time (csq1d, g_value_get_double (value));
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
    case PROP_VACUUM_TYPE:
      g_value_set_enum (value, ncm_csq1d_get_initial_condition_type (csq1d));
      break;
    case PROP_VACUUM_RELTOL:
      g_value_set_double (value, ncm_csq1d_get_vacuum_reltol (csq1d));
      break;
    case PROP_VACUUM_MAX_TIME:
      g_value_set_double (value, ncm_csq1d_get_vacuum_max_time (csq1d));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _ncm_csq1d_eval_xi         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_nu         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_nu2        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_m          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_int_1_m    (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_int_mnu2   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_int_qmnu2  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_F1         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _ncm_csq1d_eval_F2         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);

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
                                                        0.0, G_MAXDOUBLE, 1.0e0,
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
                                   PROP_VACUUM_TYPE,
                                   g_param_spec_enum ("vacuum-type",
                                                      NULL,
                                                      "The vacuum type",
                                                      NCM_TYPE_CSQ1D_VACUUM_TYPE,
                                                      NCM_CSQ1D_VACUUM_TYPE_ADIABATIC4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_VACUUM_RELTOL,
                                   g_param_spec_double ("vacuum-reltol",
                                                        NULL,
                                                        "The vacuum relative tolerance",
                                                        DBL_EPSILON, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_VACUUM_MAX_TIME,
                                   g_param_spec_double ("vacuum-max-time",
                                                        NULL,
                                                        "The vacuum maximum time",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->eval_xi         = &_ncm_csq1d_eval_xi;
  klass->eval_nu         = &_ncm_csq1d_eval_nu;
  klass->eval_nu2        = &_ncm_csq1d_eval_nu2;
  klass->eval_m          = &_ncm_csq1d_eval_m;
  klass->eval_int_1_m    = &_ncm_csq1d_eval_int_1_m;
  klass->eval_int_mnu2   = &_ncm_csq1d_eval_int_mnu2;
  klass->eval_int_qmnu2  = &_ncm_csq1d_eval_int_qmnu2;
  klass->eval_int_q2mnu2 = &_ncm_csq1d_eval_int_q2mnu2;
  klass->eval_F1         = &_ncm_csq1d_eval_F1;
  klass->eval_F2         = &_ncm_csq1d_eval_F2;
}

static gdouble
_ncm_csq1d_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_xi: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_nu: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return gsl_pow_2 (NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t));
}

static gdouble
_ncm_csq1d_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return exp (NCM_CSQ1D_GET_CLASS (csq1d)->eval_xi (csq1d, model, t)) / NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t);
}

static gdouble
_ncm_csq1d_eval_int_1_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_int_1_m: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_int_mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_int_mnu2: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_int_qmnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_int_qmnu2: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_int_q2mnu2: not implemented.");

  return 0.0;
}

static gdouble
_ncm_csq1d_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  g_error ("_ncm_csq1d_eval_F1: not implemented.");

  return 0.0;
}

typedef struct _NcmCSQ1DWS
{
  NcmCSQ1D *csq1d;
  NcmModel *model;
  gdouble reltol;
  gdouble F1_min;
} NcmCSQ1DWS;

static gdouble _ncm_csq1d_F1_func (const gdouble t, gpointer user_data);

static gdouble
_ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0, 0.0};
  const gdouble nu             = ncm_csq1d_eval_nu (csq1d, model, t);
  const gdouble twonu          = 2.0 * nu;

  gdouble err;

  return ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F1_func, &ws, &err) / twonu;
}

/* State related functions */

/**
 * ncm_csq1d_state_new:
 *
 * Creates a new and uninitialized #NcmCSQ1DState.
 *
 * Returns: (transfer full): a new #NcmCSQ1DState.
 */
NcmCSQ1DState *
ncm_csq1d_state_new (void)
{
  return g_slice_new0 (NcmCSQ1DState);
}

/**
 * ncm_csq1d_state_copy:
 * @state: a #NcmCSQ1DState
 *
 * Creates a copy of @state.
 *
 * Returns: (transfer full): a new #NcmCSQ1DState with the same contents as @state.
 */
NcmCSQ1DState *
ncm_csq1d_state_copy (NcmCSQ1DState *state)
{
  return g_slice_dup (NcmCSQ1DState, state);
}

/**
 * ncm_csq1d_state_free:
 * @state: a #NcmCSQ1DState
 *
 * Frees @state.
 *
 */
void
ncm_csq1d_state_free (NcmCSQ1DState *state)
{
  g_slice_free (NcmCSQ1DState, state);
}

/**
 * ncm_csq1d_state_set_ag:
 * @state: a #NcmCSQ1DState
 * @frame: the frame of @state
 * @t: the time of @state
 * @alpha: the alpha of @state
 * @gamma: the gamma of @state
 *
 * Sets the state using the $(\alpha, \gamma)$ parametrization.
 *
 */
void
ncm_csq1d_state_set_ag (NcmCSQ1DState *state, const NcmCSQ1DFrame frame, const gdouble t, const gdouble alpha, const gdouble gamma)
{
  state->frame = frame;
  state->t     = t;
  state->alpha = alpha;
  state->gamma = gamma;
}

/**
 * ncm_csq1d_state_set_up:
 * @state: a #NcmCSQ1DState
 * @frame: the frame of @state
 * @t: the time of @state
 * @chi: the chi of @state
 * @Up: the Up of @state
 *
 * Sets the state using the $(\chi, U_+)$ parametrization.
 *
 */
void
ncm_csq1d_state_set_up (NcmCSQ1DState *state, const NcmCSQ1DFrame frame, const gdouble t, const gdouble chi, const gdouble Up)
{
  state->frame = frame;
  state->t     = t;
  state->alpha = asinh (chi);
  state->gamma = Up - gsl_sf_lncosh (state->alpha);
}

/**
 * ncm_csq1d_state_set_um:
 * @state: a #NcmCSQ1DState
 * @frame: the frame of @state
 * @t: the time of @state
 * @chi: the chi of @state
 * @Um: the Um of @state
 *
 * Sets the state using the $(\chi, U_-)$ parametrization.
 *
 */
void
ncm_csq1d_state_set_um (NcmCSQ1DState *state, const NcmCSQ1DFrame frame, const gdouble t, const gdouble chi, const gdouble Um)
{
  state->frame = frame;
  state->t     = t;
  state->alpha = asinh (chi);
  state->gamma = -Um + gsl_sf_lncosh (state->alpha);
}

/**
 * ncm_csq1d_state_get_time:
 * @state: a #NcmCSQ1DState
 *
 * Returns: the time of @state.
 */
gdouble
ncm_csq1d_state_get_time (NcmCSQ1DState *state)
{
  return state->t;
}

/**
 * ncm_csq1d_state_get_frame:
 * @state: a #NcmCSQ1DState
 *
 * Returns: the frame of @state.
 */
NcmCSQ1DFrame
ncm_csq1d_state_get_frame (NcmCSQ1DState *state)
{
  return state->frame;
}

/**
 * ncm_csq1d_state_get_ag:
 * @state: a #NcmCSQ1DState
 * @alpha: (out): the alpha of @state
 * @gamma: (out): the gamma of @state
 *
 * Computes the $(\alpha, \gamma)$ parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_ag (NcmCSQ1DState *state, gdouble *alpha, gdouble *gamma)
{
  *alpha = state->alpha;
  *gamma = state->gamma;
}

/**
 * ncm_csq1d_state_get_up:
 * @state: a #NcmCSQ1DState
 * @chi: (out): the chi of @state
 * @Up: (out): the Up of @state
 *
 * Computes the $(\chi, U_+)$ parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_up (NcmCSQ1DState *state, gdouble *chi, gdouble *Up)
{
  *chi = sinh (state->alpha);
  *Up  = state->gamma + gsl_sf_lncosh (state->alpha);
}

/**
 * ncm_csq1d_state_get_um:
 * @state: a #NcmCSQ1DState
 * @chi: (out): the chi of @state
 * @Um: (out): the Um of @state
 *
 * Computes the $(\chi, U_-)$ parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_um (NcmCSQ1DState *state, gdouble *chi, gdouble *Um)
{
  *chi = sinh (state->alpha);
  *Um  = -state->gamma + gsl_sf_lncosh (state->alpha);
}

/**
 * ncm_csq1d_state_get_J:
 * @state: a #NcmCSQ1DState
 * @J11: (out): the J11 of @state
 * @J12: (out): the J12 of @state
 * @J22: (out): the J22 of @state
 *
 * Computes the covariant metric of @state.
 *
 */
void
ncm_csq1d_state_get_J (NcmCSQ1DState *state, gdouble *J11, gdouble *J12, gdouble *J22)
{
  const gdouble alpha = state->alpha;
  const gdouble gamma = state->gamma;

  *J11 = cosh (alpha) * exp (-gamma);
  *J22 = cosh (alpha) * exp (+gamma);
  *J12 = -sinh (alpha);
}

/**
 * ncm_csq1d_state_get_phi_Pphi:
 * @state: a #NcmCSQ1DState
 * @phi: (out caller-allocates) (array fixed-size=2): the $\phi$ of @state
 * @Pphi: (out caller-allocates) (array fixed-size=2): the $P_\phi$ of @state
 *
 * Computes the $(\phi, P_\phi)$ parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_phi_Pphi (NcmCSQ1DState *state, gdouble *phi, gdouble *Pphi)
{
  const gdouble alpha               = state->alpha;
  const gdouble gamma               = state->gamma;
  const gdouble exp_gamma_p_alpha_2 = exp (0.5 * (gamma + alpha));
  const gdouble exp_gamma_m_alpha_2 = exp (0.5 * (gamma - alpha));

  phi[0] = +0.5 / exp_gamma_m_alpha_2;
  phi[1] = -0.5 / exp_gamma_p_alpha_2;

  Pphi[0] = -0.5 * exp_gamma_p_alpha_2;
  Pphi[1] = -0.5 * exp_gamma_m_alpha_2;
}

/**
 * ncm_csq1d_state_get_poincare_half_plane:
 * @state: a #NcmCSQ1DState
 * @x: (out): the $x$ of @state
 * @lny: (out): the $\ln y$ of @state
 *
 * Computes the Poincaré half-plane parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_poincare_half_plane (NcmCSQ1DState *state, gdouble *x, gdouble *lny)
{
  const gdouble alpha = state->alpha;
  const gdouble gamma = state->gamma;

  *x   = exp (-gamma) * tanh (alpha);
  *lny = -gamma - gsl_sf_lncosh (alpha);
}

/**
 * ncm_csq1d_state_get_poincare_disc:
 * @state: a #NcmCSQ1DState
 * @x: (out): the $x$ of @state
 * @y: (out): the $y$ of @state
 *
 * Computes the Poincaré disc parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_poincare_disc (NcmCSQ1DState *state, gdouble *x, gdouble *y)
{
  const gdouble alpha = state->alpha;
  const gdouble gamma = state->gamma;

  *x = +tanh (alpha) / (1.0 / cosh (alpha) + cosh (gamma));
  *y = -tanh (gamma) / (1.0 / (cosh (alpha) * cosh (gamma)) + 1.0);
}

/**
 * ncm_csq1d_state_get_minkowski:
 * @state: a #NcmCSQ1DState
 * @x1: (out): the $x_1$ of @state
 * @x2: (out): the $x_2$ of @state
 *
 * Computes the Minkowski parametrization of @state.
 *
 */
void
ncm_csq1d_state_get_minkowski (NcmCSQ1DState *state, gdouble *x1, gdouble *x2)
{
  const gdouble alpha = state->alpha;
  const gdouble gamma = state->gamma;

  *x1 = +sinh (alpha);
  *x2 = -cosh (alpha) * sinh (gamma);
}

static gdouble
_arcsinh_exp_x (const gdouble x)
{
  return x + log1p (sqrt (1.0 + exp (-2.0 * x)));
}

static gdouble
_arccosh_exp_x (const gdouble x)
{
  return x + log1p (sqrt (1.0 - exp (-2.0 * x)));
}

/**
 * ncm_csq1d_state_get_circle:
 * @state: a #NcmCSQ1DState
 * @r: radius
 * @theta: angle
 * @cstate: (out caller-allocates): the new state
 *
 * Computes the complex structure matrix parameters for a circle
 * around the point @state with radius $r$ and angle
 * $\theta$ and stores the result in @cstate.
 *
 */
void
ncm_csq1d_state_get_circle (NcmCSQ1DState *state, const gdouble r, const gdouble theta, NcmCSQ1DState *cstate)
{
  const gdouble alpha   = state->alpha;
  const gdouble gamma   = state->gamma;
  const gdouble ct      = cos (theta);
  const gdouble st      = sin (theta);
  const gdouble tr      = tanh (r);
  const gdouble ca      = cosh (alpha);
  const gdouble ta      = tanh (alpha);
  const gdouble lncr    = gsl_sf_lncosh (r);
  const gdouble lnca    = gsl_sf_lncosh (alpha);
  const gdouble f       = (ta + ct * tr);
  const gdouble absf    = fabs (f);
  const gdouble ln_absf = log (absf);
  const gdouble signf   = GSL_SIGN (f);
  const gdouble t1      = signf * _arcsinh_exp_x (lncr + lnca + ln_absf);
  const gdouble t2      = -(2.0 * st * tr / (ca * (1.0 + tr * ct * ta) + tr * st));

  ncm_csq1d_state_set_ag (cstate, state->frame, state->t, t1, gamma + 0.5 * log1p (t2));
}

/**
 * ncm_csq1d_state_compute_distance:
 * @state: a #NcmCSQ1DState
 * @state1: a #NcmCSQ1DState
 *
 * Computes the distance between @state and @state1.
 *
 * Returns: the distance between @state and @state1.
 */
gdouble
ncm_csq1d_state_compute_distance (NcmCSQ1DState *state, NcmCSQ1DState *state1)
{
  const gdouble dgamma01 = state->gamma - state1->gamma;

  g_assert (state->frame == state1->frame);
  g_assert (state->t == state1->t);

  if ((fabs (state->alpha) > 1.0) || (fabs (state1->alpha) > 1.0) || (fabs (dgamma01) > 1.0))
  {
    const gdouble lncoshalpha0 = gsl_sf_lncosh (state->alpha);
    const gdouble lncoshalpha1 = gsl_sf_lncosh (state1->alpha);
    const gdouble lncoshdgamma = gsl_sf_lncosh (dgamma01);
    const gdouble tanhalpha0   = tanh (state->alpha);
    const gdouble tanhalpha1   = tanh (state1->alpha);
    const gdouble f            = log1p (-tanhalpha0 * tanhalpha1 * exp (-lncoshdgamma));

    return _arccosh_exp_x (lncoshalpha0 + lncoshalpha1 + lncoshdgamma + f);
  }
  else
  {
    const gdouble a        = gsl_sf_lncosh (state->alpha + state1->alpha);
    const gdouble b        = gsl_sf_lncosh (state->alpha - state1->alpha);
    const gdouble c        = gsl_sf_lncosh (dgamma01);
    const gdouble expm1a   = expm1 (a);
    const gdouble expm1b   = expm1 (b);
    const gdouble expm1apc = expm1 (a + c);
    const gdouble expm1bpc = expm1 (b + c);
    const gdouble M12_m1   = 0.5 * (expm1apc + expm1bpc + expm1b - expm1a);
    const gdouble dist     =  asinh (sqrt (M12_m1 * (2.0 + M12_m1)));

    return dist;
  }
}

/* CSQ1D methods */

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

  if (self->reltol != reltol)
  {
    self->reltol = reltol;
    ncm_model_ctrl_force_update (self->ctrl);
  }
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

  if (self->abstol != abstol)
  {
    self->abstol = abstol;
    ncm_model_ctrl_force_update (self->ctrl);
  }
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
    ncm_model_ctrl_force_update (self->ctrl);
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

  if (self->tf != tf)
  {
    self->tf = tf;
    ncm_model_ctrl_force_update (self->ctrl);
  }
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

  if (self->adiab_threshold != adiab_threshold)
  {
    ncm_model_ctrl_force_update (self->ctrl);
    self->adiab_threshold = adiab_threshold;
  }
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
 * ncm_csq1d_set_init_cond:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @evol_state: a #NcmCSQ1DEvolState
 * @initial_state: a #NcmCSQ1DState
 *
 * Sets the values of the initial conditions to @initial_state.
 * Depending on the value of @evol_state, the initial conditions
 * are set in the adiabatic frame 1 if @evol_state is #NCM_CSQ1D_EVOL_STATE_ADIABATIC,
 * or in the original frame when using the $U_+$ or $U_-$ parametrization.
 *
 */
void
ncm_csq1d_set_init_cond (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DEvolState evol_state, NcmCSQ1DState *initial_state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (evol_state)
  {
    case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
      ncm_csq1d_change_frame (csq1d, model, initial_state, NCM_CSQ1D_FRAME_ADIAB1);
      ncm_csq1d_state_get_ag (initial_state, &NV_Ith_S (self->y, 0), &NV_Ith_S (self->y, 1));
      break;
    case NCM_CSQ1D_EVOL_STATE_UP:
      ncm_csq1d_change_frame (csq1d, model, initial_state, NCM_CSQ1D_FRAME_ORIG);
      ncm_csq1d_state_get_up (initial_state, &NV_Ith_S (self->y_Up, 0), &NV_Ith_S (self->y_Up, 1));
      break;
    case NCM_CSQ1D_EVOL_STATE_UM:
      ncm_csq1d_change_frame (csq1d, model, initial_state, NCM_CSQ1D_FRAME_ORIG);
      ncm_csq1d_state_get_um (initial_state, &NV_Ith_S (self->y_Um, 0), &NV_Ith_S (self->y_Um, 1));
      break;
    default:
      g_error ("ncm_csq1d_set_init_cond: state %d not supported", evol_state);
      break;
  }

  self->t     = ncm_csq1d_state_get_time (initial_state);
  self->state = evol_state;

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
  NcmCSQ1DState *state         = ncm_csq1d_compute_adiab (csq1d, model, ti, self->cur_state, NULL, NULL);

  g_assert (state->frame == NCM_CSQ1D_FRAME_ADIAB1);

  if ((fabs (state->gamma) > self->adiab_threshold) || (fabs (state->alpha) > self->adiab_threshold))
    g_error ("ncm_csq1d_set_init_cond_adiab: time ti == % 22.15g is not a valid adiabatic time alpha, dgamma == (% 22.15g, % 22.15g)",
             ti, state->alpha, state->gamma);
  else
    ncm_csq1d_set_init_cond (csq1d, model, NCM_CSQ1D_EVOL_STATE_ADIABATIC, state);
}

/**
 * ncm_csq1d_set_initial_condition_type:
 * @csq1d: a #NcmCSQ1D
 * @initial_condition_type: the vacuum type
 *
 * Sets the vacuum type to @initial_condition_type. The vacuum type is used
 * to determine the vacuum state when preparing the object using
 * ncm_csq1d_prepare_vacuum().
 *
 */
void
ncm_csq1d_set_initial_condition_type (NcmCSQ1D *csq1d, NcmCSQ1DVacuumType initial_condition_type)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  if (self->initial_condition_type != initial_condition_type)
  {
    ncm_model_ctrl_force_update (self->ctrl);
    self->initial_condition_type = initial_condition_type;
  }
}

/**
 * ncm_csq1d_set_vacuum_reltol:
 * @csq1d: a #NcmCSQ1D
 * @vacuum_reltol: relative tolerance
 *
 * Sets the relative tolerance for the vacuum definition. This tolerance
 * is used to determine the vacuum state when preparing the object using
 * ncm_csq1d_prepare_vacuum().
 *
 */
void
ncm_csq1d_set_vacuum_reltol (NcmCSQ1D *csq1d, const gdouble vacuum_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  if (self->vacuum_reltol != vacuum_reltol)
  {
    ncm_model_ctrl_force_update (self->ctrl);
    self->vacuum_reltol = vacuum_reltol;
  }
}

/**
 * ncm_csq1d_set_vacuum_max_time:
 * @csq1d: a #NcmCSQ1D
 * @vacuum_max_time: maximum time
 *
 * Sets the maximum time for the vacuum search. This time is used
 * to determine the vacuum state when preparing the object using
 * ncm_csq1d_prepare_vacuum().
 *
 */
void
ncm_csq1d_set_vacuum_max_time (NcmCSQ1D *csq1d, const gdouble vacuum_max_time)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  if (self->vacuum_max_time != vacuum_max_time)
  {
    ncm_model_ctrl_force_update (self->ctrl);
    self->vacuum_max_time = vacuum_max_time;
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
 * ncm_csq1d_get_initial_condition_type:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the vacuum type.
 */
NcmCSQ1DVacuumType
ncm_csq1d_get_initial_condition_type (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->initial_condition_type;
}

/**
 * ncm_csq1d_get_vacuum_reltol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the relative tolerance for the vacuum definition.
 */
gdouble
ncm_csq1d_get_vacuum_reltol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->vacuum_reltol;
}

/**
 * ncm_csq1d_get_vacuum_max_time:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the maximum time for the vacuum search.
 */
gdouble
ncm_csq1d_get_vacuum_max_time (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  return self->vacuum_max_time;
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

  const gdouble nu    = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t);
  const gdouble F1    = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, t);
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

  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
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

  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
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

  const gdouble nu    = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t);
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

  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
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

  const gdouble m         = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2       = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
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
 *
 * Returns: $\xi$
 */
/**
 * ncm_csq1d_eval_nu: (virtual eval_nu)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $\nu$
 */
/**
 * ncm_csq1d_eval_nu2: (virtual eval_nu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $\nu^2$
 */
/**
 * ncm_csq1d_eval_m: (virtual eval_m)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $m$
 */
/**
 * ncm_csq1d_eval_int_1_m: (virtual eval_int_1_m)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $\int 1/m \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_int_mnu2: (virtual eval_int_mnu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $\int m\nu^2 \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_int_qmnu2: (virtual eval_int_qmnu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $\int \left(\int 1/m \mathrm{d}t\right) m\nu^2 \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_int_q2mnu2: (virtual eval_int_q2mnu2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $\int \left(\int 1/m \mathrm{d}t\right)^2 m\nu^2 \mathrm{d}t$.
 */
/**
 * ncm_csq1d_eval_F1: (virtual eval_F1)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $F_1$
 */
/**
 * ncm_csq1d_eval_F2: (virtual eval_F2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 *
 * Returns: $F_2$
 */

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_adiabatic (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a, GArray *gamma_a)
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
      const gdouble xi         = ncm_csq1d_eval_xi (csq1d, model, self->t);
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
      const gdouble xi    = ncm_csq1d_eval_xi (csq1d, model, self->t);
      const gdouble gamma = dgamma + xi;

      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      g_array_append_val (gamma_a,   gamma);
      last_asinh_t = asinh_t;
    }

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }

    if ((fabs (dgamma) > self->adiab_threshold) && (fabs (alpha) > self->adiab_threshold))
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
_ncm_csq1d_evol_Up (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a, GArray *gamma_a)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DEvolStop reason      = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t         = asinh (self->t);
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
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t);
    m       = ncm_csq1d_eval_m (csq1d, model, self->t);
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
      g_array_append_val (gamma_a,   gamma);
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

    if (dgamma < 0.0)
    {
      reason = NCM_CSQ1D_EVOL_STOP_UM_START;
      break;
    }
  }

  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Um (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a, GArray *gamma_a)
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
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t);
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
      g_array_append_val (gamma_a,   gamma);
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
_ncm_csq1d_evol_save (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DWS *ws, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a, GArray *gamma_a)
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
      stop = _ncm_csq1d_evol_adiabatic (csq1d, model, asinh_t_a, alpha_a, dgamma_a, gamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STATE_UP:
    {
      _ncm_csq1d_prepare_integrator_Up (csq1d, ws);
      stop = _ncm_csq1d_evol_Up (csq1d, ws, model, asinh_t_a, alpha_a, dgamma_a, gamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STATE_UM:
    {
      _ncm_csq1d_prepare_integrator_Um (csq1d, ws);
      stop = _ncm_csq1d_evol_Um (csq1d, model, asinh_t_a, alpha_a, dgamma_a, gamma_a);
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
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t);
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
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a, gamma_a);
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
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t);
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
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a, gamma_a);
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
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t);
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
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t);
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
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a, gamma_a);
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

/**
 * ncm_csq1d_prepare: (virtual prepare)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 *
 * Prepares the object using @model. This method expects that the initial
 * conditions have been set. It computes the evolution of the system and
 * saves the results in the internal splines.
 *
 */
void
ncm_csq1d_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0, 0.0};

  if (NCM_CSQ1D_GET_CLASS (csq1d)->prepare != NULL)
    NCM_CSQ1D_GET_CLASS (csq1d)->prepare (csq1d, model);

  g_assert (self->init_cond_set);
  g_assert_cmpfloat (self->tf, >, self->ti);

  if (self->save_evol)
  {
    GArray *asinh_t_a = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *alpha_a   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *dgamma_a  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *gamma_a   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);

    _ncm_csq1d_evol_save (csq1d, model, &ws, asinh_t_a, alpha_a, dgamma_a, gamma_a);

    ncm_spline_set_array (self->alpha_s,  asinh_t_a, alpha_a,  TRUE);
    ncm_spline_set_array (self->dgamma_s, asinh_t_a, dgamma_a, TRUE);
    ncm_spline_set_array (self->gamma_s,  asinh_t_a, gamma_a,  TRUE);

    g_array_unref (asinh_t_a);
    g_array_unref (alpha_a);
    g_array_unref (dgamma_a);
    g_array_unref (gamma_a);
  }
  else
  {
    g_assert_not_reached ();
  }
}

/**
 * ncm_csq1d_prepare_vacuum:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 *
 * Prepares the object using @model for the vacuum case. This method will set the
 * initial conditions to the vacuum case and compute the evolution of the system. The
 * results are saved in the internal splines. See ncm_csq1d_set_initial_condition_type(),
 * ncm_csq1d_set_vacuum_reltol() and ncm_csq1d_set_vacuum_max_time() for more
 * information on how to set the vacuum parameters.
 *
 * If the vacuum cannot be computed, this method will return %FALSE.
 *
 * Returns: %TRUE if the vacuum was found and the solution was computed.
 */
gboolean
ncm_csq1d_prepare_vacuum (NcmCSQ1D *csq1d, NcmModel *model)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  gdouble vacuum_final_time;

  if (!ncm_csq1d_find_adiab_time_limit (csq1d, model, self->ti, self->vacuum_max_time, self->vacuum_reltol, &vacuum_final_time))
  {
    return FALSE;
  }
  else
  {
    NcmCSQ1DState state;
    gdouble zeta_k_tau, alpha_reltol, dgamma_reltol;

    /* If the final time is greater than the vacuum time, we need to compute the numerical solution
     * from the vacuum time to the final time. Otherwise, we just need to compute the adiabatic
     * solution from the initial time to the final time.
     */
    if (self->tf > vacuum_final_time)
    {
      ncm_csq1d_compute_adiab (csq1d, model, vacuum_final_time, &state, &alpha_reltol, &dgamma_reltol);
      ncm_csq1d_set_init_cond_adiab (csq1d, model, vacuum_final_time);
      ncm_csq1d_prepare (csq1d, model);
    }
  }

  return TRUE;
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

static void _ncm_csq1d_eval_adiab_at_no_test (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol);

gdouble
_ncm_csq1d_find_adiab_time_limit_f (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) params;
  /*NcmCSQ1DPrivate * const self = ws->ncm_csq1d_get_instance_private (csq1d);*/

  gdouble alpha, dgamma, cmp, alpha_reltol, dgamma_reltol;

  _ncm_csq1d_eval_adiab_at_no_test (ws->csq1d, ws->model, t, &alpha, &dgamma, &alpha_reltol, &dgamma_reltol);

  cmp = MAX (alpha_reltol, dgamma_reltol);

  if (cmp == 0.0)
    return -GSL_DBL_MAX;
  else
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

  _ncm_csq1d_eval_adiab_at_no_test (csq1d, model, t0, &alpha0, &dgamma0, &alpha_reltol0, &dgamma_reltol0);
  _ncm_csq1d_eval_adiab_at_no_test (csq1d, model, t1, &alpha1, &dgamma1, &alpha_reltol1, &dgamma_reltol1);

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
      g_warning ("# Impossible to find the adiabatic limit: \n"
                 "\tt0 % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n"
                 "\tt1 % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n",
                 t0, alpha0, alpha_reltol0, dgamma0, dgamma_reltol0,
                 t1, alpha1, alpha_reltol1, dgamma1, dgamma_reltol1);

    return FALSE;
  }
  else
  {
    NcmCSQ1DWS ws = {csq1d, model, reltol, 0.0};
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
  const gdouble F1             = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, t);

  return F1;
}

static gdouble
_ncm_csq1d_F2_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F2             = ncm_csq1d_eval_F2 (ws->csq1d, ws->model, t);

  return F2;
}

static gdouble
_ncm_csq1d_lnnu_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble nu             = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t);

  return log (nu);
}

static gdouble _ncm_csq1d_abs_F1_asinht (gdouble at, gpointer user_data);
static gdouble _ncm_csq1d_ln_abs_F1_eps_asinht (gdouble at, gpointer user_data);

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
  NcmCSQ1DWS ws                = {csq1d, model, border_eps, 0.0};

  gsl_min_fminimizer *fmin = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  const gdouble atl        = asinh (t0);
  const gdouble atu        = asinh (t1);
  gdouble at0              = atl;
  gdouble at1              = atu;
  gdouble atm              = (at0 + at1) * 0.5;
  const guint linsearch    = 100;
  guint iter               = 0;
  gint status;
  gsl_function F;
  gint ret;
  guint i;

  F.params   = &ws;
  F.function = &_ncm_csq1d_abs_F1_asinht;

  {
    gdouble test = GSL_POSINF;

    for (i = 0; i < linsearch; i++)
    {
      const gdouble at    = at0 + (at1 - at0) * i / (linsearch - 1);
      const gdouble test0 = GSL_FN_EVAL (&F, at);

      if (test0 < test)
      {
        test = test0;
        atm  = at;
      }
    }
  }

  /* Testing for minimum on the edge */
  if (ncm_cmp (atm, atl, 1.0e-15, 0.0) == 0)
  {
    atm = atl;
    at0 = atl;
    at1 = atl;
  }
  else if (ncm_cmp (atm, atu, 1.0e-15, 0.0) == 0)
  {
    atm = atu;
    at0 = atu;
    at1 = atu;
  }
  else
  {
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
  }

  gsl_min_fminimizer_free (fmin);
  ws.F1_min = ncm_csq1d_eval_F1 (csq1d, model, sinh (atm));

  {
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    guint max_iter = 1000;

    iter = 0;

    F.function = &_ncm_csq1d_ln_abs_F1_eps_asinht;
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

      /* ncm_message ("Bl: [%d] % 22.15e % 22.15e % 22.15e\n", status, sinh (t_Bl[0]), sinh (at0), sinh (at1)); */
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_set (s, &F, atm, atu);

    do {
      iter++;
      status  = gsl_root_fsolver_iterate (s);
      t_Bu[0] = gsl_root_fsolver_root (s);
      at0     = gsl_root_fsolver_x_lower (s);
      at1     = gsl_root_fsolver_x_upper (s);
      status  = gsl_root_test_interval (at0, at1, 0.0, 1.0e-7);

      /* ncm_message ("Bu: [%d] % 22.15e % 22.15e % 22.15e\n", status, sinh (t_Bu[0]), sinh (at0), sinh (at1)); */
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
  }

  {
    const gdouble tm = sinh (atm);

    F1_min[0] = ncm_csq1d_eval_F1 (csq1d, model, tm);
    t_Bl[0]   = sinh (t_Bl[0]);
    t_Bu[0]   = sinh (t_Bu[0]);

    return tm;
  }
}

static gdouble
_ncm_csq1d_abs_F1_asinht (gdouble at, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, sinh (at));

  return fabs (F1);
}

static gdouble
_ncm_csq1d_ln_abs_F1_eps_asinht (gdouble at, gpointer user_data)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, sinh (at));

  return fabs ((F1 - ws->F1_min) / ws->reltol) - 1.0;
}

static void
_ncm_csq1d_compute_adiab2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t);
  const gdouble F2             = ncm_csq1d_eval_F2 (csq1d, model, t);

  state->frame = NCM_CSQ1D_FRAME_ADIAB1;
  state->t     = t;
  state->alpha = +F1;
  state->gamma = -F2;

  if (alpha_reltol != NULL)
    alpha_reltol[0] = fabs (F1);

  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = fabs (F2);
}

static void
_ncm_csq1d_compute_adiab4 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0, 0.0};
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t);
  const gdouble F2             = ncm_csq1d_eval_F2 (csq1d, model, t);
  const gdouble F1_2           = F1 * F1;
  const gdouble F1_3           = F1_2 * F1;
  const gdouble nu             = ncm_csq1d_eval_nu (csq1d, model, t);
  const gdouble twonu          = 2.0 * nu;
  gdouble err, F3, d2F2, dlnnu, F4, alpha_reltol0, dgamma_reltol0;

  F3    = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err) / twonu;
  d2F2  = ncm_diff_rc_d2_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err);
  dlnnu = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_lnnu_func, &ws, &err);
  F4    = d2F2 / gsl_pow_2 (twonu) - dlnnu * F3 / twonu;

  state->frame = NCM_CSQ1D_FRAME_ADIAB1;
  state->t     = t;

  if ((fabs (F3) > fabs (F2)) || (fabs (F4) > fabs (F3)))
  {
    g_warning ("WKB series with |F3| > |F2| or |F4| > |F3|, "
               "|F3/F2| = % 22.15g and |F4/F3| = % 22.15g "
               "at t = % 22.15e. Truncating series at F2.\n",
               fabs (F3 / F2), fabs (F4 / F3), t);

    alpha_reltol0  = fabs (F1);
    dgamma_reltol0 = fabs (F2);
    state->alpha   = +F1;
    state->gamma   = -F2;
  }
  else if ((F1 == 0.0) || (F2 == 0.0))
  {
    g_warning ("WKB series with F1 = 0 or F2 = 0 at t = % 22.15e. "
               "Truncating series at F2, cannot estimate error.\n", t);

    alpha_reltol0  = 0.0;
    dgamma_reltol0 = 0.0;
    state->alpha   = +F1;
    state->gamma   = -F2;
  }
  else
  {
    alpha_reltol0  = gsl_pow_2 ((F1_3 / 3.0 - F3) / F1);
    dgamma_reltol0 = gsl_pow_2 ((F4 - F1_2 * F2) / F2);
    state->alpha   = +F1 + F1_3 / 3.0 - F3;
    state->gamma   = -(1.0 + F1_2) * F2 + F4;
  }

  if (alpha_reltol != NULL)
    alpha_reltol[0] = alpha_reltol0;

  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = dgamma_reltol0;
}

static void
_ncm_csq1d_compute_adiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (self->initial_condition_type)
  {
    case NCM_CSQ1D_VACUUM_TYPE_ADIABATIC2:
      _ncm_csq1d_compute_adiab2 (csq1d, model, t, state, alpha_reltol, dgamma_reltol);
      break;
    case NCM_CSQ1D_VACUUM_TYPE_ADIABATIC4:
      _ncm_csq1d_compute_adiab4 (csq1d, model, t, state, alpha_reltol, dgamma_reltol);
      break;
    default:
      g_error ("_ncm_csq1d_compute_adiab: vacuum type '%d' not supported.", self->initial_condition_type);
      break;
  }
}

/**
 * ncm_csq1d_compute_adiab:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 * @alpha_reltol: (out) (allow-none): estimated error on $\alpha(t)$
 * @dgamma_reltol: (out) (allow-none): estimated error on $\Delta\gamma(t)$
 *
 * Computes the value of the adiabatic approximation of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 * This method computes the adiabatic approximation using the adiabatic series up to the order 2 or 4,
 * depending on the value of the property max-order-2. The result is stored in the state object in the
 * frame NCM_CSQ1D_FRAME_ADIAB1. Use ncm_csq1d_change_frame() to change the frame.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_compute_adiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  _ncm_csq1d_compute_adiab (csq1d, model, t, state, alpha_reltol, dgamma_reltol);

  return state;
}

/**
 * ncm_csq1d_compute_adiab_frame:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @frame: the frame to change
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 * @alpha_reltol: (out) (allow-none): estimated error on $\alpha(t)$
 * @dgamma_reltol: (out) (allow-none): estimated error on $\Delta\gamma(t)$
 *
 * As ncm_csq1d_compute_adiab(), but changes the frame of the result to @frame.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_compute_adiab_frame (NcmCSQ1D *csq1d, NcmModel *model, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  _ncm_csq1d_compute_adiab (csq1d, model, t, state, alpha_reltol, dgamma_reltol);
  ncm_csq1d_change_frame (csq1d, model, state, frame);

  return state;
}

static void
_ncm_csq1d_eval_adiab2_at_no_test (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t);
  const gdouble F2             = ncm_csq1d_eval_F2 (csq1d, model, t);

  alpha[0]  = +F1;
  dgamma[0] = -F2;

  if (alpha_reltol != NULL)
    alpha_reltol[0] = fabs (F1);

  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = fabs (F2);
}

static void
_ncm_csq1d_eval_adiab4_at_no_test (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, 0.0, 0.0};
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t);
  const gdouble F2             = ncm_csq1d_eval_F2 (csq1d, model, t);
  const gdouble F1_2           = F1 * F1;
  const gdouble F1_3           = F1_2 * F1;
  const gdouble nu             = ncm_csq1d_eval_nu (csq1d, model, t);
  const gdouble twonu          = 2.0 * nu;
  gdouble err, F3, d2F2, dlnnu, F4, alpha_reltol0, dgamma_reltol0;

  F3    = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err) / twonu;
  d2F2  = ncm_diff_rc_d2_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err);
  dlnnu = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_lnnu_func, &ws, &err);
  F4    = d2F2 / gsl_pow_2 (twonu) - dlnnu * F3 / twonu;

  if (F1 != 0.0)
    alpha_reltol0 = gsl_pow_2 ((F1_3 / 3.0 - F3) / F1);
  else
    alpha_reltol0 = gsl_pow_2 (F1_3 / 3.0 - F3);

  if (F2 != 0.0)
    dgamma_reltol0 = gsl_pow_2 ((F4 - F1_2 * F2) / F2);
  else
    dgamma_reltol0 = gsl_pow_2 (F4 - F1_2 * F2);

  alpha[0]  = +F1 + F1_3 / 3.0 - F3;
  dgamma[0] = -(1.0 + F1_2) * F2 + F4;

  if (alpha_reltol != NULL)
    alpha_reltol[0] = alpha_reltol0;

  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = dgamma_reltol0;
}

static void
_ncm_csq1d_eval_adiab_at_no_test (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (self->initial_condition_type)
  {
    case NCM_CSQ1D_VACUUM_TYPE_ADIABATIC2:
      _ncm_csq1d_eval_adiab2_at_no_test (csq1d, model, t, alpha, dgamma, alpha_reltol, dgamma_reltol);
      break;
    case NCM_CSQ1D_VACUUM_TYPE_ADIABATIC4:
      _ncm_csq1d_eval_adiab4_at_no_test (csq1d, model, t, alpha, dgamma, alpha_reltol, dgamma_reltol);
      break;
    default:
      g_error ("_ncm_csq1d_eval_adiab_at_no_test: vacuum type '%d' not supported.", self->initial_condition_type);
      break;
  }
}

static void
_ncm_csq1d_compute_nonadiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble p1             = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t);
  const gdouble r1             = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t);

  ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_NONADIAB1, t, +2.0 * r1, -2.0 * p1);
}

/*  OLD IMPLEMENTATION
 *  const gdouble q0             = ncm_csq1d_eval_int_1_m (csq1d, model, t);
 *  const gdouble q1             = ncm_csq1d_eval_int_q2mnu2 (csq1d, model, t);
 *  const gdouble p1             = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t);
 *  const gdouble r1             = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t);
 *
 *  switch (frame)
 *  {
 *   case NCM_CSQ1D_FRAME_ORIG:
 *     chi[0] = (2.0 * p1 - 1.0) * (q0 + q1) + 2.0 * r1;
 *     Up[0]  = -2.0 * p1;
 *     break;
 *   case NCM_CSQ1D_FRAME_NONADIAB1:
 *     chi[0] = +2.0 * r1;
 *     Up[0]  = -2.0 * p1;
 *     break;
 *   case NCM_CSQ1D_FRAME_NONADIAB2:
 *     chi[0] = 0.0;
 *     Up[0]  = 0.0;
 *     break;
 *   default:
 *     g_assert_not_reached ();
 *     break;
 *  }
 */

/**
 * ncm_csq1d_compute_nonadiab:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 *
 * Computes the value of the non-adiabatic VDC order two in the variables $\chi$ and
 * $U$ at $t$. The result is stored in the state object in the frame NCM_CSQ1D_FRAME_NONADIAB1.
 * Use ncm_csq1d_change_frame() to change the frame.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_compute_nonadiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state)
{
  _ncm_csq1d_compute_nonadiab (csq1d, model, t, state);

  return state;
}

/**
 * ncm_csq1d_compute_H:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 *
 * Computes the Hamiltonian vector state of the original frame at $t$. The result is
 * stored in the state object in the frame NCM_CSQ1D_FRAME_ORIG.
 * Use ncm_csq1d_change_frame() to change the frame.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_compute_H (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state)
{
  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ORIG, t, 0.0, ncm_csq1d_eval_xi (csq1d, model, t));

  return state;
}

static void
_ncm_csq1d_eval_state (NcmCSQ1D *csq1d, const gdouble t, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble a_t            = asinh (t);
  const gdouble alpha          = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble gamma          = ncm_spline_eval (self->gamma_s, a_t);

  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ORIG, t, alpha, gamma);
}

/**
 * ncm_csq1d_eval_at:
 * @csq1d: a #NcmCSQ1D
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 *
 * Computes the system state at $t$, the result is stored in the state object
 * in the frame NCM_CSQ1D_FRAME_ORIG. Use ncm_csq1d_change_frame() to change the frame.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_eval_at (NcmCSQ1D *csq1d, const gdouble t, NcmCSQ1DState *state)
{
  _ncm_csq1d_eval_state (csq1d, t, state);

  return state;
}

/**
 * ncm_csq1d_eval_at_frame:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @frame: a #NcmCSQ1DFrame
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 *
 * Computes the system state at $t$, the result is stored in the state object
 * in the frame @frame.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_eval_at_frame (NcmCSQ1D *csq1d, NcmModel *model, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state)
{
  _ncm_csq1d_eval_state (csq1d, t, state);

  ncm_csq1d_change_frame (csq1d, model, state, frame);

  return state;
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
  const gdouble l_theta = _ALPHA_TO_THETA (alpha);
  const gdouble l_gamma = gamma;
  gdouble thetap        = 0.0;

  ncm_util_mln_1mIexpzA_1pIexpmzA (l_gamma, l_theta, tanh (0.5 * p), gammap, &thetap);
  alphap[0] = _THETA_TO_ALPHA (thetap);
}

static void
_ncm_csq1d_change_adiab1_to_orig (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;
  const gdouble xi             = ncm_csq1d_eval_xi (csq1d, model, t);

  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ORIG, t, state->alpha, state->gamma + xi);
}

static void
_ncm_csq1d_change_orig_to_adiab1 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;
  const gdouble xi             = ncm_csq1d_eval_xi (csq1d, model, t);

  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ADIAB1, t, state->alpha, state->gamma - xi);
}

static void
_ncm_csq1d_change_adiab2_to_adiab1 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t);
  const gdouble xi1            = atanh (-F1);
  gdouble alphap               = 0.0;
  gdouble gammap               = 0.0;

  if (fabs (F1) >= 1.0)
    g_error ("ncm_csq1d_change_adiab2_to_adiab1: |F1| >= 1 at t = % 22.15e", t);

  _ncm_csq1d_ct_g0g1 (state->alpha, state->gamma, -xi1, &alphap, &gammap);
  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ADIAB1, t, alphap, gammap);
}

static void
_ncm_csq1d_change_adiab1_to_adiab2 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;
  const gdouble F1             = ncm_csq1d_eval_F1 (csq1d, model, t);
  const gdouble xi1            = atanh (-F1);
  gdouble alphap               = 0.0;
  gdouble gammap               = 0.0;

  if (fabs (F1) >= 1.0)
    g_error ("ncm_csq1d_change_adiab1_to_adiab2: |F1| >= 1 at t = % 22.15e", t);

  _ncm_csq1d_ct_g0g1 (state->alpha, state->gamma, xi1, &alphap, &gammap);
  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ADIAB2, t, alphap, gammap);
}

static void
_ncm_csq1d_change_nonadiab1_to_orig (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;
  const gdouble int_1_m        = ncm_csq1d_eval_int_1_m (csq1d, model, t);
  const gdouble q0             = int_1_m;
  const gdouble q1             = ncm_csq1d_eval_int_q2mnu2 (csq1d, model, t);
  gdouble chi, Up;

  ncm_csq1d_state_get_up (state, &chi, &Up);
  ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_ORIG, t, chi - exp (Up) * (q0 + q1), Up);
}

static void
_ncm_csq1d_change_orig_to_nonadiab1 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;
  const gdouble int_1_m        = ncm_csq1d_eval_int_1_m (csq1d, model, t);
  const gdouble q0             = int_1_m;
  const gdouble q1             = ncm_csq1d_eval_int_q2mnu2 (csq1d, model, t);
  gdouble chi, Up;

  ncm_csq1d_state_get_up (state, &chi, &Up);
  ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_NONADIAB1, t, chi + exp (Up) * (q0 + q1), Up);
}

static void
_ncm_csq1d_change_nonadiab1_to_nonadiab2 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;

  /* p1 g1 */
  {
    const gdouble p1 = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t);
    gdouble chi, Up;

    ncm_csq1d_state_get_up (state, &chi, &Up);

    Up = Up + 2.0 * p1;

    ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_NONADIAB2, t, chi, Up);
  }

  /* r1 g2 */
  {
    const gdouble r1 = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t);
    gdouble alpha, gamma, alphap, gammap;

    ncm_csq1d_state_get_ag (state, &alpha, &gamma);

    _ncm_csq1d_ct_g0g1 (alpha, gamma, -2.0 * r1, &alphap, &gammap);

    ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_NONADIAB2, t, alphap, gammap);
  }
}

static void
_ncm_csq1d_change_nonadiab2_to_nonadiab1 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble t              = state->t;

  /* r1 g2 */
  {
    const gdouble r1 = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t);
    gdouble alpha, gamma, alphap, gammap;

    ncm_csq1d_state_get_ag (state, &alpha, &gamma);

    _ncm_csq1d_ct_g0g1 (alpha, gamma, +2.0 * r1, &alphap, &gammap);

    ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_NONADIAB1, t, alphap, gammap);
  }

  {
    const gdouble p1 = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t);
    gdouble chi, Up;

    ncm_csq1d_state_get_up (state, &chi, &Up);

    /* p1 g1 */
    Up = Up - 2.0 * p1;

    ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_NONADIAB1, t, chi, Up);
  }
}

static void
_ncm_csq1d_change_frame_to_orig (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (state->frame)
  {
    case NCM_CSQ1D_FRAME_ORIG:
      break;
    case NCM_CSQ1D_FRAME_ADIAB1:
      _ncm_csq1d_change_adiab1_to_orig (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB2:
      _ncm_csq1d_change_adiab2_to_adiab1 (csq1d, model, state);
      _ncm_csq1d_change_adiab1_to_orig (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB1:
      _ncm_csq1d_change_nonadiab1_to_orig (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB2:
      _ncm_csq1d_change_nonadiab2_to_nonadiab1 (csq1d, model, state);
      _ncm_csq1d_change_nonadiab1_to_orig (csq1d, model, state);
      break;
    default:
      g_error ("_ncm_csq1d_change_frame_to_orig: Invalid frame %d", state->frame);
      break;
  }
}

static void
_ncm_csq1d_change_frame_to_adiab1 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (state->frame)
  {
    case NCM_CSQ1D_FRAME_ORIG:
      _ncm_csq1d_change_orig_to_adiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB1:
      break;
    case NCM_CSQ1D_FRAME_ADIAB2:
      _ncm_csq1d_change_adiab2_to_adiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB1:
      _ncm_csq1d_change_nonadiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_adiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB2:
      _ncm_csq1d_change_nonadiab2_to_nonadiab1 (csq1d, model, state);
      _ncm_csq1d_change_nonadiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_adiab1 (csq1d, model, state);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static void
_ncm_csq1d_change_frame_to_adiab2 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (state->frame)
  {
    case NCM_CSQ1D_FRAME_ORIG:
      _ncm_csq1d_change_orig_to_adiab1 (csq1d, model, state);
      _ncm_csq1d_change_adiab1_to_adiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB1:
      _ncm_csq1d_change_adiab1_to_adiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB2:
      break;
    case NCM_CSQ1D_FRAME_NONADIAB1:
      _ncm_csq1d_change_nonadiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_adiab1 (csq1d, model, state);
      _ncm_csq1d_change_adiab1_to_adiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB2:
      _ncm_csq1d_change_nonadiab2_to_nonadiab1 (csq1d, model, state);
      _ncm_csq1d_change_nonadiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_adiab1 (csq1d, model, state);
      _ncm_csq1d_change_adiab1_to_adiab2 (csq1d, model, state);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static void
_ncm_csq1d_change_frame_to_nonadiab1 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (state->frame)
  {
    case NCM_CSQ1D_FRAME_ORIG:
      _ncm_csq1d_change_orig_to_nonadiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB1:
      _ncm_csq1d_change_adiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_nonadiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB2:
      _ncm_csq1d_change_adiab2_to_adiab1 (csq1d, model, state);
      _ncm_csq1d_change_adiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_nonadiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB1:
      break;
    case NCM_CSQ1D_FRAME_NONADIAB2:
      _ncm_csq1d_change_nonadiab2_to_nonadiab1 (csq1d, model, state);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static void
_ncm_csq1d_change_frame_to_nonadiab2 (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (state->frame)
  {
    case NCM_CSQ1D_FRAME_ORIG:
      _ncm_csq1d_change_orig_to_nonadiab1 (csq1d, model, state);
      _ncm_csq1d_change_nonadiab1_to_nonadiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB1:
      _ncm_csq1d_change_adiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_nonadiab1 (csq1d, model, state);
      _ncm_csq1d_change_nonadiab1_to_nonadiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB2:
      _ncm_csq1d_change_adiab2_to_adiab1 (csq1d, model, state);
      _ncm_csq1d_change_adiab1_to_orig (csq1d, model, state);
      _ncm_csq1d_change_orig_to_nonadiab1 (csq1d, model, state);
      _ncm_csq1d_change_nonadiab1_to_nonadiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB1:
      _ncm_csq1d_change_nonadiab1_to_nonadiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB2:
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_csq1d_change_frame:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @state: a #NcmCSQ1DState
 * @frame: which frame to use
 *
 * Changes the frame of the @state object to the given @frame. The state object
 * must be a valid state object, it cannot be NULL. The state object is updated
 * in place.
 *
 * Returns: (transfer none): the state object in the new frame.
 */
NcmCSQ1DState *
ncm_csq1d_change_frame (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state, NcmCSQ1DFrame frame)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);

  switch (frame)
  {
    case NCM_CSQ1D_FRAME_ORIG:
      _ncm_csq1d_change_frame_to_orig (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB1:
      _ncm_csq1d_change_frame_to_adiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_ADIAB2:
      _ncm_csq1d_change_frame_to_adiab2 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB1:
      _ncm_csq1d_change_frame_to_nonadiab1 (csq1d, model, state);
      break;
    case NCM_CSQ1D_FRAME_NONADIAB2:
      _ncm_csq1d_change_frame_to_nonadiab2 (csq1d, model, state);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  return state;
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

  const gdouble m       = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2     = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
  const gdouble mnu2    = m * nu2;
  const gdouble mnu2_2  = 0.5 * mnu2;
  const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t);
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

  const gdouble m       = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2     = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
  const gdouble mnu2    = m * nu2;
  const gdouble mnu2_2  = 0.5 * mnu2;
  const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t);
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
  NcmCSQ1DWS ws                = {csq1d, model, 0.0, 0.0};

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
 * Computes the propagator for the given @csq1d and @model from @ti to @tf. The
 * propagator is computed using the integral approximation time @tii. The
 * propagator is stored in the @csq1d object and can be used to compute the
 * propagator a state from @ti to any time between @tii and @tf.
 *
 */
void
ncm_csq1d_prepare_prop (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti, const gdouble tii, const gdouble tf)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  NcmCSQ1DWS ws                = {csq1d, model, self->prop_threshold, 0.0};
  GArray *t_a                  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  gdouble t                    = 0.0;
  gboolean tf_Prop_found       = FALSE;
  GArray *R_a[4];
  gdouble u1[3];
  gint i;

  _ncm_csq1d_prepare_prop_eval_u1 (csq1d, model, ti, tii, u1);
  self->ti_Prop = ti;

/*
 *  ncm_message ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", ti, tii,
 *     u1[0], +ncm_csq1d_eval_int_qmnu2 (csq1d, model, tii),
 *     u1[1], +ncm_csq1d_eval_int_q2mnu2 (csq1d, model, tii),
 *     u1[2], -ncm_csq1d_eval_int_mnu2 (csq1d, model, tii)
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

  const gdouble m   = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2 = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);

  return m * nu2;
}

static gdouble
_ncm_csq1d_prepare_prop_qmnu2 (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble m    = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2  = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
  const gdouble mnu2 = m * nu2;

  if (mnu2 == 0)
  {
    return 0.0;
  }
  else
  {
    const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t);
    const gdouble q       = -int_1_m;

    return q * mnu2;
  }
}

static gdouble
_ncm_csq1d_prepare_prop_q2mnu2 (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws               = (NcmCSQ1DWS *) params;
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (ws->csq1d);

  const gdouble m    = ncm_csq1d_eval_m (ws->csq1d, ws->model, t);
  const gdouble nu2  = ncm_csq1d_eval_nu2 (ws->csq1d, ws->model, t);
  const gdouble mnu2 = m * nu2;

  if (mnu2 == 0)
  {
    return 0.0;
  }
  else
  {
    const gdouble int_1_m = ncm_csq1d_eval_int_1_m (ws->csq1d, ws->model, t);
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
 * ncm_csq1d_get_prop_vector:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @state: a #NcmCSQ1DState to store the result
 *
 * Computes the state vector associated with the propagator at time $t$.
 *
 * Returns: (transfer none): the @state object with the result.
 */
NcmCSQ1DState *
ncm_csq1d_compute_prop_vector (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble b              = ncm_spline_eval (self->R[1], t);
  const gdouble c              = ncm_spline_eval (self->R[2], t);
  const gdouble h              = ncm_spline_eval (self->R[3], t);
  const gdouble n0             = sqrt (-(b * c + h * h));

  ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_ORIG, t, h / n0, log (-c / n0));

  return state;
}

static void
_ncm_csq1d_evolve_prop_vector (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *initial_state, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state)
{
  NcmCSQ1DPrivate * const self = ncm_csq1d_get_instance_private (csq1d);
  const gdouble int_1_m        = ncm_csq1d_eval_int_1_m (csq1d, model, t);
  const gdouble q0             = int_1_m;
  const gdouble q1             = ncm_csq1d_eval_int_q2mnu2 (csq1d, model, t);
  const gdouble a              = ncm_spline_eval (self->R[0], t);
  const gdouble b              = ncm_spline_eval (self->R[1], t);
  const gdouble c              = ncm_spline_eval (self->R[2], t);
  const gdouble h              = ncm_spline_eval (self->R[3], t);
  const gdouble a11            = a + h;
  const gdouble a22            = a - h;
  const gdouble a12            = b;
  const gdouble a21            = c;
  const gdouble init_ti        = ncm_csq1d_state_get_time (initial_state);
  const gdouble q1_ti          = ncm_csq1d_eval_int_qmnu2 (csq1d, model, init_ti);
  gdouble chi_i, Up_i;

  if (init_ti != self->ti_Prop)
    g_error ("Initial time does not match the propagator initial time");

  if (initial_state->frame != NCM_CSQ1D_FRAME_NONADIAB1)
    g_error ("Initial state must be in the non-adiabatic 1 frame");

  ncm_csq1d_state_get_up (initial_state, &chi_i, &Up_i);

  /*
   * The propagator is computed in a intermediate frame, we need to change from the
   * non-adibatic 1 to the intermediate frame first. This transformation does not
   * change Up_i, so we can use it to compute the new chi_i.
   */
  chi_i =  (chi_i - exp (Up_i) * q1_ti);

  {
    const gdouble onepchi2 = 1.0 + chi_i * chi_i;
    const gdouble exp_Up_i = exp (Up_i);
    gdouble chi, Up;

    chi = a12 * a21 * chi_i + a22 * (a11 * chi_i - a12 * exp_Up_i) - a11 * a21 * onepchi2 / exp_Up_i;
    Up  = -(2.0 * a21 * a22 * chi_i - a22 * a22 * exp_Up_i - a21 * a21 * onepchi2 / exp_Up_i);

    /*
     * The transformation above result in a new chi and Up in the intermediate frame.
     * We can now use these values and transform to the final frame.
     */

    switch (frame)
    {
      case NCM_CSQ1D_FRAME_ORIG:
        chi = (chi - Up * q0);
        Up  = log (Up);

        ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_ORIG, t, chi, Up);
        break;
      case NCM_CSQ1D_FRAME_NONADIAB1:
        /* q1 L2p */
        chi = (chi + Up * q1);
        Up  = log (Up);

        ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_NONADIAB1, t, chi, Up);
        break;
      case NCM_CSQ1D_FRAME_NONADIAB2:
      {
        const gdouble p1 = ncm_csq1d_eval_int_qmnu2 (csq1d, model, t);

        /* q1 L2p */
        chi = (chi + Up * q1);
        Up  = log (Up);

        /* p1 g1 */
        Up = Up + 2.0 * p1;

        /* r1 g2 */
        {
          const gdouble alpha = asinh (chi);
          const gdouble gamma = Up - 0.5 * log1p (chi * chi);
          const gdouble r1    = 0.5 * ncm_csq1d_eval_int_mnu2 (csq1d, model, t);

          gdouble alphap, gammap;

          _ncm_csq1d_ct_g0g1 (alpha, gamma, -2.0 * r1, &alphap, &gammap);

          chi = sinh (alphap);
          Up  = gammap + 0.5 * log1p (chi * chi);
        }
        ncm_csq1d_state_set_up (state, NCM_CSQ1D_FRAME_NONADIAB2, t, chi, Up);
        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }
  }
}

/**
 * ncm_csq1d_evolve_prop_vector:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @initial_state: initial state
 * @frame: frame
 * @t: time $t$
 *
 * Uses the propagator to evolve the state vector @initial_state to time $t$ and
 * at frame @frame.
 *
 * Returns: (transfer none): the state vector.
 */
NcmCSQ1DState *
ncm_csq1d_evolve_prop_vector (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *initial_state, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state)
{
  _ncm_csq1d_evolve_prop_vector (csq1d, model, initial_state, frame, t, state);

  return state;
}

