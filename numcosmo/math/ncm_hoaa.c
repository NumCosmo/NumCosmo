/***************************************************************************
 *            ncm_hoaa.c
 *
 *  Fri November 04 13:27:52 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_hoaa.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_hoaa
 * @title: NcmHOAA
 * @short_description: Abstract class for Harmonic Oscillator calculation through AA variables.
 *
 * This object represents a generic time dependent harmonic oscillator for the variables
 * $q$ and its momentum $p$. The Hamiltonian of the system is given by 
 * \begin{equation}\label{eq:H}
 * H = \frac{P_q^2}{2m} + \frac{m(\nu^2 - V)q^2}{2},
 * \end{equation}
 * where the mass $m$, frequency $\nu$ and the potential $V$ are functions of the time $t$
 * and mode $k$. The mass $m$ and frequency $\nu$ are assumed to be positive definite functions.
 * For details about this system see Appendix A of [Celani et al (2016)][XCelani2016].
 * 
 * Each child should implement the functions ncm_hoaa_eval_mnu(), ncm_hoaa_eval_nu() and
 * ncm_hoaa_eval_system(). The potential term is optional, if not implemented then $V$ is assumed
 * to be zero.
 * 
 * The Action Angle variables are defined through
 * \begin{align}
 * q &= \sqrt{\frac{2I}{m\nu}}\sin\theta, \\\\
 * P_q &= \sqrt{2Im\nu}\cos\theta.
 * \end{align}
 * Therefore, the new Hamiltonian is
 * \begin{equation}\label{eq:HAA}
 * H = I\nu + \frac{\mathrm{d}\ln(\sqrt{m\nu})}{\mathrm{d}t}I\sin2\theta - \frac{V}{\nu}I\sin^2\theta.
 * \end{equation}
 * The Hamilton equations for these variables are:
 * \begin{align}
 * \frac{\mathrm{d}\theta}{\mathrm{d}t} &= \nu - \frac{V}{\nu}\sin^2\theta + \frac{\mathrm{d}\ln(\sqrt{m\nu})}{\mathrm{d}t}\sin2\theta, \\\\
 * \frac{\mathrm{d}\ln(I)}{\mathrm{d}t} &= - 2\frac{\mathrm{d}\ln(\sqrt{m\nu})}{\mathrm{d}t}\cos2\theta + \frac{V}{\nu}I\sin2\theta.
 * \end{align}
 * 
 * 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include <complex.h>

#include "math/ncm_hoaa.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_ode_spline.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_rng.h"
#include "math/ncm_diff.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#undef HAVE_SUNDIALS_ARKODE

/*
#undef NCM_HOAA_DEBUG_EVOL
#define NCM_HOAA_DEBUG_EVOL TRUE
*/

#include <nvector/nvector_serial.h>

#ifndef HAVE_SUNDIALS_ARKODE
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#else /* HAVE_SUNDIALS_ARKODE */
#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>
#endif /* HAVE_SUNDIALS_ARKODE */

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D

#include <sundials/sundials_types.h> 

#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_fit.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmHOAAPrivate
{
#ifndef HAVE_SUNDIALS_ARKODE
  gpointer cvode;
  gpointer cvode_sing;
  gpointer cvode_phase;
  gboolean cvode_init;
  gboolean cvode_sing_init;
#else
  gpointer arkode;
  gpointer arkode_sing;
  gpointer arkode_phase;
  gboolean arkode_init;
  gboolean arkode_sing_init;
#endif /* HAVE_SUNDIALS_ARKODE */
  NcmHOAAOpt opt;
  N_Vector sigma;
  N_Vector y;
  N_Vector abstol_v;
  SUNMatrix A;
  SUNLinearSolver LS;
  gdouble reltol;
  gdouble abstol;
  gdouble ti;
  gdouble tf;
  gdouble t_cur;
  gboolean save_evol;
  NcmModelCtrl *ctrl;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  NcmRNG *rng;
  gdouble sigma_c;
  gdouble tc;
  gdouble t_ad_0, t_ad_1;
  gdouble t_na_0, t_na_1;
  GArray *t, *t_m_ts, *sing_qbar, *sing_pbar, *upsilon, *gamma, *qbar, *pbar;
  NcmSpline *upsilon_s, *gamma_s, *qbar_s, *pbar_s;
  gdouble sigma0;
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
  PROP_SAVE_EVOL,
  PROP_OPT,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmHOAA, ncm_hoaa, G_TYPE_OBJECT);

static void
ncm_hoaa_init (NcmHOAA *hoaa)
{
  NcmHOAAPrivate * const self = hoaa->priv = ncm_hoaa_get_instance_private (hoaa);
  hoaa->k                      = 0.0;
  /* private properties */
#ifndef HAVE_SUNDIALS_ARKODE
  self->cvode            = CVodeCreate (CV_BDF);   /*CVodeCreate (CV_BDF, CV_NEWTON);*/
  self->cvode_sing       = CVodeCreate (CV_BDF);   /*CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);*/
  self->cvode_phase      = CVodeCreate (CV_ADAMS); /*CVodeCreate (CV_BDF, CV_NEWTON);*/
  self->cvode_init       = FALSE;
  self->cvode_sing_init  = FALSE;
#else
  self->arkode           = ARKodeCreate ();
  self->arkode_sing      = ARKodeCreate ();
  self->arkode_phase     = ARKodeCreate ();
  self->arkode_sing_init = FALSE;
#endif /* HAVE_SUNDIALS_ARKODE */
  self->opt       = NCM_HOAA_OPT_INVALID;
  self->reltol    = 0.0;
  self->abstol    = 0.0;
  self->ti        = 0.0;
  self->tf        = 0.0;
  self->t_cur     = 0.0;
  self->save_evol = FALSE;
  self->sigma     = N_VNew_Serial (1);
  self->y         = N_VNew_Serial (NCM_HOAA_VAR_SYS_SIZE);
  self->abstol_v  = N_VNew_Serial (NCM_HOAA_VAR_SYS_SIZE);
  self->A         = SUNDenseMatrix (NCM_HOAA_VAR_SYS_SIZE, NCM_HOAA_VAR_SYS_SIZE);
  self->LS        = SUNDenseLinearSolver (self->y, self->A);
  
  NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
  self->ctrl      = ncm_model_ctrl_new (NULL);
  self->T         = gsl_root_fsolver_brent;
  self->s         = gsl_root_fsolver_alloc (self->T);
  self->rng       = ncm_rng_new (NULL);

  self->t         = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->t_m_ts    = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->sing_qbar = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->sing_pbar = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->upsilon   = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->gamma     = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->qbar      = g_array_new (TRUE, TRUE, sizeof (gdouble));
  self->pbar      = g_array_new (TRUE, TRUE, sizeof (gdouble));

  self->upsilon_s = ncm_spline_cubic_notaknot_new ();
  self->gamma_s   = ncm_spline_cubic_notaknot_new ();
  self->qbar_s    = ncm_spline_cubic_notaknot_new ();
  self->pbar_s    = ncm_spline_cubic_notaknot_new ();

  self->sigma0    = 0.0;

  self->diff      = ncm_diff_new ();
  
  ncm_rng_set_random_seed (self->rng, TRUE);
}

static void
_ncm_hoaa_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmHOAA *hoaa = NCM_HOAA (object);
  NcmHOAAPrivate * const self = hoaa->priv;
  g_return_if_fail (NCM_IS_HOAA (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      ncm_hoaa_set_reltol (hoaa, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_hoaa_set_abstol (hoaa, g_value_get_double (value));
      break;
    case PROP_K:
      ncm_hoaa_set_k (hoaa, g_value_get_double (value));
      break;
    case PROP_TI:
      ncm_hoaa_set_ti (hoaa, g_value_get_double (value));
      break;
    case PROP_TF:
      ncm_hoaa_set_tf (hoaa, g_value_get_double (value));
      break;
    case PROP_SAVE_EVOL:
      ncm_hoaa_save_evol (hoaa, g_value_get_boolean (value));
      break;
    case PROP_OPT:
      self->opt = g_value_get_enum (value);
      if (self->opt == NCM_HOAA_OPT_INVALID)
        g_error ("_ncm_hoaa_set_property: `opt' property must be set.");
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_hoaa_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmHOAA *hoaa = NCM_HOAA (object);
  NcmHOAAPrivate * const self = hoaa->priv;
  g_return_if_fail (NCM_IS_HOAA (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, self->abstol);
      break;
    case PROP_K:
      g_value_set_double (value, hoaa->k);
      break;
    case PROP_TI:
      g_value_set_double (value, self->ti);
      break;
    case PROP_TF:
      g_value_set_double (value, self->tf);
      break;
    case PROP_SAVE_EVOL:
      g_value_set_boolean (value, self->save_evol);
      break;
    case PROP_OPT:
      g_value_set_enum (value, self->opt);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_hoaa_dispose (GObject *object)
{
  NcmHOAA *hoaa = NCM_HOAA (object);
  NcmHOAAPrivate * const self = hoaa->priv;

  ncm_model_ctrl_clear (&self->ctrl);
  ncm_rng_clear (&self->rng);

  ncm_spline_clear (&self->upsilon_s);
  ncm_spline_clear (&self->gamma_s);
  ncm_spline_clear (&self->qbar_s);
  ncm_spline_clear (&self->pbar_s);
  
  g_clear_pointer (&self->t,         g_array_unref);
  g_clear_pointer (&self->t_m_ts,    g_array_unref);
  g_clear_pointer (&self->sing_qbar, g_array_unref);
  g_clear_pointer (&self->sing_pbar, g_array_unref);
  g_clear_pointer (&self->upsilon,   g_array_unref);
  g_clear_pointer (&self->gamma,     g_array_unref);
  g_clear_pointer (&self->qbar,      g_array_unref);
  g_clear_pointer (&self->pbar,      g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_hoaa_parent_class)->dispose (object);
}

static void
_ncm_hoaa_finalize (GObject *object)
{
  NcmHOAA *hoaa = NCM_HOAA (object);
  NcmHOAAPrivate * const self = hoaa->priv;

#ifndef HAVE_SUNDIALS_ARKODE
  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode = NULL;
  }
  if (self->cvode_sing != NULL)
  {
    CVodeFree (&self->cvode_sing);
    self->cvode_sing = NULL;
  }
  if (self->cvode_phase != NULL)
  {
    CVodeFree (&self->cvode_phase);
    self->cvode_phase = NULL;
  }
#else
  if (self->arkode != NULL)
  {
    ARKodeFree (&self->arkode);
    self->arkode = NULL;
  }
  if (self->arkode_sing != NULL)
  {
    ARKodeFree (&self->arkode_sing);
    self->arkode_sing = NULL;
  }
  if (self->arkode_phase != NULL)
  {
    ARKodeFree (&self->arkode_phase);
    self->arkode_phase = NULL;
  }
#endif /* HAVE_SUNDIALS_ARKODE */  
  
  if (self->sigma != NULL)
  {
    N_VDestroy (self->sigma);
    self->sigma = NULL;
  }

  if (self->y != NULL)
  {
    N_VDestroy (self->y);
    self->y = NULL;
  }

  if (self->abstol_v != NULL)
  {
    N_VDestroy (self->abstol_v);
    self->abstol_v = NULL;
  }

  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }
  
  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
  }

  g_clear_pointer (&self->s, gsl_root_fsolver_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_hoaa_parent_class)->finalize (object);
}

static gdouble _ncm_hoaa_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model);

static gdouble _ncm_hoaa_eval_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)    { g_error ("_ncm_hoaa_eval_mnu: not implemented."); return 0.0; }
static gdouble _ncm_hoaa_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)     { g_error ("_ncm_hoaa_eval_nu: not implemented."); return 0.0; }
static gdouble _ncm_hoaa_eval_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k) { g_error ("_ncm_hoaa_eval_dlnmnu: not implemented."); return 0.0; }
static gdouble _ncm_hoaa_eval_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)      { g_error ("_ncm_hoaa_eval_V: not implemented."); return 0.0; }
static void _ncm_hoaa_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu) { g_error ("_ncm_hoaa_eval_system: not implemented."); }
static guint _ncm_hoaa_nsing (NcmHOAA *hoaa, NcmModel *model, const gdouble k)                          { g_error ("_ncm_hoaa_nsing: not implemented."); return 0; }
static void _ncm_hoaa_get_sing_info (NcmHOAA *hoaa, NcmModel *model, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st) { g_error ("_ncm_hoaa_get_sing_info: not implemented."); }
static gdouble _ncm_hoaa_eval_sing_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing)    { g_error ("_ncm_hoaa_eval_sing_mnu: not implemented."); return 0.0; }
static gdouble _ncm_hoaa_eval_sing_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing) { g_error ("_ncm_hoaa_eval_sing_dlnmnu: not implemented."); return 0.0; }
static gdouble _ncm_hoaa_eval_sing_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing)      { g_error ("_ncm_hoaa_eval_sing_V: not implemented."); return 0.0; }
static void _ncm_hoaa_eval_sing_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu) { { g_error ("_ncm_hoaa_eval_sing_system: not implemented."); } }
static void _ncm_hoaa_prepare (NcmHOAA *hoaa, NcmModel *model) { g_error ("_ncm_hoaa_prepare: not implemented."); }

static void
ncm_hoaa_class_init (NcmHOAAClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_hoaa_set_property;
  object_class->get_property = &_ncm_hoaa_get_property;
  object_class->dispose      = &_ncm_hoaa_dispose;
  object_class->finalize     = &_ncm_hoaa_finalize;

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
                                                        "The mode k",
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
                                   PROP_SAVE_EVOL,
                                   g_param_spec_boolean ("save-evol",
                                                         NULL,
                                                         "Save the system evolution",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_OPT,
                                   g_param_spec_enum ("opt",
                                                      NULL,
                                                      "Evolution options",
                                                      NCM_TYPE_HOAA_OPT, NCM_HOAA_OPT_INVALID,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  klass->eval_mnu         = &_ncm_hoaa_eval_mnu;
  klass->eval_nu          = &_ncm_hoaa_eval_nu;
  klass->eval_dlnmnu      = &_ncm_hoaa_eval_dlnmnu;
  klass->eval_V           = &_ncm_hoaa_eval_V;
  klass->eval_system      = &_ncm_hoaa_eval_system;

  klass->nsing            = &_ncm_hoaa_nsing;
  klass->get_sing_info    = &_ncm_hoaa_get_sing_info;
  klass->eval_sing_mnu    = &_ncm_hoaa_eval_sing_mnu;
  klass->eval_sing_dlnmnu = &_ncm_hoaa_eval_sing_dlnmnu;
  klass->eval_sing_V      = &_ncm_hoaa_eval_sing_V;
  klass->eval_sing_system = &_ncm_hoaa_eval_sing_system;

  klass->prepare          = &_ncm_hoaa_prepare;
  
	klass->eval_powspec_factor = &_ncm_hoaa_eval_powspec_factor;
}

static gdouble 
_ncm_hoaa_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model)
{
	const gdouble one_2pi2 = 1.0 / (2.0 * gsl_pow_2 (M_PI));
  return one_2pi2;
}

/**
 * ncm_hoaa_ref:
 * @hoaa: a #NcmHOAA
 *
 * Increases the reference count of @hoaa.
 *
 * Returns: (transfer full): @hoaa.
 */
NcmHOAA *
ncm_hoaa_ref (NcmHOAA *hoaa)
{
  return g_object_ref (hoaa);
}

/**
 * ncm_hoaa_free:
 * @hoaa: a #NcmHOAA
 *
 * Decreases the reference count of @hoaa.
 *
 */
void
ncm_hoaa_free (NcmHOAA *hoaa)
{
  g_object_unref (hoaa);
}

/**
 * ncm_hoaa_clear:
 * @hoaa: a #NcmHOAA
 *
 * Decreases the reference count of *@hoaa and sets the pointer *@hoaa to NULL.
 *
 */
void
ncm_hoaa_clear (NcmHOAA **hoaa)
{
  g_clear_object (hoaa);
}

/**
 * ncm_hoaa_set_reltol:
 * @hoaa: a #NcmHOAA
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance to @reltol.
 *
 */
void 
ncm_hoaa_set_reltol (NcmHOAA *hoaa, const gdouble reltol)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  self->reltol = reltol;
}

/**
 * ncm_hoaa_set_abstol:
 * @hoaa: a #NcmHOAA
 * @abstol: absolute tolerance
 *
 * Sets the absolute tolerance to @abstol.
 *
 */
void 
ncm_hoaa_set_abstol (NcmHOAA *hoaa, const gdouble abstol)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  self->abstol = abstol;
}

/**
 * ncm_hoaa_set_k:
 * @hoaa: a #NcmHOAA
 * @k: mode $k$
 *
 * Sets the mode $k$ to @k.
 *
 */
void 
ncm_hoaa_set_k (NcmHOAA *hoaa, const gdouble k)
{
  hoaa->k = k;
}

/**
 * ncm_hoaa_set_ti:
 * @hoaa: a #NcmHOAA
 * @ti: mode $t_i$
 *
 * Sets the initial time $t_i$ to @ti.
 *
 */
void 
ncm_hoaa_set_ti (NcmHOAA *hoaa, const gdouble ti)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  self->ti = ti;
}

/**
 * ncm_hoaa_set_tf:
 * @hoaa: a #NcmHOAA
 * @tf: mode $t_f$
 *
 * Sets the initial time $t_f$ to @tf.
 *
 */
void 
ncm_hoaa_set_tf (NcmHOAA *hoaa, const gdouble tf)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  self->tf = tf;
}

/**
 * ncm_hoaa_save_evol:
 * @hoaa: a #NcmHOAA
 * @save_evol: whether to save all evolution
 *
 * If true saves all evolution to be evaluted later through FIXME
 *
 */
void 
ncm_hoaa_save_evol (NcmHOAA *hoaa, gboolean save_evol)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  if (self->save_evol != save_evol)
  {
    ncm_model_ctrl_force_update (self->ctrl);    
    self->save_evol = save_evol;
  }
}

typedef struct _NcmHOAAArg
{
  NcmHOAA *hoaa;
  NcmModel *model;
  gdouble prec;
  guint sing;
  NcmHOAASingType st;
} NcmHOAAArg;

static gdouble 
_ncm_hoaa_F (gdouble t, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  return 0.5 * Vnu / nu;
}

static gdouble 
_ncm_hoaa_G (gdouble t, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);
  
  return 0.5 * dlnmnu / nu;
}

static gint
_ncm_hoaa_phase_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) f_data;
  const gdouble nu = ncm_hoaa_eval_nu (arg->hoaa, arg->model, t, arg->hoaa->k);

  NV_Ith_S (ydot, 0) = nu;

  return 0;
}

/******************************************************************************************************/
/***************************************       V and dlnmnu       *************************************/
/***************************************          START           *************************************/
/******************************************************************************************************/
#if 0
static gint
_ncm_hoaa_full_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;

  const gdouble theta      = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi        = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);

  const gdouble sin_theta  = sin (theta);
  const gdouble sin_psi    = sin (psi);
  const gdouble sin2_theta = sin_theta * sin_theta;
  const gdouble sin2_psi   = sin_psi * sin_psi;

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);

  NV_Ith_S (ydot, NCM_HOAA_VAR_THETAB)  = nu - sin2_theta * Vnu + 0.5 * dlnmnu * sin_2theta;
  NV_Ith_S (ydot, NCM_HOAA_VAR_UPSILON) = nu - sin2_psi * Vnu   + 0.5 * dlnmnu * sin_2psi;

  NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = Vnu * sin_2theta - dlnmnu * cos_2theta;
  NV_Ith_S (ydot, NCM_HOAA_VAR_LNJ)     = Vnu * sin_2psi   - dlnmnu * cos_2psi;
  
  return 0;
}

static gint
_ncm_hoaa_full_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble theta = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi   = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi, &sin_2psi, &cos_2psi);

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_THETAB) = dlnmnu * cos_2theta - Vnu * sin_2theta;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_GAMMA)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB) = 2.0 * (dlnmnu * sin_2theta + Vnu * cos_2theta);
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_GAMMA)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_UPSILON)   = dlnmnu * cos_2psi - Vnu * sin_2psi;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_LNJ)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_UPSILON)   = 2.0 * (dlnmnu * sin_2psi + Vnu * cos_2psi);
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_LNJ)   = 0.0;

  return 0;
}
/******************************************************************************************************/
/***************************************       V and dlnmnu       *************************************/
/***************************************           END            *************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/***************************************            V             *************************************/
/***************************************          START           *************************************/
/******************************************************************************************************/
static gint
_ncm_hoaa_V_only_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;

  const gdouble theta      = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi        = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);

  const gdouble sin_theta  = sin (theta);
  const gdouble sin_psi    = sin (psi);
  const gdouble sin2_theta = sin_theta * sin_theta;
  const gdouble sin2_psi   = sin_psi * sin_psi;
  const gdouble sin_2theta = sin (2.0 * theta);
  const gdouble sin_2psi   = sin (2.0 * psi);
  
  gdouble nu, dlnmnu, Vnu;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  NV_Ith_S (ydot, NCM_HOAA_VAR_THETAB) = nu - sin2_theta * Vnu;
  NV_Ith_S (ydot, NCM_HOAA_VAR_UPSILON)   = nu - sin2_psi * Vnu;

  NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = Vnu * sin_2theta;
  NV_Ith_S (ydot, NCM_HOAA_VAR_LNJ)   = Vnu * sin_2psi;
  
  return 0;
}

static gint
_ncm_hoaa_V_only_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble theta = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi   = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_THETAB) = - Vnu * sin_2theta;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_GAMMA)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB) = 2.0 * Vnu * cos_2theta;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_GAMMA)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_UPSILON)   = - Vnu * sin_2psi;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_LNJ)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_UPSILON)   = 2.0 * Vnu * cos_2psi;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_LNJ)   = 0.0;

  return 0;
}
/******************************************************************************************************/
/***************************************            V             *************************************/
/***************************************           END            *************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/***************************************         dlnmnu           *************************************/
/***************************************          START           *************************************/
/******************************************************************************************************/
#endif /* 0 */
static void
_ncm_hoaa_dlnmnu_calc_u_v (const gdouble upsilon, const gdouble lnmnu, gdouble *u, gdouble *v, gdouble *epsilon)
{
  const gdouble ln_abs_upsilon  = log (fabs (upsilon));
  const gdouble abs_lnmnu = fabs (lnmnu);
  const gdouble Tex       = ln_abs_upsilon + abs_lnmnu;
  const gdouble one_mnu2  = exp (-2.0 * abs_lnmnu);
  const gdouble lnmnu_cor = 1.0 + one_mnu2;
  const gdouble s_arg     = 0.5 * lnmnu_cor * exp (Tex);
  
  if (Tex < 0.0)
  {
    epsilon[0] = GSL_SIGN (upsilon) * asinh (s_arg);
    u[0]       = epsilon[0] + lnmnu;
    v[0]       = epsilon[0] - lnmnu;
  }
  else
  {
    const gdouble s_arg2       = s_arg * s_arg;
    const gdouble sqrt_1_1s2   = sqrt (1.0 + 1.0 / s_arg2);
    const gdouble sign_upsilon = GSL_SIGN (upsilon);
    const gdouble sign_lnmnu   = GSL_SIGN (lnmnu);

    if (sign_lnmnu * sign_upsilon > 0)
    {
      v[0]       = sign_upsilon * (ln_abs_upsilon - M_LN2 + log1p (sqrt_1_1s2 + one_mnu2 + one_mnu2 * sqrt_1_1s2));
      u[0]       = v[0] + 2.0 * lnmnu;
      epsilon[0] = v[0] + lnmnu;
    }
    else
    {
      u[0]       = sign_upsilon * (ln_abs_upsilon - M_LN2 + log1p (sqrt_1_1s2 + one_mnu2 + one_mnu2 * sqrt_1_1s2));
      v[0]       = u[0] - 2.0 * lnmnu;
      epsilon[0] = u[0] - lnmnu;
    }
  }
}

static gint
_ncm_hoaa_dlnmnu_only_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;

  const gdouble qbar    = NV_Ith_S (y, NCM_HOAA_VAR_QBAR);
  const gdouble pbar    = NV_Ith_S (y, NCM_HOAA_VAR_PBAR);
  const gdouble upsilon = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);
  const gdouble mnu     = ncm_hoaa_eval_mnu (arg->hoaa, arg->model, t, arg->hoaa->k);
  const gdouble lnmnu   = log (mnu);
  
  gdouble nu, dlnmnu, Vnu, u, v, epsilon;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);
  _ncm_hoaa_dlnmnu_calc_u_v (upsilon, lnmnu, &u, &v, &epsilon);
  
  {
    const gdouble tanh_u          = tanh (u);
    const gdouble tanh_v          = tanh (v);
    const gdouble lnch_u          = gsl_sf_lncosh (u);
    const gdouble lnch_v          = gsl_sf_lncosh (v);
    const gdouble ch_u_ch_v       = exp (- lnch_v + lnch_u);
    const gdouble qbar2           = qbar * qbar;
    const gdouble pbar2           = pbar * pbar;
    const gdouble lnch_epsilon    = gsl_sf_lncosh (epsilon);
    const gdouble lnch_lnmnu      = gsl_sf_lncosh (lnmnu);
    const gdouble one_p_ch_u_ch_v = 1.0 + ch_u_ch_v;
    const gdouble one_p_ch_v_ch_u = 1.0 + 1.0 / ch_u_ch_v;

    NV_Ith_S (ydot, NCM_HOAA_VAR_QBAR)    = + nu * pbar + dlnmnu * tanh_v * qbar / one_p_ch_u_ch_v;
    NV_Ith_S (ydot, NCM_HOAA_VAR_PBAR)    = - nu * qbar - dlnmnu * tanh_u * pbar / one_p_ch_v_ch_u;

    NV_Ith_S (ydot, NCM_HOAA_VAR_UPSILON) = - 2.0 * dlnmnu * (pbar2 / one_p_ch_v_ch_u - qbar2 / one_p_ch_u_ch_v);
    NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = - 2.0 * dlnmnu * qbar * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon);
	}

  return 0;
}

static gint
_ncm_hoaa_dlnmnu_only_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg          = (NcmHOAAArg *) jac_data;

  const gdouble qbar    = NV_Ith_S (y, NCM_HOAA_VAR_QBAR);
  const gdouble pbar    = NV_Ith_S (y, NCM_HOAA_VAR_PBAR);
  const gdouble upsilon = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);
  const gdouble mnu     = ncm_hoaa_eval_mnu (arg->hoaa, arg->model, t, arg->hoaa->k);
  const gdouble lnmnu   = log (mnu);

  gdouble nu, dlnmnu, Vnu, u, v, epsilon;
  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);
  _ncm_hoaa_dlnmnu_calc_u_v (upsilon, lnmnu, &u, &v, &epsilon);
  
  {
    const gdouble tanh_u          = tanh (u);
    const gdouble tanh_v          = tanh (v);
    const gdouble lnch_u          = gsl_sf_lncosh (u);
    const gdouble lnch_v          = gsl_sf_lncosh (v);
    const gdouble ch_u_ch_v       = exp (- lnch_v + lnch_u);
    const gdouble lnch_epsilon    = gsl_sf_lncosh (epsilon);
    const gdouble lnch_lnmnu      = gsl_sf_lncosh (lnmnu);
    const gdouble tanh_epsilon    = tanh (epsilon);
    const gdouble tanh_lnmnu      = tanh (lnmnu);
    const gdouble one_p_ch_u_ch_v = 1.0 + ch_u_ch_v;
    const gdouble one_p_ch_v_ch_u = 1.0 + 1.0 / ch_u_ch_v;

    const gdouble dtanh_epsilon_dupsilon = exp (lnch_lnmnu - 3.0 * lnch_epsilon);
    
    /* + nu * pbar + dlnmnu * tanh_v * qbar / one_p_ch_u_ch_v */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_QBAR,      NCM_HOAA_VAR_QBAR)    = + dlnmnu * tanh_v / (1.0 + ch_u_ch_v);
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_QBAR,      NCM_HOAA_VAR_PBAR)    = + nu;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_QBAR,      NCM_HOAA_VAR_UPSILON) = + 0.5 * dlnmnu * qbar * dtanh_epsilon_dupsilon;

    /* - nu * qbar - dlnmnu * tanh_u * pbar / one_p_ch_v_ch_u */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_PBAR,      NCM_HOAA_VAR_QBAR)    = - nu;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_PBAR,      NCM_HOAA_VAR_PBAR)    = - dlnmnu * tanh_u / (1.0 + 1.0 / ch_u_ch_v);
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_PBAR,      NCM_HOAA_VAR_UPSILON) = - 0.5 * dlnmnu * pbar * dtanh_epsilon_dupsilon;

    /* - 2.0 * dlnmnu * (pbar2 / one_p_ch_v_ch_u - qbar2 / one_p_ch_u_ch_v) */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_QBAR)    = + 4.0 * dlnmnu * qbar / one_p_ch_u_ch_v;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_PBAR)    = - 4.0 * dlnmnu * pbar / one_p_ch_v_ch_u;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_UPSILON) = - dlnmnu * tanh_lnmnu * exp (- 2.0 * lnch_epsilon);

    /* - 2.0 * dlnmnu * qbar * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon) */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,     NCM_HOAA_VAR_QBAR)    = - 2.0 * dlnmnu * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon);
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,     NCM_HOAA_VAR_PBAR)    = - 2.0 * dlnmnu * qbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon);
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,     NCM_HOAA_VAR_UPSILON) = + 4.0 * dlnmnu * qbar * pbar * tanh_epsilon * exp (2.0 * lnch_lnmnu - 3.0 * lnch_epsilon);    
  }

  return 0;
}

static gint
_ncm_hoaa_dlnmnu_only_sing_f (realtype t_m_ts, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;
  
  const gdouble qbar    = NV_Ith_S (y, NCM_HOAA_VAR_QBAR);
  const gdouble pbar    = NV_Ith_S (y, NCM_HOAA_VAR_PBAR);
  const gdouble upsilon = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);
  const gdouble mnu     = ncm_hoaa_eval_sing_mnu (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing);
  const gdouble lnmnu   = log (mnu);

  gdouble nu, dlnmnu, Vnu, u, v, epsilon;

  ncm_hoaa_eval_sing_system (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing, &nu, &dlnmnu, &Vnu);
  _ncm_hoaa_dlnmnu_calc_u_v (upsilon, lnmnu, &u, &v, &epsilon);

  {
    const gdouble tanh_u          = tanh (u);
    const gdouble tanh_v          = tanh (v);
    const gdouble lnch_u          = gsl_sf_lncosh (u);
    const gdouble lnch_v          = gsl_sf_lncosh (v);
    const gdouble ch_u_ch_v       = exp (- lnch_v + lnch_u);
    const gdouble qbar2           = qbar * qbar;
    const gdouble pbar2           = pbar * pbar;
    const gdouble lnch_epsilon    = gsl_sf_lncosh (epsilon);
    const gdouble lnch_lnmnu      = gsl_sf_lncosh (lnmnu);
    const gdouble one_p_ch_u_ch_v = 1.0 + ch_u_ch_v;
    const gdouble one_p_ch_v_ch_u = 1.0 + 1.0 / ch_u_ch_v;

    NV_Ith_S (ydot, NCM_HOAA_VAR_QBAR)    = + nu * pbar + dlnmnu * tanh_v * qbar / one_p_ch_u_ch_v;
    NV_Ith_S (ydot, NCM_HOAA_VAR_PBAR)    = - nu * qbar - dlnmnu * tanh_u * pbar / one_p_ch_v_ch_u;

    NV_Ith_S (ydot, NCM_HOAA_VAR_UPSILON) = - 2.0 * dlnmnu * (pbar2 / one_p_ch_v_ch_u - qbar2 / one_p_ch_u_ch_v);
    NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = - 2.0 * dlnmnu * qbar * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon); 

    if (FALSE)
    {
      printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n",
              t_m_ts,
              + nu * pbar, 
              + dlnmnu * tanh_v * qbar / one_p_ch_u_ch_v,
              - nu * qbar, 
              - dlnmnu * tanh_u * pbar / one_p_ch_v_ch_u, 
              - 2.0 * dlnmnu * pbar2 / one_p_ch_v_ch_u, 
              + 2.0 * dlnmnu * qbar2 / one_p_ch_u_ch_v, 
              - 2.0 * dlnmnu * qbar * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon),
              qbar, 
              pbar, 
              upsilon,
              NV_Ith_S (y, NCM_HOAA_VAR_GAMMA),
              nu,
              dlnmnu,
              lnmnu,
              u,
              v,
              ch_u_ch_v
              );
    }

  }

  return 0;
}

static gint
_ncm_hoaa_dlnmnu_only_sing_J (realtype t_m_ts, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg   = (NcmHOAAArg *) jac_data;

  const gdouble qbar    = NV_Ith_S (y, NCM_HOAA_VAR_QBAR);
  const gdouble pbar    = NV_Ith_S (y, NCM_HOAA_VAR_PBAR);
  const gdouble upsilon = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);
  const gdouble mnu     = ncm_hoaa_eval_sing_mnu (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing);
  const gdouble lnmnu   = log (mnu);

  gdouble nu, dlnmnu, Vnu, u, v, epsilon;

  ncm_hoaa_eval_sing_system (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing, &nu, &dlnmnu, &Vnu);
  _ncm_hoaa_dlnmnu_calc_u_v (upsilon, lnmnu, &u, &v, &epsilon);
  
  {
    const gdouble tanh_u          = tanh (u);
    const gdouble tanh_v          = tanh (v);
    const gdouble lnch_u          = gsl_sf_lncosh (u);
    const gdouble lnch_v          = gsl_sf_lncosh (v);
    const gdouble ch_u_ch_v       = exp (- lnch_v + lnch_u);
    const gdouble lnch_epsilon    = gsl_sf_lncosh (epsilon);
    const gdouble lnch_lnmnu      = gsl_sf_lncosh (lnmnu);
    const gdouble tanh_epsilon    = tanh (epsilon);
    const gdouble tanh_lnmnu      = tanh (lnmnu);
    const gdouble one_p_ch_u_ch_v = 1.0 + ch_u_ch_v;
    const gdouble one_p_ch_v_ch_u = 1.0 + 1.0 / ch_u_ch_v;

    const gdouble dtanh_epsilon_dupsilon = exp (lnch_lnmnu - 3.0 * lnch_epsilon);
    
    /* + nu * pbar + dlnmnu * tanh_v * qbar / one_p_ch_u_ch_v */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_QBAR,      NCM_HOAA_VAR_QBAR)    = + dlnmnu * tanh_v / one_p_ch_u_ch_v;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_QBAR,      NCM_HOAA_VAR_PBAR)    = + nu;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_QBAR,      NCM_HOAA_VAR_UPSILON) = + 0.5 * dlnmnu * qbar * dtanh_epsilon_dupsilon;

    /* - nu * qbar - dlnmnu * tanh_u * pbar / one_p_ch_v_ch_u */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_PBAR,      NCM_HOAA_VAR_QBAR)    = - nu;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_PBAR,      NCM_HOAA_VAR_PBAR)    = - dlnmnu * tanh_u / one_p_ch_v_ch_u;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_PBAR,      NCM_HOAA_VAR_UPSILON) = - 0.5 * dlnmnu * pbar * dtanh_epsilon_dupsilon;

    /* - 2.0 * dlnmnu * (pbar2 / one_p_ch_v_ch_u - qbar2 / one_p_ch_u_ch_v) */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_QBAR)    = + 4.0 * dlnmnu * qbar / one_p_ch_u_ch_v;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_PBAR)    = - 4.0 * dlnmnu * pbar / one_p_ch_v_ch_u;
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_UPSILON) = - dlnmnu * tanh_lnmnu * exp (- 2.0 * lnch_epsilon);

    /* - 2.0 * dlnmnu * qbar * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon) */
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,     NCM_HOAA_VAR_QBAR)    = - 2.0 * dlnmnu * pbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon);
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,     NCM_HOAA_VAR_PBAR)    = - 2.0 * dlnmnu * qbar * exp (lnch_lnmnu - 2.0 * lnch_epsilon);
    SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,     NCM_HOAA_VAR_UPSILON) = + 4.0 * dlnmnu * qbar * pbar * tanh_epsilon * exp (2.0 * lnch_lnmnu - 3.0 * lnch_epsilon);    
  }
  return 0;
}

#if 0
static gint
_ncm_hoaa_dlnmnu_only_J (glong N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble theta = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi   = NV_Ith_S (y, NCM_HOAA_VAR_UPSILON);

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_THETAB) = dlnmnu * cos_2theta;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_GAMMA)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB) = 2.0 * dlnmnu * sin_2theta;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_GAMMA)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_UPSILON)   = dlnmnu * cos_2psi;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_UPSILON,   NCM_HOAA_VAR_LNJ)   = 0.0;

  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_UPSILON)   = 2.0 * dlnmnu * sin_2psi;
  SUN_DENSE_ACCESS (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_LNJ)   = 0.0;

  return 0;
}
#endif
/******************************************************************************************************/
/***************************************         dlnmnu           *************************************/
/***************************************           END            *************************************/
/******************************************************************************************************/

void 
_ncm_hoaa_set_init_cond (NcmHOAA *hoaa, NcmModel *model, const gdouble t0)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  const gdouble Rsigma_t0 = 0.125 * M_PI; /* pi / 8 */
  gdouble Athetab, Aupsilon, Agamma;

  self->tc      = t0;
  self->sigma_c = Rsigma_t0;

  NV_Ith_S (self->sigma, 0) = Rsigma_t0;

  /* t0 must match self->tc at this point! */
  ncm_hoaa_eval_adiabatic_approx (hoaa, model, t0, &Athetab, &Aupsilon, &Agamma);

	{
    const gdouble mnu             = ncm_hoaa_eval_mnu (hoaa, model, t0, hoaa->k);
    const gdouble lnmnu           = log (mnu);
    const gdouble ch_lnmnu        = cosh (lnmnu);
    const gdouble pbar2_qbar2     = hypot (Aupsilon, 1.0 / ch_lnmnu);
    const gdouble hypot_pbar_qbar = sqrt (pbar2_qbar2);

		NV_Ith_S (self->y, NCM_HOAA_VAR_QBAR)    = hypot_pbar_qbar * sin (Athetab);
		NV_Ith_S (self->y, NCM_HOAA_VAR_PBAR)    = hypot_pbar_qbar * cos (Athetab);
		NV_Ith_S (self->y, NCM_HOAA_VAR_UPSILON) = Aupsilon;
		NV_Ith_S (self->y, NCM_HOAA_VAR_GAMMA)   = Agamma;
	}
}



void 
_ncm_hoaa_prepare_integrator (NcmHOAA *hoaa, NcmModel *model, const gdouble t0)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  gint flag;
#ifndef HAVE_SUNDIALS_ARKODE
  CVRhsFn f;
  CVLsJacFn J;
#else
  ARKRhsFn f;
  ARKLsJacFn J;
#endif /* HAVE_SUNDIALS_ARKODE */
  
  switch (self->opt)
  {
    case NCM_HOAA_OPT_FULL:
      /*f = _ncm_hoaa_full_f;*/
      /*J = _ncm_hoaa_full_J;*/
      g_assert_not_reached ();
      break;
    case NCM_HOAA_OPT_V_ONLY:
      /*f = _ncm_hoaa_V_only_f;*/
      /*J = _ncm_hoaa_V_only_J;*/
      g_assert_not_reached ();
      break;
    case NCM_HOAA_OPT_DLNMNU_ONLY:
      f = _ncm_hoaa_dlnmnu_only_f;
      J = _ncm_hoaa_dlnmnu_only_J;
      break;
    default:
      f = NULL;
      J = NULL;
      g_assert_not_reached ();
      break;
  }

  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_QBAR)    = 0.0;
  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_PBAR)    = 0.0;
  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_UPSILON) = 0.0;
  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_GAMMA)   = 0.0;

#ifndef HAVE_SUNDIALS_ARKODE
  if (!self->cvode_init)
  {
    flag = CVodeInit (self->cvode, f, t0, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->t_cur = t0;

    flag = CVodeSVtolerances (self->cvode, self->reltol, self->abstol_v);
    NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

    flag = CVodeSetMaxNumSteps (self->cvode, 0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVDlsSetLinearSolver", 1, );

    flag = CVodeSetJacFn (self->cvode, J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetJacFn", 1, );

    flag = CVodeSetInitStep (self->cvode, fabs (t0) * self->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = CVodeInit (self->cvode_phase, _ncm_hoaa_phase_f, t0, self->sigma);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSStolerances (self->cvode_phase, self->reltol, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (self->cvode_phase, 0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
    
    self->cvode_init = TRUE;
  }
  else
  {    
    flag = CVodeReInit (self->cvode, t0, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->t_cur = t0;

    flag = CVodeSetInitStep (self->cvode, fabs (t0) * self->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = CVodeReInit (self->cvode_phase, t0, self->sigma);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }  
#else
#define INTTYPE f, NULL
  if (!self->arkode_init)
  {
    flag = ARKodeInit (self->arkode, INTTYPE, t0, self->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    self->t_cur = t0;

    flag = ARKodeSVtolerances (self->arkode, self->reltol, self->abstol_v);
    NCM_CVODE_CHECK (&flag, "ARKodeSVtolerances", 1, );

    flag = ARKodeSetMaxNumSteps (self->arkode, 0);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    flag = ARKDlsSetLinearSolver (self->arkode, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetLinearSolver", 1, );

    flag = ARKDlsSetJacFn (self->arkode, J);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetJacFn", 1, );
    
    //flag = ARKodeSetLinear (self->arkode, 1);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetLinear", 1, );
    
    flag = ARKodeSetOrder (self->arkode, 7);
    NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );

    //flag = ARKodeSetERKTableNum (self->arkode, FEHLBERG_13_7_8);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetERKTableNum", 1, );
    
    flag = ARKodeSetInitStep (self->arkode, fabs (t0) * self->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = ARKodeInit (self->arkode_phase, _ncm_hoaa_phase_f,  NULL, t0, self->sigma);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSStolerances (self->arkode_phase, self->reltol, 0.0);
    NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );
    
    flag = ARKodeSetMaxNumSteps (self->arkode_phase, 1.0e9);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    //flag = ARKodeSetOrder (self->arkode_phase, 7);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );
    
    self->arkode_init = TRUE;
  }
  else
  {
    flag = ARKodeReInit (self->arkode, INTTYPE, t0, self->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    self->t_cur = t0;

    flag = ARKodeSetInitStep (self->arkode, fabs (t0) * self->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = ARKodeReInit (self->arkode_phase, _ncm_hoaa_phase_f,  NULL, t0, self->sigma);
    NCM_CVODE_CHECK (&flag, "ARKodeReInit", 1, );
  }
#endif /* HAVE_SUNDIALS_ARKODE */
}

void 
_ncm_hoaa_prepare_integrator_sing (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const guint sing, const gdouble ts, const NcmHOAASingType st)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  gint flag;
#ifndef HAVE_SUNDIALS_ARKODE
  CVRhsFn f;
  CVLsJacFn J;
#else
  ARKRhsFn f;
  ARKLsJacFn J;
#endif /* HAVE_SUNDIALS_ARKODE */
  
  /*g_assert_cmpfloat (t_m_ts + ts, ==, self->t_cur);*/
  g_assert_cmpuint (sing, <, ncm_hoaa_nsing (hoaa, model, hoaa->k));

  switch (self->opt)
  {
    case NCM_HOAA_OPT_FULL:
//      f = _ncm_hoaa_full_sing_f;
//      J = _ncm_hoaa_full_sing_J;
      break;
    case NCM_HOAA_OPT_V_ONLY:
//      f = _ncm_hoaa_V_only_sing_f;
//      J = _ncm_hoaa_V_only_sing_J;
      break;
    case NCM_HOAA_OPT_DLNMNU_ONLY:
      f = &_ncm_hoaa_dlnmnu_only_sing_f;
      J = &_ncm_hoaa_dlnmnu_only_sing_J;
      break;
    default:
      f = NULL;
      J = NULL;
      g_assert_not_reached ();
      break;
  }

  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_QBAR)    = 0.0;
  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_PBAR)    = 0.0;
  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_UPSILON) = 0.0;
  NV_Ith_S (self->abstol_v, NCM_HOAA_VAR_GAMMA)   = 0.0;
  
#ifndef HAVE_SUNDIALS_ARKODE
  if (!self->cvode_sing_init)
  {
    flag = CVodeInit (self->cvode_sing, f, t_m_ts, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSVtolerances (self->cvode_sing, self->reltol, self->abstol_v);
    NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

    flag = CVodeSetMaxNumSteps (self->cvode_sing, 0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVodeSetLinearSolver (self->cvode_sing, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVDlsSetLinearSolver", 1, );

    flag = CVodeSetJacFn (self->cvode_sing, J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetJacFn", 1, );

    flag = CVodeSetMaxErrTestFails (self->cvode_sing, 100);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxErrTestFails", 1, );

    self->cvode_sing_init = TRUE;
  }
  else
  {    
    flag = CVodeReInit (self->cvode_sing, t_m_ts, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetInitStep (self->cvode_sing, fabs (t_m_ts) * self->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
  }
#else
  if (!self->arkode_sing_init)
  {
    flag = ARKodeInit (self->arkode_sing, NULL, f, t_m_ts, self->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSVtolerances (self->arkode_sing, self->reltol, self->abstol_v);
    NCM_CVODE_CHECK (&flag, "ARKodeSVtolerances", 1, );

    flag = ARKodeSetMaxNumSteps (self->arkode_sing, 0);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    flag = ARKDlsSetLinearSolver (self->arkode_sing, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetLinearSolver", 1, );

    flag = ARKDlsSetJacFn (self->arkode_sing, J);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetJacFn", 1, );
    
    flag = ARKodeSetOrder (self->arkode_sing, 4);
    NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );

    //flag = ARKodeSetDiagnostics (self->arkode_sing, stdout);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetDiagnostics", 1, );
    
    //flag = ARKodeSetIRKTableNum (self->arkode_sing, ARK548L2SA_DIRK_8_4_5);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetIRKTableNum", 1, );

    //flag = ARKodeSetAdaptivityMethod (self->arkode_sing, 4, 1, 0, NULL);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetAdaptivityMethod", 1, );
    
    flag = ARKodeSetInitStep (self->arkode_sing, fabs (t_m_ts) * self->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

    self->arkode_sing_init = TRUE;
  }
  else
  {
    flag = ARKodeReInit (self->arkode_sing, NULL, f, t_m_ts, self->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSetInitStep (self->arkode_sing, fabs (t_m_ts) * self->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );
  }
#endif /* HAVE_SUNDIALS_ARKODE */
}

static gdouble 
_ncm_hoaa_test_tol_dlnmnu (gdouble at, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, test;
  const gdouble t = sinh (at);

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  test = dlnmnu / nu;

  return log (fabs (test / arg->prec));
}

static gdouble 
_ncm_hoaa_test_tol_Vnu (gdouble at, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, test;
  const gdouble t = sinh (at);

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  /*printf ("t=% 22.15g k=% 22.15g nu=% 22.15g dlnmnu=% 22.15g Vnu=% 22.15g\n", t, arg->hoaa->k, nu, dlnmnu, Vnu);*/
  
  test = Vnu / nu;
  
  return log (fabs (test / arg->prec));
}

static gdouble
_ncm_hoaa_search_initial_time_by_func (NcmHOAA *hoaa, NcmModel *model, gdouble tol, gdouble (*test_tol) (gdouble, gpointer))
{
  NcmHOAAPrivate * const self = hoaa->priv;
  gdouble at_hi           = asinh (self->tf);
  gdouble at_lo           = asinh (self->ti);
  gint iter               = 0;
  gint max_iter           = 100000;
  const gdouble pass_step = (at_hi - at_lo) * 1.0e-4;
  NcmHOAAArg arg          = {hoaa, model, tol, -1, NCM_HOAA_SING_TYPE_INVALID};
  
  gsl_function F;
  gint status;
  gdouble at0, test_ep;
  
  F.function = test_tol;
  F.params   = &arg;

  if ((test_ep = test_tol (at_lo, F.params)) > 0.0)
  {
    g_warning ("_ncm_hoaa_search_initial_time_by_func: system is not adiabatic at the initial time setting initial conditions at t = % 21.15g, test / tol = % 21.15g",
               sinh (at_lo), test_ep);
    return sinh (at_lo);
  }

  at_hi = at_lo + pass_step;
  while (((test_ep = test_tol (at_hi, F.params)) < 0.0) && (iter < max_iter))
  {
    at_lo  = at_hi;
    at_hi += pass_step;
    iter++;
  }

  if (iter >= max_iter)
  {
    g_warning ("_ncm_hoaa_search_initial_time_by_func: cannot find non-adiabtic regime, setting t0 to ti.");
    return self->ti;
  }

  iter = 0;
  gsl_root_fsolver_set (self->s, &F, at_lo, at_hi);
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (self->s);

    at0    = gsl_root_fsolver_root (self->s);
    at_lo  = gsl_root_fsolver_x_lower (self->s);
    at_hi  = gsl_root_fsolver_x_upper (self->s);
    status = gsl_root_test_interval (at_lo, at_hi, 0.0, self->reltol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if ((iter >= max_iter) || (status != GSL_SUCCESS))
  {
    g_warning ("_ncm_hoaa_search_initial_time_by_func: cannot find intial time with the required precision.");
  }
  
  return sinh (at0);  
}

static gdouble
_ncm_hoaa_search_initial_time (NcmHOAA *hoaa, NcmModel *model, gdouble tol)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  switch (self->opt)
  {
    case NCM_HOAA_OPT_FULL:
    {
      const gdouble t0_dlnmnu = _ncm_hoaa_search_initial_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_dlnmnu);
      const gdouble t0_Vnu    = _ncm_hoaa_search_initial_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_Vnu);
      return GSL_MIN (t0_dlnmnu, t0_Vnu);
      break;
    }
    case NCM_HOAA_OPT_V_ONLY:
    {
      return _ncm_hoaa_search_initial_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_Vnu);
      break;
    }      
    case NCM_HOAA_OPT_DLNMNU_ONLY:
    {
      return _ncm_hoaa_search_initial_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_dlnmnu);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
}

static gdouble
_ncm_hoaa_search_final_time_by_func (NcmHOAA *hoaa, NcmModel *model, gdouble tol, gdouble (*test_tol) (gdouble, gpointer))
{
  NcmHOAAPrivate * const self = hoaa->priv;
  gdouble at_hi           = asinh (self->tf);
  gdouble at_lo           = asinh (self->ti);
  gint iter               = 0;
  gint max_iter           = 100000;
  const gdouble pass_step = (at_hi - at_lo) * 1.0e-4;
  NcmHOAAArg arg          = {hoaa, model, tol, -1, NCM_HOAA_SING_TYPE_INVALID};

  
  gsl_function F;
  gint status;
  gdouble at0, test_ep;
  
  F.function = test_tol;
  F.params   = &arg;

  if ((test_ep = test_tol (at_hi, F.params)) > 0.0)
    return sinh (at_hi);

  at_lo = at_hi - pass_step;
  while (((test_ep = test_tol (at_lo, F.params)) < 0.0) && (iter < max_iter))
  {
    iter++;
    at_hi  = at_lo;
    at_lo -= pass_step;
  }

  if (iter >= max_iter)
  {
    g_warning ("_ncm_hoaa_search_final_time_by_func: cannot find non-adiabtic regime, setting t1 to tf.");
    return self->tf;
  }

  iter = 0;
  gsl_root_fsolver_set (self->s, &F, at_lo, at_hi);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (self->s);

    at0    = gsl_root_fsolver_root (self->s);
    at_lo  = gsl_root_fsolver_x_lower (self->s);
    at_hi  = gsl_root_fsolver_x_upper (self->s);
    status = gsl_root_test_interval (at_lo, at_hi, 0.0, self->reltol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if ((iter >= max_iter) || (status != GSL_SUCCESS))
  {
    g_warning ("_ncm_hoaa_search_final_time: cannot intial time with the required precision.");
  }
  
  return sinh (at0);  
}

static gdouble
_ncm_hoaa_search_final_time (NcmHOAA *hoaa, NcmModel *model, gdouble tol)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  switch (self->opt)
  {
    case NCM_HOAA_OPT_FULL:
    {
      const gdouble t1_dlnmnu = _ncm_hoaa_search_final_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_dlnmnu);
      const gdouble t1_Vnu    = _ncm_hoaa_search_final_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_Vnu);

      return GSL_MAX (t1_dlnmnu, t1_Vnu);
      break;
    }
    case NCM_HOAA_OPT_V_ONLY:
    {
      return _ncm_hoaa_search_final_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_Vnu);
      break;
    }      
    case NCM_HOAA_OPT_DLNMNU_ONLY:
    {
      return _ncm_hoaa_search_final_time_by_func (hoaa, model, tol, &_ncm_hoaa_test_tol_dlnmnu);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
}

#ifdef HAVE_SUNDIALS_ARKODE
static gint 
_ncm_hoaa_arkode_resize (N_Vector y, N_Vector ytemplate, void *user_data)
{
  return 0;
}
#endif /* HAVE_SUNDIALS_ARKODE */

void
_ncm_hoaa_evol_save (NcmHOAA *hoaa, NcmModel *model, const gdouble tf)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  NcmHOAAArg arg = {hoaa, model, 0.0, -1, NCM_HOAA_SING_TYPE_INVALID};
  gdouble last_t = -1.0e100;
  gint flag;

  if (self->sigma0 == 0.0)
    self->sigma0 = 0.125 * M_PI - hoaa->k * self->t_cur;

  g_assert_cmpfloat (tf, >, self->t_cur);

#ifndef HAVE_SUNDIALS_ARKODE
  flag = CVodeSetStopTime (self->cvode, tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetUserData (self->cvode, &arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
#else
  flag = ARKodeSetStopTime (self->arkode, tf);
  NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

  flag = ARKodeSetUserData (self->arkode, &arg);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
  
  while (TRUE)
  {
#ifndef HAVE_SUNDIALS_ARKODE
    flag = CVode (self->cvode, tf, self->y, &self->t_cur, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_hoaa_evol_save]", 1, );
#else
    flag = ARKode (self->arkode, tf, self->y, &self->t_cur, ARK_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "ARKode[_ncm_hoaa_evol_save]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
    {
      const gdouble qbar    = NV_Ith_S (self->y, NCM_HOAA_VAR_QBAR);
      const gdouble pbar    = NV_Ith_S (self->y, NCM_HOAA_VAR_PBAR);
      const gdouble upsilon = NV_Ith_S (self->y, NCM_HOAA_VAR_UPSILON);
      const gdouble gamma   = NV_Ith_S (self->y, NCM_HOAA_VAR_GAMMA);

      if (NCM_HOAA_DEBUG_EVOL)
      {
        const gdouble mnu         = ncm_hoaa_eval_mnu (hoaa, model, self->t_cur, hoaa->k);
        const gdouble lnmnu       = log (mnu);
        const gdouble ch_lnmnu    = cosh (lnmnu);
        const gdouble pbar2_qbar2 = hypot (upsilon, 1.0 / ch_lnmnu);
        gdouble q, v, Pq, Pv;

        ncm_hoaa_eval_AA2QV (hoaa, model, self->t_cur, upsilon, gamma, qbar, pbar, &q, &v, &Pq, &Pv);

        printf ("% 22.15g % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n", 
                self->t_cur, 
                upsilon, 
                gamma,
                qbar, 
                pbar,
                lnmnu,
                hypot (qbar, pbar) / sqrt (pbar2_qbar2) - 1.0,
                qbar / sqrt (pbar2_qbar2),
                pbar / sqrt (pbar2_qbar2),
                sqrt (gsl_pow_2 (pbar * pbar + qbar * qbar) - upsilon * upsilon) * ch_lnmnu
                );
      }
      
      if (self->t_cur > last_t + fabs (last_t) * NCM_HOAA_TIME_FRAC)
      {
        g_array_append_val (self->t,       self->t_cur);
        g_array_append_val (self->upsilon, upsilon);
        g_array_append_val (self->gamma,   gamma);
        g_array_append_val (self->qbar,    qbar);
        g_array_append_val (self->pbar,    pbar);

        last_t = self->t_cur;
      }          

      if (self->t_cur == tf)
        break;
    }
  }
}

static gdouble _ncm_hoaa_find_parabolic_cut_sigma (NcmVector *x_v, NcmVector *y_v, const gint n, gint i, gdouble *c0, gdouble *c1, gdouble *cov00, gdouble *cov01, gdouble *cov11, gdouble *sumsq);

static gboolean
_ncm_hoaa_find_parabolic_cut (NcmVector *x_v, NcmVector *y_v, const gdouble tol, gint *best_i, gdouble *bc0, gdouble *bc1)
{
  gdouble c0, c1, cov00, cov01, cov11, sumsq, sigma;
  const gint n         = ncm_vector_len (x_v);
  gint lb              = 0;
  gint ub              = n - NCM_HOAA_PARABOLIC_MIN_POINTS;
  const gdouble sigmab = _ncm_hoaa_find_parabolic_cut_sigma (x_v, y_v, n, ub, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

  best_i[0] = ub;
  bc0[0]    = c0;
  bc1[0]    = c1;

  if (sigmab > tol)
    return FALSE;
  
  while (TRUE)
  {
    gint i = (ub + lb) / 2;

    sigma = _ncm_hoaa_find_parabolic_cut_sigma (x_v, y_v, n, i, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    if (sigma < tol)
    {
			if (NCM_HOAA_DEBUG_SING)
				printf ("#  PASS [%d %d] => [%d %d] (% 22.15g % 22.15g)\n", lb, ub, lb, i, sigma, sigmab);

			ub = i;

      best_i[0] = i;
      bc0[0]    = c0;
      bc1[0]    = c1;
    }
    else
    {
			if (NCM_HOAA_DEBUG_SING)
				printf ("# !PASS [%d %d] => [%d %d] (% 22.15g % 22.15g)\n", lb, ub, i, ub, sigma, sigmab);
			
      lb = i;
    }

    if (ub == lb + 1)
      break;
  }

  return TRUE;
}

static gdouble
_ncm_hoaa_find_parabolic_cut_sigma (NcmVector *x_v, NcmVector *y_v, const gint n, gint i, gdouble *c0, gdouble *c1, gdouble *cov00, gdouble *cov01, gdouble *cov11, gdouble *sumsq)
{
  gdouble sigma2 = 0.0;
  gint neff      = n - i;
  gint j;

  gsl_fit_linear (ncm_vector_ptr (x_v, i), 1, ncm_vector_ptr (y_v, i), 1, neff, c0, c1, cov00, cov01, cov11, sumsq);

  for (j = 0; j < neff; j++)
  {
    sigma2 += gsl_pow_2 ((c0[0] + ncm_vector_get (x_v, i + j) * c1[0]) / ncm_vector_get (y_v, i + j) - 1.0);
  }
  
  if (NCM_HOAA_DEBUG_SING)
  {
    const gdouble ti = ncm_vector_get (x_v, i);
    printf ("# FIT: %d % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15e % 22.15e % 22.15e\n", 
            neff, ti, c0[0], c1[0], cov00[0], cov01[0], cov11[0], sumsq[0], 
            sqrt (sigma2 / neff), 
            sqrt (cov00[0]) / c0[0], sqrt (cov11[0]) / c1[0]);
  }  

  return  sqrt (sigma2 / neff);
}


void
_ncm_hoaa_evol_sing_save (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const guint sing, const gdouble ts, const NcmHOAASingType st)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  NcmHOAAArg arg         = {hoaa, model, 0.0, sing, st};
  gdouble last_t         = self->t_cur;
  gboolean check_turning = TRUE;
  gdouble tstep_m_ts;
  gint flag;

  g_assert_cmpfloat (t_m_ts + ts, >, self->t_cur);

#ifndef HAVE_SUNDIALS_ARKODE
  flag = CVodeSetStopTime (self->cvode_sing, t_m_ts);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );    

  flag = CVodeSetUserData (self->cvode_sing, &arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
#else
  flag = ARKodeSetStopTime (self->arkode_sing, t_m_ts);
  NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

  flag = ARKodeSetUserData (self->arkode_sing, &arg);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */

  g_array_set_size (self->t_m_ts,    0);
	g_array_set_size (self->sing_qbar, 0);
  g_array_set_size (self->sing_pbar, 0);
  
  while (TRUE)
  {
#ifndef HAVE_SUNDIALS_ARKODE
    flag = CVode (self->cvode_sing, t_m_ts, self->y, &tstep_m_ts, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_hoaa_evol_sing_save]", 1, );
#else
    flag = ARKode (self->arkode_sing, t_m_ts, self->y, &tstep_m_ts, ARK_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "ARKode[_ncm_hoaa_evol_sing_save]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */

    self->t_cur = ts + tstep_m_ts;

    {
      const gdouble mnu     = ncm_hoaa_eval_sing_mnu (hoaa, model, tstep_m_ts, hoaa->k, sing);
      const gdouble lnmnu   = log (mnu);
      const gdouble qbar    = NV_Ith_S (self->y, NCM_HOAA_VAR_QBAR);
      const gdouble pbar    = NV_Ith_S (self->y, NCM_HOAA_VAR_PBAR);
      const gdouble upsilon = NV_Ith_S (self->y, NCM_HOAA_VAR_UPSILON);
      const gdouble gamma   = NV_Ith_S (self->y, NCM_HOAA_VAR_GAMMA);

      if ((tstep_m_ts > 0.0) && check_turning)
      {
        gint n = self->t_m_ts->len;

        if (n > NCM_HOAA_PARABOLIC_MIN_POINTS)
        {
          const gdouble hypot_qbar_pbar = hypot (qbar, pbar);
          const gdouble sin_thetab      = qbar / hypot_qbar_pbar;
          const gdouble cos_thetab      = pbar / hypot_qbar_pbar;
          GArray *x                     = self->t_m_ts;
          GArray *y                     = NULL;
          gboolean cos_dom              = FALSE;

          if (fabs (cos_thetab) > NCM_HOAA_PARABOLIC_TRIG_ONE)
          {
            y       = self->sing_qbar;
            cos_dom = TRUE;
          }
          else if (fabs (sin_thetab) > NCM_HOAA_PARABOLIC_TRIG_ONE)
          {
            y       = self->sing_pbar;
            cos_dom = FALSE;
          }

          if (y != NULL)
          {
            NcmVector *x_v = ncm_vector_new_array (x); 
            NcmVector *y_v = ncm_vector_new_array (y); 
            gint best_i;
            gdouble bc0, bc1;
            
            ncm_vector_div (y_v, x_v);

            if (_ncm_hoaa_find_parabolic_cut (x_v, y_v, self->reltol * 1.0e3, &best_i, &bc0, &bc1))
            {
              const gdouble tfp_m_ts = GSL_MIN (-g_array_index (self->t_m_ts, gdouble, best_i), t_m_ts);

							if (NCM_HOAA_DEBUG_SING)
								printf ("# Skipping from % 22.15g to % 22.15g\n", tstep_m_ts, tfp_m_ts);
              
              if (cos_dom)
                NV_Ith_S (self->y, NCM_HOAA_VAR_QBAR) = tfp_m_ts * (tfp_m_ts * bc1 + bc0);
              else
                NV_Ith_S (self->y, NCM_HOAA_VAR_PBAR) = tfp_m_ts * (tfp_m_ts * bc1 + bc0);

              _ncm_hoaa_prepare_integrator_sing (hoaa, model, tfp_m_ts, sing, ts, st);

            }
            ncm_vector_free (y_v);
            ncm_vector_free (x_v);
          }
        }
        check_turning = FALSE;
      }

			g_array_append_val (self->t_m_ts,    tstep_m_ts);
			g_array_append_val (self->sing_qbar, qbar);
			g_array_append_val (self->sing_pbar, pbar);

      if (self->t_cur > last_t + fabs (last_t) * NCM_HOAA_TIME_FRAC)
      {
        g_array_append_val (self->t,       self->t_cur);
        g_array_append_val (self->upsilon, upsilon);
        g_array_append_val (self->gamma,   gamma);
        g_array_append_val (self->qbar,    qbar);
        g_array_append_val (self->pbar,    pbar);
        
        last_t = self->t_cur;
      }          

			if (NCM_HOAA_DEBUG_EVOL_SING)
			{
        const gdouble ch_lnmnu    = cosh (lnmnu);
        const gdouble pbar2_qbar2 = hypot (upsilon, 1.0 / ch_lnmnu);
        gdouble q, v, Pq, Pv;

        gdouble nu, dlnmnu, Vnu;
        ncm_hoaa_eval_sing_system (hoaa, model, t_m_ts, hoaa->k, sing, &nu, &dlnmnu, &Vnu);

        ncm_hoaa_eval_AA2QV (hoaa, model, self->t_cur, upsilon, gamma, qbar, pbar, &q, &v, &Pq, &Pv);

        printf ("% 22.15g[% 22.15g] % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n", 
                self->t_cur, 
                tstep_m_ts,
                upsilon, 
                gamma,
                qbar, 
                pbar,
                lnmnu,
                hypot (qbar, pbar) / sqrt (pbar2_qbar2) - 1.0,
                qbar / sqrt (pbar2_qbar2),
                pbar / sqrt (pbar2_qbar2),
                sqrt (gsl_pow_2 (pbar * pbar + qbar * qbar) - upsilon * upsilon) * ch_lnmnu
                );
      }
			
      if (tstep_m_ts == t_m_ts)
        break;
    }
  }

	g_array_set_size (self->t_m_ts,    0);
	g_array_set_size (self->sing_qbar, 0);
  g_array_set_size (self->sing_pbar, 0);
}

gdouble 
_ncm_hoaa_phase_nu (gdouble y, gdouble x, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  const gdouble nu = ncm_hoaa_eval_nu (arg->hoaa, arg->model, x, arg->hoaa->k);
  return nu;
}

gdouble 
_ncm_hoaa_phase_mnu2 (gdouble y, gdouble x, gpointer userdata)
{
  NcmHOAAArg *arg   = (NcmHOAAArg *) userdata;
  const gdouble nu  = ncm_hoaa_eval_nu (arg->hoaa, arg->model, x, arg->hoaa->k);
  const gdouble mnu = ncm_hoaa_eval_mnu (arg->hoaa, arg->model, x, arg->hoaa->k);

  return GSL_MAX (mnu * nu, 1.0e100);
}

/**
 * ncm_hoaa_prepare: (virtual prepare)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 *
 * Prepares the object using @model.
 *
 */
void 
ncm_hoaa_prepare (NcmHOAA *hoaa, NcmModel *model)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  
  if (NCM_HOAA_GET_CLASS (hoaa)->prepare != NULL)
  {
    NCM_HOAA_GET_CLASS (hoaa)->prepare (hoaa, model);
  }

  g_assert_cmpfloat (self->tf, >, self->ti);

  {
    const gdouble tol     = GSL_MAX (self->reltol, 1.0e-4);
    const gdouble t_ad_0  = _ncm_hoaa_search_initial_time (hoaa, model, tol);
    const gdouble t_ad_1  = _ncm_hoaa_search_final_time (hoaa, model, tol);
    const gdouble t_na_0  = _ncm_hoaa_search_initial_time (hoaa, model, 1.0);
    const gdouble t_na_1  = _ncm_hoaa_search_final_time (hoaa, model, 1.0);

    self->t_ad_0 = t_ad_0;
    self->t_ad_1 = t_ad_1;
    self->t_na_0 = t_na_0;
    self->t_na_1 = t_na_1;

    if (FALSE)
    {
      NcmHOAAArg arg     = {hoaa, model, 0.0, -1, NCM_HOAA_SING_TYPE_INVALID};
      NcmSpline *s       = ncm_spline_cubic_notaknot_new ();
      NcmOdeSpline *nu_s = ncm_ode_spline_new (s, _ncm_hoaa_phase_nu);
      gint i;

      ncm_ode_spline_set_xi (nu_s, -1.0);
      ncm_ode_spline_set_xf (nu_s, +2.0);

      ncm_ode_spline_set_yi (nu_s, 0.0);

      ncm_ode_spline_set_reltol (nu_s, 1.0e-13);
      ncm_ode_spline_set_abstol (nu_s, 1.0e-60);
      
      ncm_ode_spline_prepare (nu_s, &arg);

      printf ("# GO!\n");
      for (i = 0; i < 10000; i++)      
      {
        const gdouble t = -1.0 + 2.0 / (10000.0 - 1.0) * i;
        printf ("% 22.15g % 22.15g\n", t, ncm_spline_eval (ncm_ode_spline_peek_spline (nu_s), t));
      } 
      printf ("\n\n");
    }

    if (FALSE)
    {
      NcmHOAAArg arg       = {hoaa, model, 0.0, -1, NCM_HOAA_SING_TYPE_INVALID};
      NcmSpline *s         = ncm_spline_cubic_notaknot_new ();
      NcmOdeSpline *mnu2_s = ncm_ode_spline_new (s, _ncm_hoaa_phase_mnu2);
      gint i;

      ncm_ode_spline_set_xi (mnu2_s, -1.0);
      ncm_ode_spline_set_xf (mnu2_s, +2.0);

      ncm_ode_spline_set_yi (mnu2_s, 0.0);

      ncm_ode_spline_set_reltol (mnu2_s, 1.0e-13);
      ncm_ode_spline_set_abstol (mnu2_s, 1.0e-70);

      ncm_ode_spline_prepare (mnu2_s, &arg);

      printf ("# GO!\n");
      for (i = 0; i < 10000; i++)      
      {
        const gdouble t = -1.0 + 2.0 / (10000.0 - 1.0) * i;
        printf ("% 22.15g % 22.15g\n", t, ncm_spline_eval (ncm_ode_spline_peek_spline (mnu2_s), t));
      } 
      printf ("\n\n");
    }
    
    _ncm_hoaa_set_init_cond (hoaa, model, t_ad_0);
    
    _ncm_hoaa_prepare_integrator (hoaa, model, t_ad_0);

    if (self->save_evol)
    {
      const guint nsing = ncm_hoaa_nsing (hoaa, model, hoaa->k); 
      if (nsing == 0)
      {
        _ncm_hoaa_evol_save (hoaa, model, t_ad_1);
      }
      else
      {
        guint i;
        gdouble *ts_a = g_new (gdouble, nsing);
        gsize *sing_a = g_new (gsize, nsing);
        
        for (i = 0; i < nsing; i++)
        {
          gdouble dts_i, dts_f;
          NcmHOAASingType st;
          
          ncm_hoaa_get_sing_info (hoaa, model, hoaa->k, i, &ts_a[i], &dts_i, &dts_f, &st);
        }

        gsl_sort_index (sing_a, ts_a, 1, nsing);
        
        for (i = 0; i < nsing; i++)
        {
          gdouble ts, dts_i, dts_f;
          NcmHOAASingType st;

          ncm_hoaa_get_sing_info (hoaa, model, hoaa->k, sing_a[i], &ts, &dts_i, &dts_f, &st);

          if (ts < self->t_cur)
            continue;

					if (ts > t_ad_1)
						break;

          _ncm_hoaa_evol_save (hoaa, model, dts_i + ts);
          _ncm_hoaa_prepare_integrator_sing (hoaa, model, dts_i, sing_a[i], ts, st);

          _ncm_hoaa_evol_sing_save (hoaa, model, dts_f, sing_a[i], ts, st);

          _ncm_hoaa_prepare_integrator (hoaa, model, dts_f + ts);
        }
        _ncm_hoaa_evol_save (hoaa, model, t_ad_1);

        g_free (ts_a);
        g_free (sing_a);
      }

      ncm_spline_set_array (self->upsilon_s, self->t, self->upsilon, TRUE);
      ncm_spline_set_array (self->gamma_s,   self->t, self->gamma,   TRUE);
      ncm_spline_set_array (self->qbar_s,    self->t, self->qbar,    TRUE);
      ncm_spline_set_array (self->pbar_s,    self->t, self->pbar,    TRUE);
    }
  }
}

/**
 * ncm_hoaa_get_t0_t1:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t0: (out): $t_0$
 * @t1: (out): $t_1$
 * 
 * Gets the time interval where the numerical evolution was calculated;
 *
 */
void 
ncm_hoaa_get_t0_t1 (NcmHOAA *hoaa, NcmModel *model, gdouble *t0, gdouble *t1)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  t0[0] = self->t_ad_0;
  t1[0] = self->t_ad_1;
}

typedef struct _NcmHOAADAngle
{
  gdouble d1;
  gdouble d2;
  gdouble d3;
  gdouble d4;
} NcmHOAADAngle;

static gdouble
_ncm_hoaa_dangle_f (gdouble dt, gpointer userdata)
{
  NcmHOAADAngle *da = (NcmHOAADAngle *) userdata;
  gdouble sin_2dt, cos_2dt, sin_4dt, cos_4dt;
  
  sincos (2.0 * dt, &sin_2dt, &cos_2dt);
  sincos (4.0 * dt, &sin_4dt, &cos_4dt);

  return (da->d1 * sin_2dt + da->d2 * cos_2dt + da->d3 * sin_4dt + da->d4 * cos_4dt - dt);
}

gdouble 
_ncm_hoaa_dangle (NcmHOAA *hoaa, const gdouble R, const gdouble F, const gdouble G, const gdouble dF, const gdouble dG, const gdouble nu)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  const gdouble sin_2R = sin (2.0 * R);
  const gdouble cos_2R = cos (2.0 * R);
  const gdouble sin_4R = sin (4.0 * R);
  const gdouble cos_4R = cos (4.0 * R);
  const gdouble F2     = F * F;
  const gdouble G2     = G * G;
  const gdouble c1 = + 0.5 * (F + F2 + 0.5 * dG / nu);
  const gdouble c2 = - 0.5 * (G + F * G - 0.5 * dF / nu);
  const gdouble c3 = + 0.125 * (G2 - F2);
  const gdouble c4 = 0.250 * F * G;

  const gdouble d1 = (c1 * cos_2R - c2 * sin_2R);
  const gdouble d2 = (c1 * sin_2R + c2 * cos_2R);

  const gdouble d3 = (c3 * cos_4R - c4 * sin_4R);
  const gdouble d4 = (c3 * sin_4R + c4 * cos_4R);

  NcmHOAADAngle da = {d1, d2, d3, d4};

  gdouble dt0   = -0.25 * M_PI, dt1 = 0.25 * M_PI, dt = 0.0;
  gint max_iter = 100000;
  gint iter;

  gsl_function gslF;
  gint status;
    
  gslF.function = _ncm_hoaa_dangle_f;
  gslF.params   = &da;

  iter = 0;
  gsl_root_fsolver_set (self->s, &gslF, dt0, dt1);
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (self->s);

    dt     = gsl_root_fsolver_root (self->s);
    dt0    = gsl_root_fsolver_x_lower (self->s);
    dt1    = gsl_root_fsolver_x_upper (self->s);
    status = gsl_root_test_interval (dt0, dt1, 0.0, self->reltol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if ((iter >= max_iter) || (status != GSL_SUCCESS))
  {
    g_warning ("_ncm_hoaa_dangle: cannot adiabatic angle.");
  }

  return dt;
}

/**
 * ncm_hoaa_eval_adiabatic_approx:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @thetab: (out): $\theta_b$
 * @upsilon: (out): $\upsilon$
 * @gamma: (out): $\ln(I)$
 * 
 * Calculates the adiabatic approximation at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_adiabatic_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *thetab, gdouble *upsilon, gdouble *gamma)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  if (ncm_cmp (t, self->tc, self->reltol, 0.0) != 0)
  {
    gint flag;
    gdouble t_step;
    guint restarts = 0;
    NcmHOAAArg arg = {hoaa, model, 0.0, -1, NCM_HOAA_SING_TYPE_INVALID};

#ifndef HAVE_SUNDIALS_ARKODE
    flag = CVodeSetStopTime (self->cvode_phase, t);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );    

    flag = CVodeSetUserData (self->cvode_phase, &arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    while (TRUE) 
    {
      flag = CVode (self->cvode_phase, t, self->sigma, &t_step, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode[phase]", 1, );

      if (flag == CV_ROOT_RETURN)
      {
        NV_Ith_S (self->sigma, 0) = gsl_sf_angle_restrict_pos (NV_Ith_S (self->sigma, 0));

        flag = CVodeReInit (self->cvode_phase, t_step, self->sigma);
        NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
        restarts++;
      }

      if (t == t_step)
        break;
    }
#else
    flag = ARKodeSetStopTime (self->arkode_phase, t);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );    

    flag = ARKodeSetUserData (self->arkode_phase, &arg);
    NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

    while (TRUE) 
    {
      flag = ARKode (self->arkode_phase, t, self->sigma, &t_step, ARK_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "ARKode[phase]", 1, );

      if (flag == ARK_ROOT_RETURN)
      {
        NV_Ith_S (self->sigma, 0) = gsl_sf_angle_restrict_symm (NV_Ith_S (self->sigma, 0));

        flag = ARKodeResize (self->arkode_phase, self->sigma, 1.0e-5, t_step, _ncm_hoaa_arkode_resize, NULL);
        NCM_CVODE_CHECK (&flag, "ARKodeResize", 1, );
        restarts++;
      }

      if (t == t_step)
        break;
    }
#endif /* HAVE_SUNDIALS_ARKODE */
    self->sigma_c = gsl_sf_angle_restrict_symm (NV_Ith_S (self->sigma, 0));
  }
  
  {
    const gdouble Rsigma_t0 = self->sigma_c;
    const gdouble Isigma_t0 = self->sigma_c + 0.5 * M_PI;
    const gdouble mnu       = ncm_hoaa_eval_mnu (hoaa, model, t, hoaa->k);
    const gdouble lnmnu     = log (mnu);
    gdouble nu, dlnmnu, Vnu;
    NcmHOAAArg arg;

    ncm_hoaa_eval_system (hoaa, model, t, hoaa->k, &nu, &dlnmnu, &Vnu);

    arg.model = model;
    arg.hoaa  = hoaa;
    arg.prec  = 0.0;

    {
      gdouble err;
      const gdouble F  = 0.5 * Vnu / nu;
      const gdouble G  = 0.5 * dlnmnu / nu;
      const gdouble dF = ncm_diff_rf_d1_1_to_1 (self->diff, t, _ncm_hoaa_F, &arg, &err);
      const gdouble dG = ncm_diff_rf_d1_1_to_1 (self->diff, t, _ncm_hoaa_G, &arg, &err); 
      const gdouble F2 = F * F;
      const gdouble G2 = G * G;

      gdouble dtheta = _ncm_hoaa_dangle (hoaa, Rsigma_t0, F, G, dF, dG, nu);
      gdouble dpsi   = _ncm_hoaa_dangle (hoaa, Isigma_t0, F, G, dF, dG, nu);
      gdouble sin_2theta, cos_2theta, sin_2psi, cos_2psi;
      gdouble sin_4theta, cos_4theta, sin_4psi, cos_4psi;
      gdouble theta, psi, LnI, LnJ;

      if (fabs (dtheta) > 0.1 || fabs (dpsi) > 0.1)
      {
        sincos (2.0 * Rsigma_t0, &sin_2theta, &cos_2theta);
        sincos (4.0 * Rsigma_t0, &sin_4theta, &cos_4theta);
        sincos (2.0 * Isigma_t0, &sin_2psi, &cos_2psi);
        sincos (4.0 * Isigma_t0, &sin_4psi, &cos_4psi);
        
        theta  = Rsigma_t0 + 0.5 * (sin_2theta * (F + F2 + 0.5 * dG / nu) - cos_2theta * (G + F * G - 0.5 * dF / nu));
        psi    = Isigma_t0 + 0.5 * (sin_2psi *   (F + F2 + 0.5 * dG / nu) - cos_2psi *   (G + F * G - 0.5 * dF / nu));
        theta += 0.125 * (sin_4theta * (G2 - F2) + cos_4theta * 2.0 * F * G); 
        psi   += 0.125 * (sin_4psi *   (G2 - F2) + cos_4psi *   2.0 * F * G);

        sincos (2.0 * theta, &sin_2theta, &cos_2theta);
        sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);
        sincos (4.0 * theta, &sin_4theta, &cos_4theta);
        sincos (4.0 * psi,   &sin_4psi,   &cos_4psi);
        
        dtheta  = 0.5 * (sin_2theta * (F + F2 + 0.5 * dG / nu) - cos_2theta * (G + F * G - 0.5 * dF / nu));
        dtheta += 0.125 * (sin_4theta * (G2 - F2) + cos_4theta * 2.0 * F * G);

        dpsi    = 0.5 * (sin_2psi *   (F + F2 + 0.5 * dG / nu) - cos_2psi *   (G + F * G - 0.5 * dF / nu));
        dpsi   += 0.125 * (sin_4psi *   (G2 - F2) + cos_4psi *   2.0 * F * G);

        theta = Rsigma_t0 + dtheta;
        psi   = Isigma_t0 + dpsi;
      }
      else
      {
        theta = Rsigma_t0 + dtheta;
        psi   = Isigma_t0 + dpsi;
        
        sincos (2.0 * theta, &sin_2theta, &cos_2theta);
        sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);
        sincos (4.0 * theta, &sin_4theta, &cos_4theta);
        sincos (4.0 * psi,   &sin_4psi,   &cos_4psi);
      }
      
      LnI  = - 0.25 * (F2 + G2) - (sin_2theta * (G + F * G - 0.5 * dF / nu) + cos_2theta * (F + F2 + 0.5 * dG / nu));
      LnJ  = - 0.25 * (F2 + G2) - (sin_2psi *   (G + F * G - 0.5 * dF / nu) + cos_2psi *   (F + F2 + 0.5 * dG / nu));
      LnI += 0.25 * (sin_4theta * 2.0 * F * G + cos_4theta * (F2 - G2));
      LnJ += 0.25 * (sin_4psi *   2.0 * F * G + cos_4psi *   (F2 - G2));

      thetab[0]  = 0.5 * (Rsigma_t0 + Isigma_t0 + dtheta + dpsi);
      upsilon[0] = tan (dtheta - dpsi) / cosh (lnmnu);
      gamma[0]   = 0.5 * (LnI - LnJ);
    }
  }
}

/**
 * ncm_hoaa_eval_adiabatic_LnI_approx:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @theta: $\theta$
 * @psi: $\psi$
 * @LnI: (out): $\ln(I)$
 * @LnJ: (out): $\ln(J)$
 * 
 * Calculates the adiabatic approximation at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_adiabatic_LnI_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble theta, const gdouble psi, gdouble *LnI, gdouble *LnJ)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  gdouble nu, dlnmnu, Vnu;
  NcmHOAAArg arg;

  ncm_hoaa_eval_system (hoaa, model, t, hoaa->k, &nu, &dlnmnu, &Vnu);

  arg.model = model;
  arg.hoaa  = hoaa;
  arg.prec  = 0.0;

  {
    gdouble err;
    const gdouble F  = 0.5 * Vnu / nu;
    const gdouble G  = 0.5 * dlnmnu / nu;
    const gdouble dF = ncm_diff_rf_d1_1_to_1 (self->diff, t, _ncm_hoaa_F, &arg, &err);
    const gdouble dG = ncm_diff_rf_d1_1_to_1 (self->diff, t, _ncm_hoaa_G, &arg, &err);  
    const gdouble F2 = F * F;
    const gdouble G2 = G * G;
    gdouble sin_2theta, cos_2theta, sin_2psi, cos_2psi;
    gdouble sin_4theta, cos_4theta, sin_4psi, cos_4psi;

    sincos (2.0 * theta, &sin_2theta, &cos_2theta);
    sincos (4.0 * theta, &sin_4theta, &cos_4theta);
    sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);
    sincos (4.0 * psi,   &sin_4psi,   &cos_4psi);

    LnI[0]  = - 0.25 * (F2 + G2) - (sin_2theta * (G + F * G - 0.5 * dF / nu) + cos_2theta * (F + F2 + 0.5 * dG / nu));
    LnJ[0]  = - 0.25 * (F2 + G2) - (sin_2psi *   (G + F * G - 0.5 * dF / nu) + cos_2psi *   (F + F2 + 0.5 * dG / nu));
    LnI[0] += 0.25 * (sin_4theta * 2.0 * F * G + cos_4theta * (F2 - G2));
    LnJ[0] += 0.25 * (sin_4psi *   2.0 * F * G + cos_4psi *   (F2 - G2));
  }
}

/**
 * ncm_hoaa_eval_AA:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @upsilon: (out): $\theta$
 * @gamma: (out): $\psi$
 * @qbar: (out): $\ln(I)$
 * @pbar: (out): $\ln(J)$
 * 
 * Calculates the AA variables at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *upsilon, gdouble *gamma, gdouble *qbar, gdouble *pbar)
{
  NcmHOAAPrivate * const self = hoaa->priv;
  if (!self->save_evol)
  {
    g_error ("ncm_hoaa_eval_AA: save_evol must be enables.");
  }
  else if (t < self->t_ad_0)
  {
    gdouble Athetab, Aupsilon, Agamma;

    ncm_hoaa_eval_adiabatic_approx (hoaa, model, t, &Athetab, &Aupsilon, &Agamma);
    sincos (Athetab, qbar, pbar);

    {
      const gdouble mnu             = ncm_hoaa_eval_mnu (hoaa, model, t, hoaa->k);
      const gdouble lnmnu           = log (mnu);
      const gdouble ch_lnmnu        = cosh (lnmnu);
      const gdouble pbar2_qbar2     = hypot (Aupsilon, 1.0 / ch_lnmnu);
      const gdouble hypot_qbar_pbar = sqrt (pbar2_qbar2);

      qbar[0]   *= hypot_qbar_pbar;
      pbar[0]   *= hypot_qbar_pbar;

      upsilon[0] = Aupsilon;
      gamma[0]   = Agamma;
    }
  }
  else
  {
    if (t > self->t_ad_1)
      g_warning ("ncm_hoaa_eval_AA: extrapolating solution, max time = % 22.15g, t = % 22.15g.", self->t_ad_1, t);
    
    upsilon[0] = ncm_spline_eval (self->upsilon_s, t);
    gamma[0]   = ncm_spline_eval (self->gamma_s, t);
    qbar[0]    = ncm_spline_eval (self->qbar_s, t);
    pbar[0]    = ncm_spline_eval (self->pbar_s, t);
  }
}

/**
 * ncm_hoaa_eval_AA2QV:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @upsilon: $\upsilon$
 * @gamma: $\gamma$
 * @qbar: $\sin (\theta_b)$
 * @pbar: $\cos (\theta_b)$
 * @q: (out): $u$
 * @v: (out): $v$
 * @Pq: (out): $P_q$
 * @Pv: (out): $P_v$
 * 
 * Change the variables from AA to QV.
 * 
 */
void 
ncm_hoaa_eval_AA2QV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble upsilon, const gdouble gamma, const gdouble qbar, const gdouble pbar, gdouble *q, gdouble *v, gdouble *Pq, gdouble *Pv)
{
  const gdouble mnu              = ncm_hoaa_eval_mnu (hoaa, model, t, hoaa->k);
  const gdouble lnmnu            = log (mnu);
  const gdouble ch_lnmnu         = cosh (lnmnu);
  const gdouble pbar2_qbar2      = hypot (upsilon, 1.0 / ch_lnmnu);
  const gdouble sqrt_pbar2_qbar2 = sqrt (pbar2_qbar2);
  gdouble U, V, epsilon;

  _ncm_hoaa_dlnmnu_calc_u_v (upsilon, lnmnu, &U, &V, &epsilon);

  q[0]  = (- exp (0.5 * (+ gamma - U)) * pbar + exp (0.5 * (+ gamma + V)) * qbar) / sqrt_pbar2_qbar2;
  v[0]  = (+ exp (0.5 * (- gamma - U)) * pbar + exp (0.5 * (- gamma + V)) * qbar) / sqrt_pbar2_qbar2;

  Pq[0] = (+ exp (0.5 * (+ gamma + U)) * pbar + exp (0.5 * (+ gamma - V)) * qbar) / sqrt_pbar2_qbar2;
  Pv[0] = (+ exp (0.5 * (- gamma + U)) * pbar - exp (0.5 * (- gamma - V)) * qbar) / sqrt_pbar2_qbar2;
}

/**
 * ncm_hoaa_eval_QV2AA:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @q: $q$
 * @v: $v$
 * @Pq: $P_q$
 * @Pv: $P_v$
 * @upsilon: (out): $\upsilon$
 * @gamma: (out): $\gamma$
 * @qbar: (out): $\sin (\theta_b)$
 * @pbar: (out): $\cos (\theta_b)$
 * 
 * Change the variables from complex to AA.
 * 
 */
void 
ncm_hoaa_eval_QV2AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble q, const gdouble v, const gdouble Pq, const gdouble Pv, gdouble *upsilon, gdouble *gamma, gdouble *qbar, gdouble *pbar)
{
  const gdouble mnu      = ncm_hoaa_eval_mnu (hoaa, model, t, hoaa->k);
  const gdouble lnmnu    = log (mnu);

  const gdouble theta    = atan2 (q, Pq);
  const gdouble psi      = atan2 (v, Pv);
  const gdouble LnI      = log (0.5 * (mnu * q * q + Pq * Pq / mnu));
  const gdouble LnJ      = log (0.5 * (mnu * v * v + Pv * Pv / mnu));
  const gdouble thetab   = 0.5 * (psi + theta);
  const gdouble ch_lnmnu = cosh (lnmnu);

  upsilon[0] = GSL_SIGN (psi - theta) * sinh (acosh (exp (0.5 * (LnI + LnJ)))) / ch_lnmnu;
  gamma[0]   = 0.5 * (LnI - LnJ);
  
  sincos (thetab, qbar, pbar);

  {
    const gdouble pbar2_qbar2     = hypot (upsilon[0], 1.0 / ch_lnmnu);
    const gdouble hypot_qbar_pbar = sqrt (pbar2_qbar2);

    qbar[0] *= hypot_qbar_pbar;
    pbar[0] *= hypot_qbar_pbar;
  }
}

/**
 * ncm_hoaa_eval_QV:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @q: (out): $q$
 * @v: (out): $v$
 * @Pq: (out): $P_q$
 * @Pv: (out): $P_v$
 * 
 * Calculates the complex variables at $t$ and $k$.
 * 
 */
void 
ncm_hoaa_eval_QV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *q, gdouble *v, gdouble *Pq, gdouble *Pv)
{
  gdouble upsilon, gamma, qbar, pbar;
  ncm_hoaa_eval_AA (hoaa, model, t, &upsilon, &gamma, &qbar, &pbar);
  ncm_hoaa_eval_AA2QV (hoaa, model, t, upsilon, gamma, qbar, pbar, q, v, Pq, Pv);
}

/**
 * ncm_hoaa_eval_Delta:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @Delta_phi: (out): $\Delta_\phi$
 * @Delta_Pphi: (out): $\Delta_{P_\phi}$
 * 
 * Calculates the power spectra.
 *
 */
void 
ncm_hoaa_eval_Delta (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *Delta_phi, gdouble *Delta_Pphi)
{
  const gdouble factor  = ncm_hoaa_eval_powspec_factor (hoaa, model);
  const gdouble k3_2pi2 = gsl_pow_3 (hoaa->k) * factor;
  gdouble q = 0.0, v = 0.0, Pq = 0.0, Pv = 0.0;
    
  ncm_hoaa_eval_QV (hoaa, model, t, &q, &v, &Pq, &Pv);
  
  Delta_phi[0]  = k3_2pi2 * 0.25 * (q  * q  + v  * v);
  Delta_Pphi[0] = k3_2pi2 * 0.25 * (Pq * Pq + Pv * Pv);
}

/**
 * ncm_hoaa_eval_solution:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @S: $S$
 * @PS: $P_S$
 * @Aq: (out): $A_q$
 * @Av: (out): $A_v$
 * 
 * Calculates the coefficients $A_q$ and $A_v$ of the solution with 
 * initial conditions $S,\;P_S$ at $t$.
 *
 */
void 
ncm_hoaa_eval_solution (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble S, const gdouble PS, gdouble *Aq, gdouble *Av)
{
  gdouble q = 0.0, v = 0.0, Pq = 0.0, Pv = 0.0;
  ncm_hoaa_eval_QV (hoaa, model, t, &q, &v, &Pq, &Pv);

  Aq[0] = 0.5 * (PS * v - Pv * S);
  Av[0] = 0.5 * (Pq * S - PS * q);
}

/**
 * ncm_hoaa_eval_nu: (virtual eval_nu)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the frequency term at $t$ and $k$.
 *
 * Returns: the frequency $\nu$.
 */
/**
 * ncm_hoaa_eval_mnu: (virtual eval_mnu)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the mass term at $t$ and $k$.
 *
 * Returns: the mass $m$.
 */
/**
 * ncm_hoaa_eval_dlnmnu: (virtual eval_dlnmnu)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$ term at $t$ and $k$.
 *
 * Returns: the potential $\mathrm{d}(m\nu)/\mathrm{d}t$.
 */
/**
 * ncm_hoaa_eval_V: (virtual eval_V)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the potential term at $t$ and $k$.
 *
 * Returns: the potential $V$.
 */
/**
 * ncm_hoaa_eval_system: (virtual eval_system)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 * @nu: (out): the frequency $\nu$
 * @dlnmnu: (out): the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$
 * @Vnu: (out): the potential over frequency term $V/\nu$
 *
 * Evaluates the system functions at $t$ and $k$.
 *
 */

/**
 * ncm_hoaa_nsing: (virtual nsing)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @k: mode $k$
 *
 * Gets the number of singular points $m(t_s) = 0$ for the problem in hand.
 *
 * Returns: the number of singular points.
 */
/**
 * ncm_hoaa_get_sing_info: (virtual get_sing_info)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @k: mode $k$
 * @sing: singularity index
 * @ts: (out): singularity time $t_s$
 * @dts_i: (out): FIXME
 * @dts_f: (out): FIXME
 * @st: (out): FIXME
 *
 * Gets the time $t_s$ where the @sing-th singularity occour.
 *
 */
/**
 * ncm_hoaa_eval_sing_mnu: (virtual eval_sing_mnu)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 * @sing: singularity index
 *
 * Evaluates the mass term at $t$ and $k$.
 *
 * Returns: the mass $m$.
 */
/**
 * ncm_hoaa_eval_sing_V: (virtual eval_sing_V)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 * @sing: singularity index
 *
 * Evaluates the potential term at $t$ and $k$.
 *
 * Returns: the potential $V$.
 */
/**
 * ncm_hoaa_eval_sing_dlnmnu: (virtual eval_sing_dlnmnu)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 * @sing: singularity index
 *
 * Evaluates the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$ term at $t$ and $k$.
 *
 * Returns: the potential $\mathrm{d}(m\nu)/\mathrm{d}t$.
 */
/**
 * ncm_hoaa_eval_sing_system: (virtual eval_sing_system)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 * @sing: singularity index
 * @nu: (out): the frequency $\nu$
 * @dlnmnu: (out): the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$
 * @Vnu: (out): the potential over frequency term $V/\nu$
 *
 * Evaluates the system functions at $t$ and $k$.
 *
 */
/**
 * ncm_hoaa_eval_powspec_factor: (virtual eval_powspec_factor)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * 
 * Evaluates the power spectrum constant factor.
 * Default is: $1/(2\pi^2)$.
 * 
 */

#if 0
  if (t_m_ts > -6.0e-9 && TRUE)
  {
    gint i, j;

    //printf ("EXAM: % 22.15g % 22.15g % 22.15g % 22.15g\n", log(fabs(dlnmnu)), log(fabs(qbar)), dtanh_epsilon_dupsilon, lnch_lnmnu - 3.0 * lnch_epsilon);
    printf ("# (% 22.15g)\n", t_m_ts);
    printf ("#------------------------ANALITIC----------------------------------------\n");
    for (i = 0; i < 4; i++)
    {
      printf ("J:    ");
      for (j = 0; j < 4; j++)
      {
        printf ("% 22.15e ", SUN_DENSE_ACCESS (J, i, j));
      }
      printf ("\n");
    }

{
  guint i;
  const guint fparam_len = 4;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p       = NV_Ith_S (y, i);
    const gdouble p_scale = fabs (p);
    const gdouble h       = p_scale * GSL_ROOT3_DBL_EPSILON;
    const gdouble pph     = p + h;
    const gdouble pmh     = p - h;
    const gdouble twoh    = pph - pmh;
    const gdouble one_2h  = 1.0 / twoh;
    guint j;

    for (j = 0; j < 4; j++)
    {
      NV_Ith_S (tmp1, j) = NV_Ith_S (y, j);
    }

    NV_Ith_S (tmp1, i) = pph;
    _ncm_hoaa_dlnmnu_only_sing_f (t_m_ts, tmp1, tmp2, jac_data);
    
    NV_Ith_S (tmp1, i) = pmh;
    _ncm_hoaa_dlnmnu_only_sing_f (t_m_ts, tmp1, tmp3, jac_data);

    for (j = 0; j < 4; j++)
    {
      SUN_DENSE_ACCESS (J, j, i) = (NV_Ith_S (tmp2, j) - NV_Ith_S (tmp3, j)) * one_2h;
    }
  }
}
printf ("#------------------------FINITEDIFF----------------------------------------\n");
    for (i = 0; i < 4; i++)
    {
      printf ("J:    ");
      for (j = 0; j < 4; j++)
      {
        printf ("% 22.15e ", SUN_DENSE_ACCESS (J, i, j));
      }
      printf ("\n");
    }

    
  }
#endif
