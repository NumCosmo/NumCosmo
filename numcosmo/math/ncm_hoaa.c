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
 * H = \frac{p^2}{2m} + \frac{m(\nu^2 - V)q^2}{2},
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
 * p &= \sqrt{2Im\nu}\cos\theta.
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
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_rng.h"
#include "ncm_enum_types.h"

/*#undef HAVE_SUNDIALS_ARKODE*/

#ifndef HAVE_SUNDIALS_ARKODE
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#else
#include <arkode/arkode.h>
#include <arkode/arkode_dense.h>
#endif /* HAVE_SUNDIALS_ARKODE */

#include <nvector/nvector_serial.h>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sf_bessel.h>

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
  gdouble sin_shift;
  gdouble cos_shift;
  gint shift;
  glong shift_t;
  GArray *theta, *psi, *LnI, *LnJ, *tarray;
  NcmSpline *theta_s, *psi_s, *LnI_s, *LnJ_s;
  gdouble sigma0;
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

G_DEFINE_ABSTRACT_TYPE (NcmHOAA, ncm_hoaa, G_TYPE_OBJECT);

static gdouble
_ncm_hoaa_sin_thetabf (NcmHOAA *hoaa, const gdouble sin_thetab, const gdouble cos_thetab)
{
  return sin_thetab * hoaa->priv->cos_shift + cos_thetab * hoaa->priv->sin_shift;
}

static gdouble
_ncm_hoaa_cos_thetabf (NcmHOAA *hoaa, const gdouble sin_thetab, const gdouble cos_thetab)
{
  return cos_thetab * hoaa->priv->cos_shift - sin_thetab * hoaa->priv->sin_shift;
}

static void
_ncm_hoaa_update_shift (NcmHOAA *hoaa, const gint shift_d)
{
  gint shift_a = (hoaa->priv->shift + shift_d) % 4;

  hoaa->priv->shift = shift_a < 0 ? shift_a + 4 : shift_a;
 
  switch (hoaa->priv->shift)
  {
    case 0:
      hoaa->priv->cos_shift = +1.0;
      hoaa->priv->sin_shift = +0.0;
      break;
    case 1:
      hoaa->priv->cos_shift = +0.0;
      hoaa->priv->sin_shift = +1.0;
      break;
    case 2:
      hoaa->priv->cos_shift = -1.0;
      hoaa->priv->sin_shift = +0.0;
      break;
    case 3:
      hoaa->priv->cos_shift = +0.0;
      hoaa->priv->sin_shift = -1.0;
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static void
ncm_hoaa_init (NcmHOAA *hoaa)
{
  hoaa->k                      = 0.0;
  /* private properties */
  hoaa->priv                   = G_TYPE_INSTANCE_GET_PRIVATE (hoaa, NCM_TYPE_HOAA, NcmHOAAPrivate);
#ifndef HAVE_SUNDIALS_ARKODE
  hoaa->priv->cvode            = CVodeCreate (CV_BDF, CV_NEWTON); /*CVodeCreate (CV_BDF, CV_NEWTON);*/
  hoaa->priv->cvode_sing       = CVodeCreate (CV_BDF, CV_NEWTON);       /*CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);*/
  hoaa->priv->cvode_phase      = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL); /*CVodeCreate (CV_BDF, CV_NEWTON);*/
  hoaa->priv->cvode_init       = FALSE;
  hoaa->priv->cvode_sing_init  = FALSE;
#else
  hoaa->priv->arkode           = ARKodeCreate ();
  hoaa->priv->arkode_sing      = ARKodeCreate ();
  hoaa->priv->arkode_phase     = ARKodeCreate ();
  hoaa->priv->arkode_sing_init = FALSE;
#endif /* HAVE_SUNDIALS_ARKODE */
  hoaa->priv->opt       = NCM_HOAA_OPT_INVALID;
  hoaa->priv->reltol    = 0.0;
  hoaa->priv->abstol    = 0.0;
  hoaa->priv->ti        = 0.0;
  hoaa->priv->tf        = 0.0;
  hoaa->priv->t_cur     = 0.0;
  hoaa->priv->sin_shift = 0.0;
  hoaa->priv->cos_shift = 0.0;
  hoaa->priv->shift     = 0;
  hoaa->priv->shift_t   = 0;
  hoaa->priv->save_evol = FALSE;
  hoaa->priv->sigma     = N_VNew_Serial (1);
  hoaa->priv->y         = N_VNew_Serial (NCM_HOAA_VAR_SYS_SIZE);
  hoaa->priv->abstol_v  = N_VNew_Serial (NCM_HOAA_VAR_SYS_SIZE);
  hoaa->priv->ctrl      = ncm_model_ctrl_new (NULL);
  hoaa->priv->T         = gsl_root_fsolver_brent;
  hoaa->priv->s         = gsl_root_fsolver_alloc (hoaa->priv->T);
  hoaa->priv->rng       = ncm_rng_new (NULL);

  hoaa->priv->theta     = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->psi       = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->LnI       = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->LnJ       = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->tarray    = g_array_new (TRUE, TRUE, sizeof (gdouble));

  hoaa->priv->theta_s   = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->psi_s     = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->LnI_s     = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->LnJ_s     = ncm_spline_cubic_notaknot_new ();

  hoaa->priv->sigma0 = 0.0;
  
  ncm_rng_set_random_seed (hoaa->priv->rng, TRUE);
}

static void
_ncm_hoaa_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmHOAA *hoaa = NCM_HOAA (object);
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
      hoaa->priv->opt = g_value_get_enum (value);
      if (hoaa->priv->opt == NCM_HOAA_OPT_INVALID)
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
  g_return_if_fail (NCM_IS_HOAA (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, hoaa->priv->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, hoaa->priv->abstol);
      break;
    case PROP_K:
      g_value_set_double (value, hoaa->k);
      break;
    case PROP_TI:
      g_value_set_double (value, hoaa->priv->ti);
      break;
    case PROP_TF:
      g_value_set_double (value, hoaa->priv->tf);
      break;
    case PROP_SAVE_EVOL:
      g_value_set_boolean (value, hoaa->priv->save_evol);
      break;
    case PROP_OPT:
      g_value_set_enum (value, hoaa->priv->opt);
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

  ncm_model_ctrl_clear (&hoaa->priv->ctrl);
  ncm_rng_clear (&hoaa->priv->rng);

  ncm_spline_clear (&hoaa->priv->theta_s);
  ncm_spline_clear (&hoaa->priv->psi_s);
  ncm_spline_clear (&hoaa->priv->LnI_s);
  ncm_spline_clear (&hoaa->priv->LnJ_s);
  
  g_clear_pointer (&hoaa->priv->theta,  (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->psi,    (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->LnI,    (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->LnJ,    (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->tarray, (GDestroyNotify) g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_hoaa_parent_class)->dispose (object);
}

static void
_ncm_hoaa_finalize (GObject *object)
{
  NcmHOAA *hoaa = NCM_HOAA (object);

#ifndef HAVE_SUNDIALS_ARKODE
  if (hoaa->priv->cvode != NULL)
  {
    CVodeFree (&hoaa->priv->cvode);
    hoaa->priv->cvode = NULL;
  }
  if (hoaa->priv->cvode_sing != NULL)
  {
    CVodeFree (&hoaa->priv->cvode_sing);
    hoaa->priv->cvode_sing = NULL;
  }
  if (hoaa->priv->cvode_phase != NULL)
  {
    CVodeFree (&hoaa->priv->cvode_phase);
    hoaa->priv->cvode_phase = NULL;
  }
#else
  if (hoaa->priv->arkode != NULL)
  {
    ARKodeFree (&hoaa->priv->arkode);
    hoaa->priv->arkode = NULL;
  }
  if (hoaa->priv->arkode_sing != NULL)
  {
    ARKodeFree (&hoaa->priv->arkode_sing);
    hoaa->priv->arkode_sing = NULL;
  }
  if (hoaa->priv->arkode_phase != NULL)
  {
    ARKodeFree (&hoaa->priv->arkode_phase);
    hoaa->priv->arkode_phase = NULL;
  }
#endif /* HAVE_SUNDIALS_ARKODE */  
  
  if (hoaa->priv->sigma != NULL)
  {
    N_VDestroy (hoaa->priv->sigma);
    hoaa->priv->sigma = NULL;
  }

  if (hoaa->priv->y != NULL)
  {
    N_VDestroy (hoaa->priv->y);
    hoaa->priv->y = NULL;
  }

  if (hoaa->priv->abstol_v != NULL)
  {
    N_VDestroy (hoaa->priv->abstol_v);
    hoaa->priv->abstol_v = NULL;
  }

  g_clear_pointer (&hoaa->priv->s, (GDestroyNotify) gsl_root_fsolver_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_hoaa_parent_class)->finalize (object);
}

static void
ncm_hoaa_class_init (NcmHOAAClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcmHOAAPrivate));

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
  
  klass->eval_mnu         = NULL;
  klass->eval_nu          = NULL;
  klass->eval_V           = NULL;
  klass->eval_system      = NULL;

  klass->nsing            = NULL;
  klass->get_sing_info    = NULL;
  klass->eval_sing_mnu    = NULL;
  klass->eval_sing_dlnmnu = NULL;
  klass->eval_sing_V      = NULL;

  klass->prepare          = NULL;
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
  hoaa->priv->reltol = reltol;
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
  hoaa->priv->abstol = abstol;
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
  hoaa->priv->ti = ti;
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
  hoaa->priv->tf = tf;
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
  if (hoaa->priv->save_evol != save_evol)
  {
    ncm_model_ctrl_force_update (hoaa->priv->ctrl);    
    hoaa->priv->save_evol = save_evol;
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

static gint 
_ncm_hoaa_phase_root (realtype t, N_Vector y, realtype *gout, gpointer user_data)
{
  gout[0] = gsl_pow_2 (NV_Ith_S (y, 0) / (NCM_HOAA_MAX_ANGLE * M_PI + 0.25 * M_PI)) - 1.0;
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
  const gdouble psi        = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  const gdouble sin_theta  = sin (theta);
  const gdouble sin_psi    = sin (psi);
  const gdouble sin2_theta = sin_theta * sin_theta;
  const gdouble sin2_psi   = sin_psi * sin_psi;

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);

  NV_Ith_S (ydot, NCM_HOAA_VAR_THETAB) = nu - sin2_theta * Vnu + 0.5 * dlnmnu * sin_2theta;
  NV_Ith_S (ydot, NCM_HOAA_VAR_EPSILON)   = nu - sin2_psi * Vnu   + 0.5 * dlnmnu * sin_2psi;

  NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = Vnu * sin_2theta - dlnmnu * cos_2theta;
  NV_Ith_S (ydot, NCM_HOAA_VAR_LNJ)   = Vnu * sin_2psi   - dlnmnu * cos_2psi;
  
  return 0;
}

static gint
_ncm_hoaa_full_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble theta = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi   = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi, &sin_2psi, &cos_2psi);

  DENSE_ELEM (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_THETAB) = dlnmnu * cos_2theta - Vnu * sin_2theta;
  DENSE_ELEM (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_GAMMA)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB) = 2.0 * (dlnmnu * sin_2theta + Vnu * cos_2theta);
  DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_GAMMA)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON,   NCM_HOAA_VAR_EPSILON)   = dlnmnu * cos_2psi - Vnu * sin_2psi;
  DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON,   NCM_HOAA_VAR_LNJ)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_EPSILON)   = 2.0 * (dlnmnu * sin_2psi + Vnu * cos_2psi);
  DENSE_ELEM (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_LNJ)   = 0.0;

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
  const gdouble psi        = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  const gdouble sin_theta  = sin (theta);
  const gdouble sin_psi    = sin (psi);
  const gdouble sin2_theta = sin_theta * sin_theta;
  const gdouble sin2_psi   = sin_psi * sin_psi;
  const gdouble sin_2theta = sin (2.0 * theta);
  const gdouble sin_2psi   = sin (2.0 * psi);
  
  gdouble nu, dlnmnu, Vnu;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  NV_Ith_S (ydot, NCM_HOAA_VAR_THETAB) = nu - sin2_theta * Vnu;
  NV_Ith_S (ydot, NCM_HOAA_VAR_EPSILON)   = nu - sin2_psi * Vnu;

  NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = Vnu * sin_2theta;
  NV_Ith_S (ydot, NCM_HOAA_VAR_LNJ)   = Vnu * sin_2psi;
  
  return 0;
}

static gint
_ncm_hoaa_V_only_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble theta = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi   = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);

  DENSE_ELEM (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_THETAB) = - Vnu * sin_2theta;
  DENSE_ELEM (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_GAMMA)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB) = 2.0 * Vnu * cos_2theta;
  DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_GAMMA)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON,   NCM_HOAA_VAR_EPSILON)   = - Vnu * sin_2psi;
  DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON,   NCM_HOAA_VAR_LNJ)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_EPSILON)   = 2.0 * Vnu * cos_2psi;
  DENSE_ELEM (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_LNJ)   = 0.0;

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
#endif
static gint
_ncm_hoaa_dlnmnu_only_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;

  const gdouble thetab  = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble epsilon = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  gdouble nu, dlnmnu, Vnu, sin_thetab, cos_thetab;
  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);
  
  sincos (thetab, &sin_thetab, &cos_thetab);

  {
    const gdouble sin_thetabf  = _ncm_hoaa_sin_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble cos_thetabf  = _ncm_hoaa_cos_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble tanh_epsilon = tanh (epsilon);
    const gdouble cosh_epsilon = cosh (epsilon);

    NV_Ith_S (ydot, NCM_HOAA_VAR_THETAB)  = + 2.0 * nu + dlnmnu * tanh_epsilon * sin_thetabf;
    NV_Ith_S (ydot, NCM_HOAA_VAR_EPSILON) = - dlnmnu * cos_thetabf;
    NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)   = - dlnmnu * sin_thetabf / cosh_epsilon;
  }

  return 0;
}

static gint
_ncm_hoaa_dlnmnu_only_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg          = (NcmHOAAArg *) jac_data;

  const gdouble thetab  = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble epsilon = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  gdouble nu, dlnmnu, Vnu, sin_thetab, cos_thetab;
  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);
  
  sincos (thetab, &sin_thetab, &cos_thetab);

  {
    const gdouble sin_thetabf  = _ncm_hoaa_sin_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble cos_thetabf  = _ncm_hoaa_cos_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble tanh_epsilon = tanh (epsilon);
    const gdouble cosh_epsilon = cosh (epsilon);

    /* 2.0 * nu + dlnmnu * tanh_epsilon * sin_thetabf */
    DENSE_ELEM (J, NCM_HOAA_VAR_THETAB,  NCM_HOAA_VAR_THETAB)  = + dlnmnu * tanh_epsilon * cos_thetabf;
    DENSE_ELEM (J, NCM_HOAA_VAR_THETAB,  NCM_HOAA_VAR_EPSILON) = + dlnmnu * sin_thetabf / gsl_pow_2 (cosh_epsilon);

    /* - dlnmnu * cos_thetabf */
    DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON, NCM_HOAA_VAR_THETAB)  = + dlnmnu * sin_thetabf;

    /* - dlnmnu * sin_thetabf / cosh_epsilon */
    DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB)  = - dlnmnu * cos_thetabf / cosh_epsilon;
    DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_EPSILON) = + dlnmnu * tanh_epsilon * sin_thetabf / cosh_epsilon;
  }

  return 0;
}

static gint
_ncm_hoaa_dlnmnu_only_sing_f (realtype t_m_ts, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;
  
  gdouble nu, dlnmnu, Vnu, sin_thetab, cos_thetab;
  ncm_hoaa_eval_sing_system (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing, &nu, &dlnmnu, &Vnu);

  const gdouble mnu     = ncm_hoaa_eval_sing_mnu (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing);
  const gdouble lnmnu   = log (mnu);
  const gdouble thetab  = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  gdouble epsilon;

  switch (arg->st)
  {
    case NCM_HOAA_SING_TYPE_ZERO:
      epsilon = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON) - lnmnu;
      break;
    case NCM_HOAA_SING_TYPE_INF:
      epsilon = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON) + lnmnu;
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  
  sincos (thetab, &sin_thetab, &cos_thetab);

  {
    const gdouble sin_thetabf  = _ncm_hoaa_sin_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble cos_thetabf  = _ncm_hoaa_cos_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble tanh_epsilon = tanh (epsilon);
    const gdouble cosh_epsilon = cosh (epsilon);

    NV_Ith_S (ydot, NCM_HOAA_VAR_THETAB) = 2.0 * nu + dlnmnu * tanh_epsilon * sin_thetabf;
    NV_Ith_S (ydot, NCM_HOAA_VAR_GAMMA)  = - dlnmnu * sin_thetabf / cosh_epsilon;

    switch (arg->st)
    {
      case NCM_HOAA_SING_TYPE_ZERO:
        printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %d\n", 
                t_m_ts, 
                thetab, 1.0 - cos_thetabf, ncm_util_1mcosx (sin_thetabf, cos_thetabf), 
                dlnmnu, lnmnu, NV_Ith_S (y, NCM_HOAA_VAR_EPSILON),
                dlnmnu * ncm_util_1mcosx (sin_thetabf, cos_thetabf),
                2.0 * nu + dlnmnu * tanh_epsilon * sin_thetabf,
                arg->hoaa->priv->shift);
        NV_Ith_S (ydot, NCM_HOAA_VAR_EPSILON) = + dlnmnu * ncm_util_1mcosx (sin_thetabf, cos_thetabf);
        break;
      case NCM_HOAA_SING_TYPE_INF:
        NV_Ith_S (ydot, NCM_HOAA_VAR_EPSILON) = - dlnmnu * ncm_util_1pcosx (sin_thetabf, cos_thetabf);
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

  return 0;
}

static gint
_ncm_hoaa_dlnmnu_only_sing_J (_NCM_SUNDIALS_INT_TYPE N, realtype t_m_ts, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg   = (NcmHOAAArg *) jac_data;

  gdouble nu, dlnmnu, Vnu, sin_thetab, cos_thetab;
  ncm_hoaa_eval_sing_system (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing, &nu, &dlnmnu, &Vnu);

  const gdouble mnu     = ncm_hoaa_eval_sing_mnu (arg->hoaa, arg->model, t_m_ts, arg->hoaa->k, arg->sing);
  const gdouble lnmnu   = log (mnu);
  const gdouble thetab  = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  gdouble epsilon;

  switch (arg->st)
  {
    case NCM_HOAA_SING_TYPE_ZERO:
      epsilon = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON) - lnmnu;
      break;
    case NCM_HOAA_SING_TYPE_INF:
      epsilon = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON) + lnmnu;
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  
  sincos (thetab, &sin_thetab, &cos_thetab);

  {
    const gdouble sin_thetabf  = _ncm_hoaa_sin_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble cos_thetabf  = _ncm_hoaa_cos_thetabf (arg->hoaa, sin_thetab, cos_thetab);
    const gdouble tanh_epsilon = tanh (epsilon);
    const gdouble cosh_epsilon = cosh (epsilon);

    /* 2.0 * nu + dlnmnu * tanh_epsilon * sin_thetabf */
    DENSE_ELEM (J, NCM_HOAA_VAR_THETAB,  NCM_HOAA_VAR_THETAB)  = + dlnmnu * tanh_epsilon * cos_thetabf;
    DENSE_ELEM (J, NCM_HOAA_VAR_THETAB,  NCM_HOAA_VAR_EPSILON) = + dlnmnu * sin_thetabf / gsl_pow_2 (cosh_epsilon);

    /* - dlnmnu * (cos_thetabf +/- 1) */
    DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON, NCM_HOAA_VAR_THETAB)  = + dlnmnu * sin_thetabf;

    /* - dlnmnu * sin_thetabf / cosh_epsilon */
    DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB)  = - dlnmnu * cos_thetabf / cosh_epsilon;
    DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_EPSILON) = + dlnmnu * tanh_epsilon * sin_thetabf / cosh_epsilon;
  }

  return 0;
}

#if 0
static gint
_ncm_hoaa_dlnmnu_only_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble theta = NV_Ith_S (y, NCM_HOAA_VAR_THETAB);
  const gdouble psi   = NV_Ith_S (y, NCM_HOAA_VAR_EPSILON);

  gdouble nu, dlnmnu, Vnu, sin_2theta, cos_2theta, sin_2psi, cos_2psi;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  sincos (2.0 * theta, &sin_2theta, &cos_2theta);
  sincos (2.0 * psi,   &sin_2psi,   &cos_2psi);

  DENSE_ELEM (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_THETAB) = dlnmnu * cos_2theta;
  DENSE_ELEM (J, NCM_HOAA_VAR_THETAB, NCM_HOAA_VAR_GAMMA)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_THETAB) = 2.0 * dlnmnu * sin_2theta;
  DENSE_ELEM (J, NCM_HOAA_VAR_GAMMA,   NCM_HOAA_VAR_GAMMA)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON,   NCM_HOAA_VAR_EPSILON)   = dlnmnu * cos_2psi;
  DENSE_ELEM (J, NCM_HOAA_VAR_EPSILON,   NCM_HOAA_VAR_LNJ)   = 0.0;

  DENSE_ELEM (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_EPSILON)   = 2.0 * dlnmnu * sin_2psi;
  DENSE_ELEM (J, NCM_HOAA_VAR_LNJ,   NCM_HOAA_VAR_LNJ)   = 0.0;

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
  const gdouble Rsigma_t0 = 0.125 * M_PI; /* pi / 8 */
  gdouble Athetab, Aepsilon, Agamma;

  hoaa->priv->tc      = t0;
  hoaa->priv->sigma_c = Rsigma_t0;

  NV_Ith_S (hoaa->priv->sigma, 0) = Rsigma_t0;

  /* t0 must match hoaa->priv->tc at this point! */
  ncm_hoaa_eval_adiabatic_approx (hoaa, model, t0, &Athetab, &Aepsilon, &Agamma);

  NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB)  = Athetab;
  NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) = Aepsilon;
  NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_GAMMA)   = Agamma;

  printf ("# COND INI: % 21.15g % 21.15g % 21.15g % 21.15g\n", t0, 
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB), 
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON), 
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_GAMMA));
}

void 
_ncm_hoaa_prepare_integrator (NcmHOAA *hoaa, NcmModel *model, const gdouble t0)
{
  gint flag;
#ifndef HAVE_SUNDIALS_ARKODE
  CVRhsFn f;
  CVDlsDenseJacFn J;  
#else
  ARKRhsFn f;
  ARKDlsDenseJacFn J;    
#endif /* HAVE_SUNDIALS_ARKODE */
  
  switch (hoaa->priv->opt)
  {
    case NCM_HOAA_OPT_FULL:
//      f = _ncm_hoaa_full_f;
//      J = _ncm_hoaa_full_J;
      break;
    case NCM_HOAA_OPT_V_ONLY:
//      f = _ncm_hoaa_V_only_f;
//      J = _ncm_hoaa_V_only_J;
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

  NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_VAR_THETAB)  = 0.0;
  NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_VAR_EPSILON) = 0.0;
  NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_VAR_GAMMA)   = 0.0;

#ifndef HAVE_SUNDIALS_ARKODE
  if (!hoaa->priv->cvode_init)
  {
    flag = CVodeInit (hoaa->priv->cvode, f, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    hoaa->priv->t_cur = t0;

    flag = CVodeSVtolerances (hoaa->priv->cvode, hoaa->priv->reltol, hoaa->priv->abstol_v);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (hoaa->priv->cvode, 0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVDense (hoaa->priv->cvode, NCM_HOAA_VAR_SYS_SIZE);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVDlsSetDenseJacFn (hoaa->priv->cvode, J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );

    flag = CVodeSetInitStep (hoaa->priv->cvode, fabs (t0) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = CVodeInit (hoaa->priv->cvode_phase, _ncm_hoaa_phase_f, t0, hoaa->priv->sigma);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSStolerances (hoaa->priv->cvode_phase, hoaa->priv->reltol, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (hoaa->priv->cvode_phase, 0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVodeRootInit (hoaa->priv->cvode_phase, 1, &_ncm_hoaa_phase_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );
    
    hoaa->priv->cvode_init = TRUE;
  }
  else
  {    
    flag = CVodeReInit (hoaa->priv->cvode, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    hoaa->priv->t_cur = t0;

    flag = CVodeSetInitStep (hoaa->priv->cvode, fabs (t0) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = CVodeReInit (hoaa->priv->cvode_phase, t0, hoaa->priv->sigma);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }  
#else
  if (!hoaa->priv->arkode_init)
  {
    flag = ARKodeInit (hoaa->priv->arkode, f, NULL, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    hoaa->priv->t_cur = t0;

    flag = ARKodeSVtolerances (hoaa->priv->arkode, hoaa->priv->reltol, hoaa->priv->abstol_v);
    NCM_CVODE_CHECK (&flag, "ARKodeSVtolerances", 1, );

    flag = ARKodeSetMaxNumSteps (hoaa->priv->arkode, 0);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    flag = ARKDense (hoaa->priv->arkode, NCM_HOAA_VAR_SYS_SIZE);
    NCM_CVODE_CHECK (&flag, "ARKDense", 1, );

    flag = ARKDlsSetDenseJacFn (hoaa->priv->arkode, J);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetDenseJacFn", 1, );
    
    flag = ARKodeSetLinear (hoaa->priv->arkode, 1);
    NCM_CVODE_CHECK (&flag, "ARKodeSetLinear", 1, );
    
    flag = ARKodeSetOrder (hoaa->priv->arkode, 8);
    NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );

    //flag = ARKodeSetERKTableNum (hoaa->priv->arkode, FEHLBERG_13_7_8);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetERKTableNum", 1, );
    
    flag = ARKodeSetInitStep (hoaa->priv->arkode, fabs (t0) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = ARKodeInit (hoaa->priv->arkode_phase, _ncm_hoaa_phase_f,  NULL, t0, hoaa->priv->sigma);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSStolerances (hoaa->priv->arkode_phase, hoaa->priv->reltol, 0.0);
    NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );
    
    flag = ARKodeSetMaxNumSteps (hoaa->priv->arkode_phase, 0);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    flag = ARKodeSetOrder (hoaa->priv->arkode_phase, 7);
    NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );

    flag = ARKodeRootInit (hoaa->priv->arkode_phase, 1, &_ncm_hoaa_phase_root);
    NCM_CVODE_CHECK (&flag, "ARKodeRootInit", 1, );
    
    hoaa->priv->arkode_init = TRUE;
  }
  else
  {
    flag = ARKodeReInit (hoaa->priv->arkode, f, NULL, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    hoaa->priv->t_cur = t0;

    flag = ARKodeSetInitStep (hoaa->priv->arkode, fabs (t0) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

    /* PHASE integrator */
    flag = ARKodeReInit (hoaa->priv->arkode_phase, _ncm_hoaa_phase_f,  NULL, t0, hoaa->priv->sigma);
    NCM_CVODE_CHECK (&flag, "ARKodeReInit", 1, );
  }
#endif /* HAVE_SUNDIALS_ARKODE */
}

void 
_ncm_hoaa_prepare_integrator_sing (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const guint sing, const gdouble ts, const NcmHOAASingType st)
{
  gint flag;
#ifndef HAVE_SUNDIALS_ARKODE
  CVRhsFn f;
  CVDlsDenseJacFn J;  
#else
  ARKRhsFn f;
  ARKDlsDenseJacFn J;    
#endif /* HAVE_SUNDIALS_ARKODE */
  
  g_assert_cmpfloat (t_m_ts + ts, ==, hoaa->priv->t_cur);
  g_assert_cmpuint (sing, <, ncm_hoaa_nsing (hoaa, model, hoaa->k));

  switch (hoaa->priv->opt)
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

  {
    const gdouble mnu   = ncm_hoaa_eval_sing_mnu (hoaa, model, t_m_ts, hoaa->k, sing);
    const gdouble lnmnu = log (mnu);

    switch (st)
    {
      case NCM_HOAA_SING_TYPE_ZERO:
        NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) += lnmnu;
        break;
      case NCM_HOAA_SING_TYPE_INF:
        NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) -= lnmnu;
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

#ifndef HAVE_SUNDIALS_ARKODE
  if (!hoaa->priv->cvode_sing_init)
  {
    flag = CVodeInit (hoaa->priv->cvode_sing, f, t_m_ts, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSVtolerances (hoaa->priv->cvode_sing, hoaa->priv->reltol, hoaa->priv->abstol_v);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (hoaa->priv->cvode_sing, 0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVDense (hoaa->priv->cvode_sing, NCM_HOAA_VAR_SYS_SIZE);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVDlsSetDenseJacFn (hoaa->priv->cvode_sing, J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );

    flag = CVodeSetInitStep (hoaa->priv->cvode_sing, fabs (t_m_ts) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    flag = CVodeSetMaxErrTestFails (hoaa->priv->cvode_sing, 100);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxErrTestFails", 1, );

    hoaa->priv->cvode_sing_init = TRUE;
  }
  else
  {    
    flag = CVodeReInit (hoaa->priv->cvode_sing, t_m_ts, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetInitStep (hoaa->priv->cvode_sing, fabs (t_m_ts) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    flag = CVDlsSetDenseJacFn (hoaa->priv->cvode_sing, J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );
  }
#else
  if (!hoaa->priv->arkode_sing_init)
  {
    flag = ARKodeInit (hoaa->priv->arkode_sing, NULL, f, t_m_ts, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSVtolerances (hoaa->priv->arkode_sing, hoaa->priv->reltol, hoaa->priv->abstol_v);
    NCM_CVODE_CHECK (&flag, "ARKodeSVtolerances", 1, );

    flag = ARKodeSetMaxNumSteps (hoaa->priv->arkode_sing, 0);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    flag = ARKDense (hoaa->priv->arkode_sing, NCM_HOAA_VAR_SYS_SIZE);
    NCM_CVODE_CHECK (&flag, "ARKDense", 1, );

    //flag = ARKDlsSetDenseJacFn (hoaa->priv->arkode_sing, J);
    //NCM_CVODE_CHECK (&flag, "ARKDlsSetDenseJacFn", 1, );
    
    //flag = ARKodeSetOrder (hoaa->priv->arkode_sing, 5);
    //NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );

    flag = ARKodeSetIRKTableNum (hoaa->priv->arkode_sing, ARK548L2SA_DIRK_8_4_5);
    NCM_CVODE_CHECK (&flag, "ARKodeSetIRKTableNum", 1, );

    flag = ARKodeSetAdaptivityMethod (hoaa->priv->arkode_sing, 2, 1, 0, NULL);
    NCM_CVODE_CHECK (&flag, "ARKodeSetAdaptivityMethod", 1, );
    
    flag = ARKodeSetInitStep (hoaa->priv->arkode_sing, fabs (t_m_ts) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

    hoaa->priv->arkode_sing_init = TRUE;
  }
  else
  {
    flag = ARKodeReInit (hoaa->priv->arkode_sing, NULL, f, t_m_ts, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );
  }
#endif /* HAVE_SUNDIALS_ARKODE */
}

static gdouble 
_ncm_hoaa_test_epsilon_dlnmnu (gdouble at, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, test;
  const gdouble t = sinh (at);

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  test = dlnmnu / nu;

  return log (fabs (test / arg->prec));
}

static gdouble 
_ncm_hoaa_test_epsilon_Vnu (gdouble at, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, test;
  const gdouble t = sinh (at);

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu);

  test = Vnu / nu;
  
  return log (fabs (test / arg->prec));
}

static gdouble
_ncm_hoaa_search_initial_time_by_func (NcmHOAA *hoaa, NcmModel *model, gdouble epsilon, gdouble (*test_epsilon) (gdouble, gpointer))
{
  gdouble at_hi           = asinh (hoaa->priv->tf);
  gdouble at_lo           = asinh (hoaa->priv->ti);
  gint iter               = 0;
  gint max_iter           = 100000;
  const gdouble pass_step = (at_hi - at_lo) * 1.0e-4;
  NcmHOAAArg arg          = {hoaa, model, epsilon, -1, NCM_HOAA_SING_TYPE_INVALID};
  
  gsl_function F;
  gint status;
  gdouble at0, test_ep;
  
  F.function = test_epsilon;
  F.params   = &arg;

  if ((test_ep = test_epsilon (at_lo, F.params)) > 0.0)
  {
    g_warning ("_ncm_hoaa_search_initial_time_by_func: system is not adiabatic at the initial time setting initial conditions at t = % 21.15g, test / epsilon = % 21.15g",
               sinh (at_lo), test_ep);
    return sinh (at_lo);
  }

  at_hi = at_lo + pass_step;
  while (((test_ep = test_epsilon (at_hi, F.params)) < 0.0) && (iter < max_iter))
  {
    at_lo  = at_hi;
    at_hi += pass_step;
    iter++;
  }

  if (iter >= max_iter)
  {
    g_warning ("_ncm_hoaa_search_initial_time_by_func: cannot find non-adiabtic regime, setting t0 to ti.");
    return hoaa->priv->ti;
  }

  iter = 0;
  gsl_root_fsolver_set (hoaa->priv->s, &F, at_lo, at_hi);
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (hoaa->priv->s);

    at0    = gsl_root_fsolver_root (hoaa->priv->s);
    at_lo  = gsl_root_fsolver_x_lower (hoaa->priv->s);
    at_hi  = gsl_root_fsolver_x_upper (hoaa->priv->s);
    status = gsl_root_test_interval (at_lo, at_hi, 0.0, hoaa->priv->reltol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if ((iter >= max_iter) || (status != GSL_SUCCESS))
  {
    g_warning ("_ncm_hoaa_search_initial_time_by_func: cannot intial time with the required precision.");
  }
  
  return sinh (at0);  
}

static gdouble
_ncm_hoaa_search_initial_time (NcmHOAA *hoaa, NcmModel *model, gdouble epsilon)
{
  switch (hoaa->priv->opt)
  {
    case NCM_HOAA_OPT_FULL:
    {
      const gdouble t0_dlnmnu = _ncm_hoaa_search_initial_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_dlnmnu);
      const gdouble t0_Vnu    = _ncm_hoaa_search_initial_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_Vnu);
      return GSL_MIN (t0_dlnmnu, t0_Vnu);
      break;
    }
    case NCM_HOAA_OPT_V_ONLY:
    {
      return _ncm_hoaa_search_initial_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_Vnu);
      break;
    }      
    case NCM_HOAA_OPT_DLNMNU_ONLY:
    {
      return _ncm_hoaa_search_initial_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_dlnmnu);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
}

static gdouble
_ncm_hoaa_search_final_time_by_func (NcmHOAA *hoaa, NcmModel *model, gdouble epsilon, gdouble (*test_epsilon) (gdouble, gpointer))
{
  gdouble at_hi           = asinh (hoaa->priv->tf);
  gdouble at_lo           = asinh (hoaa->priv->ti);
  gint iter               = 0;
  gint max_iter           = 100000;
  const gdouble pass_step = (at_hi - at_lo) * 1.0e-4;
  NcmHOAAArg arg          = {hoaa, model, epsilon, -1, NCM_HOAA_SING_TYPE_INVALID};

  
  gsl_function F;
  gint status;
  gdouble at0, test_ep;
  
  F.function = test_epsilon;
  F.params   = &arg;

  if ((test_ep = test_epsilon (at_hi, F.params)) > 0.0)
    return sinh (at_hi);

  at_lo = at_hi - pass_step;
  while (((test_ep = test_epsilon (at_lo, F.params)) < 0.0) && (iter < max_iter))
  {
    iter++;
    at_hi  = at_lo;
    at_lo -= pass_step;
  }

  if (iter >= max_iter)
  {
    g_warning ("_ncm_hoaa_search_final_time_by_func: cannot find non-adiabtic regime, setting t1 to tf.");
    return hoaa->priv->tf;
  }

  iter = 0;
  gsl_root_fsolver_set (hoaa->priv->s, &F, at_lo, at_hi);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (hoaa->priv->s);

    at0    = gsl_root_fsolver_root (hoaa->priv->s);
    at_lo  = gsl_root_fsolver_x_lower (hoaa->priv->s);
    at_hi  = gsl_root_fsolver_x_upper (hoaa->priv->s);
    status = gsl_root_test_interval (at_lo, at_hi, 0.0, hoaa->priv->reltol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if ((iter >= max_iter) || (status != GSL_SUCCESS))
  {
    g_warning ("_ncm_hoaa_search_final_time: cannot intial time with the required precision.");
  }
  
  return sinh (at0);  
}

static gdouble
_ncm_hoaa_search_final_time (NcmHOAA *hoaa, NcmModel *model, gdouble epsilon)
{
  switch (hoaa->priv->opt)
  {
    case NCM_HOAA_OPT_FULL:
    {
      const gdouble t1_dlnmnu = _ncm_hoaa_search_final_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_dlnmnu);
      const gdouble t1_Vnu    = _ncm_hoaa_search_final_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_Vnu);

      return GSL_MAX (t1_dlnmnu, t1_Vnu);
      break;
    }
    case NCM_HOAA_OPT_V_ONLY:
    {
      return _ncm_hoaa_search_final_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_Vnu);
      break;
    }      
    case NCM_HOAA_OPT_DLNMNU_ONLY:
    {
      return _ncm_hoaa_search_final_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_dlnmnu);
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
  gint flag;
  gdouble delta_thetab       = 0.0;
  gdouble delta_thetab_i     = 0.0;
  guint nrestarts            = 0;
  gdouble last_t             = -1.0e100;
  NcmHOAAArg arg             = {hoaa, model, 0.0, -1, NCM_HOAA_SING_TYPE_INVALID};
  const gdouble sin_table[8] = {0.0, 1.0 / M_SQRT2, 1.0, +1.0 / M_SQRT2, +0.0, -1.0 / M_SQRT2, -1.0, -1.0 / M_SQRT2};
  const gdouble cos_table[8] = {1.0, 1.0 / M_SQRT2, 0.0, -1.0 / M_SQRT2, -1.0, -1.0 / M_SQRT2, +0.0, +1.0 / M_SQRT2};

  if (hoaa->priv->sigma0 == 0.0)
    hoaa->priv->sigma0 = 0.125 * M_PI - hoaa->k * hoaa->priv->t_cur;

  g_assert_cmpfloat (tf, >, hoaa->priv->t_cur);

#ifndef HAVE_SUNDIALS_ARKODE
  flag = CVodeSetStopTime (hoaa->priv->cvode, tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetUserData (hoaa->priv->cvode, &arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
#else
  flag = ARKodeSetStopTime (hoaa->priv->arkode, tf);
  NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

  flag = ARKodeSetUserData (hoaa->priv->arkode, &arg);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */

  while (TRUE)
  {
    gint restart = 0;
#ifndef HAVE_SUNDIALS_ARKODE
    flag = CVode (hoaa->priv->cvode, tf, hoaa->priv->y, &hoaa->priv->t_cur, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_hoaa_evol_save]", 1, );
#else
    flag = ARKode (hoaa->priv->arkode, tf, hoaa->priv->y, &hoaa->priv->t_cur, ARK_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "ARKode[_ncm_hoaa_evol_save]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
    {
      const gdouble thetab  = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB);
      const gdouble epsilon = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON);
      const gdouble gamma   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_GAMMA);
      
      if (hoaa->priv->t_cur > last_t + fabs (last_t) * 0.001)
      {
        const gdouble mnu          = ncm_hoaa_eval_mnu (hoaa, model, hoaa->priv->t_cur, hoaa->k);
        const gdouble lnmnu        = log (mnu);
        const glong n              = hoaa->priv->shift_t % 8;
        const gdouble cos_thetab_2 = cos (0.5 * thetab) * cos_table[n] - sin (0.5 * thetab) * sin_table[n];
        const gdouble sin_thetab_2 = sin (0.5 * thetab) * cos_table[n] + cos (0.5 * thetab) * sin_table[n];
        const gdouble q            = - exp (0.5 * (+ gamma - epsilon - lnmnu)) * cos_thetab_2 + exp (0.5 * (+ gamma + epsilon - lnmnu)) * sin_thetab_2;
        const gdouble v            = + exp (0.5 * (- gamma - epsilon - lnmnu)) * cos_thetab_2 + exp (0.5 * (- gamma + epsilon - lnmnu)) * sin_thetab_2;
        const gdouble t            = hoaa->priv->t_cur;
        const gdouble arg          = hoaa->k * t + hoaa->priv->sigma0;
        const gdouble hypot1       = hypot (q, v) * sqrt (mnu / 2.0);
        const gdouble Aq           = - M_SQRT2 * sin (arg) / (sqrt (hoaa->k) * t);
        const gdouble Av           = - M_SQRT2 * cos (arg) / (sqrt (hoaa->k) * t);
        const gdouble hypot2       = hypot (Aq, Av) * sqrt (mnu / 2.0);

        printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %8ld %d\n",
                hoaa->priv->t_cur,
                thetab / M_PI,
                epsilon,
                gamma,
                q,
                v,
                hypot1,
                Aq,
                Av,
                hypot2,
                hypot1 / hypot2 - 1.0,
                hoaa->priv->shift_t % 8,
                hoaa->priv->shift
                );
/*
        gdouble Athetab, Aepsilon, Agamma;
        ncm_hoaa_eval_adiabatic_approx (hoaa, model, hoaa->priv->t_cur, &Athetab, &Aepsilon, &Agamma);
        const gdouble thetab_f  = gsl_sf_angle_restrict_symm (thetab + hoaa->priv->shift * 0.5 * M_PI);
        const gdouble Athetab_f = gsl_sf_angle_restrict_symm (Athetab);
        
        printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %8ld %d\n",
                hoaa->priv->t_cur,
                Athetab_f / thetab_f - 1.0,
                Aepsilon / epsilon - 1.0,
                Agamma / gamma - 1.0,
                thetab_f,
                epsilon,
                gamma,                
                Athetab_f,
                Aepsilon,
                Agamma,                
                q,
                v,
                0.0,
                hoaa->priv->shift_t % 8,
                hoaa->priv->shift
                );
*/
        last_t = hoaa->priv->t_cur;
      }          

      if (fabs (thetab) > 0.99 * 0.5 * M_PI)
      {
        const gdouble thetab_s = gsl_sf_angle_restrict_symm (thetab);
        const gint shift_d = thetab_s >= 0.0 ? 1 : -1;

        if (fabs (thetab_s) > 0.99 * 0.5 * M_PI)
        {
          _ncm_hoaa_update_shift (hoaa, shift_d);
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB) = thetab_s - shift_d * 0.5 * M_PI;
        }
        else
        {
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB) = thetab_s;
        }

        delta_thetab_i       = (thetab - NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB));
        delta_thetab        += delta_thetab_i;
        hoaa->priv->shift_t += round (delta_thetab_i / (0.5 * M_PI));

        /*printf ("SHIFT %ld %d %ld | %d\n", hoaa->priv->shift_t, hoaa->priv->shift, hoaa->priv->shift_t % 4, shift_d);*/

        restart = 1;
      }

      if (hoaa->priv->t_cur == tf)
        break;

      if (restart)
      {
#ifndef HAVE_SUNDIALS_ARKODE
        flag = CVodeReInit (hoaa->priv->cvode, hoaa->priv->t_cur, hoaa->priv->y);
        NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

        flag = CVodeSetInitStep (hoaa->priv->cvode, fabs (hoaa->priv->t_cur) * hoaa->priv->reltol);
        NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
#else
        flag = ARKodeResize (hoaa->priv->arkode, hoaa->priv->y, 1.0, hoaa->priv->t_cur, _ncm_hoaa_arkode_resize, NULL);
        NCM_CVODE_CHECK (&flag, "ARKodeResize", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
        nrestarts++;
      }
    }
  }
}

void
_ncm_hoaa_evol_sing_save (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const guint sing, const gdouble ts, const NcmHOAASingType st)
{
  gint flag;
  gdouble delta_thetab       = 0.0;
  gdouble delta_thetab_i     = 0.0;
  guint nrestarts            = 0;
  gdouble last_t             = -1.0e100;
  NcmHOAAArg arg             = {hoaa, model, 0.0, sing, st};
  const gdouble sin_table[8] = {0.0, 1.0 / M_SQRT2, 1.0, +1.0 / M_SQRT2, +0.0, -1.0 / M_SQRT2, -1.0, -1.0 / M_SQRT2};
  const gdouble cos_table[8] = {1.0, 1.0 / M_SQRT2, 0.0, -1.0 / M_SQRT2, -1.0, -1.0 / M_SQRT2, +0.0, +1.0 / M_SQRT2};
  gdouble tstep_m_ts;

  g_assert_cmpfloat (t_m_ts + ts, >, hoaa->priv->t_cur);

#ifndef HAVE_SUNDIALS_ARKODE
  flag = CVodeSetStopTime (hoaa->priv->cvode_sing, t_m_ts);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );    

  flag = CVodeSetUserData (hoaa->priv->cvode_sing, &arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
#else
  flag = ARKodeSetStopTime (hoaa->priv->arkode_sing, t_m_ts);
  NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

  flag = ARKodeSetUserData (hoaa->priv->arkode_sing, &arg);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */

  while (TRUE)
  {
    gint restart = 0;
#ifndef HAVE_SUNDIALS_ARKODE
    flag = CVode (hoaa->priv->cvode_sing, t_m_ts, hoaa->priv->y, &tstep_m_ts, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_hoaa_evol_sing_save]", 1, );
#else
    flag = ARKode (hoaa->priv->arkode_sing, t_m_ts, hoaa->priv->y, &tstep_m_ts, ARK_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "ARKode[_ncm_hoaa_evol_sing_save]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */

    hoaa->priv->t_cur = ts + tstep_m_ts;

    {
      const gdouble mnu    = ncm_hoaa_eval_sing_mnu (hoaa, model, tstep_m_ts, hoaa->k, sing);
      const gdouble lnmnu  = log (mnu);
      const gdouble thetab = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB);
      const gdouble gamma  = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_GAMMA);
      gdouble epsilon;

      switch (st)
      {
        case NCM_HOAA_SING_TYPE_ZERO:
          epsilon = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) - lnmnu;
          break;
        case NCM_HOAA_SING_TYPE_INF:
          epsilon = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) + lnmnu;
          break;
        default:
          g_assert_not_reached ();
          break;
      }
      
      if (tstep_m_ts > last_t + fabs (last_t) * 0.01)
      {
        const glong n              = hoaa->priv->shift_t % 8;
        const gdouble cos_thetab_2 = cos (0.5 * thetab) * cos_table[n] - sin (0.5 * thetab) * sin_table[n];
        const gdouble sin_thetab_2 = sin (0.5 * thetab) * cos_table[n] + cos (0.5 * thetab) * sin_table[n];
        
        switch (st)
        {
          case NCM_HOAA_SING_TYPE_ZERO:
          {
            const gdouble q = - exp (0.5 * (+ gamma - epsilon - lnmnu)) * cos_thetab_2 + exp (0.5 * (+ gamma + epsilon - lnmnu)) * sin_thetab_2;
            const gdouble v = + exp (0.5 * (- gamma - epsilon - lnmnu)) * cos_thetab_2 + exp (0.5 * (- gamma + epsilon - lnmnu)) * sin_thetab_2;
           
/*            printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %8ld %d\n",
                    hoaa->priv->t_cur,
                    thetab / M_PI,
                    epsilon - lnmnu,
                    gamma,
                    q,
                    v,
                    hypot (q, v),
                    tstep_m_ts,
                    hoaa->priv->shift_t % 8,
                    hoaa->priv->shift
                    );
*/
        const gdouble t            = hoaa->priv->t_cur;
        const gdouble arg          = hoaa->k * t + hoaa->priv->sigma0;
        const gdouble hypot1       = hypot (q, v) * sqrt (mnu / 2.0);
        const gdouble Aq           = - M_SQRT2 * sin (arg) / (sqrt (hoaa->k) * t);
        const gdouble Av           = - M_SQRT2 * cos (arg) / (sqrt (hoaa->k) * t);
        const gdouble hypot2       = hypot (Aq, Av) * sqrt (mnu / 2.0);

        printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %8ld %d\n",
                hoaa->priv->t_cur,
                thetab / M_PI,
                epsilon,
                gamma,
                q,
                v,
                hypot1,
                Aq,
                Av,
                hypot2,
                hypot1 / hypot2 - 1.0,
                hoaa->priv->shift_t % 8,
                hoaa->priv->shift
                );


          break;
           }
          case NCM_HOAA_SING_TYPE_INF:
          {
            const gdouble q = - exp (0.5 * (+ gamma - epsilon - lnmnu)) * cos_thetab_2 + exp (0.5 * (+ gamma + epsilon - lnmnu)) * sin_thetab_2;
            const gdouble v = + exp (0.5 * (- gamma - epsilon - lnmnu)) * cos_thetab_2 + exp (0.5 * (- gamma + epsilon - lnmnu)) * sin_thetab_2;

            printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %8ld %d\n",
                    hoaa->priv->t_cur,
                    thetab / M_PI,
                    epsilon + lnmnu,
                    gamma,
                    q,
                    v,
                    hypot (q, v),
                    tstep_m_ts,
                    hoaa->priv->shift_t % 8,
                    hoaa->priv->shift
                    );
            break;
          }
          default:
            g_assert_not_reached ();
            break;
        }
        last_t = tstep_m_ts;
      }          

      if (fabs (thetab) > 0.99 * 0.5 * M_PI)
      {
        const gdouble thetab_s = gsl_sf_angle_restrict_symm (thetab);
        const gint shift_d = thetab_s >= 0.0 ? 1 : -1;

        if (fabs (thetab_s) > 0.99 * 0.5 * M_PI)
        {
          _ncm_hoaa_update_shift (hoaa, shift_d);
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB) = thetab_s - shift_d * 0.5 * M_PI;
        }
        else
        {
          NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB) = thetab_s;
        }

        delta_thetab_i       = (thetab - NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB));
        delta_thetab        += delta_thetab_i;
        hoaa->priv->shift_t += round (delta_thetab_i / (0.5 * M_PI));

        /*printf ("SHIFT %ld %d %ld | %d\n", hoaa->priv->shift_t, hoaa->priv->shift, hoaa->priv->shift_t % 4, shift_d);*/
        
        restart = 1;
      }
      if (tstep_m_ts == t_m_ts)
        break;

      if (restart)
      {
#ifndef HAVE_SUNDIALS_ARKODE
        flag = CVodeReInit (hoaa->priv->cvode_sing, tstep_m_ts, hoaa->priv->y);
        NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

        flag = CVodeSetInitStep (hoaa->priv->cvode_sing, fabs (tstep_m_ts) * hoaa->priv->reltol);
        NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
#else
        flag = ARKodeResize (hoaa->priv->arkode_sing, hoaa->priv->y, 1.0, tstep_m_ts, _ncm_hoaa_arkode_resize, NULL);
        NCM_CVODE_CHECK (&flag, "ARKodeResize", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
        nrestarts++;
      }
    }
  }

  {
    const gdouble mnu = ncm_hoaa_eval_sing_mnu (hoaa, model, t_m_ts, hoaa->k, sing);
    
    switch (st)
    {
      case NCM_HOAA_SING_TYPE_ZERO:
        NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) -= log (mnu);
        break;
      case NCM_HOAA_SING_TYPE_INF:
        NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) += log (mnu);
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }
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
  if (NCM_HOAA_GET_CLASS (hoaa)->prepare != NULL)
  {
    NCM_HOAA_GET_CLASS (hoaa)->prepare (hoaa, model);
  }

  g_assert_cmpfloat (hoaa->priv->tf, >, hoaa->priv->ti);

  {
    const gdouble epsilon = GSL_MAX (hoaa->priv->reltol, 1.0e-4);
    const gdouble t_ad_0  = _ncm_hoaa_search_initial_time (hoaa, model, epsilon);
    const gdouble t_ad_1  = _ncm_hoaa_search_final_time (hoaa, model, epsilon);
    const gdouble t_na_0  = _ncm_hoaa_search_initial_time (hoaa, model, 1.0);
    const gdouble t_na_1  = _ncm_hoaa_search_final_time (hoaa, model, 1.0);
    
    hoaa->priv->t_ad_0 = t_ad_0;
    hoaa->priv->t_ad_1 = t_ad_1;
    hoaa->priv->t_na_0 = t_na_0;
    hoaa->priv->t_na_1 = t_na_1;

    _ncm_hoaa_set_init_cond (hoaa, model, t_ad_0);
    hoaa->priv->cos_shift = 1.0;
    hoaa->priv->sin_shift = 0.0;
    hoaa->priv->shift     = 0;
    hoaa->priv->shift_t   = 0;
    
    _ncm_hoaa_prepare_integrator (hoaa, model, t_ad_0);

    if (hoaa->priv->save_evol)
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
          printf ("# SING % 21.15g % 21.15g % 21.15g\n", ts_a[i], dts_i, dts_f);
        }

        gsl_sort_index (sing_a, ts_a, 1, nsing);
        
        for (i = 0; i < nsing; i++)
        {
          gdouble ts, dts_i, dts_f;
          NcmHOAASingType st;

          ncm_hoaa_get_sing_info (hoaa, model, hoaa->k, sing_a[i], &ts, &dts_i, &dts_f, &st);

          if (ts < hoaa->priv->t_cur)
            continue;

          _ncm_hoaa_evol_save (hoaa, model, dts_i + ts);

          _ncm_hoaa_prepare_integrator_sing (hoaa, model, dts_i, sing_a[i], ts, st);
          _ncm_hoaa_evol_sing_save (hoaa, model, -1.0e-45, sing_a[i], ts, st);
          _ncm_hoaa_evol_sing_save (hoaa, model, +1.0e-45, sing_a[i], ts, st);
          _ncm_hoaa_evol_sing_save (hoaa, model, dts_f, sing_a[i], ts, st);

          _ncm_hoaa_prepare_integrator (hoaa, model, dts_f + ts);
        }
        _ncm_hoaa_evol_save (hoaa, model, t_ad_1);

        g_free (ts_a);
        g_free (sing_a);
      }
    }
    else
    {
#if 0      
      gdouble t_step;
      guint nrestarts = 0;
      while (TRUE)
      {
        gint restart = 0;
#ifndef HAVE_SUNDIALS_ARKODE
        gint flag = CVode (hoaa->priv->cvode, t_ad_1, hoaa->priv->y, &t_step, CV_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "CVode[ncm_hoaa_prepare]", 1, );
#else
        gint flag = ARKode (hoaa->priv->arkode, t_ad_1, hoaa->priv->y, &t_step, ARK_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "ARKode[ncm_hoaa_prepare]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
        {
          const gdouble theta   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB);
          const gdouble psi   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON);

          if (theta >= NCM_HOAA_MAX_ANGLE * M_PI)
          {
            NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_THETAB) = fmod (theta, M_PI);
            restart = 1;
          }

          if (psi >= NCM_HOAA_MAX_ANGLE * M_PI)
          {
            NV_Ith_S (hoaa->priv->y, NCM_HOAA_VAR_EPSILON) = fmod (psi, M_PI);
            restart = 1;
          }
          
          if (t_step == t_ad_1)
            break;
                  
          if (restart)
          {
#ifndef HAVE_SUNDIALS_ARKODE
            flag = CVodeReInit (hoaa->priv->cvode, t_step, hoaa->priv->y);
            NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
#else
            flag = ARKodeResize (hoaa->priv->arkode, hoaa->priv->y, 1.0, t_step, _ncm_hoaa_arkode_resize, NULL);
            NCM_CVODE_CHECK (&flag, "ARKodeResize", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
            nrestarts++;
          }
        }
      }
#endif
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
  t0[0] = hoaa->priv->t_ad_0;
  t1[0] = hoaa->priv->t_ad_1;
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
  gsl_root_fsolver_set (hoaa->priv->s, &gslF, dt0, dt1);
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (hoaa->priv->s);

    dt     = gsl_root_fsolver_root (hoaa->priv->s);
    dt0    = gsl_root_fsolver_x_lower (hoaa->priv->s);
    dt1    = gsl_root_fsolver_x_upper (hoaa->priv->s);
    status = gsl_root_test_interval (dt0, dt1, 0.0, hoaa->priv->reltol);
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
 * @thetab: (out): $\theta$
 * @epsilon: (out): $\psi$
 * @gamma: (out): $\ln(I)$
 * 
 * Calculates the adiabatic approximation at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_adiabatic_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *thetab, gdouble *epsilon, gdouble *gamma)
{
  if (ncm_cmp (t, hoaa->priv->tc, hoaa->priv->reltol) != 0)
  {
    gint flag;
    gdouble t_step;
    guint restarts = 0;
    NcmHOAAArg arg = {hoaa, model, 0.0, -1, NCM_HOAA_SING_TYPE_INVALID};

#ifndef HAVE_SUNDIALS_ARKODE
    flag = CVodeSetStopTime (hoaa->priv->cvode_phase, t);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );    

    flag = CVodeSetUserData (hoaa->priv->cvode_phase, &arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    while (TRUE) 
    {
      flag = CVode (hoaa->priv->cvode_phase, t, hoaa->priv->sigma, &t_step, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode[phase]", 1, );

      if (flag == CV_ROOT_RETURN)
      {
        NV_Ith_S (hoaa->priv->sigma, 0) = gsl_sf_angle_restrict_pos (NV_Ith_S (hoaa->priv->sigma, 0));

        flag = CVodeReInit (hoaa->priv->cvode_phase, t_step, hoaa->priv->sigma);
        NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
        restarts++;
      }

      if (t == t_step)
        break;
    }
#else
    flag = ARKodeSetStopTime (hoaa->priv->arkode_phase, t);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );    

    flag = ARKodeSetUserData (hoaa->priv->arkode_phase, &arg);
    NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

    while (TRUE) 
    {
      flag = ARKode (hoaa->priv->arkode_phase, t, hoaa->priv->sigma, &t_step, ARK_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "ARKode[phase]", 1, );

      if (flag == ARK_ROOT_RETURN)
      {
        NV_Ith_S (hoaa->priv->sigma, 0) = gsl_sf_angle_restrict_symm (NV_Ith_S (hoaa->priv->sigma, 0));

        flag = ARKodeResize (hoaa->priv->arkode_phase, hoaa->priv->sigma, 1.0e-5, t_step, _ncm_hoaa_arkode_resize, NULL);
        NCM_CVODE_CHECK (&flag, "ARKodeResize", 1, );
        restarts++;
      }

      if (t == t_step)
        break;
    }
#endif /* HAVE_SUNDIALS_ARKODE */
    hoaa->priv->sigma_c = gsl_sf_angle_restrict_symm (NV_Ith_S (hoaa->priv->sigma, 0));
  }
  
  {
    const gdouble Rsigma_t0 = hoaa->priv->sigma_c;
    const gdouble Isigma_t0 = hoaa->priv->sigma_c + 0.5 * M_PI;
    gdouble nu, dlnmnu, Vnu;
    gsl_function F_F, F_G;
    NcmHOAAArg arg;

    ncm_hoaa_eval_system (hoaa, model, t, hoaa->k, &nu, &dlnmnu, &Vnu);

    arg.model = model;
    arg.hoaa  = hoaa;
    arg.prec  = 0.0;

    F_F.params   = &arg;
    F_G.params   = &arg;

    F_F.function = &_ncm_hoaa_F;
    F_G.function = &_ncm_hoaa_G;
    
    {
      gdouble err;
      const gdouble F  = 0.5 * Vnu / nu;
      const gdouble G  = 0.5 * dlnmnu / nu;
      const gdouble dF = ncm_numdiff_1 (&F_F, t, fabs (t) * 0.1, &err);
      const gdouble dG = ncm_numdiff_1 (&F_G, t, fabs (t) * 0.1, &err); 
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

      thetab[0]  = Rsigma_t0 + Isigma_t0 + dtheta + dpsi;
      epsilon[0] = asinh (tan (dtheta - dpsi));
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
  gdouble nu, dlnmnu, Vnu;
  gsl_function F_F, F_G;
  NcmHOAAArg arg;

  ncm_hoaa_eval_system (hoaa, model, t, hoaa->k, &nu, &dlnmnu, &Vnu);

  arg.model = model;
  arg.hoaa  = hoaa;
  arg.prec  = 0.0;

  F_F.params   = &arg;
  F_G.params   = &arg;

  F_F.function = &_ncm_hoaa_F;
  F_G.function = &_ncm_hoaa_G;

  {
    gdouble err;
    const gdouble F  = 0.5 * Vnu / nu;
    const gdouble G  = 0.5 * dlnmnu / nu;
    const gdouble dF = ncm_numdiff_1 (&F_F, t, fabs (t) * 0.1, &err);
    const gdouble dG = ncm_numdiff_1 (&F_G, t, fabs (t) * 0.1, &err); 
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
 * @theta: (out): $\theta$
 * @psi: (out): $\psi$
 * @LnI: (out): $\ln(I)$
 * @LnJ: (out): $\ln(J)$
 * 
 * Calculates the AA variables at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *theta, gdouble *psi, gdouble *LnI, gdouble *LnJ)
{
  if (!hoaa->priv->save_evol)
  {
    g_error ("ncm_hoaa_eval_adiabatic_approx: save_evol must be enables.");
  }
  else if ((t < hoaa->priv->t_ad_0) || (t > hoaa->priv->t_ad_1))
  {
    gdouble Athetab, Aepsilon, Agamma;
    
    ncm_hoaa_eval_adiabatic_approx (hoaa, model, t, &Athetab, &Aepsilon, &Agamma);
    g_assert_not_reached ();
  }
  else
  {
    theta[0] = ncm_spline_eval (hoaa->priv->theta_s, t);
    psi[0]   = ncm_spline_eval (hoaa->priv->psi_s, t);
    LnI[0]   = ncm_spline_eval (hoaa->priv->LnI_s, t);
    LnJ[0]   = ncm_spline_eval (hoaa->priv->LnJ_s, t);
  }
}

/**
 * ncm_hoaa_eval_AA2CV:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @theta: $\theta$
 * @psi: $\psi$
 * @LnI: $\ln(I)$
 * @LnJ: $\ln(J)$
 * @Rphi: (out): $\Re(\phi)$
 * @Iphi: (out): $\Im(\phi)$
 * @RPphi: (out): $\Re(P_\phi)$
 * @IPphi: (out): $\Im(P_\phi)$
 * 
 * Change the variables from AA to complex.
 * 
 */
void 
ncm_hoaa_eval_AA2CV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble theta, const gdouble psi, const gdouble LnI, const gdouble LnJ, gdouble *Rphi, gdouble *Iphi, gdouble *RPphi, gdouble *IPphi)
{
  const gdouble mnu  = ncm_hoaa_eval_mnu (hoaa, model, t, hoaa->k);
  const gdouble smnu = sqrt (mnu);

  gdouble sin_theta, cos_theta, sin_psi, cos_psi;

  sincos (theta, &sin_theta, &cos_theta);
  sincos (psi,   &sin_psi,   &cos_psi);

  {
    const gdouble q  = M_SQRT2 / smnu * exp (0.5 * LnI) * sin_theta;
    const gdouble Pq = M_SQRT2 * smnu * exp (0.5 * LnI) * cos_theta;

    const gdouble v  = M_SQRT2 / smnu * exp (0.5 * LnJ) * sin_psi;
    const gdouble Pv = M_SQRT2 * smnu * exp (0.5 * LnJ) * cos_psi;

    const complex double phi  = (q  + I * v ) / (2.0 * I);
    const complex double Pphi = (Pq + I * Pv) / (2.0 * I);

    Rphi[0]  = creal (phi);
    Iphi[0]  = cimag (phi);
    RPphi[0] = creal (Pphi);
    IPphi[0] = cimag (Pphi);
  }
}

/**
 * ncm_hoaa_eval_CV2AA:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @Rphi: $\Re(\phi)$
 * @Iphi: $\Im(\phi)$
 * @RPphi: $\Re(P_\phi)$
 * @IPphi: $\Im(P_\phi)$
 * @theta: (out): $\theta$
 * @psi: (out): $\psi$
 * @LnI: (out): $\ln(I)$
 * @LnJ: (out): $\ln(J)$
 * 
 * Change the variables from complex to AA.
 * 
 */
void 
ncm_hoaa_eval_CV2AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble Rphi, const gdouble Iphi, const gdouble RPphi, const gdouble IPphi, gdouble *theta, gdouble *psi, gdouble *LnI, gdouble *LnJ)
{
  const gdouble mnu  = ncm_hoaa_eval_mnu (hoaa, model, t, hoaa->k);

  const complex double twoI_phi  = 2.0 * I * (Rphi  + I * Iphi);
  const complex double twoI_Pphi = 2.0 * I * (RPphi + I * IPphi);

  const gdouble q  = creal (twoI_phi);
  const gdouble Pq = creal (twoI_Pphi);;

  const gdouble v  = cimag (twoI_phi);
  const gdouble Pv = cimag (twoI_Pphi);

  theta[0] = atan2 (q, Pq);
  psi[0]   = atan2 (v, Pv);
  LnI[0]   = log (0.5 * (mnu * q * q + Pq * Pq / mnu));
  LnJ[0]   = log (0.5 * (mnu * v * v + Pv * Pv / mnu));
}

/**
 * ncm_hoaa_eval_CV:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @Rphi: (out): $\Re(\phi)$
 * @Iphi: (out): $\Im(\phi)$
 * @RPphi: (out): $\Re(P_\phi)$
 * @IPphi: (out): $\Im(P_\phi)$
 * 
 * Calculates the complex variables at $t$ and $k$.
 * 
 */
void 
ncm_hoaa_eval_CV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *Rphi, gdouble *Iphi, gdouble *RPphi, gdouble *IPphi)
{
  gdouble theta, psi, LnI, LnJ;
  ncm_hoaa_eval_AA (hoaa, model, t, &theta, &psi, &LnI, &LnJ);
  ncm_hoaa_eval_AA2CV (hoaa, model, t, theta, psi, LnI, LnJ, Rphi, Iphi, RPphi, IPphi);
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
  const gdouble twopi2  = 2.0 * gsl_pow_2 (M_PI);
  const gdouble k3_2pi2 = gsl_pow_3 (hoaa->k) / twopi2;
  gdouble Rphi, Iphi, RPphi, IPphi;
    
  ncm_hoaa_eval_CV (hoaa, model, t, &Rphi, &Iphi, &RPphi, &IPphi);

  Delta_phi[0]  = k3_2pi2 * (Rphi  * Rphi  + Iphi  * Iphi);
  Delta_Pphi[0] = k3_2pi2 * (RPphi * RPphi + IPphi * IPphi);
}

/**
 * ncm_hoaa_eval_mnu: (virtual eval_mnu)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the mass term at $t$ and $k$.
 *
 * Returns: the mass $m$.
 */
/**
 * ncm_hoaa_eval_nu: (virtual eval_nu)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the frequency term at $t$ and $k$.
 *
 * Returns: the frequency $\nu$.
 */
/**
 * ncm_hoaa_eval_V: (virtual eval_V)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the potential term at $t$ and $k$.
 *
 * Returns: the potential $V$.
 */
/**
 * ncm_hoaa_eval_dlnmnu: (virtual eval_dlnmnu)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Evaluates the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$ term at $t$ and $k$.
 *
 * Returns: the potential $\mathrm{d}(m\nu)/\mathrm{d}t$.
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
 * @model: a #NcmModel
 * @k: mode $k$
 *
 * Gets the number of singular points $m(t_s) = 0$ for the problem in hand.
 *
 * Returns: the number of singular points.
 */
/**
 * ncm_hoaa_get_sing_info: (virtual get_sing_info)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @k: mode $k$
 * @sing: singularity number
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
 * @model: a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 *
 * Evaluates the mass term at $t$ and $k$.
 *
 * Returns: the mass $m$.
 */
/**
 * ncm_hoaa_eval_sing_nu: (virtual eval_sing_nu)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 *
 * Evaluates the frequency term at $t$ and $k$.
 *
 * Returns: the frequency $\nu$.
 */
/**
 * ncm_hoaa_eval_sing_V: (virtual eval_sing_V)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
 *
 * Evaluates the potential term at $t$ and $k$.
 *
 * Returns: the potential $V$.
 */
/**
 * ncm_hoaa_eval_sing_dlnmnu: (virtual eval_sing_dlnmnu)
 * @hoaa: a #NcmHOAA
 * @model: a #NcmModel
 * @t_m_ts: time $t - t_s$
 * @k: mode $k$
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
 * @nu: (out): the frequency $\nu$
 * @dlnmnu: (out): the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$
 * @Vnu: (out): the potential over frequency term $V/\nu$
 *
 * Evaluates the system functions at $t$ and $k$.
 *
 */
