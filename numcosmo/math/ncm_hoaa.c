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
 * Each child should implement the functions ncm_hoaa_eval_m(), ncm_hoaa_eval_nu() and
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

#ifndef HAVE_SUNDIALS_ARKODE
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#else
#include <arkode/arkode.h>
#include <arkode/arkode_dense.h>
#endif /* HAVE_SUNDIALS_ARKODE */

#include <nvector/nvector_serial.h>

#include <gsl/gsl_roots.h>

typedef struct _NcmHOAAArg
{
  NcmModel *model;
  NcmHOAA *hoaa;
  gdouble prec;
} NcmHOAAArg;

struct _NcmHOAAPrivate
{
#ifndef HAVE_SUNDIALS_ARKODE
  gpointer cvode;
  gboolean cvode_init;
#else
  gpointer arkode;
  gboolean arkode_init;
#endif /* HAVE_SUNDIALS_ARKODE */
  N_Vector y;
  N_Vector abstol_v;
  gdouble reltol;
  gdouble abstol;
  gdouble ti;
  gdouble tf;
  gboolean save_evol;
  NcmModelCtrl *ctrl;
  NcmSpline *nu_s;
  NcmHOAAArg arg;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  NcmRNG *rng;
  gdouble sigma_c;
  gdouble tc;
  gdouble t0, t1;
  GArray *RQ, *IQ, *RLnI, *ILnI, *tarray;
  NcmSpline *RQ_s, *IQ_s, *RLnI_s, *ILnI_s;
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
};

G_DEFINE_ABSTRACT_TYPE (NcmHOAA, ncm_hoaa, G_TYPE_OBJECT);

static void
ncm_hoaa_init (NcmHOAA *hoaa)
{
  hoaa->k                = 0.0;
  hoaa->priv             = G_TYPE_INSTANCE_GET_PRIVATE (hoaa, NCM_TYPE_HOAA, NcmHOAAPrivate);
#ifndef HAVE_SUNDIALS_ARKODE
  hoaa->priv->cvode      = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL); /*CVodeCreate (CV_BDF, CV_NEWTON);*/
  hoaa->priv->cvode_init = FALSE;
#else
  hoaa->priv->arkode      = ARKodeCreate ();
  hoaa->priv->arkode_init = FALSE;
#endif /* HAVE_SUNDIALS_ARKODE */
  hoaa->priv->reltol     = 0.0;
  hoaa->priv->abstol     = 0.0;
  hoaa->priv->ti         = 0.0;
  hoaa->priv->tf         = 0.0;
  hoaa->priv->save_evol  = FALSE;
  hoaa->priv->y          = N_VNew_Serial (NCM_HOAA_SYS_SIZE);
  hoaa->priv->abstol_v   = N_VNew_Serial (NCM_HOAA_SYS_SIZE);
  hoaa->priv->ctrl       = ncm_model_ctrl_new (NULL);
  hoaa->priv->nu_s       = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->T          = gsl_root_fsolver_brent;
  hoaa->priv->s          = gsl_root_fsolver_alloc (hoaa->priv->T);
  hoaa->priv->rng        = ncm_rng_new (NULL);

  hoaa->priv->RQ         = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->IQ         = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->RLnI       = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->ILnI       = g_array_new (TRUE, TRUE, sizeof (gdouble));
  hoaa->priv->tarray     = g_array_new (TRUE, TRUE, sizeof (gdouble));

  hoaa->priv->RQ_s       = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->IQ_s       = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->RLnI_s     = ncm_spline_cubic_notaknot_new ();
  hoaa->priv->ILnI_s     = ncm_spline_cubic_notaknot_new ();
  
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

  ncm_spline_clear (&hoaa->priv->nu_s);

  ncm_spline_clear (&hoaa->priv->RQ_s);
  ncm_spline_clear (&hoaa->priv->IQ_s);
  ncm_spline_clear (&hoaa->priv->RLnI_s);
  ncm_spline_clear (&hoaa->priv->ILnI_s);
  
  g_clear_pointer (&hoaa->priv->RQ,     (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->IQ,     (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->RLnI,   (GDestroyNotify) g_array_unref);
  g_clear_pointer (&hoaa->priv->ILnI,   (GDestroyNotify) g_array_unref);
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
#else
  if (hoaa->priv->arkode != NULL)
  {
    ARKodeFree (&hoaa->priv->arkode);
    hoaa->priv->arkode = NULL;
  }
#endif /* HAVE_SUNDIALS_ARKODE */  
  
  if (hoaa->priv->y != NULL)
  {
    N_VDestroy (hoaa->priv->y);
    hoaa->priv->y = NULL;
  }

  if (hoaa->priv->abstol_v != NULL)
  {
    N_VDestroy (hoaa->priv->y);
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

  klass->eval_m      = NULL;
  klass->eval_nu     = NULL;
  klass->eval_V      = NULL;
  klass->eval_system = NULL;
  klass->prepare     = NULL;
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

static gdouble 
_ncm_hoaa_F (gdouble t, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, Nt;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

  /*printf ("F % 21.15e % 21.15g\n", t, 0.5 * Vnu / nu);*/

  return 0.5 * Vnu / nu;
}

static gdouble 
_ncm_hoaa_G (gdouble t, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, Nt;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

  /*printf ("G % 21.15e % 21.15g\n", t, 0.5 * dlnmnu / nu);*/
  
  return 0.5 * dlnmnu / nu;
}

static gdouble 
_ncm_hoaa_phase (gdouble t, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu = ncm_hoaa_eval_nu (arg->hoaa, arg->model, t, arg->hoaa->k);

  return nu;
}

void 
_ncm_hoaa_prepare_phase (NcmHOAA *hoaa, NcmModel *model)
{
  const gdouble prec = hoaa->priv->reltol;
  NcmHOAAArg arg;
  gsl_function F;

  arg.model = model;
  arg.hoaa  = hoaa;
  arg.prec  = 0.0;

  F.params   = &arg;

  F.function = &_ncm_hoaa_phase;
  ncm_spline_set_func (hoaa->priv->nu_s, NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, &F, hoaa->priv->ti, hoaa->priv->tf, 0, prec);
}

/******************************************************************************************************/
/*************************************** Time gauge, V and dlnmnu *************************************/
/***************************************          START           *************************************/
/******************************************************************************************************/
static gint
_ncm_hoaa_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmHOAAArg *arg = (NcmHOAAArg *) f_data;

  const gdouble RQ     = NV_Ith_S (y, NCM_HOAA_RQi);
  const gdouble IQ     = NV_Ith_S (y, NCM_HOAA_IQi);

  const gdouble sinRQ  = sin (RQ);
  const gdouble sinIQ  = sin (IQ);
  const gdouble sinRQ2 = sinRQ * sinRQ;
  const gdouble sinIQ2 = sinIQ * sinIQ;

  gdouble nu, dlnmnu, Vnu, Nt, sin2RQ, cos2RQ, sin2IQ, cos2IQ;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

  sincos (2.0 * RQ, &sin2RQ, &cos2RQ);
  sincos (2.0 * IQ, &sin2IQ, &cos2IQ);

  NV_Ith_S (ydot, NCM_HOAA_RQi)   = Nt * (nu - sinRQ2 * Vnu + 0.5 * dlnmnu * sin2RQ);
  NV_Ith_S (ydot, NCM_HOAA_IQi)   = Nt * (nu - sinIQ2 * Vnu + 0.5 * dlnmnu * sin2IQ);

  NV_Ith_S (ydot, NCM_HOAA_RLnIi) = Nt * (- (dlnmnu * cos2RQ - Vnu * sin2RQ));
  NV_Ith_S (ydot, NCM_HOAA_ILnIi) = Nt * (- (dlnmnu * cos2IQ - Vnu * sin2IQ));

/*printf ("% 21.15e % 21.15e % 21.15e % 21.15e\n", t, nu, sinRQ2 * Vnu / nu, 0.5 * dlnmnu * sin2RQ / nu);*/
  
#if 0
  printf ("Y  % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", 
          t,
          NV_Ith_S (y, NCM_HOAA_RQi),
          NV_Ith_S (y, NCM_HOAA_IQi),
          NV_Ith_S (y, NCM_HOAA_RLnIi),
          NV_Ith_S (y, NCM_HOAA_ILnIi)        
          );
  printf ("DY % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", 
          t,
          NV_Ith_S (ydot, NCM_HOAA_RQi),
          NV_Ith_S (ydot, NCM_HOAA_IQi),
          NV_Ith_S (ydot, NCM_HOAA_RLnIi),
          NV_Ith_S (ydot, NCM_HOAA_ILnIi)        
          );
#endif
  
  return 0;
}

static gint
_ncm_hoaa_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) jac_data;

  const gdouble RQ = NV_Ith_S (y, NCM_HOAA_RQi);
  const gdouble IQ = NV_Ith_S (y, NCM_HOAA_IQi);

  gdouble nu, dlnmnu, Vnu, Nt, sin2RQ, cos2RQ, sin2IQ, cos2IQ;

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

  sincos (2.0 * RQ, &sin2RQ, &cos2RQ);
  sincos (2.0 * IQ, &sin2IQ, &cos2IQ);

  DENSE_ELEM (J, NCM_HOAA_RQi,   NCM_HOAA_RQi)   = Nt * (dlnmnu * cos2RQ - Vnu * sin2RQ);
  DENSE_ELEM (J, NCM_HOAA_RQi,   NCM_HOAA_RLnIi) = 0.0;

  DENSE_ELEM (J, NCM_HOAA_RLnIi, NCM_HOAA_RQi)   = Nt * (2.0 * (dlnmnu * sin2RQ + Vnu * cos2RQ));
  DENSE_ELEM (J, NCM_HOAA_RLnIi, NCM_HOAA_RLnIi) = 0.0;

  DENSE_ELEM (J, NCM_HOAA_IQi,   NCM_HOAA_IQi)   = Nt * (dlnmnu * cos2IQ - Vnu * sin2IQ);
  DENSE_ELEM (J, NCM_HOAA_IQi,   NCM_HOAA_ILnIi) = 0.0;

  DENSE_ELEM (J, NCM_HOAA_ILnIi, NCM_HOAA_IQi)   = Nt * (2.0 * (dlnmnu * sin2IQ + Vnu * cos2IQ));
  DENSE_ELEM (J, NCM_HOAA_ILnIi, NCM_HOAA_ILnIi) = 0.0;

  return 0;
}
/******************************************************************************************************/
/*************************************** Time gauge, V and dlnmnu *************************************/
/***************************************           END            *************************************/
/******************************************************************************************************/

void 
_ncm_hoaa_prepare_integrator (NcmHOAA *hoaa, NcmModel *model, const gdouble t0, const gdouble t1)
{
  gint flag;

  hoaa->priv->arg.hoaa  = hoaa;
  hoaa->priv->arg.model = model;

  {
    const gdouble Rsigma_t0 = 0.125 * M_PI; /* pi / 8 */
    gdouble ARQ, AIQ, ARLnI, AILnI;
  
    hoaa->priv->tc      = t0;
    hoaa->priv->sigma_c = Rsigma_t0;

    ncm_hoaa_eval_adiabatic_approx (hoaa, model, t0, &ARQ, &AIQ, &ARLnI, &AILnI);

    NV_Ith_S (hoaa->priv->y, NCM_HOAA_RQi)   = ARQ;
    NV_Ith_S (hoaa->priv->y, NCM_HOAA_IQi)   = AIQ;
    NV_Ith_S (hoaa->priv->y, NCM_HOAA_RLnIi) = ARLnI;
    NV_Ith_S (hoaa->priv->y, NCM_HOAA_ILnIi) = AILnI;

    NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_RQi)   = 0.0;
    NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_IQi)   = 0.0;
    NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_RLnIi) = hoaa->priv->abstol;
    NV_Ith_S (hoaa->priv->abstol_v, NCM_HOAA_ILnIi) = hoaa->priv->abstol;
  }

#ifndef HAVE_SUNDIALS_ARKODE
  if (!hoaa->priv->cvode_init)
  {
    flag = CVodeInit (hoaa->priv->cvode, &_ncm_hoaa_f, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSVtolerances (hoaa->priv->cvode, hoaa->priv->reltol, hoaa->priv->abstol_v);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (hoaa->priv->cvode, -1);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVDense (hoaa->priv->cvode, NCM_HOAA_SYS_SIZE);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVDlsSetDenseJacFn (hoaa->priv->cvode, &_ncm_hoaa_J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );

    flag = CVodeSetUserData (hoaa->priv->cvode, &hoaa->priv->arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetInitStep (hoaa->priv->cvode, fabs (t0) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    flag = CVodeSetStopTime (hoaa->priv->cvode, t1);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );    
    
    hoaa->priv->cvode_init = TRUE;
  }
  else
  {    
    flag = CVodeReInit (hoaa->priv->cvode, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetInitStep (hoaa->priv->cvode, fabs (t0) * hoaa->priv->reltol);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    flag = CVodeSetStopTime (hoaa->priv->cvode, t1);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );
  }  
#else
  if (!hoaa->priv->arkode_init)
  {
    flag = ARKodeInit (hoaa->priv->arkode, &_ncm_hoaa_f,  NULL, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSVtolerances (hoaa->priv->arkode, hoaa->priv->reltol, hoaa->priv->abstol_v);
    NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );

    flag = ARKodeSetMaxNumSteps (hoaa->priv->arkode, -1);
    NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

    flag = ARKDense (hoaa->priv->arkode, NCM_HOAA_SYS_SIZE);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = ARKDlsSetDenseJacFn (hoaa->priv->arkode, &_ncm_hoaa_J);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetDenseJacFn", 1, );
    
    flag = ARKodeSetUserData (hoaa->priv->arkode, &hoaa->priv->arg);
    NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

    flag = ARKodeSetStopTime (hoaa->priv->arkode, t1);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );    

    flag = ARKodeSetOrder (hoaa->priv->arkode, 8);
    NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );  
    
    hoaa->priv->arkode_init = TRUE;
  }
  else
  {    
    flag = ARKodeReInit (hoaa->priv->arkode, &_ncm_hoaa_f,  NULL, t0, hoaa->priv->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

    flag = ARKodeSetStopTime (hoaa->priv->arkode, t1);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );
  }  
#endif /* HAVE_SUNDIALS_ARKODE */
}

static gdouble 
_ncm_hoaa_test_epsilon_dlnmnu (gdouble at, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, Nt, test;
  const gdouble t = sinh (at);

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

  test = dlnmnu / nu;

  return log (fabs (test / arg->prec));
}

static gdouble 
_ncm_hoaa_test_epsilon_Vnu (gdouble at, gpointer userdata)
{
  NcmHOAAArg *arg  = (NcmHOAAArg *) userdata;
  gdouble nu, dlnmnu, Vnu, Nt, test;
  const gdouble t = sinh (at);

  ncm_hoaa_eval_system (arg->hoaa, arg->model, t, arg->hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

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
  
  gsl_function F;
  gint status;
  gdouble at0, test_ep;
  
  F.function = test_epsilon;
  F.params   = &hoaa->priv->arg;

  hoaa->priv->arg.hoaa  = hoaa;
  hoaa->priv->arg.model = model;
  hoaa->priv->arg.prec  = epsilon;

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
    /*printf ("# SL % 21.15g % 21.15g % 21.15g % 21.15g\n", sinh (at_lo), sinh (at_hi), sinh (at0), test_epsilon (at0, F.params));*/
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
  const gdouble t0_dlnmnu = _ncm_hoaa_search_initial_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_dlnmnu);
  const gdouble t0_Vnu    = _ncm_hoaa_search_initial_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_Vnu);

  /*printf ("t0_dlnmnu = % 21.15g, t0_Vnu = % 21.15g\n", t0_dlnmnu, t0_Vnu);*/
  
  return GSL_MIN (t0_dlnmnu, t0_Vnu);
}

static gdouble
_ncm_hoaa_search_final_time_by_func (NcmHOAA *hoaa, NcmModel *model, gdouble epsilon, gdouble (*test_epsilon) (gdouble, gpointer))
{
  gdouble at_hi           = asinh (hoaa->priv->tf);
  gdouble at_lo           = asinh (hoaa->priv->ti);
  gint iter               = 0;
  gint max_iter           = 100000;
  const gdouble pass_step = (at_hi - at_lo) * 1.0e-4;
  
  gsl_function F;
  gint status;
  gdouble at0, test_ep;
  
  F.function = test_epsilon;
  F.params   = &hoaa->priv->arg;

  hoaa->priv->arg.hoaa  = hoaa;
  hoaa->priv->arg.model = model;
  hoaa->priv->arg.prec  = epsilon;
  
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
  /*printf ("% 21.15g % 21.15g % 21.15g % 21.15g\n", sinh (at_lo), sinh (at_hi), sinh (at0), _ncm_hoaa_test_epsilon (at0, F.params));*/
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (hoaa->priv->s);

    at0    = gsl_root_fsolver_root (hoaa->priv->s);
    at_lo  = gsl_root_fsolver_x_lower (hoaa->priv->s);
    at_hi  = gsl_root_fsolver_x_upper (hoaa->priv->s);
    status = gsl_root_test_interval (at_lo, at_hi, 0.0, hoaa->priv->reltol);

    /*printf ("% 21.15g % 21.15g % 21.15g % 21.15g\n", sinh (at_lo), sinh (at_hi), sinh (at0), _ncm_hoaa_test_epsilon (at0, F.params));*/
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
  const gdouble t1_dlnmnu = _ncm_hoaa_search_final_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_dlnmnu);
  const gdouble t1_Vnu    = _ncm_hoaa_search_final_time_by_func (hoaa, model, epsilon, &_ncm_hoaa_test_epsilon_Vnu);

  /*printf ("t1_dlnmnu = % 21.15g, t1_Vnu = % 21.15g\n", t1_dlnmnu, t1_Vnu);*/

  return GSL_MAX (t1_dlnmnu, t1_Vnu);
}

#ifdef HAVE_SUNDIALS_ARKODE
static gint 
_ncm_hoaa_arkode_resize (N_Vector y, N_Vector ytemplate, void *user_data)
{
  return 0;
}
#endif /* HAVE_SUNDIALS_ARKODE */

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
  
  _ncm_hoaa_prepare_phase (hoaa, model);

  {
    const gdouble epsilon = GSL_MAX (hoaa->priv->reltol, 1.0e-4);
    const gdouble t0      = _ncm_hoaa_search_initial_time (hoaa, model, epsilon);
    const gdouble t1      = _ncm_hoaa_search_final_time (hoaa, model, epsilon);

    /*printf ("t0 % 21.15e t1 % 21.15e\n", t0, t1);*/
    
    hoaa->priv->t0 = t0;
    hoaa->priv->t1 = t1;
    
    _ncm_hoaa_prepare_integrator (hoaa, model, t0, t1);

    if (hoaa->priv->save_evol)
    {
      gdouble t_step;
      gdouble delta_RQ = 0.0;
      gdouble delta_IQ = 0.0;
      guint nrestarts = 0;
      
      while (TRUE)
      {
        gint restart = 0;
#ifndef HAVE_SUNDIALS_ARKODE
        gint flag = CVode (hoaa->priv->cvode, t1, hoaa->priv->y, &t_step, CV_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "CVode[ncm_hoaa_prepare]", 1, );
#else
        gint flag = ARKode (hoaa->priv->arkode, t1, hoaa->priv->y, &t_step, ARK_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "ARKode[ncm_hoaa_prepare]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
        {
          const gdouble RQ   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_RQi);
          const gdouble IQ   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_IQi);
          const gdouble RLnI = NV_Ith_S (hoaa->priv->y, NCM_HOAA_RLnIi);
          const gdouble ILnI = NV_Ith_S (hoaa->priv->y, NCM_HOAA_ILnIi);
          const gdouble FRQ  = RQ + delta_RQ;
          const gdouble FIQ  = IQ + delta_IQ;          

          if (RQ >= NCM_HOAA_MAX_ANGLE * M_PI)
          {
            NV_Ith_S (hoaa->priv->y, NCM_HOAA_RQi) = fmod (RQ, M_PI);

            delta_RQ += (RQ - NV_Ith_S (hoaa->priv->y, NCM_HOAA_RQi));

            restart = 1;
          }

          if (IQ >= NCM_HOAA_MAX_ANGLE * M_PI)
          {
            NV_Ith_S (hoaa->priv->y, NCM_HOAA_IQi) = fmod (IQ, M_PI);

            delta_IQ += (IQ - NV_Ith_S (hoaa->priv->y, NCM_HOAA_IQi));
            
            restart = 1;
          }
          
          if (FALSE)
          {
            gdouble nu, dlnmnu, Vnu, Nt;
            gdouble ARQ, AIQ, ARLnI, AILnI;
            static gdouble merr_RQ = 0.0, merr_IQ = 0.0, merr_RLnI = 0.0, merr_ILnI = 0.0;
            gdouble err_RQ = 0.0, err_IQ = 0.0, err_RLnI = 0.0, err_ILnI = 0.0;
            
            ncm_hoaa_eval_adiabatic_approx (hoaa, model, t_step, &ARQ, &AIQ, &ARLnI, &AILnI);
            ncm_hoaa_eval_adiabatic_LnI_approx (hoaa, model, t_step, RQ, IQ, &ARLnI, &AILnI);
            ncm_hoaa_eval_system (hoaa, model, t_step, hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

            err_RQ   = fabs (fmod (RQ, M_PI) / fmod (ARQ, M_PI) - 1.0);
            err_IQ   = fabs (fmod (IQ, M_PI) / fmod (AIQ, M_PI) - 1.0);
            err_RLnI = fabs (RLnI / ARLnI - 1.0);
            err_ILnI = fabs (ILnI / AILnI - 1.0);
            
            merr_RQ   = GSL_MAX (merr_RQ,   err_RQ);
            merr_IQ   = GSL_MAX (merr_IQ,   err_IQ);
            merr_RLnI = GSL_MAX (merr_RLnI, err_RLnI);
            merr_ILnI = GSL_MAX (merr_ILnI, err_ILnI);
            
            printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                    t_step, 
                    err_RQ, 
                    err_IQ, 
                    err_RLnI, 
                    err_ILnI, 
                    merr_RQ, 
                    merr_IQ, 
                    merr_RLnI, 
                    merr_ILnI, 
                    fabs (dlnmnu) / nu, 
                    fabs (Vnu) / nu);
          }
          
          g_array_append_val (hoaa->priv->tarray, t_step);
          g_array_append_val (hoaa->priv->RQ,     FRQ);
          g_array_append_val (hoaa->priv->IQ,     FIQ);
          g_array_append_val (hoaa->priv->RLnI,   RLnI);
          g_array_append_val (hoaa->priv->ILnI,   ILnI);
          
          if (t_step == t1)
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
      
      /*printf ("# nrestarts: %u.\n", nrestarts);*/

      {
        NcmVector *t_v    = ncm_vector_new_array (hoaa->priv->tarray);
        NcmVector *RQ_v   = ncm_vector_new_array (hoaa->priv->RQ);
        NcmVector *IQ_v   = ncm_vector_new_array (hoaa->priv->IQ);
        NcmVector *RLnI_v = ncm_vector_new_array (hoaa->priv->RLnI);
        NcmVector *ILnI_v = ncm_vector_new_array (hoaa->priv->ILnI);

        ncm_spline_set (hoaa->priv->RQ_s,   t_v, RQ_v, TRUE);
        ncm_spline_set (hoaa->priv->IQ_s,   t_v, IQ_v, TRUE);
        ncm_spline_set (hoaa->priv->RLnI_s, t_v, RLnI_v, TRUE);
        ncm_spline_set (hoaa->priv->ILnI_s, t_v, ILnI_v, TRUE);

        ncm_vector_free (t_v);
        ncm_vector_free (RQ_v);
        ncm_vector_free (IQ_v);
        ncm_vector_free (RLnI_v);
        ncm_vector_free (ILnI_v);
      }
    }
    else
    {
      gdouble t_step;
      guint nrestarts = 0;
      
      while (TRUE)
      {
        gint restart = 0;
#ifndef HAVE_SUNDIALS_ARKODE
        gint flag = CVode (hoaa->priv->cvode, t1, hoaa->priv->y, &t_step, CV_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "CVode[ncm_hoaa_prepare]", 1, );
#else
        gint flag = ARKode (hoaa->priv->arkode, t1, hoaa->priv->y, &t_step, ARK_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "ARKode[ncm_hoaa_prepare]", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
        {
          const gdouble RQ   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_RQi);
          const gdouble IQ   = NV_Ith_S (hoaa->priv->y, NCM_HOAA_IQi);

          if (RQ >= NCM_HOAA_MAX_ANGLE * M_PI)
          {
            NV_Ith_S (hoaa->priv->y, NCM_HOAA_RQi) = fmod (RQ, M_PI);
            restart = 1;
          }

          if (IQ >= NCM_HOAA_MAX_ANGLE * M_PI)
          {
            NV_Ith_S (hoaa->priv->y, NCM_HOAA_IQi) = fmod (IQ, M_PI);
            restart = 1;
          }
          
          if (t_step == t1)
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
  t0[0] = hoaa->priv->t0;
  t1[0] = hoaa->priv->t1;
}

/**
 * ncm_hoaa_eval_adiabatic_approx:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @RQ: (out): $\theta_R$
 * @IQ: (out): $\theta_I$
 * @RLnI: (out): $I_R$
 * @ILnI: (out): $I_I$
 * 
 * Calculates the adiabatic approximation at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_adiabatic_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *RQ, gdouble *IQ, gdouble *RLnI, gdouble *ILnI)
{
  if (t != hoaa->priv->tc)
  {
    gdouble t1;
    do
    {
      gdouble dsigma, sigma;
      t1 = t;
      /*printf ("% 21.15g % 21.15g % 21.15g\n", t, t1, hoaa->priv->tc);*/

      while (fabs (dsigma = ncm_spline_eval_integ (hoaa->priv->nu_s, hoaa->priv->tc, t1)) > NCM_HOAA_MAX_ANGLE * M_PI)
      {
        t1 = 0.5 * (hoaa->priv->tc + t1);
/*        printf ("==> % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", t, t1, hoaa->priv->tc, dsigma,
                ncm_spline_eval_integ (hoaa->priv->nu_s, hoaa->priv->tc, t1));
*/
      }

      sigma = fmod (hoaa->priv->sigma_c + dsigma, M_PI);

      g_assert_cmpfloat (hoaa->priv->tc, !=, t1);
    
      hoaa->priv->tc      = t1;
      hoaa->priv->sigma_c = sigma;

    } while (t1 != t);
  }

  {
    const gdouble Rsigma_t0 = hoaa->priv->sigma_c;
    const gdouble Isigma_t0 = hoaa->priv->sigma_c + 0.5 * M_PI;
    gdouble nu, dlnmnu, Vnu, Nt;
    gsl_function F_F, F_G;
    NcmHOAAArg arg;

    ncm_hoaa_eval_system (hoaa, model, t, hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

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

/*
      printf ("=== \t\t % 21.15e % 21.15e\n",
              ncm_numdiff_1 (&F_F, t, fabs (t) * 0.1, &err) / dF - 1.0,
              ncm_numdiff_1 (&F_G, t, fabs (t) * 0.1, &err) / dG - 1.0
              );
  */    
      gdouble sin2RQ, cos2RQ, sin2IQ, cos2IQ;
      gdouble sin4RQ, cos4RQ, sin4IQ, cos4IQ;

      sincos (2.0 * Rsigma_t0, &sin2RQ, &cos2RQ);
      sincos (2.0 * Isigma_t0, &sin2IQ, &cos2IQ);
      sincos (4.0 * Rsigma_t0, &sin4RQ, &cos4RQ);
      sincos (4.0 * Isigma_t0, &sin4IQ, &cos4IQ);

      {
        RQ[0]  = Rsigma_t0 + 0.5 * (sin2RQ * (F + F2 + 0.5 * dG / nu) - cos2RQ * (G + F * G - 0.5 * dF / nu));
        IQ[0]  = Isigma_t0 + 0.5 * (sin2IQ * (F + F2 + 0.5 * dG / nu) - cos2IQ * (G + F * G - 0.5 * dF / nu));
        RQ[0] += 0.125 * (sin4RQ * (G2 - F2) + cos4RQ * 2.0 * F * G); 
        IQ[0] += 0.125 * (sin4IQ * (G2 - F2) + cos4IQ * 2.0 * F * G);

        sincos (2.0 * RQ[0], &sin2RQ, &cos2RQ);
        sincos (2.0 * IQ[0], &sin2IQ, &cos2IQ);
        sincos (4.0 * RQ[0], &sin4RQ, &cos4RQ);
        sincos (4.0 * IQ[0], &sin4IQ, &cos4IQ);

        RQ[0]  = Rsigma_t0 + 0.5 * (sin2RQ * (F + F2 + 0.5 * dG / nu) - cos2RQ * (G + F * G - 0.5 * dF / nu));
        IQ[0]  = Isigma_t0 + 0.5 * (sin2IQ * (F + F2 + 0.5 * dG / nu) - cos2IQ * (G + F * G - 0.5 * dF / nu));
        RQ[0] += 0.125 * (sin4RQ * (G2 - F2) + cos4RQ * 2.0 * F * G); 
        IQ[0] += 0.125 * (sin4IQ * (G2 - F2) + cos4IQ * 2.0 * F * G);
        
        RLnI[0]  = - 0.25 * (F2 + G2) - (sin2RQ * (G + F * G - 0.5 * dF / nu) + cos2RQ * (F + F2 + 0.5 * dG / nu));
        ILnI[0]  = - 0.25 * (F2 + G2) - (sin2IQ * (G + F * G - 0.5 * dF / nu) + cos2IQ * (F + F2 + 0.5 * dG / nu));
        RLnI[0] += 0.25 * (sin4RQ * 2.0 * F * G + cos4RQ * (F2 - G2));
        ILnI[0] += 0.25 * (sin4IQ * 2.0 * F * G + cos4IQ * (F2 - G2));

        /*printf ("F = % 21.15g, G = % 21.15g\n", F, G);*/
        //RLnI[0] = - (dlnmnu * sin2RQ + Vnu * cos2RQ) / (2.0 * nu);
        //ILnI[0] = - (dlnmnu * sin2IQ + Vnu * cos2IQ) / (2.0 * nu);
      }
    }
  }
}

/**
 * ncm_hoaa_eval_adiabatic_LnI_approx:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @RQ: $\theta_R$
 * @IQ: $\theta_I$
 * @RLnI: (out): $I_R$
 * @ILnI: (out): $I_I$
 * 
 * Calculates the adiabatic approximation at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_adiabatic_LnI_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble RQ, const gdouble IQ, gdouble *RLnI, gdouble *ILnI)
{
  gdouble nu, dlnmnu, Vnu, Nt;
  gsl_function F_F, F_G;
  NcmHOAAArg arg;

  ncm_hoaa_eval_system (hoaa, model, t, hoaa->k, &nu, &dlnmnu, &Vnu, &Nt);

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
    gdouble sin2RQ, cos2RQ, sin2IQ, cos2IQ;
    gdouble sin4RQ, cos4RQ, sin4IQ, cos4IQ;

    sincos (2.0 * RQ, &sin2RQ, &cos2RQ);
    sincos (2.0 * IQ, &sin2IQ, &cos2IQ);
    sincos (4.0 * RQ, &sin4RQ, &cos4RQ);
    sincos (4.0 * IQ, &sin4IQ, &cos4IQ);

    RLnI[0]  = - 0.25 * (F2 + G2) - (sin2RQ * (G + F * G - 0.5 * dF / nu) + cos2RQ * (F + F2 + 0.5 * dG / nu));
    ILnI[0]  = - 0.25 * (F2 + G2) - (sin2IQ * (G + F * G - 0.5 * dF / nu) + cos2IQ * (F + F2 + 0.5 * dG / nu));
    RLnI[0] += 0.25 * (sin4RQ * 2.0 * F * G + cos4RQ * (F2 - G2));
    ILnI[0] += 0.25 * (sin4IQ * 2.0 * F * G + cos4IQ * (F2 - G2));
  }
}

/**
 * ncm_hoaa_eval_AA:
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @RQ: (out): $\theta_R$
 * @IQ: (out): $\theta_I$
 * @RLnI: (out): $I_R$
 * @ILnI: (out): $I_I$
 * 
 * Calculates the AA variables at $t$ and $k$.
 *
 */
void 
ncm_hoaa_eval_AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *RQ, gdouble *IQ, gdouble *RLnI, gdouble *ILnI)
{
  if (!hoaa->priv->save_evol)
  {
    g_error ("ncm_hoaa_eval_adiabatic_approx: save_evol must be enables.");
  }
  else if ((t < hoaa->priv->t0) || (t > hoaa->priv->t1))
  {
    ncm_hoaa_eval_adiabatic_approx (hoaa, model, t, RQ, IQ, RLnI, ILnI);
  }
  else
  {
    RQ[0]   = ncm_spline_eval (hoaa->priv->RQ_s, t);
    IQ[0]   = ncm_spline_eval (hoaa->priv->IQ_s, t);
    RLnI[0] = ncm_spline_eval (hoaa->priv->RLnI_s, t);
    ILnI[0] = ncm_spline_eval (hoaa->priv->ILnI_s, t);
  }
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
  gdouble RQ, IQ, RLnI, ILnI;

  ncm_hoaa_eval_AA (hoaa, model, t, &RQ, &IQ, &RLnI, &ILnI);

  {
    const gdouble m    = ncm_hoaa_eval_m (hoaa, model, t, hoaa->k);
    const gdouble nu   = ncm_hoaa_eval_nu (hoaa, model, t, hoaa->k);
    const gdouble mnu  = m * nu;
    const gdouble smnu = sqrt (mnu);
    
    const gdouble q    = M_SQRT2 / smnu * exp (0.5 * RLnI) * sin (RQ);
    const gdouble Pq   = M_SQRT2 * smnu * exp (0.5 * RLnI) * cos (RQ);

    const gdouble v    = M_SQRT2 / smnu * exp (0.5 * ILnI) * sin (IQ);
    const gdouble Pv   = M_SQRT2 * smnu * exp (0.5 * ILnI) * cos (IQ);

    const complex double phi  = (q  + I * v ) / (2.0 * I);
    const complex double Pphi = (Pq + I * Pv) / (2.0 * I);

    Rphi[0]  = creal (phi);
    Iphi[0]  = cimag (phi);
    RPphi[0] = creal (Pphi);
    IPphi[0] = cimag (Pphi);
  }
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
  gdouble RQ, IQ, RLnI, ILnI;

  ncm_hoaa_eval_AA (hoaa, model, t, &RQ, &IQ, &RLnI, &ILnI);

  {
    const gdouble m    = ncm_hoaa_eval_m (hoaa, model, t, hoaa->k);
    const gdouble nu   = ncm_hoaa_eval_nu (hoaa, model, t, hoaa->k);
    const gdouble mnu  = m * nu;
    const gdouble smnu = sqrt (mnu);
    
    const gdouble q    = M_SQRT2 / smnu * exp (0.5 * RLnI) * sin (RQ);
    const gdouble Pq   = M_SQRT2 * smnu * exp (0.5 * RLnI) * cos (RQ);

    const gdouble v    = M_SQRT2 / smnu * exp (0.5 * ILnI) * sin (IQ);
    const gdouble Pv   = M_SQRT2 * smnu * exp (0.5 * ILnI) * cos (IQ);

    const complex double phi  = (q  + I * v ) / (2.0 * I);
    const complex double Pphi = (Pq + I * Pv) / (2.0 * I);
    const gdouble k3_2pi2     = gsl_pow_3 (hoaa->k) / (2.0 * gsl_pow_2 (M_PI));

    Delta_phi[0]  = k3_2pi2 * conj (phi) * phi;
    Delta_Pphi[0] = k3_2pi2 * conj (Pphi) * Pphi;
  }
}

/**
 * ncm_hoaa_eval_m: (virtual eval_m)
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
 * ncm_hoaa_eval_system: (virtual eval_system)
 * @hoaa: a #NcmHOAA
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 * @nu: (out): the frequency $\nu$
 * @dlnmnu: (out): the derivative $\mathrm{d}(m\nu)/\mathrm{d}t$
 * @Vnu: (out): the potential over frequency term $V/\nu$
 * @Nt: (out): time gauge $N_t$
 *
 * Evaluates the system functions at $t$ and $k$.
 *
 */
