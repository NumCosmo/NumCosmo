/***************************************************************************
 *            nc_matter_var.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_matter_var
 * @title: Matter Fluctuation Variance
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_matter_var.h"
#include "lss/nc_window_tophat.h"
#include "lss/nc_window_gaussian.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"
#include "math/integral.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline_cubic_notaknot.h"

#include <complex.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_expint.h>

#define NC_MATTER_VAR_SIGMA2_NP 150

G_DEFINE_TYPE (NcMatterVar, nc_matter_var, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_STRAT,
  PROP_WINDOW,
  PROP_TRANSFER
};

/**
 * nc_matter_var_new:
 * @vs: a #NcMatterVarStrategy.
 * @wp: a #NcWindow.
 * @tf: a #NcTransferFunc.
 *
 * This function allocates memory for a new #NcMatterVar object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: A new #NcMatterVar.
 */
NcMatterVar *
nc_matter_var_new (NcMatterVarStrategy vs, NcWindow *wp, NcTransferFunc *tf)
{
  NcMatterVar *vp = g_object_new (NC_TYPE_MATTER_VAR,
                                  "strategy", vs,
                                  "window", wp,
                                  "transfer", tf,
                                  NULL);
  return vp;
}

/**
 * nc_matter_var_copy:
 * @vp: a #NcMatterVar.
 *
 * This function duplicates the #NcMatterVar object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcMatterVar.
   */
NcMatterVar *
nc_matter_var_copy (NcMatterVar * vp)
{
  return nc_matter_var_new (vp->vs, vp->wp, vp->tf);
}

/**
 * nc_matter_var_free:
 * @vp: a #NcMatterVar.
 *
 * Atomically decrements the reference count of @vp by one. If the reference count drops to 0,
 * all memory allocated by @vp is released.
 *
 */
void
nc_matter_var_free (NcMatterVar *vp)
{
  g_object_unref (vp);
}

/**
 * nc_matter_var_clear:
 * @vp: a #NcMatterVar.
 *
 * Atomically decrements the reference count of @vp by one. If the reference count drops to 0,
 * all memory allocated by @vp is released. Set pointer to NULL.
 *
 */
void
nc_matter_var_clear (NcMatterVar **vp)
{
  g_clear_object (vp);
}

static void
nc_matter_var_init (NcMatterVar *vp)
{
  NcmVector *sigma2_over_growth_xv = ncm_vector_new (NC_MATTER_VAR_SIGMA2_NP);
  NcmVector *sigma2_over_growth_yv = ncm_vector_new (NC_MATTER_VAR_SIGMA2_NP);
  NcmVector *dsigma2_over_growth_xv = ncm_vector_new (NC_MATTER_VAR_SIGMA2_NP);
  NcmVector *dsigma2_over_growth_yv = ncm_vector_new (NC_MATTER_VAR_SIGMA2_NP);

  vp->r = 0.0;
  vp->spline_init = 0;
  vp->n_points_spline = 0;

  vp->integrand_overw2_spline = NULL;
  vp->sigma2_over_growth = ncm_spline_cubic_notaknot_new_full (sigma2_over_growth_xv, sigma2_over_growth_yv, FALSE);
  vp->deriv_sigma2_over_growth = ncm_spline_cubic_notaknot_new_full (dsigma2_over_growth_xv, dsigma2_over_growth_yv, FALSE);

  ncm_vector_free (sigma2_over_growth_xv);
  ncm_vector_free (sigma2_over_growth_yv);
  ncm_vector_free (dsigma2_over_growth_xv);
  ncm_vector_free (dsigma2_over_growth_yv);

  vp->ctrl = ncm_model_ctrl_new (NULL);
}

static void
_nc_matter_var_dispose (GObject * object)
{
  NcMatterVar *vp = NC_MATTER_VAR (object);

  nc_window_clear (&vp->wp);
  nc_transfer_func_clear (&vp->tf);
  ncm_model_ctrl_clear (&vp->ctrl);
  ncm_spline_clear (&vp->integrand_overw2_spline);
  ncm_spline_clear (&vp->sigma2_over_growth);
  ncm_spline_clear (&vp->deriv_sigma2_over_growth);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_matter_var_parent_class)->dispose (object);
}

static void
_nc_matter_var_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_matter_var_parent_class)->finalize (object);
}

static void
_nc_matter_var_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMatterVar *vp = NC_MATTER_VAR (object);
  g_return_if_fail (NC_IS_MATTER_VAR (object));

  switch (prop_id)
  {
    case PROP_STRAT:
      vp->vs = g_value_get_enum (value);
      break;
    case PROP_WINDOW:
      vp->wp = g_value_dup_object (value);
      break;
    case PROP_TRANSFER:
      vp->tf = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_matter_var_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMatterVar *vp = NC_MATTER_VAR (object);
  g_return_if_fail (NC_IS_MATTER_VAR (object));

  switch (prop_id)
  {
    case PROP_STRAT:
      g_value_set_enum (value, vp->vs);
      break;
    case PROP_WINDOW:
      g_value_set_object (value, vp->wp);
      break;
    case PROP_TRANSFER:
      g_value_set_object (value, vp->tf);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_matter_var_class_init (NcMatterVarClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_matter_var_dispose;
  object_class->finalize = _nc_matter_var_finalize;
  object_class->set_property = _nc_matter_var_set_property;
  object_class->get_property = _nc_matter_var_get_property;

  /**
   * NcMatterVar:strategy:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_STRAT,
                                   g_param_spec_enum ("strategy",
                                                      NULL,
                                                      "Strategy",
                                                      NC_TYPE_MATTER_VAR_STRATEGY, NC_MATTER_VAR_FFT,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                      | G_PARAM_STATIC_BLURB));

  /**
   * NcMatterVar:window:
   *
   * This property keeps the window object.
   */
  g_object_class_install_property (object_class,
                                   PROP_WINDOW,
                                   g_param_spec_object ("window",
                                                        NULL,
                                                        "Window Function.",
                                                        NC_TYPE_WINDOW,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));

  /**
   * NcMatterVar:transfer:
   *
   * This property keeps the transferfunc object.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER,
                                   g_param_spec_object ("transfer",
                                                        NULL,
                                                        "Transfer Function.",
                                                        NC_TYPE_TRANSFER_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));
}

static void _nc_matter_var_prepare_numint (NcMatterVar *vp, NcHICosmo *model);
static void _nc_matter_var_prepare_splineint (NcMatterVar *vp, NcHICosmo *model);
static void _nc_matter_var_prepare_fft (NcMatterVar *vp, NcHICosmo *model);

/**
 * nc_matter_var_prepare:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 *
 * FIXME
 */
void
nc_matter_var_prepare (NcMatterVar *vp, NcHICosmo *model)
{
  switch (vp->vs)
  {
    case NC_MATTER_VAR_NUMINT:
      _nc_matter_var_prepare_numint (vp, model);
      break;
    case NC_MATTER_VAR_SPLINEINT:
      _nc_matter_var_prepare_splineint (vp, model);
      break;
    case NC_MATTER_VAR_FFT:
#ifdef NUMCOSMO_HAVE_FFTW3
      _nc_matter_var_prepare_fft (vp, model);
#else
      g_error ("nc_matter_var_prepare: Cannot use NC_MATTER_VAR_FFT: fftw not installed.");
#endif
      break;
  }
}

/* Estrutura criada para podermos calcular todos os momentos espectrais, a variancia eh o momento para j=0. */

typedef struct _int_temp_stc
{
  int j;
  gdouble R;
  NcMatterVar *vp;
  NcHICosmo *model;
} int_temp_stc;

static gdouble
_nc_matter_var_integrand_gaussian (gdouble kR, gpointer params)     /* integrando de \sigma^2 sem a funcao crescimento */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp; //(NcVariance *) params;
  gdouble k = kR / ts->R;        /* por causa da integracao, para altos valores de R temos problema */
  const gdouble c1 = 2.0 * M_PI * M_PI;    /* 2\pi^2 */
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble W = nc_window_eval_fourier (vp->wp, k, ts->R);
  gdouble k2 = k*k;
  gdouble v_integrand = matter_P * k2 * gsl_pow_int (k2, ts->j) * W * W / (c1 * ts->R);  /* O termo 1/R é para a transformacao de variavel k -> kR. A integral eh feita em kR. */

  //printf ("%g %g %g %g %g %g\n", P, k2, gsl_pow_int (k2, ts->j), T, W, vp->wp->R);

  return v_integrand;
}

static gdouble
_nc_matter_var_over_growth2_gaussian (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)     /* \frac{\sigma^2}{D^2} */
{
  gdouble v_over_growth2_gaussian, error;
  static gsl_integration_workspace *w = NULL; //gsl_integration_workspace_alloc (INT_PARTITION);
  gsl_function F;
  int_temp_stc ts;

  ts.j = 0;
  ts.R = exp (lnR);
  ts.vp = vp;
  ts.model = model;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

  F.function = &_nc_matter_var_integrand_gaussian;
  F.params = &ts;

  gsl_integration_qagiu (&F, 0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &v_over_growth2_gaussian, &error);

  return v_over_growth2_gaussian;

}

static gdouble
_nc_matter_var_integrand_tophat_ksmall (gdouble kR, gpointer params)    /* integrando de \sigma^2 (sem a funcao crescimento) completo para se fazer a integral de k=0 a 1  */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  const gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble W = nc_window_eval_fourier (vp->wp, k, ts->R);
  gdouble k2 = k*k;
  gdouble v_integrand = matter_P * k2 * gsl_pow_int (k2, ts->j) * W * W / (c1 * ts->R);  /* O termo 1/R X para a transformacao de variavel k -> kR. */

  return v_integrand;
}

static gdouble
_nc_matter_var_integrand_tophat_one (gdouble kR, gpointer params)     /* integrando de \sigma^2 sem a funcao crescimento */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  const gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble k2 = k * k;
  gdouble k4 = k2 * k2;
  gdouble v_integrand_one = 4.5 * matter_P * gsl_pow_int (k2, ts->j) * (1.0 + k2 * ts->R * ts->R) / (gsl_pow_6(ts->R) * k4 * c1 * ts->R);

  return v_integrand_one;
}

static gdouble
_nc_matter_var_integrand_tophat_cosine (gdouble kR, gpointer params)     /* integrando de \sigma^2 sem a funcao crescimento */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  const gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble k2 = k*k;
  gdouble k4 = k2*k2;
  gdouble v_integrand_2 = 4.5 * matter_P * gsl_pow_int (k2, ts->j) * (-1.0 + k2*ts->R*ts->R)/ (gsl_pow_6(ts->R) * k4 * c1 * ts->R);  /* *cos(2kR) */

  return v_integrand_2;
}

static gdouble
_nc_matter_var_integrand_tophat_sine (gdouble kR, gpointer params)     /* integrando de \sigma^2 sem a funcao crescimento */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  const gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble k2 = k*k;
  gdouble k3 = k2*k;
  gdouble v_integrand_3 = -9.0 * matter_P * gsl_pow_int (k2, ts->j) / (gsl_pow_5(ts->R) * k3 * c1 * ts->R);   /* *sin(2kR) */

  return v_integrand_3;
}

#define _NC_LSTEP2 5.0

static gdouble
_nc_matter_var_over_growth2_tophat_old (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)     /* \frac{\sigma^2}{D^2} */
{
  gdouble v_over_growth2_ksmall = 0.0, v_over_growth2_one = 0.0, v_over_growth2_osc = 0.0;
  gdouble error, error_1, error_cos, error_sin, step, pres, result = 0.0;
  static gsl_integration_workspace *w = NULL; //gsl_integration_workspace_alloc (INT_PARTITION);
  static gsl_integration_qawo_table *wosc_cos = NULL;
  static gsl_integration_qawo_table *wosc_sin = NULL;
  gsl_function F, F_cos, F_sin, Ft;
  int_temp_stc ts;

  //  int i;
  //  gsl_spline *var_spline;

  ts.j = 0;
  ts.R = exp (lnR);
  ts.vp = vp;
  ts.model = model;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);
  if (wosc_cos == NULL)
    wosc_cos = gsl_integration_qawo_table_alloc (2.0, _NC_LSTEP2, GSL_INTEG_COSINE, NCM_INTEGRAL_PARTITION);
  if (wosc_sin == NULL)
    wosc_sin = gsl_integration_qawo_table_alloc (2.0, _NC_LSTEP2, GSL_INTEG_SINE, NCM_INTEGRAL_PARTITION);

  F.function = &_nc_matter_var_integrand_tophat_one;
  F_cos.function = &_nc_matter_var_integrand_tophat_cosine;
  F_sin.function = &_nc_matter_var_integrand_tophat_sine;
  Ft.function = &_nc_matter_var_integrand_tophat_ksmall;

  F.params = &ts;
  F_cos.params = &ts;
  F_sin.params = &ts;
  Ft.params = &ts;
  //  printf ("# Old code-------------\n");

  //  nc_matter_var_integrand_prepare_spline (vp);
  //  var_spline =  vp->integrand_overw2_spline;

  //  for (i = 0; i < var_spline->size-1; i++)
  //  {
  //    gsl_integration_qag (&Ft, var_spline->x[i] * vp->wp->R, var_spline->x[i+1] * vp->wp->R, 0, 1e-11, INT_PARTITION, 6, w, &v_over_growth2_ksmall, &error);
  gsl_integration_qag (&Ft, 0.0, 1.0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &v_over_growth2_ksmall, &error);   /* integral completa de 0 a 1 */

  result += v_over_growth2_ksmall;

  //  }

  gsl_integration_qag (&F, 1.0, 1000.0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &v_over_growth2_one, &error_1);
  result += v_over_growth2_one;

  for (step = 1.0; 1; step += _NC_LSTEP2)
  {
    gdouble pres_t = 0.0;
    gsl_integration_qawo (&F_cos, step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, wosc_cos, &pres, &error_cos);
    pres_t = pres;

    gsl_integration_qawo (&F_sin, step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, wosc_sin, &pres, &error_sin);
    pres_t += pres;

    if (fabs(pres_t/v_over_growth2_osc) < NCM_DEFAULT_PRECISION * 1e-3)
      break;
    v_over_growth2_osc += pres_t;
    if (pres_t == 0.0)
      break;
  }

  result += v_over_growth2_osc;

  return result;
}

typedef gdouble (* _NcMatterVarVar0) (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
typedef gdouble (* _NcMatterVardVar0dR) (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
static gdouble _nc_matter_var_over_growth2_tophat_old (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
static gdouble _nc_matter_var_dvariance_over_growth2_dR_tophat (NcMatterVar *vp, NcHICosmo *model, gdouble lnR); /* Check if it is R or lnR*/
static gdouble _nc_matter_var_over_growth2_gaussian (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
static gdouble _nc_matter_var_dvariance_over_growth2_dR_gaussian (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);

static void
_nc_matter_var_prepare_numint (NcMatterVar *vp, NcHICosmo *model)
{
  gint i;

  nc_transfer_func_prepare (vp->tf, model);

  {
    _NcMatterVarVar0 var0 = NULL;
    _NcMatterVardVar0dR dvar0dR = NULL;

    if (NC_IS_WINDOW_TOPHAT (vp->wp))
    {
      var0 = &_nc_matter_var_over_growth2_tophat_old;
      dvar0dR = &_nc_matter_var_dvariance_over_growth2_dR_tophat;
    }
    else if (NC_IS_WINDOW_GAUSSIAN (vp->wp))
    {
      var0 = &_nc_matter_var_over_growth2_gaussian;
      dvar0dR = &_nc_matter_var_dvariance_over_growth2_dR_gaussian;
    }
    else
      g_assert_not_reached ();


    for (i = 0; i < NC_MATTER_VAR_SIGMA2_NP; i++)
    {
      gdouble lnR = log (0.1) + log (10000.0) / (NC_MATTER_VAR_SIGMA2_NP - 1.0) * i;
      //gdouble R = exp (lnR);
      gdouble sigma2_0, deriv_sigma2_0; /* variance over growth2 */
      sigma2_0 = var0 (vp, model, lnR);
      deriv_sigma2_0 = dvar0dR (vp, model, lnR);

      ncm_vector_set (vp->sigma2_over_growth->xv, i, lnR);
      ncm_vector_set (vp->sigma2_over_growth->yv, i, log(sigma2_0));
      ncm_vector_set (vp->deriv_sigma2_over_growth->xv, i, lnR);
      ncm_vector_set (vp->deriv_sigma2_over_growth->yv, i, (exp(lnR) / sigma2_0) * deriv_sigma2_0);
    }
    ncm_spline_prepare (vp->sigma2_over_growth);
    ncm_spline_prepare (vp->deriv_sigma2_over_growth);
  }
  return;
}

static void _top_hat_cubic_spline_integration_rule (gdouble x, gdouble *rules);
static void _derivative_top_hat_cubic_spline_integration_rule (gdouble x, gdouble *rules);
static void _top_hat_cubic_spline_integration_rule_taylor (gdouble xa, gdouble xb, gdouble R, gdouble *rules);
//static void _derivative_top_hat_cubic_spline_integration_rule_taylor (gdouble xa, gdouble xb, gdouble R, gdouble *rules);

/**
 * nc_matter_var_var0:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnR: logarithm base e of the radius.
 *
 * This function returns the variance of the density contrast at redshift \f$ z = 0 \f$ computed at scale R FIXME
 *
 * Returns: a gdouble which is the variance \f$ \sigma^2 (R, z = 0) \f$.
 */
gdouble
nc_matter_var_var0 (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)
{
  //const gdouble lnR = log (R);
  if (ncm_model_ctrl_update (vp->ctrl, NCM_MODEL(model)))
    nc_matter_var_prepare (vp, model);
  return exp (ncm_spline_eval (vp->sigma2_over_growth, lnR));
}

/**
 * nc_matter_var_dlnvar0_dR:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnR: logarithm base e of the radius.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_dlnvar0_dR (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)
{
  //const gdouble lnR = log (R);
  if (ncm_model_ctrl_update (vp->ctrl, NCM_MODEL(model)))
    nc_matter_var_prepare (vp, model);
  return ncm_spline_eval (vp->deriv_sigma2_over_growth, lnR) / exp(lnR);
}

/**
 * nc_matter_var_dlnvar0_dlnR:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnR: logarithm base e of the radius.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_dlnvar0_dlnR (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)
{
  if (ncm_model_ctrl_update (vp->ctrl, NCM_MODEL(model)))
    nc_matter_var_prepare (vp, model);
  return ncm_spline_eval (vp->deriv_sigma2_over_growth, lnR);
}

/**
 * nc_matter_var_mass_to_R:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @M: mass enclosed in the volume specified by the window function.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_mass_to_R (NcMatterVar *vp, NcHICosmo *model, gdouble M)
{
  const gdouble Omega_m = nc_hicosmo_Omega_m (model);
  return cbrt (M / (Omega_m * nc_window_volume(vp->wp) * ncm_c_crit_mass_density_solar_Mpc ()));
}

/**
 * nc_matter_var_R_to_mass:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @R: radius.
 *
 * FIXME mass enclosed in the volume specified by the window function
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_R_to_mass (NcMatterVar *vp, NcHICosmo *model, gdouble R)
{
  const gdouble Omega_m = nc_hicosmo_Omega_m (model);
  return gsl_pow_3(R) * Omega_m * nc_window_volume(vp->wp) * ncm_c_crit_mass_density_solar_Mpc ();
}

/**
 * nc_matter_var_lnM_to_lnR:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of the mass enclosed in the volume specified by the window function.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_lnM_to_lnR (NcMatterVar *vp, NcHICosmo *model, gdouble lnM)
{
  const gdouble Omega_m = nc_hicosmo_Omega_m (model);
  return (lnM - log (Omega_m * nc_window_volume(vp->wp) * ncm_c_crit_mass_density_solar_Mpc ())) / 3.0;
}

/**
 * nc_matter_var_lnR_to_lnM:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnR: logarithm base e of the radius.
 *
 * FIXME mass enclosed in the volume specified by the window function
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_lnR_to_lnM (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)
{
  const gdouble Omega_m = nc_hicosmo_Omega_m (model);
  return (3.0 * lnR + log (Omega_m * nc_window_volume(vp->wp) * ncm_c_crit_mass_density_solar_Mpc ()));
}

/**
 * nc_matter_var_integrand_over_window2:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @k: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_integrand_over_window2 (NcMatterVar *vp, NcHICosmo *model, gdouble k)
{
  gdouble k2, integrand_overw2;
  k2 = k * k;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, model, k);
  integrand_overw2 = k2 * matter_P / (2.0 * M_PI * M_PI);

  return integrand_overw2;
}

static gdouble nc_matter_var_integrate_spline (NcMatterVar *vp, NcmSpline *s, void (*calc_rule)(gdouble, gdouble *), void (*shift_rule)(gdouble *, gdouble *, gdouble, gdouble, gdouble *), void (*calc_rule_taylor) (gdouble, gdouble, gdouble, gdouble *), gdouble R);
static void nc_matter_var_shift_scale_rule (gdouble *rules_a, gdouble *rules_b, gdouble R, gdouble a, gdouble *rule);
static void nc_matter_var_deriv_shift_scale_rule (gdouble *rules_a, gdouble *rules_b, gdouble R, gdouble a, gdouble *rule);

typedef struct __NcMatterVarIntegOverW2
{
  NcMatterVar *vp;
  NcHICosmo *model;
} _NcMatterVarIntegOverW2;

static gdouble
_nc_matter_var_integrand_over_window2 (gdouble k, gpointer data)
{
  _NcMatterVarIntegOverW2 *vpd = (_NcMatterVarIntegOverW2 *) data;
  gdouble k2, integrand_overw2;
  k2 = k * k;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vpd->vp->tf, vpd->model, k);
  integrand_overw2 = k2 * matter_P / (2.0 * M_PI * M_PI);

  return integrand_overw2;
}

/**
 * nc_matter_var_integrand_prepare_spline:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 *
 * This function prepares the spline for \f$ k^2 P(k) T^2(k)\f$, where P(k) is the power spectrum and T(k) is the
   * transfer function. Given this spline, we optimize the computation of the variance of the density contrast (comparing
                                                                                                                * to nc_matter_var_over_growth2_tophat, for example). The variance with the gaussian window funtion is given by
     * \f$ \sigma^2(R, z) = A \left(\frac{D(z)}{D(0)}\right)^2 \int_0^\infty \frac{dk}{2\pi^2} k^2 P(k) T^2(k) \exp(-(kR)^2) \f$,
     * where \f$ A \f$ is the power spectrum normalization and \f$ D(z) \f$ is the growth function at redshift \f$ z \f$.
       * The variance for the top hat window function takes the form
 * \f$ \sigma^2(R, z) = A \left(\frac{D(z)}{D(0)}\right)^2 \int_0^\infty dk \frac{9}{2\pi^2(kR)^2} k^2 P(k) T^2(k) (j_1(kR))^2 \f$,
 * where \f$ j_1(kR) \f$ is the spherical Bessel function of the first kind.
   */
static void
_nc_matter_var_prepare_splineint (NcMatterVar *vp, NcHICosmo *model)
{
  gint i;
  //  GTimer *bench = g_timer_new ();

  if (TRUE)
  {
    gint last_i;
    gdouble first_k[] = {0.0, 1.0e-20, 1.0e-19, 1.0e-18, 1.0e-17, 1.0e-16, 1.0e-15, 1.0e-14, 1.0e-13, 1.0e-12, 1.0e-11,
      1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5};
    vp->n_points_spline = 230;

    if (vp->integrand_overw2_spline == NULL)
    {
      NcmVector *integrand_overw2_xv = ncm_vector_new (vp->n_points_spline);
      NcmVector *integrand_overw2_yv = ncm_vector_new (vp->n_points_spline);

      vp->integrand_overw2_spline =
        ncm_spline_cubic_notaknot_new_full (integrand_overw2_xv,
                                            integrand_overw2_yv,
                                            FALSE);
      ncm_vector_free (integrand_overw2_xv);
      ncm_vector_free (integrand_overw2_yv);
    }

    vp->spline_init = TRUE;

    nc_transfer_func_prepare (vp->tf, model);

    for (i = 0; i < 17; i++)
    {
      const gdouble k = first_k[i];
      ncm_vector_set (vp->integrand_overw2_spline->xv, i, k);
      ncm_vector_set (vp->integrand_overw2_spline->yv, i, nc_matter_var_integrand_over_window2 (vp, model, k));
    }
    last_i = i;

    for (; i < last_i + 66; i++)
    {
      const gdouble k = 4e-2 / (66.0 - 1.0) * (i - last_i) + 2.0e-5;
      ncm_vector_set (vp->integrand_overw2_spline->xv, i, k);
      ncm_vector_set (vp->integrand_overw2_spline->yv, i, nc_matter_var_integrand_over_window2 (vp, model, k));
    }
    last_i = i;

    for (; i < last_i + 76; i++)
    {
      const gdouble k = 0.75 / (76.0 - 1.0) * (i - last_i) + 5e-2;
      ncm_vector_set (vp->integrand_overw2_spline->xv, i, k);
      ncm_vector_set (vp->integrand_overw2_spline->yv, i, nc_matter_var_integrand_over_window2 (vp, model, k));
    }
    last_i = i;

    for (; i < last_i + 51; i++)
    {
      const gdouble k = exp (log(0.85) + (log(400.0) - log(0.85)) * (i - last_i) / 50.0);
      ncm_vector_set (vp->integrand_overw2_spline->xv, i, k);
      ncm_vector_set (vp->integrand_overw2_spline->yv, i, nc_matter_var_integrand_over_window2 (vp, model, k));
    }
    last_i = i;

    for (; i < last_i + 10; i++)
    {
      const gdouble k = 400.0 / (10.0 - 1.0) * (i - last_i) + 450.0;
      ncm_vector_set (vp->integrand_overw2_spline->xv, i, k);
      ncm_vector_set (vp->integrand_overw2_spline->yv, i, nc_matter_var_integrand_over_window2 (vp, model, k));
    }
    last_i = i;

    for (; i < last_i + 10; i++)
    {
      const gdouble k = 100.0 / (10.0 - 1.0) * (i - last_i) + 900.0;
      ncm_vector_set (vp->integrand_overw2_spline->xv, i, k);
      ncm_vector_set (vp->integrand_overw2_spline->yv, i, nc_matter_var_integrand_over_window2 (vp, model, k));
    }

    /* last_i = i; */
    ncm_spline_prepare (vp->integrand_overw2_spline);
  }
  else if (FALSE)
  {
    gsl_function F;
    _NcMatterVarIntegOverW2 vpd = {vp, model};
    F.function = &_nc_matter_var_integrand_over_window2;
    F.params = &vpd;
    vp->integrand_overw2_spline = ncm_spline_cubic_notaknot_new ();
    ncm_spline_set_func (vp->integrand_overw2_spline, NCM_SPLINE_FUNCTION_SPLINE_LNKNOT, &F, 1e-20, 1000.0, 0, 1e-5);
    printf ("Ok Spline %zd\n", vp->integrand_overw2_spline->len);fflush (stdout);
  }

  if (FALSE)
  {
    /* Testing precision */
    gdouble max_err = 0.0;
    for (i = 0; i <= 10000; i++)
    {
      gdouble k = pow (10.0, -3.0 + 6.0 / (10000.0) * i);
      gdouble f = nc_matter_var_integrand_over_window2 (vp, model, k);
      gdouble fs = ncm_spline_eval (vp->integrand_overw2_spline, k);
      gdouble err = fabs ((f-fs) / f);
      max_err = GSL_MAX (err, max_err);
      printf ("% 20.15g % 20.15g % 20.15g % 8.5e % 8.5e\n", k, f, fs, err, max_err);
    }
  }

  for (i = 0; i < NC_MATTER_VAR_SIGMA2_NP; i++)
  {
    gdouble lnR = log (0.1) + log (1.0e4) / (NC_MATTER_VAR_SIGMA2_NP - 1.0) * i;
    gdouble R = exp (lnR);
    gdouble sigma2_0, deriv_sigma2_0; /* variance over growth2 */

    //gdouble sigma2_0_exact = nc_matter_var_over_growth2_tophat_old (vp, model, lnR);
    //gdouble deriv_sigma2_0_exact = _nc_matter_var_dvariance_over_growth2_dR_tophat (vp, model, R);

    sigma2_0 = nc_matter_var_integrate_spline (vp, vp->integrand_overw2_spline, &_top_hat_cubic_spline_integration_rule, nc_matter_var_shift_scale_rule, &_top_hat_cubic_spline_integration_rule_taylor, R);
    deriv_sigma2_0 = nc_matter_var_integrate_spline (vp, vp->integrand_overw2_spline, &_derivative_top_hat_cubic_spline_integration_rule, &nc_matter_var_deriv_shift_scale_rule, NULL, R);

    ncm_vector_set (vp->sigma2_over_growth->xv, i, lnR);
    ncm_vector_set (vp->sigma2_over_growth->yv, i, log(sigma2_0));
    ncm_vector_set (vp->deriv_sigma2_over_growth->xv, i, lnR);
    ncm_vector_set (vp->deriv_sigma2_over_growth->yv, i, (R / sigma2_0) * deriv_sigma2_0); /* dln(var0)/dlnR */

    //printf ("%d % 20.15g % 20.15g % 20.15g\n", i, lnR, exp(lnR) * sigma2_0, 0.0);
    //printf("% 20.15g % 20.15g % 20.15g\n", sigma2_0_exact, sigma2_0, (sigma2_0_exact - sigma2_0) / sigma2_0_exact);
    //printf("% 20.5e % 20.15g % 20.15g % 20.15g\n", nc_matter_var_R_to_mass (vp, model, R), deriv_sigma2_0_exact, deriv_sigma2_0, (deriv_sigma2_0_exact - deriv_sigma2_0) / deriv_sigma2_0_exact);
  }
  ncm_spline_prepare (vp->sigma2_over_growth);
  ncm_spline_prepare (vp->deriv_sigma2_over_growth);

  return;
}

#ifdef NUMCOSMO_HAVE_FFTW3
static void
_nc_matter_var_prepare_fft (NcMatterVar *vp, NcHICosmo *model)
{
  const gint N = 500 + 1;
  const gint N_2 = 250;
  const gdouble L = 17.0 * M_LN10;
  const gdouble dr = L / (N * 1.0);
  const gdouble r0 = 3.0 * M_LN10;
  const gdouble lnk0 = -3.0 * M_LN10;
  static fftw_complex *in = NULL;
  static fftw_complex *out = NULL;
  static fftw_complex *u = NULL;
  static fftw_complex *fftRdsigma2_dr = NULL;
  static fftw_complex *Rdsigma2_dr = NULL;
  static fftw_complex *du = NULL;
  static fftw_plan p, p_Rsigma2, p_Rdsigma2_dr;
  static gdouble planned = FALSE;
  static gdouble calc_u = FALSE;
  gint i, j;

  //fprintf (stderr, "# %d % 20.15g [% 20.15g]\n", N, dr, L);
  if (!planned)
  {
    in             = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * N);
    out            = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * N);
    fftRdsigma2_dr = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * N);
    Rdsigma2_dr    = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * N);

    u  = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * N);
    du = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * N);

    ncm_cfg_load_fftw_wisdom ("nc_matter_var_wisdown.fftw3");
    p = fftw_plan_dft_1d (N, in, out, FFTW_FORWARD, FFTW_PATIENT | FFTW_DESTROY_INPUT);
    p_Rsigma2 = fftw_plan_dft_1d (N, out, in, FFTW_FORWARD, FFTW_PATIENT | FFTW_DESTROY_INPUT);
    p_Rdsigma2_dr = fftw_plan_dft_1d (N, fftRdsigma2_dr, Rdsigma2_dr, FFTW_FORWARD, FFTW_PATIENT | FFTW_DESTROY_INPUT);
    ncm_cfg_save_fftw_wisdom ("nc_matter_var_wisdown.fftw3");
    
    planned = TRUE;
  }

  for (i = -N_2; i <= N_2; i++)
  {
    const gint ii = (i < 0) ? i + N : i;
    const gdouble c1 = 2.0 * M_PI * M_PI;
    const gdouble mu = dr * i;
    const gdouble k = exp (lnk0 + mu);
    const gdouble k2 = k * k;
    const gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, model, k);
    const gdouble f = matter_P * k2 / c1;
    in[ii] = f;
  }
  fftw_execute (p);

  if (!calc_u)
  {
    for (i = -N_2; i <= N_2; i++)
    {
      const gint ii = (i < 0) ? i + N : i;
      const gdouble a = 2.0 * M_PI / L * i;
      const gdouble abs_a = fabs (a);
      const gdouble sign_a = a < 0 ? -1.0 : 1.0;
      complex double U = 36.0 * (a + I) / (cexp (M_LN2 * I * a) * (a * I - 5.0));
      gsl_sf_result lnr;
      gsl_sf_result arg;
      gsl_sf_lngamma_complex_e (-3.0, a, &lnr, &arg);
      U = (a != 0) ? U * sign_a * cexp (lnr.val + arg.val * I + gsl_sf_lnsinh (M_PI * abs_a * 0.5)) : 3.0 * M_PI / 5.0;
      u[ii] = U * cexp (-2.0 * M_PI / L * i * I * (lnk0 + r0));
      du[ii] = -(1.0 + I * a) * u[ii];
      fftRdsigma2_dr[ii] = out[ii] * du[ii];
      out[ii] *= u[ii];
    }
    calc_u = TRUE;
  }
  else
  {
    gint ii;
    for (ii = 0; ii < N; ii++)
    {
      fftRdsigma2_dr[ii] = out[ii] * du[ii];
      out[ii] *= u[ii];
    }
  }
  fftw_execute (p_Rsigma2);
  fftw_execute (p_Rdsigma2_dr);

  j = 0;
  for (i = -120; i < 30; i++)
  {
    const gint ii = (i < 0) ? i + N : i;
    const gdouble r = r0 + i * dr;
    const complex double Rsigma2 = in[ii] / N;
    const complex double lnsigma2 = log (Rsigma2) - r;
    const complex double dlnsigma2_dr = Rdsigma2_dr[ii] / (Rsigma2 * N);
    /*
     const gdouble R = exp (r);
     const gdouble Rsigma2_old = R * nc_matter_var_over_growth2_tophat_old (vp, model, r);
     const gdouble Rsigma2_spline = R * exp(ncm_spline_eval (vp->sigma2_over_growth, r));
     const gdouble err_spline = fabs((Rsigma2_old - Rsigma2_spline) / Rsigma2_old);
     const gdouble err_fft = fabs((Rsigma2_old - creal(Rsigma2)) / Rsigma2_old);
     printf ("%d %d % 20.15g % 20.15g % 20.15g % 20.15g % 8.5e % 8.5e\n", i, ii, r / M_LN10, creal(Rsigma2), Rsigma2_spline, Rsigma2_old, err_fft, err_spline);
     */
    //printf ("%d %d % 20.15g % 20.15g % 20.15g\n", i, ii, r / M_LN10, creal(Rsigma2), creal(dlnsigma2_dr));
    ncm_vector_set (vp->sigma2_over_growth->xv, j, r);
    ncm_vector_set (vp->sigma2_over_growth->yv, j, lnsigma2);
    ncm_vector_set (vp->deriv_sigma2_over_growth->xv, j, r);
    ncm_vector_set (vp->deriv_sigma2_over_growth->yv, j, dlnsigma2_dr);
    j++;
  }
  ncm_spline_prepare (vp->sigma2_over_growth);
  ncm_spline_prepare (vp->deriv_sigma2_over_growth);
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

/**
 * _2f3:
 * @n: FIXME
 * @x: FIXME
 * @prec: is the precision required to compute this function.
 *
 * Hypergeometric of type /$ {}_2F_3 /$.
 *
 * Returns: FIXME
 */
static gdouble
_2f3 (int n, gdouble x, gdouble prec)
{
  gdouble fact = 1.0;
  gdouble sum = 1.0;
  gdouble mx2 = - x * x;
  unsigned long int i = 0;

  while (TRUE)
  {
    fact *= (2.0 + i) * ((n + 1.0) / 2.0 + i) * mx2;
    fact /= (5.0 / 2.0 + i) * (4.0 + i) * ((n + 3.0)/2.0 + i) * (i + 1.0);

    if (fabs(fact/sum) < prec)
      break;

    sum += fact;

    i++;
  }

  return gsl_pow_int (x, n + 1) * sum / (9.0 * (n + 1.0));
}

/**
 * _3f4:
 * @n: FIXME
 * @x: FIXME
 * @prec: is the precision required to compute this function.
 *
 * Hypergeometric of type /$ {}_3F_4 /$.
 *
 * Returns: FIXME
 */
static gdouble
_3f4 (int n, gdouble x, gdouble prec)
{
  gdouble fact = 1.0;
  gdouble sum = 1.0;
  gdouble mx2 = - x * x;
  unsigned long int i = 0;

  while (TRUE)
  {
    fact *= (5.0 / 2.0 + i) * (3.0 + i) * ((n + 3.0) / 2.0 + i) * mx2;
    fact /= (7.0 / 2.0 + i) * (5.0 / 2.0 + i) * (5.0 + i) * ((n + 5.0)/2.0 + i) * (i + 1.0);

    if (fabs(fact/sum) < prec)
      break;

    sum += fact;

    i++;
  }

  return gsl_pow_int (x, n + 3) * sum / (45.0 * (n + 3.0));
}

static void
_top_hat_cubic_spline_integration_rule (gdouble x, gdouble *rules)
{
  if (x < 1.0)
  {
    rules[0] = 9.0 * _2f3 (0, x, GSL_DBL_EPSILON); // limit when x-> 0: 0.0;
    rules[1] = 9.0 * _2f3 (1, x, GSL_DBL_EPSILON); // limit when x-> 0: 9.0*(-0.25);
    rules[2] = 9.0 * _2f3 (2, x, GSL_DBL_EPSILON); // limit when x-> 0: 0.0;
    rules[3] = 9.0 * _2f3 (3, x, GSL_DBL_EPSILON); // limit when x-> 0: 9.0*(1.0 - M_LN2 - M_EULER)/2.0;
  }
  else
  {
    const gdouble x2 = x * x;   /* x = k * R */
    const gdouble x3 = x2 * x;
    const gdouble x4 = x3 * x;
    const gdouble x5 = x4 * x;
    const gdouble x_2 = 2.0 * x;
    const gdouble cos_x = cos(x);
    const gdouble cos_x_2 = cos_x * cos_x;
    const gdouble sin_x = sin(x);
    const gdouble sin_x_2 = sin_x * sin_x;
    const gdouble cos_2x = cos_x_2 - sin_x_2;
    const gdouble sin_2x = 2.0 * cos_x * sin_x;
    const gdouble sin_integral_2x = gsl_sf_Si (x_2);
    const gdouble cos_integral_2x = gsl_sf_Ci (x_2);
    rules[0] = 9.0 * (-3.0 -5.0 * x2 + (3.0 - x2 + 2.0 * x4) * cos_2x + x * (6.0 + x2) * sin_2x + 4.0 * x5 * sin_integral_2x) / (30.0 * x5);
    rules[1] = 9.0 * ((-1.0 -2.0 * x2 + cos_2x + 2.0 * x * sin_2x) / (8.0 * x4) + 0.25);
    rules[2] = 9.0 * (-1.0 -3.0 * x2 + (1.0 + x2) * cos_2x + 2.0 * x * sin_2x + 2.0 * x3 * sin_integral_2x) / (6.0 * x3);
    rules[3] = 9.0 * ((-sin_x_2 + x2 * log(x) + x * (sin_2x - x * cos_integral_2x))/ (2.0 * x2) - (1.0 - M_LN2 - M_EULER)/2.0);
  }
}

static void
_top_hat_cubic_spline_integration_derivatives_rule (gdouble x, gdouble *dnj12_x2)
{
  if (x == 0.0)
  {
    dnj12_x2[0] = 1.0;
    dnj12_x2[1] = 0.0;
    dnj12_x2[2] = -2.0 / 5.0;
    dnj12_x2[3] = 0.0;
    dnj12_x2[4] = 72.0 / 175.0;
  }
  else
  {
    const gdouble x2 = x * x;   /* x = k * R */
    const gdouble x3 = x2 * x;
    const gdouble x4 = x2 * x2;
    const gdouble x5 = x3 * x2;
    const gdouble j1 = gsl_sf_bessel_j1 (x);
    const gdouble j12 = j1 * j1;
    const gdouble j2 = gsl_sf_bessel_j2 (x);
    const gdouble j22 = j2 * j2;
    const gdouble j1j2 = j1 * j2;

    dnj12_x2[0] = 9.0 * j12 / x2;
    dnj12_x2[1] = -18.0 * j1j2 / x2;
    dnj12_x2[2] = 18.0 * ((- j12 / x2) + (4.0 * j1j2 / x3) + (j22 / x2));
    dnj12_x2[3] = 72.0 * (x * j12 + (x2 - 5.0) * j1j2 - 3.0 * x * j22) / x4;
    dnj12_x2[4] = 72.0 * ((x3 - 6.0 * x) * j12 + (30.0 - 12.0 * x2) * j1j2 + (32.0 * x - x3) * j22) / x5;
  }
}

/* To compute dsigma^2/dR
 static void
 _derivative_top_hat_cubic_spline_integration_derivatives_rule (gdouble x, gdouble *drule, gdouble *d2rule, gdouble *d3rule)
 {
   const gdouble x2 = x * x;   // x = k * R
   const gdouble j1 = gsl_sf_bessel_j1 (x);
   const gdouble j12 = j1 * j1;
   const gdouble j2 = gsl_sf_bessel_j2 (x);
   const gdouble j22 = j2 * j2;
   const gdouble j1j2 = j1 * j2;

   drule[0] = -18.0 * j12 / x;
   d2rule[0] = -18.0 * (x * (j12 - j22) - 3.0 * j1j2) / x2;
   d3rule[0] = -36.0 * (-j12 / x - 2.0 * j1j2 + 6.0 * j1j2 / x2 + 5.0 * j22 / x) / x;

   drule[1] = -18.0 * j1j2;
   d2rule[1] = -18.0 * (j12 - j22 - 2.0 * j1j2 / x);
   d3rule[1] = -18.0 * (8.0 * j22 / x - 4.0 * j1j2 + 6.0 * j1j2 / x2);

   drule[2] = -18.0 * j1j2 * x;
   d2rule[2] = -18.0 * ((j12 - j22) * x - j1j2);
   d3rule[2] = -36.0 * (j12 + 3.0 * j22 - 2.0 * x * j1j2 + j1j2 / x);

   drule[3] = -18.0 * j1j2 * x2;
   d2rule[3] = -18.0 * x2 * (j12 - j22);
   d3rule[3] = -72.0 * x * (j12 + j22 - x * j1j2);
   }
   */

static void
_derivative_top_hat_cubic_spline_integration_rule (gdouble x, gdouble *rules)
{
  if (x < 1.0)
  {
    rules[0] = -18.0 * _3f4 (0, x, GSL_DBL_EPSILON); // limit when x-> 0: 0.0;
    rules[1] = -18.0 * _3f4 (1, x, GSL_DBL_EPSILON); // limit when x-> 0: -18.0*(-0.25);
    rules[2] = -18.0 * _3f4 (2, x, GSL_DBL_EPSILON); // limit when x-> 0: 0.0;
    rules[3] = -18.0 * _3f4 (3, x, GSL_DBL_EPSILON); // limit when x-> 0: -18.0*(1.0 - M_LN2 - M_EULER);
  }
  else
  {
    const gdouble x2 = x * x;    /* x = k * R */
    const gdouble x3 = x2 * x;
    const gdouble x4 = x3 * x;
    const gdouble x5 = x4 * x;
    const gdouble x_2 = 2.0 * x;
    const gdouble cos_x = cos(x);
    const gdouble cos_x_2 = cos_x * cos_x;
    const gdouble sin_x = sin(x);
    const gdouble sin_x_2 = sin_x * sin_x;
    const gdouble cos_2x = cos_x_2 - sin_x_2;
    const gdouble sin_2x = 2.0 * cos_x * sin_x;
    const gdouble sin_integral_2x = gsl_sf_Si (x_2);
    const gdouble cos_integral_2x = gsl_sf_Ci (x_2);

    rules[0] = -18.0 * (-18.0 -20.0 * x2 + 2.0 * (9.0 - 8.0 * x2 + x4) * cos_2x + x * (36.0 + x2) * sin_2x + 4.0 * x5 * sin_integral_2x) / (60.0 * x5);
    rules[1] = -18.0 * ((-3.0 -4.0 * x2 + (3.0 - 2.0 * x2) * cos_2x + 6.0 * x * sin_2x) / (8.0 * x4) + 0.25);
    rules[2] = -18.0 * (-1.0 -2.0 * x2 + cos_2x + 2.0 * x * sin_2x + x3 * sin_integral_2x) / (2.0 * x3);
    rules[3] = -18.0 * (log(x) - cos_integral_2x + (-3.0 * sin_x_2 - x2 * cos_x_2 + 6.0 * x * cos_x * sin_x)/ (2.0 * x2) - (1.0 - M_LN2 - M_EULER));
  }
}

static void
_top_hat_cubic_spline_integration_rule_taylor (gdouble xa, gdouble xb, gdouble R, gdouble *rules)
{
  const gdouble Rxa = R * xa;
  const gdouble delta = (xb - xa);
  const gdouble delta2 = delta * delta;
  const gdouble delta3 = delta2 * delta;
  const gdouble delta4 = delta2 * delta2;
  const gdouble Rdelta = delta * R;
  const gdouble Rdelta2 = Rdelta * Rdelta;
  const gdouble Rdelta3 = Rdelta2 * Rdelta;
  const gdouble Rdelta4 = Rdelta2 * Rdelta2;
  gdouble dnj12_x2[5];

  _top_hat_cubic_spline_integration_derivatives_rule (Rxa, dnj12_x2);

  rules[0] = delta  * (dnj12_x2[0] / 1.0 + dnj12_x2[1] * Rdelta / 2.0 + dnj12_x2[2] * Rdelta2 /  6.0 + dnj12_x2[3] * Rdelta3 / 24.0 + dnj12_x2[4] * Rdelta4 / 120.0);
  rules[1] = delta2 * (dnj12_x2[0] / 2.0 + dnj12_x2[1] * Rdelta / 3.0 + dnj12_x2[2] * Rdelta2 /  8.0 + dnj12_x2[3] * Rdelta3 / 30.0 + dnj12_x2[4] * Rdelta4 / 144.0);
  rules[2] = delta3 * (dnj12_x2[0] / 3.0 + dnj12_x2[1] * Rdelta / 4.0 + dnj12_x2[2] * Rdelta2 / 10.0 + dnj12_x2[3] * Rdelta3 / 36.0 + dnj12_x2[4] * Rdelta4 / 168.0);
  rules[3] = delta4 * (dnj12_x2[0] / 4.0 + dnj12_x2[1] * Rdelta / 5.0 + dnj12_x2[2] * Rdelta2 / 12.0 + dnj12_x2[3] * Rdelta3 / 42.0 + dnj12_x2[4] * Rdelta4 / 192.0);
}

/*
 * nc_matter_var_gaussian_cubic_spline_integration_rule:
 * @x: x = k * R FIXME
 * @rules: FIXME
 * FIXME
 *
 *
 static void
 nc_matter_var_gaussian_cubic_spline_integration_rule (gdouble x, gdouble *rules)
 {
   gdouble x2 = x * x;   // x = k * R
   gdouble erf_x_sqrt_pi = M_SQRTPI * gsl_sf_erf (x);
   gdouble exp_mx2 = exp(-x2);

   rules[0] = erf_x_sqrt_pi / 2.0;

   rules[1] = -exp_mx2 / 2.0;

   rules[2] = (-2.0 * x * exp_mx2 + erf_x_sqrt_pi) / 4.0;

   rules[3] = -exp_mx2 * (1.0 + x2)/ 2.0;
   }
   */

/*
 * nc_matter_var_derivative_gaussian_cubic_spline_integration_rule:
 * @x: x = k * R FIXME
 * @rules: FIXME
 * FIXME
 *
 *
 static void
 nc_matter_var_derivative_gaussian_cubic_spline_integration_rule (gdouble x, gdouble *rules)
 {
   gdouble x2 = x * x;    // x = k * R
   gdouble x4 = x2 * x2;
   gdouble erf_x_sqrt_pi = M_SQRTPI * gsl_sf_erf (x);
   gdouble exp_mx2 = exp(-x2);

   rules[0] = -2.0 * (-2.0 * x * exp_mx2 + erf_x_sqrt_pi) / 4.0;

   rules[1] = 2.0 * exp_mx2 * (1.0 + x2)/ 2.0;

   rules[2] = -2.0 * (-exp_mx2 * x * (3.0 + 2.0 * x2)/ 4.0 + 3.0 * erf_x_sqrt_pi / 8.0);

   rules[3] = 2.0 * exp_mx2 * (1.0 + x2 + x4/2.0);
   }
   */

/* The structure below corresponds to the spline' s structure. Valid for version 1.13 of gsl. Verify if it is the same for older and future versions. */
typedef struct
{
  gdouble * c;
  gdouble * g;
  gdouble * diag;
  gdouble * offdiag;
} cspline_state_t;

/**
 * nc_matter_var_shift_scale_rule:
 * @rules_a: are the rules computed at the point k*R = a.
 * @rules_b: are the rules computed at the point k*R = b.
 * @R: is the radius of the cluster which is related to its mass.
 * @a: is the lower point (k*R) of the spline.
 * @rule: is the vector with the four components whose sum gives the integration of the
 * spline with nodes a (lower) and b (upper).
 *
 * FIXME
 * 
 */
static void
nc_matter_var_shift_scale_rule (gdouble *rules_a, gdouble *rules_b, gdouble R, gdouble a, gdouble *rule)
{
  gdouble a2 = a * a;
  gdouble a3 = a2 * a;
  gdouble R2 = R * R;
  gdouble R3 = R2 * R;
  gdouble R4 = R3 * R;
  gdouble r0 = (rules_b[0] - rules_a[0]) / R;
  gdouble r1 = (rules_b[1] - rules_a[1]) / R2;
  gdouble r2 = (rules_b[2] - rules_a[2]) / R3;
  gdouble r3 = (rules_b[3] - rules_a[3]) / R4;

  rule[0] = r0;
  rule[1] = (r1 - a * r0);
  rule[2] = (r2 - 2.0 * a * r1 + a2 * r0);
  rule[3] = (r3 - 3.0 * a * r2 + 3.0 * a2 * r1 - a3 * r0);
}

/**
 * nc_matter_var_deriv_shift_scale_rule:
 * @rules_a: are the rules computed at the point k*R = a.
 * @rules_b: are the rules computed at the point k*R = b.
 * @R: is the radius of the cluster which is related to its mass.
 * @a: is the lower point (k*R) of the spline.
   * @rule: is the vector with the four components whose sum gives the integration of the
 * spline with nodes a (lower) and b (upper).
   *
 * FIXME
 */
static void
nc_matter_var_deriv_shift_scale_rule (gdouble *rules_a, gdouble *rules_b, gdouble R, gdouble a, gdouble *rule)
{
  gdouble a2 = a * a;
  gdouble a3 = a2 * a;
  gdouble R2 = R * R;
  gdouble R3 = R2 * R;
  gdouble R4 = R3 * R;
  gdouble R5 = R4 * R;
  gdouble r0 = (rules_b[0] - rules_a[0]) / R2;
  gdouble r1 = (rules_b[1] - rules_a[1]) / R3;
  gdouble r2 = (rules_b[2] - rules_a[2]) / R4;
  gdouble r3 = (rules_b[3] - rules_a[3]) / R5;

  rule[0] = r0;
  rule[1] = (r1 - a * r0);
  rule[2] = (r2 - 2.0 * a * r1 + a2 * r0);
  rule[3] = (r3 - 3.0 * a * r2 + 3.0 * a2 * r1 - a3 * r0);
}

/**
 * nc_matter_var_integrate_node:
 * @s: FIXME
 * @n: indicates in which node the rules will be computed.
 * @rules: are the rules computed at the point FIXME.
 *
 * FIXME
 *
 * Returns: FIXME
 */
static gdouble
nc_matter_var_integrate_node (NcmSpline *s, int n, gdouble *rules)
{
  gdouble res = 0.0;
  gdouble y_lo;

  y_lo = ncm_vector_get (s->yv, n);

  if (NCM_IS_SPLINE_CUBIC (s))
  {
    const NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
    const gdouble b_i = ncm_vector_get (sc->b, n);
    const gdouble c_i = ncm_vector_get (sc->c, n);
    const gdouble d_i = ncm_vector_get (sc->d, n);
    res += y_lo * rules[0] + b_i * rules[1] + c_i * rules[2] + d_i * rules[3];
  }
  else
    g_error ("The spline integration algorithm must be used with a cubic spline.");

  return res;
}

static gdouble
nc_matter_var_integrate_spline (NcMatterVar *vp, NcmSpline *s, void (*calc_rule)(gdouble, gdouble *), void (*shift_rule)(gdouble *, gdouble *, gdouble, gdouble, gdouble *), void (*calc_rule_taylor) (gdouble, gdouble, gdouble, gdouble *), gdouble R)
{
  gdouble rules_a[4];
  gdouble rules_b[4];
  gdouble rules[4];
  gdouble *ra = rules_a;
  gdouble *rb = rules_b;
  gdouble *temp;
  gboolean taylor = FALSE;
  guint i;
  gdouble res = 0.0;
  gdouble partial;

  NCM_UNUSED (vp);

  calc_rule(ncm_vector_get (s->xv, 0) * R, ra);
  for(i = 0; i < s->len - 1; i++)
  {
    gdouble delta = (ncm_vector_get (s->xv, i + 1) - ncm_vector_get (s->xv, i)) * R;
    if (delta < 1e-2 && (calc_rule_taylor != NULL))
    {
      calc_rule_taylor (ncm_vector_get (s->xv, i), ncm_vector_get (s->xv, i + 1), R, rules);
      taylor = TRUE;
    }
    else
    {
      if (taylor)
        calc_rule (ncm_vector_get (s->xv, i) * R, ra);
      taylor = FALSE;
      calc_rule (ncm_vector_get (s->xv, i + 1) * R, rb);
      shift_rule (ra, rb, R, ncm_vector_get (s->xv, i), rules);
    }

    partial = nc_matter_var_integrate_node (s, i, rules);
    res += partial;

    temp = ra;
    ra = rb;
    rb = temp;
  }

  return res;
}

/******************* Spectral Moment ************************************************************************************/

/**
 * nc_matter_var_spectral_moment_over_growth2:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @n: FIXME
 *
 * FIXME
 * \frac{\sigma^2}{D^2}
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_spectral_moment_over_growth2 (NcMatterVar *vp, NcHICosmo *model, gint n)
{
  if (NC_IS_WINDOW_TOPHAT (vp->wp))
  {
    return nc_matter_var_spectral_moment_over_growth2_tophat (vp, model, n);
  }
  else if (NC_IS_WINDOW_GAUSSIAN (vp->wp))
  {
    return nc_matter_var_spectral_moment_over_growth2_gaussian (vp, model, n);
  }
  else
    g_assert_not_reached ();

  return 0.0;
}

/**
 * nc_matter_var_spectral_moment_over_growth2_gaussian:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_spectral_moment_over_growth2_gaussian (NcMatterVar *vp, NcHICosmo *model, gint n)     /* \frac{\sigma_j^2}{D^2} */
{
  gdouble sm_over_growth2, error;
  static gsl_integration_workspace *w = NULL; //gsl_integration_workspace_alloc (INT_PARTITION);
  gsl_function F;
  int_temp_stc ts;

  ts.j = n;
  ts.vp = vp;
  ts.model = model;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

  F.function = &_nc_matter_var_integrand_gaussian;
  F.params = &ts;
  gsl_integration_qagiu (&F, 0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &sm_over_growth2, &error);

  // gsl_integration_workspace_free (w);
  return sm_over_growth2;
}

/**
 * nc_matter_var_spectral_moment_over_growth2_tophat:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_spectral_moment_over_growth2_tophat (NcMatterVar *vp, NcHICosmo *model, gint n)     /* \frac{\sigma_j^2}{D^2} */
{
  gdouble sm_over_growth2_ksmall = 0.0, sm_over_growth2_one = 0.0, sm_over_growth2_osc = 0.0;
  gdouble error, error_1, error_cos, error_sin, step, pres, result = 0.0;
  static gsl_integration_workspace *w = NULL; //gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);
  static gsl_integration_qawo_table *wosc_cos = NULL;
  static gsl_integration_qawo_table *wosc_sin = NULL;
  gsl_function F, F_cos, F_sin, Ft;
  int_temp_stc ts;

  ts.j = n;
  ts.vp = vp;
  ts.model = model;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);
  if (wosc_cos == NULL)
    wosc_cos = gsl_integration_qawo_table_alloc (2.0, _NC_LSTEP2, GSL_INTEG_COSINE, NCM_INTEGRAL_PARTITION);
  if (wosc_sin == NULL)
    wosc_sin = gsl_integration_qawo_table_alloc (2.0, _NC_LSTEP2, GSL_INTEG_SINE, NCM_INTEGRAL_PARTITION);

  F.function = &_nc_matter_var_integrand_tophat_one;
  F_cos.function = &_nc_matter_var_integrand_tophat_cosine;
  F_sin.function = &_nc_matter_var_integrand_tophat_sine;
  Ft.function = &_nc_matter_var_integrand_tophat_ksmall;

  F.params = &ts;
  F_cos.params = &ts;
  F_sin.params = &ts;
  Ft.params = &ts;


  gsl_integration_qag (&Ft, 0.0, 1.0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &sm_over_growth2_ksmall, &error);   /* integral completa de 0 a 1 */
  result += sm_over_growth2_ksmall;

  gsl_integration_qag (&F, 1.0, 1000.0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &sm_over_growth2_one, &error_1);
  result += sm_over_growth2_one;

  for (step = 1.0; 1; step += _NC_LSTEP2)
  {
    gdouble pres_t = 0.0;
    gsl_integration_qawo (&F_cos, step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, wosc_cos, &pres, &error_cos);
    pres_t = pres;
    gsl_integration_qawo (&F_sin, step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, wosc_sin, &pres, &error_sin);
    pres_t += pres;

    if (fabs(pres_t/sm_over_growth2_osc) < NCM_DEFAULT_PRECISION * 1e-3)
      break;
    sm_over_growth2_osc += pres_t;
    if (pres_t == 0.0)
      break;
  }

  result += sm_over_growth2_osc;
  return result;
}

/* Fizemos ums troca de variavel k = exp(l) - 1, para facilitar a convergencia, jah que para
alguns valores de R a convergencia eh lenta.
*/

/******************* Variance derivative with respect to R ****************************************************************/

static gdouble
nc_matter_var_deriv_variance_over_growth2_integrand_gaussian (gdouble kR, gpointer params)     /* (1/D^2)*d(integrando de \sigma^2)/dR, integrando com W * dW/dR */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble W = nc_window_eval_fourier (vp->wp, k, ts->R);
  gdouble dW = nc_window_deriv_fourier (vp->wp, k, ts->R);
  gdouble k2 = k*k;
  gdouble dv_integrand = matter_P * k2 * 2.0 * W * dW / (c1 * ts->R);

  return dv_integrand;
}

static gdouble
_nc_matter_var_dvariance_over_growth2_dR_gaussian (NcMatterVar *vp, NcHICosmo *model, gdouble R)
{
  static gsl_integration_workspace *w = NULL;
  gsl_function F;
  gdouble dv_D2dR, error;
  int_temp_stc ts;

  ts.j = 0;
  ts.R = R;
  ts.vp = vp;
  ts.model = model;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

  F.function = &nc_matter_var_deriv_variance_over_growth2_integrand_gaussian;
  F.params = &ts;

  gsl_integration_qagiu (&F, 0, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &dv_D2dR, &error);

  return dv_D2dR;
}

static gdouble
_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_ksmall (gdouble kR, gpointer params)     /* (1/D^2)*d(integrando de \sigma^2)/dR, integrando com W * dW/dR */
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble W = nc_window_eval_fourier (vp->wp, k, ts->R);
  gdouble dW = nc_window_deriv_fourier (vp->wp, k, ts->R);
  gdouble k2 = k*k;
  gdouble dv_integrand = matter_P * k2 * 2.0 * W * dW / (c1 * ts->R);

  return dv_integrand;
}

static gdouble
_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_one (gdouble kR, gpointer params)
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble k2 = k*k;
  gdouble dv_integrand = - matter_P * (18.0 + 27.0/(k2 * ts->R * ts->R)) / (c1 * k2 * gsl_pow_5(ts->R) * ts->R);

  return dv_integrand;
}

static gdouble
_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_cosine (gdouble kR, gpointer params)
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble k2 = k*k;
  gdouble dv_integrand = matter_P * (27.0/(k2 * ts->R * ts->R) - 36.0) / (c1 * k2 * gsl_pow_5(ts->R) * ts->R);

  return dv_integrand;
}

static gdouble
_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_sine (gdouble kR, gpointer params)
{
  int_temp_stc *ts = (int_temp_stc *)params;
  NcMatterVar *vp = ts->vp;
  gdouble k = kR / ts->R;
  gdouble c1 = 2.0 * M_PI * M_PI;
  gdouble matter_P = nc_transfer_func_matter_powerspectrum (vp->tf, ts->model, k);
  gdouble k2 = k*k;
  gdouble dv_integrand = matter_P * (54.0/(k2 * ts->R * ts->R) - 9.0) / (c1 * k * gsl_pow_4(ts->R) * ts->R);

  return dv_integrand;
}

#define _NC_LSTEP2 5.0

static gdouble
_nc_matter_var_dvariance_over_growth2_dR_tophat (NcMatterVar *vp, NcHICosmo *model, gdouble R)     /* \frac{d\sigma^2}{D^2dR} */
{
  gdouble dv_over_growth2_ksmall = 0.0, dv_over_growth2_one = 0.0, dv_over_growth2_osc = 0.0;
  gdouble error, error_1, error_cos, error_sin, step, pres, result = 0.0;
  static gsl_integration_workspace *w = NULL;
  static gsl_integration_qawo_table *wosc_cos = NULL;
  static gsl_integration_qawo_table *wosc_sin = NULL;
  gsl_function F, F_cos, F_sin, Ft;
  int_temp_stc ts;

  ts.j = 0;
  ts.R = R;
  ts.vp = vp;
  ts.model = model;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);
  if (wosc_cos == NULL)
    wosc_cos = gsl_integration_qawo_table_alloc (2.0, _NC_LSTEP2, GSL_INTEG_COSINE, NCM_INTEGRAL_PARTITION);
  if (wosc_sin == NULL)
    wosc_sin = gsl_integration_qawo_table_alloc (2.0, _NC_LSTEP2, GSL_INTEG_SINE, NCM_INTEGRAL_PARTITION);

  F.function = &_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_one;
  F_cos.function = &_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_cosine;
  F_sin.function = &_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_sine;
  Ft.function = &_nc_matter_var_deriv_variance_over_growth2_integrand_tophat_ksmall;

  F.params = &ts;
  F_cos.params = &ts;
  F_sin.params = &ts;
  Ft.params = &ts;

  //  for (i = 0; i < var_spline->size-1; i++)
  //  {
  //    gsl_integration_qag (&Ft, var_spline->x[i] * vp->wp->R, var_spline->x[i+1] * vp->wp->R, 0, 1e-9, INT_PARTITION, 6, w, &dv_over_growth2_ksmall, &error);   /* integral completa de 0 a 1 */

  gsl_integration_qag (&Ft, 1.0e-4, 1.0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &dv_over_growth2_ksmall, &error);

  result += dv_over_growth2_ksmall;
  //  }

  gsl_integration_qag (&F, 1.0, 1000.0, 0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &dv_over_growth2_one, &error_1);
  result += dv_over_growth2_one;

  for (step = 1.0; 1; step += _NC_LSTEP2)
    //  for (step = 1.0; step < 11.0; step++)
  {
    gdouble pres_t = 0.0;
    gsl_integration_qawo (&F_cos, step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, wosc_cos, &pres, &error_cos);
    pres_t = pres;
    gsl_integration_qawo (&F_sin, step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, wosc_sin, &pres, &error_sin);
    pres_t += pres;
    if (fabs(pres_t/dv_over_growth2_osc) < NCM_DEFAULT_PRECISION * 1e-3)
      break;
    dv_over_growth2_osc += pres_t;
    if (pres_t == 0.0)
      break;
  }

  result += dv_over_growth2_osc;

  return result;
}

/**
 * nc_matter_var_dsigma_over_growth_dR:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnR: logarithm base e of the radius.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_dsigma0_dR (NcMatterVar *vp, NcHICosmo *model, gdouble lnR)
{
  gdouble dlnvar0_dR = nc_matter_var_dlnvar0_dR (vp, model, lnR);
  gdouble var0 = nc_matter_var_var0 (vp, model, lnR);
  gdouble sigma0 = sqrt(var0);
  gdouble dsigma0_dR = sigma0 * dlnvar0_dR / 2.0;

  return dsigma0_dR;
}

/****************** Variance normalization with sigma8 ******************************************************************/

/**
 * nc_matter_var_sigma_norma:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_matter_var_sigma8_sqrtvar0 (NcMatterVar *vp, NcHICosmo *model)    /* eh calculada usando a integracao para funcoes oscilatorias */
{
  gdouble sigma0;
  gdouble sigma0_norma;
  const gdouble ln8 = log (8.0);

  if (!NC_IS_WINDOW_TOPHAT (vp->wp))
  {
    NcWindow *temp = vp->wp;
    vp->wp = nc_window_tophat_new (); /* Normalizacao sempre com tophat */
    nc_matter_var_prepare (vp, model);
    sigma0 = sqrt(nc_matter_var_var0 (vp, model, ln8));
    nc_window_free (vp->wp);
    vp->wp = temp;
    nc_matter_var_prepare (vp, model);
  }
  else
    sigma0 = sqrt(nc_matter_var_var0 (vp, model, ln8));

  sigma0_norma = nc_hicosmo_sigma_8 (model) / sigma0;

  return sigma0_norma;
}
