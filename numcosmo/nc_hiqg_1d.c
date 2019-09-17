/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hiqg_1d.c
 *
 *  Thu February 15 14:44:56 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiqg_1d.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hiqg_1d
 * @title: NcHIQG1D
 * @short_description: Minisuperspace 1D quantum gravity models
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hiqg_1d.h"
#include "math/ncm_matrix.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_integral1d.h"
#include "math/ncm_integral1d_ptr.h"
#include "math/ncm_util.h"
#include "math/ncm_spline_rbf.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_lapack.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_poly.h>

#ifdef HAVE_ACB_H
#include <acb.h>
#include <acb_hypgeom.h>
#endif /* HAVE_ACB_H  */

#include <nvector/nvector_serial.h>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <arkode/arkode_bandpre.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIQG1DPrivate
{
  gdouble lambda;
  gdouble acs_a;
  gdouble basis_a;
  gdouble nu;
  gdouble mu;
  gdouble abstol;
  gdouble reltol;
  guint nknots;
  NcmVector *fknots;
  NcmVector *knots;
  NcmMatrix *fpsi;
  NcmMatrix *psi;
  NcmVector *fRePsi;
  NcmVector *fImPsi;
  NcmVector *RePsi;
  NcmVector *ImPsi;
  NcmVector *frho;
  NcmSpline *RePsi_s;
  NcmSpline *ImPsi_s;
  NcmSpline *rho_s;
  NcmMatrix *C0;
  NcmVector *ReC0;
  NcmVector *ImC0;
  NcmVector *A0;
  NcmMatrix *C;
  NcmVector *ReC;
  NcmVector *ImC;
  gboolean up_splines;
  gboolean noboundary;
  SUNLinearSolver LS;
  SUNMatrix A;
  gdouble h;
  gdouble ti;
  gdouble xi;
  gdouble xf;
  gpointer bohm;
  GArray *work;
  GArray *ipiv;
  GArray *iwork;
  NcmMatrix *IM;
  NcmMatrix *IM_fact;
  NcmMatrix *KM;
  NcmMatrix *JM;
  NcmMatrix *err_norm;
  NcmMatrix *err_comp;
  NcmVector *scaling;
  NcmVector *Escaling;
  NcmVector *berr;
  NcmVector *wr;
  NcmVector *wi;
  NcmMatrix *vl;
  NcmMatrix *vr;
  NcmVector *vlvr;
  gint nBohm;
  N_Vector yBohm;
 };

enum
{
  PROP_0,
  PROP_LAMBDA,
  PROP_ABSTOL,
  PROP_RELTOL,
  PROP_NKNOTS,
  PROP_NOBOUNDARY,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHIQG1D, nc_hiqg_1d, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcHIQG1DGauss, nc_hiqg_1d_gauss, nc_hiqg_1d_gauss_dup, nc_hiqg_1d_gauss_free);
G_DEFINE_BOXED_TYPE (NcHIQG1DExp,   nc_hiqg_1d_exp,   nc_hiqg_1d_exp_dup,   nc_hiqg_1d_exp_free);

static void
nc_hiqg_1d_init (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv = nc_hiqg_1d_get_instance_private (qg1d);

  self->lambda     = 0.0;
  self->acs_a      = 0.0;
  self->basis_a    = 0.0;
  self->nu         = 0.0;
  self->mu         = 0.0;

  self->abstol     = 0.0;
  self->reltol     = 0.0;

  /* PDE */
  self->nknots     = 0;
  self->fknots     = NULL;
  self->knots      = NULL;
  self->psi        = NULL;
  self->fpsi       = NULL;
  self->fRePsi     = NULL;
  self->fImPsi     = NULL;
  self->RePsi      = NULL;
  self->ImPsi      = NULL;
  self->frho       = NULL;
  self->RePsi_s    = NULL;
  self->ImPsi_s    = NULL;
  self->rho_s      = NULL;
  
  self->C0         = NULL;
  self->ReC0       = NULL;
  self->ImC0       = NULL;
  self->A0         = NULL;
  self->C          = NULL;
  self->ReC        = NULL;
  self->ImC        = NULL;

  self->up_splines = FALSE;
  self->noboundary = FALSE;
  self->LS         = NULL;
  self->A          = NULL;
  self->h          = 0.0;
  self->ti         = 0.0;
  self->xi         = 0.0;
  self->xf         = 0.0;
  
  self->bohm       = NULL;

  self->work       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->ipiv       = g_array_new (FALSE, FALSE, sizeof (gint));
  self->iwork      = g_array_new (FALSE, FALSE, sizeof (gint));

  self->IM         = NULL;
  self->IM_fact    = NULL;
  self->KM         = NULL;
  self->JM         = NULL;
  self->err_norm   = NULL;
  self->err_comp   = NULL;
  self->scaling    = NULL;
  self->Escaling   = NULL;
  self->berr       = NULL;
  self->wr         = NULL;
  self->wi         = NULL;
  self->vl         = NULL;
  self->vr         = NULL;
  self->vlvr       = NULL;

  self->nBohm      = 0;
  self->yBohm      = NULL;
}

#define _XI(i,n)   (i)
#define _LNRI(i,n) ((n) + (i) * 2 + 0)
#define _SI(i,n)   ((n) + (i) * 2 + 1)

#define _XI_STRIDE   (1)
#define _LNRI_STRIDE (2)
#define _SI_STRIDE   (2)

static void
_nc_hiqg_1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIQG1D *qg1d = NC_HIQG_1D (object);
  NcHIQG1DPrivate * const self = qg1d->priv;
  g_return_if_fail (NC_IS_HIQG_1D (object));

  switch (prop_id)
  {
    case PROP_LAMBDA:
      self->lambda  = g_value_get_double (value);
      self->acs_a   = 1.0 + 2.0 * self->lambda + 2.0 * sqrt (self->lambda * (self->lambda - 1.0));
      self->mu      = 0.25 * (self->acs_a - 1.0);
      self->basis_a = 0.5 * ncm_util_sqrt1px_m1 (4.0 * self->lambda);
      self->nu      = sqrt (self->lambda + 0.25);
      break;
    case PROP_ABSTOL:
      self->abstol = g_value_get_double (value);
      break;
    case PROP_RELTOL:
      self->reltol = g_value_get_double (value);
      break;
    case PROP_NKNOTS:
      nc_hiqg_1d_set_nknots (qg1d, g_value_get_uint (value));
      break;
    case PROP_NOBOUNDARY:
      self->noboundary = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hiqg_1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIQG1D *qg1d = NC_HIQG_1D (object);
  NcHIQG1DPrivate * const self = qg1d->priv;
  g_return_if_fail (NC_IS_HIQG_1D (object));

  switch (prop_id)
  {
    case PROP_LAMBDA:
      g_value_set_double (value, self->lambda);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, self->abstol);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_NKNOTS:
      g_value_set_uint (value, nc_hiqg_1d_get_nknots (qg1d));
      break;
    case PROP_NOBOUNDARY:
      g_value_set_boolean (value, self->noboundary);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hiqg_1d_dispose (GObject *object)
{
  NcHIQG1D *qg1d = NC_HIQG_1D (object);
  NcHIQG1DPrivate * const self = qg1d->priv;

  g_clear_pointer (&self->LS, SUNLinSolFree);
  g_clear_pointer (&self->A, SUNMatDestroy);

  ncm_vector_clear (&self->fknots);
  ncm_vector_clear (&self->knots);
  ncm_matrix_clear (&self->fpsi);
  ncm_matrix_clear (&self->psi);
  ncm_vector_clear (&self->fRePsi);
  ncm_vector_clear (&self->fImPsi);
  ncm_vector_clear (&self->RePsi);
  ncm_vector_clear (&self->ImPsi);
  ncm_vector_clear (&self->frho);
  ncm_spline_clear (&self->RePsi_s);
  ncm_spline_clear (&self->ImPsi_s);
  ncm_spline_clear (&self->rho_s);

  ncm_matrix_clear (&self->C0);
  ncm_vector_clear (&self->ReC0);
  ncm_vector_clear (&self->ImC0);
  ncm_vector_clear (&self->A0);
  ncm_matrix_clear (&self->C);
  ncm_vector_clear (&self->ReC);
  ncm_vector_clear (&self->ImC);

  g_clear_pointer (&self->work, g_array_unref);
  g_clear_pointer (&self->ipiv, g_array_unref);
  g_clear_pointer (&self->iwork, g_array_unref);

  ncm_matrix_clear (&self->IM);
  ncm_matrix_clear (&self->IM_fact);
  ncm_matrix_clear (&self->KM);
  ncm_matrix_clear (&self->JM);
  ncm_matrix_clear (&self->err_norm);
  ncm_matrix_clear (&self->err_comp);
  ncm_vector_clear (&self->scaling);
  ncm_vector_clear (&self->Escaling);
  ncm_vector_clear (&self->berr);
  ncm_vector_clear (&self->wr);
  ncm_vector_clear (&self->wi);
  ncm_matrix_clear (&self->vl);
  ncm_matrix_clear (&self->vr);
  ncm_vector_clear (&self->vlvr);

  g_clear_pointer (&self->yBohm, N_VDestroy);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiqg_1d_parent_class)->dispose (object);
}

static void
_nc_hiqg_1d_finalize (GObject *object)
{
  NcHIQG1D *qg1d = NC_HIQG_1D (object);
  NcHIQG1DPrivate * const self = qg1d->priv;

  if (self->bohm != NULL)
  {
    ARKStepFree (&self->bohm);
    self->bohm = NULL;
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiqg_1d_parent_class)->finalize (object);
}

static void
nc_hiqg_1d_class_init (NcHIQG1DClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hiqg_1d_set_property;
  object_class->get_property = &_nc_hiqg_1d_get_property;
  object_class->dispose      = &_nc_hiqg_1d_dispose;
  object_class->finalize     = &_nc_hiqg_1d_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_LAMBDA,
                                   g_param_spec_double ("lambda",
                                                        NULL,
                                                        "\\lambda",
                                                        0.0, G_MAXDOUBLE, 0.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "absolute tolerance",
                                                        0.0, G_MAXDOUBLE, 1.0e-50,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "relative tolerance",
                                                        0.0, 1.0, 1.0e-9,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NKNOTS,
                                   g_param_spec_uint ("nknots",
                                                      NULL,
                                                      "n_k",
                                                      6, G_MAXUINT, 50,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NOBOUNDARY,
                                   g_param_spec_boolean ("noboundary",
                                                         NULL,
                                                         "no boundary condition at x_f",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
}

/**
 * nc_hiqg_1d_gauss_new:
 * @mean: gaussian mean
 * @alpha: power-law 
 * @sigma: Gaussian width 
 * @Hi: linear imaginary phase
 * 
 * Creates a new Gaussian wave function.
 * 
 * Returns: (transfer full): a new #NcHIQG1DGauss
 */
NcHIQG1DGauss *
nc_hiqg_1d_gauss_new (const gdouble mean, const gdouble alpha, const gdouble sigma, const gdouble Hi)
{
  NcHIQG1DGauss *qm_gauss = g_new (NcHIQG1DGauss, 1);

  qm_gauss->mean  = mean;
  qm_gauss->alpha = alpha;
  qm_gauss->sigma = sigma;
  qm_gauss->Hi    = Hi;

  g_assert_cmpfloat (alpha, >, -0.5);
  
  {
    gint signp             = 0;
    const gdouble f1       = (alpha - 0.5) * M_LN2 + (2.0 * alpha + 1.0) * log (sigma);
    const gdouble ap12     = (alpha + 0.5);
    const gdouble lng_ap12 = lgamma_r (ap12, &signp);
    if (mean != 0.0)
    {
      const gdouble ap1     = (alpha + 1.0);
      const gdouble arg     = mean / (M_SQRT2 * sigma);
      const gdouble arg2    = arg * arg;
      const gdouble lng_ap1 = lgamma_r (ap1,  &signp);
      const gdouble t1      = log (2.0 * fabs (arg)) + lng_ap1  + log (gsl_sf_hyperg_1F1 (0.5 - alpha, 1.5, - arg2));
      const gdouble t2      =                          lng_ap12 + log (gsl_sf_hyperg_1F1 (    - alpha, 0.5, - arg2));

      qm_gauss->lnNorm = log (exp (t2) + GSL_SIGN (arg) * exp (t1)) + f1;
    }
    else
    {
      qm_gauss->lnNorm = lng_ap12 + f1;
    }
  }
  
  return qm_gauss;
}

/**
 * nc_hiqg_1d_gauss_dup:
 * @qm_gauss: a #NcHIQG1DGauss
 * 
 * Duplicates @qm_gauss.
 * 
 * Returns: (transfer full): a duplicate of @qm_gauss.
 */
NcHIQG1DGauss *
nc_hiqg_1d_gauss_dup (NcHIQG1DGauss *qm_gauss)
{
  NcHIQG1DGauss *qm_gauss_dup = g_new (NcHIQG1DGauss, 1);
  qm_gauss_dup[0] = qm_gauss[0];

  return qm_gauss_dup;
}

/**
 * nc_hiqg_1d_gauss_free:
 * @qm_gauss: a #NcHIQG1DGauss
 * 
 * Frees @qm_gauss.
 * 
 */
void 
nc_hiqg_1d_gauss_free (NcHIQG1DGauss *qm_gauss)
{
  g_free (qm_gauss);
}

/**
 * nc_hiqg_1d_gauss_eval:
 * @qm_gauss: a #NcHIQG1DGauss
 * @x: the point where to evaluate $\psi(x)$
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi$
 * 
 * Evaluates @qm_gauss at @x.
 * 
 */
void 
nc_hiqg_1d_gauss_eval (NcHIQG1DGauss *qm_gauss, const gdouble x, gdouble *psi)
{  
  const gdouble lnx     = log (x);
  const gdouble xmean   = x - qm_gauss->mean;
  const gdouble xmean2  = xmean * xmean;
  const gdouble sigma2  = qm_gauss->sigma * qm_gauss->sigma;
  complex double psi0c  = cexp (- 0.5 * qm_gauss->lnNorm + qm_gauss->alpha * lnx - 0.25 * xmean2 / sigma2 + 0.5 * xmean2 * I * qm_gauss->Hi);
  //complex double psi0c  = cexp (- 0.5 * qm_gauss->lnNorm + qm_gauss->alpha * lnx - 0.25 * xmean2 / sigma2 + xmean * I * qm_gauss->Hi);
  
  psi[0] = creal (psi0c);
  psi[1] = cimag (psi0c);
}

/**
 * nc_hiqg_1d_gauss_eval_hermit:
 * @qm_gauss: a #NcHIQG1DGauss
 * @x: the point where to evaluate $\psi(x)$
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi$
 * 
 * Evaluates @qm_gauss at @x removing the Hemite weight.
 * 
 */
void 
nc_hiqg_1d_gauss_eval_hermit (NcHIQG1DGauss *qm_gauss, const gdouble x, gdouble *psi)
{  
  const gdouble xmean   = x - qm_gauss->mean;
  const gdouble xmean2  = xmean * xmean;
  complex double psi0c  = cexp (- 0.5 * qm_gauss->lnNorm + 0.5 * xmean2 * I * qm_gauss->Hi);

  psi[0] = creal (psi0c);
  psi[1] = cimag (psi0c);
}

/**
 * nc_hiqg_1d_gauss_eval_lnRS:
 * @qm_gauss: a #NcHIQG1DGauss
 * @x: the point where to evaluate $\psi(x)$
 * @lnRS: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\ln(R)$ and $S$ in $\psi = e^{\ln(R) + iS}$
 * 
 * Evaluates @qm_gauss at @x.
 * 
 */
void 
nc_hiqg_1d_gauss_eval_lnRS (NcHIQG1DGauss *qm_gauss, const gdouble x, gdouble *lnRS)
{  
  const gdouble lnx     = log (x);
  const gdouble xmean   = x - qm_gauss->mean;
  const gdouble xmean2  = xmean * xmean;
  const gdouble sigma2  = qm_gauss->sigma * qm_gauss->sigma;
  const gdouble lnR     = - 0.5 * qm_gauss->lnNorm + qm_gauss->alpha * lnx - 0.25 * xmean2 / sigma2;
  const gdouble S       = 0.5 * x * x * qm_gauss->Hi;
  //const gdouble S       = x * qm_gauss->Hi;

  lnRS[0] = lnR;
  lnRS[1] = S;
}

/**
 * nc_hiqg_1d_exp_new:
 * @n: power-law 
 * @V: Volume 
 * @pV: Volume momentum
 * 
 * Creates a new Exponential wave function.
 * 
 * Returns: (transfer full): a new #NcHIQG1DExp
 */
NcHIQG1DExp *
nc_hiqg_1d_exp_new (const gdouble n, const gdouble V, const gdouble pV)
{
  NcHIQG1DExp *qm_exp = g_new (NcHIQG1DExp, 1);

  qm_exp->n  = n;
  qm_exp->V  = V;
  qm_exp->pV = pV;

  g_assert_cmpfloat (n, >, 2.0);
  g_assert_cmpfloat (V, >, 0.0);
  
  {
    gint signp     = 0;
    qm_exp->lnNorm = - n * log (n - 1.0) + lgamma_r (n, &signp) + log (V);
  }
  
  return qm_exp;
}

/**
 * nc_hiqg_1d_exp_dup:
 * @qm_exp: a #NcHIQG1DExp
 * 
 * Duplicates @qm_exp.
 * 
 * Returns: (transfer full): a duplicate of @qm_exp.
 */
NcHIQG1DExp *
nc_hiqg_1d_exp_dup (NcHIQG1DExp *qm_exp)
{
  NcHIQG1DExp *qm_exp_dup = g_new (NcHIQG1DExp, 1);
  qm_exp_dup[0] = qm_exp[0];

  return qm_exp_dup;
}

/**
 * nc_hiqg_1d_exp_free:
 * @qm_exp: a #NcHIQG1DExp
 * 
 * Frees @qm_exp.
 * 
 */
void 
nc_hiqg_1d_exp_free (NcHIQG1DExp *qm_exp)
{
  g_free (qm_exp);
}

/**
 * nc_hiqg_1d_exp_eval:
 * @qm_exp: a #NcHIQG1DExp
 * @x: the point where to evaluate $\psi(x)$
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $psi$
 * 
 * Evaluates @qm_exp at @x.
 * 
 */
void 
nc_hiqg_1d_exp_eval (NcHIQG1DExp *qm_exp, const gdouble x, gdouble *psi)
{  
  const gdouble xV      = x / qm_exp->V;
  const gdouble lnxV    = log (xV);
  complex double psi0c  = cexp (-0.5 * qm_exp->lnNorm + 0.5 * (qm_exp->n - 1.0) * (lnxV - xV) + I * qm_exp->pV * x);
  
  psi[0] = creal (psi0c);
  psi[1] = cimag (psi0c);
}

/**
 * nc_hiqg_1d_exp_eval_lnRS:
 * @qm_exp: a #NcHIQG1DExp
 * @x: the point where to evaluate $\psi(x)$
 * @lnRS: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\ln(R)$ and $S$ in $\psi = e^{\ln(R) + iS}$
 * 
 * Evaluates @qm_exp at @x.
 * 
 */
void 
nc_hiqg_1d_exp_eval_lnRS (NcHIQG1DExp *qm_exp, const gdouble x, gdouble *lnRS)
{  
  const gdouble xV   = x / qm_exp->V;
  const gdouble lnxV = log (xV);
  const gdouble lnR  = -0.5 * qm_exp->lnNorm + 0.5 * (qm_exp->n - 1.0) * (lnxV - xV);
  const gdouble S    = qm_exp->pV * x;

  lnRS[0] = lnR;
  lnRS[1] = S;
}

/**
 * nc_hiqg_1d_new:
 * 
 * Creates a new #NcHIQG1D object.
 * 
 * Returns: (transfer full): a new #NcHIQG1D.
 */
NcHIQG1D *
nc_hiqg_1d_new (void)
{
  NcHIQG1D *qg1d = g_object_new (NC_TYPE_HIQG_1D,
                                 NULL);
  return qg1d;
}

/**
 * nc_hiqg_1d_new_full:
 * @nknots: number of knots
 * @lambda: $\lambda$
 * 
 * Creates a new #NcHIQG1D object.
 * 
 * Returns: (transfer full): a new #NcHIQG1D.
 */
NcHIQG1D *
nc_hiqg_1d_new_full (guint nknots, gdouble lambda)
{
  NcHIQG1D *qg1d = g_object_new (NC_TYPE_HIQG_1D,
                                 "nknots", nknots,
                                 "lambda", lambda,
                                  NULL);
  return qg1d;
}

/**
 * nc_hiqg_1d_ref:
 * @qg1d: a #NcHIQG1D
 *
 * Increase the reference of @qg1d by one.
 *
 * Returns: (transfer full): @qg1d.
 */
NcHIQG1D *
nc_hiqg_1d_ref (NcHIQG1D *qg1d)
{
  return g_object_ref (qg1d);
}

/**
 * nc_hiqg_1d_free:
 * @qg1d: a #NcHIQG1D
 *
 * Decrease the reference count of @qg1d by one.
 *
 */
void
nc_hiqg_1d_free (NcHIQG1D *qg1d)
{
  g_object_unref (qg1d);
}

/**
 * nc_hiqg_1d_clear:
 * @qg1d: a #NcHIQG1D
 *
 * Decrease the reference count of @qg1d by one, and sets the pointer *qg1d to
 * NULL.
 *
 */
void
nc_hiqg_1d_clear (NcHIQG1D **qg1d)
{
  g_clear_object (qg1d);
}

/**
 * nc_hiqg_1d_set_nknots:
 * @qg1d: a #NcHIQG1D
 * @nknots: number of knots
 *
 * Sets the initial number of knots to be used in the wave function mesh.
 * 
 */
void 
nc_hiqg_1d_set_nknots (NcHIQG1D *qg1d, const guint nknots)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  if (self->nknots != nknots)
  {
    if (self->nknots > 0)
    {
      ncm_vector_clear (&self->fknots);
      ncm_vector_clear (&self->knots);
      ncm_matrix_clear (&self->fpsi);
      ncm_matrix_clear (&self->psi);
      ncm_matrix_clear (&self->C0);
      ncm_vector_clear (&self->ReC0);
      ncm_vector_clear (&self->ImC0);
      ncm_vector_clear (&self->A0);
      ncm_matrix_clear (&self->C);
      ncm_vector_clear (&self->ReC);
      ncm_vector_clear (&self->ImC);
      
      ncm_matrix_clear (&self->IM);
      ncm_matrix_clear (&self->IM_fact);
      ncm_matrix_clear (&self->KM);
      ncm_matrix_clear (&self->JM);
      ncm_matrix_clear (&self->err_norm);
      ncm_matrix_clear (&self->err_comp);
      ncm_vector_clear (&self->scaling);
      ncm_vector_clear (&self->Escaling);
      ncm_vector_clear (&self->berr);
      ncm_vector_clear (&self->wr);
      ncm_vector_clear (&self->wi);
      ncm_matrix_clear (&self->vl);
      ncm_matrix_clear (&self->vr);
      ncm_vector_clear (&self->vlvr);
    }

    self->nknots = nknots;
    if (self->nknots > 0)
    {
      self->fknots   = ncm_vector_new (self->nknots + 1);
      self->knots    = ncm_vector_get_subvector_stride (self->fknots, 1, self->nknots, 1);
      self->fpsi     = ncm_matrix_new (2, self->nknots + 1);
      self->psi      = ncm_matrix_get_submatrix (self->fpsi, 0, 1, 2, self->nknots);
      self->fRePsi   = ncm_matrix_get_row (self->fpsi, 0);
      self->fImPsi   = ncm_matrix_get_row (self->fpsi, 1);
      self->RePsi    = ncm_matrix_get_row (self->psi, 0);
      self->ImPsi    = ncm_matrix_get_row (self->psi, 1);
      self->frho     = ncm_vector_new (self->nknots + 1);
      self->RePsi_s  = ncm_spline_cubic_notaknot_new_full (self->fknots, self->fRePsi, FALSE); 
      self->ImPsi_s  = ncm_spline_cubic_notaknot_new_full (self->fknots, self->fImPsi, FALSE); 
      self->rho_s    = ncm_spline_cubic_notaknot_new_full (self->fknots, self->frho, FALSE);       
      self->C0       = ncm_matrix_new (2, self->nknots);
      self->A0       = ncm_vector_new (self->nknots * 2);
      self->ReC0     = ncm_matrix_get_row (self->C0, 0);
      self->ImC0     = ncm_matrix_get_row (self->C0, 1);
      self->C        = ncm_matrix_new (2, self->nknots);
      self->ReC      = ncm_matrix_get_row (self->C, 0);
      self->ImC      = ncm_matrix_get_row (self->C, 1);
      self->IM       = ncm_matrix_new (self->nknots, self->nknots);
      self->IM_fact  = ncm_matrix_new (self->nknots, self->nknots);
      self->KM       = ncm_matrix_new (self->nknots, self->nknots);
      self->JM       = ncm_matrix_new (self->nknots, self->nknots);
      self->err_norm = ncm_matrix_new (3, self->nknots);
      self->err_comp = ncm_matrix_new (3, self->nknots);
      self->scaling  = ncm_vector_new (self->nknots);
      self->Escaling = ncm_vector_new (self->nknots);
      self->berr     = ncm_vector_new (self->nknots);
      self->wr       = ncm_vector_new (self->nknots);
      self->wi       = ncm_vector_new (self->nknots);
      self->vl       = ncm_matrix_new (self->nknots, self->nknots);
      self->vr       = ncm_matrix_new (self->nknots, self->nknots);
      self->vlvr     = ncm_vector_new (self->nknots);
    }
  }
}

/**
 * nc_hiqg_1d_get_nknots:
 * @qg1d: a #NcHIQG1D
 *
 * Gets the current number of knots used in the wave function mesh.
 * 
 * Returns: the current number of knots used in the wave function mesh.
 */
guint 
nc_hiqg_1d_get_nknots (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->nknots;
}

typedef struct _NcHIQG1DInitCond
{
  NcHIQG1DPsi psi0_lnRS;
  gpointer psi_data;
} NcHIQG1DInitCond;

/**
 * nc_hiqg_1d_set_init_cond:
 * @qg1d: a #NcHIQG1D
 * @psi0_lnRS: (scope call): Initial wave-function in polar form
 * @psi_data: Initial wave-function data
 * @xi: initial point
 * @xf: final point
 * 
 * Sets the initial condition using @psi0, it calculates the best
 * mesh for the initial condition using the real part of @psi0.
 * 
 */
void
nc_hiqg_1d_set_init_cond (NcHIQG1D *qg1d, NcHIQG1DPsi psi0_lnRS, gpointer psi_data, const gdouble xi, const gdouble xf)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i;

  g_assert_cmpfloat (xi, <, xf);
  self->ti = 0.0;
  self->xi = xi;
  self->xf = xf;

  ncm_vector_set (self->fRePsi, 0, 0.0);
  ncm_vector_set (self->fImPsi, 0, 0.0);
  ncm_vector_set (self->frho,   0, 0.0);
  
  for (i = 0; i < self->nknots; i++)
  {
    const gdouble x = xi + (xf - xi) / (self->nknots * 1.0) * (i + 1);
    complex double ln_psi;
    complex double psi;

    ncm_vector_fast_set (self->knots, i, x);

    psi0_lnRS (psi_data, x, (gdouble *)&ln_psi);

    psi = cexp (ln_psi);

    ncm_vector_set (self->RePsi, i, creal (psi));
    ncm_vector_set (self->ImPsi, i, cimag (psi));

    ncm_vector_set (self->frho, i + 1, gsl_pow_2 (cabs (psi)));
  }

  ncm_spline_prepare (self->RePsi_s);
  ncm_spline_prepare (self->ImPsi_s);
  ncm_spline_prepare (self->rho_s);
}

/**
 * nc_hiqg_1d_set_init_cond_gauss:
 * @qg1d: a #NcHIQG1D
 * @qm_gauss: Initial wave-function data
 * @xi: initial point
 * @xf: final point
 * 
 * Sets the initial condition using @psi0 and nc_hiqg_1d_gauss_eval(), 
 * it calculates the best mesh for the initial condition using the real 
 * part of @psi0.
 * 
 */
void
nc_hiqg_1d_set_init_cond_gauss (NcHIQG1D *qg1d, NcHIQG1DGauss *qm_gauss, const gdouble xi, const gdouble xf)
{
  nc_hiqg_1d_set_init_cond (qg1d, (NcHIQG1DPsi) &nc_hiqg_1d_gauss_eval_lnRS, qm_gauss, xi, xf);
}

/**
 * nc_hiqg_1d_set_init_cond_exp:
 * @qg1d: a #NcHIQG1D
 * @qm_exp: Initial wave-function data
 * @xi: initial point
 * @xf: final point
 * 
 * Sets the initial condition using @psi0 and nc_hiqg_1d_exp_eval(), 
 * it calculates the best mesh for the initial condition using the real 
 * part of @psi0.
 * 
 */
void
nc_hiqg_1d_set_init_cond_exp (NcHIQG1D *qg1d, NcHIQG1DExp *qm_exp, const gdouble xi, const gdouble xf)
{
  nc_hiqg_1d_set_init_cond (qg1d, (NcHIQG1DPsi) &nc_hiqg_1d_exp_eval_lnRS, qm_exp, xi, xf);
}

gdouble
_nc_hiqg_1d_basis (const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  const gdouble h2   = h * h;
  const gdouble xy   = x * y;
  const gdouble xyh2 = xy * h2;
  const gdouble x2   = x * x;
  const gdouble y2   = y * y;
  const gdouble a1   = 0.5 * h2 * (x2 + y2);
  const gdouble xya  = pow (xy, a);

  if (xyh2 < 1.0)
  {
    return 2.0 * xya * sinh (xyh2) * exp (-a1);
  }
  else
  {
    return 2.0 * xya * exp (-a1 + gsl_sf_lnsinh (xyh2));
  }
}

/**
 * nc_hiqg_1d_basis:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @y: FIXME
 * @h: FIXME
 * @a: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_basis (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  return _nc_hiqg_1d_basis (x, y, h, a);
}

/*
 * READ HERE
 */


/*
 * Computes H . K(x, y)
 */
static gdouble
_nc_hiqg_1d_Hbasis (const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  const gdouble h2   = h * h;
  const gdouble xy   = x * y;
  const gdouble xyh2 = xy * h2;
  const gdouble x2   = x * x;
  const gdouble y2   = y * y;
  const gdouble a1   = 0.5 * h2 * (x2 + y2);
  const gdouble xya  = pow (xy, a);

  if (xyh2 < 1.0)
  {
    return 2.0 * h2 * xya * (2.0 * xyh2 * cosh (xyh2) + (1.0 + 2.0 * a - 2.0 * a1) * sinh (xyh2) + 2.0 * h2 * xyh2 * y2 * a * ncm_util_sinhx_m_xcoshx_x3 (xyh2)) * exp (-a1);
  }
  else
  {
    return 2.0 * h2 * xya * (2.0 * xyh2 * (1.0 - a / (h2 * x2)) * exp (-a1 + gsl_sf_lncosh (xyh2)) + (1.0 + 2.0 * a - 2.0 * a1 + 2.0 * a / (h2 * x2)) * exp (-a1 + gsl_sf_lnsinh (xyh2)));
  }
}

/**
 * nc_hiqg_1d_Hbasis:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @y: FIXME
 * @h: FIXME
 * @a: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_Hbasis (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  return _nc_hiqg_1d_Hbasis (x, y, h, a);
}

/*
 * Computes: K(x, y) / (x x^a)
 * 
 * It uses the function sinh1(x) = sinh(x)/x, such that when computing close
 * to the border x = 0 it avoids computing 0/0.
 * 
 */
static gdouble
_nc_hiqg_1d_basis_xxa (const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  const gdouble h2   = h * h;
  const gdouble xy   = x * y;
  const gdouble xyh2 = xy * h2;
  const gdouble x2   = x * x;
  const gdouble y2   = y * y;
  const gdouble a1   = 0.5 * h2 * (x2 + y2);
  const gdouble ya   = pow (y, a);

  if (xyh2 < 1.0)
  {
    return 2.0 * ya * y * h2 * ncm_util_sinh1 (xyh2) * exp (-a1);
  }
  else
  {
    return 2.0 * ya * exp (-a1 + gsl_sf_lnsinh (xyh2)) / x;
  }
}

/*
 * Computes: (K(x,y1)dK(x,y2)/dx - K(x,y2)dK(x,y1)/dx) / x^3
 * 
 * Computes the subtraction above (necessary to compute the derivative of the phase)
 * taking care of all the cancellations.
 * 
 */
static gdouble
_nc_hiqg_1d_Sbasis_x3 (const gdouble x, const gdouble y1, const gdouble y2, const gdouble h, const gdouble a)
{
  const gdouble h2      = h * h;
  const gdouble h2x     = h * h * x;
  const gdouble h8      = gsl_pow_4 (h2);
  const gdouble y12     = y1 * y1;
  const gdouble y22     = y2 * y2;
  const gdouble x2      = x * x;
  const gdouble y1my2   = y1 - y2;
  const gdouble y1py2   = y1 + y2;
  const gdouble y1my2_2 = y1my2 * y1my2;
  const gdouble y1py2_2 = y1py2 * y1py2;
  const gdouble y12my22 = y1py2 * y1my2;
  const gdouble a1      = 0.5 * h2 * (y12 + y22 + 2.0 * x2);
  const gdouble y1y2    = y1 * y2;
  const gdouble y1y2a   = pow (y1y2, a);

  if ((fabs (h2x * y1my2) < 1.0) && (fabs (h2x * y1py2) < 1.0))
  {
    return exp (-a1) * (h8 / 6.0) * y1y2a * y12my22 * (ncm_util_sinh3 (h2x * y1my2) * y1my2_2 - ncm_util_sinh3 (h2x * y1py2) * y1py2_2);
  }
  else
  {
    const gdouble x3 = x2 * x;
    if (y1 > y2)
      return (h2 / x3) * y1y2a * (+exp (-a1 + gsl_sf_lnsinh (+h2x * y1my2)) * y1py2 - exp (-a1 + gsl_sf_lnsinh (h2x * y1py2)) * y1my2);
    else
      return (h2 / x3) * y1y2a * (-exp (-a1 + gsl_sf_lnsinh (-h2x * y1my2)) * y1py2 - exp (-a1 + gsl_sf_lnsinh (h2x * y1py2)) * y1my2);
  }
}

/**
 * nc_hiqg_1d_Sbasis_x3:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @y1: FIXME
 * @y2: FIXME
 * @h: FIXME
 * @a: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_Sbasis_x3 (NcHIQG1D *qg1d, const gdouble x, const gdouble y1, const gdouble y2, const gdouble h, const gdouble a)
{
  return _nc_hiqg_1d_Sbasis_x3 (x, y1, y2, h, a);
}

/**
 * nc_hiqg_1d_get_lambda:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_get_lambda (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->lambda;
}

/**
 * nc_hiqg_1d_get_basis_a:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_get_basis_a (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->basis_a;
}

/**
 * nc_hiqg_1d_get_acs_a:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_get_acs_a (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->acs_a;
}

/**
 * nc_hiqg_1d_get_nu:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_get_nu (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->nu;
}

/**
 * nc_hiqg_1d_get_mu:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_get_mu (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->mu;
}

static gdouble
_nc_hiqg_1d_If2 (const gdouble y1, const gdouble y2, const gdouble h, const gdouble a)
{
  const gdouble h2      = h * h;
  const gdouble y12     = y1 * y1;
  const gdouble y22     = y2 * y2;
  const gdouble y1y2    = y1 * y2;
  const gdouble y1y2ha  = pow (y1y2 / h2, a);
  const gdouble a1      = 0.25 * h2 * (y12 + y22);
  gint signp;

  return exp (-a1 + lgamma_r (0.5 + a, &signp)) * y1y2ha / h * (gsl_sf_hyperg_1F1 (0.5 + a, 0.5, gsl_pow_2 (y1 + y2)) - gsl_sf_hyperg_1F1 (0.5 + a, 0.5, gsl_pow_2 (y1 - y2)));
}

static gdouble
_nc_hiqg_1d_Ixf2 (const gdouble y1, const gdouble y2, const gdouble h)
{
  const gdouble h2     = h * h;
  const gdouble ym12   = 0.5 * h * (y1 - y2);
  const gdouble yp12   = 0.5 * h * (y1 + y2);
  const gdouble ym12_2 = ym12 * ym12;
  const gdouble yp12_2 = yp12 * yp12;

  return (exp (-ym12_2) * yp12 * erf (yp12) - exp (-yp12_2) * ym12 * erf (ym12)) * ncm_c_sqrt_pi () / h2;
}

static void _nc_hiqg_1d_evol_C (NcHIQG1D *qg1d, const gdouble t);
static void _nc_hiqg_1d_prepare_splines (NcHIQG1D *qg1d);

static gint
_nc_hiqg_1d_bohm_f (gdouble t, N_Vector y, N_Vector ydot, gpointer user_data)
{
  NcHIQG1D *qg1d = NC_HIQG_1D (user_data);
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble *psi  = N_VGetArrayPointer (y);
  gdouble *dpsi = N_VGetArrayPointer (ydot);
  gint i;

  _nc_hiqg_1d_evol_C (qg1d, t);
  /*_nc_hiqg_1d_prepare_splines (qg1d);*/

  for (i = 0; i < self->nBohm; i++)
  {
    dpsi[i] = 2.0 * nc_hiqg_1d_eval_dS (qg1d, psi[i]);
  }

  return 0;
}

void
_nc_hiqg_1d_init_solver (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  const gdouble t0    = 0.0;
  const gdouble ini_q = nc_hiqg_1d_int_xrho_0_inf (qg1d);
  gint flag;
  gint i;

  if (self->bohm != NULL)
    ARKStepFree (&self->bohm);

  self->nBohm = 1;

  g_clear_pointer (&self->yBohm, N_VDestroy);
  self->yBohm = N_VNew_Serial (self->nBohm);

  for (i = 0; i < self->nBohm; i++)
  {
    NV_Ith_S (self->yBohm, i) = ini_q / (1.0 * self->nBohm) * (i + 1.0);
  }

  self->bohm = ARKStepCreate (&_nc_hiqg_1d_bohm_f, NULL, t0, self->yBohm);
  NCM_CVODE_CHECK (&self->bohm, "ARKodeCreate", 0, );

  flag = ARKStepSetUserData (self->bohm, (void *) qg1d);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

  flag = ARKStepSetMaxNumSteps (self->bohm, 10000);
  NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

  flag = ARKStepSStolerances (self->bohm, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );
}

/**
 * nc_hiqg_1d_prepare:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 * 
 * READ HERE
 *
 */
void
nc_hiqg_1d_prepare (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i, ret;

  ncm_vector_set (self->fknots, 0, 0.0);
  ncm_matrix_set (self->fpsi,   0, 0, 0.0);
  ncm_matrix_set (self->fpsi,   1, 0, 0.0);
  ncm_vector_set (self->frho,   0, 0.0);

  /* WARNING self->h must not be a integer multiple of y-x */
  self->h = 1.0 / (1.99 * (self->xf - self->xi) / (1.0 * self->nknots));
  ncm_matrix_set_zero (self->IM);

	for (i = 0; i < self->nknots; i++)
	{
		const gdouble xi    = ncm_vector_get (self->knots, i);
    const gdouble Ixixi = _nc_hiqg_1d_basis (xi, xi, self->h, self->basis_a);
    const gdouble Kxixi = /*Ixixi;//*/_nc_hiqg_1d_Hbasis (xi, xi, self->h, self->basis_a);
    guint j;

    ncm_matrix_set (self->IM, i, i, Ixixi);
    ncm_matrix_set (self->KM, i, i, Kxixi);
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g\n", h, xi, xi, Ixixi);*/
		for (j = i + 1; j < self->nknots; j++)
		{
			const gdouble xj    = ncm_vector_get (self->knots, j);
			const gdouble Ixixj = _nc_hiqg_1d_basis (xi, xj, self->h, self->basis_a);
			const gdouble Kxixj = /*Ixixj;//*/_nc_hiqg_1d_Hbasis (xi, xj, self->h, self->basis_a);
			const gdouble Kxjxi = /*Ixixj;//*/_nc_hiqg_1d_Hbasis (xj, xi, self->h, self->basis_a);

      /*printf ("% 22.15g % 22.15g\n", Kxixj, Kxjxi);*/

			ncm_matrix_set (self->IM, j, i, Ixixj);
      ncm_matrix_set (self->KM, i, j, Kxjxi);
      ncm_matrix_set (self->KM, j, i, Kxixj);
    }
	}

  /*ncm_vector_log_vals (self->knots, "#    KNOTS:  ", "% .5e", TRUE);*/
  /*ncm_vector_log_vals (self->hv,    "#       HV:  ", "% .5e", TRUE);*/
  /*printf ("#        H:  % 22.15g\n", self->h);*/

  {
    gint lwork  = -1;
    gchar equed = 'Y';
    gdouble rcond;
#ifdef HAVE_DSYSVXX_
    gint nparams = 3;
    gdouble params[3] = {1.0, 99.0, 0.0};
    gdouble rpvgrw;
#endif /* HAVE_DSYSVXX_ */
    gint ilo, ihi;
    gdouble abnrm;

    g_array_set_size (self->ipiv, self->nknots);
    g_array_set_size (self->iwork, 2 * self->nknots);

    if (self->work->len < 1)
      g_array_set_size (self->work, 1);

    /*ncm_matrix_log_vals (self->IM, "#     IM: ", "% .5e");*/
    /*ncm_matrix_log_vals (self->KM, "#     KM: ", "% .5e");*/

    ncm_matrix_set_zero (self->err_norm);
    ncm_matrix_set_zero (self->err_comp);

#ifdef HAVE_DSYSVXX_
    g_array_set_size (self->work, 4 * self->nknots);
    ret = ncm_lapack_dsysvxx ('E', 'L', self->nknots, self->nknots,
                              ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                              ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                              (gint *)self->ipiv->data, &equed,
                              ncm_vector_data (self->scaling),
                              ncm_matrix_data (self->KM), ncm_matrix_tda (self->KM),
                              ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                              &rcond, &rpvgrw,
                              ncm_vector_data (self->berr),
                              3, ncm_matrix_data (self->err_norm), ncm_matrix_data (self->err_comp), 0, NULL,
                              (gdouble *)self->work->data, (gint *)self->iwork->data);
#elif defined (HAVE_DSYSVX_)
    lwork = -1;
    equed = 'N';
    ret = ncm_lapack_dsysvx ('N', 'L', self->nknots, self->nknots,
                             ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                             ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                             (gint *)self->ipiv->data, 
                             ncm_matrix_data (self->KM), ncm_matrix_tda (self->KM),
                             ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                             &rcond, ncm_matrix_data (self->err_norm),
                             ncm_vector_data (self->berr),
                             (gdouble *)self->work->data, lwork, (gint *)self->iwork->data); 
    g_assert_cmpint (ret, ==, 0);

    lwork = g_array_index (self->work, gdouble, 0);
    if (lwork > self->work->len)
      g_array_set_size (self->work, lwork);

    ret = ncm_lapack_dsysvx ('N', 'L', self->nknots, self->nknots,
                             ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                             ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                             (gint *)self->ipiv->data, 
                             ncm_matrix_data (self->KM), ncm_matrix_tda (self->KM),
                             ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                             &rcond, ncm_matrix_data (self->err_norm),
                             ncm_vector_data (self->berr),
                             (gdouble *)self->work->data, lwork, (gint *)self->iwork->data); 
#endif /* HAVE_DSYSVXX_ */
    /*printf ("#    RCOND:  % 22.15g\n", rcond);*/
    /*printf ("#   RPVGRW:  % 22.15g\n", rpvgrw);*/

    /*ncm_vector_log_vals (self->scaling,  "#  SCALING:  ", "% .5e", TRUE);*/
    /*ncm_vector_log_vals (self->berr,     "#     BERR:  ", "% .5e", TRUE);*/
    /*ncm_matrix_log_vals (self->err_norm, "# ERR_NORM: ",  "% .5e");*/
    /*ncm_matrix_log_vals (self->err_comp, "# ERR_COMP: ",  "% .5e");*/
    /*ncm_matrix_log_vals (self->JM,       "#       JM: ",  "% .5e");*/

    g_assert_cmpint (ret, ==, 0);

    lwork = -1;
    ret = ncm_lapack_dgeevx ('B', 'V', 'V', 'B', self->nknots,
                             ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                             ncm_vector_data (self->wr), ncm_vector_data (self->wi),
                             ncm_matrix_data (self->vr), ncm_matrix_tda (self->vr),
                             ncm_matrix_data (self->vl), ncm_matrix_tda (self->vl),
                             &ilo, &ihi,
                             ncm_vector_data (self->Escaling),
                             &abnrm, ncm_matrix_data (self->err_norm), ncm_matrix_data (self->err_comp),
                             (gdouble *)self->work->data, lwork, (gint *)self->iwork->data);
    g_assert_cmpint (ret, ==, 0);

    lwork = g_array_index (self->work, gdouble, 0);
    if (lwork > self->work->len)
      g_array_set_size (self->work, lwork);

    ret = ncm_lapack_dgeevx ('B', 'V', 'V', 'B', self->nknots,
                             ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                             ncm_vector_data (self->wr), ncm_vector_data (self->wi),
                             ncm_matrix_data (self->vr), ncm_matrix_tda (self->vr),
                             ncm_matrix_data (self->vl), ncm_matrix_tda (self->vl),
                             &ilo, &ihi,
                             ncm_vector_data (self->Escaling),
                             &abnrm, ncm_matrix_data (self->err_norm), ncm_matrix_data (self->err_comp),
                             (gdouble *)self->work->data, lwork, (gint *)self->iwork->data);
    g_assert_cmpint (ret, ==, 0);

    /*ncm_vector_log_vals (self->wr, "#       WR:  ", "% .5e", TRUE);*/
    /*ncm_vector_log_vals (self->wi, "#       WI:  ", "% .5e", TRUE);*/

    /*ncm_matrix_log_vals (self->vl, "#       VL: ", "% .5e");*/
    /*ncm_matrix_log_vals (self->vr, "#       VR: ", "% .5e");*/

    /*ncm_vector_log_vals (self->ReC0, "# RePsi:  ", "% .5e", TRUE);*/
    /*ncm_vector_log_vals (self->ImC0, "# ImPsi:  ", "% .5e", TRUE);*/

    ncm_matrix_set_zero (self->err_norm);
    ncm_matrix_set_zero (self->err_comp);

    /*ncm_matrix_log_vals (self->psi, "#      PSI: ", "% .5e");*/
    /*printf ("#  PSI TDA:  %d\n", ncm_matrix_tda (self->psi));*/

#ifdef HAVE_DSYSVXX_
    ret = ncm_lapack_dsysvxx ('F', 'L', self->nknots, 2,
                              ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                              ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                              (gint *)self->ipiv->data, &equed,
                              ncm_vector_data (self->scaling),
                              ncm_matrix_data (self->psi), ncm_matrix_tda (self->psi),
                              ncm_matrix_data (self->C0), ncm_matrix_tda (self->C0),
                              &rcond, &rpvgrw,
                              ncm_vector_data (self->berr),
                              3, ncm_matrix_data (self->err_norm), ncm_matrix_data (self->err_comp), nparams, params,
                              (gdouble *)self->work->data, (gint *)self->iwork->data);
#elif defined (HAVE_DSYSVX_)
    lwork = -1;
    ret = ncm_lapack_dsysvx ('F', 'L', self->nknots, 2,
                             ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                             ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                             (gint *)self->ipiv->data, 
                             ncm_matrix_data (self->psi), ncm_matrix_tda (self->psi),
                             ncm_matrix_data (self->C0), ncm_matrix_tda (self->C0),
                             &rcond, ncm_matrix_data (self->err_norm),
                             ncm_vector_data (self->berr),
                             (gdouble *)self->work->data, lwork, (gint *)self->iwork->data);
    g_assert_cmpint (ret, ==, 0);
    
    lwork = g_array_index (self->work, gdouble, 0);
    if (lwork > self->work->len)
      g_array_set_size (self->work, lwork);
    
    ret = ncm_lapack_dsysvx ('F', 'L', self->nknots, 2,
                             ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                             ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                             (gint *)self->ipiv->data, 
                             ncm_matrix_data (self->psi), ncm_matrix_tda (self->psi),
                             ncm_matrix_data (self->C0), ncm_matrix_tda (self->C0),
                             &rcond, ncm_matrix_data (self->err_norm),
                             ncm_vector_data (self->berr),
                             (gdouble *)self->work->data, lwork, (gint *)self->iwork->data);
#endif /* HAVE_DSYSVXX_ */

    /*printf ("#    RCOND:  % 22.15g\n", rcond);*/
    /*printf ("#   RPVGRW:  % 22.15g\n", rpvgrw);*/

    /*ncm_vector_log_vals (self->scaling,  "#  SCALING:",  "% .5e", TRUE);*/
    /*ncm_vector_log_vals (self->berr,     "#     BERR:",  "% .5e", TRUE);*/
    /*ncm_matrix_log_vals (self->err_norm, "# ERR_NORM: ", "% .5e");*/
    /*ncm_matrix_log_vals (self->err_comp, "# ERR_COMP: ", "% .5e");*/

    g_assert_cmpint (ret, ==, 0);

    /*ncm_matrix_log_vals (self->C0, "#     C0: ", "% .5e");*/

    if (equed == 'Y')
    {
      for (i = 0; i < self->nknots; i++)
      {
        const gdouble si_m1 = 1.0 / ncm_vector_get (self->scaling, i);
        gint j;

        for (j = i; j < self->nknots; j++)
        {
          const gdouble sj_m1    = 1.0 / ncm_vector_get (self->scaling, j);
          const gdouble sisjIMij = ncm_matrix_get (self->IM, j, i);
          ncm_matrix_set (self->IM, j, i, sisjIMij * si_m1 * sj_m1);
        }
      }
    }

    for (i = 0; i < self->nknots; i++)
    {
      NcmVector *vr_i      = ncm_matrix_get_row (self->vr, i);
      NcmVector *vl_i      = ncm_matrix_get_row (self->vl, i);
      const gdouble n_i    = ncm_vector_dot (vl_i, vr_i);
      const gdouble ReC0_i = ncm_matrix_get (self->C0, 0, i);
      const gdouble ImC0_i = ncm_matrix_get (self->C0, 1, i);
      const gdouble ReA0_i = ncm_vector_dot (vl_i, self->ReC0) / n_i;
      const gdouble ImA0_i = ncm_vector_dot (vl_i, self->ImC0) / n_i;

      ncm_vector_set (self->vlvr, i, n_i);

      ncm_vector_set (self->A0, 2 * i + 0, ReA0_i);
      ncm_vector_set (self->A0, 2 * i + 1, ImA0_i);

      ncm_vector_set (self->ReC, i, ReC0_i);
      ncm_vector_set (self->ImC, i, ImC0_i);

      /*printf ("#      DOT:  % 22.15g % 22.15g % 22.15g\n", ncm_vector_dot (vl_i, vl_i), ncm_vector_dot (vr_i, vr_i), ncm_vector_dot (vl_i, vr_i));*/

      ncm_vector_free (vr_i);
      ncm_vector_free (vl_i);
    }
  }

  _nc_hiqg_1d_init_solver (qg1d);
}

/**
 * nc_hiqg_1d_peek_knots:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmVector *
nc_hiqg_1d_peek_knots (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->knots;
}

/**
 * nc_hiqg_1d_eval_ev:
 * @qg1d: a #NcHIQG1D
 * @i: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_eval_ev (NcHIQG1D *qg1d, const gint i, const gdouble x)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  g_assert_cmpint (i, <, self->nknots);
  {
    NcmVector *vr_i = ncm_matrix_get_row (self->vr, i);
    gdouble res = 0.0;
    gint j;

    for (j = 0; j < self->nknots; j++)
    {
      const gdouble xj = ncm_vector_get (self->knots, j);
      const gdouble cj = ncm_vector_get (vr_i, j);

      res += cj * _nc_hiqg_1d_basis (x, xj, self->h, self->basis_a);
    }

    ncm_vector_free (vr_i);

    return res;
  }
}

/**
 * nc_hiqg_1d_eval_psi0:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @psi0: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi_0$
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_eval_psi0 (NcHIQG1D *qg1d, const gdouble x, gdouble *psi0)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i;

  psi0[0] = 0.0;
  psi0[1] = 0.0;
  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi     = ncm_vector_get (self->knots, i);
    const gdouble ReC0_i = ncm_vector_get (self->ReC0, i);
    const gdouble ImC0_i = ncm_vector_get (self->ImC0, i);
    const gdouble Kxxi   = _nc_hiqg_1d_basis (x, xi, self->h, self->basis_a);

    psi0[0] += ReC0_i * Kxxi;
    psi0[1] += ImC0_i * Kxxi;
  }
}

static void
_nc_hiqg_1d_evol_C (NcHIQG1D *qg1d, const gdouble t)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  complex double *A0 = (complex double *) ncm_vector_data (self->A0);
  gint i;

  ncm_matrix_set_zero (self->C);
  for (i = 0; i < self->nknots; i++)
  {
    const gdouble wr = ncm_vector_get (self->wr, i);
    complex double A = A0[i] * cexp (-I * wr * t);
    cblas_daxpy (self->nknots, creal (A), ncm_matrix_ptr (self->vr, i, 0), 1, ncm_vector_data (self->ReC), 1);
    cblas_daxpy (self->nknots, cimag (A), ncm_matrix_ptr (self->vr, i, 0), 1, ncm_vector_data (self->ImC), 1);
  }
}

static void
_nc_hiqg_1d_prepare_splines (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint ret;
  gint i;

  ret = gsl_blas_dsymv (CblasLower, 1.0, ncm_matrix_gsl (self->IM), ncm_vector_gsl (self->ReC), 0.0, ncm_vector_gsl (self->RePsi));
  NCM_TEST_GSL_RESULT ("_nc_hiqg_1d_prepare_splines", ret);

  ret = gsl_blas_dsymv (CblasLower, 1.0, ncm_matrix_gsl (self->IM), ncm_vector_gsl (self->ImC), 0.0, ncm_vector_gsl (self->ImPsi));
  NCM_TEST_GSL_RESULT ("_nc_hiqg_1d_prepare_splines", ret);

  for (i = 0; i < self->nknots; i++)
  {
    const gdouble RePsi_i = ncm_vector_get (self->RePsi, i);
    const gdouble ImPsi_i = ncm_vector_get (self->ImPsi, i);

    ncm_vector_set (self->frho, i + 1, RePsi_i * RePsi_i + ImPsi_i * ImPsi_i);
  }

  ncm_spline_prepare (self->RePsi_s);
  ncm_spline_prepare (self->ImPsi_s);
  ncm_spline_prepare (self->rho_s);
}

/**
 * nc_hiqg_1d_evol:
 * @qg1d: a #NcHIQG1D
 * @t: FIXME
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_evol (NcHIQG1D *qg1d, const gdouble t)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  /*gdouble *psi = N_VGetArrayPointer (self->yBohm);*/
  gdouble ts = 0.0;
  gint flag;

  if (t > 0.0)
  {
    flag = ARKStepSetStopTime (self->bohm, t);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

    flag = ARKStepEvolve (self->bohm, t, self->yBohm, &ts, ARK_NORMAL);
    NCM_CVODE_CHECK (&flag, "ARKode", 1, );
  }

  _nc_hiqg_1d_evol_C (qg1d, t);
  _nc_hiqg_1d_prepare_splines (qg1d);
}

/**
 * nc_hiqg_1d_eval_psi:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi_0$
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_eval_psi (NcHIQG1D *qg1d, const gdouble x, gdouble *psi)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i;

  if (FALSE)
  {
    psi[0] = 0.0;
    psi[1] = 0.0;
    for (i = 0; i < self->nknots; i++)
    {
      const gdouble xi    = ncm_vector_get (self->knots, i);
      const gdouble ReC_i = ncm_vector_get (self->ReC, i);
      const gdouble ImC_i = ncm_vector_get (self->ImC, i);
      const gdouble Kxxi  = _nc_hiqg_1d_basis (x, xi, self->h, self->basis_a);

      psi[0] += ReC_i * Kxxi;
      psi[1] += ImC_i * Kxxi;
    }
  }
  else
  {
    psi[0] = ncm_spline_eval (self->RePsi_s, x);
    psi[1] = ncm_spline_eval (self->ImPsi_s, x);
  }
}

/**
 * nc_hiqg_1d_eval_dS:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_eval_dS (NcHIQG1D *qg1d, const gdouble x)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble dS_rho_x3 = 0.0;
  gdouble Re_psi_x  = 0.0;
  gdouble Im_psi_x  = 0.0;
  gint i, j;

  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi     = ncm_vector_get (self->knots, i);
    const gdouble ReC_i  = ncm_vector_get (self->ReC, i);
    const gdouble ImC_i  = ncm_vector_get (self->ImC, i);
    const gdouble Kxxi_x = _nc_hiqg_1d_basis_xxa (x, xi, self->h, self->basis_a);

    Re_psi_x += ReC_i * Kxxi_x;
    Im_psi_x += ImC_i * Kxxi_x;

    /*printf ("PSI_x % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", self->h, x, xi, Kxxi_x, Re_psi_x, Im_psi_x);*/

    for (j = i + 1; j < self->nknots; j++)
    {
      const gdouble xj    = ncm_vector_get (self->knots, j);
      const gdouble ReC_j = ncm_vector_get (self->ReC, j);
      const gdouble ImC_j = ncm_vector_get (self->ImC, j);
      const gdouble Sij   = _nc_hiqg_1d_Sbasis_x3 (x, xi, xj, self->h, self->basis_a);
      const gdouble Cij   = (ImC_j * ReC_i - ReC_j * ImC_i);

      /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", self->h, x, xi, xj, Sij, Cij, dS_rho_x3);*/
      dS_rho_x3 += Sij * Cij;
    }
  }
  /*printf ("# BLA % 22.15g % 22.15g % 22.15g\n", x * x * x * dS_rho_x3, x * x * (Re_psi_x * Re_psi_x + Im_psi_x * Im_psi_x), x * dS_rho_x3 / (Re_psi_x * Re_psi_x + Im_psi_x * Im_psi_x));*/

  if (FALSE)
  {
    const gdouble dS_rho   = 2.0 * x * x * x * dS_rho_x3;
    const gdouble dS_rho_s = ncm_spline_eval (self->RePsi_s, x) * ncm_spline_eval_deriv (self->ImPsi_s, x) - ncm_spline_eval_deriv (self->RePsi_s, x) * ncm_spline_eval (self->ImPsi_s, x);
    printf ("% 22.15g % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g\n", 
            x, 
            dS_rho, 
            dS_rho_s, 
            dS_rho_s / dS_rho, 
            ncm_spline_eval (self->RePsi_s, x) * ncm_spline_eval_deriv (self->ImPsi_s, x), 
            ncm_spline_eval_deriv (self->RePsi_s, x) * ncm_spline_eval (self->ImPsi_s, x)
            );
  }
  
  return 2.0 * x * dS_rho_x3 / (Re_psi_x * Re_psi_x + Im_psi_x * Im_psi_x);
}

/**
 * nc_hiqg_1d_int_rho_0_inf:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_int_rho_0_inf (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble int_rho = 0.0;

  if (FALSE)
  {
    gint i, j;

    for (i = 0; i < self->nknots; i++)
    {
      const gdouble xi       = ncm_vector_get (self->knots, i);
      const gdouble ReC_i    = ncm_vector_get (self->ReC, i);
      const gdouble ImC_i    = ncm_vector_get (self->ImC, i);
      const gdouble IKxixi_x = _nc_hiqg_1d_If2 (xi, xi, self->h, self->basis_a);

      int_rho += (ReC_i * ReC_i + ImC_i * ImC_i) * IKxixi_x;
      /*printf ("% 22.15g % 22.15g % 22.15g\n", xi, xi, IKxixi_x);*/

      for (j = i + 1; j < self->nknots; j++)
      {
        const gdouble xj       = ncm_vector_get (self->knots, j);
        const gdouble ReC_j    = ncm_vector_get (self->ReC, j);
        const gdouble ImC_j    = ncm_vector_get (self->ImC, j);
        const gdouble IKxixj_x = _nc_hiqg_1d_If2 (xi, xj, self->h, self->basis_a);
        const gdouble CRij     = 2.0 * (ReC_i * ReC_j + ImC_i * ImC_j);

        /*printf ("% 22.15g % 22.15g % 22.15g\n", xi, xj, IKxixj_x);*/

        int_rho += IKxixj_x * CRij;
      }
    }
  }
  else
    int_rho = ncm_spline_eval_integ (self->rho_s, 0.0, self->xf);

  return int_rho;
}

static gdouble
_nc_hiqg_1d_xrho (gdouble x, NcHIQG1DPrivate * const self)
{
  return x * ncm_spline_eval (self->rho_s, x);
}

static gdouble
_nc_hiqg_1d_mean_p_int (gdouble x, NcHIQG1DPrivate * const self)
{
  return ncm_spline_eval (self->RePsi_s, x) * ncm_spline_eval_deriv (self->ImPsi_s, x) - ncm_spline_eval_deriv (self->RePsi_s, x) * ncm_spline_eval (self->ImPsi_s, x);
}

/**
 * nc_hiqg_1d_int_xrho_0_inf:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_int_xrho_0_inf (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble int_xrho = 0.0;

  if (FALSE)
  {
    gint i, j;

    for (i = 0; i < self->nknots; i++)
    {
      const gdouble xi       = ncm_vector_get (self->knots, i);
      const gdouble ReC_i    = ncm_vector_get (self->ReC, i);
      const gdouble ImC_i    = ncm_vector_get (self->ImC, i);
      const gdouble IKxixi_x = _nc_hiqg_1d_Ixf2 (xi, xi, self->h);

      int_xrho += (ReC_i * ReC_i + ImC_i * ImC_i) * IKxixi_x;

      for (j = i + 1; j < self->nknots; j++)
      {
        const gdouble xj       = ncm_vector_get (self->knots, j);
        const gdouble ReC_j    = ncm_vector_get (self->ReC, j);
        const gdouble ImC_j    = ncm_vector_get (self->ImC, j);
        const gdouble IKxixi_x = _nc_hiqg_1d_Ixf2 (xi, xj, self->h);
        const gdouble CRij     = 2.0 * (ReC_i * ReC_j + ImC_i * ImC_j);

        int_xrho += IKxixi_x * CRij;
      }
    }
  }
  else
  {
    gsl_integration_workspace **w = ncm_integral_get_workspace ();
    gsl_function F;
    gdouble abserr;

    F.function = (gdouble (*)(gdouble, gpointer))_nc_hiqg_1d_xrho;
    F.params   = (gpointer) self;

    gsl_integration_qag (&F, 0.0, self->xf, self->abstol, self->reltol, NCM_INTEGRAL_PARTITION, 0, *w, &int_xrho, &abserr);
  }

  return int_xrho;
}

/**
 * nc_hiqg_1d_expect_p:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_expect_p (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble mean_p = 0.0;

  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gsl_function F;
  gdouble abserr;

  F.function = (gdouble (*)(gdouble, gpointer))_nc_hiqg_1d_mean_p_int;
  F.params   = (gpointer) self;

  gsl_integration_qag (&F, 0.0, self->xf, self->abstol, self->reltol, NCM_INTEGRAL_PARTITION, 0, *w, &mean_p, &abserr);

  return mean_p;
}

/**
 * nc_hiqg_1d_nBohm:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint
nc_hiqg_1d_nBohm (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;

  return self->nBohm;
}

/**
 * nc_hiqg_1d_Bohm:
 * @qg1d: a #NcHIQG1D
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_Bohm (NcHIQG1D *qg1d, gint i)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return NV_Ith_S (self->yBohm, i);
}

/**
 * nc_hiqg_1d_Bohm_p:
 * @qg1d: a #NcHIQG1D
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_Bohm_p (NcHIQG1D *qg1d, gint i)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  const gdouble a = NV_Ith_S (self->yBohm, i);
  return nc_hiqg_1d_eval_dS (qg1d, a);
}
