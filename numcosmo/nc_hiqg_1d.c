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
#include "math/memory_pool.h"
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
/*
#include <nvector/nvector_pthreads.h>
#include <nvector/nvector_openmp.h>
*/
#if HAVE_SUNDIALS_MAJOR == 3
#include <arkode/arkode.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <arkode/arkode_direct.h>
#include <arkode/arkode_spils.h>
#include <arkode/arkode_bandpre.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#endif /* HAVE_SUNDIALS_MAJOR == 3 */

#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIQG1DPrivate
{
  gdouble lambda;
  gdouble acs_a;
  gdouble basis_a;
  gdouble nu;
  gsize n;
  gsize np;
  gdouble abstol;
  gdouble reltol;
  complex double *weights;
  complex double *G;
  gdouble *nodes;
  gdouble *Jnu;
  NcmIntegral1d *Re_int;
  NcmIntegral1d *Im_int;
  /* PDE */
  NcmSpline *psi0_s;
  guint nknots;
  NcmVector *knots;
  NcmVector *spec_knots;
  NcmVector *spec_Y;
  NcmMatrix *spec_C0;
  NcmVector *spec_ReC0;
  NcmVector *spec_ImC0;
  NcmVector *spec_A0;
  NcmVector *spec_C;
  NcmVector *dS_v;
  NcmVector *rho_v;
  NcmVector *Re_psi_v;
  NcmVector *Im_psi_v;
  NcmSpline *dS_s;
  NcmSpline *rho_s;
  NcmSpline *Re_psi_s;
  NcmSpline *Im_psi_s;
  gdouble gamma;
  gboolean up_splines;
  gboolean noboundary;
  NcmSpline *YNp1_Re_R;
  NcmSpline *YNp1_Im_R;
  GArray *YNp1;
  gdouble h;
  gdouble aN;
  gdouble bN;
  gdouble cN;
  gdouble sqrt_aN_cN;
  gdouble ti;
  gdouble xi, xf;
  N_Vector y;
  N_Vector y2;
  N_Vector yt;
#if HAVE_SUNDIALS_MAJOR == 3
  SUNLinearSolver LS;
  SUNMatrix A;
#endif /* HAVE_SUNDIALS_MAJOR == 3 */
  gpointer arkode;
  gpointer bohm;
  gsl_dht *dht;
  NcmVector *spec_ReY;
  NcmVector *spec_ImY;
  NcmVector *spec_ReYt;
  NcmVector *spec_ImYt;
  NcmVector *specReW;
  NcmVector *specImW;
  GArray *work;
  GArray *ipiv;
  GArray *iwork;
  NcmVector *hv;
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
  PROP_NP,
  PROP_ABSTOL,
  PROP_RELTOL,
  PROP_NKNOTS,
  PROP_NOBOUNDARY,
};

G_DEFINE_TYPE (NcHIQG1D, nc_hiqg_1d, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcHIQG1DGauss, nc_hiqg_1d_gauss, nc_hiqg_1d_gauss_dup, nc_hiqg_1d_gauss_free);
G_DEFINE_BOXED_TYPE (NcHIQG1DExp,   nc_hiqg_1d_exp,   nc_hiqg_1d_exp_dup,   nc_hiqg_1d_exp_free);

static void
nc_hiqg_1d_init (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv = G_TYPE_INSTANCE_GET_PRIVATE (qg1d, NC_TYPE_HIQG_1D, NcHIQG1DPrivate);

  self->lambda  = 0.0;
  self->acs_a   = 0.0;
  self->basis_a = 0.0;
  self->nu      = 0.0;

  self->np      = 0;
  self->n       = 0;
  self->abstol  = 0.0;
  self->reltol  = 0.0;
  self->weights = NULL;
  self->G       = NULL;
  self->nodes   = NULL;
  self->Jnu     = NULL;
  self->Re_int  = NULL;
  self->Im_int  = NULL;

  /* PDE */
  self->psi0_s     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->nknots     = 0;
  self->knots      = NULL;
  self->spec_knots = NULL;
  self->spec_Y     = NULL;
  self->spec_C0    = NULL;
  self->spec_ReC0  = NULL;
  self->spec_ImC0  = NULL;
  self->spec_A0    = NULL;
  self->spec_C     = NULL;
  self->dS_v       = NULL;
  self->rho_v      = NULL;
  self->Re_psi_v   = NULL;
  self->Im_psi_v   = NULL;
  self->dS_s       = NULL;
  self->rho_s      = NULL;
  self->Re_psi_s   = NULL;
  self->Im_psi_s   = NULL;
  self->gamma      = 0.0;

  self->up_splines = FALSE;
  self->noboundary = FALSE;
  self->YNp1_Re_R  = ncm_spline_cubic_notaknot_new ();
  self->YNp1_Im_R  = ncm_spline_cubic_notaknot_new ();
  self->YNp1       = g_array_new (FALSE, FALSE, sizeof (complex double));
  self->ti         = 0.0;
  self->y          = NULL;
  self->y2         = NULL;
  self->yt         = NULL;
#if HAVE_SUNDIALS_MAJOR == 3
  self->LS         = NULL;
  self->A          = NULL;
#endif /* HAVE_SUNDIALS_MAJOR == 3 */
  self->arkode     = NULL;
  self->bohm       = NULL;

  self->dht        = NULL;
  self->spec_ReY   = NULL;
  self->spec_ImY   = NULL;
  self->spec_ReYt  = NULL;
  self->spec_ImYt  = NULL;
  self->specReW    = NULL;
  self->specImW    = NULL;

  self->work       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->ipiv       = g_array_new (FALSE, FALSE, sizeof (gint));
  self->iwork      = g_array_new (FALSE, FALSE, sizeof (gint));

  self->hv         = NULL;
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
      self->basis_a = 0.5 * ncm_util_sqrt1px_m1 (4.0 * self->lambda);
      self->nu      = sqrt (self->lambda + 0.25);
      break;
    case PROP_NP:
      self->np = g_value_get_uint (value);
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
    case PROP_NP:
      g_value_set_uint (value, self->np);
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

  g_clear_pointer (&self->nodes,   g_free);
  g_clear_pointer (&self->weights, g_free);
  g_clear_pointer (&self->G,       g_free);
  g_clear_pointer (&self->Jnu,     g_free);

  ncm_integral1d_clear (&self->Re_int);
  ncm_integral1d_clear (&self->Im_int);

  g_clear_pointer (&self->YNp1,    g_array_unref);

  g_clear_pointer (&self->y, N_VDestroy);
#if HAVE_SUNDIALS_MAJOR == 3
  g_clear_pointer (&self->LS, SUNLinSolFree);
  g_clear_pointer (&self->A, SUNMatDestroy);
#endif /* HAVE_SUNDIALS_MAJOR == 3 */

  ncm_vector_clear (&self->knots);
  ncm_vector_clear (&self->spec_knots);
  ncm_vector_clear (&self->spec_Y);
  ncm_matrix_clear (&self->spec_C0);
  ncm_vector_clear (&self->spec_ReC0);
  ncm_vector_clear (&self->spec_ImC0);
  ncm_vector_clear (&self->spec_A0);
  ncm_vector_clear (&self->spec_C);
  ncm_vector_clear (&self->dS_v);
  ncm_vector_clear (&self->rho_v);
  ncm_vector_clear (&self->Re_psi_v);
  ncm_vector_clear (&self->Im_psi_v);
  
  ncm_spline_clear (&self->psi0_s);

  ncm_spline_clear (&self->dS_s);
  ncm_spline_clear (&self->rho_s);

  ncm_spline_clear (&self->Re_psi_s);
  ncm_spline_clear (&self->Im_psi_s);

  ncm_spline_clear (&self->YNp1_Re_R);
  ncm_spline_clear (&self->YNp1_Im_R);

  ncm_vector_clear (&self->spec_ReY);
  ncm_vector_clear (&self->spec_ImY);

  ncm_vector_clear (&self->spec_ReYt);
  ncm_vector_clear (&self->spec_ImYt);

  ncm_vector_clear (&self->specReW);
  ncm_vector_clear (&self->specImW);

  g_clear_pointer (&self->work, g_array_unref);
  g_clear_pointer (&self->ipiv, g_array_unref);
  g_clear_pointer (&self->iwork, g_array_unref);

  ncm_vector_clear (&self->hv);
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

#if HAVE_SUNDIALS_MAJOR == 3
  if (self->arkode != NULL)
  {
    ARKodeFree (&self->arkode);
    self->arkode = NULL;
  }
  if (self->bohm != NULL)
  {
    ARKodeFree (&self->bohm);
    self->bohm = NULL;
  }
#endif /* HAVE_SUNDIALS_MAJOR == 3 */

  g_clear_pointer (&self->dht, gsl_dht_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiqg_1d_parent_class)->finalize (object);
}

static void
nc_hiqg_1d_class_init (NcHIQG1DClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIQG1DPrivate));

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
                                   PROP_NP,
                                   g_param_spec_uint ("np",
                                                      NULL,
                                                      "n_p",
                                                      10, G_MAXUINT, 20000,
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
      ncm_vector_clear (&self->hv);
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
      self->hv       = ncm_vector_new (self->nknots);
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

static gdouble 
_nc_hiqg_1d_Jnu (const gdouble nu, const gdouble z)
{
  if (nu >= 0.0)
    return gsl_sf_bessel_Jnu (nu, z);
  else
    return gsl_sf_bessel_Jnu (-nu, z) * cos (- nu * M_PI) - gsl_sf_bessel_Ynu (-nu, z) * sin (- nu * M_PI);
}

/**
 * nc_hiqg_1d_eval:
 * @qg1d: a #NcHIQG1D
 * @t: time difference
 * @x: $x$ coordinate
 * @y: $y$ coordinate
 * @G: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $G$
 * 
 * Calculates the propagator at $G (x,\;y;\;t)$.
 * 
 */
void
nc_hiqg_1d_eval (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble t, gdouble *G)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  const gdouble x2  = x * x;
  const gdouble y2  = y * y;
  const gdouble f1  = 0.5 * sqrt (x * y) / t;
  const gdouble f2  = _nc_hiqg_1d_Jnu (self->nu, 0.5 * x * y / t);
  complex double Gc = -I * f1 * f2 * cexp (I * (x2 + y2) * 0.25 / t - I * M_PI * 0.5 * self->nu);

  /*printf ("% 22.15g\n", f2);*/
  
  G[0] = creal (Gc);
  G[1] = cimag (Gc);
}

/**
 * nc_hiqg_1d_eval_array:
 * @qg1d: a #NcHIQG1D
 * @x: $x$ coordinate
 * @ya: $y$ coordinate
 * @n: number of elements 
 * @t: time difference
 * @G: $G$
 * 
 * Calculates the propagator at $G (x,\;y;\;t)$.
 * 
 */
void
nc_hiqg_1d_eval_array (NcHIQG1D *qg1d, const gdouble x, const gdouble *ya, gsize n, const gdouble t, gdouble *G)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  const gdouble x2   = x * x;
  complex double *Gc = (complex double *) G;

  gint i;

  if (FALSE)
  {
#ifdef HAVE_ACB_H
    acb_t nu, a, res;

    acb_init (nu);
    acb_init (a);
    acb_init (res);

    acb_set_d (nu, self->nu);

    for (i = 0; i < n; i++)
    {
      const gdouble arg = 0.5 * x * ya[i] / t;

      acb_set_d (a, arg);
      acb_hypgeom_bessel_j (res, nu, a, 120);

      self->Jnu[i] = arf_get_d (arb_midref (acb_realref (res)), ARF_RND_NEAR);

      /*printf ("% 22.15g % 22.15g\n", arg, _nc_hiqg_1d_Jnu (self->nu, arg) / self->Jnu[i] - 1.0);*/
    }

    acb_clear (nu);
    acb_clear (a);
    acb_clear (res);
#endif /* HAVE_ACB_H  */
  }
  else if (FALSE)
  {
    for (i = 0; i < n; i++)
    {
      const gdouble arg = 0.5 * x * ya[i] / t;
      self->Jnu[i] = arg;
    }
    
    gsl_sf_bessel_sequence_Jnu_e (self->nu, GSL_PREC_DOUBLE, n, self->Jnu);
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      const gdouble arg = 0.5 * x * ya[i] / t;
      self->Jnu[i] = _nc_hiqg_1d_Jnu (self->nu, arg);
    }
  }
  
  //gsl_sf_bessel_sequence_Jnu_e (self->nu, GSL_PREC_DOUBLE, n, self->Jnu); GSL_PREC_APPROX
  //gsl_sf_bessel_sequence_Jnu_e (self->nu, GSL_PREC_APPROX, n, self->Jnu); 
  
  for (i = 0; i < n; i++)
  {
    const gdouble y   = ya[i];
    const gdouble y2  = y * y;
    const gdouble f1  = 0.5 * sqrt (x * y) / t;

    Gc[i] = -I * f1 * self->Jnu[i] * cexp (I * (x2 + y2) * 0.25 / t - I * M_PI * 0.5 * self->nu); 
  }  
}

/**
 * nc_hiqg_1d_gauss_ini:
 * @qg1d: a #NcHIQG1D
 * @mean: Gaussian mean
 * @alpha: power-law
 * @sigma: standard deviation
 * @Hi: Initial $H(t_i) = H_i$
 * 
 * Calculates the propagator at $G (x,\;y;\;t)$.
 * 
 */
void
nc_hiqg_1d_gauss_ini (NcHIQG1D *qg1d, const gdouble mean, const gdouble alpha, const gdouble sigma, const gdouble Hi)
{
#ifdef HAVE_GSL_2_4
  NcHIQG1DPrivate * const self = qg1d->priv;
  NcHIQG1DGauss *qm_gauss = nc_hiqg_1d_gauss_new (mean, alpha, sigma, Hi);

  gsl_integration_fixed_workspace * ws = 
    gsl_integration_fixed_alloc (gsl_integration_fixed_hermite, self->np, 0.0, 0.25 / (sigma * sigma), alpha, 0.0);
  const gdouble *nodes   = gsl_integration_fixed_nodes (ws);
  const gdouble *weights = gsl_integration_fixed_weights (ws);
  gint i;

  self->n = 0;
  g_clear_pointer (&self->nodes,   g_free);
  g_clear_pointer (&self->weights, g_free);
  g_clear_pointer (&self->G,       g_free);
  g_clear_pointer (&self->Jnu,     g_free);

  self->nodes   = g_new0 (gdouble, self->np);
  self->weights = g_new0 (complex double, self->np);
  self->Jnu     = g_new0 (gdouble, self->np);
  self->G       = g_new0 (complex double, self->np);

  for (i = self->np / 2; i < self->np; i++)
  {
    if (weights[i] < self->abstol)
      break;
    else
    {
      const gdouble x = nodes[i];
      complex double psi0c;

      nc_hiqg_1d_gauss_eval_hermit (qm_gauss, x, (gdouble *)&psi0c);

      self->nodes[self->n]   = x;
      self->weights[self->n] = psi0c * weights[i];

      self->n++;
      /*printf ("%d % 22.15g % 22.15g | % 22.15g % 22.15g\n", i, nodes[i], weights[i], creal (psic), cimag (psic));*/
    }
  }

  nc_hiqg_1d_gauss_free (qm_gauss);
  gsl_integration_fixed_free (ws);
#endif /* HAVE_GSL_2_4 */
}

typedef struct _NcHIQG1DInt
{
  NcHIQG1D *qg1d;
  NcHIQG1DPrivate * const self;
  const gdouble t;
  const gdouble x;
  const gdouble alpha;
  const gdouble sigma;
  const gdouble Hi;
} NcHIQG1DInt;

static gdouble 
_nc_hiqg_1d_propto_Re_integ (gpointer userdata, const gdouble y, const gdouble w)
{
  NcHIQG1DInt *integ   = (NcHIQG1DInt *) userdata;
  gint signp            = 0;
  const gdouble x       = integ->x;
  const gdouble t       = integ->t;
  const gdouble alpha   = integ->alpha;
  const gdouble sigma   = integ->sigma;
  const gdouble Hi      = integ->Hi;
  const gdouble lny     = log (y);
  const gdouble x2      = x * x;
  const gdouble y2      = y * y;
  const gdouble ap12    = alpha + 0.5;
  const gdouble f1      = (1.0 - 2.0 * alpha) * 0.25 * M_LN2;
  const gdouble lnsigma = log (sigma);
  const gdouble lng     = lgamma_r (ap12, &signp);
  complex double psi0   = cexp (f1 - ap12 * lnsigma + (alpha - 1.0) * lny - 0.5 * lng + 0.5 * y2 * I * Hi);
  const gdouble f2      = 0.5 * sqrt (x * y) / t;
  complex double G      = -I * f2 *  _nc_hiqg_1d_Jnu (integ->self->nu, 0.5 * x * y / t) * cexp (I * (x2 + y2) * 0.25 / t - I * M_PI * 0.5 * integ->self->nu); 

  /*printf ("% 22.15g % 22.15g % 22.15g\n", y, creal (psi0 * G), cimag (psi0 * G));*/
  
  return creal (psi0 * G);
}

static gdouble 
_nc_hiqg_1d_propto_Im_integ (gpointer userdata, const gdouble y, const gdouble w)
{
  NcHIQG1DInt *integ   = (NcHIQG1DInt *) userdata;
  gint signp            = 0;
  const gdouble x       = integ->x;
  const gdouble t       = integ->t;
  const gdouble alpha   = integ->alpha;
  const gdouble sigma   = integ->sigma;
  const gdouble Hi      = integ->Hi;
  const gdouble lny     = log (y);
  const gdouble x2      = x * x;
  const gdouble y2      = y * y;
  const gdouble ap12    = alpha + 0.5;
  const gdouble f1      = (1.0 - 2.0 * alpha) * 0.25 * M_LN2;
  const gdouble lnsigma = log (sigma);
  const gdouble lng     = lgamma_r (ap12, &signp);
  complex double psi0   = cexp (f1 - ap12 * lnsigma + (alpha - 1.0) * lny - 0.5 * lng + 0.5 * y2 * I * Hi);
  const gdouble f2      = 0.5 * sqrt (x * y) / t;
  complex double G      = -I * f2 *  _nc_hiqg_1d_Jnu (integ->self->nu, 0.5 * x * y / t) * cexp (I * (x2 + y2) * 0.25 / t - I * M_PI * 0.5 * integ->self->nu); 

  return cimag (psi0 * G);
}

/**
 * nc_hiqg_1d_propto:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @t: FIXME
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $G$
 * 
 * Calculates the propagator at $G (x,\;y;\;t)$.
 * 
 */
void
nc_hiqg_1d_propto (NcHIQG1D *qg1d, const gdouble x, const gdouble t, gdouble *psi)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  complex double psic = 0.0;
  gint i;

  if (TRUE)
  {
    nc_hiqg_1d_eval_array (qg1d, x, self->nodes, self->n, t, (gdouble *)self->G);

    for (i = 0; i < self->n; i++)
    {
      psic += self->weights[i] * self->G[i];
      //printf ("%d % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g\n", i, self->nodes[i], creal (self->weights[i]), cimag (self->weights[i]), creal (psic), cimag (psic));
    }

    psi[0] = creal (psic);
    psi[1] = cimag (psic);
  }
  else if (FALSE)
  {
    const gdouble alpha = 1.0;
    const gdouble sigma = 1.0;
    const gdouble Hi    = 1.0;
    NcHIQG1DInt integ = {qg1d, self, t, x, alpha, sigma, Hi};
    gdouble err = 0.0;
    
    ncm_integral1d_clear (&self->Re_int);
    ncm_integral1d_clear (&self->Im_int);

    self->Re_int = NCM_INTEGRAL1D (ncm_integral1d_ptr_new (&_nc_hiqg_1d_propto_Re_integ, NULL));
    self->Im_int = NCM_INTEGRAL1D (ncm_integral1d_ptr_new (&_nc_hiqg_1d_propto_Im_integ, NULL));

    ncm_integral1d_set_reltol (self->Re_int, 1.0e-5);
    ncm_integral1d_set_reltol (self->Im_int, 1.0e-5);
    
    ncm_integral1d_ptr_set_userdata (NCM_INTEGRAL1D_PTR (self->Re_int), &integ);
    ncm_integral1d_ptr_set_userdata (NCM_INTEGRAL1D_PTR (self->Im_int), &integ);

    psi[0] = ncm_integral1d_eval_gauss_hermite1_r_p (self->Re_int, 0.25 / (sigma * sigma), &err);
    psi[1] = ncm_integral1d_eval_gauss_hermite1_r_p (self->Im_int, 0.25 / (sigma * sigma), &err);
  }
}

static gdouble 
_nc_hiqg_1d_propto_norm (const gdouble x, gpointer p)
{
  NcHIQG1DInt *integ = (NcHIQG1DInt *) p;
  complex double psi;

  nc_hiqg_1d_propto (integ->qg1d, x, integ->t, (gdouble *)&psi);

  return creal (psi) * creal (psi) + cimag (psi) * cimag (psi);
}

/**
 * nc_hiqg_1d_propto_norm:
 * @qg1d: a #NcHIQG1D
 * @t: FIXME
 * 
 * Calculates the propagator at $G (x,\;y;\;t)$.
 * 
 */
gdouble
nc_hiqg_1d_propto_norm (NcHIQG1D *qg1d, const gdouble t)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  NcHIQG1DInt integ = {qg1d, self, t};
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gdouble result, abserr;
  gsl_function F;

  F.params   = &integ;
  F.function = &_nc_hiqg_1d_propto_norm;

  gsl_integration_qagiu (&F, 0.0, 1.0e-7, 1.0e-7, NCM_INTEGRAL_PARTITION, *w, &result, &abserr);

  ncm_memory_pool_return (w);

  return result;
}

typedef struct _NcHIQG1DInitCond
{
  NcHIQG1DPsi psi0_lnRS;
  gpointer psi_data;
} NcHIQG1DInitCond;

static gdouble
_nc_hiqg_1d_set_init_cond_real (gdouble x, gpointer p)
{
  NcHIQG1DInitCond *ic = (NcHIQG1DInitCond *) p;
  gdouble lnRS[2];

  ic->psi0_lnRS (ic->psi_data, x, lnRS);
  
  return creal (cexp (lnRS[0] + I * lnRS[1]));
}

static complex double
_nc_hiqg_1d_set_init_cond_complex (gdouble x, gpointer p)
{
  NcHIQG1DInitCond *ic = (NcHIQG1DInitCond *) p;
  gdouble lnRS[2];

  ic->psi0_lnRS (ic->psi_data, x, lnRS);

  return cexp (lnRS[0] + I * lnRS[1]);
}

void _nc_hiqg_1d_init_solver (NcHIQG1D *qg1d);
void _nc_hiqg_1d_spec_init_solver (NcHIQG1D *qg1d);
void _nc_hiqg_1d_init_spec_solver (NcHIQG1D *qg1d, NcHIQG1DInitCond *ic);

complex double
_nc_hiqg_1d_get_YN (NcHIQG1D *qg1d, const gdouble t)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  if (self->YNp1->len == 0)
  {
    return 0.0;
  }

  if (t == 0.0)
  {
    return g_array_index (self->YNp1, complex double, 0);
  }
  else
  {
    complex double A      = 1.0 / self->sqrt_aN_cN;
    complex double T      = -I / self->sqrt_aN_cN;
    complex double YN     = 0.0;
    const gdouble arg     = t * self->cN * self->sqrt_aN_cN;
    /*const gdouble abstol  = GSL_DBL_EPSILON;*/
    const gdouble epsilon = 1.0;
    const gint last_m     = self->YNp1->len - 1;
    const gint asym_m     = GSL_MAX (gsl_pow_2 (2.0 * arg / epsilon), 1);
    gdouble Jn, Jnm1, Jnp1;
    gint n, m;

    if (asym_m > last_m)
    {
      m = last_m;
    }
    else
    {
      const gdouble logabstol = log (self->abstol);
      const gdouble logarg    = log (arg);
      gint newm, oldm = asym_m;
      gint c = 0;
      gint signp;

      while (TRUE)
      {
        newm =  logabstol / (logarg - lgamma_r (oldm, &signp) / oldm);
        c++;
        if (abs (oldm - newm) <= 1)
        {
          newm = GSL_MAX (oldm, newm);
          m    = GSL_MIN (last_m, newm);
          break;
        }
        else
        {
          oldm = newm;
        }
      }
      
      printf ("[%3d] %3d %3d %3d | %3d %3d % 22.15g % 22.15g\n", c, m, asym_m, last_m, newm, oldm, jn (newm, 2.0 * arg), jn (oldm, 2.0 * arg));
    }

    Jn   = jn (m + 0, 2.0 * arg);
    Jnp1 = jn (m + 1, 2.0 * arg);
    
    for (n = m; n > 0; n--)
    {
      const gdouble t1        = n * Jn / (t * self->cN);
      const complex double dY = g_array_index (self->YNp1, complex double, n - 1) * A * t1;

      YN += dY;
/*
      if (
          (fabs (creal (dY) / creal (YN)) < abstol) && 
          (fabs (cimag (dY) / cimag (YN)) < abstol)
          )
        break;
*/
      /*printf ("%d % 22.15g % 22.15g\n", n, cabs (YN), cabs (g_array_index (self->YNp1, complex double, i) * A * t1));*/
      /*printf ("%d % 22.15g % 22.15g %e | % 22.15g\n", n, Jn, jn (n, 2.0 * arg), Jn / jn (n, 2.0 * arg) - 1.0, cabs (YN));*/

      Jnm1 = n * 1.0 / arg * Jn - Jnp1;
      Jnp1 = Jn;
      Jn   = Jnm1;
      
      A *= T;
    }

    return YN * cexp (-2.0 * I * t * self->bN);
  }
}

gdouble
_nc_hiqg_1d_get_Re_YN (const gdouble t, gpointer qg1d)
{
  /*printf ("% 22.15g % 22.15g\n", t, creal (_nc_hiqg_1d_get_YN (qg1d, t)));*/
  return creal (_nc_hiqg_1d_get_YN (qg1d, t));
}

gdouble
_nc_hiqg_1d_get_Im_YN (const gdouble t, gpointer qg1d)
{
  return cimag (_nc_hiqg_1d_get_YN (qg1d, t));
}

static void _nc_hiqg_1d_prepare_splines (NcHIQG1DPrivate * const self, const gdouble t, NcmVector *knots, gdouble *Y);

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
  NcHIQG1DInitCond ic = {psi0_lnRS, psi_data};
  gboolean use_spline  = FALSE;
  gsl_function F;
  gint i;

  g_assert_cmpfloat (xi, <, xf);

  self->ti = 0.0;
  self->xi = xi;
  self->xf = xf;

  ncm_vector_clear (&self->knots);
  
  if (use_spline)
  {
    F.params   = &ic;
    F.function = _nc_hiqg_1d_set_init_cond_real;
    ncm_spline_set_func (self->psi0_s, NCM_SPLINE_FUNCTION_SPLINE, &F, xi, xf, 0, 1.0e-2/*self->reltol*/);

    self->knots  = ncm_spline_get_xv (self->psi0_s);
    self->nknots = ncm_vector_len (self->knots);
  }
  else
  {
    self->knots         = ncm_vector_new (self->nknots);
    const guint lp      = (self->nknots / 10) > 1 ? (self->nknots / 10) : 2;
    const gint np       = self->nknots - lp;
    const gdouble fxi   = xf / (np * 1.0);
    const gdouble lnfxi = log (fxi);
    const gdouble lnxi  = (xi == 0.0) ? (lnfxi - 4.0 * M_LN10) : log (xi);

    g_assert_cmpfloat (lnfxi, >, lnxi);
    g_assert_cmpfloat (xi, >=, 0.0);
    g_assert_cmpfloat (xf, >, xi);

    if (FALSE)
    {
      ncm_vector_fast_set (self->knots, 0, exp (lnxi));
      for (i = 0; i < lp; i++)
      {
        const gdouble lnx = lnxi + (lnfxi - lnxi) / (1.0 * lp) * i;
        const gdouble x   = exp (lnx);

        ncm_vector_fast_set (self->knots, i, x);
      }

      for (i = 0; i < np; i++)
      {
        const gdouble x = fxi * (i + 1.0);

        ncm_vector_fast_set (self->knots, i + lp, x);
      }

      /*ncm_vector_log_vals (self->knots, "# KNOTS: ", "% 22.15g", TRUE);*/
    }
    else
    {
      for (i = 0; i < self->nknots; i++)
      {
        const gdouble x = xi + (xf - xi) / (self->nknots * 1.0) * (i + 1);
        ncm_vector_fast_set (self->knots, i, x);
      }
    }
  }
  
  {
    ncm_vector_clear (&self->spec_knots);
    ncm_vector_clear (&self->spec_Y);
    ncm_matrix_clear (&self->spec_C0);
    ncm_vector_clear (&self->spec_ReC0);
    ncm_vector_clear (&self->spec_ImC0);
    ncm_vector_clear (&self->spec_A0);
    ncm_vector_clear (&self->spec_C);

    ncm_vector_clear (&self->dS_v);
    ncm_vector_clear (&self->rho_v);

    ncm_vector_clear (&self->Re_psi_v);
    ncm_vector_clear (&self->Im_psi_v);

    ncm_spline_clear (&self->dS_s);
    ncm_spline_clear (&self->rho_s);

    ncm_spline_clear (&self->Re_psi_s);
    ncm_spline_clear (&self->Im_psi_s);
    
    self->rho_s      = ncm_spline_cubic_notaknot_new ();
    self->dS_s       = ncm_spline_cubic_notaknot_new ();

    self->Re_psi_s   = ncm_spline_cubic_notaknot_new ();
    self->Im_psi_s   = ncm_spline_cubic_notaknot_new ();

    self->spec_knots = ncm_vector_new (self->nknots);
    self->spec_Y     = ncm_vector_new (self->nknots * 2);
    self->spec_C0    = ncm_matrix_new (2, self->nknots);
    self->spec_A0    = ncm_vector_new (self->nknots * 2);
    self->spec_C     = ncm_vector_new (self->nknots * 2);
    self->spec_ReY   = ncm_vector_get_subvector_stride (self->spec_Y, 0, self->nknots, 2);
    self->spec_ImY   = ncm_vector_get_subvector_stride (self->spec_Y, 1, self->nknots, 2);
    self->spec_ReC0  = ncm_matrix_get_row (self->spec_C0, 0);
    self->spec_ImC0  = ncm_matrix_get_row (self->spec_C0, 1);

    self->rho_v      = ncm_vector_new (self->nknots);
    self->dS_v       = ncm_vector_new (self->nknots);

    self->Re_psi_v   = ncm_vector_new (self->nknots);
    self->Im_psi_v   = ncm_vector_new (self->nknots);
    
    ncm_spline_set (self->rho_s, self->spec_knots, self->rho_v, FALSE);
    ncm_spline_set (self->dS_s,  self->spec_knots, self->dS_v,  FALSE);

    ncm_spline_set (self->Re_psi_s, self->spec_knots, self->Re_psi_v, FALSE);
    ncm_spline_set (self->Im_psi_s, self->spec_knots, self->Im_psi_v, FALSE);

    self->up_splines = FALSE;
  }

  g_clear_pointer (&self->y, N_VDestroy);
  self->y = N_VNew_Serial (3 * self->nknots);
  /*self->y = N_VNew_OpenMP (3 * self->nknots, 20);*/
  NCM_CVODE_CHECK (&self->y, "N_VNew_Serial", 0, );

  if (FALSE)
  {
    complex double *Y = (complex double *) N_VGetArrayPointer (self->y);

    printf ("# USING %d knots!\n", self->nknots);
    
    for (i = 0; i < self->nknots - 1; i++)
    {
      const gdouble x = ncm_vector_fast_get (self->knots, i);
      Y[i] = _nc_hiqg_1d_set_init_cond_complex (x, &ic);
    }

    if (!self->noboundary)
      Y[self->nknots - 1] = 0.0;
    else
    {
      const gdouble xNm1 = ncm_vector_fast_get (self->knots, self->nknots - 1);
      const gdouble xNm2 = ncm_vector_fast_get (self->knots, self->nknots - 2);
      const gdouble h    = xNm1 - xNm2;
      const gdouble xN   = xNm1 + h;
      gdouble x          = xN;
      gdouble absY;

      g_array_set_size (self->YNp1, 0);
      self->h          = h;
      self->aN         = -1.0 / (h * h);
      self->cN         = -1.0 / (h * h);
      self->bN         = +2.0 / (h * h) + self->lambda / gsl_pow_2 (xNm1);
      self->sqrt_aN_cN = sqrt (self->aN / self->cN);

      Y[self->nknots - 1] = _nc_hiqg_1d_set_init_cond_complex (xNm1, &ic);
      
      while (TRUE)
      {
        const complex double Yx = _nc_hiqg_1d_set_init_cond_complex (x, &ic);
        absY = cabs (Yx);

        if (absY > self->abstol)
        {
          g_array_append_val (self->YNp1, Yx);
          x += h;
        }
        else
          break;
      }

      if (FALSE)
      {
        gsl_function F;

        F.params   = qg1d;
        
        F.function = &_nc_hiqg_1d_get_Re_YN;
        ncm_spline_set_func (self->YNp1_Re_R, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, 10.0, 0, 1.0e-2);

        F.function = &_nc_hiqg_1d_get_Im_YN;
        ncm_spline_set_func (self->YNp1_Im_R, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, 10.0, 0, 1.0e-2);      
      }
    }
  }
  else
  {
    gdouble *Y = (gdouble *) N_VGetArrayPointer (self->y);
    
    for (i = 0; i < self->nknots; i++)
    {
      const gdouble x = ncm_vector_fast_get (self->knots, i);
      complex double psi;

      Y[_XI (i, self->nknots)] = x;
      psi0_lnRS (psi_data, x, &Y[_LNRI (i, self->nknots)]);

      psi = cexp (*((complex double *)&Y[_LNRI (i, self->nknots)]));

      ncm_vector_set (self->spec_Y, 2 * i + 0, creal (psi));
      ncm_vector_set (self->spec_Y, 2 * i + 1, cimag (psi));

      /*printf ("% 22.15g % 22.15g % 22.15g\n", Y[_XI (i, self->nknots)], Y[_LNRI (i, self->nknots)], Y[_SI (i, self->nknots)]);*/
    }

    /* _nc_hiqg_1d_prepare_splines (self, self->ti, self->knots, Y); */
  }

  if (FALSE)
    _nc_hiqg_1d_init_spec_solver (qg1d, &ic);
  _nc_hiqg_1d_init_solver (qg1d);
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

static gint _nc_hiqg_1d_f_impl (gdouble t, N_Vector y, N_Vector ydot, gpointer user_data);
static gint _nc_hiqg_1d_f_expl (gdouble t, N_Vector y, N_Vector ydot, gpointer user_data);
/*static gint _nc_hiqg_1d_J (N_Vector v, N_Vector Jv, gdouble t, N_Vector y, N_Vector fy, gpointer user_data, N_Vector tmp);*/

void
_nc_hiqg_1d_init_solver (NcHIQG1D *qg1d)
{
#if HAVE_SUNDIALS_MAJOR == 3
  NcHIQG1DPrivate * const self = qg1d->priv;
  const gdouble t0 = 0.0;
  gint flag;

  g_clear_pointer (&self->LS, SUNLinSolFree);

  if (self->arkode != NULL)
    ARKodeFree (&self->arkode);

  self->arkode = ARKodeCreate ();
  NCM_CVODE_CHECK (&self->arkode, "ARKodeCreate", 0, );

  flag = ARKodeInit (self->arkode, &_nc_hiqg_1d_f_expl, &_nc_hiqg_1d_f_impl, t0, self->y);
  NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

  flag = ARKodeSetUserData (self->arkode, (void *) qg1d);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

  flag = ARKodeSetMaxNumSteps (self->arkode, 10000);
  NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

  flag = ARKodeSStolerances (self->arkode, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );

  //flag = ARKodeSetARKTableNum (self->arkode, ARK324L2SA_DIRK_4_2_3, ARK324L2SA_ERK_4_2_3);
  //flag = ARKodeSetARKTableNum (self->arkode, ARK548L2SA_DIRK_8_4_5, ARK548L2SA_ERK_8_4_5);
  //NCM_CVODE_CHECK (&flag, "ARKodeSetARKTableNum", 1, );

  //flag = ARKodeSetOrder (self->arkode, 3);
  //NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );
  
  //flag = ARKodeSetAdaptivityMethod (self->arkode, 2, 1, 0, NULL);
  //NCM_CVODE_CHECK (&flag, "ARKodeSetAdaptivityMethod", 1, );

  //flag = ARKodeSetPredictorMethod (self->arkode, 4);
  //NCM_CVODE_CHECK (&flag, "ARKodeSetPredictorMethod", 1, );

  //flag = ARKodeSetLinear (self->arkode, 0);
  //NCM_CVODE_CHECK (&flag, "ARKodeSetLinear", 1, );

  if (TRUE)
  {
    self->LS = SUNSPGMR (self->y, PREC_NONE, self->nknots);
    NCM_CVODE_CHECK (&flag, "SUNSPGMR", 1, );  

/*    self->LS = SUNPCG (self->y, 0, self->nknots);
    NCM_CVODE_CHECK (&flag, "SUNPCG", 1, );
*/
    flag = ARKSpilsSetLinearSolver (self->arkode, self->LS);
    NCM_CVODE_CHECK (&flag, "ARKSpilsSetLinearSolver", 1, );

/*
     flag = ARKSpilsSetJacTimes (self->arkode, NULL, _nc_hiqg_1d_J);
     NCM_CVODE_CHECK (&flag, "ARKSpilsSetJacTimes", 1, );
*/

    flag = ARKBandPrecInit (self->arkode, 3 * self->nknots, 8, 8);
    NCM_CVODE_CHECK (&flag, "ARKBandPrecInit", 1, );

  }
  else if (TRUE)
  {
    self->A = SUNBandMatrix (3 * self->nknots, 12, 12, 24);
    NCM_CVODE_CHECK ((gpointer)self->A, "SUNBandMatrix", 0, );

    self->LS = SUNBandLinearSolver (self->y, self->A);
    NCM_CVODE_CHECK ((gpointer)self->LS, "SUNBandLinearSolver", 0, );

    flag = ARKDlsSetLinearSolver (self->arkode, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetLinearSolver", 1, );

    /* Use a difference quotient Jacobian */
    flag = ARKDlsSetJacFn (self->arkode, NULL);
    NCM_CVODE_CHECK (&flag, "ARKDlsSetJacFn", 1, );
    
  }
  else
  {
    flag = ARKodeSetFixedPoint (self->arkode, 100);
    NCM_CVODE_CHECK (&flag, "ARKodeSetFixedPoint", 1, );
  }
#endif /* HAVE_SUNDIALS_MAJOR == 3 */
}

void
_nc_hiqg_1d_init_spec_solver (NcHIQG1D *qg1d, NcHIQG1DInitCond *ic)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  complex double *Y = (complex double *) ncm_vector_data (self->spec_Y);
  gint i;

  g_clear_pointer (&self->dht, gsl_dht_free);
  
  ncm_vector_clear (&self->spec_ReY);
  ncm_vector_clear (&self->spec_ImY);

  ncm_vector_clear (&self->spec_ReYt);
  ncm_vector_clear (&self->spec_ImYt);

  ncm_vector_clear (&self->specReW);
  ncm_vector_clear (&self->specImW);

  self->dht = gsl_dht_alloc (self->nknots - 1);
  gsl_dht_init (self->dht, self->nu, self->xf);

  self->spec_ReY    = ncm_vector_new (self->nknots - 1);
  self->spec_ImY    = ncm_vector_new (self->nknots - 1);

  self->spec_ReYt   = ncm_vector_new (self->nknots - 1);
  self->spec_ImYt   = ncm_vector_new (self->nknots - 1);

  self->specReW    = ncm_vector_new (self->nknots - 1);
  self->specImW    = ncm_vector_new (self->nknots - 1);

  ncm_vector_fast_set (self->spec_knots, 0, 0.0);
  ncm_vector_set      (self->Re_psi_v,   0, 0.0);
  ncm_vector_set      (self->Im_psi_v,   0, 0.0);

  Y[0] = 0.0;
  for (i = 0; i < self->nknots - 1; i++)
  {
    const gdouble x   = gsl_dht_x_sample (self->dht, i);
    complex double Yi = _nc_hiqg_1d_set_init_cond_complex (x, ic);
    complex double Yn = Yi / sqrt (x);

    ncm_vector_fast_set (self->spec_knots, i + 1, x);
    Y[i + 1] = Yi;

    ncm_vector_fast_set (self->spec_ReY, i, creal (Yn));
    ncm_vector_fast_set (self->spec_ImY, i, cimag (Yn));
  }
  
  gsl_dht_apply (self->dht, ncm_vector_data (self->spec_ReY), ncm_vector_data (self->specReW));
  gsl_dht_apply (self->dht, ncm_vector_data (self->spec_ImY), ncm_vector_data (self->specImW));

  /*_nc_hiqg_1d_prepare_splines (self, self->ti, self->spec_knots, ncm_vector_data (self->spec_Y));*/
}

/*
static gdouble
_nc_hiqg_1d_diff (const gdouble fp2, const gdouble fp1, const gdouble f, const gdouble dx2, const gdouble dx1)
{
  const gdouble ddx = dx2 - dx1;
  const gdouble a   = -dx1 / (ddx * dx2);
  const gdouble c   = +dx2 / (ddx * dx1);
  const gdouble b   = -(dx2 + dx1) / (dx1 * dx2);

  return a * fp2 + b * f + c * fp1;
}
*/

/*static complex double
_nc_hiqg_1d_cdiff2 (const complex double fp2, const complex double fp1, const complex double f, const gdouble dx2, const gdouble dx1)
{
  const gdouble ddx = dx2 - dx1;
  const gdouble a   = +2.0 / (ddx * dx2);
  const gdouble c   = -2.0 / (ddx * dx1);
  const gdouble b   = +2.0 / (dx1 * dx2);
  
  return a * fp2 + b * f + c * fp1;
}
*/
#if HAVE_SUNDIALS_MAJOR == 3

#define LOCAL_STENCIL 6

static gint 
_nc_hiqg_1d_f_expl (gdouble t, N_Vector y, N_Vector ydot, gpointer user_data) 
{
  NcHIQG1D *qg1d = NC_HIQG_1D (user_data);
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble *Yt   = N_VGetArrayPointer (y);
  gdouble *dYt  = N_VGetArrayPointer (ydot);
  const gint np = LOCAL_STENCIL;
  const gint sp = np / 2;
  gdouble x[LOCAL_STENCIL], S[LOCAL_STENCIL];
  NcmVector *xv      = ncm_vector_new_data_static (x, LOCAL_STENCIL, 1);
  NcmVector *Sv      = ncm_vector_new_data_static (S, LOCAL_STENCIL, 1);
  /*NcmSplineRBF *Srbf = ncm_spline_rbf_new (NCM_SPLINE_RBF_TYPE_GAUSS);*/
  /*NcmSpline *Ss      = NCM_SPLINE (Srbf);*/
  NcmSpline *Ss      = ncm_spline_gsl_new (gsl_interp_polynomial);
  /*NcmSpline *Ss      = ncm_spline_cubic_notaknot_new ();*/
  gint i;

  N_VConst (0.0, ydot);

  ncm_spline_set (Ss, xv, Sv, FALSE);

  for (i = 0; i < self->nknots; i++)
  {
    gdouble xi = Yt[_XI (i, self->nknots)];
    gdouble dSi;
    gint fi, j;

    if (i < sp)
      fi = 0;
    else if (i + sp + 1 > self->nknots)
      fi = self->nknots - np;
    else
      fi = i - sp;

    for (j = 0; j < np; j++)
    {
      gint k = fi + j;
      x[j]   = Yt[_XI (k, self->nknots)];
      S[j]   = Yt[_SI (k, self->nknots)];
    }

    ncm_spline_prepare (Ss);

    dSi = ncm_spline_eval_deriv (Ss, xi);

    dYt[_XI (i, self->nknots)] = 2.0 * dSi;
  }

  ncm_vector_free (xv);
  ncm_vector_free (Sv);
  ncm_spline_free (Ss);

  return 0;
}

static gint 
_nc_hiqg_1d_f_impl (gdouble t, N_Vector y, N_Vector ydot, gpointer user_data) 
{
  NcHIQG1D *qg1d = NC_HIQG_1D (user_data);
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble *Yt   = N_VGetArrayPointer (y);
  gdouble *dYt  = N_VGetArrayPointer (ydot);
  const gint np = LOCAL_STENCIL;
  const gint sp = np / 2;
  gdouble x[LOCAL_STENCIL], lnR[LOCAL_STENCIL], S[LOCAL_STENCIL];
  NcmVector *xv        = ncm_vector_new_data_static (x, LOCAL_STENCIL, 1);
  NcmVector *Sv        = ncm_vector_new_data_static (S, LOCAL_STENCIL, 1);
  NcmVector *lnRv      = ncm_vector_new_data_static (lnR, LOCAL_STENCIL, 1);
  /*NcmSplineRBF *Srbf   = ncm_spline_rbf_new (NCM_SPLINE_RBF_TYPE_GAUSS);*/
  /*NcmSplineRBF *lnRrbf = ncm_spline_rbf_new (NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS);*/
  /*NcmSpline *Ss        = NCM_SPLINE (Srbf);*/
  /*NcmSpline *lnRs      = NCM_SPLINE (lnRrbf);*/
  NcmSpline *Ss        = ncm_spline_gsl_new (gsl_interp_polynomial);
  NcmSpline *lnRs      = ncm_spline_gsl_new (gsl_interp_polynomial);
  /*NcmSpline *Ss        = ncm_spline_cubic_notaknot_new ();*/
  /*NcmSpline *lnRs      = ncm_spline_cubic_notaknot_new ();*/
  gint i;

  N_VConst (0.0, ydot);

  ncm_spline_set (Ss,   xv, Sv,   FALSE);
  ncm_spline_set (lnRs, xv, lnRv, FALSE);

  for (i = 0; i < self->nknots; i++)
  {
    const gdouble l0 = 1.0;
    gdouble xi       = Yt[_XI (i, self->nknots)];
    gdouble dSi, d2Si, dlnRi, d2lnRi;
    gint fi, j;


    if (i < sp)
      fi = 0;
    else if (i + sp + 1 > self->nknots)
      fi = self->nknots - np;
    else
      fi = i - sp;

    for (j = 0; j < np; j++)
    {
      gint k = fi + j;
      x[j]   = Yt[_XI   (k, self->nknots)];
      lnR[j] = Yt[_LNRI (k, self->nknots)] - l0 * log (x[j]);
      S[j]   = Yt[_SI   (k, self->nknots)];
    }

    ncm_spline_prepare (Ss);
    ncm_spline_prepare (lnRs);

    dlnRi  = ncm_spline_eval_deriv  (lnRs, xi) + l0 / xi;
    d2lnRi = ncm_spline_eval_deriv2 (lnRs, xi) - l0 / (xi * xi);
    dSi    = ncm_spline_eval_deriv  (Ss,   xi);
    d2Si   = ncm_spline_eval_deriv2 (Ss,   xi);

    dYt[_LNRI (i, self->nknots)] = - d2Si;
    dYt[_SI   (i, self->nknots)] = dSi * dSi - self->lambda / (xi * xi) + d2lnRi + dlnRi * dlnRi;

    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g\n", t, dSi * dSi, d2lnRi, dlnRi * dlnRi);*/
  }

  ncm_vector_free (xv);
  ncm_vector_free (lnRv);
  ncm_vector_free (Sv);

  ncm_spline_free (lnRs);
  ncm_spline_free (Ss);

  return 0;
}

/*
static gint 
_nc_hiqg_1d_J (N_Vector v, N_Vector Jv, gdouble t, N_Vector y, N_Vector fy, gpointer user_data, N_Vector tmp) 
{ 
  NcHIQG1D *qg1d = NC_HIQG_1D (user_data);
  NcHIQG1DPrivate * const self = qg1d->priv;
  complex double *V  = NULL;
  complex double *JV = NULL;
  sunindextype i;
  
  V  = (complex double *) N_VGetArrayPointer (v);
  JV = (complex double *) N_VGetArrayPointer (Jv);

  N_VConst (0.0, Jv);

  JV[0] = 0.0;
  for (i = 1; i < self->nknots - 1; i++)
  {
    const gdouble xi          = ncm_vector_fast_get (self->knots, i);
    const gdouble xim1        = ncm_vector_fast_get (self->knots, i - 1);
    const gdouble xip1        = ncm_vector_fast_get (self->knots, i + 1);
    const complex double d2Vi = _nc_hiqg_1d_cdiff2 (V[i + 1], V[i - 1], V[i], xip1 - xi, xim1 - xi); 

    JV[i] = -I * (- d2Vi + self->lambda * V[i] / (xi * xi));
  }
  i = self->nknots - 1;
  if (self->noboundary)
  {
    const gdouble xi          = ncm_vector_fast_get (self->knots, i);
    const gdouble xim1        = ncm_vector_fast_get (self->knots, i - 1);
    const complex double d2Vi = _nc_hiqg_1d_cdiff2 (0.0, V[i - 1], V[i], self->h, xim1 - xi); 

    JV[i] = -I * (- d2Vi + self->lambda * V[i] / (xi * xi));
  }
  else
    JV[i] = 0.0;
  
  return 0; 
}
*/
#endif /* HAVE_SUNDIALS_MAJOR == 3 */

static gdouble
_nc_hiqg_1d_calc_dS (gdouble x, NcHIQG1DPrivate * const self)
{
  const gdouble x0 = ncm_vector_fast_get (self->knots, 1);

  if (x < x0)
  {
    NcmSplineCubic *Re_sc = NCM_SPLINE_CUBIC (self->Re_psi_s);
    NcmSplineCubic *Im_sc = NCM_SPLINE_CUBIC (self->Im_psi_s);
    const gdouble br      = ncm_vector_get (Re_sc->b, 0);
    const gdouble bi      = ncm_vector_get (Im_sc->b, 0);
    const gdouble cr      = ncm_vector_get (Re_sc->c, 0);
    const gdouble ci      = ncm_vector_get (Im_sc->c, 0);
    const gdouble dr      = ncm_vector_get (Re_sc->d, 0);
    const gdouble di      = ncm_vector_get (Im_sc->d, 0);
    const gdouble rho_x2  = gsl_pow_2 (br + (cr + dr * x) * x) + gsl_pow_2 (bi + (ci + di * x) * x);
    const gdouble dS      = ((br * ci - bi * cr) + (2.0 * (br * di - bi * dr) + (cr * di - ci * dr) * x) * x) / rho_x2;

    return dS;
  }
  else
  {
    const gdouble Re_psi  = ncm_spline_eval (self->Re_psi_s, x);
    const gdouble Im_psi  = ncm_spline_eval (self->Im_psi_s, x);
    const gdouble Re_dpsi = ncm_spline_eval_deriv (self->Re_psi_s, x);
    const gdouble Im_dpsi = ncm_spline_eval_deriv (self->Im_psi_s, x);
    const gdouble rho     = gsl_pow_2 (Re_psi) + gsl_pow_2 (Im_psi); // ncm_spline_eval (self->rho_s, x); //
    const gdouble dS      = (Re_psi * Im_dpsi - Re_dpsi * Im_psi) / rho;

    return dS;
  }
}

static void
_nc_hiqg_1d_prepare_splines (NcHIQG1DPrivate * const self, const gdouble t, NcmVector *knots, gdouble *Y)
{
  if (!self->up_splines)
  {
    if (FALSE)
    {
    gdouble maxrho = GSL_NEGINF, minrho = GSL_POSINF, xf = 0.0;
    gdouble gamma, lna;
    gdouble Rlna, Ilna;
    gdouble Rs, Is;
    gint i;

    if (TRUE)
    {
      const complex double *Yc  = (const complex double *)Y;
      const gdouble x1          = ncm_vector_fast_get (knots, 1);
      const gdouble x2          = ncm_vector_fast_get (knots, 2);
      const gdouble lnx1        = log (x1);
      const gdouble lnx2        = log (x2);
      const complex double psi1 = Yc[1];
      const complex double psi2 = Yc[2];
      const gdouble R1          = cabs (psi1);
      const gdouble R2          = cabs (psi2);
      gdouble lnxa[2]           = {lnx1, lnx2};
      gdouble lnrhoa[2]         = {log (R1), log (R2)};
      gdouble dd[2], c[2] = {0.0, }, w[2];
      gdouble b1, cmp_R, cmp_I, cmp_B;
      /*gdouble Rgamma, Igamma;*/

      gsl_poly_dd_init (dd, lnxa, lnrhoa, 2);
      gsl_poly_dd_taylor (c, 0.0, dd, lnxa, 2, w);

      lna   = c[0];
      gamma = c[1];

      gsl_poly_dd_init (dd, lnxa, lnrhoa, 2);
      gsl_poly_dd_taylor (c, 0.0, dd, lnxa, 2, w);

      lnrhoa[0] = log (fabs (creal (psi1)));
      lnrhoa[1] = log (fabs (creal (psi2)));

      gsl_poly_dd_init (dd, lnxa, lnrhoa, 2);
      gsl_poly_dd_taylor (c, 0.0, dd, lnxa, 2, w);

      Rlna   = c[0];
      /*Rgamma = c[1];*/

      lnrhoa[0] = log (fabs (cimag (psi1)));
      lnrhoa[1] = log (fabs (cimag (psi2)));

      gsl_poly_dd_init (dd, lnxa, lnrhoa, 2);
      gsl_poly_dd_taylor (c, 0.0, dd, lnxa, 2, w);

      Ilna   = c[0];
      /*Igamma = c[1];*/

      Rs = GSL_SIGN (creal (psi1));
      Is = GSL_SIGN (cimag (psi2));

      b1 = log (hypot (exp (Rlna), exp (Ilna)));

      cmp_R = fabs (Rlna / lna - 1.0);
      cmp_I = fabs (Ilna / lna - 1.0);
      cmp_B = fabs (b1   / lna - 1.0);
/*
      printf ("%.6f | % 17.10g % 17.10g % 17.10g % 17.10g | % 17.10g % 17.10g % 17.10g | ", 
              t,
              lna,   Rlna,   Ilna, b1,
              cmp_R, cmp_I, cmp_B
              );
*/
      if (cmp_I < cmp_R)
      {
        if (ncm_cmp (cmp_I, cmp_B, 1.0e-3, 0.0) <= 0)
        {
          /*printf ("I dominates % 17.10g % 17.10g % 17.10g!\n", gamma, Rgamma, Igamma);*/
          ncm_vector_fast_set (self->Re_psi_v, 0, 0.0);
          ncm_vector_fast_set (self->Im_psi_v, 0, Is * exp (Ilna));
        }
        else
        {
          /*printf ("B dominates % 17.10g % 17.10g % 17.10g!\n", gamma, Rgamma, Igamma);*/
          ncm_vector_fast_set (self->Re_psi_v, 0, Rs * exp (Rlna));
          ncm_vector_fast_set (self->Im_psi_v, 0, Is * exp (Ilna));
        }
      }
      else
      {
        if (ncm_cmp (cmp_R, cmp_B, 1.0e-3, 0.0) <= 0)
        {
          /*printf ("R dominates % 17.10g % 17.10g % 17.10g!\n", gamma, Rgamma, Igamma);*/
          ncm_vector_fast_set (self->Re_psi_v, 0, Rs * exp (Rlna));
          ncm_vector_fast_set (self->Im_psi_v, 0, 0.0);
        }
        else
        {
          /*printf ("B dominates % 17.10g % 17.10g % 17.10g!\n", gamma, Rgamma, Igamma);*/
          ncm_vector_fast_set (self->Re_psi_v, 0, Rs * exp (Rlna));
          ncm_vector_fast_set (self->Im_psi_v, 0, Is * exp (Ilna));
        }
      }
    }

    if (gamma <= 0.0 || TRUE)
    {
      self->gamma = gamma = 0.0;
      ncm_vector_fast_set (self->Re_psi_v, 0, 0.0);
      ncm_vector_fast_set (self->Im_psi_v, 0, 0.0);
      ncm_vector_fast_set (self->rho_v,    0, 0.0);
    }
    else
    {
      self->gamma = gamma;
      ncm_vector_fast_set (self->rho_v, 0, exp (2.0 * lna));
    }

    for (i = 1; i < self->nknots; i++)
    {
      const gdouble x_i   = ncm_vector_fast_get (knots, i);
      const gdouble b     = hypot (Y[2 * i], Y[2 * i + 1]);
      const gdouble rho_i = b * b;
      const gdouble f1    = pow (x_i, -gamma);
      
      ncm_vector_fast_set (self->Re_psi_v, i, Y[2 * i]     * f1);
      ncm_vector_fast_set (self->Im_psi_v, i, Y[2 * i + 1] * f1);
      ncm_vector_fast_set (self->rho_v,    i, rho_i        * f1 * f1);

      /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x_i, Y[2 * i] * f1, Y[2 * i + 1] * f1, rho_i * f1 * f1, f1);*/
      
      minrho = MIN (minrho, rho_i);
      maxrho = MAX (maxrho, rho_i);
    }

    /*
    ncm_spline_prepare (self->Re_psi_s);
    ncm_spline_prepare (self->Im_psi_s);
    ncm_spline_prepare (self->rho_s);
    */

    ncm_spline_set (self->Re_psi_s, knots, self->Re_psi_v, TRUE);
    ncm_spline_set (self->Im_psi_s, knots, self->Im_psi_v, TRUE);
    ncm_spline_set (self->rho_s,    knots, self->rho_v,    TRUE);
    
    {
      gsl_function F;

      /*printf ("A minrho % 22.15g % 22.15g\n", minrho, maxrho);*/

      if (minrho == 0.0)
        minrho = maxrho * self->abstol;
      else
        minrho *= 1.0e4;

      /*printf ("D minrho % 22.15g % 22.15g\n", minrho, maxrho);*/
      
      if (TRUE)
      {
        gdouble last_rho = 0.0;
        for (i = 1; i < self->nknots; i++)
        {
          const gdouble x_i   = ncm_vector_fast_get (self->spec_knots, i);
          const gdouble rho_i = ncm_vector_fast_get (self->rho_v, i) * pow (x_i, 2.0 * gamma);
          
          /*printf ("# x_i = % 22.15g, MIN % 22.15g % 22.15g\n", x_i, minrho, rho_i);*/

          if ((rho_i < last_rho) && (rho_i < minrho))
          {
            xf = x_i;
            break;
          }
          last_rho = rho_i;
        }
      }
      printf ("# xf = % 22.15g, gamma % 22.15g\n", xf, gamma);
      
      F.params   = self;
      F.function = (gdouble (*) (gdouble, gpointer)) &_nc_hiqg_1d_calc_dS;

      ncm_spline_set_func (self->dS_s, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, xf, 10000, self->reltol);
    }
    }
    else
    {
      NcmVector *x   = ncm_vector_new_data_static (&Y[_XI   (0, self->nknots)], self->nknots, _XI_STRIDE);
      NcmVector *lnR = ncm_vector_new_data_static (&Y[_LNRI (0, self->nknots)], self->nknots, _LNRI_STRIDE);
      NcmVector *S   = ncm_vector_new_data_static (&Y[_SI   (0, self->nknots)], self->nknots, _SI_STRIDE);
 
      ncm_spline_set (self->dS_s,     x, S,   TRUE);
      ncm_spline_set (self->rho_s,    x, lnR, TRUE);
      ncm_spline_set (self->Re_psi_s, x, S,   TRUE);
      ncm_spline_set (self->Im_psi_s, x, S,   TRUE);

      ncm_spline_prepare (self->dS_s);
      ncm_spline_prepare (self->rho_s);
      ncm_spline_prepare (self->Re_psi_s);
      ncm_spline_prepare (self->Im_psi_s);      

      if (FALSE)
      {
        gint i;
        for (i = 0; i < self->nknots; i++)
        {
          printf ("T % 22.15g % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g\n", t, 
                  Y[_XI (i, self->nknots)], Y[_LNRI (i, self->nknots)], Y[_SI (i, self->nknots)],
                  ncm_vector_get (x, i), ncm_vector_get (lnR, i), ncm_vector_get (S, i),
                  ncm_spline_eval (self->rho_s, ncm_vector_get (x, i)), ncm_spline_eval (self->dS_s, ncm_vector_get (x, i)));
        }
      }
      
      ncm_vector_free (x);
      ncm_vector_free (lnR);
      ncm_vector_free (S);
    }
    self->up_splines = TRUE;
  }
}

/**
 * nc_hiqg_1d_evolve:
 * @qg1d: a #NcHIQG1D
 * @tf: final time
 * 
 * Evolve the wave-function to @tf.
 * 
 */
void
nc_hiqg_1d_evolve (NcHIQG1D *qg1d, const gdouble tf)
{
#if HAVE_SUNDIALS_MAJOR == 3
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble t  = self->ti;
  gdouble *Y = N_VGetArrayPointer (self->y);
  gint flag;

  flag = ARKodeSetStopTime (self->arkode, tf);
  NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

  /*printf ("# Evolving! % 22.15g => % 22.15g\n", t, tf);*/
#if 0
  glong iout = 0;
  glong nni, nni_cur = 0; 
  glong nli, nli_cur = 0;
  glong nncfails, nncfails_cur = 0;
  glong njvevals, njvevals_cur = 0;
  gdouble olddt, newdt;  
#endif
  while (t < tf) 
  {
    flag = ARKode (self->arkode, tf, self->y, &t, ARK_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "ARKode", 1, );

    Y = N_VGetArrayPointer (self->y);

    /*_nc_hiqg_1d_prepare_splines (self, t, Y);*/
    /*printf ("# STEP % 22.15g\n", t);*/
/*
    printf ("# STEP % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g\n",
            t,
            Y[_XI   (0, self->nknots)], Y[_XI   (1, self->nknots)], Y[_XI   (2, self->nknots)],
            Y[_LNRI (0, self->nknots)], Y[_LNRI (1, self->nknots)], Y[_LNRI (2, self->nknots)],
            Y[_SI   (0, self->nknots)], Y[_SI   (1, self->nknots)], Y[_SI   (2, self->nknots)]
            );
*/

#if 0
    flag = ARKodeGetLastStep (self->arkode, &olddt);
    NCM_CVODE_CHECK (&flag, "ARKodeGetLastStep", 1, );
    
    flag = ARKodeGetCurrentStep (self->arkode, &newdt);
    NCM_CVODE_CHECK (&flag, "ARKodeGetCurrentStep", 1, );
    
    flag = ARKodeGetNumNonlinSolvIters (self->arkode, &nni);
    NCM_CVODE_CHECK (&flag, "ARKodeGetNumNonlinSolvIters", 1, );
    
    flag = ARKSpilsGetNumLinIters (self->arkode, &nli);
    NCM_CVODE_CHECK (&flag, "ARKSpilsGetNumLinIters", 1, );

    flag = ARKodeGetNumNonlinSolvConvFails (self->arkode, &nncfails);
    NCM_CVODE_CHECK (&flag, "ARKodeGetNumNonlinSolvConvFails", 1, );

    flag = ARKodeGetNumNonlinSolvConvFails (self->arkode, &nncfails);
    NCM_CVODE_CHECK (&flag, "ARKodeGetNumNonlinSolvConvFails", 1, );

    flag = ARKSpilsGetNumJtimesEvals (self->arkode, &njvevals);
    NCM_CVODE_CHECK (&flag, "ARKSpilsGetNumJtimesEvals", 1, );

    // print current solution stats 
    iout++;
    printf(" %4ld  %22.15g  %22.15g nni %3ld nli %3ld nncfails %3ld njvevals %3ld\n", 
           iout, olddt, newdt, 
           nni - nni_cur, 
           nli - nli_cur, 
           nncfails - nncfails_cur, 
           njvevals - njvevals_cur);

    nni_cur      = nni;
    nli_cur      = nli;
    nncfails_cur = nncfails;
    njvevals_cur = njvevals;
#endif
#if 0
    {
      gint i;
      printf ("% 22.15g", t);
      for (i = 0; i < self->nknots; i++)
      {
        printf (" % 22.15g % 22.15g", Y[2 * i], Y[2 * i + 1]);
      }
      printf ("\n");
    }
#endif

  }

  self->ti = t;
  self->up_splines = FALSE;
  /* _nc_hiqg_1d_prepare_splines (self, self->ti, self->knots, Y); */
#endif /* HAVE_SUNDIALS_MAJOR == 3 */
}

/**
 * nc_hiqg_1d_evolve_spec:
 * @qg1d: a #NcHIQG1D
 * @t: final time
 * 
 * Evolves the system using spectral methods.
 * 
 */
void
nc_hiqg_1d_evolve_spec (NcHIQG1D *qg1d, const gdouble t)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  complex double *Y   = (complex double *) ncm_vector_data (self->spec_Y);
  const gdouble norma = gsl_pow_2 (gsl_dht_k_sample (self->dht, 0) / gsl_dht_x_sample (self->dht, 0));
  gint i;

  for (i = 0; i < self->nknots - 1; i++)
  {
    const gdouble li  = gsl_dht_k_sample (self->dht, i);
    complex double W  = ncm_vector_fast_get (self->specReW, i) + I * ncm_vector_fast_get (self->specImW, i);
    complex double Wt = cexp (-I * li * li * t) * W;

    ncm_vector_fast_set (self->spec_ReY, i, creal (Wt));
    ncm_vector_fast_set (self->spec_ImY, i, cimag (Wt));
  }

  gsl_dht_apply (self->dht, ncm_vector_data (self->spec_ReY), ncm_vector_data (self->spec_ReYt));
  gsl_dht_apply (self->dht, ncm_vector_data (self->spec_ImY), ncm_vector_data (self->spec_ImYt));

  Y[0] = 0.0;
  for (i = 0; i < self->nknots - 1; i++)
  {
    const gdouble x   = gsl_dht_x_sample (self->dht, i);    
    complex double Yt = norma * sqrt (x) * (ncm_vector_fast_get (self->spec_ReYt, i) + I * ncm_vector_fast_get (self->spec_ImYt, i));

    Y[i + 1] = Yt;
  }

  self->ti = t;
  self->up_splines = FALSE;

  /* _nc_hiqg_1d_prepare_splines (self, self->ti, self->spec_knots, ncm_vector_data (self->spec_Y)); */
}

/**
 * nc_hiqg_1d_eval_psi:
 * @qg1d: a #NcHIQG1D
 * @x: (array length=len) (element-type gdouble): array of points
 * @len: array length
 * 
 * Evaluates $\psi$ in the array @x.
 * 
 * Returns: (transfer full) (array) (element-type gdouble): array of length 2*@len containing the real and imaginary parts of $\psi$
 */
GArray *
nc_hiqg_1d_eval_psi (NcHIQG1D *qg1d, const gdouble *x, const guint len)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  GArray *psi = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gint i;

  g_array_set_size (psi, 2 * len);
  
  for (i = 0; i < len; i++)
  {
    const gdouble lnR_i        = ncm_spline_eval (self->rho_s, x[i]);
    const gdouble S_i          = ncm_spline_eval (self->dS_s, x[i]);
    const complex double psi_i = cexp (lnR_i + I * S_i);
    /*const gdouble Re_psi_i = ncm_spline_eval (self->Re_psi_s, x[i]);*/
    /*const gdouble Im_psi_i = ncm_spline_eval (self->Im_psi_s, x[i]);*/
    /*const gdouble f1       = pow (x[i], self->gamma);*/

    g_array_index (psi, gdouble, 2 * i + 0) = creal (psi_i);
    g_array_index (psi, gdouble, 2 * i + 1) = cimag (psi_i);
  }

  return psi;
}

/**
 * nc_hiqg_1d_eval_rho:
 * @qg1d: a #NcHIQG1D
 * @x: (array length=len) (element-type gdouble): array of points
 * @len: array length
 * 
 * Evaluates $\rho$ in the array @x.
 * 
 * Returns: (transfer full) (array) (element-type gdouble): array of length @len containing $\rho$
 */
GArray *
nc_hiqg_1d_eval_rho (NcHIQG1D *qg1d, const gdouble *x, const guint len)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  GArray *rho = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gint i;

  g_array_set_size (rho, len);

  for (i = 0; i < len; i++)
  {
    /*const gdouble rho_i = ncm_spline_eval (self->rho_s, x[i]);*/
    const gdouble lnR_i = ncm_spline_eval (self->rho_s, x[i]);
    /*const gdouble f1    = pow (x[i], 2.0 * self->gamma);*/

    /*printf ("rho(% 22.15g) = % 22.15g, % 22.15g\n", x[i], rho_i, sqrt (rho_i));*/
    g_array_index (rho, gdouble, i) = exp (2.0 * lnR_i);
  }

  return rho;
}

/**
 * nc_hiqg_1d_eval_dS:
 * @qg1d: a #NcHIQG1D
 * @x: (array length=len) (element-type gdouble): array of points
 * @len: array length
 * 
 * Evaluates $\partial_x S$ in the array @x.
 * 
 * Returns: (transfer full) (array) (element-type gdouble): array of length @len containing $\partial_x S$
 */
GArray *
nc_hiqg_1d_eval_dS (NcHIQG1D *qg1d, const gdouble *x, const guint len)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  GArray *dS = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gint i;

  g_array_set_size (dS, len);

  for (i = 0; i < len; i++)
  {
    const gdouble dS_i  = ncm_spline_eval_deriv (self->dS_s, x[i]);
    g_array_index (dS, gdouble, i) = dS_i;
  }

  return dS;
}

/**
 * nc_hiqg_1d_peek_rho_s:
 * @qg1d: a #NcHIQG1D
 * 
 * Peeks the current $\rho(x)$ #NcmSpline.
 * 
 * Returns: (transfer none): $\rho(x)$ #NcmSpline
 */
NcmSpline *
nc_hiqg_1d_peek_rho_s (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->rho_s;
}

static gdouble
_nc_hiqg_1d_eval_rho (const gdouble x, gpointer user_data)
{
  NcHIQG1DPrivate * const self = user_data;
  const gdouble lnR_i = ncm_spline_eval (self->rho_s, x);

  return exp (2.0 * lnR_i);
}

static gdouble
_nc_hiqg_1d_eval_xrho (const gdouble x, gpointer user_data)
{
  NcHIQG1DPrivate * const self = user_data;
  const gdouble lnR_i = ncm_spline_eval (self->rho_s, x);

  return x * exp (2.0 * lnR_i);
}

/**
 * nc_hiqg_1d_eval_int_rho:
 * @qg1d: a #NcHIQG1D
 * 
 * Computes $\int\mathrm{d}x\rho(x)$.
 * 
 * Returns: $\int\mathrm{d}x\rho(x)$.
 */
gdouble
nc_hiqg_1d_eval_int_rho (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gsl_function F;
  gdouble int_rho = 0.0;
  gdouble int_err = 0.0;
  NcmVector *x = ncm_spline_get_xv (self->rho_s);

  F.function = &_nc_hiqg_1d_eval_rho;
  F.params   = self;

  gsl_integration_qag (&F, 
                       ncm_vector_get (x, 0), 
                       ncm_vector_get (x, self->nknots - 1), 
                       0.0, self->reltol, 
                       NCM_INTEGRAL_PARTITION, 6, *w, &int_rho, &int_err);

  ncm_memory_pool_return (w);
  ncm_vector_free (x);
  
  return int_rho;
}

/**
 * nc_hiqg_1d_eval_int_xrho:
 * @qg1d: a #NcHIQG1D
 * 
 * Computes $\int\mathrm{d}x\;x\rho(x)$.
 * 
 * Returns: $\int\mathrm{d}x\;x\rho(x)$.
 */
gdouble
nc_hiqg_1d_eval_int_xrho (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gsl_function F;
  gdouble int_rho = 0.0;
  gdouble int_err = 0.0;
  NcmVector *x = ncm_spline_get_xv (self->rho_s);

  F.function = &_nc_hiqg_1d_eval_xrho;
  F.params   = self;

  gsl_integration_qag (&F, 
                       ncm_vector_get (x, 0), 
                       ncm_vector_get (x, self->nknots - 1), 
                       0.0, self->reltol, 
                       NCM_INTEGRAL_PARTITION, 6, *w, &int_rho, &int_err);

  ncm_memory_pool_return (w);
  ncm_vector_free (x);
  
  return int_rho;
}

gdouble
_nc_hiqg_1d_spec_basis (const gdouble x, const gdouble y, const gdouble h, const gdouble a)
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
 * nc_hiqg_1d_spec_basis:
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
nc_hiqg_1d_spec_basis (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  return _nc_hiqg_1d_spec_basis (x, y, h, a);
}

static gdouble
_nc_hiqg_1d_spec_Hbasis (const gdouble x, const gdouble y, const gdouble h, const gdouble a)
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
 * nc_hiqg_1d_spec_Hbasis:
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
nc_hiqg_1d_spec_Hbasis (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble h, const gdouble a)
{
  return _nc_hiqg_1d_spec_Hbasis (x, y, h, a);
}

static gdouble
_nc_hiqg_1d_spec_basis_xxa (const gdouble x, const gdouble y, const gdouble h, const gdouble a)
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

static gdouble
_nc_hiqg_1d_spec_Sbasis_x3 (const gdouble x, const gdouble y1, const gdouble y2, const gdouble h, const gdouble a)
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
 * nc_hiqg_1d_spec_Sbasis_x3:
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
nc_hiqg_1d_spec_Sbasis_x3 (NcHIQG1D *qg1d, const gdouble x, const gdouble y1, const gdouble y2, const gdouble h, const gdouble a)
{
  return _nc_hiqg_1d_spec_Sbasis_x3 (x, y1, y2, h, a);
}

/**
 * nc_hiqg_1d_spec_get_lambda:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_get_lambda (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->lambda;
}

/**
 * nc_hiqg_1d_spec_get_basis_a:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_get_basis_a (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->basis_a;
}

/**
 * nc_hiqg_1d_spec_get_acs_a:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_get_acs_a (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->acs_a;
}

/**
 * nc_hiqg_1d_spec_get_nu:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_get_nu (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->nu;
}

static gdouble
_nc_hiqg_1d_spec_If2 (const gdouble y1, const gdouble y2, const gdouble h)
{
  const gdouble h2      = h * h;
  const gdouble y12     = y1 * y1;
  const gdouble y22     = y2 * y2;
  const gdouble b1      = 0.5 * h2 * y1 * y2;
  const gdouble a1      = 0.25 * h2 * (y12 + y22);

  if (b1 < 1.0)
  {
    return exp (-a1) * sinh (b1) * 2.0 * ncm_c_sqrt_pi () / h;
  }
  else
  {
    return exp (-a1 + gsl_sf_lnsinh (b1)) * 2.0 * ncm_c_sqrt_pi () / h;
  }
}

static gdouble
_nc_hiqg_1d_spec_Ixf2 (const gdouble y1, const gdouble y2, const gdouble h)
{
  const gdouble h2     = h * h;
  const gdouble ym12   = 0.5 * h * (y1 - y2);
  const gdouble yp12   = 0.5 * h * (y1 + y2);
  const gdouble ym12_2 = ym12 * ym12;
  const gdouble yp12_2 = yp12 * yp12;

  return (exp (-ym12_2) * yp12 * erf (yp12) - exp (-yp12_2) * ym12 * erf (ym12)) * ncm_c_sqrt_pi () / h2;
}

static void _nc_hiqg_1d_spec_evol_C (NcHIQG1D *qg1d, const gdouble t);

static gint
_nc_hiqg_1d_spec_bohm_f (gdouble t, N_Vector y, N_Vector ydot, gpointer user_data)
{
  NcHIQG1D *qg1d = NC_HIQG_1D (user_data);
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble *Y  = N_VGetArrayPointer (y);
  gdouble *dY = N_VGetArrayPointer (ydot);
  gint i;

  _nc_hiqg_1d_spec_evol_C (qg1d, t);

  for (i = 0; i < self->nBohm; i++)
  {
    dY[i] = 2.0 * nc_hiqg_1d_spec_eval_dS (qg1d, Y[i]);
  }

  return 0;
}

void
_nc_hiqg_1d_spec_init_solver (NcHIQG1D *qg1d)
{
#if HAVE_SUNDIALS_MAJOR == 3
  NcHIQG1DPrivate * const self = qg1d->priv;
  const gdouble t0 = 0.0;
  gint flag;

  if (self->bohm != NULL)
    ARKodeFree (&self->bohm);

  self->nBohm = 10;

  g_clear_pointer (&self->y, N_VDestroy);
  self->yBohm = N_VNew_Serial (self->nBohm);

  NV_Ith_S (self->yBohm, 0) = 0.01;
  NV_Ith_S (self->yBohm, 1) = 0.21;
  NV_Ith_S (self->yBohm, 2) = 0.41;
  NV_Ith_S (self->yBohm, 3) = 0.61;
  NV_Ith_S (self->yBohm, 4) = 0.81;
  NV_Ith_S (self->yBohm, 5) = 1.01;
  NV_Ith_S (self->yBohm, 6) = 1.21;
  NV_Ith_S (self->yBohm, 7) = 1.41;
  NV_Ith_S (self->yBohm, 8) = 1.61;
  NV_Ith_S (self->yBohm, 9) = 1.81;

  self->bohm = ARKodeCreate ();
  NCM_CVODE_CHECK (&self->bohm, "ARKodeCreate", 0, );

  flag = ARKodeInit (self->bohm, &_nc_hiqg_1d_spec_bohm_f, NULL, t0, self->yBohm);
  NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );

  flag = ARKodeSetUserData (self->bohm, (void *) qg1d);
  NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

  flag = ARKodeSetMaxNumSteps (self->bohm, 10000);
  NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

  flag = ARKodeSStolerances (self->bohm, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );

#endif /* HAVE_SUNDIALS_MAJOR == 3 */
}

/**
 * nc_hiqg_1d_spec_prepare:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_spec_prepare (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i, ret;

  /* WARNING self->h must not be a integer multiple of y-x */
  self->h = 1.0 / (1.99 * (self->xf - self->xi) / (1.0 * self->nknots));
  ncm_matrix_set_zero (self->IM);

	for (i = 0; i < self->nknots; i++)
	{
		const gdouble x1 = ncm_vector_get (self->knots, i);
    const gdouble f1 = 2.49;
    if (i == 0)
    {
      const gdouble si = x1;
      const gdouble hi = 1.0 / (f1 * si);

      ncm_vector_set (self->hv, i, hi);
    }
    else if (i + 1 == self->nknots)
    {
      const gdouble x1m1 = ncm_vector_get (self->knots, i - 1);
      const gdouble si   = x1 - x1m1;
      const gdouble hi   = 1.0 / (f1 * si);

      ncm_vector_set (self->hv, i, hi);
    }
    else
    {
      const gdouble x1m1 = ncm_vector_get (self->knots, i - 1);
      const gdouble x1p1 = ncm_vector_get (self->knots, i + 1);
      const gdouble si   = MAX (x1 - x1m1, x1p1 - x1);
      const gdouble hi   = 1.0 / (f1 * si);

      ncm_vector_set (self->hv, i, hi);
    }
    ncm_vector_set (self->hv, i, self->h);
  }

	for (i = 0; i < self->nknots; i++)
	{
		const gdouble x1    = ncm_vector_get (self->knots, i);
    const gdouble hi    = ncm_vector_get (self->hv, i);
    const gdouble Ix1x1 = _nc_hiqg_1d_spec_basis (x1, x1, hi, self->basis_a);
    const gdouble Kx1x1 = /*Ix1x1;//*/_nc_hiqg_1d_spec_Hbasis (x1, x1, hi, self->basis_a);
    guint j;

    ncm_matrix_set (self->IM, i, i, Ix1x1);
    ncm_matrix_set (self->KM, i, i, Kx1x1);
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g\n", h, x1, x1, Ix1x1);*/
		for (j = i + 1; j < self->nknots; j++)
		{
			const gdouble x2    = ncm_vector_get (self->knots, j);
      const gdouble hj    = ncm_vector_get (self->hv, j);
			const gdouble Ix1x2 = _nc_hiqg_1d_spec_basis (x1, x2, hj, self->basis_a);
			const gdouble Kx1x2 = /*Ix1x2;//*/_nc_hiqg_1d_spec_Hbasis (x1, x2, hj, self->basis_a);
			const gdouble Kx2x1 = /*Ix1x2;//*/_nc_hiqg_1d_spec_Hbasis (x2, x1, hj, self->basis_a);

			ncm_matrix_set (self->IM, j, i, Ix1x2);
      ncm_matrix_set (self->KM, i, j, Kx2x1);
      ncm_matrix_set (self->KM, j, i, Kx1x2);
    }
	}
  /*ncm_vector_log_vals (self->knots, "#    KNOTS:  ", "% .5e", TRUE);*/
  /*ncm_vector_log_vals (self->hv,    "#       HV:  ", "% .5e", TRUE);*/
  printf ("#        H:  % 22.15g\n", self->h);

  {
    gint lwork  = -1;
    gchar equed = 'Y';

    g_array_set_size (self->ipiv, self->nknots);
    g_array_set_size (self->iwork, 2 * self->nknots);

    if (self->work->len < 1)
      g_array_set_size (self->work, 1);

    /*ncm_matrix_log_vals (self->IM, "#     IM: ", "% .5e");*/
    /*ncm_matrix_log_vals (self->KM, "#     KM: ", "% .5e");*/

    if (FALSE)
    {
      ret = ncm_lapack_dsytrf ('L', self->nknots,
                               ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                               (gint *)self->ipiv->data, (gdouble *)self->work->data, lwork);
      g_assert_cmpint (ret, ==, 0);

      lwork = g_array_index (self->work, gdouble, 0);
      if (lwork > self->work->len)
        g_array_set_size (self->work, lwork);

      ret = ncm_lapack_dsytrf ('L', self->nknots,
                               ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                               (gint *)self->ipiv->data, (gdouble *)self->work->data, lwork);
      g_assert_cmpint (ret, ==, 0);

      ncm_matrix_memcpy (self->IM_fact, self->IM);
      /*ncm_matrix_log_vals (self->IM, "# DEC IM: ", "% .5e");*/

      ret = ncm_lapack_dsytrs ('L', self->nknots, self->nknots,
                               ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                               (gint *)self->ipiv->data,
                               ncm_matrix_data (self->KM), ncm_matrix_tda (self->KM));
      g_assert_cmpint (ret, ==, 0);

      ncm_matrix_memcpy (self->JM, self->KM);
      /*ncm_matrix_log_vals (self->JM, "#     JM: ", "% .5e");*/
    }
    else
    {
      gdouble rcond, rpvgrw;

      ncm_matrix_set_zero (self->err_norm);
      ncm_matrix_set_zero (self->err_comp);

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

      printf ("#    RCOND:  % 22.15g\n", rcond);
      printf ("#   RPVGRW:  % 22.15g\n", rpvgrw);

      /*ncm_vector_log_vals (self->scaling,  "#  SCALING:  ", "% .5e", TRUE);*/
      /*ncm_vector_log_vals (self->berr,     "#     BERR:  ", "% .5e", TRUE);*/
      /*ncm_matrix_log_vals (self->err_norm, "# ERR_NORM: ",  "% .5e");*/
      /*ncm_matrix_log_vals (self->err_comp, "# ERR_COMP: ",  "% .5e");*/
      /*ncm_matrix_log_vals (self->JM,       "#       JM: ",  "% .5e");*/

      g_assert_cmpint (ret, ==, 0);
    }

    if (FALSE)
    {
      lwork = -1;
      ret = ncm_lapack_dgeev ('V', 'V', self->nknots,
                              ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                              ncm_vector_data (self->wr),
                              ncm_vector_data (self->wi),
                              ncm_matrix_data (self->vl), ncm_matrix_tda (self->vl),
                              ncm_matrix_data (self->vr), ncm_matrix_tda (self->vr),
                              (gdouble *)self->work->data, lwork);
      g_assert_cmpint (ret, ==, 0);

      lwork = g_array_index (self->work, gdouble, 0);
      if (lwork > self->work->len)
        g_array_set_size (self->work, lwork);

      ret = ncm_lapack_dgeev ('V', 'V', self->nknots,
                              ncm_matrix_data (self->JM), ncm_matrix_tda (self->JM),
                              ncm_vector_data (self->wr), ncm_vector_data (self->wi),
                              ncm_matrix_data (self->vl), ncm_matrix_tda (self->vl),
                              ncm_matrix_data (self->vr), ncm_matrix_tda (self->vr),
                              (gdouble *)self->work->data, lwork);
      g_assert_cmpint (ret, ==, 0);
    }
    else
    {
      gint ilo, ihi;
      gdouble abnrm;

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
    }

    /*ncm_vector_log_vals (self->wr, "#       WR:  ", "% .5e", TRUE);*/
    /*ncm_vector_log_vals (self->wi, "#       WI:  ", "% .5e", TRUE);*/

    /*ncm_matrix_log_vals (self->vl, "#       VL: ", "% .5e");*/
    /*ncm_matrix_log_vals (self->vr, "#       VR: ", "% .5e");*/

    ncm_vector_memcpy (self->spec_ReC0, self->spec_ReY);
    ncm_vector_memcpy (self->spec_ImC0, self->spec_ImY);

    memcpy (ncm_vector_data (self->spec_Y), ncm_matrix_data (self->spec_C0), sizeof (gdouble) * self->nknots * 2);

    /*ncm_vector_log_vals (self->spec_ReC0, "# ReY:  ", "% .5e", TRUE);*/
    /*ncm_vector_log_vals (self->spec_ImC0, "# ImY:  ", "% .5e", TRUE);*/

    if (FALSE)
    {
      ret = ncm_lapack_dsytrs ('L', self->nknots, 2,
                               ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                               (gint *)self->ipiv->data,
                               ncm_matrix_data (self->spec_C0), ncm_matrix_tda (self->spec_C0));
      g_assert_cmpint (ret, ==, 0);
    }
    else
    {
      gdouble rcond, rpvgrw;

      ncm_matrix_set_zero (self->err_norm);
      ncm_matrix_set_zero (self->err_comp);

      ret = ncm_lapack_dsysvxx ('F', 'L', self->nknots, 2,
                                ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM),
                                ncm_matrix_data (self->IM_fact), ncm_matrix_tda (self->IM_fact),
                                (gint *)self->ipiv->data, &equed,
                                ncm_vector_data (self->scaling),
                                ncm_vector_data (self->spec_Y), ncm_matrix_tda (self->spec_C0),
                                ncm_matrix_data (self->spec_C0), ncm_matrix_tda (self->spec_C0),
                                &rcond, &rpvgrw,
                                ncm_vector_data (self->berr),
                                3, ncm_matrix_data (self->err_norm), ncm_matrix_data (self->err_comp), 0, NULL,
                                (gdouble *)self->work->data, (gint *)self->iwork->data);

      printf ("#    RCOND:  % 22.15g\n", rcond);
      printf ("#   RPVGRW:  % 22.15g\n", rpvgrw);

      /*ncm_vector_log_vals (self->scaling,  "#  SCALING:",  "% .5e", TRUE);*/
      /*ncm_vector_log_vals (self->berr,     "#     BERR:",  "% .5e", TRUE);*/
      /*ncm_matrix_log_vals (self->err_norm, "# ERR_NORM: ", "% .5e");*/
      /*ncm_matrix_log_vals (self->err_comp, "# ERR_COMP: ", "% .5e");*/

      g_assert_cmpint (ret, ==, 0);
    }

    /*ncm_matrix_log_vals (self->spec_C0, "#     C0: ", "% .5e");*/

    for (i = 0; i < self->nknots; i++)
    {
      NcmVector *vr_i      = ncm_matrix_get_row (self->vr, i);
      NcmVector *vl_i      = ncm_matrix_get_row (self->vl, i);
      const gdouble n_i    = ncm_vector_dot (vl_i, vr_i);
      const gdouble ReC0_i = ncm_matrix_get (self->spec_C0, 0, i);
      const gdouble ImC0_i = ncm_matrix_get (self->spec_C0, 1, i);
      const gdouble ReA0_i = ncm_vector_dot (vl_i, self->spec_ReC0) / n_i;
      const gdouble ImA0_i = ncm_vector_dot (vl_i, self->spec_ImC0) / n_i;

      ncm_vector_set (self->vlvr, i, n_i);

      ncm_vector_set (self->spec_A0, 2 * i + 0, ReA0_i);
      ncm_vector_set (self->spec_A0, 2 * i + 1, ImA0_i);

      ncm_vector_set (self->spec_C,  2 * i + 0, ReC0_i);
      ncm_vector_set (self->spec_C,  2 * i + 1, ImC0_i);

      /*printf ("#      DOT:  % 22.15g % 22.15g % 22.15g\n", ncm_vector_dot (vl_i, vl_i), ncm_vector_dot (vr_i, vr_i), ncm_vector_dot (vl_i, vr_i));*/

      ncm_vector_free (vr_i);
      ncm_vector_free (vl_i);
    }
  }

  _nc_hiqg_1d_spec_init_solver (qg1d);
}

/**
 * nc_hiqg_1d_spec_peek_knots:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmVector *
nc_hiqg_1d_spec_peek_knots (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return self->knots;
}

/**
 * nc_hiqg_1d_spec_eval_ev:
 * @qg1d: a #NcHIQG1D
 * @i: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_eval_ev (NcHIQG1D *qg1d, const gint i, const gdouble x)
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

      res += cj * _nc_hiqg_1d_spec_basis (x, xj, self->h, self->basis_a);
    }

    ncm_vector_free (vr_i);

    return res;
  }
}

/**
 * nc_hiqg_1d_spec_eval_psi0:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @psi0: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi_0$
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_spec_eval_psi0 (NcHIQG1D *qg1d, const gdouble x, gdouble *psi0)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i;

  psi0[0] = 0.0;
  psi0[1] = 0.0;
  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi     = ncm_vector_get (self->knots, i);
    const gdouble ReC0_i = ncm_vector_get (self->spec_ReC0, i);
    const gdouble ImC0_i = ncm_vector_get (self->spec_ImC0, i);
    const gdouble Kxxi   = _nc_hiqg_1d_spec_basis (x, xi, self->h, self->basis_a);

    psi0[0] += ReC0_i * Kxxi;
    psi0[1] += ImC0_i * Kxxi;
  }
}

static void
_nc_hiqg_1d_spec_evol_C (NcHIQG1D *qg1d, const gdouble t)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  complex double *A0 = (complex double *) ncm_vector_data (self->spec_A0);
  gint i;

  ncm_vector_set_zero (self->spec_C);
  for (i = 0; i < self->nknots; i++)
  {
    const gdouble wr = ncm_vector_get (self->wr, i);
    complex double A = A0[i] * cexp (-I * wr * t);
    cblas_daxpy (self->nknots, creal (A), ncm_matrix_ptr (self->vr, i, 0), 1, ncm_vector_ptr (self->spec_C, 0), 2);
    cblas_daxpy (self->nknots, cimag (A), ncm_matrix_ptr (self->vr, i, 0), 1, ncm_vector_ptr (self->spec_C, 1), 2);
  }
}


/**
 * nc_hiqg_1d_spec_evol:
 * @qg1d: a #NcHIQG1D
 * @t: FIXME
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_spec_evol (NcHIQG1D *qg1d, const gdouble t)
{
#if HAVE_SUNDIALS_MAJOR == 3
  NcHIQG1DPrivate * const self = qg1d->priv;
  /*gdouble *Y = N_VGetArrayPointer (self->yBohm);*/
  gdouble ts = 0.0;
  gint flag;

  if (t > 0.0)
  {
    flag = ARKodeSetStopTime (self->bohm, t);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

    flag = ARKode (self->bohm, t, self->yBohm, &ts, ARK_NORMAL);
    NCM_CVODE_CHECK (&flag, "ARKode", 1, );
#endif /* HAVE_SUNDIALS_MAJOR == 3 */
  }

  /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", t, Y[0], Y[1], Y[2], Y[3]);*/

  _nc_hiqg_1d_spec_evol_C (qg1d, t);
}

/**
 * nc_hiqg_1d_spec_eval_psi:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi_0$
 *
 * FIXME
 *
 */
void
nc_hiqg_1d_spec_eval_psi (NcHIQG1D *qg1d, const gdouble x, gdouble *psi)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gint i;

  psi[0] = 0.0;
  psi[1] = 0.0;
  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi    = ncm_vector_get (self->knots, i);
    const gdouble ReC_i = ncm_vector_get (self->spec_C, 2 * i + 0);
    const gdouble ImC_i = ncm_vector_get (self->spec_C, 2 * i + 1);
    const gdouble Kxxi  = _nc_hiqg_1d_spec_basis (x, xi, self->h, self->basis_a);

    psi[0] += ReC_i * Kxxi;
    psi[1] += ImC_i * Kxxi;
  }
}

/**
 * nc_hiqg_1d_spec_eval_dS:
 * @qg1d: a #NcHIQG1D
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_eval_dS (NcHIQG1D *qg1d, const gdouble x)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble dS_rho_x3 = 0.0;
  gdouble Re_psi_x  = 0.0;
  gdouble Im_psi_x  = 0.0;
  gint i, j;

  /*return 1.0;*/

  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi     = ncm_vector_get (self->knots, i);
    const gdouble ReC_i  = ncm_vector_get (self->spec_C, 2 * i + 0);
    const gdouble ImC_i  = ncm_vector_get (self->spec_C, 2 * i + 1);
    const gdouble Kxxi_x = _nc_hiqg_1d_spec_basis_xxa (x, xi, self->h, self->basis_a);

    Re_psi_x += ReC_i * Kxxi_x;
    Im_psi_x += ImC_i * Kxxi_x;

    /*printf ("PSI_x % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", self->h, x, xi, Kxxi_x, Re_psi_x, Im_psi_x);*/

    for (j = i + 1; j < self->nknots; j++)
    {
      const gdouble xj    = ncm_vector_get (self->knots, j);
      const gdouble ReC_j = ncm_vector_get (self->spec_C, 2 * j + 0);
      const gdouble ImC_j = ncm_vector_get (self->spec_C, 2 * j + 1);
      const gdouble Sij   = _nc_hiqg_1d_spec_Sbasis_x3 (x, xi, xj, self->h, self->basis_a);
      const gdouble Cij   = (ImC_j * ReC_i - ReC_j * ImC_i);

      /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", self->h, x, xi, xj, Sij, Cij, dS_rho_x3);*/
      dS_rho_x3 += Sij * Cij;
    }
  }
  /*printf ("# BLA % 22.15g % 22.15g % 22.15g\n", x * x * x * dS_rho_x3, x * x * (Re_psi_x * Re_psi_x + Im_psi_x * Im_psi_x), x * dS_rho_x3 / (Re_psi_x * Re_psi_x + Im_psi_x * Im_psi_x));*/

  return 2.0 * x * dS_rho_x3 / (Re_psi_x * Re_psi_x + Im_psi_x * Im_psi_x);
}

/**
 * nc_hiqg_1d_spec_int_rho_0_inf:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_int_rho_0_inf (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble int_rho = 0.0;
  gint i, j;

  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi       = ncm_vector_get (self->knots, i);
    const gdouble ReC_i    = ncm_vector_get (self->spec_C, 2 * i + 0);
    const gdouble ImC_i    = ncm_vector_get (self->spec_C, 2 * i + 1);
    const gdouble IKxixi_x = _nc_hiqg_1d_spec_If2 (xi, xi, self->h);

    int_rho += (ReC_i * ReC_i + ImC_i * ImC_i) * IKxixi_x;

    for (j = i + 1; j < self->nknots; j++)
    {
      const gdouble xj       = ncm_vector_get (self->knots, j);
      const gdouble ReC_j    = ncm_vector_get (self->spec_C, 2 * j + 0);
      const gdouble ImC_j    = ncm_vector_get (self->spec_C, 2 * j + 1);
      const gdouble IKxixi_x = _nc_hiqg_1d_spec_If2 (xi, xj, self->h);
      const gdouble CRij     = 2.0 * (ReC_i * ReC_j + ImC_i * ImC_j);

      int_rho += IKxixi_x * CRij;
    }
  }

  return int_rho;
}

/**
 * nc_hiqg_1d_spec_int_xrho_0_inf:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_int_xrho_0_inf (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  gdouble int_xrho = 0.0;
  gint i, j;

  for (i = 0; i < self->nknots; i++)
  {
    const gdouble xi       = ncm_vector_get (self->knots, i);
    const gdouble ReC_i    = ncm_vector_get (self->spec_C, 2 * i + 0);
    const gdouble ImC_i    = ncm_vector_get (self->spec_C, 2 * i + 1);
    const gdouble IKxixi_x = _nc_hiqg_1d_spec_Ixf2 (xi, xi, self->h);

    int_xrho += (ReC_i * ReC_i + ImC_i * ImC_i) * IKxixi_x;

    for (j = i + 1; j < self->nknots; j++)
    {
      const gdouble xj       = ncm_vector_get (self->knots, j);
      const gdouble ReC_j    = ncm_vector_get (self->spec_C, 2 * j + 0);
      const gdouble ImC_j    = ncm_vector_get (self->spec_C, 2 * j + 1);
      const gdouble IKxixi_x = _nc_hiqg_1d_spec_Ixf2 (xi, xj, self->h);
      const gdouble CRij     = 2.0 * (ReC_i * ReC_j + ImC_i * ImC_j);

      int_xrho += IKxixi_x * CRij;
    }
  }

  return int_xrho;
}

/**
 * nc_hiqg_1d_spec_nBohm:
 * @qg1d: a #NcHIQG1D
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint
nc_hiqg_1d_spec_nBohm (NcHIQG1D *qg1d)
{
  NcHIQG1DPrivate * const self = qg1d->priv;

  return self->nBohm;
}

/**
 * nc_hiqg_1d_spec_Bohm:
 * @qg1d: a #NcHIQG1D
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hiqg_1d_spec_Bohm (NcHIQG1D *qg1d, gint i)
{
  NcHIQG1DPrivate * const self = qg1d->priv;
  return NV_Ith_S (self->yBohm, i);
}

