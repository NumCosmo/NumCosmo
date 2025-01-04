/***************************************************************************
 *            ncm_fftlog_sbessel_jljm.c
 *
 *  Sun March 24 16:53:52 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fftlog_sbessel_jljm.c
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmFftlogSBesselJLJM:
 *
 * Logarithm fast fourier transform for the base kernel for angular projections.
 *
 * This object computes the function (see #NcmFftlog) $$Y_n = \int_0^\infty
 * t^{\frac{2\pi i n}{L}} K(t) dt,$$ where the kernel are the product of spherical
 * bessel function of the first kind $K(t) = t^q j_{\ell}(t r) j_{\ell+\delta\ell}(t /
 * r)$, where $\delta\ell = m - l$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog_sbessel_jljm.h"
#include "math/ncm_cfg.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>
#include <math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif /* HAVE_FFTW3 */

#ifdef HAVE_ACB_H
#ifdef HAVE_FLINT_ACB_H
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#else /* HAVE_FLINT_ACB_H */
#include <acb.h>
#include <acb_hypgeom.h>
#endif /* HAVE_FLINT_ACB_H */
#endif /* HAVE_ACB_H */
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmFftlogSBesselJLJMInt
{
  gdouble c;
  gdouble _2ar_p_1;
  gdouble _2ar_p_3;
  gdouble ar2_ai2;
  gdouble arp12_ai2;
} NcmFftlogSBesselJLJMInt;

typedef struct _NcmFftlogSBesselJLJMPrivate
{
  gint ell;
  gint dell;
  gdouble q;
  gdouble lnw;
  gdouble w;
  N_Vector y;
  SUNMatrix J;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
  gpointer cvode;
  gboolean cvode_init;
  NcmFftlogSBesselJLJMInt int_data;
} NcmFftlogSBesselJLJMPrivate;

struct _NcmFftlogSBesselJLJM
{
  NcmFftlog parent_instance;
};

enum
{
  PROP_0,
  PROP_ELL,
  PROP_DELL,
  PROP_LNW,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFftlogSBesselJLJM, ncm_fftlog_sbessel_jljm, NCM_TYPE_FFTLOG)

static void
ncm_fftlog_sbessel_jljm_init (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  self->ell        = 0;
  self->dell       = 0;
  self->lnw        = G_MAXDOUBLE;
  self->w          = 0.0;
  self->y          = N_VNew_Serial (4);
  self->J          = SUNDenseMatrix (4, 4);
  self->LS         = SUNLinSol_Dense (self->y, self->J);
  self->NLS        = SUNNonlinSol_Newton (self->y);
  self->cvode      = CVodeCreate (CV_BDF);
  self->cvode_init = FALSE;
}

static void
_ncm_fftlog_sbessel_jljm_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFftlogSBesselJLJM *fftlog_jljm = NCM_FFTLOG_SBESSEL_JLJM (object);

  g_return_if_fail (NCM_IS_FFTLOG_SBESSEL_JLJM (object));

  switch (prop_id)
  {
    case PROP_ELL:
      ncm_fftlog_sbessel_jljm_set_ell (fftlog_jljm, g_value_get_int (value));
      break;
    case PROP_DELL:
      ncm_fftlog_sbessel_jljm_set_dell (fftlog_jljm, g_value_get_int (value));
      break;
    case PROP_LNW:
      ncm_fftlog_sbessel_jljm_set_lnw (fftlog_jljm, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fftlog_sbessel_jljm_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFftlogSBesselJLJM *fftlog_jljm = NCM_FFTLOG_SBESSEL_JLJM (object);

  g_return_if_fail (NCM_IS_FFTLOG_SBESSEL_JLJM (object));

  switch (prop_id)
  {
    case PROP_ELL:
      g_value_set_int (value, ncm_fftlog_sbessel_jljm_get_ell (fftlog_jljm));
      break;
    case PROP_DELL:
      g_value_set_int (value, ncm_fftlog_sbessel_jljm_get_dell (fftlog_jljm));
      break;
    case PROP_LNW:
      g_value_set_double (value, ncm_fftlog_sbessel_jljm_get_lnw (fftlog_jljm));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fftlog_sbessel_jljm_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_sbessel_jljm_parent_class)->dispose (object);
}

static void
_ncm_fftlog_sbessel_jljm_finalize (GObject *object)
{
  NcmFftlogSBesselJLJM *fftlog_jljm        = NCM_FFTLOG_SBESSEL_JLJM (object);
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode = NULL;
  }

  if (self->y != NULL)
  {
    N_VDestroy (self->y);
    self->y = NULL;
  }

  if (self->NLS != NULL)
  {
    SUNNonlinSolFree (self->NLS);
    self->NLS = NULL;
  }

  if (self->J != NULL)
  {
    SUNMatDestroy (self->J);
    self->J = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_sbessel_jljm_parent_class)->finalize (object);
}

static void _ncm_fftlog_sbessel_jljm_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0);

static void
ncm_fftlog_sbessel_jljm_class_init (NcmFftlogSBesselJLJMClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->set_property = &_ncm_fftlog_sbessel_jljm_set_property;
  object_class->get_property = &_ncm_fftlog_sbessel_jljm_get_property;
  object_class->dispose      = &_ncm_fftlog_sbessel_jljm_dispose;
  object_class->finalize     = &_ncm_fftlog_sbessel_jljm_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ELL,
                                   g_param_spec_int ("ell",
                                                     NULL,
                                                     "Spherical Bessel integer order j_\\ell j_{\\ell+d\\ell}",
                                                     0, G_MAXINT32, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DELL,
                                   g_param_spec_int ("dell",
                                                     NULL,
                                                     "Spherical Bessel integer order difference $j_\\ell j_{\\ell+d\\ell}$",
                                                     G_MININT32, G_MAXINT32, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNW,
                                   g_param_spec_double ("lnw",
                                                        NULL,
                                                        "Spherical Bessel scale difference log(w)",
                                                        GSL_LOG_DBL_MIN, 0.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  fftlog_class->name       = "sbessel_jljm";
  fftlog_class->compute_Ym = &_ncm_fftlog_sbessel_jljm_compute_Ym;
}

#if defined (HAVE_FFTW3) && defined (HAVE_ACB_H)

static gint _ncm_fftlog_sbessel_jljm_cpu_integrate_2f1_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_fftlog_sbessel_jljm_cpu_integrate_2f1_jac (gdouble x, N_Vector y, N_Vector fy, SUNMatrix Jac, gpointer user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint
_ncm_fftlog_sbessel_jljm_cpu_integrate_2f1_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmFftlogSBesselJLJMInt *data = (NcmFftlogSBesselJLJMInt *) f_data;
  const gdouble w1              = NV_Ith_S (y, 0);
  const gdouble w1p             = NV_Ith_S (y, 1);
  const gdouble w2              = NV_Ith_S (y, 2);
  const gdouble w2p             = NV_Ith_S (y, 3);
  const gdouble denom           = x * (1.0 + x);

  NV_Ith_S (ydot, 0) = w1p;
  NV_Ith_S (ydot, 1) = -(data->ar2_ai2 * w1   + (data->_2ar_p_1 * x + data->c) * w1p) / denom;

  NV_Ith_S (ydot, 2) = w2p;
  NV_Ith_S (ydot, 3) = -(data->arp12_ai2 * w2 + (data->_2ar_p_3 * x + data->c) * w2p) / denom;

  /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x, data->ar2_ai2 * w1, data->_2ar_p_1 * x * w1p, data->c * w1p, NV_Ith_S (ydot, 1));*/

  return 0;
}

static gint
_ncm_fftlog_sbessel_jljm_cpu_integrate_2f1_jac (gdouble x, N_Vector y, N_Vector fy, SUNMatrix Jac, gpointer user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmFftlogSBesselJLJMInt *data = (NcmFftlogSBesselJLJMInt *) user_data;
  const gdouble denom           = x * (1.0 + x);

  SM_ELEMENT_D (Jac, 0, 0) = 1.0;
  SM_ELEMENT_D (Jac, 0, 1) = 0.0;
  SM_ELEMENT_D (Jac, 0, 2) = 0.0;
  SM_ELEMENT_D (Jac, 0, 3) = 0.0;

  SM_ELEMENT_D (Jac, 1, 0) = -data->ar2_ai2 / denom;
  SM_ELEMENT_D (Jac, 1, 1) = -(data->_2ar_p_1 * x + data->c) / denom;
  SM_ELEMENT_D (Jac, 1, 2) = 0.0;
  SM_ELEMENT_D (Jac, 1, 3) = 0.0;

  SM_ELEMENT_D (Jac, 2, 0) = 0.0;
  SM_ELEMENT_D (Jac, 2, 1) = 0.0;
  SM_ELEMENT_D (Jac, 2, 2) = 1.0;
  SM_ELEMENT_D (Jac, 2, 3) = 0.0;

  SM_ELEMENT_D (Jac, 3, 0) = 0.0;
  SM_ELEMENT_D (Jac, 3, 1) = 0.0;
  SM_ELEMENT_D (Jac, 3, 2) = -data->arp12_ai2 / denom;
  SM_ELEMENT_D (Jac, 3, 3) = -(data->_2ar_p_3 * x + data->c) / denom;

  return 0;
}

static gdouble
hyperg_2F1_conj_luke (const double aR, const double aI, const double c, const double xin)
{
  const double RECUR_BIG = 1.0e+50;
  const int nmax         = 10000;
  int n                  = 3;
  const double x         = -xin;
  const double x3        = x * x * x;
  const double atimesb   = aR * aR + aI * aI;
  const double apb       = 2.0 * aR;
  const double t0        = atimesb / c;
  const double t1        = (atimesb +     apb + 1.0) / (2.0 * c);
  const double t2        = (atimesb + 2.0 * apb + 4.0) / (2.0 * (c + 1.0));
  double F               = 1.0;
  double prec;

  double Bnm3 = 1.0;                                 /* B0 */
  double Bnm2 = 1.0 + t1 * x;                        /* B1 */
  double Bnm1 = 1.0 + t2 * x * (1.0 + t1 / 3.0 * x); /* B2 */

  double Anm3 = 1.0;                                                                /* A0 */
  double Anm2 = Bnm2 - t0 * x;                                                      /* A1 */
  double Anm1 = Bnm1 - t0 * (1.0 + t2 * x) * x + t0 * t1 * (c / (c + 1.0)) * x * x; /* A2 */

  while (1)
  {
    double nm1         = n - 1;
    double nm2         = n - 2;
    double npam1_npbm1 = atimesb + nm1 * apb + nm1 * nm1;
    double npam2_npbm2 = atimesb + nm2 * apb + nm2 * nm2;
    double npcm1       = nm1 + c;
    double npcm2       = nm2 + c;
    double tnm1        = 2 * n - 1;
    double tnm3        = 2 * n - 3;
    double tnm5        = 2 * n - 5;
    double n2          = n * n;
    double F1          =  (3.0 * n2 + (apb - 6) * n + 2 - atimesb - 2 * apb) / (2 * tnm3 * npcm1);
    double F2          = -(3.0 * n2 - (apb + 6) * n + 2 - atimesb) * npam1_npbm1 / (4 * tnm1 * tnm3 * npcm2 * npcm1);
    double F3          = (npam2_npbm2 * npam1_npbm1 * (nm2 * nm2 - nm2 * apb + atimesb)) / (8 * tnm3 * tnm3 * tnm5 * (n + c - 3) * npcm2 * npcm1);
    double E           = -npam1_npbm1 * (n - c - 1) / (2 * tnm3 * npcm2 * npcm1);

    double An = (1.0 + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
    double Bn = (1.0 + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
    double r  = An / Bn;

    prec = fabs (F - r) / fabs (F);
    F    = r;

    if ((prec < GSL_DBL_EPSILON) || (n > nmax))
      break;

    if ((fabs (An) > RECUR_BIG) || (fabs (Bn) > RECUR_BIG))
    {
      An   /= RECUR_BIG;
      Bn   /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
      Anm3 /= RECUR_BIG;
      Bnm3 /= RECUR_BIG;
    }
    else if ((fabs (An) < 1.0 / RECUR_BIG) || (fabs (Bn) < 1.0 / RECUR_BIG))
    {
      An   *= RECUR_BIG;
      Bn   *= RECUR_BIG;
      Anm1 *= RECUR_BIG;
      Bnm1 *= RECUR_BIG;
      Anm2 *= RECUR_BIG;
      Bnm2 *= RECUR_BIG;
      Anm3 *= RECUR_BIG;
      Bnm3 *= RECUR_BIG;
    }

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }

  return F;
}

static complex double
_ncm_fftlog_sbessel_jljm_cpu_integrate_2f1 (NcmFftlogSBesselJLJMPrivate * const self, const gdouble ar, const gdouble ai, const gdouble c, const gdouble l, const gdouble xref, const gdouble x)
{
  const gdouble x_ini        = (xref > 1.0) ? 0.99 : xref;
  const gdouble _2f1_ini     = gsl_sf_hyperg_2F1_conj (ar, ai, c, -x_ini);
  const gdouble _d2f1_ini    = -(ar * ar + ai * ai) / c * gsl_sf_hyperg_2F1_conj (ar + 1.0, ai, c + 1.0, -x_ini);
  const gdouble _2f1_p1_ini  = gsl_sf_hyperg_2F1_conj (ar + 1.0, ai, c, -x_ini);
  const gdouble _d2f1_p1_ini = -(gsl_pow_2 (ar + 1.0) + ai * ai) / c * gsl_sf_hyperg_2F1_conj (ar + 2.0, ai, c + 1.0, -x_ini);
  gdouble xf                 = 0.0;
  gint flag;

  printf ("# INI % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", x_ini, _2f1_ini, _d2f1_ini, _2f1_p1_ini, _d2f1_p1_ini, ar, ai, c);

  NV_Ith_S (self->y, 0) = _2f1_ini;
  NV_Ith_S (self->y, 1) = _d2f1_ini;
  NV_Ith_S (self->y, 2) = _2f1_p1_ini;
  NV_Ith_S (self->y, 3) = _d2f1_p1_ini;

  self->int_data.c         = 1.5 + self->ell;
  self->int_data._2ar_p_1  = 2.0 * ar + 1.0;
  self->int_data._2ar_p_3  = 2.0 * ar + 3.0;
  self->int_data.ar2_ai2   = ar * ar + ai * ai;
  self->int_data.arp12_ai2 = gsl_pow_2 (ar + 1.0) + ai * ai;

  if (!self->cvode_init)
  {
    flag = CVodeInit (self->cvode, &_ncm_fftlog_sbessel_jljm_cpu_integrate_2f1_f, x_ini, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, 0.0);
    self->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode, x_ini, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, 0.0);
  }

  flag = CVodeSetNonlinearSolver (self->cvode, self->NLS);
  NCM_CVODE_CHECK (&flag, "CVodeSetNonlinearSolver", 1, 0.0);

  flag = CVodeSetLinearSolver (self->cvode, self->LS, self->J);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, 0.0);

  flag = CVodeSStolerances (self->cvode, 1.0e-13, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, 0.0);

  flag = CVodeSetMaxNumSteps (self->cvode, 10000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, 0.0);

  flag = CVodeSetJacFn (self->cvode, &_ncm_fftlog_sbessel_jljm_cpu_integrate_2f1_jac);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, 0.0);

  flag = CVodeSetUserData (self->cvode, &self->int_data);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, 0.0);

  flag = CVodeSetStopTime (self->cvode, x);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, 0.0);

  do {
    flag = CVode (self->cvode, x, self->y, &xf, CV_NORMAL);
    NCM_CVODE_CHECK (&flag, "CVode", 1, 0.0);

    printf ("# STEP % 22.15g % 22.15g % 22.15g % 11.5e % 22.15g % 22.15g % 22.15g % 11.5e % 22.15g\n", xf,
            NV_Ith_S (self->y, 0), hyperg_2F1_conj_luke (ar, ai, c, -xf), fabs (hyperg_2F1_conj_luke (ar, ai, c, -xf) / NV_Ith_S (self->y, 0) - 1.0),
            NV_Ith_S (self->y, 1),
            NV_Ith_S (self->y, 2), hyperg_2F1_conj_luke (ar + 1.0, ai, c, -xf), fabs (hyperg_2F1_conj_luke (ar + 1.0, ai, c, -xf) / NV_Ith_S (self->y, 2) - 1.0),
            NV_Ith_S (self->y, 3));
  } while (x != xf);

  {
    const gdouble _2lp1     = 2.0 * l + 1.0;
    const gdouble xp1       = 1.0 + x;
    const complex double a  = ar + I * ai;
    const complex double ac = ar - I * ai;
    const gdouble _2f1      = NV_Ith_S (self->y, 0);
    const gdouble _2f1_p1   = NV_Ith_S (self->y, 2);

    return cpow (xp1, a) * ((2.0 * ac - _2lp1) * _2f1 + 2.0 * a * xp1 * _2f1_p1) / (4.0 * ar - _2lp1);
  }
}

static complex double
_ncm_fftlog_sbessel_jljm_cpu_direct_2f1 (const gdouble ar, const gdouble ai, const gdouble c, const gdouble l, const gdouble x)
{
  const gdouble _2lp1     = 2.0 * l + 1.0;
  const gdouble xp1       = 1.0 + x;
  const complex double a  = ar + I * ai;
  const complex double ac = ar - I * ai;
  const gdouble _2f1      = gsl_sf_hyperg_2F1_conj (ar,       ai, c, -x);
  const gdouble _2f1_p1   = gsl_sf_hyperg_2F1_conj (ar + 1.0, ai, c, -x);

  return cpow (xp1, a) * ((2.0 * ac - _2lp1) * _2f1 + 2.0 * a * xp1 * _2f1_p1) / (4.0 * ar - _2lp1);
}

static complex double
_ncm_fftlog_sbessel_jljm_cpu_2f1 (NcmFftlogSBesselJLJMPrivate * const self, const gdouble L, const gint n, const gdouble x)
{
  const gdouble p    = 1.0;
  const gdouble ar   = -0.5 * self->dell;
  const gdouble ai   = M_PI * n / L;
  const gdouble c    = 1.5 + self->ell;
  const gdouble xref = (n == 0) ? 1000.0 : c * p  / (gsl_pow_2 (ai) + gsl_pow_2 (ar));

  if (x < xref)
    return _ncm_fftlog_sbessel_jljm_cpu_direct_2f1 (ar, ai, c, self->ell, x);
  else
    return _ncm_fftlog_sbessel_jljm_cpu_integrate_2f1 (self, ar, ai, c, self->ell, xref, x);
}

static void
_ncm_fftlog_sbessel_jljm_compute_2f1 (NcmFftlog *fftlog)
{
  NcmFftlogSBesselJLJM *fftlog_jljm        = NCM_FFTLOG_SBESSEL_JLJM (fftlog);
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  complex double res = _ncm_fftlog_sbessel_jljm_cpu_2f1 (self, ncm_fftlog_get_full_length (fftlog), 1000, 2.0);

  printf ("# RES % 22.15g % 22.15g\n", creal (res), cimag (res));
}

#endif /* defined (HAVE_FFTW3) && defined (HAVE_ACB_H) */

#define _NCM_FFTLOG_SBESSEL_JLJM_CACHE_FILE "ncm_fftlog_sbessel_j%dj%d_lnw%.14g_L%.14g.fftlog", self->ell, self->ell + self->dell, self->lnw, ncm_fftlog_get_full_length (fftlog)

static void
_ncm_fftlog_sbessel_jljm_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
#if defined (HAVE_FFTW3) && defined (HAVE_ACB_H)
  NcmFftlogSBesselJLJM *fftlog_jljm        = NCM_FFTLOG_SBESSEL_JLJM (fftlog);
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);
  const guint maxprec                      = 64 * 64;
  guint prec                               = 64 * 2;
  const gint Nf                            = ncm_fftlog_get_full_size (fftlog);
  const gint Nf_2                          = Nf / 2;
  fftw_complex *Ym_base                    = (fftw_complex *) Ym_0;
  fftw_complex *Ym_cache                   = NULL;
  NcmVector *Ym_v_cache                    = NULL;
  gint Nf_cached                           = 0;
  gint Nf_cached_2                         = 0;
  acb_t two, onehalf, twopi_Lt, Lt, a_n, A_n, B_n, C_n, D_n, E_n, F_n, N_n, w4, x, q, lnw, res, pi;
  gint i, i_ini = 0, i_size = Nf;
  glong prec_sf;

  if (FALSE)
    _ncm_fftlog_sbessel_jljm_compute_2f1 (fftlog);

  /* printf ("# Computing coefs for ell = %4d, delta ell = %4d, Nf = %6d, w = % 22.15g, w4 = % 22.15g\n", self->ell, self->dell, Nf, exp (self->lnw), exp (4.0 * self->lnw)); */

  if (ncm_cfg_exists (_NCM_FFTLOG_SBESSEL_JLJM_CACHE_FILE))
  {
    gchar *cache_file = ncm_cfg_get_fullpath (_NCM_FFTLOG_SBESSEL_JLJM_CACHE_FILE);
    NcmVector *Ym_v   = Ym_v_cache = NCM_VECTOR (ncm_serialize_global_from_binfile (cache_file));
    gint Ym_v_len     = ncm_vector_len (Ym_v);

    Nf_cached   = Ym_v_len / 2;
    Nf_cached_2 = Nf_cached / 2;

    g_assert_cmpint (Ym_v_len % 2, ==, 0);

    /*printf ("# Found cached coefs, Nf %u!\n", Nf_cached);*/
    Ym_cache = (fftw_complex *) ncm_vector_data (Ym_v);

    if (Nf_cached >= Nf)
    {
      memcpy (&Ym_base[0], &Ym_cache[0], sizeof (fftw_complex) * (Nf_2 + 1));
      memcpy (&Ym_base[Nf_2 + 1], &Ym_cache[Nf_cached - (Nf_2 - 1)], sizeof (fftw_complex) * (Nf_2 - 1));

      ncm_vector_clear (&Ym_v_cache);
      g_free (cache_file);

      return;
    }
    else
    {
      memcpy (&Ym_base[0], &Ym_cache[0], sizeof (fftw_complex) * (Nf_cached_2 + 1));
      memcpy (&Ym_base[Nf - (Nf_cached_2 - 1)], &Ym_cache[Nf_cached - (Nf_cached_2 - 1)], sizeof (fftw_complex) * (Nf_cached_2 - 1));

      /*printf ("# INITIAL I %d I SIZE %d\n", Nf_cached_2 + 1, Nf - Nf_cached);*/
      i_ini  = Nf_cached_2 + 1;
      i_size = i_ini + (Nf - Nf_cached);
    }

    g_free (cache_file);
  }

  acb_init (pi);
  acb_init (Lt);
  acb_init (two);
  acb_init (onehalf);
  acb_init (twopi_Lt);
  acb_init (a_n);
  acb_init (A_n);
  acb_init (B_n);
  acb_init (C_n);
  acb_init (D_n);
  acb_init (E_n);
  acb_init (F_n);
  acb_init (N_n);
  acb_init (w4);
  acb_init (x);
  acb_init (q);
  acb_init (lnw);
  acb_init (res);

  acb_set_d (Lt,  ncm_fftlog_get_full_length (fftlog));
  acb_set_d (q,   self->q);
  acb_set_d (lnw, self->lnw);

  acb_const_pi (pi, maxprec);

  acb_set (twopi_Lt, pi);                      /* twopi_Lt = pi */
  acb_mul_ui (twopi_Lt, twopi_Lt, 2, maxprec); /* twopi_Lt = 2 * twopi_Lt */
  acb_div (twopi_Lt, twopi_Lt, Lt, maxprec);   /* twopi_Lt = twopi_Lt / Lt */

  acb_set_ui (two, 2);                       /* two = 2 */
  acb_set_ui (onehalf, 1);                   /* onehalf = 1 */
  acb_div_ui (onehalf, onehalf, 2, maxprec); /* onehalf = onehalf / 2 */

  acb_mul_ui (w4, lnw, 4, maxprec); /* w4 = lnw * 4 */
  acb_exp (w4, w4, maxprec);        /* w4 = exp (w4) */

  acb_mul_si (x, w4, -2, maxprec); /* x = -2 * w4 */
  acb_add_ui (x, x, 1, maxprec);   /* x = x + 1 */

  for (i = i_ini; i < i_size; i++)
  {
    const gint phys_i = ncm_fftlog_get_mode_index (fftlog, i);

    acb_mul_si (a_n, twopi_Lt, phys_i, prec); /* a_n = twopi_Lt * phys_i */
    acb_mul_onei (a_n, a_n);                  /* a_n = a_n * I */

    if (FALSE)
    {
      acb_set (A_n, q);                        /* A_n = q */
      acb_add (A_n, A_n, a_n, prec);           /* A_n = A_n + a_n */
      acb_add_ui (A_n, A_n, 1, prec);          /* A_n = A_n + 1 */
      acb_add_si (A_n, A_n, self->ell,  prec); /* A_n = A_n + ell */
      acb_add_si (A_n, A_n, self->ell,  prec); /* A_n = A_n + ell */
      acb_add_si (A_n, A_n, self->dell, prec); /* A_n = A_n + dell */
      acb_div_ui (A_n, A_n, 2, prec);          /* A_n = A_n / 2 */

      acb_set_si (B_n, -1);                    /* B_n = -1 */
      acb_div_ui (B_n, B_n, 2, prec);          /* B_n = B_n / 2 */
      acb_add (B_n, B_n, A_n, prec);           /* B_n = B_n + A_n */
      acb_sub_si (B_n, B_n, self->ell,  prec); /* B_n = B_n - ell */
      acb_sub_si (B_n, B_n, self->dell, prec); /* B_n = B_n - dell */

      acb_set_ui (C_n, 3);                     /* C_n = 3 */
      acb_div_ui (C_n, C_n, 2, prec);          /* C_n = C_n / 2 */
      acb_sub (C_n, C_n, A_n, prec);           /* C_n = C_n - A_n */
      acb_add_si (C_n, C_n, self->ell,  prec); /* C_n = C_n + ell */
      acb_add_si (C_n, C_n, self->dell, prec); /* C_n = C_n + dell */
      acb_gamma (C_n, C_n, prec);              /* C_n = Gamma (C_n) */

      acb_set_ui (D_n, 3);                    /* D_n = 3 */
      acb_div_ui (D_n, D_n, 2, prec);         /* D_n = D_n / 2 */
      acb_add_ui (D_n, D_n, self->ell, prec); /* D_n = D_n + ell */

      acb_set (E_n, A_n);                      /* E_n = A_n */
      acb_mul_ui (E_n, E_n, 2, prec);          /* E_n = E_n * 2 */
      acb_sub_si (E_n, E_n, self->dell, prec); /* E_n = E_n - dell */
      acb_mul (E_n, E_n, lnw, prec);           /* E_n = E_n * lnw */
      acb_exp (E_n, E_n, prec);                /* E_n = exp (E_n) */

      acb_set (F_n, A_n);                      /* F_n = A_n */
      acb_mul_ui (F_n, F_n, 2, prec);          /* F_n = F_n * 2 */
      acb_sub_si (F_n, F_n, self->ell, prec);  /* F_n = F_n - ell */
      acb_sub_si (F_n, F_n, self->ell, prec);  /* F_n = F_n - ell */
      acb_sub_si (F_n, F_n, self->dell, prec); /* F_n = F_n - dell */
      acb_sub_ui (F_n, F_n, 3, prec);          /* F_n = F_n - 3 */
      acb_pow (F_n, two, F_n, prec);           /* F_n = 2^(E_n) */

      acb_hypgeom_2f1 (res, A_n, B_n, D_n, w4, ACB_HYPGEOM_2F1_REGULARIZED, prec); /* res = 2f1 (A_n, B_n, D_n; w4) */
      /*acb_hypgeom_2f1_transform (res, A_n, B_n, D_n, w4, ACB_HYPGEOM_2F1_REGULARIZED, 4, prec); */
      /*acb_hypgeom_2f1_direct (res, A_n, B_n, D_n, w4, 1, prec); */

      acb_mul (res, res, pi, prec);  /* res = res * pi */
      acb_mul (res, res, E_n, prec); /* res = res * E_n */
      acb_mul (res, res, F_n, prec); /* res = res * F_n */

      acb_gamma (A_n, A_n, prec);    /* A_n = Gamma (A_n) */
      acb_mul (res, res, A_n, prec); /* res = res * A_n */
      acb_div (res, res, C_n, prec); /* res = res / C_n */
    }
    else
    {
      acb_set (A_n, onehalf);                 /* A_n = 1 / 2 */
      acb_add_si (A_n, A_n, self->ell, prec); /* A_n = A_n + ell */

      acb_add (B_n, a_n, q, prec);    /* B_n = a_n + q */
      acb_sub_ui (B_n, B_n, 1, prec); /* B_n = B_n - 1 */

      acb_add (C_n, a_n, q, prec);                                 /* C_n = a_n + q */
      acb_add_si (C_n, C_n, 2 * self->ell + self->dell + 1, prec); /* C_n = C_n + 2ell + dell + 1 */
      acb_div_ui (C_n, C_n, 2, prec);                              /* C_n = C_n / 2 */

      acb_conj (D_n, C_n);            /* D_n = C_n^* */
      acb_sub (D_n, D_n, q, prec);    /* D_n = D_n - q */
      acb_add_ui (D_n, D_n, 1, prec); /* D_n = D_n + 1 */

      acb_neg (N_n, a_n);                      /* N_n = -a_n */
      acb_add_si (N_n, N_n, self->dell, prec); /* N_n = N_n + dell */
      acb_sub (N_n, N_n, q, prec);             /* N_n = N_n - q */
      acb_div_ui (N_n, N_n, 2, prec);          /* N_n = N_n / 2 */

      acb_add (E_n, a_n, q, prec);                    /* E_n = a_n + q */
      acb_set (F_n, E_n);                             /* F_n = E_n */
      acb_sub_ui (E_n, E_n, 2, prec);                 /* E_n = E_n - 2 */
      acb_add_si (F_n, F_n, 2 * self->ell + 1, prec); /* F_n = F_n + 2 ell + 1 */
      acb_mul (F_n, F_n, lnw, prec);                  /* F_n = F_n * ln (w) */

      acb_lgamma (C_n, C_n, prec); /* C_n = ln (Gamma (C_n)) */
      acb_lgamma (D_n, D_n, prec); /* D_n = ln (Gamma (D_n)) */

      acb_add (F_n, F_n, C_n, prec); /* F_n = F_n + C_n */
      acb_sub (F_n, F_n, D_n, prec); /* F_n = F_n - D_n */

      acb_pow (E_n, two, E_n, prec); /* E_n = 2^E_n */
      acb_exp (F_n, F_n, prec);      /* F_n = e^F_n */

      acb_hypgeom_jacobi_p (res, N_n, A_n, B_n, x, prec); /* res = P_{N_n}^{A_n,B_n}(x) */
      prec_sf = acb_rel_accuracy_bits (res);

      /*F2f1 = ncm_acb_get_complex (res);*/

      acb_mul (res, res, pi, prec);  /* res = res * pi */
      acb_mul (res, res, E_n, prec); /* res = res * E_n */
      acb_mul (res, res, F_n, prec); /* res = res * F_n */
    }

    prec_sf = acb_rel_accuracy_bits (res);
    /*prec_sf = prec_sf < 0 ? 0 : prec_sf;*/

    if (prec_sf < 53)
    {
      /*printf ("INCREASING PREC: %6d %4ld %4u => %4u\n", i, prec_sf, prec, prec + 16);*/
      if (prec + 16 > maxprec)
        g_error ("_ncm_fftlog_sbessel_jljm_compute_Ym: max precision was not enough, giving up.");

      i    -= 1;
      prec += 16;
      continue;
    }
    else if (prec_sf > 96)
    {
      /*printf ("DECREASING PREC: %6d %4ld %4u => %4u\n", i, prec_sf, prec, prec - 2);*/
      prec -= 2;
    }

    /*printf ("\r% 22d", phys_i);*/
    if (FALSE)
    {
      if ((phys_i <= Nf_cached_2) && (phys_i > -Nf_cached_2))
      {
        gint cache_array_index = (phys_i < 0) ? phys_i + Nf_cached : phys_i;

        printf ("# CMP CACHE %d %d %d %d %e % 22.15g % 22.15g % 22.15g % 22.15g\n",
                i, phys_i, cache_array_index, Nf_cached, cabs (Ym_base[i] / Ym_cache[cache_array_index] - 1.0), creal (Ym_base[i]), cimag (Ym_base[i]), creal (Ym_cache[cache_array_index]), cimag (Ym_cache[cache_array_index]));
      }
    }

    Ym_base[i] = ncm_acb_get_complex (res);

/*
 *   printf ("%d %d %d %u %ld % 22.15g % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g | % 22.15g % 22.15g | % 22.15g % 22.15g | % 22.15g % 22.15g | % 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g % 22.15g\n",
 *           phys_i, self->ell, self->dell,
 *           prec, prec_sf,
 *           self->q, self->w, gsl_pow_4 (self->w),
 *           ncm_fftlog_get_full_length (fftlog),
 *           creal (Ym_base[i]), cimag (Ym_base[i]),
 *           creal (ncm_acb_get_complex (N_n)),
 *           cimag (ncm_acb_get_complex (N_n)),
 *           creal (ncm_acb_get_complex (A_n)),
 *           cimag (ncm_acb_get_complex (A_n)),
 *           creal (ncm_acb_get_complex (B_n)),
 *           cimag (ncm_acb_get_complex (B_n)),
 *           creal (ncm_acb_get_complex (x)),
 *           cimag (ncm_acb_get_complex (x)),
 *           creal (ncm_acb_get_complex (C_n)),
 *           cimag (ncm_acb_get_complex (C_n)),
 *           creal (ncm_acb_get_complex (D_n)),
 *           cimag (ncm_acb_get_complex (D_n))
 *           );
 */
  }

  acb_clear (pi);
  acb_clear (Lt);
  acb_clear (two);
  acb_clear (onehalf);
  acb_clear (twopi_Lt);
  acb_clear (a_n);
  acb_clear (A_n);
  acb_clear (B_n);
  acb_clear (C_n);
  acb_clear (D_n);
  acb_clear (E_n);
  acb_clear (F_n);
  acb_clear (N_n);
  acb_clear (w4);
  acb_clear (x);
  acb_clear (q);
  acb_clear (lnw);
  acb_clear (res);

  ncm_vector_clear (&Ym_v_cache);

  {
    NcmVector *Ym_v   = ncm_vector_new_data_static (Ym_0, 2 * Nf, 1);
    gchar *cache_file = ncm_cfg_get_fullpath (_NCM_FFTLOG_SBESSEL_JLJM_CACHE_FILE);

    ncm_serialize_global_to_binfile (G_OBJECT (Ym_v), cache_file);

    g_free (cache_file);
    ncm_vector_free (Ym_v);
  }

#else
  g_error ("_ncm_fftlog_sbessel_jljm_compute_Ym: this object requires both FFTW3 and ARB dependencies.");
#endif /* HAVE_FFTW3 */
}

/**
 * ncm_fftlog_sbessel_jljm_new:
 * @ell: Spherical Bessel Integer order
 * @dell: Spherical Bessel Integer order
 * @lnw: log-scale difference $\ln(r)$
 * @lnr0: output center $\ln(r_0)$
 * @lnk0: input center $\ln(k_0)$
 * @Lk: input/output interval size
 * @N: number of knots
 *
 * Creates a new fftlog Spherical Bessel $j_\ell(xr) j_{\ell+\delta\ell}(x/r)$ object.
 *
 * Returns: (transfer full): a new #NcmFftlogSBesselJLJM
 */
NcmFftlogSBesselJLJM *
ncm_fftlog_sbessel_jljm_new (gint ell, gint dell, gdouble lnw, gdouble lnr0, gdouble lnk0, gdouble Lk, guint N)
{
  NcmFftlogSBesselJLJM *fftlog_jljm = g_object_new (NCM_TYPE_FFTLOG_SBESSEL_JLJM,
                                                    "ell",  ell,
                                                    "dell", dell,
                                                    "lnw",  lnw,
                                                    "lnr0", lnr0,
                                                    "lnk0", lnk0,
                                                    "Lk",   Lk,
                                                    "N",    N,
                                                    NULL);

  return fftlog_jljm;
}

/**
 * ncm_fftlog_sbessel_jljm_set_ell:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 * @ell: Spherical Bessel integer order $\ell$
 *
 * Sets @ell as the Spherical Bessel integer order $\ell$.
 *
 */
void
ncm_fftlog_sbessel_jljm_set_ell (NcmFftlogSBesselJLJM *fftlog_jljm, const gint ell)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  if (self->ell != ell)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jljm);

    g_assert_cmpint (ell, >=, 0);
    g_assert_cmpint (ell + self->dell, >=, 0);

    self->ell = ell;
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_sbessel_jljm_get_ell:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 *
 * Returns: the current Spherical Bessel integer order $\ell$.
 */
gint
ncm_fftlog_sbessel_jljm_get_ell (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  return self->ell;
}

/**
 * ncm_fftlog_sbessel_jljm_set_dell:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 * @dell: Spherical Bessel integer order $\dell$
 *
 * Sets @dell as the Spherical Bessel integer order $\delta\ell$.
 *
 */
void
ncm_fftlog_sbessel_jljm_set_dell (NcmFftlogSBesselJLJM *fftlog_jljm, const gint dell)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  if (self->dell != dell)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jljm);

    if (self->ell >= 0)
      g_assert_cmpint (self->ell + dell, >=, 0);

    self->dell = dell;

    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_sbessel_jljm_get_dell:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 *
 * Returns: the current Spherical Bessel integer order $\delta\ell$.
 */
gint
ncm_fftlog_sbessel_jljm_get_dell (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  return self->dell;
}

/**
 * ncm_fftlog_sbessel_jljm_set_q:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 * @q: Spherical Bessel power factor $q$
 *
 * Sets @q as the Spherical Bessel power $q$.
 *
 */
void
ncm_fftlog_sbessel_jljm_set_q (NcmFftlogSBesselJLJM *fftlog_jljm, const gdouble q)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  if (self->q != q)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jljm);

    self->q = q;
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_sbessel_jljm_get_q:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 *
 * Returns: the current Spherical Bessel power $q$.
 */
gdouble
ncm_fftlog_sbessel_jljm_get_q (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  return self->q;
}

/**
 * ncm_fftlog_sbessel_jljm_set_lnw:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 * @lnw: Spherical Bessel log-scale difference $\ln(w)$
 *
 * Sets @lnw as the Spherical Bessel log-scale difference $\ln(w)$.
 *
 */
void
ncm_fftlog_sbessel_jljm_set_lnw (NcmFftlogSBesselJLJM *fftlog_jljm, const gdouble lnw)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  if (self->lnw != lnw)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jljm);

    g_assert_cmpfloat (GSL_LOG_DBL_MIN, <, lnw);
    g_assert_cmpfloat (lnw, <=, 0.0);

    self->lnw = lnw;
    self->w   = exp (lnw);

    /*printf ("%p % 22.15g % 22.15g\n", self, self->lnw, self->w);*/
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_sbessel_jljm_get_lnw:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 *
 * Returns: the current Spherical Bessel power $\ln(w)$.
 */
gdouble
ncm_fftlog_sbessel_jljm_get_lnw (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);

  return self->lnw;
}

/**
 * ncm_fftlog_sbessel_jljm_set_best_lnr0:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 *
 * Sets the value of $\ln(r_0)$ which gives the best results for
 * the transformation based on the current value of $\ln(k_0)$,
 * this is based in the rule of thumb $\mathrm{max}_{x^*}(j_l)$
 * where $ x^* \approx l + 1$.
 *
 */
void
ncm_fftlog_sbessel_jljm_set_best_lnr0 (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);
  NcmFftlog *fftlog                        = NCM_FFTLOG (fftlog_jljm);

  const gdouble lnk0      = ncm_fftlog_get_lnk0 (fftlog);
  const gdouble r1        = gsl_sf_bessel_zero_Jnu (0.5 + self->ell, 1);
  const gdouble r2        = gsl_sf_bessel_zero_Jnu (0.5 + self->ell + self->dell, 1);
  const gdouble r1m       = 0.5 * (r1 * self->w + r1 / self->w);
  const gdouble r2m       = 0.5 * (r2 * self->w + r2 / self->w);
  const gdouble lnc0      = log (0.5 * (r1m + r2m));
  const gdouble best_lnr0 = -lnk0 + lnc0;

  ncm_fftlog_set_lnr0 (fftlog, best_lnr0);
}

/**
 * ncm_fftlog_sbessel_jljm_set_best_lnk0:
 * @fftlog_jljm: a #NcmFftlogSBesselJLJM
 *
 * Sets the value of $\ln(k_0)$ which gives the best results for
 * the transformation based on the current value of $\ln(r_0)$,
 * this is based in the rule of thumb $\mathrm{max}_{x^*}(j_l)$
 * where $ x^* \approx l + 1$.
 *
 */
void
ncm_fftlog_sbessel_jljm_set_best_lnk0 (NcmFftlogSBesselJLJM *fftlog_jljm)
{
  NcmFftlogSBesselJLJMPrivate * const self = ncm_fftlog_sbessel_jljm_get_instance_private (fftlog_jljm);
  NcmFftlog *fftlog                        = NCM_FFTLOG (fftlog_jljm);

  const gdouble lnr0 = ncm_fftlog_get_lnr0 (fftlog);
  const gdouble r1   = gsl_sf_bessel_zero_Jnu (0.5 + self->ell, 1);
  const gdouble r2   = gsl_sf_bessel_zero_Jnu (0.5 + self->ell + self->dell, 1);
  const gdouble r1m  = 0.5 * (r1 * self->w + r1 / self->w);
  const gdouble r2m  = 0.5 * (r2 * self->w + r2 / self->w);
  const gdouble lnc0 = log (0.5 * (r1m + r2m));

  const gdouble best_lnk0 = -lnr0 + lnc0;

  ncm_fftlog_set_lnk0 (fftlog, best_lnk0);
}

