/***************************************************************************
 *            nc_hipert_adiab.c
 *
 *  Tue June 03 17:20:42 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_adiab
 * @title: Adiabatic Perturbation Object
 * @short_description: Perturbation object for adiabatic mode only 
 *
 * This object provides the computation of the adiabatic mode for the cosmological 
 * perturbations. It solves the equation of motion for the gauge invariant variable
 * (see <link linkend="XVitenti2013">Vitenti (2013)</link> for notation and details)
 * $$\zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + H\mathcal{V}.$$
 * Its conjugated momentum is give by
 * \begin{split}
 * P_\zeta &= \frac{2\bar{D}^2_\bar{K}\Psi}{x^3H},
 * \end{split}
 * 
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{P_\zeta}{m_\zeta}, \\\\
 * P_\zeta^\prime &= -m_\zeta\mu_\zeta^2\zeta.
 * \end{align}
 * The mass $m_\zeta$ and the frequency $\mu_\zeta$ are defined by
 * \begin{align}
 * m_\zeta     &= \frac{3\Delta_\bar{K}(\bar{\rho} + \bar{p})}{\rho_\text{crit0} N x^3 c_s^2 E^2}, \\\\
 * \mu_\zeta^2 &= x^2N^2c_s^2k^2,
 * \end{align}
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure defined by nc_hicosmo_rhopp(),
 * $E^2 = H^2/H_0^2$ is the adimensional Hubble function squared (nc_hicosmo_E2()), $c_s^2$ the speed of sound (nc_hicosmo_cs2()),
 * $N$ is the lapse function that in this case (using $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$
 * is the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and $$\Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_k}.$$
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_iadiab.h"
#include "perturbations/nc_hipert_adiab.h"

#include "math/cvode_util.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 
#include <gsl/gsl_roots.h>

G_DEFINE_TYPE (NcHIPertAdiab, nc_hipert_adiab, NC_TYPE_HIPERT);

enum {
  PROP_0,
  PROP_SIZE,
};

static gdouble _nc_hipert_adiab_wkb_phase (gdouble y, gdouble x, gpointer userdata);

static void
nc_hipert_adiab_init (NcHIPertAdiab *pa)
{
  NcmSpline *s     = ncm_spline_cubic_notaknot_new ();
  pa->wkb_phase    = ncm_ode_spline_new (s, &_nc_hipert_adiab_wkb_phase);

  ncm_ode_spline_set_abstol (pa->wkb_phase, 0.0);
  ncm_spline_free (s);
}

static void
nc_hipert_adiab_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_adiab_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_adiab_finalize (GObject *object)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);

  ncm_ode_spline_clear (&pa->wkb_phase);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->finalize (object);
}

static void
nc_hipert_adiab_class_init (NcHIPertAdiabClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &nc_hipert_adiab_set_property;
  object_class->get_property = &nc_hipert_adiab_get_property;
  object_class->finalize     = &nc_hipert_adiab_finalize;
}

/**
 * nc_hipert_adiab_new:
 * 
 * Creates a new #NcHIPertAdiab object.
 * 
 * Returns: (transfer full): a new #NcHIPertAdiab.
 */
NcHIPertAdiab *
nc_hipert_adiab_new (void)
{
  NcHIPertAdiab *pa = g_object_new (NC_TYPE_HIPERT_ADIAB, 
                                        NULL);

  return pa;
}

/**
 * nc_hipert_adiab_ref:
 * @pa: a #NcHIPertAdiab.
 * 
 * Increases the reference count of @pa.
 * 
 * Returns: (transfer full): @pa. 
 */
NcHIPertAdiab *
nc_hipert_adiab_ref (NcHIPertAdiab *pa)
{
  return g_object_ref (pa);
}

/**
 * nc_hipert_adiab_free:
 * @pa: a #NcHIPertAdiab.
 * 
 * Decreases the reference count of @pa.
 * 
 */
void 
nc_hipert_adiab_free (NcHIPertAdiab *pa)
{
  g_object_unref (pa);
}

/**
 * nc_hipert_adiab_clear:
 * @pa: a #NcHIPertAdiab.
 * 
 * Decreases the reference count of *@pa and sets *@pa to NULL.
 * 
 */
void 
nc_hipert_adiab_clear (NcHIPertAdiab **pa)
{
  g_clear_object (pa);
}

typedef struct _NcHIPertAdiabWKBPhaseArg
{
  NcHICosmo *cosmo;
  gdouble k;
} NcHIPertAdiabWKBPhaseArg;

static gdouble 
_nc_hipert_adiab_wkb_phase (gdouble y, gdouble x, gpointer userdata)
{
  NcHIPertAdiabWKBPhaseArg *arg = (NcHIPertAdiabWKBPhaseArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), x, arg->k);
  return sqrt (nuA2);
}

/**
 * nc_hipert_adiab_prepare_wkb:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_adiab_prepare_wkb (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (pa);
  if (!pert->prepared)
  {
    NcHIPertAdiabWKBPhaseArg arg;

    arg.cosmo = cosmo;
    arg.k     = pert->k;

    ncm_ode_spline_set_interval (pa->wkb_phase, 2.0 * M_PI, alpha_i, alpha_f);

    ncm_ode_spline_prepare (pa->wkb_phase, &arg);

    pert->prepared = TRUE;
  }
}

/**
 * nc_hipert_adiab_wkb_zeta:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiab_wkb_zeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double zeta;
  gdouble int_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2); 

  g_assert (pert->prepared);

  int_nuA = ncm_spline_eval (pa->wkb_phase->s, alpha) + M_PI * 0.25; 
  
  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->m * nuA);

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);
}

/**
 * nc_hipert_adiab_wkb_zeta_Pzeta:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * @Re_Pzeta: (out caller-allocates): Real part of the wkb solution momentum.
 * @Im_Pzeta: (out caller-allocates): Imaginary part of the wkb solution momentum.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ and its momentum for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiab_wkb_zeta_Pzeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double zeta, Pzeta;
  gdouble int_nuA;
  gdouble dmzetanuA_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2); 

  g_assert (pert->prepared);

  int_nuA = ncm_spline_eval (pa->wkb_phase->s, alpha) + M_PI * 0.25;

  dmzetanuA_nuA = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);

  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->m * nuA);

  Pzeta = -I * cexp (-I * int_nuA) * sqrt (0.5 * eom->m * nuA) - 0.5 * dmzetanuA_nuA * zeta;

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);

  *Re_Pzeta = creal (Pzeta);
  *Im_Pzeta = cimag (Pzeta);  
}

static gdouble 
_nc_hipert_adiab_wkb_nuA2 (gdouble x, gpointer userdata)
{
  NcHIPertAdiabWKBPhaseArg *arg = (NcHIPertAdiabWKBPhaseArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), x, arg->k);
  return nuA2;
}

static gint
_nc_hipert_adiab_wkb_nu2_root (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha0, gdouble *alpha)
{
  NcHIPert *pert = NC_HIPERT (pa);
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble prec = pert->reltol, alpha1 = *alpha;
  NcHIPertAdiabWKBPhaseArg arg;

  arg.cosmo = cosmo;
  arg.k     = pert->k;

  F.function = &_nc_hipert_adiab_wkb_nuA2;
  F.params   = &arg;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, alpha0, alpha1);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      gsl_root_fsolver_free (s);
      return status;
    }

    *alpha = gsl_root_fsolver_root (s);
    alpha0 = gsl_root_fsolver_x_lower (s);
    alpha1 = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (alpha0, alpha1, 0, prec);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return 0;
}

/**
 * nc_hipert_adiab_wkb_maxtime:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$. 
 * 
 * Returns: the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_adiab_wkb_maxtime (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  guint ret = _nc_hipert_adiab_wkb_nu2_root (pa, cosmo, alpha0, &alpha);

  if (ret == 0)
    return alpha;
  else
    return GSL_NAN;
}

/**
 * nc_hipert_adiab_wkb_v:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_v: (out caller-allocates): Real part of the wkb solution.
 * @Im_v: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $v_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 * Returns: the value of $v_\text{WKB}(k, \alpha)$.
 */
void 
nc_hipert_adiab_wkb_v (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_v, gdouble *Im_v)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double v;
  gdouble int_theta;

  g_assert (pert->prepared);

  int_theta = ncm_spline_eval (pa->wkb_phase->s, alpha) + M_PI * 0.25; 
  
  v = cexp (-I * int_theta) / sqrt (2.0 * sqrt (eom->nu2));

  *Re_v = creal (v);
  *Im_v = cimag (v);
}

typedef struct _NcHIPertAdiabArg
{
  NcHICosmo *cosmo;
  NcHIPertAdiab *pa;
} NcHIPertAdiabArg;

static gint
_nc_hipert_adiab_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) f_data;
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), x, NC_HIPERT (arg->pa)->k);
  
  NV_Ith_S (ydot, NC_HIPERT_ADIAB_RE_ZETA) = NV_Ith_S (y, NC_HIPERT_ADIAB_RE_PZETA) / eom->m;
  NV_Ith_S (ydot, NC_HIPERT_ADIAB_IM_ZETA) = NV_Ith_S (y, NC_HIPERT_ADIAB_IM_PZETA) / eom->m;

  NV_Ith_S (ydot, NC_HIPERT_ADIAB_RE_PZETA) = - eom->m * eom->nu2 * NV_Ith_S (y, NC_HIPERT_ADIAB_RE_ZETA);
  NV_Ith_S (ydot, NC_HIPERT_ADIAB_IM_PZETA) = - eom->m * eom->nu2 * NV_Ith_S (y, NC_HIPERT_ADIAB_IM_ZETA);

  return 0;
}

static gint
_nc_hipert_adiab_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) jac_data;
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), t, NC_HIPERT (arg->pa)->k);

  DENSE_ELEM (J, NC_HIPERT_ADIAB_RE_ZETA, NC_HIPERT_ADIAB_RE_PZETA) = 1.0 / eom->m;
  DENSE_ELEM (J, NC_HIPERT_ADIAB_IM_ZETA, NC_HIPERT_ADIAB_IM_PZETA) = 1.0 / eom->m;
  DENSE_ELEM (J, NC_HIPERT_ADIAB_RE_PZETA, NC_HIPERT_ADIAB_RE_ZETA) = - eom->m * eom->nu2;
  DENSE_ELEM (J, NC_HIPERT_ADIAB_IM_PZETA, NC_HIPERT_ADIAB_IM_ZETA) = - eom->m * eom->nu2;

  return 0;
}

/**
 * nc_hipert_adiab_set_init_cond:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alphai: the log-redshift time.
 * @Re_zeta: Real part of the wkb solution.
 * @Im_zeta: Imaginary part of the wkb solution.
 * @Re_Pzeta: Real part of the wkb solution momentum.
 * @Im_Pzeta: Imaginary part of the wkb solution momentum.
 * 
 * Sets the initial conditions for the zeta evolution. 
 * 
 */
void 
nc_hipert_adiab_set_init_cond (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphai, gdouble Re_zeta, gdouble Im_zeta, gdouble Re_Pzeta, gdouble Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  gint flag;

  pert->alpha0 = alphai;

  if (pert->y == NULL)
    pert->y = N_VNew_Serial (4);

  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_ZETA) = Re_zeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_ZETA) = Im_zeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_PZETA) = Re_Pzeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_PZETA) = Im_Pzeta;

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_adiab_f, alphai, pert->y);
    CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alphai, pert->y);
    CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSStolerances (pert->cvode, pert->reltol, pert->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, 1000000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, 4);
  CVODE_CHECK(&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_adiab_J);
  CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, );  
}

/**
 * nc_hipert_adiab_set_init_cond_wkb:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alphai: the log-redshift time.
 * 
 * Sets the initial conditions for the zeta evolution using the value of the WKB solution at @alphai. 
 * 
 */
void 
nc_hipert_adiab_set_init_cond_wkb (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphai)
{
  gdouble Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta;

  nc_hipert_adiab_wkb_zeta_Pzeta (pa, cosmo, alphai, &Re_zeta, &Im_zeta, &Re_Pzeta, &Im_Pzeta);
  nc_hipert_adiab_set_init_cond (pa, cosmo, alphai, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta);
}

/**
 * nc_hipert_adiab_evolve:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alphaf: the final log-redshift time.
 * 
 * Evolve the system until @alphaf.
 * 
 */
void 
nc_hipert_adiab_evolve (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphaf)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertAdiabArg arg;
  gdouble alpha_i = 0.0;
  gint flag;
  
  arg.cosmo = cosmo;
  arg.pa = pa;

  flag = CVodeSetUserData (pert->cvode, &arg);
  CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  flag = CVodeSetUserData (pert->cvode, &arg);
  CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  flag = CVode (pert->cvode, alphaf, pert->y, &alpha_i, CV_NORMAL);
  CVODE_CHECK (&flag, "CVode[nc_hipert_adiab_evolve]", 1, );

  pert->alpha0 = alpha_i;
}

/**
 * nc_hipert_adiab_get_values:
 * @pa: a #NcHIPertAdiab.
 * @alphai: (out caller-allocates): Current time.
 * @Re_zeta: (out caller-allocates): Real part of the solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the solution.
 * @Re_Pzeta: (out caller-allocates): Real part of the solution momentum.
 * @Im_Pzeta: (out caller-allocates): Imaginary part of the solution momentum.
 * 
 * Get the current time and values of the numerical solution.
 * 
 */
void
nc_hipert_adiab_get_values (NcHIPertAdiab *pa, gdouble *alphai, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  
  *alphai = pert->alpha0;
  *Re_zeta  = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_ZETA);
  *Im_zeta  = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_ZETA);
  *Re_Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_PZETA);
  *Im_Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_PZETA);
}
