/***************************************************************************
 *            nc_hipert_adiabatic.c
 *
 *  Tue June 03 17:20:42 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiabatic.c
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
 * SECTION:nc_hipert_adiabatic
 * @title: Adiabatic Perturbation Object
 * @short_description: Perturbation object for adiabatic mode only 
 *
 * This object provides the computation of the adiabatic mode for the cosmological 
 * perturbations. It solves the equation of motion for the gauge invariant variable
 * (see <link linkend="XVitenti2013">Vitenti (2013)</link> for notation and details)
 * $$\zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + \frac{\bar{\Theta}}{3}\hat{\mathcal{V}}.$$
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{\Pi_\zeta}{m_\zeta}, \\\\
 * \Pi_\zeta^\prime &= m_\zeta\mu_\zeta^2\zeta.
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
#include "nc_hipert_adiabatic.h"

#include "math/cvode_util.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 

G_DEFINE_TYPE (NcHIPertAdiabatic, nc_hipert_adiabatic, NC_TYPE_HIPERT);

enum {
  PROP_0,
  PROP_K,
  PROP_RELTOL,
  PROP_SIZE,
};

static gdouble _nc_hipert_adiabatic_wkb_phase (gdouble y, gdouble x, gpointer userdata);

static void
nc_hipert_adiabatic_init (NcHIPertAdiabatic *pa)
{
  NcmSpline *s     = ncm_spline_cubic_notaknot_new ();
  pa->wkb_phase    = ncm_ode_spline_new (s, &_nc_hipert_adiabatic_wkb_phase);

  ncm_ode_spline_set_abstol (pa->wkb_phase, 0.0);
  ncm_spline_free (s);
}

static void
nc_hipert_adiabatic_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiabatic *pa = NC_HIPERT_ADIABATIC (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIABATIC (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_adiabatic_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiabatic *pa = NC_HIPERT_ADIABATIC (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIABATIC (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_adiabatic_finalize (GObject *object)
{
  /* NcHIPertAdiabatic *pa = NC_HIPERT_ADIABATIC (object); */
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiabatic_parent_class)->finalize (object);
}

static void
nc_hipert_adiabatic_class_init (NcHIPertAdiabaticClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &nc_hipert_adiabatic_set_property;
  object_class->get_property = &nc_hipert_adiabatic_get_property;
  object_class->finalize     = &nc_hipert_adiabatic_finalize;
}

/**
 * nc_hipert_adiabatic_new:
 * 
 * Creates a new #NcHIPertAdiabatic object.
 * 
 * Returns: (transfer full): a new #NcHIPertAdiabatic.
 */
NcHIPertAdiabatic *
nc_hipert_adiabatic_new (void)
{
  NcHIPertAdiabatic *pa = g_object_new (NC_TYPE_HIPERT_ADIABATIC, 
                                        NULL);

  return pa;
}

/**
 * nc_hipert_adiabatic_ref:
 * @pa: a #NcHIPertAdiabatic.
 * 
 * Increases the reference count of @pa.
 * 
 * Returns: (transfer full): @pa. 
 */
NcHIPertAdiabatic *
nc_hipert_adiabatic_ref (NcHIPertAdiabatic *pa)
{
  return g_object_ref (pa);
}

/**
 * nc_hipert_adiabatic_free:
 * @pa: a #NcHIPertAdiabatic.
 * 
 * Decreases the reference count of @pa.
 * 
 */
void 
nc_hipert_adiabatic_free (NcHIPertAdiabatic *pa)
{
  g_object_unref (pa);
}

/**
 * nc_hipert_adiabatic_clear:
 * @pa: a #NcHIPertAdiabatic.
 * 
 * Decreases the reference count of *@pa and sets *@pa to NULL.
 * 
 */
void 
nc_hipert_adiabatic_clear (NcHIPertAdiabatic **pa)
{
  g_clear_object (pa);
}

typedef struct _NcHIPertAdiabaticWKBPhaseArg
{
  NcHICosmo *cosmo;
  gdouble k;
} NcHIPertAdiabaticWKBPhaseArg;

static gdouble 
_nc_hipert_adiabatic_wkb_phase (gdouble y, gdouble x, gpointer userdata)
{
  NcHIPertAdiabaticWKBPhaseArg *arg = (NcHIPertAdiabaticWKBPhaseArg *) userdata;
  const gdouble theta = nc_hicosmo_wkb_adiab_theta (arg->cosmo, x, arg->k);

  /* printf ("# % 20.15g % 20.15g % 20.15g % 20.15g\n", arg->k, x, nc_hicosmo_x_alpha (arg->cosmo, x), theta); */

  return theta;
}

/**
 * nc_hipert_adiabatic_prepare_wkb:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_adiabatic_prepare_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (pa);
  if (!pert->prepared)
  {
    NcHIPertAdiabaticWKBPhaseArg arg;

    arg.cosmo = cosmo;
    arg.k     = pert->k;

    ncm_ode_spline_set_interval (pa->wkb_phase, 2.0 * M_PI, alpha_i, alpha_f);

    ncm_ode_spline_prepare (pa->wkb_phase, &arg);

    pert->prepared = TRUE;
  }
}

/**
 * nc_hipert_adiabatic_mass_zeta:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * 
 * Computes the effective mass of the $\zeta$ equation of motion. 
 * 
 * Returns: $m_\zeta(\alpha)$.
 */
gdouble 
nc_hipert_adiabatic_mass_zeta (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble alpha)
{
  
  return 0.0;
}

/**
 * nc_hipert_adiabatic_freq2_zeta:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @k: the mode $k$.
 * @alpha: the log-redshift time.
 * 
 * Computes the effective frequency squared of the $\zeta$ equation of motion. 
 * 
 * Returns: $\mu_\zeta(\alpha)$.
 */
gdouble 
nc_hipert_adiabatic_freq2_zeta (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha)
{
  
  return 0.0;
}

/**
 * nc_hipert_adiabatic_zeta_wkb:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @k: the mode $k$.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiabatic_zeta_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHICosmoEOMAdiabZeta *eom = nc_hicosmo_adiabatic_zeta (cosmo, alpha, k);
  complex double zeta;
  gdouble int_theta;

  g_assert_cmpfloat (k, ==, pert->k);
  g_assert (pert->prepared);

  int_theta = ncm_spline_eval (pa->wkb_phase->s, alpha); 
  
  zeta = cexp (-I * int_theta) / sqrt (2.0 * eom->m * sqrt (eom->mu2));

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);
}

/**
 * nc_hipert_adiabatic_zeta_Pzeta_wkb:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @k: the mode $k$.
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
nc_hipert_adiabatic_zeta_Pzeta_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHICosmoEOMAdiabZeta *eom = nc_hicosmo_adiabatic_zeta (cosmo, alpha, k);
  complex double zeta, Pzeta;
  gdouble int_theta;
  gdouble dmtheta_theta;

  g_assert_cmpfloat (k, ==, pert->k);
  g_assert (pert->prepared);

  int_theta = ncm_spline_eval (pa->wkb_phase->s, alpha); 

  dmtheta_theta = nc_hicosmo_wkb_adiab_dmtheta (cosmo, alpha, k);
  
  zeta = cexp (-I * int_theta) / sqrt (2.0 * eom->m * sqrt (eom->mu2));

  Pzeta = -I * cexp (-I * int_theta) * sqrt (0.5 * eom->m * sqrt (eom->mu2)) - 0.5 * dmtheta_theta * zeta;
  
  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);

  *Re_Pzeta = creal (Pzeta);
  *Im_Pzeta = cimag (Pzeta);  
}

/**
 * nc_hipert_adiabatic_v_wkb:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @k: the mode $k$.
 * @alpha: the log-redshift time.
 * @Re_v: (out caller-allocates): Real part of the wkb solution.
 * @Im_v: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $v_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 * Returns: the value of $v_\text{WKB}(k, \alpha)$.
 */
void 
nc_hipert_adiabatic_v_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha, gdouble *Re_v, gdouble *Im_v)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHICosmoEOMAdiabZeta *eom = nc_hicosmo_adiabatic_zeta (cosmo, alpha, k);
  complex double v;
  gdouble int_theta;

  g_assert_cmpfloat (k, ==, pert->k);
  g_assert (pert->prepared);

  int_theta = ncm_spline_eval (pa->wkb_phase->s, alpha); 
  
  v = cexp (-I * int_theta) / sqrt (2.0 * sqrt (eom->mu2));

  *Re_v = creal (v);
  *Im_v = cimag (v);
}

typedef struct _NcHIPertAdiabaticArg
{
  NcHICosmo *cosmo;
  NcHIPertAdiabatic *pa;
} NcHIPertAdiabaticArg;

static gint
_nc_hipert_adiabatic_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertAdiabaticArg *arg = (NcHIPertAdiabaticArg *) f_data;
  NcHICosmoEOMAdiabZeta *eom = nc_hicosmo_adiabatic_zeta (arg->cosmo, x, NC_HIPERT (arg->pa)->k);
  
  NV_Ith_S (ydot, NC_HIPERT_ADIABATIC_RE_ZETA) = NV_Ith_S (y, NC_HIPERT_ADIABATIC_RE_PZETA) / eom->m;
  NV_Ith_S (ydot, NC_HIPERT_ADIABATIC_IM_ZETA) = NV_Ith_S (y, NC_HIPERT_ADIABATIC_IM_PZETA) / eom->m;

  NV_Ith_S (ydot, NC_HIPERT_ADIABATIC_RE_PZETA) = - eom->m * eom->mu2 * NV_Ith_S (y, NC_HIPERT_ADIABATIC_RE_ZETA);
  NV_Ith_S (ydot, NC_HIPERT_ADIABATIC_IM_PZETA) = - eom->m * eom->mu2 * NV_Ith_S (y, NC_HIPERT_ADIABATIC_IM_ZETA);
  
  return 0;
}

static gint
_nc_hipert_adiabatic_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertAdiabaticArg *arg = (NcHIPertAdiabaticArg *) jac_data;
  NcHICosmoEOMAdiabZeta *eom = nc_hicosmo_adiabatic_zeta (arg->cosmo, t, NC_HIPERT (arg->pa)->k);

  DENSE_ELEM (J, NC_HIPERT_ADIABATIC_RE_ZETA, NC_HIPERT_ADIABATIC_RE_PZETA) = 1.0 / eom->m;
  DENSE_ELEM (J, NC_HIPERT_ADIABATIC_IM_ZETA, NC_HIPERT_ADIABATIC_IM_PZETA) = 1.0 / eom->m;
  DENSE_ELEM (J, NC_HIPERT_ADIABATIC_RE_PZETA, NC_HIPERT_ADIABATIC_RE_ZETA) = - eom->m * eom->mu2;
  DENSE_ELEM (J, NC_HIPERT_ADIABATIC_IM_PZETA, NC_HIPERT_ADIABATIC_IM_ZETA) = - eom->m * eom->mu2;

  return 0;
}

/**
 * nc_hipert_adiabatic_set_init_cond:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @k: the mode $k$.
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
nc_hipert_adiabatic_set_init_cond (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alphai, gdouble Re_zeta, gdouble Im_zeta, gdouble Re_Pzeta, gdouble Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  gint flag;

  nc_hipert_set_mode_k (pert, k);

  pert->alpha0 = alphai;

  if (pert->y == NULL)
    pert->y = N_VNew_Serial (4);

  NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_RE_ZETA) = Re_zeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_IM_ZETA) = Im_zeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_RE_PZETA) = Re_Pzeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_IM_PZETA) = Im_Pzeta;

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_adiabatic_f, alphai, pert->y);
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

  flag = CVDense (pert->cvode, 2);
  CVODE_CHECK(&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_adiabatic_J);
  CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, );  
}

/**
 * nc_hipert_adiabatic_set_init_cond_wkb:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @k: the mode $k$.
 * @alphai: the log-redshift time.
 * 
 * Sets the initial conditions for the zeta evolution using the value of the WKB solution at @alphai. 
 * 
 */
void 
nc_hipert_adiabatic_set_init_cond_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alphai)
{
  gdouble Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta;
  nc_hipert_adiabatic_zeta_Pzeta_wkb (pa, cosmo, k, alphai, &Re_zeta, &Im_zeta, &Re_Pzeta, &Im_Pzeta);
  nc_hipert_adiabatic_set_init_cond (pa, cosmo, k, alphai, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta);
}

/**
 * nc_hipert_adiabatic_evolve:
 * @pa: a #NcHIPertAdiabatic.
 * @cosmo: a #NcHICosmo.
 * @alphaf: the final log-redshift time.
 * 
 * Evolve the system until @alphaf.
 * 
 */
void 
nc_hipert_adiabatic_evolve (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble alphaf)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertAdiabaticArg arg;
  gdouble alpha_i = 0.0;
  gint flag;
  
  arg.cosmo = cosmo;
  arg.pa = pa;

  flag = CVodeSetUserData (pert->cvode, &arg);
  CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  flag = CVodeSetUserData (pert->cvode, &arg);
  CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  flag = CVode (pert->cvode, alphaf, pert->y, &alpha_i, CV_NORMAL);
  CVODE_CHECK (&flag, "CVode[nc_hipert_adiabatic_evolve]", 1, );

  pert->alpha0 = alpha_i;
}

/**
 * nc_hipert_adiabatic_get_values:
 * @pa: a #NcHIPertAdiabatic.
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
nc_hipert_adiabatic_get_values (NcHIPertAdiabatic *pa, gdouble *alphai, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  
  *alphai = pert->alpha0;
  *Re_zeta  = NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_RE_ZETA);
  *Im_zeta  = NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_IM_ZETA);
  *Re_Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_RE_PZETA);
  *Im_Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ADIABATIC_IM_PZETA);
}
