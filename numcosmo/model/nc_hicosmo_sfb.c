/***************************************************************************
 *            nc_hicosmo_sfb.c
 *
 *  Thu August 04 08:28:24 2022
 *  Copyright  2022  Eduardo José Barroso and Sandro Dias Pinto Vitenti
 *  <eduardo.jsbarroso@uel.br> <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hicosmo_sfb.c
 * Copyright (C) 2022 Eduardo José Barroso and Sandro Dias Pinto Vitenti  <
 * eduardo.jsbarroso@uel.br> <sandro@isoftware.com.br> 
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
 * SECTION:nc_hicosmo_sfb
 * @title: NcHICosmoSFB
 * @short_description:One $w$-fluid model with a bounce phase model.
 * 
 * This model implements the methods from the interface #NcHIPertIAdiab. It is a cosmological model to solve the complex harmonic
 * oscillator problem defined in #NcmCSQ1D, such that the action that described this problem is given by
 * 
 * \begin{align}
 * S &= \int d^3x d\tau \left[ \Pi_\zeta \frac{\partial\zeta}{\partial \tau} - \frac{\Pi_\zeta^2}{4 z^2} + 2 c_s^2 z^2 N \zeta D^2 \zeta   \right]
 * \end{align}
 * ,where $\tau$ is the time defined by $x(\tau) = x_b e^{-|\tau|}$, $x = \frac{a_0}{a}, 
 * $N$ is the lapse function given by $E^{-1}$ for our time choice, $\zeta$ is the curvature perturbation,
 * $\Pi_\zeta$ is its conjugated momentum and $D^2$ is the laplacian operator. Bear in mind that this is
 * the dimensionless action and all the variables are also dimensionless.
 *
 * In this model the adiabatic mode $\zeta$ has its mass, speed of sound square $c_s^2$ and frequency squared $\nu_\zeta^2$ given by
 * \begin{align}
 * m_\zeta &= 2 z^2 = \frac{(1 + w)}{c_s^2 N x^3}, \\\\
 * c_s^2 &= w, \\\\
 * \nu_\zeta^2 &= m_\zeta c_s^2 k^2 N^2 x^2.
 * ,\end{align}
 * where $k$ is from the Fourier expansion of the field.
 *
 * The model is implemented for a universe with null curvature and containing only cold dark matter,
 * such that $w = 0.3$. This model shall be used as the argument of the #NcmCSQ1D class methods.
 * For an example, check NumCosmo/examples/example\_adiab\_spec.py.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_sfb.h"

static void nc_hipert_adiab_interface_init (NcHIPertIAdiabInterface *iface);

G_DEFINE_TYPE_WITH_CODE (NcHICosmoSFB, nc_hicosmo_sfb, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IADIAB,
                                                nc_hipert_adiab_interface_init)
                        );

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_sfb_init (NcHICosmoSFB *sfb)
{
}

static void
nc_hicosmo_sfb_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_sfb_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_sfb_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_sfb_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_sfb_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_sfb_bgp_cs2 (NcHICosmo *cosmo, gdouble tau);

static gdouble _nc_hicosmo_sfb_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_sfb_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_sfb_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_sfb_xb (NcHICosmo *cosmo);

static void
nc_hicosmo_sfb_class_init (NcHICosmoSFBClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);
  
  
  ncm_model_class_set_name_nick (model_class, "SFB", "SFB");
  ncm_model_class_add_params (model_class, NC_HICOSMO_SFB_SPARAM_LEN, 0, PROP_SIZE);
  
  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_H0, "H_0", "H0",
                              10.0, 500.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_OMEGA_W, "\\Omega_{w0}", "Omegaw",
                              1e-8,  10.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_OMEGA_W,
                              NCM_PARAM_TYPE_FREE);
  /* Set w param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_W, "w", "w",
                              1e-50,  10.0, 1.0e-8,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_W,
                              NCM_PARAM_TYPE_FIXED);
  /* Set xb param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_X_B, "x_b", "xb",
                              1.0e-50,  1.0, 1.0e25,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_X_B,
                              NCM_PARAM_TYPE_FIXED);
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
  
  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_sfb_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_sfb_E2);
  nc_hicosmo_set_Omega_c0_impl  (parent_class, &_nc_hicosmo_sfb_Omega_c0);
  nc_hicosmo_set_Omega_t0_impl  (parent_class, &_nc_hicosmo_sfb_Omega_t0);
  nc_hicosmo_set_xb_impl        (parent_class, &_nc_hicosmo_sfb_xb);
  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_sfb_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_sfb_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl   (parent_class, &_nc_hicosmo_sfb_bgp_cs2);
  
  object_class->finalize = nc_hicosmo_sfb_finalize;
}

static gdouble _nc_hicosmo_sfb_iadiab_eval_m   (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_nu  (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_xi  (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_F1  (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_F2  (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_H   (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_sfb_iadiab_eval_x   (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static void _nc_hicosmo_sfb_iadiab_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *xi, gdouble *F1);

static void
nc_hipert_adiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_m      = &_nc_hicosmo_sfb_iadiab_eval_m;
  iface->eval_mnu    = &_nc_hicosmo_sfb_iadiab_eval_mnu;
  iface->eval_nu     = &_nc_hicosmo_sfb_iadiab_eval_nu;
  iface->eval_xi     = &_nc_hicosmo_sfb_iadiab_eval_xi;
  iface->eval_F1     = &_nc_hicosmo_sfb_iadiab_eval_F1;
  iface->eval_F2     = &_nc_hicosmo_sfb_iadiab_eval_F2;
  iface->eval_H      = &_nc_hicosmo_sfb_iadiab_eval_H;
  iface->eval_x      = &_nc_hicosmo_sfb_iadiab_eval_x;
  iface->eval_system = &_nc_hicosmo_sfb_iadiab_eval_system;
}

#define VECTOR   (NCM_MODEL (cosmo)->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_H0))
#define OMEGA_W  (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_OMEGA_W))
#define W        (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_W))
#define X_B      (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_X_B))

/****************************************************************************
 * Future Implementation if Needed to compute Z and rho+p for more than one fluid
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble lnX_B = log (X_B);
  const gdouble x = z + 1.0;
  const gdouble lnx = log(x);
  const gdouble w = W;
  const gdouble x_3_1pw = pow (x, 3.0 * (1.0 + w));
  const gdouble x_2_p3w = pow (x, 2.0 + 3.0 * w);
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * (lnx - lnX_B));
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt(E2);
  const gdouble e_3t1w = exp (-3.0 *(1.0 - w) * (lnx - lnX_B));

  return 2.0 * E * Omega_w * x_2_p3w * (3.0 * (1.0 + w) - 6.0 * w * e_3t1w);

}

static gdouble
_nc_hicosmo_sfb_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble lnX_B = log (X_B);
  const gdouble x = z + 1.0;
  const gdouble lnx = log(x);
  const gdouble w = W;
  const gdouble x_3_1pw = pow (x, 3.0 * (1.0 + w));
  const gdouble x_2_p3w = pow (x, 2.0 + 3.0 * w);
  const gdouble x_1_p3w = pow (x, 1.0 + 3.0 * w);
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * (lnx - lnX_B));
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt(E2);
  const gdouble e_3t1w = exp (-3.0 *(1.0 - w) * (lnx - lnX_B));
  const gdouble dEdz = Omega_w * x_2_p3w * (3.0 * (1.0 + w) - 6.0 * w * e_3t1w);
  const gdouble dEdz2 = Omega_w * x_1_p3w * (3.0 * (1.0 + w) * (3.0 * w + 2.0) - 6.0 * w * (6.0 * w + 1.0) * e_3t1w);
  
  return 2.0 * dEdz * dEdz + 2.0 * E * dEdz2; 
}

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_sfb_E2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble lnX_B = log (X_B);
  const gdouble x = z + 1.0;
  const gdouble lnx = log(x);
  const gdouble w = W;
  const gdouble x_3_1pw = pow (x, 3.0 * (1.0 + w));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * (lnx - lnX_B));

  return -Omega_w *x_3_1pw * e_3t_1mw;
}


/****************************************************************************
 * Speed of sound squared
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_bgp_cs2 (NcHICosmo *cosmo, gdouble tau)
{
  gdouble w = W;
  
  
  return w;
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0 * (OMEGA_W);
}

static gdouble
_nc_hicosmo_sfb_Omega_t0 (NcHICosmo *cosmo)
{
  return 1.0;
}

static gdouble
_nc_hicosmo_sfb_Omega_c0 (NcHICosmo *cosmo)
{
  return OMEGA_W;
}

static gdouble
_nc_hicosmo_sfb_xb (NcHICosmo *cosmo)
{
  return X_B;
}

/*****************************************************************************
 *  Interface functions
 ******************************************************************************/
static gdouble
_nc_hicosmo_sfb_iadiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble tabs = fabs (tau);
  const gdouble lnX_B = log (X_B);
  const gdouble x = exp(lnX_B - tabs);

  return x;
}

static gdouble
_nc_hicosmo_sfb_iadiab_eval_H (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);

  return E;
}

static gdouble
_nc_hicosmo_sfb_iadiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble mult = 3.0 * (1.0 + w) / (2.0 * w);
  const gdouble x = exp(lnX_B - tabs);
  const gdouble x3   = x * x * x;
  const gdouble z2 = mult * 1.0 / (x3 * N);

  return z2 * 2.0;
 } 
 
static gdouble
_nc_hicosmo_sfb_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble x = exp(lnX_B - tabs);

  return sqrt (w) * k * N * x;
 }

static gdouble
_nc_hicosmo_sfb_iadiab_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
 {
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);  
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble x = exp(lnX_B - tabs);
  const gdouble mult = 3.0 * (1.0 + w) / (2.0 * w);
  const gdouble x3   = x * x * x;
  const gdouble z2 = mult * 1.0 / (x3 * N);
  const gdouble m = 2.0 * z2;
  const gdouble nu = sqrt (w) * k * N * x;

  return m * nu;
 }
 
static gdouble
_nc_hicosmo_sfb_iadiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);  
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble x = exp(lnX_B - tabs);
  const gdouble mult = 3.0 * (1.0 + w) / (2.0 * w);
  const gdouble x3   = x * x * x;
  const gdouble z2 = mult * 1.0 / (x3 * N);
  const gdouble m = 2.0 * z2;
  const gdouble nu = sqrt (w) * k * N * x;
  const gdouble mnu =  m * nu;

  return log (mnu);
}

static gdouble
_nc_hicosmo_sfb_iadiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble sign = GSL_SIGN(tau);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs); 
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble x = exp(lnX_B - tabs);
  const gdouble nu = sqrt (w) * k * N * x;

  return  1.0 * sign / nu;
 }
static gdouble
_nc_hicosmo_sfb_iadiab_eval_F2 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble sign = GSL_SIGN(tau);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_t3 = exp (-3.0 *(1.0 - w) * tabs);
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble x = exp(lnX_B - tabs);
  const gdouble nu = sqrt (w) * k * N * x;
  const gdouble dN_dtau = 1.0 / 2.0 * gsl_pow_3 (N) * Omega_w * x_3_1pw * sign * (-3.0 * (1.0 + w) + 6.0 * w * e_t3);
  const gdouble N2 = N * N;

  return 1.0 /  (2.0 * nu * k * sqrt(w)) * (1.0  /(x * N) - dN_dtau / (x * N2));
}

static void
_nc_hicosmo_sfb_iadiab_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *xi, gdouble *F1)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  const gdouble sign = GSL_SIGN(tau);
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (3.0 * (1.0 + w) * (-tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble e_3t_1mw = expm1 (-3.0 *(1.0 - w) * tabs);
  const gdouble E2 = -Omega_w *x_3_1pw * e_3t_1mw;
  const gdouble E = sqrt (E2);
  const gdouble N = 1.0 / E;
  const gdouble x = exp(lnX_B - tabs);
  const gdouble mult = 3.0 * (1.0 + w) / (2.0 * w);
  const gdouble x3   = x * x * x;
  const gdouble z2 = mult * 1.0 / (x3 * N);
  const gdouble m = 2.0 * z2;
  const gdouble _nu = sqrt (w) * k * N * x;
  const gdouble mnu =  m * _nu;
  
  xi[0] = log(mnu); 
  F1[0] = 1.0 * sign / _nu;
  nu[0] = _nu;
}

/**
 * nc_hicosmo_sfb_new:
 *
 * Creates a new #NcHICosmoSFB object.
 *
 * Returns: (transfer full): a #NcHICosmoSFB model object.
 */
NcHICosmoSFB *
nc_hicosmo_sfb_new (void)
{
  NcHICosmoSFB *sfb = g_object_new (NC_TYPE_HICOSMO_SFB, NULL);
  
  
  return sfb;
}

