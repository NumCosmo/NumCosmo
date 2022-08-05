/***************************************************************************
 *            nc_hicosmo_sfb.c
 *
 *  Wed June 04 10:04:24 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hicosmo_sfb.c
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
 * SECTION:nc_hicosmo_sfb
 * @title: NcHICosmoSFB
 * @short_description:One $w$-fluid model with a bounce phase model.
 *
 * In this model the adiabatic mode $\zeta$ has its mass, speed of sound square $c_s^2$ and frequency square $\nu_\zeta^2$ given by
 * \begin{align}
 * m_\zeta &= 2 z^2 = 2\sqrt{\frac{3(1+w)}{2w}}a, \\\\
 * c_s^2 &= w, \\\\
 * \nu_\zeta^2 &= \frac {m_\zeta c_s^2 k^2}{a^2}.
 * \end{align}
 *
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

static gdouble _nc_hicosmo_sfb_eval_z (NcHICosmo *cosmo, gdouble tau);
static gdouble _nc_hicosmo_sfb_x (NcHICosmo *cosmo, gdouble tau);
static gdouble _nc_hicosmo_sfb_N_dtau (NcHICosmo *cosmo, gdouble tau);
static gdouble _nc_hicosmo_sfb_N_dtau2 (NcHICosmo *cosmo, gdouble tau);
static gdouble _nc_hicosmo_sfb_E2 (NcHICosmo *cosmo, gdouble tau);
static gdouble _nc_hicosmo_sfb_dE2_dz (NcHICosmo *cosmo, gdouble tau);
static gdouble _nc_hicosmo_sfb_d2E2_dz2 (NcHICosmo *cosmo, gdouble tau);
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
  /* Set Omega_r0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_OMEGA_R, "\\Omega_{r0}", "Omegar",
                              1e-8,  10.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_OMEGA_R,
                              NCM_PARAM_TYPE_FREE);
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
  /* Set taub param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_TAU_B, "tau_b", "taub",
                              1.0e-10,  1.0e50, 1.0e-25,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_TAU_B,
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

static gdouble _nc_hipert_iadiab_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hipert_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hipert_iadiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hipert_iadiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hipert_iadiab_eval_F2 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static void _nc_hipert_iadiab_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu);

static void
nc_hipert_adiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_mnu    = &_nc_hipert_iadiab_eval_mnu;
  iface->eval_nu     = &_nc_hipert_iadiab_eval_nu;
  iface->eval_xi     = &_nc_hipert_iadiab_eval_xi;
  iface->eval_F1     = &_nc_hipert_iadiab_eval_F1;
  iface->eval_F2     = &_nc_hipert_iadiab_eval_F2;
  iface->eval_system = &_nc_hipert_iadiab_eval_system;
}

#define VECTOR   (NCM_MODEL (cosmo)->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_H0))
#define OMEGA_R  (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_OMEGA_R))
#define OMEGA_W  (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_OMEGA_W))
#define W        (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_W))
#define X_B      (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_X_B))
#define TAU_B    (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_TAU_B))

/****************************************************************************
 * Future Implementation if Needed to compute Z and rho+p for more than one fluid
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_dE2_dz (NcHICosmo *cosmo, gdouble tau)
{
  return 0.0;
}

static gdouble
_nc_hicosmo_sfb_d2E2_dz2 (NcHICosmo *cosmo, gdouble tau)
{
  return 0.0;
}

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_sfb_E2 (NcHICosmo *cosmo, gdouble tau)
{
  const gdouble lnX_B = log (X_B);
  const gdouble tabs = fabs (tau);
  const gdouble x = exp (-tabs + lnX_B);
  const gdouble w = W;
  const gdouble x_3_1pw = exp (-3 * (1.0 + w) * (tabs + lnX_B));
  const gdouble Omega_w = OMEGA_W;
  const gdouble exp_t, et_3w;
  
  
  return -Omega_w *x_3_1pw *expm1 (-3.0 *(1.0 - w) *(tabs + lnX_B));
}

/****************************************************************************
 * d/dtau of the lapse function N
 ****************************************************************************/

static gdouble
_nc_hicosmo_sfb_N_dtau (NcHICosmo *cosmo, gdouble tau)
{
  const gdouble sign = GSL_SIGN (tau);
  const gdouble tabs = fabs (tau);
  const gdouble E = sqrt (_nc_hicosmo_sfb_E2 (cosmo, tau));
  const gdouble N = 1.0 / E;
  const gdouble x = _nc_hicosmo_sfb_x (cosmo, tau);
  const gdouble w = W;
  const gdouble x_3w = pow (x, 3 * (1 + w));
  const gdouble omega_w = OMEGA_W;
  const gdouble exp_t, e_3w;

  if (tabs <= 1e-5)
  {
    exp_t = expm1 (tabs) + 1.0;
  }
  else
  {
    exp_t = exp (tabs);
  }
  
  e_3w = pow (exp_t, -3 * (1 - w));
  
  /*printf("dNdtau %.20f", 1/2 * pow(N, 3) * omega_w * x_3w * sign * (- 3 * (1 + w) + 6 * w * e_3w));*/
  return 1.0 / 2.0 * gsl_pow_3 (N) * omega_w * x_3w * sign * (-3.0 * (1.0 + w) + 6.0 * w * e_3w);
}

/****************************************************************************
 * d2/dtau2 of the lapse function N
 ****************************************************************************/

static gdouble
_nc_hicosmo_sfb_N_dtau2 (NcHICosmo *cosmo, gdouble tau)
{
  gdouble tabs = fabs (tau);
  gdouble N = pow (_nc_hicosmo_sfb_E2 (cosmo, tau), -1);
  gdouble dNdtau = _nc_hicosmo_sfb_N_dtau (cosmo, tau);
  gdouble x = _nc_hicosmo_sfb_x (cosmo, tau);
  gdouble w = W;
  gdouble x_3w = pow (x, 3 * (1 + w));
  gdouble omega_w = OMEGA_W;
  gdouble exp_t, e_3w;
  
  
  if (tabs <= 1e-5)
  {
    exp_t = expm1 (tabs) + 1.0;
  }
  else
  {
    exp_t = exp (tabs);
  }
  
  e_3w = pow (exp_t, -3 * (1 - w));
  
  /*printf("dndtau2 %.20f", 3/2 * pow(N, 2) * dNdtau + 1/2 * pow(N,3) * omega_w * x_3w * (9 * pow(1 + w, 2) - e_3w * 36.0));*/
  return 3 / 2 * pow (N, 2) * dNdtau + 1 / 2 * pow (N, 3) * omega_w * x_3w * (9 * pow (1 + w, 2) - e_3w * 36.0);
}

/****************************************************************************
 * x function, the scale factor today divided by the scale factor at a time \tau
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_x (NcHICosmo *cosmo, gdouble tau)
{
  gdouble exp_t;
  gdouble tabs = fabs (tau);
  gdouble result;
  
  
  if (tabs <= 1.0e-6)
  {
    exp_t  = expm1 (-tabs) + 1.0;
    result = exp_t * X_B;
  }
  else
  {
    exp_t  = log (X_B);
    result = exp (3.0 * (1.0 + W) * (log (X_B) - tabs));
    printf (" %.20f    ", result);
  }
  
  return result;
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_dE2_dt (NcHICosmo *cosmo, gdouble tau)
{
  gdouble sign;
  
  
  if (tau > 0)
  {
    sign = 1.0;
  }
  else
  {
    sign = -1.0;
  }
  
  gdouble tabs = fabs (tau);
  gdouble x = _nc_hicosmo_sfb_x (cosmo, tau);
  gdouble w = W;
  gdouble x_3w = pow (x, 3 * (1 + w));
  gdouble omega_w = OMEGA_W;
  gdouble exp_t, et_3w;
  gdouble H_tau = sign;
  
  
  if (tabs <= 1e-5)
  {
    exp_t = -expm1 (tau) + 1.0;
  }
  else
  {
    exp_t = exp (tabs);
  }
  
  et_3w = 1.0 - pow (exp_t, -3 * (1 - w));
  
  return 3.0  * omega_w * x_3w * H_tau * et_3w - 3 / 2 * omega_w * (1.0 + w) * x_3w * H_tau;
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

static gdouble
_nc_hicosmo_sfb_eval_z (NcHICosmo *cosmo, gdouble tau)
{
  gdouble x    = _nc_hicosmo_sfb_x (cosmo, tau);
  gdouble w    = W;
  gdouble mult = pow (3 * (1 + w) / (2 * w), 0.5);
  
  
  return mult * pow (x, -1);
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0 * (OMEGA_R + OMEGA_W);
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
_nc_hipert_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  gdouble w        = W;
  gdouble N        = pow (_nc_hicosmo_sfb_E2 (cosmo, tau), -1.0);
  gdouble x        = _nc_hicosmo_sfb_x (cosmo, tau);
  
  
  return sqrt (w) * k * N * x;
}

static gdouble
_nc_hipert_iadiab_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  gdouble nu       = _nc_hipert_iadiab_eval_nu (iad, tau, k);
  gdouble z        = _nc_hicosmo_sfb_eval_z (cosmo, tau);
  gdouble m        = 2 * pow (z, 2);
  
  
  return m * nu;
}

static gdouble
_nc_hipert_iadiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  gdouble mnu = _nc_hipert_iadiab_eval_mnu (iad, tau, k);
  
  
  return log (mnu);
}

static gdouble
_nc_hipert_iadiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  gdouble H        = _nc_hicosmo_sfb_E2 (cosmo, tau);
  gdouble nu       = _nc_hipert_iadiab_eval_nu (iad, tau, k);
  gdouble dNdtau   = _nc_hicosmo_sfb_N_dtau (cosmo, tau);
  
  
  return H / nu + (H * dNdtau) / nu;
}

static gdouble
_nc_hipert_iadiab_eval_F2 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (iad);
  gdouble H        = _nc_hicosmo_sfb_E2 (cosmo, tau);
  gdouble nu       = _nc_hipert_iadiab_eval_nu (iad, tau, k);
  gdouble dNdtau   = _nc_hicosmo_sfb_N_dtau (cosmo, tau);
  gdouble dNdtau2  = _nc_hicosmo_sfb_N_dtau2 (cosmo, tau);
  gdouble nu2_2    = 2.0 * pow (nu, 2);
  
  
  return -((dNdtau * pow (H, 2)) / nu2_2) - ((pow (dNdtau, 2) * pow (H, 2)) / nu2_2) + ((dNdtau2 * pow (H, 2)) / nu2_2);
}

static void
_nc_hipert_iadiab_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu)
{
  printf ("Not Implemented.");
}

/**
 * nc_hicosmo_sfb_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoSFB *
nc_hicosmo_sfb_new (void)
{
  NcHICosmoSFB *sfb = g_object_new (NC_TYPE_HICOSMO_SFB, NULL);
  
  
  return sfb;
}

