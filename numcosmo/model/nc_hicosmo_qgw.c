/***************************************************************************
 *            nc_hicosmo_qgw.c
 *
 *  Tue March 19 10:45:17 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>>
 ****************************************************************************/
/*
 * nc_hicosmo_qgw.c
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_hicosmo_qgw
 * @title: NcHICosmoQGW
 * @short_description: $w$-fluid model with a quantum generated bounce phase model.
 *
 * The NcHICosmoQGW class implements the $w$-fluid model with a quantum generated
 * bounce phase model.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qgw.h"
#include "perturbations/nc_hipert_adiab.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_trig.h>
#endif /* NUMCOSMO_GIR_SCAN */

static void nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface);

struct _NcHICosmoQGWPrivate
{
  guint place_holder;
};

G_DEFINE_TYPE_WITH_CODE (NcHICosmoQGW, nc_hicosmo_qgw, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IADIAB,
                                                nc_hipert_iadiab_interface_init)
                         G_ADD_PRIVATE (NcHICosmoQGW)
                        );

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_qgw_init (NcHICosmoQGW *cosmo_qgw)
{
  NcHICosmoQGWPrivate * const self = cosmo_qgw->priv = nc_hicosmo_qgw_get_instance_private (cosmo_qgw);

  self->place_holder = 0;
}

static void
_nc_hicosmo_qgw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoQGW *qgw = NC_HICOSMO_QGW (object);

  g_return_if_fail (NC_IS_HICOSMO_QGW (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qgw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoQGW *qgw = NC_HICOSMO_QGW (object);

  g_return_if_fail (NC_IS_HICOSMO_QGW (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hicosmo_qgw_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qgw_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_qgw_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgw_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgw_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgw_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgw_bgp_cs2 (NcHICosmo *cosmo, gdouble z);

static gdouble _nc_hicosmo_qgw_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qgw_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qgw_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qgw_xb (NcHICosmo *cosmo);

static void
nc_hicosmo_qgw_class_init (NcHICosmoQGWClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_hicosmo_qgw_set_property;
  model_class->get_property = &_nc_hicosmo_qgw_get_property;
  object_class->finalize    = nc_hicosmo_qgw_finalize;

  ncm_model_class_set_name_nick (model_class, "QGW", "QGW");
  ncm_model_class_add_params (model_class, NC_HICOSMO_QGW_SPARAM_LEN, 0, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGW_H0, "H_0", "H0",
                              10.0, 500.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGW_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGW_OMEGA_W, "\\Omega_{w0}", "Omegaw",
                              1e-8,  10.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGW_DEFAULT_OMEGA_W,
                              NCM_PARAM_TYPE_FREE);
  /* Set w param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGW_W, "w", "w",
                              1e-50,  1.0, 1.0e-8,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGW_DEFAULT_W,
                              NCM_PARAM_TYPE_FIXED);
  /* Set xb param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGW_X_B, "x_b", "xb",
                              1.0e10,  1.0e40, 1.0e25,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGW_DEFAULT_X_B,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_qgw_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_qgw_E2);
  nc_hicosmo_set_Omega_c0_impl  (parent_class, &_nc_hicosmo_qgw_Omega_c0);
  nc_hicosmo_set_Omega_t0_impl  (parent_class, &_nc_hicosmo_qgw_Omega_t0);
  nc_hicosmo_set_xb_impl        (parent_class, &_nc_hicosmo_qgw_xb);
  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_qgw_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_qgw_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl   (parent_class, &_nc_hicosmo_qgw_bgp_cs2);
}

#define VECTOR   (NCM_MODEL (cosmo))
#define MACRO_H0 (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_QGW_H0))
#define OMEGA_W  (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_QGW_OMEGA_W))
#define W        (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_QGW_W))
#define X_B      (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_QGW_X_B))

static gdouble _nc_hicosmo_qgw_adiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_qgw_adiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_qgw_adiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_qgw_adiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_qgw_adiab_eval_unit (NcHIPertIAdiab *iad);
static gdouble _nc_hicosmo_qgw_adiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau);
static gdouble _nc_hicosmo_qgw_adiab_eval_p2Psi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_qgw_adiab_eval_p2drho (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);

static void
nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_xi     = &_nc_hicosmo_qgw_adiab_eval_xi;
  iface->eval_F1     = &_nc_hicosmo_qgw_adiab_eval_F1;
  iface->eval_nu     = &_nc_hicosmo_qgw_adiab_eval_nu;
  iface->eval_m      = &_nc_hicosmo_qgw_adiab_eval_m;
  iface->eval_unit   = &_nc_hicosmo_qgw_adiab_eval_unit;
  iface->eval_x      = &_nc_hicosmo_qgw_adiab_eval_x;
  iface->eval_p2Psi  = &_nc_hicosmo_qgw_adiab_eval_p2Psi;
  iface->eval_p2drho = &_nc_hicosmo_qgw_adiab_eval_p2drho;
}

/*
 * Interface implementation -- NcHIPertIAdiab
 */

static gdouble
_nc_hicosmo_qgw_adiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoQGW *qgw       = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo        = NC_HICOSMO (iad);
  const gdouble w         = W;
  const gdouble xb        = X_B;
  const gdouble Omega_w   = OMEGA_W;
  const gdouble xb2       = xb * xb;
  const gdouble cs        = sqrt (w);
  const gdouble tanh_tau  = tanh (tau);
  const gdouble tanh_tau2 = tanh_tau * tanh_tau;
  const gdouble lnfact    = log (3.0 * (1.0 + w) * k / (cs * tanh_tau2));
  const gdouble x         = _nc_hicosmo_qgw_adiab_eval_x (iad, tau);

  return lnfact - 2.0 * log (x);
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoQGW *qgw          = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo           = NC_HICOSMO (iad);
  const gdouble w            = W;
  const gdouble xb           = X_B;
  const gdouble Omega_w      = OMEGA_W;
  const gdouble xb_1p3w_2    = pow (xb, (1.0 + 3.0 * w) / 2.0);
  const gdouble sqrt_Omega_w = sqrt (Omega_w);
  const gdouble cs           = sqrt (w);
  const gdouble tanh_tau     = tanh (tau);
  const gdouble cosh_tau     = cosh (tau);
  const gdouble cosh_2tau    = cosh (2.0 * tau);
  const gdouble fact         = xb_1p3w_2 * sqrt_Omega_w / (2.0 * cs * k);

  return fact / pow (cosh_tau, (7.0 - 3.0 * w) / (3.0 * (1.0 - w))) / tanh_tau * (cosh_2tau + 3.0 * w - 4.0);
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoQGW *qgw          = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo           = NC_HICOSMO (iad);
  const gdouble w            = W;
  const gdouble Omega_w      = OMEGA_W;
  const gdouble w2           = w * w;
  const gdouble sqrt_Omega_w = sqrt (Omega_w);
  const gdouble cs           = sqrt (w);
  const gdouble x            = _nc_hicosmo_qgw_adiab_eval_x (iad, tau);

  return cs * 2.0 / (3.0 * (1.0 - w)) * k / sqrt_Omega_w / pow (x, 0.5 + 1.5 * w);
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoQGW *qgw          = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo           = NC_HICOSMO (iad);
  const gdouble w            = W;
  const gdouble Omega_w      = OMEGA_W;
  const gdouble w2           = w * w;
  const gdouble sqrt_Omega_w = sqrt (Omega_w);
  const gdouble cs2          = w;
  const gdouble tanh_tau     = tanh (tau);
  const gdouble tanh_tau2    = tanh_tau * tanh_tau;
  const gdouble x            = _nc_hicosmo_qgw_adiab_eval_x (iad, tau);

  return 9.0 * (1.0 - w2) * sqrt_Omega_w / (2.0 * pow (x, 1.5 * (1.0 - w)) * cs2 * tanh_tau2);
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_unit (NcHIPertIAdiab *iad)
{
  NcHICosmoQGW *qgw    = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo     = NC_HICOSMO (iad);
  const gdouble RH_lp  = nc_hicosmo_RH_planck (cosmo);
  const gdouble factor = sqrt (8.0 * ncm_c_pi () / 3.0);

  return factor / RH_lp;
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau)
{
  NcHICosmoQGW *qgw = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo  = NC_HICOSMO (iad);
  const gdouble xb  = X_B;
  const gdouble w   = W;

  return xb / pow (cosh (tau), 2.0 / (3.0 * (1 - w)));
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_p2Psi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoQGW *qgw          = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo           = NC_HICOSMO (iad);
  const gdouble w            = W;
  const gdouble Omega_w      = OMEGA_W;
  const gdouble tanh_tau     = tanh (tau);
  const gdouble x            = _nc_hicosmo_qgw_adiab_eval_x (iad, tau);
  const gdouble sqrt_Omega_w = sqrt (Omega_w);

  return -sqrt_Omega_w *pow (x, 2.5 + 1.5 *w) * tanh_tau / (2.0 * k * k);
}

static gdouble
_nc_hicosmo_qgw_adiab_eval_p2drho (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoQGW *qgw          = NC_HICOSMO_QGW (iad);
  NcHICosmo *cosmo           = NC_HICOSMO (iad);
  const gdouble w            = W;
  const gdouble Omega_w      = OMEGA_W;
  const gdouble tanh_tau     = tanh (tau);
  const gdouble x            = _nc_hicosmo_qgw_adiab_eval_x (iad, tau);
  const gdouble sqrt_Omega_w = sqrt (Omega_w);

  return -tanh_tau *pow (x, 1.5 *(1.0 - w)) / (sqrt_Omega_w * (1.0 + w));
}

/*
 * NcHICosmo virtual methods
 */

static gdouble
_nc_hicosmo_qgw_E2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQGW *qgw = NC_HICOSMO_QGW (cosmo);

  NcHICosmoQGWPrivate * const self = qgw->priv = nc_hicosmo_qgw_get_instance_private (qgw);
  const gdouble x                  = 1.0 + z;
  const gdouble x3_1pw             = pow (x, 3.0 * (1.0 + W));
  const gdouble Omega_w            = OMEGA_W;

  return Omega_w * x3_1pw;
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_qgw_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x         = 1.0 + z;
  const gdouble three_1pw = 3.0 * (1.0 + W);
  const gdouble x2p3w     = pow (x, three_1pw - 1.0);

  return three_1pw * x2p3w;
}

static gdouble
_nc_hicosmo_qgw_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x         = 1.0 + z;
  const gdouble three_1pw = 3.0 * (1.0 + W);
  const gdouble twop3w    = 2.0 + 3.0 * W;

  return three_1pw * twop3w * pow (x, twop3w - 1.0);
}

/****************************************************************************
 * Speed of sound squared
 ****************************************************************************/
static gdouble
_nc_hicosmo_qgw_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble w = W;

  return w;
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_qgw_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0 * OMEGA_W;
}

static gdouble
_nc_hicosmo_qgw_Omega_t0 (NcHICosmo *cosmo)
{
  return 1.0;
}

static gdouble
_nc_hicosmo_qgw_Omega_c0 (NcHICosmo *cosmo)
{
  return OMEGA_W;
}

static gdouble
_nc_hicosmo_qgw_xb (NcHICosmo *cosmo)
{
  return X_B;
}

/**
 * nc_hicosmo_qgw_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoQGW *
nc_hicosmo_qgw_new (void)
{
  NcHICosmoQGW *qgw = g_object_new (NC_TYPE_HICOSMO_QGW, NULL);

  return qgw;
}

