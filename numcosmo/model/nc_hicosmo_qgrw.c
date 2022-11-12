/***************************************************************************
 *            nc_hicosmo_qgrw.c
 *
 *  Wed June 04 10:04:24 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hicosmo_qgrw.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_hicosmo_qgrw
 * @title: NcHICosmoQGRW
 * @short_description: Radiation plus $w$-fluid model with a quantum generated bounce phase model.
 *
 * In this model the adiabatic mode $\zeta$ has its mass, speed of sound square $c_s^2$ and frequency square $\nu_\zeta^2$ given by
 * \begin{align}
 * m_\zeta &= 3 \Delta_\bar{K}\sqrt{\Omega_{w0}} x^{-3(1-w)/2}\frac{(1 + w) +  4R/3}{c_s^2}\frac{1}{\sqrt{(1-exp(-2\vert\alpha\vert)) + (1-exp(-3(1-w)\vert\alpha\vert))}}, \\\\
 * c_s^2 &= \frac{w (1 + w) + 4R/9}{(1+w) + 4R/3}, \\\\
 * \nu_\zeta^2 &= \frac {c_s^2 k^2}{\Omega_{w0} x^{1+3w} ((1-exp(-2\vert\alpha\vert)) + (1-exp(-3(1-w)\vert\alpha\vert)))},
 * \end{align}
 * where $$R \equiv \frac{\Omega_{r0} x}{\Omega_{w0} x^{3w}}.$$
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qgrw.h"

static void nc_hipert_itwo_fluids_interface_init (NcHIPertITwoFluidsInterface *iface);

G_DEFINE_TYPE_WITH_CODE (NcHICosmoQGRW, nc_hicosmo_qgrw, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_ITWO_FLUIDS,
                                                nc_hipert_itwo_fluids_interface_init)
                         );

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_qgrw_init (NcHICosmoQGRW *qgrw)
{

}

static void
nc_hicosmo_qgrw_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qgrw_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_qgrw_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgrw_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgrw_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgrw_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qgrw_bgp_cs2 (NcHICosmo *cosmo, gdouble z);

static gdouble _nc_hicosmo_qgrw_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qgrw_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qgrw_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qgrw_xb (NcHICosmo *cosmo);

static NcHIPertITwoFluidsEOM *_nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
static NcHIPertITwoFluidsTV *_nc_hipert_itwo_fluids_tv (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

static void
nc_hicosmo_qgrw_class_init (NcHICosmoQGRWClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  ncm_model_class_set_name_nick (model_class, "QGRW", "QGRW");
  ncm_model_class_add_params (model_class, NC_HICOSMO_QGRW_SPARAM_LEN, 0, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_r0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_OMEGA_R, "\\Omega_{r0}", "Omegar",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_OMEGA_R,
                               NCM_PARAM_TYPE_FREE);
  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_OMEGA_W, "\\Omega_{w0}", "Omegaw",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_OMEGA_W,
                               NCM_PARAM_TYPE_FREE);
  /* Set w param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_W, "w", "w",
                               1e-50,  1.0, 1.0e-8,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_W,
                               NCM_PARAM_TYPE_FIXED);
  /* Set xb param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_X_B, "x_b", "xb",
                               1.0e10,  1.0e40, 1.0e25,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_X_B,
                               NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_qgrw_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_qgrw_E2);
  nc_hicosmo_set_Omega_c0_impl  (parent_class, &_nc_hicosmo_qgrw_Omega_c0);
  nc_hicosmo_set_Omega_t0_impl  (parent_class, &_nc_hicosmo_qgrw_Omega_t0);
  nc_hicosmo_set_xb_impl        (parent_class, &_nc_hicosmo_qgrw_xb);
  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_qgrw_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_qgrw_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl   (parent_class, &_nc_hicosmo_qgrw_bgp_cs2);

  object_class->finalize = nc_hicosmo_qgrw_finalize;
}

static void
nc_hipert_itwo_fluids_interface_init (NcHIPertITwoFluidsInterface *iface)
{
  iface->eom = &_nc_hipert_itwo_fluids_eom;
  iface->tv  = &_nc_hipert_itwo_fluids_tv;
}

#define VECTOR   (NCM_MODEL (cosmo)->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_H0))
#define OMEGA_R  (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_OMEGA_R))
#define OMEGA_W  (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_OMEGA_W))
#define W        (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_W))
#define X_B      (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_X_B))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_qgrw_E2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x   = 1.0 + z;
  const gdouble x2  = x * x;
  const gdouble x3  = x2 * x;
  const gdouble x4  = x2 * x2;
  const gdouble x3w = pow (x3, W);
  const gdouble x6  = x3 * x3;
  const gdouble xb  = X_B;
  const gdouble xb2 = xb * xb;
  const gdouble xb3 = xb2 * xb;
  const gdouble xb3w = pow (xb3, W);

#ifdef NC_HICOSMO_QGRW_CHECK_INTERVAL
  if (G_UNLIKELY (ncm_cmp (x, xb, 1e-4) == 0))
    g_warning ("_nc_hicosmo_qgrw_E2: used outside if its valid interval.");
#endif /* NC_HICOSMO_QGRW_CHECK_INTERVAL */

  return (OMEGA_R * (x4 - x6 / xb2) + OMEGA_W * (x3 * x3w - x6 * xb3w / xb3));
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_qgrw_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x5 = x3 * x2;
  const gdouble x3w = pow (x3, W);
  const gdouble xb  = X_B;
  const gdouble xb2 = xb * xb;
  const gdouble xb3 = xb2 * xb;
  const gdouble xb3w = pow (xb3, W);
  const gdouble poly = OMEGA_R * (4.0 * x3 - 6.0 * x5 / xb2) + OMEGA_W * (3.0 * (1.0 + W) * x2 * x3w - 6.0 * x5 * xb3w / xb3);

#ifdef NC_HICOSMO_QGRW_CHECK_INTERVAL
  if (G_UNLIKELY (ncm_cmp (x, xb, 1e-4) == 0))
    g_warning ("_nc_hicosmo_qgrw_E2: used outside if its valid interval.");
#endif /* NC_HICOSMO_QGRW_CHECK_INTERVAL */

  return poly;
}

static gdouble
_nc_hicosmo_qgrw_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x2 * x2;
  const gdouble x3w = pow (x3, W);
  const gdouble three1p3w = 3.0 * (1.0 + W);
  const gdouble three1p3wm1 = three1p3w - 1.0;
  const gdouble xb  = X_B;
  const gdouble xb2 = xb * xb;
  const gdouble xb3 = xb2 * xb;
  const gdouble xb3w = pow (xb3, W);

  const gdouble poly = OMEGA_R * (12.0 * x2 - 30.0 * x4 / xb2) + OMEGA_W * (three1p3w * three1p3wm1 * x * x3w - 30.0 * x4  * xb3w / xb3);

#ifdef NC_HICOSMO_QGRW_CHECK_INTERVAL
  if (G_UNLIKELY (ncm_cmp (x, xb, 1e-4) == 0))
    g_warning ("_nc_hicosmo_qgrw_E2: used outside if its valid interval.");
#endif /* NC_HICOSMO_QGRW_CHECK_INTERVAL */

  return poly;
}

/****************************************************************************
 * Speed of sound squared
 ****************************************************************************/
static gdouble
_nc_hicosmo_qgrw_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble w  = W;
  const gdouble x3w = pow (x3, w);
  const gdouble R = (4.0 * OMEGA_R * x / (3.0 * (1.0 + w) * OMEGA_W * x3w));

  return (w + R * 1.0 / 3.0) / (1.0 + R);
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_qgrw_H0 (NcHICosmo *cosmo) { return MACRO_H0 * (OMEGA_R + OMEGA_W); }
static gdouble _nc_hicosmo_qgrw_Omega_t0 (NcHICosmo *cosmo) { return 1.0; }
static gdouble _nc_hicosmo_qgrw_Omega_c0 (NcHICosmo *cosmo) { return OMEGA_W; }
static gdouble _nc_hicosmo_qgrw_xb (NcHICosmo *cosmo) { return X_B; }

/****************************************************************************
 * Wkb theta
 ****************************************************************************/

#define _d1ln(pn) (d1##pn / pn)
#define _d2ln(pn) (d2##pn / pn - d1ln##pn * d1ln##pn)
#define _d3ln(pn) (d3##pn / pn - d1ln##pn * d1ln##pn * d1ln##pn - 3.0 * d1ln##pn * d2ln##pn)

#define _NC_HICOSMO_QGRW_WKB_COMMON1 \
  const gdouble w2       = 1.0 / 3.0; \
  const gdouble alpha2   = alpha * alpha; \
  const gdouble alpha2_2 = alpha2 * 0.5; \
  const gdouble w        = W; \
  const gdouble Omegar   = OMEGA_R; \
  const gdouble Omegaw   = OMEGA_W; \
 \
  const gdouble x   = nc_hicosmo_x_alpha (cosmo, alpha); \
  const gdouble x2  = x * x; \
  const gdouble x3  = x2 * x; \
  const gdouble x3w = pow (x3, w); \
 \
  const gdouble onepw       = 1.0 + w; \
  const gdouble threex1mw   = 3.0 * (1.0 - w); \
  const gdouble onepw2      = 4.0 / 3.0; \
 \
  const gdouble R    = Omegar * x / (Omegaw * x3w); \
  const gdouble p2f1 = 0.5 * threex1mw * ncm_exprel (-threex1mw * alpha2_2); \
  const gdouble p2f2 = R * ncm_exprel (-2.0 * alpha2_2); \
  const gdouble p2   = p2f1 + p2f2; \
  const gdouble p4   = onepw + onepw2 * R;

#define _NC_HICOSMO_QGRW_WKB_COMMON2 \
  const gdouble onem3w      = 1.0 - 3.0 * w; \
 \
  const gdouble d1R  = onem3w * R; \
 \
  const gdouble d1p2f1 = 0.5 * gsl_pow_2 (threex1mw) * ncm_d1exprel (-threex1mw * alpha2_2); \
  const gdouble d1p2f2 = 2.0 * R * ncm_d1exprel (-2.0 * alpha2_2) + d1R * ncm_exprel (-2.0 * alpha2_2); \
  const gdouble d1p2   = d1p2f1 + d1p2f2; \
  const gdouble d1p4   = onepw2 * d1R; \
 \
  const gdouble d1lnp2 = _d1ln (p2); \
  const gdouble d1lnp4 = _d1ln (p4);

#define _NC_HICOSMO_QGRW_WKB_COMMON22 \
  const gdouble d2R    = onem3w * d1R;

#define _NC_HICOSMO_QGRW_WKB_COMMON222 \
  const gdouble d2p2f1 = 0.5 * gsl_pow_3 (threex1mw) * ncm_d2exprel (-threex1mw * alpha2_2); \
  const gdouble d2p2f2 = 4.0 * R * ncm_d2exprel (-2.0 * alpha2_2) + 4.0 * d1R * ncm_d1exprel (-2.0 * alpha2_2) + d2R * ncm_exprel (-2.0 * alpha2_2); \
  const gdouble d2p2   = d2p2f1 + d2p2f2; \
  const gdouble d2lnp2 = _d2ln (p2);

#define _NC_HICOSMO_QGRW_WKB_COMMON2222 \
  const gdouble d2p4   = onepw2 * d2R; \
  const gdouble d2lnp4 = _d2ln (p4);

#define _NC_HICOSMO_QGRW_WKB_COMMON3 \
  const gdouble p1   = w * onepw + w2 * onepw2 * R;

#define _NC_HICOSMO_QGRW_WKB_COMMON33 \
  const gdouble p5   = sqrt (x3 / (x3w * Omegaw)) / 3.0;

#define _NC_HICOSMO_QGRW_WKB_COMMON4 \
  const gdouble k2          = k * k; \
  const gdouble p0   = k2 / (Omegaw * x * x3w);

#define _NC_HICOSMO_QGRW_WKB_COMMON44 \
  const gdouble p3       = w2 * onepw + w * onepw2 * R;

#define _NC_HICOSMO_QGRW_WKB_COMMON5 \
  const gdouble onep3w   = 1.0 + 3.0 * w; \
  const gdouble d1p0     = - onep3w * p0; \
  const gdouble d1lnp0   = _d1ln (p0);

#define _NC_HICOSMO_QGRW_WKB_COMMON51 \
  const gdouble d1p1     = w2 * onepw2 * d1R; \
  const gdouble d1lnp1   = _d1ln (p1);

#define _NC_HICOSMO_QGRW_WKB_COMMON53 \
  const gdouble d1p3     = w * onepw2 * d1R; \
  const gdouble d1lnp3   = _d1ln (p3);

#define _NC_HICOSMO_QGRW_WKB_D1LNSQRTNUS \
  const gdouble d1lnsqrtnuS = 0.25 * (d1lnp0 + d1lnp3 - d1lnp2 - d1lnp4);

#define _NC_HICOSMO_QGRW_WKB_D1LNSQRTNUZETA \
  const gdouble d1lnsqrtnuzeta = 0.25 * (d1lnp0 + d1lnp1 - d1lnp2 - d1lnp4);

#define _NC_HICOSMO_QGRW_WKB_COMMON56 \
  const gdouble onep3w_2 = onep3w * 0.5; \
  const gdouble d1lnp6   = onep3w_2;

#define _NC_HICOSMO_QGRW_WKB_COMMON55 \
  const gdouble d2p3     = w * onepw2 * d2R;

#define _NC_HICOSMO_QGRW_WKB_COMMON555 \
  const gdouble d2lnp0   = 0.0;

#define _NC_HICOSMO_QGRW_WKB_COMMON5555 \
  const gdouble d2lnp3   = _d2ln (p3); \
  const gdouble d2lnp6   = 0.0;

#define _NC_HICOSMO_QGRW_WKB_NUZETA2 \
  const gdouble nuzeta2 = p0 * p1 / (p2 * p4);

#define _NC_HICOSMO_QGRW_WKB_D1LNP5 \
  const gdouble threex1mw_2 = threex1mw * 0.5; \
  const gdouble d1p5   = threex1mw_2 * p5; \
  const gdouble d1lnp5 = _d1ln (p5);

#define _NC_HICOSMO_QGRW_WKB_DLNSQRTMZETANUZETA \
  _NC_HICOSMO_QGRW_WKB_D1LNP5; \
  const gdouble d1lnsqrtmzetanuzeta = 0.75 * d1lnp4 + 0.25 * d1lnp0 - 0.25 * d1lnp1 - 0.5 * d1lnp2 - 0.5 * d1lnp5 + 1.0 / alpha2;

#define _NC_HICOSMO_QGRW_WKB_VZETA \
  const gdouble d2p1   = w2 * onepw2 * d2R; \
  const gdouble d2lnp1 = _d2ln (p1); \
  const gdouble d2lnp5 = 0.0; \
  const gdouble d2lnsqrtmzetanuzeta = 0.75 * d2lnp4 + 0.25 * d2lnp0 - 0.25 * d2lnp1 - 0.5 * d2lnp2 - 0.5 * d2lnp5 + 2.0 / (alpha2 * alpha2); \
  const gdouble d2sqrtmzetanuzeta_sqrtmzetanuzeta = d1lnsqrtmzetanuzeta * d1lnsqrtmzetanuzeta + d2lnsqrtmzetanuzeta; \
  const gdouble Vzeta   = (alpha2 * d2sqrtmzetanuzeta_sqrtmzetanuzeta -  d1lnsqrtmzetanuzeta) - 2.0 * alpha2 * d1lnsqrtmzetanuzeta * d1lnsqrtnuzeta;

#define _NC_HICOSMO_QGRW_WKB_NUA2 \
  const gdouble nuA2    = nuzeta2 - Vzeta;

#define _NC_HICOSMO_QGRW_WKB_MZETA \
  const gdouble mzeta   = p4 * p4 / (alpha2 * p1 * p5 * sqrt (p2));

#define _NC_HICOSMO_QGRW_WKB_NUS2 \
  const gdouble nuS2 = p0 * p3 / (p2 * p4);

#define _NC_HICOSMO_QGRW_WKB_DLNSQRTMSNUS \
  const gdouble d1lnsqrtmSnuS = (0.25 * d1lnp0 + 0.5 * d1lnp6) * 0.0 + 0.75 * d1lnp4 - 0.25 * d1lnp3;

#define _NC_HICOSMO_QGRW_WKB_VS \
  const gdouble d2lnsqrtmSnuS = 0.25 * d2lnp0 + 0.75 * d2lnp4 + 0.5 * d2lnp6 - 0.25 * d2lnp3; \
  const gdouble d2sqrtmSnuS_sqrtmSnuS = d1lnsqrtmSnuS * d1lnsqrtmSnuS + d2lnsqrtmSnuS; \
  const gdouble VS = (alpha2 * d2sqrtmSnuS_sqrtmSnuS -  d1lnsqrtmSnuS) - 2.0 * alpha2 * d1lnsqrtmSnuS * d1lnsqrtnuS;

#define _NC_HICOSMO_QGRW_WKB_NUB2 \
  const gdouble nuB2    = nuS2 - VS;

#define _NC_HICOSMO_QGRW_WKB_P6 \
  const gdouble p6      = sqrt (x3 / (x3w * Omegaw)) / (3.0 * onepw * onepw2 * R);

#define _NC_HICOSMO_QGRW_WKB_MS \
  _NC_HICOSMO_QGRW_WKB_P6; \
  const gdouble mS      = sqrt (p2) * p4 * p4 * p6 / p3;

#define _NC_HICOSMO_QGRW_ONE_M_Y2 \
  const gdouble one_m_y2 = w * w2 * gsl_pow_2 (onepw + onepw2 * R) / (p1 * p3);

#define _NC_HICOSMO_QGRW_ODE_FULL_DIFF \
  /*const gdouble d1lnp1mp3 = -onepw2 * d1R / (onepw - onepw2 * R);*/ \
  const gdouble dlnmzetanuzeta2 = -2.0 / alpha - alpha * (d1lnp0 - 1.5 * d1lnp2 - d1lnp5); \
  const gdouble dlnmSnuS2 = - alpha * (d1lnp0 - 0.5 * d1lnp2 + d1lnp4 + d1lnp6); \
  const gdouble dlnmzetanuzeta2_mSnuS2 = -2.0 / alpha + alpha * (d1lnp2 + d1lnp4 + d1lnp5 + d1lnp6); \
  const gdouble dlnlambda_s2_lambda_zeta2 = - alpha * d1R / R; \
  /*const gdouble dlnnu_plus2 = - alpha * (d1lnp0 - d1lnp2);*/ \
  const gdouble dlnnu_minus2 = -alpha * (d1lnp0 - d1lnp2); \
  const gdouble dlnmSnuS2_dlnnu_plus = - alpha * (d1lnp4 + 0.0 * (0.5 * d1lnp0 + d1lnp6));

/****************************************************************************
 * Equations of motion
 ****************************************************************************/

static NcHIPertITwoFluidsEOM *
_nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (itf);
  NcHICosmoQGRW *qgrw = NC_HICOSMO_QGRW (itf);
  if (qgrw->eom_two_fluids.skey  != NCM_MODEL (itf)->pkey ||
      qgrw->eom_two_fluids.alpha != alpha ||
      qgrw->eom_two_fluids.k     != k)
  {
    const gdouble epsilon   = GSL_SIGN (alpha);
    const gdouble absalpha  = fabs (alpha);
    const gdouble w2        = W;
    const gdouble x         = X_B * exp (- absalpha);
    const gdouble x2        = x * x;
    const gdouble x3        = x2 * x;
    const gdouble x4        = x2 * x2;
    const gdouble x3_1pw    = x3 * pow (x3, w2);
    const gdouble x_m3w     = x / pow (x3, w2);
    const gdouble three_1pw = 3.0 * (1.0 + w2);
    const gdouble three_1mw = 3.0 * (1.0 - w2);
    const gdouble x_xb3_1mw = exp (- three_1mw * absalpha);
    const gdouble x_xb2     = exp (- 2.0 * absalpha);
    const gdouble oneF1_r   = ncm_exprel (-2.0 * absalpha);
    const gdouble oneF1_p   = ncm_exprel (-three_1mw * absalpha);
    const gdouble R0        = OMEGA_R / OMEGA_W;

    const gdouble F         = (R0 * x_m3w * 2.0 * oneF1_r + three_1mw * oneF1_p) * epsilon * alpha;
    const gdouble d1F       = epsilon * (three_1mw * x_xb3_1mw + x_m3w * R0 * (three_1mw * (x_xb2 - 1.0) + 2.0));
    const gdouble d2F       = - x_xb3_1mw * three_1mw * three_1mw + R0 * x_m3w * (gsl_pow_2 (1.0 - 3.0 * w2) - x_xb2 * three_1mw * three_1mw);
    
    const gdouble E2        = OMEGA_W * x3_1pw * F;

    const gdouble d1lnF     = d1F / F;
    
    const gdouble d1lnE2    = (- three_1pw * epsilon + d1lnF);
    const gdouble d2lnE2    = (d2F / F - d1lnF * d1lnF);

    const gdouble d1lnE     = 0.5 * d1lnE2;
    const gdouble d2lnE     = 0.5 * d2lnE2;
    
    const gdouble d1lnE_2   = d1lnE * d1lnE;

    const gdouble absE      = sqrt (E2);
    const gdouble c12       = 1.0 / 3.0;
    const gdouble c22       = w2;
    const gdouble c14       = c12 * c12;
    const gdouble c24       = c22 * c22;
    const gdouble c1        = sqrt (c12);
    const gdouble c2        = sqrt (c22);
    const gdouble sc1       = sqrt (c1);
    const gdouble sc2       = sqrt (c2);
    const gdouble Fnu       = x * k / absE;
    const gdouble rhopp1    = OMEGA_R * (1.0 + 1.0 / 3.0) * x4;
    const gdouble rhopp2    = OMEGA_W * (1.0 + w2) * x3_1pw;
    const gdouble phip      = - atan (epsilon * sqrt (rhopp1 / rhopp2));
    const gdouble sin_phi   = - cos (phip);
    const gdouble cos_phi   = sin (phip);
    const gdouble sin_2phi  = 2.0 * sin_phi * cos_phi;
    const gdouble sin2_phi  = sin_phi * sin_phi;
    const gdouble cos2_phi  = cos_phi * cos_phi;
    const gdouble cos_2phi  = - cos (2.0 * phip);
    const gdouble sin2_2phi = sin_2phi * sin_2phi;

    qgrw->eom_two_fluids.nu1 = c1 * Fnu;
    qgrw->eom_two_fluids.nu2 = c2 * Fnu;
    
    qgrw->eom_two_fluids.gammabar11 = 
      (
       + 0.0 * (9.0 * c14 - 1.0) / 4.0 /* c12 == 1/3 */
       - (epsilon * (6.0 * sin2_2phi * (c22 - c12) + 4.0 * (1.0 - 3.0 * c12) * 0.0) * d1lnE) / 8.0
       + (c1 - c2) * sin2_2phi * d1lnE_2 / (2.0 * c2)
       - cos2_phi * d2lnE
       ) / (c1 * Fnu);
    
    qgrw->eom_two_fluids.gammabar22 = 
      (
       + 1.0 * (9.0 * c24 - 1.0) / 4.0
       - (epsilon * (6.0 * sin2_2phi * (c12 - c22) + 4.0 * (1.0 - 3.0 * c22)) * d1lnE) / 8.0
       + (c2 - c1) * sin2_2phi * d1lnE_2 / (2.0 * c1)
       - sin2_phi * d2lnE
       ) / (c2 * Fnu);
    
    qgrw->eom_two_fluids.gammabar12 = 
      (
       + epsilon * sin_2phi * (3.0 * cos_2phi * c1 * c2 * (c12 - c22) - gsl_pow_2 (c1 - c2)) * d1lnE / (2.0 * c1 * c2)
       + sin_2phi * (c1 - c2) * (cos2_phi * c2 - sin2_phi * c1) * d1lnE_2 / (c1 * c2)
       + sin_2phi * d2lnE
       ) / (2.0 * Fnu * sc1 * sc2);
    
    qgrw->eom_two_fluids.taubar = sin_2phi * (c2 - c1) * d1lnE / (sc1 * sc2);

    {
      const gdouble rhopp = rhopp1 + rhopp2;
      const gdouble cs2   = c1 * c1 * cos2_phi + c2 * c2 * sin2_phi;
      const gdouble cm2   = c2 * c2 * cos2_phi + c1 * c1 * sin2_phi;

      qgrw->eom_two_fluids.m_zeta    = rhopp / (cs2 * absE * x3);
      qgrw->eom_two_fluids.m_s       = absE * x3 / (cm2 * rhopp * cos2_phi * sin2_phi);
      qgrw->eom_two_fluids.mnu2_zeta = rhopp * k * k / (x * gsl_pow_3 (absE));
      qgrw->eom_two_fluids.mnu2_s    = x3 * x2 * k * k / (absE * rhopp * cos2_phi * sin2_phi);
      qgrw->eom_two_fluids.y         = epsilon * (c1 * c1 - c2 * c2) * cos2_phi * sin2_phi;
    }
    
    qgrw->eom_two_fluids.skey  = NCM_MODEL (cosmo)->pkey;
    qgrw->eom_two_fluids.alpha = alpha;
    qgrw->eom_two_fluids.k     = k;

    /*printf ("% 20.15g % 20.15g % 20.15g % 20.15e\n", alpha, nuS2, VS, fabs ((nuS2 - VS)/ VS));*/
  }

  return &qgrw->eom_two_fluids;
}

static NcHIPertITwoFluidsTV *
_nc_hipert_itwo_fluids_tv (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  NcHICosmo *cosmo = NC_HICOSMO (itf);
  NcHICosmoQGRW *qgrw = NC_HICOSMO_QGRW (itf);
  if (qgrw->tv_two_fluids.skey  != NCM_MODEL (itf)->pkey ||
      qgrw->tv_two_fluids.alpha != alpha ||
      qgrw->tv_two_fluids.k     != k)
  {
    const gdouble epsilon   = GSL_SIGN (alpha);
    const gdouble absalpha  = fabs (alpha);
    const gdouble w2        = W;
    const gdouble x         = X_B * exp (- absalpha);
    const gdouble x2        = x * x;
    const gdouble x3        = x2 * x;
    const gdouble x4        = x2 * x2;
    const gdouble x3_1pw    = x3 * pow (x3, w2);
    const gdouble x_m3w     = x / pow (x3, w2);
    const gdouble three_1pw = 3.0 * (1.0 + w2);
    const gdouble three_1mw = 3.0 * (1.0 - w2);
    const gdouble x_xb3_1mw = exp (- three_1mw * absalpha);
    const gdouble x_xb2     = exp (- 2.0 * absalpha);
    const gdouble oneF1_r   = ncm_exprel (-2.0 * absalpha);
    const gdouble oneF1_p   = ncm_exprel (-three_1mw * absalpha);
    const gdouble R0        = OMEGA_R / OMEGA_W;

    const gdouble F         = (R0 * x_m3w * 2.0 * oneF1_r + three_1mw * oneF1_p) * epsilon * alpha;
    const gdouble d1F       = epsilon * (three_1mw * x_xb3_1mw + x_m3w * R0 * (three_1mw * (x_xb2 - 1.0) + 2.0));
    
    const gdouble E2        = OMEGA_W * x3_1pw * F;

    const gdouble d1lnF     = d1F / F;
    
    const gdouble d1lnE2    = (- three_1pw * epsilon + d1lnF);

    const gdouble d1lnE     = 0.5 * d1lnE2;
    
    const gdouble absE      = sqrt (E2);
    const gdouble c12       = 1.0 / 3.0;
    const gdouble c22       = w2;
    const gdouble c1        = sqrt (c12);
    const gdouble c2        = sqrt (c22);
    const gdouble sc1       = sqrt (c1);
    const gdouble sc2       = sqrt (c2);
    const gdouble Fnu       = x * k / absE;
    const gdouble rhopp1    = OMEGA_R * (1.0 + 1.0 / 3.0) * x4;
    const gdouble rhopp2    = OMEGA_W * (1.0 + w2) * x3_1pw;
    const gdouble phip      = - atan (epsilon * sqrt (rhopp1 / rhopp2));
    const gdouble sin_phi   = - cos (phip);
    const gdouble cos_phi   = sin (phip);
    const gdouble sin_2phi  = 2.0 * sin_phi * cos_phi;
    const gdouble sin2_phi  = sin_phi * sin_phi;
    const gdouble cos2_phi  = cos_phi * cos_phi;

    const gdouble rhopp     = rhopp1 + rhopp2;
    const gdouble srhopp    = sqrt (rhopp);
    const gdouble sk        = sqrt (k);
    const gdouble Az_sFnu   = srhopp * sk / (x * absE);
    const gdouble As_sFnu   = 2.0 * epsilon * x2 * sk / (sin_2phi * srhopp);

    qgrw->tv_two_fluids.zeta[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1]  = + cos_phi * sc1 / Az_sFnu;
    qgrw->tv_two_fluids.zeta[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2]  = - sin_phi * sc2 / Az_sFnu;
    qgrw->tv_two_fluids.zeta[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1]  = 0.0;
    qgrw->tv_two_fluids.zeta[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2]  = 0.0;

    qgrw->tv_two_fluids.s[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1]     = + sin_phi * sc1 / As_sFnu;
    qgrw->tv_two_fluids.s[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2]     = + cos_phi * sc2 / As_sFnu;
    qgrw->tv_two_fluids.s[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1]     = 0.0;
    qgrw->tv_two_fluids.s[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2]     = 0.0;
    
    qgrw->tv_two_fluids.Pzeta[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1] = + (cos_phi * Az_sFnu / (Fnu * sc1)) * (epsilon * (1.0 + 3.0 * c12) / (2.0 * c1) + (cos2_phi / c1 + sin2_phi / c2) * d1lnE);
    qgrw->tv_two_fluids.Pzeta[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2] = - (sin_phi * Az_sFnu / (Fnu * sc2)) * (epsilon * (1.0 + 3.0 * c22) / (2.0 * c2) + (cos2_phi / c1 + sin2_phi / c2) * d1lnE);
    qgrw->tv_two_fluids.Pzeta[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1] = + cos_phi * Az_sFnu / sc1;
    qgrw->tv_two_fluids.Pzeta[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2] = - sin_phi * Az_sFnu / sc2;

    /* 1.0 - 3.0 * c12 == 0.0 */
    qgrw->tv_two_fluids.Ps[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1]    = + (sin_phi * As_sFnu / (Fnu * sc1)) * (epsilon * 0.0 * (1.0 - 3.0 * c12) / (2.0 * c1) + cos2_phi * (1.0 / c1 - 1.0 / c2) * d1lnE);
    qgrw->tv_two_fluids.Ps[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2]    = + (cos_phi * As_sFnu / (Fnu * sc2)) * (epsilon * 1.0 * (1.0 - 3.0 * c22) / (2.0 * c2) + sin2_phi * (1.0 / c2 - 1.0 / c1) * d1lnE);
    qgrw->tv_two_fluids.Ps[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1]    = + sin_phi * As_sFnu / sc1;
    qgrw->tv_two_fluids.Ps[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2]    = + cos_phi * As_sFnu / sc2;
    
    qgrw->tv_two_fluids.skey  = NCM_MODEL (cosmo)->pkey;
    qgrw->tv_two_fluids.alpha = alpha;
    qgrw->tv_two_fluids.k     = k;
  }

  return &qgrw->tv_two_fluids;
}

/**
 * nc_hicosmo_qgrw_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoQGRW *
nc_hicosmo_qgrw_new (void)
{
  NcHICosmoQGRW *qgrw = g_object_new (NC_TYPE_HICOSMO_QGRW, NULL);
  return qgrw;
}

typedef NcHIPertWKB NcHIPertWKBQgrwZeta;
typedef NcHIPertWKBClass NcHIPertWKBQgrwZetaClass;

G_DEFINE_TYPE (NcHIPertWKBQgrwZeta, nc_hipert_wkb_qgrw_zeta, NC_TYPE_HIPERT_WKB);

static void 
nc_hipert_wkb_qgrw_zeta_init (NcHIPertWKBQgrwZeta *wkb_zeta)
{

}

static void _nc_hipert_wkb_qgrw_zeta_get_nu_V (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *nu, gdouble *V);
static void _nc_hipert_wkb_qgrw_zeta_get_mnu_dmnu (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *mnu, gdouble *dmnu);
static gdouble _nc_hipert_wkb_qgrw_zeta_get_m (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);
static gdouble _nc_hipert_wkb_qgrw_zeta_get_nu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);
static gdouble _nc_hipert_wkb_qgrw_zeta_get_dVnu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);

static void
nc_hipert_wkb_qgrw_zeta_class_init (NcHIPertWKBQgrwZetaClass *klass)
{
  klass->get_nu_V     = &_nc_hipert_wkb_qgrw_zeta_get_nu_V;
  klass->get_mnu_dmnu = &_nc_hipert_wkb_qgrw_zeta_get_mnu_dmnu;
  klass->get_m        = &_nc_hipert_wkb_qgrw_zeta_get_m;
  klass->get_nu2      = &_nc_hipert_wkb_qgrw_zeta_get_nu2;
  klass->get_dVnu2    = &_nc_hipert_wkb_qgrw_zeta_get_dVnu2;
}

static void 
_nc_hipert_wkb_qgrw_zeta_get_nu_V (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *nu, gdouble *V)
{

  
}

static void 
_nc_hipert_wkb_qgrw_zeta_get_mnu_dmnu (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *mnu, gdouble *dmnu)
{
}

static gdouble 
_nc_hipert_wkb_qgrw_zeta_get_m (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  g_assert_not_reached ();

  return 0.0;
}

static gdouble 
_nc_hipert_wkb_qgrw_zeta_get_nu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  g_assert_not_reached ();

  return 0.0;
}

static gdouble 
_nc_hipert_wkb_qgrw_zeta_get_dVnu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  g_assert_not_reached ();

  return 0.0;
}
