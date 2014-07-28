/***************************************************************************
 *            nc_hicosmo_qgrw.c
 *
 *  Wed June 04 10:04:24 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hicosmo_qgrw.c
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
 * SECTION:nc_hicosmo_qgrw
 * @title: Quantum Gravity Radiation $w$ model
 * @short_description: Radiation plus $w$-fluid model with a quantum generated bounce phase.
 *
 * In this model the adiabatic mode $\zeta$ has its mass, speed of sound square $c_s^2$ and frequency square $\nu_\zeta^2$ given by
 * \begin{align}
 * m_\zeta &= 3 \Delta_\bar{K}\sqrt{\Omega_w} x^{-3(1-w)/2}\frac{(1 + w) +  4R/3}{c_s^2}\frac{1}{\sqrt{(1-exp(-2\vert\alpha\vert)) + (1-exp(-3(1-w)\vert\alpha\vert))}}, \\\\
 * c_s^2 &= \frac{w (1 + w) + 4R/9}{(1+w) + 4R/3}, \\\\
 * \nu_\zeta^2 &= \frac {c_s^2 k^2}{\Omega_w x^{1+3w} ((1-exp(-2\vert\alpha\vert)) + (1-exp(-3(1-w)\vert\alpha\vert)))}, 
 * \end{align}
 * where $$R \equiv \frac{\Omega_r x}{\Omega_w x^{3w}}.$$
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qgrw.h"

#include <gsl/gsl_sf_hyperg.h>

static void nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface);
static void nc_hipert_itwo_fluids_interface_init (NcHIPertITwoFluidsInterface *iface);

G_DEFINE_TYPE_WITH_CODE (NcHICosmoQGRW, nc_hicosmo_qgrw, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IADIAB,
                                                nc_hipert_iadiab_interface_init)
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
  memset (&qgrw->eom_adiab_zeta, 0, sizeof (NcHIPertIAdiabEOM));
}

static void
nc_hicosmo_qgrw_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qgrw_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_qgrw_E2 (NcmModel *model, gdouble z);
static gdouble _nc_hicosmo_qgrw_dE2_dz (NcmModel *model, gdouble z);
static gdouble _nc_hicosmo_qgrw_E2 (NcmModel *model, gdouble z);
static gdouble _nc_hicosmo_qgrw_d2E2_dz2 (NcmModel *model, gdouble z);
static gdouble _nc_hicosmo_qgrw_cs2 (NcmModel *model, gdouble z);
static gdouble _nc_hicosmo_qgrw_rhopp (NcmModel *model, gdouble z);

static gdouble _nc_hicosmo_qgrw_H0 (NcmModel *model);
static gdouble _nc_hicosmo_qgrw_Omega_t (NcmModel *model);
static gdouble _nc_hicosmo_qgrw_Omega_c (NcmModel *model);
static gdouble _nc_hicosmo_qgrw_xb (NcmModel *model);

static gdouble _nc_hipert_iadiab_nuA2 (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
static gdouble _nc_hipert_iadiab_dlnmzeta (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
static gdouble _nc_hipert_iadiab_dmzetanuA_nuA (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
static NcHIPertIAdiabEOM *_nc_hipert_iadiab_eom (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);

static gdouble _nc_hipert_itwo_fluids_nuB2 (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
static gdouble _nc_hipert_itwo_fluids_dmSnuB_nuB (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
static NcHIPertITwoFluidsEOM *_nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

static void
nc_hicosmo_qgrw_class_init (NcHICosmoQGRWClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  ncm_model_class_add_params (model_class, NC_HICOSMO_QGRW_SPARAM_LEN, 0, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "QGRW", "QGRW");

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_r param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_OMEGA_R, "\\Omega_r", "Omegar",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_OMEGA_R,
                               NCM_PARAM_TYPE_FREE);
  /* Set Omega_x param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_OMEGA_W, "\\Omega_w", "Omegaw",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QGRW_DEFAULT_OMEGA_W,
                               NCM_PARAM_TYPE_FREE);
  /* Set w param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QGRW_W, "w", "w",
                               1e-25,  1.0, 1.0e-8,
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
  nc_hicosmo_set_Omega_c_impl   (parent_class, &_nc_hicosmo_qgrw_Omega_c);
  nc_hicosmo_set_Omega_t_impl   (parent_class, &_nc_hicosmo_qgrw_Omega_t);
  nc_hicosmo_set_xb_impl        (parent_class, &_nc_hicosmo_qgrw_xb);
  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_qgrw_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_qgrw_d2E2_dz2);
  nc_hicosmo_set_cs2_impl       (parent_class, &_nc_hicosmo_qgrw_cs2);
  nc_hicosmo_set_rhopp_impl     (parent_class, &_nc_hicosmo_qgrw_rhopp);

  object_class->finalize = nc_hicosmo_qgrw_finalize;
}

static void
nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->nuA2          = &_nc_hipert_iadiab_nuA2;
  iface->dlnmzeta      = &_nc_hipert_iadiab_dlnmzeta;
  iface->dmzetanuA_nuA = &_nc_hipert_iadiab_dmzetanuA_nuA;
  iface->eom           = &_nc_hipert_iadiab_eom;
}

static void
nc_hipert_itwo_fluids_interface_init (NcHIPertITwoFluidsInterface *iface)
{
  iface->nuA2          = (NcHIPertITwoFluidsFuncNuA2) &_nc_hipert_iadiab_nuA2;
  iface->nuB2          = &_nc_hipert_itwo_fluids_nuB2;
  iface->dmzetanuA_nuA = (NcHIPertITwoFluidsFuncDmzetanuAnuA) &_nc_hipert_iadiab_dmzetanuA_nuA;
  iface->dmSnuB_nuB    = &_nc_hipert_itwo_fluids_dmSnuB_nuB;
  iface->eom           = &_nc_hipert_itwo_fluids_eom;
}

#define VECTOR   (model->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_H0))
#define OMEGA_R  (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_OMEGA_R))
#define OMEGA_W  (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_OMEGA_W))
#define W        (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_W))
#define X_B      (ncm_vector_get (VECTOR, NC_HICOSMO_QGRW_X_B))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_qgrw_E2 (NcmModel *model, gdouble z)
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
_nc_hicosmo_qgrw_dE2_dz (NcmModel *model, gdouble z)
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
_nc_hicosmo_qgrw_d2E2_dz2 (NcmModel *model, gdouble z)
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
_nc_hicosmo_qgrw_cs2 (NcmModel *model, gdouble z)
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
 * rho plus p
 ****************************************************************************/
static gdouble
_nc_hicosmo_qgrw_rhopp (NcmModel *model, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x2 * x2;
  const gdouble w  = W;
  const gdouble x3w = pow (x3, w);

  return (1.0 + w) * OMEGA_W * x3 * x3w + (4.0 / 3.0) * OMEGA_R * x4;
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_qgrw_H0 (NcmModel *model) { return MACRO_H0 * (OMEGA_R + OMEGA_W); }
static gdouble _nc_hicosmo_qgrw_Omega_t (NcmModel *model) { return 1.0; }
static gdouble _nc_hicosmo_qgrw_Omega_c (NcmModel *model) { return OMEGA_W; }
static gdouble _nc_hicosmo_qgrw_xb (NcmModel *model) { return X_B; }

/****************************************************************************
 * Wkb theta
 ****************************************************************************/
#define _exprel(x) gsl_sf_hyperg_1F1_int (1, 2, (x))
#define _d1exprel(x) (0.5 * gsl_sf_hyperg_1F1_int (2, 3, (x)))
#define _d2exprel(x) (gsl_sf_hyperg_1F1_int (3, 4, (x)) / 3.0)
#define _d3exprel(x) (0.25 * gsl_sf_hyperg_1F1_int (4, 5, (x)))

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
  const gdouble p2f1 = 0.5 * threex1mw * _exprel (-threex1mw * alpha2_2); \
  const gdouble p2f2 = R * _exprel (-2.0 * alpha2_2); \
  const gdouble p2   = p2f1 + p2f2; \
  const gdouble p4   = onepw + onepw2 * R;

#define _NC_HICOSMO_QGRW_WKB_COMMON2 \
  const gdouble onem3w      = 1.0 - 3.0 * w; \
 \
  const gdouble d1R  = onem3w * R; \
 \
  const gdouble d1p2f1 = 0.5 * gsl_pow_2 (threex1mw) * _d1exprel (-threex1mw * alpha2_2); \
  const gdouble d1p2f2 = 2.0 * R * _d1exprel (-2.0 * alpha2_2) + d1R * _exprel (-2.0 * alpha2_2); \
  const gdouble d1p2   = d1p2f1 + d1p2f2; \
  const gdouble d1p4   = onepw2 * d1R; \
 \
  const gdouble d1lnp2 = _d1ln (p2); \
  const gdouble d1lnp4 = _d1ln (p4);

#define _NC_HICOSMO_QGRW_WKB_COMMON22 \
  const gdouble d2R    = onem3w * d1R; \
  const gdouble d2p2f1 = 0.5 * gsl_pow_3 (threex1mw) * _d2exprel (-threex1mw * alpha2_2); \
  const gdouble d2p2f2 = 4.0 * R * _d2exprel (-2.0 * alpha2_2) + 4.0 * d1R * _d1exprel (-2.0 * alpha2_2) + d2R * _exprel (-2.0 * alpha2_2); \
  const gdouble d2p2   = d2p2f1 + d2p2f2; \
  const gdouble d2p4   = onepw2 * d2R; \
 \
  const gdouble d2lnp2 = _d2ln (p2); \
  const gdouble d2lnp4 = _d2ln (p4);

#define _NC_HICOSMO_QGRW_WKB_COMMON3 \
  const gdouble p1   = w * onepw + w2 * onepw2 * R; \
  const gdouble p5   = sqrt (x3 / (x3w * Omegaw)) / 3.0;

#define _NC_HICOSMO_QGRW_WKB_COMMON4 \
  const gdouble p3       = w2 * onepw + w * onepw2 * R;

#define _NC_HICOSMO_QGRW_WKB_COMMON5 \
  const gdouble onep3w   = 1.0 + 3.0 * w; \
  const gdouble onep3w_2 = onep3w * 0.5; \
  const gdouble d1p3     = w * onepw2 * d1R; \
  const gdouble d2p3     = w * onepw2 * d2R; \
  const gdouble d1lnp3   = _d1ln (p3); \
  const gdouble d1lnp6   = onep3w_2; \
  const gdouble d2lnp3   = _d2ln (p3); \
  const gdouble d2lnp6   = 0.0;

#define _NC_HICOSMO_QGRW_WKB_COMMON6 \
  const gdouble k2          = k * k; \
  const gdouble p0   = k2 / (Omegaw * x * x3w);

#define _NC_HICOSMO_QGRW_WKB_NUZETA2 \
  const gdouble nuzeta2 = p0 * p1 / (p2 * p4);

#define _NC_HICOSMO_QGRW_WKB_DLNSQRTMZETA \
  const gdouble threex1mw_2 = threex1mw * 0.5; \
  const gdouble d1p1   = w2 * onepw2 * d1R; \
  const gdouble d1p5   = threex1mw_2 * p5; \
  const gdouble d1lnp1 = _d1ln (p1); \
  const gdouble d1lnp5 = _d1ln (p5); \
  const gdouble d1lnsqrtmzeta = d1lnp4 - 0.5 * d1lnp1 - 0.5 * d1lnp5 - 0.25 * d1lnp2 + 1.0 / alpha2;

#define _NC_HICOSMO_QGRW_WKB_VZETA \
  const gdouble d2p1   = w2 * onepw2 * d2R; \
  const gdouble d2lnp1 = _d2ln (p1); \
  const gdouble d2lnp5 = 0.0; \
  const gdouble d2lnsqrtmzeta = d2lnp4 - 0.5 * d2lnp1 - 0.5 * d2lnp5 - 0.25 * d2lnp2 + 2.0 / (alpha2 * alpha2); \
  const gdouble d2sqrtmzeta_sqrtmzeta = d1lnsqrtmzeta * d1lnsqrtmzeta + d2lnsqrtmzeta; \
  const gdouble Vzeta   = alpha2 * d2sqrtmzeta_sqrtmzeta -  d1lnsqrtmzeta;

#define _NC_HICOSMO_QGRW_WKB_NUA2 \
  const gdouble nuA2    = nuzeta2 - Vzeta;

#define _NC_HICOSMO_QGRW_WKB_MZETA \
  const gdouble mzeta   = p4 * p4 / (alpha2 * p1 * p5 * sqrt (p2));

#define _NC_HICOSMO_QGRW_WKB_NUS2 \
  const gdouble nuS2 = p0 * p3 / (p2 * p4);

#define _NC_HICOSMO_QGRW_WKB_VS \
  const gdouble d1lnsqrtmS = d1lnp4 + 0.5 * d1lnp6 + 0.25 * d1lnp2 - 0.5 * d1lnp3; \
  const gdouble d2lnsqrtmS = d2lnp4 + 0.5 * d2lnp6 + 0.25 * d2lnp2 - 0.5 * d2lnp3; \
  const gdouble d2sqrtmS_sqrtmS = d1lnsqrtmS * d1lnsqrtmS + d2lnsqrtmS; \
  const gdouble VS = alpha2 * d2sqrtmS_sqrtmS -  d1lnsqrtmS;

#define _NC_HICOSMO_QGRW_WKB_NUB2 \
  const gdouble nuB2    = nuS2 - VS;

#define _NC_HICOSMO_QGRW_WKB_MS \
  const gdouble p6      = sqrt (x3 / (x3w * Omegaw)) / (3.0 * onepw * onepw2 * R); \
  const gdouble mS      = sqrt (p2) * p4 * p4 * p6 / p3;

static gdouble
_nc_hipert_iadiab_nuA2 (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (iadiab);
  NcHICosmo *cosmo = NC_HICOSMO (iadiab);
  _NC_HICOSMO_QGRW_WKB_COMMON1;
  _NC_HICOSMO_QGRW_WKB_COMMON2;
  _NC_HICOSMO_QGRW_WKB_COMMON22;
  _NC_HICOSMO_QGRW_WKB_COMMON3;
  _NC_HICOSMO_QGRW_WKB_COMMON6;
  _NC_HICOSMO_QGRW_WKB_NUZETA2;
  _NC_HICOSMO_QGRW_WKB_DLNSQRTMZETA;
  _NC_HICOSMO_QGRW_WKB_VZETA;
  _NC_HICOSMO_QGRW_WKB_NUA2;
  
  return nuA2;
}

static gdouble
_nc_hipert_iadiab_dlnmzeta (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (iadiab);
  NcHICosmo *cosmo = NC_HICOSMO (iadiab);
  _NC_HICOSMO_QGRW_WKB_COMMON1;
  _NC_HICOSMO_QGRW_WKB_COMMON2;
  _NC_HICOSMO_QGRW_WKB_COMMON3;
  _NC_HICOSMO_QGRW_WKB_DLNSQRTMZETA;
  
  return - alpha * 2.0 * d1lnsqrtmzeta;
}

static gdouble
_nc_hipert_iadiab_dmzetanuA_nuA (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (iadiab);
  NcHICosmo *cosmo = NC_HICOSMO (iadiab);
  _NC_HICOSMO_QGRW_WKB_COMMON1;
  _NC_HICOSMO_QGRW_WKB_COMMON2;
  _NC_HICOSMO_QGRW_WKB_COMMON22;
  _NC_HICOSMO_QGRW_WKB_COMMON3;
  _NC_HICOSMO_QGRW_WKB_COMMON6;
  _NC_HICOSMO_QGRW_WKB_NUZETA2;
  _NC_HICOSMO_QGRW_WKB_DLNSQRTMZETA;
  _NC_HICOSMO_QGRW_WKB_VZETA;
  _NC_HICOSMO_QGRW_WKB_NUA2;
  _NC_HICOSMO_QGRW_WKB_MZETA;

  const gdouble onep3w  = 1.0 + 3.0 * w;
  const gdouble d1p0    = -onep3w * p0;
  const gdouble d1lnp0  = _d1ln (p0);

  const gdouble d1lnnuzeta2 = d1lnp0 + d1lnp1 - d1lnp2 - d1lnp4;

  const gdouble d3R  = onem3w * d2R;
  const gdouble d3p1   = w2 * onepw2 * d3R;
  const gdouble d3p2f1 = 0.5 * gsl_pow_4 (threex1mw) * _d3exprel (-threex1mw * alpha2_2);
  const gdouble d3p2f2 = 
    8.0 * R * _d3exprel (-2.0 * alpha2_2) + 
    12.0 * d1R * _d2exprel (-2.0 * alpha2_2) +
    6.0 * d2R * _d1exprel (-2.0 * alpha2_2) + 
    d3R * _exprel (-2.0 * alpha2_2);
  const gdouble d3p2   = d3p2f1 + d3p2f2;
  const gdouble d3p4   = onepw2 * d3R;

  const gdouble d3lnp1 = _d3ln (p1);
  const gdouble d3lnp2 = _d3ln (p2);
  const gdouble d3lnp4 = _d3ln (p4);
  const gdouble d3lnp5 = 0.0;
  
  const gdouble d3lnsqrtmzeta = d3lnp4 - 0.5 * d3lnp1 - 0.5 * d3lnp5 - 0.25 * d3lnp2 + 8.0 / (alpha2 * alpha2 * alpha2);
  
  const gdouble d1Vzeta       = alpha2 * (d3lnsqrtmzeta + 2.0 * d1lnsqrtmzeta * d2lnsqrtmzeta) - (3.0 * d2lnsqrtmzeta + 2.0 * d1lnsqrtmzeta * d1lnsqrtmzeta);
  
  const gdouble dmzetanuA_nuA = - alpha * mzeta * (2.0 * d1lnsqrtmzeta + (nuzeta2 * d1lnnuzeta2 - d1Vzeta) * 0.5 / nuA2);
  
  return dmzetanuA_nuA;
}

static gdouble 
_nc_hipert_itwo_fluids_nuB2 (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (itf);
  NcHICosmo *cosmo = NC_HICOSMO (itf);
  _NC_HICOSMO_QGRW_WKB_COMMON1;
  _NC_HICOSMO_QGRW_WKB_COMMON2;
  _NC_HICOSMO_QGRW_WKB_COMMON22;
  _NC_HICOSMO_QGRW_WKB_COMMON4;
  _NC_HICOSMO_QGRW_WKB_COMMON5;
  _NC_HICOSMO_QGRW_WKB_COMMON6;
  _NC_HICOSMO_QGRW_WKB_NUS2;
  _NC_HICOSMO_QGRW_WKB_VS;
  _NC_HICOSMO_QGRW_WKB_NUB2;

  return nuB2;
}

static gdouble 
_nc_hipert_itwo_fluids_dmSnuB_nuB (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (itf);
  NcHICosmo *cosmo = NC_HICOSMO (itf);
  _NC_HICOSMO_QGRW_WKB_COMMON1;
  _NC_HICOSMO_QGRW_WKB_COMMON2;
  _NC_HICOSMO_QGRW_WKB_COMMON22;
  _NC_HICOSMO_QGRW_WKB_COMMON4;
  _NC_HICOSMO_QGRW_WKB_COMMON5;
  _NC_HICOSMO_QGRW_WKB_COMMON6;
  _NC_HICOSMO_QGRW_WKB_NUS2;
  _NC_HICOSMO_QGRW_WKB_VS;
  _NC_HICOSMO_QGRW_WKB_NUB2;
  _NC_HICOSMO_QGRW_WKB_MS;

  const gdouble d1p0    = -onep3w * p0;
  const gdouble d1lnp0  = _d1ln (p0);

  const gdouble d1lnmS2nuS2 = (d1lnp0 + 2.0 * d1lnp6) * 0.0 + 3.0 * d1lnp4 - d1lnp3; /* p0p6 = cte */

  const gdouble d3R    = onem3w * d2R;
  const gdouble d3p3   = w * onepw2 * d3R;
  const gdouble d3p2f1 = 0.5 * gsl_pow_4 (threex1mw) * _d3exprel (-threex1mw * alpha2_2);
  const gdouble d3p2f2 = 
    8.0 * R * _d3exprel (-2.0 * alpha2_2) + 
    12.0 * d1R * _d2exprel (-2.0 * alpha2_2) +
    6.0 * d2R * _d1exprel (-2.0 * alpha2_2) + 
    d3R * _exprel (-2.0 * alpha2_2);
  const gdouble d3p2   = d3p2f1 + d3p2f2;
  const gdouble d3p4   = onepw2 * d3R;

  const gdouble d3lnp2 = _d3ln (p2);
  const gdouble d3lnp3 = _d3ln (p3);
  const gdouble d3lnp4 = _d3ln (p4);
  const gdouble d3lnp6 = 0.0;
  
  const gdouble d3lnsqrtmS = d3lnp4 + 0.5 * d3lnp6 + 0.25 * d3lnp2 - 0.5 * d3lnp3;
  
  const gdouble d1VS       = alpha2 * (d3lnsqrtmS + 2.0 * d1lnsqrtmS * d2lnsqrtmS) - (3.0 * d2lnsqrtmS + 2.0 * d1lnsqrtmS * d1lnsqrtmS);
  
  const gdouble dmSnuB_nuB = - alpha * mS * (nuS2 * d1lnmS2nuS2 - 4.0 * d1lnsqrtmS * VS - d1VS) / (2.0 * nuB2);

/*  printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", alpha, d1lnp0, 3.0 * d1lnp4, 2.0 * d1lnp6, - d1lnp3, dmSnuB_nuB);*/
  
  return dmSnuB_nuB;
}

/****************************************************************************
 * Equations of motion
 ****************************************************************************/

static NcHIPertIAdiabEOM *
_nc_hipert_iadiab_eom (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (iadiab);
  NcHICosmo *cosmo = NC_HICOSMO (iadiab);
  NcHICosmoQGRW *qgrw = NC_HICOSMO_QGRW (iadiab);
  if (qgrw->eom_adiab_zeta.skey  != NCM_MODEL (iadiab)->pkey ||
      qgrw->eom_adiab_zeta.alpha != alpha ||
      qgrw->eom_adiab_zeta.k     != k)
  {
    _NC_HICOSMO_QGRW_WKB_COMMON1;
    _NC_HICOSMO_QGRW_WKB_COMMON3;
    _NC_HICOSMO_QGRW_WKB_COMMON6;
    _NC_HICOSMO_QGRW_WKB_NUZETA2;
    _NC_HICOSMO_QGRW_WKB_MZETA;
    
    qgrw->eom_adiab_zeta.m     = mzeta;
    qgrw->eom_adiab_zeta.nu2   = nuzeta2;

    qgrw->eom_adiab_zeta.skey  = NCM_MODEL (cosmo)->pkey;
    qgrw->eom_adiab_zeta.alpha = alpha;
    qgrw->eom_adiab_zeta.k     = k;
  }
  
  return &qgrw->eom_adiab_zeta;
}

static NcHIPertITwoFluidsEOM *
_nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (itf);
  NcHICosmo *cosmo = NC_HICOSMO (itf);
  NcHICosmoQGRW *qgrw = NC_HICOSMO_QGRW (itf);
  if (qgrw->eom_two_fluids.skey  != NCM_MODEL (itf)->pkey ||
      qgrw->eom_two_fluids.alpha != alpha ||
      qgrw->eom_two_fluids.k     != k)
  {
    _NC_HICOSMO_QGRW_WKB_COMMON1;
    _NC_HICOSMO_QGRW_WKB_COMMON3;
    _NC_HICOSMO_QGRW_WKB_COMMON4;
    _NC_HICOSMO_QGRW_WKB_COMMON6;
    _NC_HICOSMO_QGRW_WKB_NUZETA2;
    _NC_HICOSMO_QGRW_WKB_MZETA;
    _NC_HICOSMO_QGRW_WKB_MS;
    _NC_HICOSMO_QGRW_WKB_NUS2;
/*
    _NC_HICOSMO_QGRW_WKB_COMMON2;
    _NC_HICOSMO_QGRW_WKB_COMMON5;
    _NC_HICOSMO_QGRW_WKB_VS;
*/
    
    qgrw->eom_two_fluids.mzeta   = mzeta;
    qgrw->eom_two_fluids.nuzeta2 = nuzeta2;

    qgrw->eom_two_fluids.mS      = mS;
    qgrw->eom_two_fluids.nuS2    = nuS2;

    qgrw->eom_two_fluids.Y       = alpha * (w - w2) * p5 / (p4 * p4 * p6);
    
    qgrw->eom_two_fluids.skey  = NCM_MODEL (cosmo)->pkey;
    qgrw->eom_two_fluids.alpha = alpha;
    qgrw->eom_two_fluids.k     = k;

    /*printf ("% 20.15g % 20.15g % 20.15g % 20.15e\n", alpha, nuS2, VS, fabs ((nuS2 - VS)/ VS));*/
  }
  
  return &qgrw->eom_two_fluids;
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
