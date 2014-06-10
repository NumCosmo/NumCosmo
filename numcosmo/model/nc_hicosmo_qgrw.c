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
 * In this model the adiabatic mode $\zeta$ has its mass, speed of sound square $c_s^2$ and frequency square $\mu_\zeta^2$ given by
 * \begin{align}
 * m_\zeta &= 3 \Delta_\bar{K}\sqrt{\Omega_w} x^{-3(1-w)/2}\frac{(1 + w) +  4R/3}{c_s^2}\frac{1}{\sqrt{(1-exp(-2\vert\alpha\vert)) + (1-exp(-3(1-w)\vert\alpha\vert))}}, \\\\
 * c_s^2 &= \frac{w (1 + w) + 4R/9}{(1+w) + 4R/3}, \\\\
 * \mu_\zeta^2 &= \frac {c_s^2 k^2}{\Omega_w x^{1+3w} ((1-exp(-2\vert\alpha\vert)) + (1-exp(-3(1-w)\vert\alpha\vert)))}, 
 * \end{align}
 * where $$R \equiv \frac{\Omega_r x}{\Omega_w x^{3w}}.$$
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qgrw.h"

G_DEFINE_TYPE (NcHICosmoQGRW, nc_hicosmo_qgrw, NC_TYPE_HICOSMO);

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_qgrw_init (NcHICosmoQGRW *qgrw)
{
  memset (&qgrw->eom_adiab_zeta, 0, sizeof (NcHICosmoEOMAdiabZeta));
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

static gdouble _nc_hicosmo_qgrw_wkb_adiab_theta (NcmModel *model, gdouble alpha, gdouble k);
static gdouble _nc_hicosmo_qgrw_wkb_adiab_dmtheta (NcmModel *model, gdouble alpha, gdouble k);

static NcHICosmoEOMAdiabZeta *_nc_hicosmo_qgrw_eom_adiab_zeta (NcHICosmo *cosmo, gdouble alpha, gdouble k);

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

  nc_hicosmo_set_wkb_adiab_theta_impl (parent_class, _nc_hicosmo_qgrw_wkb_adiab_theta);
  nc_hicosmo_set_wkb_adiab_dmtheta_impl (parent_class, _nc_hicosmo_qgrw_wkb_adiab_dmtheta);
  
  nc_hicosmo_set_eom_adiab_zeta_impl (parent_class, &_nc_hicosmo_qgrw_eom_adiab_zeta);

  object_class->finalize = nc_hicosmo_qgrw_finalize;
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
#define _NC_HICOSMO_QGRW_WKB_COMMON \
  const gdouble abs_alpha = fabs (alpha); \
  const gdouble xb        = X_B; \
  const gdouble w         = W; \
  const gdouble Omegar    = OMEGA_R; \
  const gdouble Omegaw    = OMEGA_W; \
 \
  const gdouble x   = xb * exp (-abs_alpha); \
  const gdouble x2  = x * x; \
  const gdouble x3  = x2 * x; \
  const gdouble x3w = pow (x3, w); \
 \
  const gdouble xb2 = xb * xb; \
  const gdouble xb3 = xb2 * xb; \
  const gdouble xb3w = pow (xb3, w); \
 \
  const gdouble k2          = k * k; \
  const gdouble onepw       = 1.0 + w; \
  const gdouble onem3w      = 1.0 - 3.0 * w; \
  const gdouble threex1mw   = 3.0 * (1.0 - w); \
  const gdouble four_thirds = 4.0 / 3.0; \
  const gdouble one_third   = 1.0 / 3.0; \
 \
  const gdouble R      = Omegar * x / (Omegaw * x3w); \
  const gdouble frhopp = onepw + four_thirds * R; \
  const gdouble cs2    = (w * onepw + one_third * four_thirds * R) / frhopp; \
  const gdouble fE2    = -(R * expm1 (-2.0 * abs_alpha) + expm1 (-threex1mw * abs_alpha)); \
 \
  /* m = 3.0 * DeltaK * sqrt (Omegaw) * A * B^2 / (C * sqrt(D)) */ \
 \
  /* const gdouble A = 1.0 / x3x1mw_2; */ \
  const gdouble B = frhopp; \
  const gdouble C = (w * onepw + one_third * four_thirds * R); \
  const gdouble D = fE2; \
   \
  const gdouble dR_dlnx    = onem3w * R; \
  const gdouble d2R_dlnx2  = onem3w * dR_dlnx; \
 \
  const gdouble dlnA_dlnx = -threex1mw * 0.5; \
 \
  const gdouble dB_dlnx     = four_thirds * dR_dlnx; \
  const gdouble d2B_dlnx2   = four_thirds * d2R_dlnx2; \
  const gdouble dlnB_dlnx   = dB_dlnx / B; \
  const gdouble d2lnB_dlnx2 = -(dlnB_dlnx * dlnB_dlnx) + d2B_dlnx2 / B; \
 \
  const gdouble dC_dlnx     = one_third * dB_dlnx; \
  const gdouble d2C_dlnx2   = one_third * d2B_dlnx2; \
  const gdouble dlnC_dlnx   = dC_dlnx / C; \
  const gdouble d2lnC_dlnx2 = -(dlnC_dlnx * dlnC_dlnx) + d2C_dlnx2 / C; \
 \
  const gdouble dD_dlnx     = dR_dlnx - 3.0 * onepw * R * x2 / xb2 - threex1mw * x3 * xb3w / (xb3 * x3w); \
  const gdouble d2D_dlnx2   = d2R_dlnx2 - 9.0 * onepw * onepw * R * x2 / xb2 - threex1mw * threex1mw * x3 * xb3w / (xb3 * x3w); \
  const gdouble dlnD_dlnx   = dD_dlnx / D; \
  const gdouble d2lnD_dlnx2 = -(dlnD_dlnx * dlnD_dlnx) + d2D_dlnx2 / D; \
 \
  const gdouble dlnm_dlnx   = dlnA_dlnx + 2.0 * dlnB_dlnx - dlnC_dlnx - 0.5 * dlnD_dlnx; \
  const gdouble d2lnm_dlnx2 = 2.0 * d2lnB_dlnx2 - d2lnC_dlnx2 - 0.5 * d2lnD_dlnx2; \
 \
  const gdouble sqrtmpp_sqrtm = 0.25 * gsl_pow_2 (dlnm_dlnx) + 0.5 * d2lnm_dlnx2;


static gdouble
_nc_hicosmo_qgrw_wkb_adiab_theta (NcmModel *model, gdouble alpha, gdouble k)
{
  _NC_HICOSMO_QGRW_WKB_COMMON;  
  const gdouble mu2 = cs2 * k2 / (Omegaw * x * x3w * fE2);  
  const gdouble theta2 = mu2 - sqrtmpp_sqrtm;

  /* printf ("# WKB alpha % 20.15g x % 20.15g cs2 % 20.15g k2 % 20.15g n2cs2x2k2 % 20.15g pot % 20.15g\n", alpha, x, cs2, k2, cs2 * k2 / (Omegaw * x * x3w * fE2), sqrtmpp_sqrtm); */
  
  return sqrt (theta2);
}

static gdouble
_nc_hicosmo_qgrw_wkb_adiab_dmtheta (NcmModel *model, gdouble alpha, gdouble k)
{
  _NC_HICOSMO_QGRW_WKB_COMMON;

  const gdouble onep3w = 1.0 + 3.0 * w;
  const gdouble dlnmu2_dlnx = (dlnC_dlnx - dlnB_dlnx - onep3w - dlnD_dlnx);
  
  const gdouble mu2 = cs2 * k2 / (Omegaw * x * x3w * fE2);  
  const gdouble theta2 = mu2 - sqrtmpp_sqrtm;

  const gdouble Omegak = nc_hicosmo_Omega_k (NC_HICOSMO (model));
  const gdouble DeltaK = k2 / (k2 + Omegak);
  const gdouble A      = sqrt (x3w / x3);
  const gdouble m      = 3.0 * DeltaK * sqrt (Omegaw) * A * B * B / (C * sqrt (D));

  const gdouble d3R_dlnx3   = onem3w * d2R_dlnx2;

  const gdouble d3B_dlnx3   = four_thirds * d3R_dlnx3;
  const gdouble d3C_dlnx3   = one_third * d3B_dlnx3;
  const gdouble d3D_dlnx3   = d3R_dlnx3 - 27.0 * onepw * onepw * onepw * R * x2 / xb2 - threex1mw * threex1mw * threex1mw * x3 * xb3w / (xb3 * x3w);

  const gdouble d3lnB_dlnx3 = d3B_dlnx3 / B - gsl_pow_3 (dlnB_dlnx) - 3.0 * dlnB_dlnx * d2lnB_dlnx2;
  const gdouble d3lnC_dlnx3 = d3C_dlnx3 / C - gsl_pow_3 (dlnC_dlnx) - 3.0 * dlnC_dlnx * d2lnC_dlnx2;
  const gdouble d3lnD_dlnx3 = d3D_dlnx3 / D - gsl_pow_3 (dlnD_dlnx) - 3.0 * dlnD_dlnx * d2lnD_dlnx2;

  const gdouble d3lnm_dlnx3 = 2.0 * d3lnB_dlnx3 - d3lnC_dlnx3 - 0.5 * d3lnD_dlnx3;
  
  const gdouble dmtheta_theta = - GSL_SIGN (alpha) * m * (dlnm_dlnx + 0.5 / theta2 * (mu2 * dlnmu2_dlnx - 0.5 * dlnm_dlnx * d2lnm_dlnx2 - 0.5 * d3lnm_dlnx3));

  /* printf ("# WKB alpha % 20.15g x % 20.15g cs2 % 20.15g k2 % 20.15g n2cs2x2k2 % 20.15g pot % 20.15g\n", alpha, x, cs2, k2, cs2 * k2 / (Omegaw * x * x3w * fE2), sqrtmpp_sqrtm); */
  
  return dmtheta_theta;
}


/****************************************************************************
 * Equations of motion
 ****************************************************************************/

static NcHICosmoEOMAdiabZeta *
_nc_hicosmo_qgrw_eom_adiab_zeta (NcHICosmo *cosmo, gdouble alpha, gdouble k)
{
  NcmModel *model = NCM_MODEL (cosmo);
  NcHICosmoQGRW *qgrw = NC_HICOSMO_QGRW (cosmo);
  if (qgrw->eom_adiab_zeta.skey  != NCM_MODEL (cosmo)->pkey ||
      qgrw->eom_adiab_zeta.alpha != alpha ||
      qgrw->eom_adiab_zeta.k     != k)
  {
    const gdouble abs_alpha = fabs (alpha);
    const gdouble xb        = X_B;
    const gdouble w         = W;
    const gdouble Omegar    = OMEGA_R;
    const gdouble Omegaw    = OMEGA_W;
    const gdouble sOmegaw   = sqrt (Omegaw);
    const gdouble Omegak    = nc_hicosmo_Omega_k (cosmo);

    const gdouble x   = xb * exp (-abs_alpha);
    const gdouble x2  = x * x;
    const gdouble x3  = x2 * x;
    const gdouble x3w = pow (x3, w);
    const gdouble x3x1mw_2 = sqrt (x3 / x3w);

    const gdouble k2          = k * k;
    const gdouble onepw       = 1.0 + w;
    const gdouble four_thirds = 4.0 / 3.0;
    const gdouble one_third   = 1.0 / 3.0;
    
    const gdouble R      = Omegar * x / (Omegaw * x3w);
    const gdouble frhopp = onepw + four_thirds * R;
    const gdouble cs2    = (w * onepw + one_third * four_thirds * R) / frhopp;
    const gdouble fE2    = -(R * expm1 (-2.0 * abs_alpha) + expm1 (-3.0 * (1.0 - w) * abs_alpha));
    const gdouble fE     = sqrt (fE2);

    const gdouble DeltaK = k2 / (k2 + Omegak);

    qgrw->eom_adiab_zeta.m     = 3.0 * DeltaK * (sOmegaw / x3x1mw_2) * frhopp / (cs2 * fE);
    qgrw->eom_adiab_zeta.mu2   = cs2 * k2 / (Omegaw * x * x3w * fE2);

    qgrw->eom_adiab_zeta.skey  = NCM_MODEL (cosmo)->pkey;
    qgrw->eom_adiab_zeta.alpha = alpha;
    qgrw->eom_adiab_zeta.k     = k;
  }
  
  return &qgrw->eom_adiab_zeta;
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
