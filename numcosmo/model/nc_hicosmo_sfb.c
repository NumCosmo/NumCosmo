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

#include "model/nc_hicosmo_sfb.h"

static void nc_hipert_adiab_interface_init (NcHIPertIAdiabInterface *iface);

G_DEFINE_TYPE_WITH_CODE (NcHICosmoSFB, nc_hicosmo_sfb, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IADIAB,
                                                nc_hipert_adiab_interface_init)
                         );

enum {
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
static gdouble _nc_hicosmo_sfb_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_sfb_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_sfb_bgp_cs2 (NcHICosmo *cosmo, gdouble z);

static gdouble _nc_hicosmo_sfb_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_sfb_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_sfb_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_sfb_xb (NcHICosmo *cosmo);

/* I must define here the functions that csq1d implements 
  static NcHIPertITwoFluidsEOM *_nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);*/
/*static NcHIPertITwoFluidsTV *_nc_hipert_itwo_fluids_tv (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);*/

static void
nc_hicosmo_sfb_class_init (NcHICosmoSFBClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
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
                               1e-50,  1.0, 1.0e-8,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_SFB_DEFAULT_W,
                               NCM_PARAM_TYPE_FIXED);
  /* Set xb param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_SFB_X_B, "x_b", "xb",
                               1.0e10,  1.0e40, 1.0e25,
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

static void
nc_hipert_adiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_mnu            = NULL;
  iface->eval_nu             = NULL;
  iface->eval_dlnmnu         = NULL;
  iface->eval_system         = NULL;
  iface->eval_powspec_factor = NULL;
}

#define VECTOR   (NCM_MODEL (cosmo)->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_H0))
#define OMEGA_R  (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_OMEGA_R))
#define OMEGA_W  (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_OMEGA_W))
#define W        (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_W))
#define X_B      (ncm_vector_get (VECTOR, NC_HICOSMO_SFB_X_B))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_sfb_E2 (NcHICosmo *cosmo, gdouble z)
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

#ifdef NC_HICOSMO_SFB_CHECK_INTERVAL
  if (G_UNLIKELY (ncm_cmp (x, xb, 1e-4) == 0))
    g_warning ("_nc_hicosmo_sfb_E2: used outside if its valid interval.");
#endif /* NC_HICOSMO_SFB_CHECK_INTERVAL */

  return (OMEGA_R * (x4 - x6 / xb2) + OMEGA_W * (x3 * x3w - x6 * xb3w / xb3));
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_sfb_dE2_dz (NcHICosmo *cosmo, gdouble z)
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

#ifdef NC_HICOSMO_SFB_CHECK_INTERVAL
  if (G_UNLIKELY (ncm_cmp (x, xb, 1e-4) == 0))
    g_warning ("_nc_hicosmo_sfb_E2: used outside if its valid interval.");
#endif /* NC_HICOSMO_SFB_CHECK_INTERVAL */

  return poly;
}

static gdouble
_nc_hicosmo_sfb_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
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

#ifdef NC_HICOSMO_SFB_CHECK_INTERVAL
  if (G_UNLIKELY (ncm_cmp (x, xb, 1e-4) == 0))
    g_warning ("_nc_hicosmo_sfb_E2: used outside if its valid interval.");
#endif /* NC_HICOSMO_SFB_CHECK_INTERVAL */

  return poly;
}

/****************************************************************************
 * Speed of sound squared
 ****************************************************************************/
static gdouble
_nc_hicosmo_sfb_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
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
static gdouble _nc_hicosmo_sfb_H0 (NcHICosmo *cosmo) { return MACRO_H0 * (OMEGA_R + OMEGA_W); }
static gdouble _nc_hicosmo_sfb_Omega_t0 (NcHICosmo *cosmo) { return 1.0; }
static gdouble _nc_hicosmo_sfb_Omega_c0 (NcHICosmo *cosmo) { return OMEGA_W; }
static gdouble _nc_hicosmo_sfb_xb (NcHICosmo *cosmo) { return X_B; }


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

