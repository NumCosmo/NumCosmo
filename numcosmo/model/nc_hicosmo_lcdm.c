/***************************************************************************
 *            nc_hicosmo_lcdm.c
 *
 *  Thu May 31 21:52:45 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * SECTION:nc_hicosmo_lcdm
 * @title: NcHICosmoLCDM
 * @short_description: $\Lambda$CDM model
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_lcdm.h"
#include "model/nc_hicosmo_de.h"

G_DEFINE_TYPE (NcHICosmoLCDM, nc_hicosmo_lcdm, NC_TYPE_HICOSMO)


enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_lcdm_init (NcHICosmoLCDM *lcdm)
{
  NCM_UNUSED (lcdm);
}

static void
nc_hicosmo_lcdm_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_lcdm_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_lcdm_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_lcdm_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_lcdm_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_lcdm_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Omega_r0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Omega_b0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Omega_g0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Omega_nu0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_T_gamma0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_Yp_4He (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_lcdm_bgp_cs2 (NcHICosmo *cosmo, gdouble z);

static void
nc_hicosmo_lcdm_class_init (NcHICosmoLCDMClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->finalize = &nc_hicosmo_lcdm_finalize;

  ncm_model_class_set_name_nick (model_class, "\\Lambda{}CDM", "LCDM");
  ncm_model_class_add_params (model_class, NC_HICOSMO_DE_SPARAM_LEN, 0, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_H0, "H_0", "H0",
                              31.0, 99.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_c0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_C, "\\Omega_{c0}", "Omegac",
                              0.01,  0.9, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_C,
                              NCM_PARAM_TYPE_FREE);
  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_X, "\\Omega_{x0}", "Omegax",
                              0.01,  2.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_X,
                              NCM_PARAM_TYPE_FREE);
  /* Set T_gamma0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_T_GAMMA0, "T_{\\gamma0}", "Tgamma0",
                              2.0, 3.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_T_GAMMA0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set He Yp param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_HE_YP, "Y_p", "Yp",
                              0.0,  1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_HE_YP,
                              NCM_PARAM_TYPE_FIXED);
  /* Set ENnu param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_ENNU, "N_\\nu", "ENnu",
                              0.0,  10.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_ENNU,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_b0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_B, "\\Omega_{b0}", "Omegab",
                              0.03,  0.05, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_B,
                              NCM_PARAM_TYPE_FIXED);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_lcdm_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_lcdm_E2);
  nc_hicosmo_set_Omega_c0_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_c0);
  nc_hicosmo_set_Omega_r0_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_r0);
  nc_hicosmo_set_Omega_b0_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_b0);
  nc_hicosmo_set_Omega_g0_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_g0);
  nc_hicosmo_set_Omega_nu0_impl  (parent_class, &_nc_hicosmo_lcdm_Omega_nu0);
  nc_hicosmo_set_Omega_t0_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_t0);
  nc_hicosmo_set_T_gamma0_impl  (parent_class, &_nc_hicosmo_lcdm_T_gamma0);
  nc_hicosmo_set_Yp_4He_impl    (parent_class, &_nc_hicosmo_lcdm_Yp_4He);

  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_lcdm_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_lcdm_d2E2_dz2);

  nc_hicosmo_set_bgp_cs2_impl   (parent_class, &_nc_hicosmo_lcdm_bgp_cs2);
}

#define VECTOR   (NCM_MODEL (cosmo))
#define MACRO_H0 (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_H0))
#define OMEGA_C  (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_OMEGA_C))
#define OMEGA_X  (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define T_GAMMA0 (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_T_GAMMA0))
#define HE_YP    (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_HE_YP))
#define ENNU     (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_ENNU))
#define OMEGA_R  nc_hicosmo_Omega_r0 (cosmo)
#define OMEGA_B  (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_OMEGA_B))
#define OMEGA_M  (OMEGA_B + OMEGA_C)
#define OMEGA_K  (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_lcdm_E2 (NcHICosmo *cosmo, gdouble z)
{
  gdouble omega_k  = OMEGA_K;
  const gdouble x  = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x3 * x;

  return (OMEGA_R * x4 + OMEGA_M * x3 + omega_k * x2 + OMEGA_X);
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_lcdm_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x       = 1.0 + z;
  const gdouble x2      = x * x;
  const gdouble x3      = x2 * x;
  const gdouble poly    = 4.0 * OMEGA_R * x3 +
                          3.0 * OMEGA_M * x2 +
                          2.0 * omega_k * x;

  return poly;
}

static gdouble
_nc_hicosmo_lcdm_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x       = 1.0 + z;
  const gdouble x2      = x * x;
  const gdouble poly    = 12.0 * OMEGA_R * x2 +
                          6.0 * OMEGA_M * x +
                          2.0 * omega_k;

  return poly;
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_lcdm_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0;
}

static gdouble
_nc_hicosmo_lcdm_Omega_t0 (NcHICosmo *cosmo)
{
  return OMEGA_M + OMEGA_X + OMEGA_R;
}

static gdouble
_nc_hicosmo_lcdm_Omega_c0 (NcHICosmo *cosmo)
{
  return OMEGA_C;
}

static gdouble
_nc_hicosmo_lcdm_T_gamma0 (NcHICosmo *cosmo)
{
  return T_GAMMA0;
}

static gdouble
_nc_hicosmo_lcdm_Yp_4He (NcHICosmo *cosmo)
{
  return HE_YP;
}

static gdouble
_nc_hicosmo_lcdm_Omega_g0 (NcHICosmo *cosmo)
{
  const gdouble h  = MACRO_H0 / 100.0;
  const gdouble h2 = h * h;

  return ncm_c_radiation_temp_to_h2Omega_r0 (T_GAMMA0) / h2;
}

static gdouble
_nc_hicosmo_lcdm_Omega_nu0 (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);

  return ENNU * conv * _nc_hicosmo_lcdm_Omega_g0 (cosmo);
}

static gdouble
_nc_hicosmo_lcdm_Omega_r0 (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);

  return (1.0 + ENNU * conv) * _nc_hicosmo_lcdm_Omega_g0 (cosmo);
}

static gdouble
_nc_hicosmo_lcdm_Omega_b0 (NcHICosmo *cosmo)
{
  return OMEGA_B;
}

static gdouble
_nc_hicosmo_lcdm_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x        = 1.0 + z;
  const gdouble Omega_g0 = _nc_hicosmo_lcdm_Omega_g0 (cosmo);
  const gdouble Omega_b0 = _nc_hicosmo_lcdm_Omega_b0 (cosmo);
  const gdouble nine_4   = 9.0 / 4.0;

  return 1.0 / (3.0 + nine_4 * Omega_b0 / (Omega_g0 * x));
}

/**
 * nc_hicosmo_lcdm_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoLCDM *
nc_hicosmo_lcdm_new (void)
{
  NcHICosmoLCDM *lcdm = g_object_new (NC_TYPE_HICOSMO_LCDM, NULL);

  return lcdm;
}

