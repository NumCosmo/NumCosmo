/***************************************************************************
 *            nc_hicosmo_lcdm.c
 *
 *  Thu May 31 21:52:45 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * @title: $\Lambda$CDM
 * @short_description: Implementation of $\Lambda$CDM model
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

G_DEFINE_TYPE (NcHICosmoLCDM, nc_hicosmo_lcdm, NC_TYPE_HICOSMO);

#define VECTOR    (model->params)
#define MACRO_H0  (ncm_vector_get (VECTOR, NC_HICOSMO_DE_H0))
#define OMEGA_C   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_C))
#define OMEGA_X   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define T_GAMMA0  (ncm_vector_get (VECTOR, NC_HICOSMO_DE_T_GAMMA0))
#define OMEGA_R   nc_hicosmo_Omega_r (NC_HICOSMO (model))
#define OMEGA_B   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_B))
#define SPECINDEX (ncm_vector_get (VECTOR, NC_HICOSMO_DE_SPECINDEX))
#define SIGMA8    (ncm_vector_get (VECTOR, NC_HICOSMO_DE_SIGMA8))

#define OMEGA_M (OMEGA_B + OMEGA_C)
#define OMEGA_K (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_lcdm_E2 (NcmModel *model, gdouble z)
{
  gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x*x;
  const gdouble x3 = x2*x;
  const gdouble x4 = x3*x;

  return (OMEGA_R * x4 + OMEGA_M * x3 + omega_k * x2 + OMEGA_X);
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_lcdm_dE2_dz (NcmModel *model, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble poly = 4.0 * OMEGA_R * x3 +
	3.0 * OMEGA_M * x2 +
	2.0 * omega_k * x;
  return poly;
}

static gdouble
_nc_hicosmo_lcdm_d2E2_dz2 (NcmModel *model, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble poly = 12.0 * OMEGA_R * x2 +
	6.0 * OMEGA_M * x +
	2.0 * omega_k;
  return poly;
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_lcdm_H0 (NcmModel *model) { return MACRO_H0; }
static gdouble _nc_hicosmo_lcdm_Omega_t (NcmModel *model) { return OMEGA_M + OMEGA_X + OMEGA_R; }
static gdouble _nc_hicosmo_lcdm_Omega_c (NcmModel *model) { return OMEGA_C; }
static gdouble _nc_hicosmo_lcdm_T_gamma0 (NcmModel *model) { return T_GAMMA0; }
static gdouble
_nc_hicosmo_lcdm_Omega_r (NcmModel *model)
{
  const gdouble h = MACRO_H0 / 100.0;
  const gdouble h2 = h * h;
  return (1.0 + 0.2271 * ncm_c_neutrino_n_eff ()) * ncm_c_radiation_temp_to_h2omega_r (T_GAMMA0) / h2;
}
static gdouble _nc_hicosmo_lcdm_Omega_b (NcmModel *model) { return OMEGA_B; }
static gdouble _nc_hicosmo_lcdm_sigma_8 (NcmModel *model) { return SIGMA8; }
static gdouble _nc_hicosmo_lcdm_powspec (NcmModel *model, gdouble k) { return pow (k, SPECINDEX); }

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

enum {
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

static void
nc_hicosmo_lcdm_class_init (NcHICosmoLCDMClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class   = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize     = &nc_hicosmo_lcdm_finalize;

  ncm_model_class_add_params (model_class, 7, 0, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "\\Lambda{}CDM", "LCDM");

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_c param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_C, "\\Omega_c", "Omegac",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_C,
                               NCM_PARAM_TYPE_FREE);
  /* Set Omega_x param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_X, "\\Omega_x", "Omegax",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_X,
                               NCM_PARAM_TYPE_FREE);
  /* Set T_gamma0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_T_GAMMA0, "T_{\\gamma0}", "Tgamma0",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_T_GAMMA0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_b param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_B, "\\Omega_b", "Omegab",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_B,
                               NCM_PARAM_TYPE_FIXED);
  /* Set n_s param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_SPECINDEX, "n_s", "ns",
                               0.5,   1.5, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_SPECINDEX,
                               NCM_PARAM_TYPE_FIXED);
  /* Set sigma_8 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_SIGMA8, "\\sigma_8", "sigma8",
                               0.2,   1.8, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_SIGMA8,
                               NCM_PARAM_TYPE_FIXED);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_lcdm_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_lcdm_E2);
  nc_hicosmo_set_Omega_c_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_c);
  nc_hicosmo_set_Omega_r_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_r);
  nc_hicosmo_set_Omega_b_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_b);
  nc_hicosmo_set_Omega_t_impl   (parent_class, &_nc_hicosmo_lcdm_Omega_t);
  nc_hicosmo_set_sigma_8_impl   (parent_class, &_nc_hicosmo_lcdm_sigma_8);
  nc_hicosmo_set_T_gamma0_impl  (parent_class, &_nc_hicosmo_lcdm_T_gamma0);

  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_lcdm_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_lcdm_d2E2_dz2);

  nc_hicosmo_set_powspec_impl   (parent_class, &_nc_hicosmo_lcdm_powspec);
}
