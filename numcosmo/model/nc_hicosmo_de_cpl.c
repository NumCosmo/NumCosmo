/***************************************************************************
 *            nc_hicosmo_de_cpl.c
 *
 *  Tue Jul 31 01:37:57 2007
 *  Copyright  2007  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_hicosmo_de_cpl
 * @title: NcHICosmoDECpl
 * @short_description: Dark Energy -- Chevallier–Polarski–Linder equation of state
 *
 * See [Chevallier (2001)][XChevallier2001] and [Linder (2003)][XLinder2003]: $ w(z) = w_0 + w_1 \frac{z}{1.0 + z}$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de_cpl.h"
#include "nc_hiprim_power_law.h"
#include "nc_hireion_camb.h"
#include "nc_powspec_ml_cbe.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcHICosmoDECpl, nc_hicosmo_de_cpl, NC_TYPE_HICOSMO_DE)

#define VECTOR (ncm_model_orig_params_peek_vector (NCM_MODEL (cosmo_de)))
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define OMEGA_0 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_CPL_W0))
#define OMEGA_1 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_CPL_W1))

static gdouble
_nc_hicosmo_de_cpl_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  gdouble x   = 1.0 + z;
  gdouble lnx = log1p (z);

  return OMEGA_X * exp (-3.0 * OMEGA_1 * z / x + 3.0 * (1.0 + OMEGA_0 + OMEGA_1) * lnx);
}

static gdouble
_nc_hicosmo_de_cpl_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x          = 1.0 + z;
  const gdouble x2         = x * x;
  const gdouble lnx        = log1p (z);
  const gdouble E2Omega_de = OMEGA_X * exp (-3.0 * OMEGA_1 * z / x + 3.0 * (1.0 + OMEGA_0 + OMEGA_1) * lnx);

  return 3.0 * ((x * (OMEGA_0 + 1.0) + z * OMEGA_1) / x2) * E2Omega_de;
}

static gdouble
_nc_hicosmo_de_cpl_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x          = 1.0 + z;
  const gdouble x2         = x * x;
  const gdouble x4         = x2 * x2;
  const gdouble lnx        = log1p (z);
  const gdouble E2Omega_de = OMEGA_X * exp (-3.0 * OMEGA_1 * z / x + 3.0 * (1.0 + OMEGA_0 + OMEGA_1) * lnx);
  const gdouble w0         = OMEGA_0;
  const gdouble w1         = OMEGA_1;
  const gdouble w12        = w1 * w1;

  return (9.0 * w12 - 6.0 * w1 * (2.0 + 3.0 * w0 + 3.0 * w1) * x + 3.0 * (1.0 + w0 + w1) * (2.0 + 3.0 * w0 + 3.0 * w1) * x2) * E2Omega_de / x4;
}

static gdouble
_nc_hicosmo_de_cpl_w_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble w0 = OMEGA_0;
  const gdouble w1 = OMEGA_1;

  return w0 + w1 * z / (1.0 + z);
}

/**
 * nc_hicosmo_de_cpl_new:
 *
 * This function instantiates a new object of type #NcHICosmoDECpl.
 *
 * Returns: A new #NcHICosmoDECpl
 */
NcHICosmoDECpl *
nc_hicosmo_de_cpl_new (void)
{
  NcHICosmoDECpl *cpl = g_object_new (NC_TYPE_HICOSMO_DE_CPL, NULL);

  return cpl;
}

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_cpl_init (NcHICosmoDECpl *cpl)
{
  NCM_UNUSED (cpl);
}

static void
nc_hicosmo_de_cpl_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_cpl_parent_class)->finalize (object);
}

static void
nc_hicosmo_de_cpl_class_init (NcHICosmoDECplClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcHICosmoDEClass *parent_class = NC_HICOSMO_DE_CLASS (klass);
  NcmModelClass *model_class     = NCM_MODEL_CLASS (klass);

  object_class->finalize = &nc_hicosmo_de_cpl_finalize;

  ncm_model_class_set_name_nick (model_class, "Chevalier-Polarski-Linder parametrization", "CPL");
  ncm_model_class_add_params (model_class, 2, 0, PROP_SIZE);
  /* Set w_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_CPL_W0, "w_0", "w0",
                              -10.0, 1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_CPL_DEFAULT_W0,
                              NCM_PARAM_TYPE_FREE);
  /* Set w_1 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_CPL_W1, "w_1", "w1",
                              -5.0, 5.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_CPL_DEFAULT_W1,
                              NCM_PARAM_TYPE_FREE);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_de_set_E2Omega_de_impl (parent_class,       &_nc_hicosmo_de_cpl_E2Omega_de);
  nc_hicosmo_de_set_dE2Omega_de_dz_impl (parent_class,   &_nc_hicosmo_de_cpl_dE2Omega_de_dz);
  nc_hicosmo_de_set_d2E2Omega_de_dz2_impl (parent_class, &_nc_hicosmo_de_cpl_d2E2Omega_de_dz2);
  nc_hicosmo_de_set_w_de_impl (parent_class,             &_nc_hicosmo_de_cpl_w_de);
}

