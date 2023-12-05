/***************************************************************************
 *            nc_hicosmo_de_jbp.c
 *
 *  Thu Dec  6 15:49:35 2007
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
 * SECTION:nc_hicosmo_de_jbp
 * @title: NcHICosmoDEJbp
 * @short_description: Dark Energy -- Jassal-Bagla-Padmanabhan equation of state
 *
 * See [Jassal et al. (2005)][XJassal2005]: $w(z) = w_0 + w_1 \frac{z}{(1.0 + z)^2}$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de_jbp.h"

G_DEFINE_TYPE (NcHICosmoDEJbp, nc_hicosmo_de_jbp, NC_TYPE_HICOSMO_DE)


enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_jbp_init (NcHICosmoDEJbp *jbp)
{
  NCM_UNUSED (jbp);
}

static void
nc_hicosmo_de_jbp_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_jbp_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_de_jbp_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_jbp_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_jbp_w_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_jbp_dw_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_jbp_ln_rho_rho0 (NcHICosmoDE *cosmo_de, gdouble z);

static void
nc_hicosmo_de_jbp_class_init (NcHICosmoDEJbpClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcHICosmoDEClass *parent_class = NC_HICOSMO_DE_CLASS (klass);
  NcmModelClass *model_class     = NCM_MODEL_CLASS (klass);

  object_class->finalize = &nc_hicosmo_de_jbp_finalize;

  ncm_model_class_set_name_nick (model_class, "JBP parametrization", "JBP");
  ncm_model_class_add_params (model_class, 2, 0, PROP_SIZE);
  /* Set w_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_JBP_W0, "w_0", "w0",
                              -10.0, 1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_JBP_DEFAULT_W0,
                              NCM_PARAM_TYPE_FREE);
  /* Set w_1 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_JBP_W1, "w_1", "w1",
                              -5.0, 5.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_JBP_DEFAULT_W1,
                              NCM_PARAM_TYPE_FREE);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_de_set_E2Omega_de_impl (parent_class, &_nc_hicosmo_de_jbp_E2Omega_de);
  nc_hicosmo_de_set_dE2Omega_de_dz_impl (parent_class, &_nc_hicosmo_de_jbp_dE2Omega_de_dz);
  nc_hicosmo_de_set_w_de_impl (parent_class, &_nc_hicosmo_de_jbp_w_de);
  nc_hicosmo_de_set_dw_de_impl (parent_class, &_nc_hicosmo_de_jbp_dw_de);
  nc_hicosmo_de_set_ln_rho_rho0_impl (parent_class, &_nc_hicosmo_de_jbp_ln_rho_rho0);
}

#define VECTOR  (NCM_MODEL (cosmo_de)->params)
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define OMEGA_0 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_JBP_W0))
#define OMEGA_1 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_JBP_W1))

static gdouble
_nc_hicosmo_de_jbp_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  gdouble x   = 1.0 + z;
  gdouble lnx = log1p (z);

  return OMEGA_X * exp (3.0 / 2.0 * OMEGA_1 * gsl_pow_2 (z / x) + 3.0 * (1.0 + OMEGA_0) * lnx);
}

static gdouble
_nc_hicosmo_de_jbp_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x          = 1.0 + z;
  const gdouble x3         = gsl_pow_3 (x);
  const gdouble lnx        = log1p (z);
  const gdouble E2Omega_de = OMEGA_X * exp (3.0 / 2.0 * OMEGA_1 * gsl_pow_2 (z / x) + 3.0 * (1.0 + OMEGA_0) * lnx);

  return 3.0 * ((1.0 + OMEGA_0) / x + z * OMEGA_1 / x3) * E2Omega_de;
}

static gdouble
_nc_hicosmo_de_jbp_w_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble w0 = OMEGA_0;
  const gdouble w1 = OMEGA_1;

  return w0 + w1 * z / gsl_pow_2 (1.0 + z);
}

static gdouble
_nc_hicosmo_de_jbp_dw_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble w1 = OMEGA_1;

  return w1 * (1.0 - z) / gsl_pow_3 (1.0 + z);
}

static gdouble
_nc_hicosmo_de_jbp_ln_rho_rho0 (NcHICosmoDE *cosmo_de, gdouble z)
{
  gdouble x   = 1.0 + z;
  gdouble lnx = log1p (z);

  return 3.0 / 2.0 * OMEGA_1 * gsl_pow_2 (z / x) + 3.0 * (1.0 + OMEGA_0) * lnx;
}

/**
 * nc_hicosmo_de_jbp_new:
 *
 * This function instantiates a new object of type #NcHICosmoDEJbp.
 *
 * Returns: A new #NcHICosmoDEJbp
 */
NcHICosmoDEJbp *
nc_hicosmo_de_jbp_new (void)
{
  NcHICosmoDEJbp *jbp = g_object_new (NC_TYPE_HICOSMO_DE_JBP, NULL);

  return jbp;
}

