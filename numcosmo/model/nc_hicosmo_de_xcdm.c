/***************************************************************************
 *            nc_hicosmo_de_xcdm.c
 *
 *  Thu May 31 21:52:59 2007
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
 * SECTION:nc_hicosmo_de_xcdm
 * @title: NcHICosmoDEXcdm
 * @short_description: Dark Energy -- constant equation of state
 *
 * Dark Energy equation of state: $w(z) = w$, where $w$ is constant.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de_xcdm.h"
#include "model/nc_hicosmo_de_cpl.h"

G_DEFINE_TYPE (NcHICosmoDEXcdm, nc_hicosmo_de_xcdm, NC_TYPE_HICOSMO_DE);

#define VECTOR  (NCM_MODEL (cosmo_de)->params)
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define W       (ncm_vector_get (VECTOR, NC_HICOSMO_DE_XCDM_W))

static gdouble
_nc_hicosmo_de_xcdm_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x = 1.0 + z;  
  return OMEGA_X * pow (x, 3.0 * ( 1.0 + W ) );
}

static gdouble
_nc_hicosmo_de_xcdm_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble E2Omega_de = OMEGA_X * pow (x, 3.0 * ( 1.0 + W ) );

  return 3.0 * ( 1.0 + W ) / x * E2Omega_de;
}

static gdouble
_nc_hicosmo_de_xcdm_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble E2Omega_de = OMEGA_X * pow (x, 3.0 * ( 1.0 + W ));

  return 3.0 * ( 1.0 + W ) * (2.0 + 3.0 * W) / x2 * E2Omega_de;
}

static gdouble _nc_hicosmo_de_xcdm_w_de (NcHICosmoDE *cosmo_de, gdouble z) { return W; }

/**
 * nc_hicosmo_de_xcdm_new:
 *
 * This function instantiates a new object of type #NcHICosmoDEXcdm.
 *
 * Returns: A new #NcHICosmoDEXcdm
 */
NcHICosmoDEXcdm *
nc_hicosmo_de_xcdm_new (void)
{
  NcHICosmoDEXcdm *xcdm = g_object_new (NC_TYPE_HICOSMO_DE_XCDM, NULL);
  return xcdm;
}

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_xcdm_init (NcHICosmoDEXcdm *xcdm)
{
  NCM_UNUSED (xcdm);
}

static void
nc_hicosmo_de_xcdm_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_xcdm_parent_class)->finalize (object);
}

static void
nc_hicosmo_de_xcdm_class_init (NcHICosmoDEXcdmClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoDEClass* parent_class = NC_HICOSMO_DE_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize     = &nc_hicosmo_de_xcdm_finalize;

  nc_hicosmo_de_set_E2Omega_de_impl (parent_class, &_nc_hicosmo_de_xcdm_E2Omega_de);
  nc_hicosmo_de_set_dE2Omega_de_dz_impl (parent_class, &_nc_hicosmo_de_xcdm_dE2Omega_de_dz);
  nc_hicosmo_de_set_d2E2Omega_de_dz2_impl (parent_class, &_nc_hicosmo_de_xcdm_d2E2Omega_de_dz2);
  nc_hicosmo_de_set_w_de_impl (parent_class, &_nc_hicosmo_de_xcdm_w_de);

  ncm_model_class_set_name_nick (model_class, "XCDM - Constant EOS", "XCDM");
  ncm_model_class_add_params (model_class, 1, 0, PROP_SIZE);
  /* Set w_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_XCDM_W, "w", "w",
                               -5.0, 0.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_CPL_DEFAULT_W0,
                               NCM_PARAM_TYPE_FREE);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}
