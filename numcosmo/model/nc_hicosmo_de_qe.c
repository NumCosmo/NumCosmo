/***************************************************************************
 *            nc_hicosmo_de_qe.c
 *
 *  Fri Dec  7 15:17:38 2007
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
 * SECTION:nc_hicosmo_de_qe
 * @title: NcHICosmoDEQe
 * @short_description: Dark Energy -- quintessence (inspired).
 *
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de_qe.h"
#include "model/nc_hicosmo_de_linder.h"

G_DEFINE_TYPE (NcHICosmoDEQe, nc_hicosmo_de_qe, NC_TYPE_HICOSMO_DE);

#define VECTOR  (NCM_MODEL (cosmo_de)->params)
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define OMEGA_0 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_QE_W0))
#define OMEGA_1 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_QE_W1))

static gdouble
_nc_hicosmo_de_qe_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  gdouble x = 1.0 + z;
  gdouble lnx = log1p (z);
  return OMEGA_X * pow (x, 3.0 * (1.0 + OMEGA_0 / (1.0 + OMEGA_1 * lnx)));
}

static gdouble
_nc_hicosmo_de_qe_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble lnx = log1p (z);
  const gdouble E2Omega_de = OMEGA_X * pow (x, 3.0 * (1.0 + OMEGA_0 / (1.0 + OMEGA_1 * lnx)));

  return 3.0 * ((1.0 + OMEGA_0 / gsl_pow_2 (1.0 + OMEGA_1 * lnx)) / x) * E2Omega_de;
}

/**
 * nc_hicosmo_de_qe_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoDEQe *
nc_hicosmo_de_qe_new (void)
{
  NcHICosmoDEQe *qe = g_object_new (NC_TYPE_HICOSMO_DE_QE, NULL);
  return qe;
}

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_qe_init (NcHICosmoDEQe *qe)
{
  NCM_UNUSED (qe);
}

static void
nc_hicosmo_de_qe_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_qe_parent_class)->finalize (object);
}

static void
nc_hicosmo_de_qe_class_init (NcHICosmoDEQeClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoDEClass* parent_class = NC_HICOSMO_DE_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize     = &nc_hicosmo_de_qe_finalize;

  nc_hicosmo_de_set_E2Omega_de_impl (parent_class, &_nc_hicosmo_de_qe_E2Omega_de);
  nc_hicosmo_de_set_dE2Omega_de_dz_impl (parent_class, &_nc_hicosmo_de_qe_dE2Omega_de_dz);

  ncm_model_class_set_name_nick (model_class, "Quintessence-inspired parametrization", "Quintessence");
  ncm_model_class_add_params (model_class, 2, 0, PROP_SIZE);
  /* Set w_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_QE_W0, "w_0", "w0",
                               -10.0, 1.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_LINDER_DEFAULT_W0,
                               NCM_PARAM_TYPE_FREE);
  /* Set w_1 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_QE_W1, "w_1", "w1",
                               -5.0, 5.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_LINDER_DEFAULT_W1,
                               NCM_PARAM_TYPE_FREE);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}
