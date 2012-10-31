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
 * @title: XCDM Darkenergy Model
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de_xcdm.h"
#include "model/nc_hicosmo_de_linder.h"

G_DEFINE_TYPE (NcHICosmoDEXcdm, nc_hicosmo_de_xcdm, NC_TYPE_HICOSMO_DE);

#define VECTOR  (model->params)
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define W       (ncm_vector_get (VECTOR, NC_HICOSMO_DE_XCDM_W))

static gdouble
_nc_hicosmo_de_xcdm_weff (NcmModel *model, gdouble z)
{
  gdouble x = 1.0 + z;
  return OMEGA_X * pow (x, 3.0 * ( 1.0 + W ) );
}

static gdouble
_nc_hicosmo_de_xcdm_dweff_dz (NcmModel *model, gdouble z)
{
  const gdouble x = 1.0 + z;
  const gdouble weff = OMEGA_X * pow (x, 3.0 * ( 1.0 + W ) );

  return 3.0 * ( 1.0 + W ) / x * weff;
}

/**
 * nc_hicosmo_de_xcdm_new:
 *
 * FIXME
 *
 * Returns: FIXME
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
  object_class->set_property = &ncm_model_class_set_property;
  object_class->get_property = &ncm_model_class_get_property;
  object_class->finalize     = &nc_hicosmo_de_xcdm_finalize;

  nc_hicosmo_de_set_weff_impl (parent_class, &_nc_hicosmo_de_xcdm_weff);
  nc_hicosmo_de_set_dweff_dz_impl (parent_class, &_nc_hicosmo_de_xcdm_dweff_dz);

  ncm_model_class_set_name_nick (model_class, "XCDM - Constant EOS", "XCDM");
  ncm_model_class_add_params (model_class, 1, 0, PROP_SIZE);
  /* Set w_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_XCDM_W, "w", "w",
                               -10.0, 1.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_LINDER_DEFAULT_W0,
                               NCM_PARAM_TYPE_FREE);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}
