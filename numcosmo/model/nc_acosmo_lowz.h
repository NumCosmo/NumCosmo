/***************************************************************************
 *            nc_acosmo_lowz.h
 *
 *  Fri Jan 26 15:37:04 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2024 <vitenti@uel.br>
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

#ifndef _NC_ACOSMO_LOWZ_H_
#define _NC_ACOSMO_LOWZ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_ACOSMO_LOWZ (nc_acosmo_lowz_get_type ())

G_DECLARE_FINAL_TYPE (NcACosmoLowz, nc_acosmo_lowz, NC, ACOSMO_LOWZ, NcmModel)

NCM_MSET_MODEL_DECLARE_ID (nc_acosmo_lowz);

/**
 * NcACosmoLowzSParams:
 * @NC_ACOSMO_LOWZ_H0: Hubble constant
 *
 * Kinematic parameters
 */
typedef enum /*< enum,underscore_name=NC_ACOSMO_LOWZ_SPARAMS >*/
{
  NC_ACOSMO_LOWZ_H0 = 0,
  /* < private > */
  NC_ACOSMO_LOWZ_SPARAM_LEN, /*< skip >*/
} NcACosmoLowzSParams;


/**
 * NcACosmoLowzVParams:
 * @NC_ACOSMO_LOWZ_ACCEL: acceleration
 * @NC_ACOSMO_LOWZ_SHEAR: shear
 *
 * Kinematic parameters
 *
 */
typedef enum /*< enum,underscore_name=NC_ACOSMO_LOWZ_VPARAMS >*/
{
  NC_ACOSMO_LOWZ_ACCEL = 0,
  NC_ACOSMO_LOWZ_SHEAR,
  /* < private > */
  NC_ACOSMO_LOWZ_VPARAM_LEN, /*< skip >*/
} NcACosmoLowzVParams;


NcACosmoLowz *nc_acosmo_lowz_new (void);
gdouble nc_acosmo_lowz_distance_modulus (NcACosmoLowz *aclz, const gdouble z, const gdouble theta, const gdouble phi);


G_END_DECLS

#endif /* _NC_ACOSMO_LOWZ_H_ */

