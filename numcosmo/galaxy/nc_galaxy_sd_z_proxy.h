/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_z_proxy.h
 *
 *  Sat May 20 23:09:07 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_z_proxy.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_GALAXY_SD_Z_PROXY_H_
#define _NC_GALAXY_SD_Z_PROXY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_Z_PROXY (nc_galaxy_sd_z_proxy_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDZProxy, nc_galaxy_sd_z_proxy, NC, GALAXY_SD_Z_PROXY, NcmModel)

struct _NcGalaxySDZProxyClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gboolean (*gen) (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z, gdouble *gen_zp);
  gdouble (*integ) (NcGalaxySDZProxy *gsdzp, const gdouble z, const gdouble zp);
  void (*get_true_z_lim) (NcGalaxySDZProxy *gsdzp, const gdouble zp, gdouble *z_min, gdouble *z_max);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_z_proxy);

NcGalaxySDZProxy *nc_galaxy_sd_z_proxy_ref (NcGalaxySDZProxy *gsdzp);

void nc_galaxy_sd_z_proxy_free (NcGalaxySDZProxy *gsdzp);
void nc_galaxy_sd_z_proxy_clear (NcGalaxySDZProxy **gsdzp);

gboolean nc_galaxy_sd_z_proxy_gen (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z, gdouble *gen_zp);
gdouble nc_galaxy_sd_z_proxy_integ (NcGalaxySDZProxy *gsdzp, const gdouble z, const gdouble zp);
void nc_galaxy_sd_z_proxy_get_true_z_lim (NcGalaxySDZProxy *gsdzp, const gdouble zp, gdouble *z_min, gdouble *z_max);

G_END_DECLS

#endif /* _NC_GALAXY_SD_Z_PROXY_H_ */

