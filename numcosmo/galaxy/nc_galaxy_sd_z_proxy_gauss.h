/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_z_proxy_gauss.h
 *
 *  Tue June 1 19:29:53 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_z_proxy_gauss.h
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

#ifndef _NC_GALAXY_SD_Z_PROXY_GAUSS_H_
#define _NC_GALAXY_SD_Z_PROXY_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_sd_z_proxy.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS             (nc_galaxy_sd_z_proxy_gauss_get_type ())
#define NC_GALAXY_SD_Z_PROXY_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS, NcGalaxySDZProxyGauss))
#define NC_GALAXY_SD_Z_PROXY_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS, NcGalaxySDZProxyGaussClass))
#define NC_IS_GALAXY_SD_Z_PROXY_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS))
#define NC_IS_GALAXY_SD_Z_PROXY_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS))
#define NC_GALAXY_SD_Z_PROXY_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS, NcGalaxySDZProxyGaussClass))

typedef struct _NcGalaxySDZProxyGaussClass NcGalaxySDZProxyGaussClass;
typedef struct _NcGalaxySDZProxyGauss NcGalaxySDZProxyGauss;
typedef struct _NcGalaxySDZProxyGaussPrivate NcGalaxySDZProxyGaussPrivate;

struct _NcGalaxySDZProxyGaussClass
{
  /*< private >*/
  NcGalaxySDZProxyClass parent_class;
};

struct _NcGalaxySDZProxyGauss
{
  /*< private >*/
  NcGalaxySDZProxy parent_instance;
  NcGalaxySDZProxyGaussPrivate *priv;
};

GType nc_galaxy_sd_z_proxy_gauss_get_type (void) G_GNUC_CONST;

NcGalaxySDZProxyGauss *nc_galaxy_sd_z_proxy_gauss_new ();
NcGalaxySDZProxyGauss *nc_galaxy_sd_z_proxy_gauss_ref (NcGalaxySDZProxyGauss *gsdzpgauss);

void nc_galaxy_sd_z_proxy_gauss_free (NcGalaxySDZProxyGauss *gsdzpgauss);
void nc_galaxy_sd_z_proxy_gauss_clear (NcGalaxySDZProxyGauss **gsdzpgauss);

void nc_galaxy_sd_z_proxy_gauss_set_z_lim (NcGalaxySDZProxyGauss *gsdzpgauss, NcmVector *lim);
NcmVector *nc_galaxy_sd_z_proxy_gauss_peek_z_lim (NcGalaxySDZProxyGauss *gsdzpgauss);
void nc_galaxy_sd_z_proxy_gauss_set_sigma (NcGalaxySDZProxyGauss *gsdzpgauss, gdouble sigma);
gdouble nc_galaxy_sd_z_proxy_gauss_get_sigma (NcGalaxySDZProxyGauss *gsdzpgauss);

G_END_DECLS

#endif /* _NC_GALAXY_SD_Z_PROXY_GAUSS_H_ */

