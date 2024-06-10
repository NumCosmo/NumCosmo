/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_z_proxy_dirac.h
 *
 *  Wed June 21 19:40:56 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_z_proxy_dirac.h
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

#ifndef _NC_GALAXY_SD_Z_PROXY_DIRAC_H_
#define _NC_GALAXY_SD_Z_PROXY_DIRAC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_sd_z_proxy.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC             (nc_galaxy_sd_z_proxy_dirac_get_type ())
#define NC_GALAXY_SD_Z_PROXY_DIRAC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC, NcGalaxySDZProxyDirac))
#define NC_GALAXY_SD_Z_PROXY_DIRAC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC, NcGalaxySDZProxyDiracClass))
#define NC_IS_GALAXY_SD_Z_PROXY_DIRAC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC))
#define NC_IS_GALAXY_SD_Z_PROXY_DIRAC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC))
#define NC_GALAXY_SD_Z_PROXY_DIRAC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC, NcGalaxySDZProxyDiracClass))

typedef struct _NcGalaxySDZProxyDiracClass NcGalaxySDZProxyDiracClass;
typedef struct _NcGalaxySDZProxyDirac NcGalaxySDZProxyDirac;
typedef struct _NcGalaxySDZProxyDiracPrivate NcGalaxySDZProxyDiracPrivate;

struct _NcGalaxySDZProxyDiracClass
{
  /*< private >*/
  NcGalaxySDZProxyClass parent_class;
};

struct _NcGalaxySDZProxyDirac
{
  /*< private >*/
  NcGalaxySDZProxy parent_instance;
  NcGalaxySDZProxyDiracPrivate *priv;
};

GType nc_galaxy_sd_z_proxy_dirac_get_type (void) G_GNUC_CONST;

NcGalaxySDZProxyDirac *nc_galaxy_sd_z_proxy_dirac_new ();
NcGalaxySDZProxyDirac *nc_galaxy_sd_z_proxy_dirac_ref (NcGalaxySDZProxyDirac *gsdzpdirac);

void nc_galaxy_sd_z_proxy_dirac_free (NcGalaxySDZProxyDirac *gsdzpdirac);
void nc_galaxy_sd_z_proxy_dirac_clear (NcGalaxySDZProxyDirac **gsdzpdirac);

G_END_DECLS

#endif /* _NC_GALAXY_SD_Z_PROXY_DIRAC_H_ */

