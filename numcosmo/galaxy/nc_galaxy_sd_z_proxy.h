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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_GALAXY_SD_Z_PROXY_H_
#define _NC_GALAXY_SD_Z_PROXY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>


G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION             (nc_galaxy_sd_z_proxy_get_type ())
#define NC_GALAXY_SD_Z_PROXY(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SD_POSITION, NcGalaxySDZProxy))
#define NC_GALAXY_SD_Z_PROXY_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SD_POSITION, NcGalaxySDZProxyClass))
#define NC_IS_GALAXY_SD_POSITION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SD_POSITION))
#define NC_IS_GALAXY_SD_POSITION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SD_POSITION))
#define NC_GALAXY_SD_Z_PROXY_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SD_POSITION, NcGalaxySDZProxyClass))

typedef struct _NcGalaxySDZProxyClass NcGalaxySDZProxyClass;
typedef struct _NcGalaxySDZProxy NcGalaxySDZProxy;
typedef struct _NcGalaxySDZProxyPrivate NcGalaxySDZProxyPrivate;

struct _NcGalaxySDZProxyClass
{
  /*< private >*/
  NcGalaxySDZProxyDistClass parent_class;
};

struct _NcGalaxySDZProxy
{
  /*< private >*/
  NcGalaxySDZProxyDist parent_instance;
  NcGalaxySDZProxyPrivate *priv;
};

GType nc_galaxy_sd_z_proxy_get_type (void) G_GNUC_CONST;

NcGalaxySDZProxy *nc_galaxy_sd_z_proxy_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDZProxy *rz_dist);
NcGalaxySDZProxy *nc_galaxy_sd_z_proxy_ref (NcGalaxySDZProxy *gsdzp);

void nc_galaxy_sd_z_proxy_free (NcGalaxySDZProxy *gsdzp);
void nc_galaxy_sd_z_proxy_clear (NcGalaxySDZProxy **gsdzp);

NCM_INLINE gdouble nc_galaxy_sd_z_proxy_gen (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z);
NCM_INLINE gdouble nc_galaxy_sd_z_proxy_integ (NcGalaxySDZProxy *gsdzp, const gdouble z);

G_END_DECLS

#endif /* _NC_GALAXY_SD_Z_PROXY_H_ */

#ifndef _NC_GALAXY_SD_Z_PROXY_INLINE_H_
#define _NC_GALAXY_SD_Z_PROXY_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
nc_galaxy_sd_z_proxy_gen (NcGalaxySDZProxy *gsdzp, const gdouble z)
{
  NcGalaxySDZProxyClass *gsdzp_class = NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp);

  if (gsdzp_class->gen != NULL)
    NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp)->gen (gsdzp, rng, z);
}

NCM_INLINE void
nc_galaxy_sd_z_proxy_integ (NcGalaxySDZProxy *gsdzp, const gdouble z)
{
  NcGalaxySDZProxyClass *gsdzp_class = NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp);

  if (gsdzp_class->integ != NULL)
    NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp)->integ (gsdzp, z);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_GALAXY_SD_Z_PROXY_INLINE_H_ */

