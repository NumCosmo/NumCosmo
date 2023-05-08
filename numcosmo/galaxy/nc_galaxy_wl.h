/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl.c
 *
 *  Mon May 08 18:13:24 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl.c
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

#ifndef _NC_GALAXY_WL_H_
#define _NC_GALAXY_WL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/galaxy/nc_gsd_shape.h>
#include <numcosmo/galaxy/nc_gsd_z_proxy.h>
#include <numcosmo/galaxy/nc_gsd_position.h>


G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL             (nc_galaxy_wl_get_type ())
#define NC_GALAXY_WL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL, NcGalaxyWL))
#define NC_GALAXY_WL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL, NcGalaxyWLClass))
#define NC_IS_GALAXY_WL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL))
#define NC_IS_GALAXY_WL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL))
#define NC_GALAXY_WL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL, NcGalaxyWLClass))

typedef struct _NcGalaxyWLClass NcGalaxyWLClass;
typedef struct _NcGalaxyWL NcGalaxyWL;
typedef struct _NcGalaxyWLPrivate NcGalaxyWLPrivate;

struct _NcGalaxyWLClass
{
  /*< private >*/
  NcGalaxyWLDistClass parent_class;
};

struct _NcGalaxyWL
{
  /*< private >*/
  NcGalaxyWLDist parent_instance;
  NcGalaxyWLPrivate *priv;
};

GType nc_galaxy_wl_get_type (void) G_GNUC_CONST;

NcGalaxyWL *nc_galaxy_wl_new (NcGSDShape *s_dist, NcGSDZProxy *zp_dist, NcGSDPosition *rz_dist);
NcGalaxyWL *nc_galaxy_wl_ref (NcGalaxyWL *gwl);

void nc_galaxy_wl_free (NcGalaxyWL *gwl);
void nc_galaxy_wl_clear (NcGalaxyWL **gwl);

G_END_DECLS

#endif /* _NC_GALAXY_WL_H_ */

