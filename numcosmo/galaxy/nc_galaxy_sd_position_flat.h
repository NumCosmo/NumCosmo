/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position_flat.h
 *
 *  Wed March 1 12:53:13 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_flat.h
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

#ifndef _NC_GALAXY_SD_POSITION_FLAT_H_
#define _NC_GALAXY_SD_POSITION_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION_FLAT             (nc_galaxy_sd_position_flat_get_type ())
#define NC_GALAXY_SD_POSITION_FLAT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SD_POSITION_FLAT, NcGalaxySDPositionFlat))
#define NC_GALAXY_SD_POSITION_FLAT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SD_POSITION_FLAT, NcGalaxySDPositionFlatClass))
#define NC_IS_GALAXY_SD_POSITION_FLAT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SD_POSITION_FLAT))
#define NC_IS_GALAXY_SD_POSITION_FLAT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SD_POSITION_FLAT))
#define NC_GALAXY_SD_POSITION_FLAT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SD_POSITION_FLAT, NcGalaxySDPositionFlatClass))

typedef struct _NcGalaxySDPositionFlatClass NcGalaxySDPositionFlatClass;
typedef struct _NcGalaxySDPositionFlat NcGalaxySDPositionFlat;
typedef struct _NcGalaxySDPositionFlatPrivate NcGalaxySDPositionFlatPrivate;

struct _NcGalaxySDPositionFlatClass
{
  /*< private >*/
  NcGalaxySDPositionClass parent_class;
};

struct _NcGalaxySDPositionFlat
{
  /*< private >*/
  NcGalaxySDPosition parent_instance;
  NcGalaxySDPositionFlatPrivate *priv;
};

GType nc_galaxy_sd_position_flat_get_type (void) G_GNUC_CONST;

NcGalaxySDPositionFlat *nc_galaxy_sd_position_flat_new ();
NcGalaxySDPositionFlat *nc_galaxy_sd_position_flat_ref (NcGalaxySDPositionFlat *gsdpflat);

void nc_galaxy_sd_position_flat_free (NcGalaxySDPositionFlat *gsdpflat);
void nc_galaxy_sd_position_flat_clear (NcGalaxySDPositionFlat **gsdpflat);

void nc_galaxy_sd_position_flat_set_z_lim (NcGalaxySDPositionFlat *gsdpflat, NcmVector *lim);
NcmVector *nc_galaxy_sd_position_flat_peek_z_lim (NcGalaxySDPositionFlat *gsdpflat);
void nc_galaxy_sd_position_flat_set_r_lim (NcGalaxySDPositionFlat *gsdpflat, NcmVector *lim);
NcmVector *nc_galaxy_sd_position_flat_peek_r_lim (NcGalaxySDPositionFlat *gsdpflat);

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_FLAT_H_ */

