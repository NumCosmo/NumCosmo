/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position.h
 *
 *  Sat May 20 18:41:47 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position.h
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
#ifndef _NC_GALAXY_SD_POSITION_H_
#define _NC_GALAXY_SD_POSITION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>


G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION             (nc_galaxy_sd_position_get_type ())
#define NC_GALAXY_SD_POSITION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SD_POSITION, NcGalaxySDPosition))
#define NC_GALAXY_SD_POSITION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SD_POSITION, NcGalaxySDPositionClass))
#define NC_IS_GALAXY_SD_POSITION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SD_POSITION))
#define NC_IS_GALAXY_SD_POSITION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SD_POSITION))
#define NC_GALAXY_SD_POSITION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SD_POSITION, NcGalaxySDPositionClass))

typedef struct _NcGalaxySDPositionClass NcGalaxySDPositionClass;
typedef struct _NcGalaxySDPosition NcGalaxySDPosition;
typedef struct _NcGalaxySDPositionPrivate NcGalaxySDPositionPrivate;

struct _NcGalaxySDPositionClass
{
  /*< private >*/
  NcGalaxySDPositionDistClass parent_class;
};

struct _NcGalaxySDPosition
{
  /*< private >*/
  NcGalaxySDPositionDist parent_instance;
  NcGalaxySDPositionPrivate *priv;
};

GType nc_galaxy_sd_position_get_type (void) G_GNUC_CONST;

NcGalaxySDPosition *nc_galaxy_sd_position_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *rz_dist);
NcGalaxySDPosition *nc_galaxy_sd_position_ref (NcGalaxySDPosition *gsdp);

void nc_galaxy_sd_position_free (NcGalaxySDPosition *gsdp);
void nc_galaxy_sd_position_clear (NcGalaxySDPosition **gsdp);

NCM_INLINE NcmVector nc_galaxy_sd_position_gen (NcGalaxySDPosition *gsdp, NcmRNG *rng);
NCM_INLINE gdouble nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, NcmVector *pos);

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_H_ */

#ifndef _NC_GALAXY_SD_POSITION_INLINE_H_
#define _NC_GALAXY_SD_POSITION_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
nc_galaxy_sd_position_gen (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionClass *gsdp_class = NC_GALAXY_SD_POSITION_GET_CLASS (gsdp);

  if (gsdp_class->gen != NULL)
    NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen (gsdp, rng);
}

NCM_INLINE void
nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, NcmVector *pos)
{
  NcGalaxySDPositionClass *gsdp_class = NC_GALAXY_SD_POSITION_GET_CLASS (gsdp);

  if (gsdp_class->integ != NULL)
    NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->integ (gsdp, pos);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_GALAXY_SD_POSITION_INLINE_H_ */

