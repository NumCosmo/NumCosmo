/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_obs_redshift.h
 *
 *  Thu Aug 1 00:48:12 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
#ifndef _NC_GALAXY_SD_OBS_REDSHIFT_H_
#define _NC_GALAXY_SD_OBS_REDSHIFT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDObsRedshift, nc_galaxy_sd_obs_redshift, NC, GALAXY_SD_OBS_REDSHIFT, GObject)

struct _NcGalaxySDObsRedshiftClass
{
  /*< private >*/
  GObjectClass parent_class;

  gdouble (*gen) (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data);
  gdouble (*integ) (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data);
  gboolean (*get_header) (NcGalaxySDObsRedshift *gsdor);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

NcGalaxySDObsRedshift *nc_galaxy_sd_obs_redshift_ref (NcGalaxySDObsRedshift *gsdor);

void nc_galaxy_sd_obs_redshift_free (NcGalaxySDObsRedshift *gsdor);
void nc_galaxy_sd_obs_redshift_clear (NcGalaxySDObsRedshift **gsdor);

gdouble nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data);
gdouble nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data);
gboolean nc_galaxy_sd_obs_redshift_get_header (NcGalaxySDObsRedshift *gsdor);

G_END_DECLS

#endif /* _NC_GALAXY_SD_OBS_REDSHIFT_H_ */

