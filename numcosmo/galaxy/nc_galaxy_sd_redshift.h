/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_redshift.h
 *
 *  Wed Jul 31 21:09:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_redshift.h
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
#ifndef _NC_GALAXY_SD_REDSHIFT_H_
#define _NC_GALAXY_SD_REDSHIFT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_REDSHIFT (nc_galaxy_sd_redshift_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDRedshift, nc_galaxy_sd_redshift, NC, GALAXY_SD_REDSHIFT, NcmModel)

struct _NcGalaxySDRedshiftClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*gen) (NcGalaxySDRedshift *gsdr, NcmRNG *rng);
  gdouble (*integ) (NcGalaxySDRedshift *gsdr, NcmVector *data);
  gboolean (*set_lim) (NcGalaxySDRedshift *gsdr, const gdouble z_min, const gdouble z_max);
  gboolean (*get_lim) (NcGalaxySDRedshift *gsdr, gdouble *z_min, gdouble *z_max);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_redshift);

NcGalaxySDRedshift *nc_galaxy_sd_redshift_ref (NcGalaxySDRedshift *gsdr);

void nc_galaxy_sd_redshift_free (NcGalaxySDRedshift *gsdr);
void nc_galaxy_sd_redshift_clear (NcGalaxySDRedshift **gsdr);

gboolean nc_galaxy_sd_redshift_set_lim (NcGalaxySDRedshift *gsdr, const gdouble z_min, const gdouble z_max);
gboolean nc_galaxy_sd_redshift_get_lim (NcGalaxySDRedshift *gsdr, gdouble *z_min, gdouble *z_max);

gdouble nc_galaxy_sd_redshift_gen (NcGalaxySDRedshift *gsdr, NcmRNG *rng);
gdouble nc_galaxy_sd_redshift_integ (NcGalaxySDRedshift *gsdr, NcmVector *data);

G_END_DECLS

#endif /* _NC_GALAXY_SD_REDSHIFT_H_ */

