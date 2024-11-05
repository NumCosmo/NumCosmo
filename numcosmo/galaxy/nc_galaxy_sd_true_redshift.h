/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_true_redshift.h
 *
 *  Wed Jul 31 21:09:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_true_redshift.h
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
#ifndef _NC_GALAXY_SD_TRUE_REDSHIFT_H_
#define _NC_GALAXY_SD_TRUE_REDSHIFT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_stats_dist1d.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDTrueRedshift, nc_galaxy_sd_true_redshift, NC, GALAXY_SD_TRUE_REDSHIFT, NcmModel)

struct _NcGalaxySDTrueRedshiftClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*gen) (NcGalaxySDTrueRedshift *gsdtr, NcmRNG *rng);
  gdouble (*integ) (NcGalaxySDTrueRedshift *gsdtr, gdouble z);
  gboolean (*set_lim) (NcGalaxySDTrueRedshift *gsdtr, const gdouble z_min, const gdouble z_max);
  gboolean (*get_lim) (NcGalaxySDTrueRedshift *gsdtr, gdouble *z_min, gdouble *z_max);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_true_redshift);

NcGalaxySDTrueRedshift *nc_galaxy_sd_true_redshift_ref (NcGalaxySDTrueRedshift *gsdtr);

void nc_galaxy_sd_true_redshift_free (NcGalaxySDTrueRedshift *gsdtr);
void nc_galaxy_sd_true_redshift_clear (NcGalaxySDTrueRedshift **gsdtr);

gboolean nc_galaxy_sd_true_redshift_set_lim (NcGalaxySDTrueRedshift *gsdtr, const gdouble z_min, const gdouble z_max);
gboolean nc_galaxy_sd_true_redshift_get_lim (NcGalaxySDTrueRedshift *gsdtr, gdouble *z_min, gdouble *z_max);

gdouble nc_galaxy_sd_true_redshift_gen (NcGalaxySDTrueRedshift *gsdtr, NcmRNG *rng);
gdouble nc_galaxy_sd_true_redshift_integ (NcGalaxySDTrueRedshift *gsdtr, gdouble z);

NcmStatsDist1d *nc_galaxy_sd_true_redshift_dist (NcGalaxySDTrueRedshift *gsdtr, const gdouble reltol, const gdouble abstol);

G_END_DECLS

#endif /* _NC_GALAXY_SD_TRUE_REDSHIFT_H_ */

