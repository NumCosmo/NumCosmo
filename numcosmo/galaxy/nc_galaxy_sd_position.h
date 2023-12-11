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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_GALAXY_SD_POSITION_H_
#define _NC_GALAXY_SD_POSITION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION (nc_galaxy_sd_position_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDPosition, nc_galaxy_sd_position, NC, GALAXY_SD_POSITION, NcmModel)

struct _NcGalaxySDPositionClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*gen_r) (NcGalaxySDPosition *gsdp, NcmRNG *rng);
  gdouble (*gen_z) (NcGalaxySDPosition *gsdp, NcmRNG *rng);
  gdouble (*integ) (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z);
  void (*set_z_lim) (NcGalaxySDPosition *gsdp, const gdouble z_min, const gdouble z_max);
  void (*set_r_lim) (NcGalaxySDPosition *gsdp, const gdouble r_min, const gdouble r_max);
  void (*get_z_lim) (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max);
  void (*get_r_lim) (NcGalaxySDPosition *gsdp, gdouble *r_min, gdouble *r_max);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[11];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_position);

NcGalaxySDPosition *nc_galaxy_sd_position_ref (NcGalaxySDPosition *gsdp);

void nc_galaxy_sd_position_free (NcGalaxySDPosition *gsdp);
void nc_galaxy_sd_position_clear (NcGalaxySDPosition **gsdp);

void nc_galaxy_sd_position_set_z_lim (NcGalaxySDPosition *gsdp, const gdouble z_min, const gdouble z_max);
void nc_galaxy_sd_position_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max);
void nc_galaxy_sd_position_set_r_lim (NcGalaxySDPosition *gsdp, const gdouble r_min, const gdouble r_max);
void nc_galaxy_sd_position_get_r_lim (NcGalaxySDPosition *gsdp, gdouble *r_min, gdouble *r_max);

gdouble nc_galaxy_sd_position_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng);
gdouble nc_galaxy_sd_position_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng);
gdouble nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z);

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_H_ */

