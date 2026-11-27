/***************************************************************************
 *            nc_galaxy_shape_pop_shear_scatter.h
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_shear_scatter.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_H_
#define _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER (nc_galaxy_shape_pop_shear_scatter_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyShapePopShearScatter, nc_galaxy_shape_pop_shear_scatter, NC, GALAXY_SHAPE_POP_SHEAR_SCATTER, NcGalaxyShapePop)

struct _NcGalaxyShapePopShearScatterClass
{
  /*< private >*/
  NcGalaxyShapePopClass parent_class;

  gdouble (*eval_ln_p) (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g);
  gdouble (*gen_delta) (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, NcmRNG *rng, const gdouble var_noise_g, const gdouble eps_to_g);

  /* Padding to allow 16 virtual functions without breaking ABI. */
  gpointer padding[14];
};

NcGalaxyShapePopShearScatter *nc_galaxy_shape_pop_shear_scatter_ref (NcGalaxyShapePopShearScatter *gspss);

void nc_galaxy_shape_pop_shear_scatter_free (NcGalaxyShapePopShearScatter *gspss);
void nc_galaxy_shape_pop_shear_scatter_clear (NcGalaxyShapePopShearScatter **gspss);

void nc_galaxy_shape_pop_shear_scatter_set_include_noise (NcGalaxyShapePopShearScatter *gspss, gboolean include_noise);
gboolean nc_galaxy_shape_pop_shear_scatter_get_include_noise (NcGalaxyShapePopShearScatter *gspss);

gdouble nc_galaxy_shape_pop_shear_scatter_eval_ln_p (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g);
gdouble nc_galaxy_shape_pop_shear_scatter_eval_p (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g);
gdouble nc_galaxy_shape_pop_shear_scatter_gen_delta (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, NcmRNG *rng, const gdouble var_noise_g, const gdouble eps_to_g);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_H_ */

