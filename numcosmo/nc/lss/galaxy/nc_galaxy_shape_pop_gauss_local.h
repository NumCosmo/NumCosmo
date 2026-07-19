/***************************************************************************
 *            nc_galaxy_shape_pop_gauss_local.h
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_gauss_local.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_H_
#define _NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_POP_GAUSS_LOCAL (nc_galaxy_shape_pop_gauss_local_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapePopGaussLocal, nc_galaxy_shape_pop_gauss_local, NC, GALAXY_SHAPE_POP_GAUSS_LOCAL, NcGalaxyShapePop)

NcGalaxyShapePopGaussLocal *nc_galaxy_shape_pop_gauss_local_new (void);
NcGalaxyShapePopGaussLocal *nc_galaxy_shape_pop_gauss_local_ref (NcGalaxyShapePopGaussLocal *gspgl);

void nc_galaxy_shape_pop_gauss_local_free (NcGalaxyShapePopGaussLocal *gspgl);
void nc_galaxy_shape_pop_gauss_local_clear (NcGalaxyShapePopGaussLocal **gspgl);

void nc_galaxy_shape_pop_gauss_local_data_set (NcGalaxyShapePopGaussLocal *gspgl, NcGalaxyShapePopData *data, const gdouble e_rms);
gdouble nc_galaxy_shape_pop_gauss_local_data_get (NcGalaxyShapePopGaussLocal *gspgl, NcGalaxyShapePopData *data);

#define NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_COL_E_RMS "e_rms"

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_H_ */

