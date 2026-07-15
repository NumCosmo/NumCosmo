/***************************************************************************
 *            nc_galaxy_shape_pop_gauss_private.h
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_gauss_private.h
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

#ifndef _NC_GALAXY_SHAPE_POP_GAUSS_PRIVATE_H_
#define _NC_GALAXY_SHAPE_POP_GAUSS_PRIVATE_H_

#include <glib.h>
#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss.h"
#include "ncm/algebra/ncm_laurent_series.h"

G_BEGIN_DECLS

/*
 * Shared with sibling models resolving the same truncated-Gaussian family
 * differently (e.g. GaussLocal, which sources sigma from a per-galaxy input
 * instead of a model parameter). NcGalaxyShapePopGauss and its siblings are
 * NOT related by inheritance: they share these plain functions directly
 * (reused as class vfunc pointers or called outright), not a common
 * instantiable base. Not part of the public API.
 */
void _nc_galaxy_shape_pop_gauss_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

gdouble _nc_galaxy_shape_pop_gauss_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x);
void _nc_galaxy_shape_pop_gauss_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
void _nc_galaxy_shape_pop_gauss_ldata_set_sigma (NcGalaxyShapePopData *data, const gdouble sigma);
gdouble _nc_galaxy_shape_pop_gauss_ldata_get_sigma (NcGalaxyShapePopData *data);
void _nc_galaxy_shape_pop_gauss_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                                      const NcmLaurentSeriesTPS *rho2_series, NcmLaurentSeriesTPS *out);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_GAUSS_PRIVATE_H_ */

