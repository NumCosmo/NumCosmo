/***************************************************************************
 *            nc_galaxy_shape_factor_private.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_private.h
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

/*
 * The internal interface of NcGalaxyShapeFactor that NcDataClusterWLFactor and the
 * tests use and the bindings do not: the evaluation of a galaxy at many shears at once,
 * the line of the shears of a galaxy (its redshift calibration) and the critical
 * densities a sampler's copies share. The functions are documented in
 * nc_galaxy_shape_factor.c.
 */

#ifndef _NC_GALAXY_SHAPE_FACTOR_PRIVATE_H_
#define _NC_GALAXY_SHAPE_FACTOR_PRIVATE_H_

#ifndef NUMCOSMO_GIR_SCAN

#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

typedef struct _NcGalaxyShapeFactorLine NcGalaxyShapeFactorLine;

void nc_galaxy_shape_factor_eval_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);
void nc_galaxy_shape_factor_eval_line (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *t, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);

NcGalaxyShapeFactorLine *nc_galaxy_shape_factor_line_new (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const gdouble z_lo, const gdouble z_hi);
gdouble nc_galaxy_shape_factor_line_eval (NcGalaxyShapeFactorLine *line, const gdouble z);
gboolean nc_galaxy_shape_factor_line_has_table (NcGalaxyShapeFactorLine *line);
void nc_galaxy_shape_factor_line_free (NcGalaxyShapeFactorLine *line);

GBytes *nc_galaxy_shape_factor_data_peek_crit (NcGalaxyShapeFactorData *data);
void nc_galaxy_shape_factor_data_share_crit (NcGalaxyShapeFactorData *data, GBytes *crit, const guint n_nodes);

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _NC_GALAXY_SHAPE_FACTOR_PRIVATE_H_ */

