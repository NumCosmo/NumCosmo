/***************************************************************************
 *            test_nc_galaxy_shape_factor_quad_common.h
 *
 *  Tue Nov 10 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

/*
 * The NcGalaxyShapeFactorQuad tests compare ln P of the default path and of
 * NcGalaxyShapeFactorQuad:use-lp with values certified by the Arb reference, in two modes
 * chosen by each executable's main(): %TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_UNIT checks the rows
 * of the table marked unit, plus the |g| = 1 closed form and the properties;
 * %TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_ACCEPTANCE checks every row.
 */

#ifndef _TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_COMMON_H
#define _TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_COMMON_H

#include <glib.h>
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/**
 * TestNcGalaxyShapeFactorQuadMode:
 * @TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_UNIT: the rows marked unit, the |g| = 1 and property tests
 * @TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_ACCEPTANCE: every row of the certified table
 */
typedef enum _TestNcGalaxyShapeFactorQuadMode
{
  TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_UNIT,
  TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_ACCEPTANCE,
} TestNcGalaxyShapeFactorQuadMode;

gint test_nc_galaxy_shape_factor_quad_main (gint argc, gchar *argv[], TestNcGalaxyShapeFactorQuadMode mode);

G_END_DECLS

#endif /* _TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_COMMON_H */

