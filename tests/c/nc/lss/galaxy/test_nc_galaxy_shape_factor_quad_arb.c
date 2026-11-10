/***************************************************************************
 *            test_nc_galaxy_shape_factor_quad_arb.c
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

/* Acceptance tier of the NcGalaxyShapeFactorQuad tests, see test_nc_galaxy_shape_factor_quad_common.h. */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */

#include "test_nc_galaxy_shape_factor_quad_common.h"

gint
main (gint argc, gchar *argv[])
{
  return test_nc_galaxy_shape_factor_quad_main (argc, argv, TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_ACCEPTANCE);
}

