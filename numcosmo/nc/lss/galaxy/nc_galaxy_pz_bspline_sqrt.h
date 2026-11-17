/***************************************************************************
 *            nc_galaxy_pz_bspline_sqrt.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_pz_bspline_sqrt.h
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_GALAXY_PZ_BSPLINE_SQRT_H_
#define _NC_GALAXY_PZ_BSPLINE_SQRT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/spline/ncm_bspline.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_pz.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_PZ_BSPLINE_SQRT (nc_galaxy_pz_bspline_sqrt_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyPzBSplineSqrt, nc_galaxy_pz_bspline_sqrt, NC, GALAXY_PZ_BSPLINE_SQRT, NcGalaxyPz)

NcGalaxyPzBSplineSqrt *nc_galaxy_pz_bspline_sqrt_new (NcmBSpline * s);
NcGalaxyPzBSplineSqrt *nc_galaxy_pz_bspline_sqrt_new_from_table (NcmVector *zv, NcmVector *pv);

NcmBSpline *nc_galaxy_pz_bspline_sqrt_peek_bspline (NcGalaxyPzBSplineSqrt *gpzs);

/**
 * NC_GALAXY_PZ_BSPLINE_SQRT_ORDER:
 *
 * The order of the B-spline of a #NcGalaxyPzBSplineSqrt, degree 5.
 */
#define NC_GALAXY_PZ_BSPLINE_SQRT_ORDER (6)

G_END_DECLS

#endif /* _NC_GALAXY_PZ_BSPLINE_SQRT_H_ */

