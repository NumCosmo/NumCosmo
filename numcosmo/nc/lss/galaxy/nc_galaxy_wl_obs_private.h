/***************************************************************************
 *            nc_galaxy_wl_obs_private.h
 *
 *  Sun Nov 15 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_private.h
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

#ifndef _NC_GALAXY_WL_OBS_PRIVATE_H_
#define _NC_GALAXY_WL_OBS_PRIVATE_H_

#include <glib.h>
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"

G_BEGIN_DECLS

/*
 * The P(z) spline of row i as the catalog holds it, possibly not prepared: for the
 * likelihood, which reads only its knots and values. Not part of the public API, whose
 * nc_galaxy_wl_obs_peek_pz() prepares the spline it returns.
 */
NcmSpline *_nc_galaxy_wl_obs_peek_pz_unprepared (NcGalaxyWLObs *obs, const guint i);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_PRIVATE_H_ */

