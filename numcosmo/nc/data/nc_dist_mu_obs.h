/***************************************************************************
 *            nc_dist_mu_obs.h
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_dist_mu_obs.h
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

#ifndef _NC_DIST_MU_OBS_H_
#define _NC_DIST_MU_OBS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_catalog.h>
#include <numcosmo/ncm/algebra/ncm_quaternion.h>

G_BEGIN_DECLS

#define NC_TYPE_DIST_MU_OBS (nc_dist_mu_obs_get_type ())

G_DECLARE_FINAL_TYPE (NcDistMuObs, nc_dist_mu_obs, NC, DIST_MU_OBS, NcmCatalog)

#define NC_DIST_MU_OBS_COL_RA "ra"
#define NC_DIST_MU_OBS_COL_DEC "dec"
#define NC_DIST_MU_OBS_COL_Z "z"
#define NC_DIST_MU_OBS_COL_MU "mu"
#define NC_DIST_MU_OBS_COL_SIGMA_MU "sigma_mu"

NcDistMuObs * nc_dist_mu_obs_new (guint nrows, GStrv col_names);
NcDistMuObs *nc_dist_mu_obs_ref (NcDistMuObs *obs);

void nc_dist_mu_obs_free (NcDistMuObs *obs);
void nc_dist_mu_obs_clear (NcDistMuObs **obs);

void nc_dist_mu_obs_set_direction (NcDistMuObs *obs, const guint i, const NcmTriVec *n);
void nc_dist_mu_obs_get_direction (NcDistMuObs *obs, const guint i, NcmTriVec *n);

G_END_DECLS

#endif /* _NC_DIST_MU_OBS_H_ */

