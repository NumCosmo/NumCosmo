/***************************************************************************
 *            nc_galaxy_wl_binned_scatter_gauss.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_scatter_gauss.h
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_GALAXY_WL_BINNED_SCATTER_GAUSS_H_
#define _NC_GALAXY_WL_BINNED_SCATTER_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_binned_scatter.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_BINNED_SCATTER_GAUSS (nc_galaxy_wl_binned_scatter_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLBinnedScatterGauss, nc_galaxy_wl_binned_scatter_gauss, NC, GALAXY_WL_BINNED_SCATTER_GAUSS, NcGalaxyWLBinnedScatter)

NcGalaxyWLBinnedScatterGauss *nc_galaxy_wl_binned_scatter_gauss_new (void);
NcGalaxyWLBinnedScatterGauss *nc_galaxy_wl_binned_scatter_gauss_ref (NcGalaxyWLBinnedScatterGauss *sg);
void nc_galaxy_wl_binned_scatter_gauss_free (NcGalaxyWLBinnedScatterGauss *sg);
void nc_galaxy_wl_binned_scatter_gauss_clear (NcGalaxyWLBinnedScatterGauss **sg);

void nc_galaxy_wl_binned_scatter_gauss_set_include_lnnorm (NcGalaxyWLBinnedScatterGauss *sg, gboolean include_lnnorm);
gboolean nc_galaxy_wl_binned_scatter_gauss_get_include_lnnorm (NcGalaxyWLBinnedScatterGauss *sg);

G_END_DECLS

#endif /* _NC_GALAXY_WL_BINNED_SCATTER_GAUSS_H_ */

