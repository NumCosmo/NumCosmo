/***************************************************************************
 *            nc_halo_density_profile_hernquist.h
 *
 *  Sat April 11 17:23:06 2020
 *  Copyright  2020
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_hernquist.h
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_HALO_DENSITY_PROFILE_HERNQUIST_H_
#define _NC_HALO_DENSITY_PROFILE_HERNQUIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST (nc_halo_density_profile_hernquist_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloDensityProfileHernquist, nc_halo_density_profile_hernquist, NC, HALO_DENSITY_PROFILE_HERNQUIST, NcHaloDensityProfile)

NcHaloDensityProfileHernquist *nc_halo_density_profile_hernquist_new (NcHaloMassSummary * hms);

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_HERNQUIST_H_ */

