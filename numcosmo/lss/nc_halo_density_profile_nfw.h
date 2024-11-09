/***************************************************************************
 *            nc_halo_density_profile_nfw.h
 *
 *  Tue June 10 16:40:06 2014
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_nfw.h
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_HALO_DENSITY_PROFILE_NFW_H_
#define _NC_HALO_DENSITY_PROFILE_NFW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_NFW (nc_halo_density_profile_nfw_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloDensityProfileNFW, nc_halo_density_profile_nfw, NC, HALO_DENSITY_PROFILE_NFW, NcHaloDensityProfile)

void nc_halo_density_profile_nfw_class_set_ni (gboolean num);
NcHaloDensityProfileNFW *nc_halo_density_profile_nfw_new (NcHaloMassSummary *hms);

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_NFW_H_ */

