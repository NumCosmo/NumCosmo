/***************************************************************************
 *            nc_halo_density_profile_einasto.h
 *
 *  Wed July 17 12:33:27 2019
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_einasto.h
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

#ifndef _NC_HALO_DENSITY_PROFILE_EINASTO_H_
#define _NC_HALO_DENSITY_PROFILE_EINASTO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_EINASTO (nc_halo_density_profile_einasto_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloDensityProfileEinasto, nc_halo_density_profile_einasto, NC, HALO_DENSITY_PROFILE_EINASTO, NcHaloDensityProfile)

/**
 * NcHaloDensityProfileEinastoParams:
 * @NC_HALO_DENSITY_PROFILE_EINASTO_ALPHA: determines how quickly the slope of the inner Einasto profile steepens
 *
 * FIXME
 */
typedef enum _NcHaloDensityProfileEinastoParams
{
  NC_HALO_DENSITY_PROFILE_EINASTO_ALPHA,
  /* < private > */
  NC_HALO_DENSITY_PROFILE_EINASTO_SPARAM_LEN, /*< skip >*/
} NcHaloDensityProfileEinastoParams;

#define NC_HALO_DENSITY_PROFILE_EINASTO_LOCAL_SPARAM_LEN (NC_HALO_DENSITY_PROFILE_EINASTO_SPARAM_LEN - 0)
#define NC_HALO_DENSITY_PROFILE_EINASTO_DEFAULT_ALPHA (0.25)
#define NC_HALO_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL (0.0)

NcHaloDensityProfileEinasto *nc_halo_density_profile_einasto_new (NcHaloMassSummary *hms);

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_EINASTO_H_ */

