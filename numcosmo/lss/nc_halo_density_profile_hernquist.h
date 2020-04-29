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

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST             (nc_halo_density_profile_hernquist_get_type ())
#define NC_HALO_DENSITY_PROFILE_HERNQUIST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST, NcHaloDensityProfileHernquist))
#define NC_HALO_DENSITY_PROFILE_HERNQUIST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST, NcHaloDensityProfileHernquistClass))
#define NC_IS_HALO_DENSITY_PROFILE_HERNQUIST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST))
#define NC_IS_HALO_DENSITY_PROFILE_HERNQUIST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST))
#define NC_HALO_DENSITY_PROFILE_HERNQUIST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST, NcHaloDensityProfileHernquistClass))

typedef struct _NcHaloDensityProfileHernquistClass NcHaloDensityProfileHernquistClass;
typedef struct _NcHaloDensityProfileHernquist NcHaloDensityProfileHernquist;

struct _NcHaloDensityProfileHernquistClass
{
  /*< private > */
  NcHaloDensityProfileClass parent_class;
};

struct _NcHaloDensityProfileHernquist
{
  /*< private > */
  NcHaloDensityProfile parent_instance;
};

GType nc_halo_density_profile_hernquist_get_type (void) G_GNUC_CONST;

NcHaloDensityProfileHernquist *nc_halo_density_profile_hernquist_new (const NcHaloDensityProfileMassDef mdef, const gdouble Delta);

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_HERNQUIST_H_ */

