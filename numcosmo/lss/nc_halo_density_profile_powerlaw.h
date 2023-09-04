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

#ifndef _NC_HALO_DENSITY_PROFILE_MEAN_H_
#define _NC_HALO_DENSITY_PROFILE_MEAN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_density_profile.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_MEAN             (nc_halo_density_profile_MEAN_get_type ())
#define NC_HALO_DENSITY_PROFILE_MEAN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_DENSITY_PROFILE_MEAN, NcHaloDensityProfileMean))
#define NC_HALO_DENSITY_PROFILE_MEAN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_DENSITY_PROFILE_MEAN, NcHaloDensityProfileMeanClass))
#define NC_IS_HALO_DENSITY_PROFILE_MEAN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_DENSITY_PROFILE_MEAN))
#define NC_IS_HALO_DENSITY_PROFILE_MEAN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_DENSITY_PROFILE_MEAN))
#define NC_HALO_DENSITY_PROFILE_MEAN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_DENSITY_PROFILE_MEAN, NcHaloDensityProfileMeanClass))

typedef struct _NcHaloDensityProfileMeanClass NcHaloDensityProfileMeanClass;
typedef struct _NcHaloDensityProfileMean NcHaloDensityProfileMean;

struct _NcHaloDensityProfileMeanClass
{
  /*< private > */
  NcHaloDensityProfileClass parent_class;
};

struct _NcHaloDensityProfileMean
{
  /*< private > */
  NcHaloDensityProfile parent_instance;
};

GType nc_halo_density_profile_mean_get_type (void) G_GNUC_CONST;

void nc_halo_density_profile_mean_class_set_ni (gboolean num);
NcHaloDensityProfileMean *nc_halo_density_profile_mean_new (const NcHaloDensityProfileMassDef mdef, const gdouble Delta);

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_MEAN_H_ */

