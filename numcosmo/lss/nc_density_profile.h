/***************************************************************************
 *            nc_density_profile.h
 *
 *  Sat June 07 19:45:55 2014
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile.h
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

#ifndef _NC_DENSITY_PROFILE_H_
#define _NC_DENSITY_PROFILE_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_DENSITY_PROFILE             (nc_density_profile_get_type ())
#define NC_DENSITY_PROFILE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DENSITY_PROFILE, NcDensityProfile))
#define NC_DENSITY_PROFILE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DENSITY_PROFILE, NcDensityProfileClass))
#define NC_IS_DENSITY_PROFILE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DENSITY_PROFILE))
#define NC_IS_DENSITY_PROFILE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DENSITY_PROFILE))
#define NC_DENSITY_PROFILE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DENSITY_PROFILE, NcDensityProfileClass))

typedef struct _NcDensityProfileClass NcDensityProfileClass;
typedef struct _NcDensityProfile NcDensityProfile;

struct _NcDensityProfileClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval_fourier) (NcDensityProfile *dp, NcHICosmo *model, const gdouble k, const gdouble M, const gdouble z);
};

struct _NcDensityProfile
{
  /*< private >*/
  GObject parent_instance;

};

GType nc_density_profile_get_type (void) G_GNUC_CONST;

NcDensityProfile *nc_density_profile_new_from_name (gchar *density_profile_name);
gdouble nc_density_profile_eval_fourier (NcDensityProfile *dp, NcHICosmo *model, const gdouble k, const gdouble M, const gdouble z); 
void nc_density_profile_free (NcDensityProfile *dp);
void nc_density_profile_clear (NcDensityProfile **dp);

G_END_DECLS

#endif /* _NC_DENSITY_PROFILE_H_ */

