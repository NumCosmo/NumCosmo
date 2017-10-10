/***************************************************************************
 *            nc_density_profile_nfw.h
 *
 *  Tue June 10 16:40:06 2014
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile_nfw.h
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

#ifndef _NC_DENSITY_PROFILE_NFW_H_
#define _NC_DENSITY_PROFILE_NFW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_density_profile.h>

G_BEGIN_DECLS

#define NC_TYPE_DENSITY_PROFILE_NFW             (nc_density_profile_nfw_get_type ())
#define NC_DENSITY_PROFILE_NFW(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DENSITY_PROFILE_NFW, NcDensityProfileNFW))
#define NC_DENSITY_PROFILE_NFW_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DENSITY_PROFILE_NFW, NcDensityProfileNFWClass))
#define NC_IS_DENSITY_PROFILE_NFW(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DENSITY_PROFILE_NFW))
#define NC_IS_DENSITY_PROFILE_NFW_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DENSITY_PROFILE_NFW))
#define NC_DENSITY_PROFILE_NFW_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DENSITY_PROFILE_NFW, NcDensityProfileNFWClass))

typedef struct _NcDensityProfileNFWClass NcDensityProfileNFWClass;
typedef struct _NcDensityProfileNFW NcDensityProfileNFW;

/**
 * NcDensityProfileNFWParams:
 * @NC_DENSITY_PROFILE_NFW_C_DELTA: concentration parameter
 * @NC_DENSITY_PROFILE_NFW_M_DELTA: halo mass
 *
 * FIXME
 */
typedef enum _NcDensityProfileNFWParams
{
  NC_DENSITY_PROFILE_NFW_C_DELTA = 0,
  NC_DENSITY_PROFILE_NFW_M_DELTA, 
  /* < private > */
  NC_DENSITY_PROFILE_NFW_SPARAM_LEN, /*< skip >*/
} NcDensityProfileNFWParams;

#define NC_DENSITY_PROFILE_NFW_DEFAULT_C_DELTA  (4.0)
#define NC_DENSITY_PROFILE_NFW_DEFAULT_M_DELTA  (2.0e14)

#define NC_DENSITY_PROFILE_NFW_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcDensityProfileNFWClass
{
  /*< private > */
  NcDensityProfileClass parent_class;
};

struct _NcDensityProfileNFW
{
  /*< private > */
  NcDensityProfile parent_instance;
  gdouble Delta;
  gdouble r_Delta;
};

GType nc_density_profile_nfw_get_type (void) G_GNUC_CONST;

NcDensityProfile *nc_density_profile_nfw_new (void);

G_END_DECLS

#endif /* _NC_DENSITY_PROFILE_NFW_H_ */

