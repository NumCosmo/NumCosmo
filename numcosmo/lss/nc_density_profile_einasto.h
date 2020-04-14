/***************************************************************************
 *            nc_density_profile_einasto.h
 *
 *  Wed July 17 12:33:27 2019
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile_einasto.h
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

#ifndef _NC_DENSITY_PROFILE_EINASTO_H_
#define _NC_DENSITY_PROFILE_EINASTO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_density_profile.h>

G_BEGIN_DECLS

#define NC_TYPE_DENSITY_PROFILE_EINASTO             (nc_density_profile_einasto_get_type ())
#define NC_DENSITY_PROFILE_EINASTO(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DENSITY_PROFILE_EINASTO, NcDensityProfileEinasto))
#define NC_DENSITY_PROFILE_EINASTO_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DENSITY_PROFILE_EINASTO, NcDensityProfileEinastoClass))
#define NC_IS_DENSITY_PROFILE_EINASTO(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DENSITY_PROFILE_EINASTO))
#define NC_IS_DENSITY_PROFILE_EINASTO_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DENSITY_PROFILE_EINASTO))
#define NC_DENSITY_PROFILE_EINASTO_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DENSITY_PROFILE_EINASTO, NcDensityProfileEinastoClass))

typedef struct _NcDensityProfileEinastoClass NcDensityProfileEinastoClass;
typedef struct _NcDensityProfileEinasto NcDensityProfileEinasto;

struct _NcDensityProfileEinastoClass
{
  /*< private > */
  NcDensityProfileClass parent_class;
};

/**
 * NcDensityProfileEinastoParams:
 * @NC_DENSITY_PROFILE_EINASTO_RHOS: the central scale density in $M_\odot h^{2}/kpc^3$.
 * @NC_DENSITY_PROFILE_EINASTO_RS: the scale radius in physical $kpc/h$
 * @NC_DENSITY_PROFILE_EINASTO_ALPHA: determines how quickly the slope of the inner Einasto profile steepens 
 *
 * FIXME
 */
typedef enum _NcDensityProfileEinastoParams
{
  NC_DENSITY_PROFILE_EINASTO_RHOS = NC_DENSITY_PROFILE_SPARAM_LEN,
  NC_DENSITY_PROFILE_EINASTO_RS,
	NC_DENSITY_PROFILE_EINASTO_ALPHA,
  /* < private > */
  NC_DENSITY_PROFILE_EINASTO_SPARAM_LEN, /*< skip >*/
} NcDensityProfileEinastoParams;

struct _NcDensityProfileEinasto
{
  /*< private > */
  NcDensityProfile parent_instance;
  gdouble Delta;
  gdouble r_Delta;
};

GType nc_density_profile_einasto_get_type (void) G_GNUC_CONST;

NcDensityProfile *nc_density_profile_einasto_new (void);

#define NC_DENSITY_PROFILE_EINASTO_DEFAULT_RHOS  (0.0)
#define NC_DENSITY_PROFILE_EINASTO_DEFAULT_RS    (0.0)
#define NC_DENSITY_PROFILE_EINASTO_DEFAULT_ALPHA (0.0) 

#define NC_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_DENSITY_PROFILE_EINASTO_H_ */

