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

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_EINASTO             (nc_halo_density_profile_einasto_get_type ())
#define NC_HALO_DENSITY_PROFILE_EINASTO(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_DENSITY_PROFILE_EINASTO, NcHaloDensityProfileEinasto))
#define NC_HALO_DENSITY_PROFILE_EINASTO_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_DENSITY_PROFILE_EINASTO, NcHaloDensityProfileEinastoClass))
#define NC_IS_HALO_DENSITY_PROFILE_EINASTO(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_DENSITY_PROFILE_EINASTO))
#define NC_IS_HALO_DENSITY_PROFILE_EINASTO_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_DENSITY_PROFILE_EINASTO))
#define NC_HALO_DENSITY_PROFILE_EINASTO_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_DENSITY_PROFILE_EINASTO, NcHaloDensityProfileEinastoClass))

typedef struct _NcHaloDensityProfileEinastoClass NcHaloDensityProfileEinastoClass;
typedef struct _NcHaloDensityProfileEinasto NcHaloDensityProfileEinasto;

struct _NcHaloDensityProfileEinastoClass
{
  /*< private > */
  NcHaloDensityProfileClass parent_class;
};

/**
 * NcHaloDensityProfileEinastoParams:
 * @NC_HALO_DENSITY_PROFILE_EINASTO_ALPHA: determines how quickly the slope of the inner Einasto profile steepens
 *
 * FIXME
 */
typedef enum _NcHaloDensityProfileEinastoParams
{
  NC_HALO_DENSITY_PROFILE_EINASTO_ALPHA = NC_HALO_DENSITY_PROFILE_SPARAM_LEN,
  /* < private > */
  NC_HALO_DENSITY_PROFILE_EINASTO_SPARAM_LEN, /*< skip >*/
} NcHaloDensityProfileEinastoParams;

#define NC_HALO_DENSITY_PROFILE_EINASTO_LOCAL_SPARAM_LEN (NC_HALO_DENSITY_PROFILE_EINASTO_SPARAM_LEN - NC_HALO_DENSITY_PROFILE_SPARAM_LEN)

struct _NcHaloDensityProfileEinasto
{
  /*< private > */
  NcHaloDensityProfile parent_instance;
};

GType nc_halo_density_profile_einasto_get_type (void) G_GNUC_CONST;

NcHaloDensityProfileEinasto *nc_halo_density_profile_einasto_new (const NcHaloDensityProfileMassDef mdef, const gdouble Delta);

#define NC_HALO_DENSITY_PROFILE_EINASTO_DEFAULT_ALPHA (0.25)
#define NC_HALO_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_EINASTO_H_ */

