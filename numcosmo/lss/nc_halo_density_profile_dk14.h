/***************************************************************************
 *            nc_halo_density_profile_dk14.h
 *
 *  Tue July 16 15:20:15 2019
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_dk14.h
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

#ifndef _NC_HALO_DENSITY_PROFILE_DK14_H_
#define _NC_HALO_DENSITY_PROFILE_DK14_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_density_profile.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE_DK14             (nc_halo_density_profile_dk14_get_type ())
#define NC_HALO_DENSITY_PROFILE_DK14(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_DENSITY_PROFILE_DK14, NcHaloDensityProfileDK14))
#define NC_HALO_DENSITY_PROFILE_DK14_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_DENSITY_PROFILE_DK14, NcHaloDensityProfileDK14Class))
#define NC_IS_HALO_DENSITY_PROFILE_DK14(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_DENSITY_PROFILE_DK14))
#define NC_IS_HALO_DENSITY_PROFILE_DK14_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_DENSITY_PROFILE_DK14))
#define NC_HALO_DENSITY_PROFILE_DK14_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_DENSITY_PROFILE_DK14, NcHaloDensityProfileDK14Class))

typedef struct _NcHaloDensityProfileDK14Class NcHaloDensityProfileDK14Class;
typedef struct _NcHaloDensityProfileDK14 NcHaloDensityProfileDK14;

struct _NcHaloDensityProfileDK14Class
{
  /*< private > */
  NcHaloDensityProfileClass parent_class;
};

/**
 * NcHaloDensityProfileDK14Params:
 * @NC_HALO_DENSITY_PROFILE_DK14_RT: the truncation radius, the radius where the profile steepens beyond the Einasto profile, in $kpc/h$
 * @NC_HALO_DENSITY_PROFILE_DK14_BETA: sharpness of the steepening
 * @NC_HALO_DENSITY_PROFILE_DK14_GAMMA: asymptotic negative slope of the steepening term
 * @NC_HALO_DENSITY_PROFILE_DK14_BE: amplitude of one factor of the outer profile.
 * @NC_HALO_DENSITY_PROFILE_DK14_SE: slope of the outer profile.
 *
 * The first three parameters, $\rho_s$, $r_s$ and $\alpha$, are the Einasto profile's parameters.
 * The transition term $f_{trans}$ is a function parametrized by $r_t$, $beta$ and $\gamma$. These two functions determine the inner
 * profile, whereas the outer profile is parametrized br $b_e$ and $s_e$.
 *
 */
typedef enum _NcHaloDensityProfileDK14Params
{
  NC_HALO_DENSITY_PROFILE_DK14_RT = NC_HALO_DENSITY_PROFILE_SPARAM_LEN,
  NC_HALO_DENSITY_PROFILE_DK14_BETA,
  NC_HALO_DENSITY_PROFILE_DK14_GAMMA,
  NC_HALO_DENSITY_PROFILE_DK14_BE,
  NC_HALO_DENSITY_PROFILE_DK14_SE,
  /* < private > */
  NC_HALO_DENSITY_PROFILE_DK14_SPARAM_LEN, /*< skip >*/
} NcHaloDensityProfileDK14Params;

/**
 * NcHaloDensityProfileDK14MethodParams:
 * @NC_HALO_DENSITY_PROFILE_DK14_MC2RHOSRS: FIXME
 * @NC_HALO_DENSITY_PROFILE_DK14_DIRECT_RHOSRS: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcHaloDensityProfileDK14MethodParams
{
  NC_HALO_DENSITY_PROFILE_DK14_MC2RHOSRS = 0,
  NC_HALO_DENSITY_PROFILE_DK14_DIRECT_RHOSRS,
  /* < private > */
} NcHaloDensityProfileDK14MethodParams;

struct _NcHaloDensityProfileDK14
{
  /*< private > */
  NcHaloDensityProfile parent_instance;
  gdouble Delta;
  gdouble r_Delta;
};

GType nc_halo_density_profile_dk14_get_type (void) G_GNUC_CONST;

NcHaloDensityProfile *nc_halo_density_profile_dk14_new (void);

#define NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_RT (1.0)
#define NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_BETA (4.0)
#define NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_GAMMA ()
#define NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_BE ()
#define NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_SE ()

#define NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_DK14_H_ */

