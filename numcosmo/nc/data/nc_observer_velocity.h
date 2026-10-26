/***************************************************************************
 *            nc_observer_velocity.h
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_observer_velocity.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_OBSERVER_VELOCITY_H_
#define _NC_OBSERVER_VELOCITY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_OBSERVER_VELOCITY (nc_observer_velocity_get_type ())

G_DECLARE_FINAL_TYPE (NcObserverVelocity, nc_observer_velocity, NC, OBSERVER_VELOCITY, NcmModel)

/**
 * NcObserverVelocitySParams:
 * @NC_OBSERVER_VELOCITY_VX: $v_x$ [km/s]
 * @NC_OBSERVER_VELOCITY_VY: $v_y$ [km/s]
 * @NC_OBSERVER_VELOCITY_VZ: $v_z$ [km/s]
 *
 * #NcObserverVelocity parameters: the velocity of the observer relative to the rest frame
 * of the cosmological model, in the equatorial Cartesian frame.
 *
 */
typedef enum /*< enum,underscore_name=NC_OBSERVER_VELOCITY_SPARAMS,prefix=NC_OBSERVER_VELOCITY >*/
{
  NC_OBSERVER_VELOCITY_VX,
  NC_OBSERVER_VELOCITY_VY,
  NC_OBSERVER_VELOCITY_VZ,
  /* < private > */
  NC_OBSERVER_VELOCITY_SPARAM_LEN, /*< skip >*/
} NcObserverVelocitySParams;

NCM_MSET_MODEL_DECLARE_ID (nc_observer_velocity);

NcObserverVelocity *nc_observer_velocity_new (void);
NcObserverVelocity *nc_observer_velocity_ref (NcObserverVelocity *ov);

void nc_observer_velocity_free (NcObserverVelocity *ov);
void nc_observer_velocity_clear (NcObserverVelocity **ov);

gdouble nc_observer_velocity_speed (NcObserverVelocity *ov);
gdouble nc_observer_velocity_beta_n (NcObserverVelocity *ov, const gdouble nx, const gdouble ny, const gdouble nz);
gdouble nc_observer_velocity_doppler (NcObserverVelocity *ov, const gdouble nx, const gdouble ny, const gdouble nz);
gdouble nc_observer_velocity_z_rest (NcObserverVelocity *ov, const gdouble z_obs, const gdouble nx, const gdouble ny, const gdouble nz);

#define NC_OBSERVER_VELOCITY_DEFAULT_V (0.0)
#define NC_OBSERVER_VELOCITY_MAX_V (2000.0)
#define NC_OBSERVER_VELOCITY_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_OBSERVER_VELOCITY_H_ */

