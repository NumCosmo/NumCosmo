/***************************************************************************
 *            nc_data_galaxy_lowz.h
 *
 *  Thu Apr 22 10:37:39 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_GALAXY_LOWZ_H_
#define _NC_DATA_GALAXY_LOWZ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_GALAXY_LOWZ             (nc_data_galaxy_lowz_get_type ())

G_DECLARE_FINAL_TYPE (NcDataGalaxyLowz, nc_data_galaxy_lowz, NC, DATA_GALAXY_LOWZ, NcmDataGaussDiag)

NcDataGalaxyLowz *nc_data_galaxy_lowz_new_empty (void);

NcmMatrix *nc_data_galaxy_lowz_peek_position_z_theta_phi (NcDataGalaxyLowz *galaxy_lowz);

G_END_DECLS

#endif /* _NC_DATA_GALAXY_LOWZ_H_ */

