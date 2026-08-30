/***************************************************************************
 *            nc_data_planck_commander.h
 *
 *  Wed July 23 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_commander.h
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

#ifndef _NC_DATA_PLANCK_COMMANDER_H_
#define _NC_DATA_PLANCK_COMMANDER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_data.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/nc/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_PLANCK_COMMANDER (nc_data_planck_commander_get_type ())

G_DECLARE_FINAL_TYPE (NcDataPlanckCommander, nc_data_planck_commander, NC, DATA_PLANCK_COMMANDER, NcmData)

NcDataPlanckCommander *nc_data_planck_commander_new (guint lmin,
                                                     guint lmax,
                                                     guint nbin,
                                                     guint delta_l,
                                                     NcmVector * cl2x_xa,
                                                     NcmVector * cl2x_ya,
                                                     NcmVector * cl2x_y2a,
                                                     NcmVector * mu,
                                                     NcmMatrix * cov,
                                                     NcmVector * mu_sigma);
NcDataPlanckCommander *nc_data_planck_commander_new_from_file (const gchar *filename);

void nc_data_planck_commander_set_hipert_boltzmann (NcDataPlanckCommander *cmd, NcHIPertBoltzmann *pb);
NcHIPertBoltzmann *nc_data_planck_commander_peek_hipert_boltzmann (NcDataPlanckCommander *cmd);

G_END_DECLS

#endif /* _NC_DATA_PLANCK_COMMANDER_H_ */

