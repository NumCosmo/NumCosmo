/***************************************************************************
 *            nc_data_planck_simall.h
 *
 *  Wed July 23 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_simall.h
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

#ifndef _NC_DATA_PLANCK_SIMALL_H_
#define _NC_DATA_PLANCK_SIMALL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_data.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_PLANCK_SIMALL (nc_data_planck_simall_get_type ())

G_DECLARE_FINAL_TYPE (NcDataPlanckSimall, nc_data_planck_simall, NC, DATA_PLANCK_SIMALL, NcmData)

NcDataPlanckSimall *nc_data_planck_simall_new (guint lmin,
                                               guint lmax,
                                               gdouble step_ee,
                                               NcmVector * prob_ee,
                                               gdouble step_bb,
                                               NcmVector * prob_bb,
                                               gdouble step_te,
                                               NcmVector * prob_te);
NcDataPlanckSimall *nc_data_planck_simall_new_from_file (const gchar *filename);

void nc_data_planck_simall_set_hipert_boltzmann (NcDataPlanckSimall *simall, NcHIPertBoltzmann *pb);
NcHIPertBoltzmann *nc_data_planck_simall_peek_hipert_boltzmann (NcDataPlanckSimall *simall);

G_END_DECLS

#endif /* _NC_DATA_PLANCK_SIMALL_H_ */

