/***************************************************************************
 *            nc_data_planck_plik_lite.h
 *
 *  Fri June 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_plik_lite.h
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

#ifndef _NC_DATA_PLANCK_PLIK_LITE_H_
#define _NC_DATA_PLANCK_PLIK_LITE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_data_gauss_cov.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_PLANCK_PLIK_LITE (nc_data_planck_plik_lite_get_type ())

G_DECLARE_FINAL_TYPE (NcDataPlanckPlikLite, nc_data_planck_plik_lite, NC, DATA_PLANCK_PLIK_LITE, NcmDataGaussCov)

NcDataPlanckPlikLite *nc_data_planck_plik_lite_new_from_file (const gchar *filename);

void nc_data_planck_plik_lite_set_hipert_boltzmann (NcDataPlanckPlikLite *plik, NcHIPertBoltzmann *pb);
NcHIPertBoltzmann *nc_data_planck_plik_lite_peek_hipert_boltzmann (NcDataPlanckPlikLite *plik);

G_END_DECLS

#endif /* _NC_DATA_PLANCK_PLIK_LITE_H_ */

