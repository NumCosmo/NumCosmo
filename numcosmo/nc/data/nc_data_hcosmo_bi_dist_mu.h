/***************************************************************************
 *            nc_data_hcosmo_bi_dist_mu.h
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_hcosmo_bi_dist_mu.h
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

#ifndef _NC_DATA_HCOSMO_BI_DIST_MU_H_
#define _NC_DATA_HCOSMO_BI_DIST_MU_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_data_gauss_diag.h>
#include <numcosmo/nc/background/nc_hcosmo_bi_dist.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_HCOSMO_BI_DIST_MU (nc_data_hcosmo_bi_dist_mu_get_type ())

G_DECLARE_FINAL_TYPE (NcDataHCosmoBIDistMu, nc_data_hcosmo_bi_dist_mu, NC, DATA_HCOSMO_BI_DIST_MU, NcmDataGaussDiag)

NcDataHCosmoBIDistMu *nc_data_hcosmo_bi_dist_mu_new (NcHCosmoBIDist * dist);
NcDataHCosmoBIDistMu *nc_data_hcosmo_bi_dist_mu_ref (NcDataHCosmoBIDistMu *dmu);

void nc_data_hcosmo_bi_dist_mu_free (NcDataHCosmoBIDistMu *dmu);
void nc_data_hcosmo_bi_dist_mu_clear (NcDataHCosmoBIDistMu **dmu);

void nc_data_hcosmo_bi_dist_mu_set_dist (NcDataHCosmoBIDistMu *dmu, NcHCosmoBIDist *dist);
void nc_data_hcosmo_bi_dist_mu_set_data (NcDataHCosmoBIDistMu *dmu, NcmVector *z, NcmMatrix *p, NcmVector *mu, NcmVector *sigma);

NcmVector *nc_data_hcosmo_bi_dist_mu_peek_z (NcDataHCosmoBIDistMu *dmu);
NcmMatrix *nc_data_hcosmo_bi_dist_mu_peek_p (NcDataHCosmoBIDistMu *dmu);

G_END_DECLS

#endif /* _NC_DATA_HCOSMO_BI_DIST_MU_H_ */

