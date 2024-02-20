/***************************************************************************
 *            nc_data_planck_lkl.h
 *
 *  Tue October 20 15:24:54 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_lkl.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_PLANCK_LKL_H_
#define _NC_DATA_PLANCK_LKL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/data/nc_data_cmb.h>
#include <numcosmo/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_PLANCK_LKL (nc_data_planck_lkl_get_type ())

G_DECLARE_FINAL_TYPE (NcDataPlanckLKL, nc_data_planck_lkl, NC, DATA_PLANCK_LKL, NcmData)

/**
 * NcDataPlanckLKLType:
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_TT: Planck 2018 baseline low-ell likelihood for TT.
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_EE: Planck 2018 baseline low-ell likelihood for EE.
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_BB: Planck 2018 baseline low-ell likelihood for BB.
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_EEBB: Planck 2018 baseline low-ell likelihood for TE.
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TT: Planck 2018 baseline high-ell likelihood for TT.
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TTTEEE: Planck 2018 baseline high-ell likelihood for TT, TE and EE.
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TT_LITE: Planck 2018 baseline high-ell likelihood for TT (lite).
 * @NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TTTEEE_LITE: Planck 2018 baseline high-ell likelihood for TT, TE and EE (lite).
 *
 * The Planck likelihood types.
 *
 */
typedef enum _NcDataPlanckLKLType
{
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_TT,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_EE,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_BB,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_EEBB,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TT,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TTTEEE,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TT_LITE,
  NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_HIGHL_TTTEEE_LITE,
  /* < private > */
  NC_DATA_PLANCK_LKL_TYPE_LENGTH,
} NcDataPlanckLKLType;

NcDataPlanckLKL *nc_data_planck_lkl_new (const gchar *filename);
NcDataPlanckLKL *nc_data_planck_lkl_full_new (const gchar *filename, NcHIPertBoltzmann *pb);
NcDataPlanckLKL *nc_data_planck_lkl_full_new_id (NcDataPlanckLKLType id, NcHIPertBoltzmann *pb);

const gchar *nc_data_planck_lkl_get_param_name (NcDataPlanckLKL *plik, guint i);
gchar **nc_data_planck_lkl_get_param_names (NcDataPlanckLKL *plik);

void nc_data_planck_lkl_set_hipert_boltzmann (NcDataPlanckLKL *plik, NcHIPertBoltzmann *pb);
void nc_data_planck_lkl_download_baseline (const gchar *dir);

G_END_DECLS

#endif /* _NC_DATA_PLANCK_LKL_H_ */

