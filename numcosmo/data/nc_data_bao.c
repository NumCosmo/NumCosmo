/***************************************************************************
 *            nc_data_bao.c
 *
 *  Thu Apr 22 15:31:32 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

/**
 * SECTION:nc_data_bao
 * @title: Baryonic Oscillation Data
 * @short_description: BAO Data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao.h"
#include "data/data_gaussian.h"

/***************************************************************************
 * BAO Percival/Eisenstein priors data (arXiv:0705.3323)
 *
 ****************************************************************************/

gdouble nc_bao_distance_priors_percival_z[]       = {    0.2,   0.35 };
gdouble nc_bao_distance_priors_percival_bestfit[] = { 0.1980, 0.1094 };
gdouble nc_bao_distance_priors_percival_inv_cov[] =
{  35059.0, -24031.0,
  -24031.0, 108300.0 };

gdouble nc_bao_distance_priors_percival_data[] =
{ 0.20, 0.1980,  35059.0, -24031.0,
  0.35, 0.1094, -24031.0, 108300.0 };

/**
 * nc_data_bao:
 * @dist: a #NcDistance
 * @bao_id: bao id #NcDataBaoId
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_bao (NcDistance *dist, NcDataBaoId bao_id)
{
  NcData *data;

  switch (bao_id)
  {
	case NC_DATA_BAO_R_DV_PERCIVAL:
	{
	  NcmMatrix *cm = ncm_matrix_new_data_static (nc_bao_distance_priors_percival_data, 2, 4);
	  NcmMSetFunc *func = nc_distance_func1_new (dist, &nc_distance_bao_r_Dv);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_COV);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "Percival BAO Sample R-Dv";
	  break;
	}
	case NC_DATA_BAO_DV_DV_PERCIVAL:
	{
	  gdouble _data[2] = { ncm_c_bao_percival_DV_DV (), ncm_c_bao_percival_sigma_DV_DV () };
	  NcmMatrix *cm = ncm_matrix_new_data_static (_data, 1, 2);
	  NcmMSetFunc *func = nc_distance_func0_new (dist, &nc_distance_dilation_scale_ratio);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_SIGMA);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "Percival BAO Sample Dv-Dv";
	  break;
	}
	case NC_DATA_BAO_A_EISENSTEIN:
	{
	  gdouble _data[3] = { ncm_c_bao_eisenstein_z (), ncm_c_bao_eisenstein_A (), ncm_c_bao_eisenstein_sigma_A () };
	  NcmMatrix *cm = ncm_matrix_new_data_static (_data, 1, 3);
	  NcmMSetFunc *func = nc_distance_func1_new (dist, &nc_distance_bao_A_scale);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_SIGMA);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "Eisenstein BAO Sample A";
	  break;
	}
	case NC_DATA_BAO_DV_EISENSTEIN:
	{
	  gdouble _data[3] = { ncm_c_bao_eisenstein_z (), ncm_c_bao_eisenstein_DV (), ncm_c_bao_eisenstein_sigma_DV () };
	  NcmMatrix *cm = ncm_matrix_new_data_static (_data, 1, 3);
	  NcmMSetFunc *func = nc_distance_func1_new (dist, &nc_distance_dilation_scale);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_SIGMA);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "Eisenstein BAO Sample Dv";
	  break;
	}
	default:
	  g_error ("Bao dataset %d do not exist", bao_id);
	  data = NULL;
  }

  return data;
}
