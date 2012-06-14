/***************************************************************************
 *            data_baryonic_oscillation.c
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
 * SECTION:data_baryonic_oscillation
 * @title: Baryonic Oscillation Data
 * @short_description: BAO Data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>

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
 * @bao_id: bao id #NcDataBaoSampleId
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_bao (NcDistance *dist, NcDataBaoSampleId bao_id)
{
  NcData *data;

  switch (bao_id)
  {
	case NC_DATA_BAO_R_DV_SAMPLE_PERCIVAL:
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
	case NC_DATA_BAO_DV_DV_SAMPLE_PERCIVAL:
	{
	  gdouble _data[2] = { NC_C_BAO_PERCIVAL_Dv_Dv, NC_C_BAO_PERCIVAL_SIGMA_Dv_Dv };
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
	case NC_DATA_BAO_A_SAMPLE_EISENSTEIN:
	{
	  gdouble _data[3] = { NC_C_BAO_EISENSTEIN_REDSHIFT, NC_C_BAO_EISENSTEIN_A, NC_C_BAO_EISENSTEIN_SIGMA_A };
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
	case NC_DATA_BAO_Dv_SAMPLE_EISENSTEIN:
	{
	  gdouble _data[3] = { NC_C_BAO_EISENSTEIN_REDSHIFT, NC_C_BAO_EISENSTEIN_Dv, NC_C_BAO_EISENSTEIN_SIGMA_Dv };
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
