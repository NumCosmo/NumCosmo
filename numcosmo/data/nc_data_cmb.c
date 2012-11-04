/***************************************************************************
 *            nc_data_cmb.c
 *
 *  Thu Apr 22 15:56:44 2010
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
 * SECTION:nc_data_cmb
 * @title: Cosmic Microwave Background Data
 * @short_description: CMB Data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cmb.h"
#include "data/data_gaussian.h"

/***************************************************************************
 * WMAP5 Distance priors data (arXiv:0803.0547), (astro-ph/0604051)
 *
 ****************************************************************************/

gdouble nc_cmb_distance_priors_wmap5_bestfit[] = { 302.1000,    1.710, 1090.04000 };
gdouble nc_cmb_distance_priors_wmap5_inv_cov[] =
{ 1.8000,   27.968,   -1.10300,
  27.968, 5667.577,  -92.26300,
  -1.103,  -92.263,    2.92300 };

/*     BF,        COV MATRIX  */
gdouble nc_cmb_distance_priors_wmap5_distpriors[] =
{ 302.10,   1.800,   27.968,  -1.103,
  1.71,    27.968, 5667.577, -92.263,
  1090.04, -1.103,  -92.263,   2.923 };

/***************************************************************************
 * WMAP7 Distance priors data (arXiv:1001.4538): tables 9 and 10
 *
 ****************************************************************************/

gdouble nc_cmb_distance_priors_wmap7_bestfit[] = { 302.0900,   1.725, 1091.30000 };
gdouble nc_cmb_distance_priors_wmap7_inv_cov[] =
{ 2.3050,   29.698,   -1.333,
  29.698, 6825.270,  -113.18,
  -1.333,  -113.18,    3.414 };

/*     BF,        COV MATRIX  */
gdouble nc_cmb_distance_priors_wmap7_distpriors[] =
{ 302.09,   2.3050,   29.698,   -1.333,
  1.725,    29.698, 6825.270,  -113.18,
  1091.30, -1.333,  -113.18,    3.414 };


/**
 * nc_data_cmb:
 * @dist: a #NcDistance
 * @cmb_id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_cmb (NcDistance *dist, NcDataCMBId cmb_id)
{
  NcData *data;

  switch (cmb_id)
  {
	case NC_DATA_CMB_SHIFT_PARAMETER_WMAP3:
	{
	  gdouble _data[3] = { ncm_c_wmap3_cmb_z (), ncm_c_wmap3_cmb_R (), ncm_c_wmap3_cmb_sigma_R () };
	  NcmMatrix *cm = ncm_matrix_new_data_static (_data, 1, 3);
	  NcmMSetFunc *func = nc_distance_func1_new (dist, &nc_distance_shift_parameter);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_SIGMA);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "WMAP3 shift parameter";
	  break;
	}
	case NC_DATA_CMB_SHIFT_PARAMETER_WMAP5:
	{
	  gdouble _data[3] = { ncm_c_wmap5_cmb_z (), ncm_c_wmap5_cmb_R (), ncm_c_wmap5_cmb_sigma_R () };
	  NcmMatrix *cm = ncm_matrix_new_data_static (_data, 1, 3);
	  NcmMSetFunc *func = nc_distance_func1_new (dist, &nc_distance_shift_parameter);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_SIGMA);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "WMAP5 shift parameter";
	  break;
	}
	case NC_DATA_CMB_DISTANCE_PRIORS_WMAP5:
	{
	  NcmMatrix *cm = ncm_matrix_new_data_static (nc_cmb_distance_priors_wmap5_distpriors, 3, 4);
	  NcmMSetFunc *func[3];
	  func[0] = nc_distance_func0_new (dist, &nc_distance_acoustic_scale);
	  func[1] = nc_distance_func0_new (dist, &nc_distance_shift_parameter_lss);
	  func[2] = nc_distance_func0_new (dist, &nc_distance_decoupling_redshift);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_COV);
	  nc_data_gaussian_set_func_array (data, func, 3);
	  ncm_mset_func_free (func[0]);
	  ncm_mset_func_free (func[1]);
	  ncm_mset_func_free (func[2]);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "WMAP5 distance priors";
	  break;
	}
	case NC_DATA_CMB_SHIFT_PARAMETER_WMAP7:
	{
	  gdouble _data[3] = { ncm_c_wmap7_cmb_z (), ncm_c_wmap7_cmb_R (), ncm_c_wmap7_cmb_sigma_R () };
	  NcmMatrix *cm = ncm_matrix_new_data_static (_data, 1, 3);
	  NcmMSetFunc *func = nc_distance_func1_new (dist, &nc_distance_shift_parameter);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_SIGMA);
	  nc_data_gaussian_set_func (data, func);
	  ncm_mset_func_free (func);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "WMAP7 shift parameter";
	  break;
	}  
	case NC_DATA_CMB_DISTANCE_PRIORS_WMAP7:
	{
	  NcmMatrix *cm = ncm_matrix_new_data_static (nc_cmb_distance_priors_wmap7_distpriors, 3, 4);
	  NcmMSetFunc *func[3];
	  func[0] = nc_distance_func0_new (dist, &nc_distance_acoustic_scale);
	  func[1] = nc_distance_func0_new (dist, &nc_distance_shift_parameter_lss);
	  func[2] = nc_distance_func0_new (dist, &nc_distance_decoupling_redshift);
	  data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_COV);
	  nc_data_gaussian_set_func_array (data, func, 3);
	  ncm_mset_func_free (func[0]);
	  ncm_mset_func_free (func[1]);
	  ncm_mset_func_free (func[2]);
	  nc_data_gaussian_init_from_matrix (data, cm, NULL);
	  ncm_matrix_free (cm);
	  data->name = "WMAP7 distance priors";
	  break;
	}
	default:
	  g_error ("CMB dataset %d do not exist", cmb_id);
	  data = NULL;
  }

  return data;
}
