/***************************************************************************
 *            test_nc_ccl_Pk.c
 *
 *  Thu February 28 10:57:18 2019
 *  Copyright  2019  Mariana Penna Lima 
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * test_nc_ccl_dist.c
 *
 * Copyright (C) 2019 - Mariana Penna Lima
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <ccl.h>

typedef struct _TestNcCCLDist
{
  ccl_cosmology * ccl_cosmo;
  NcHICosmo *cosmo;
  NcDistance *dist;
} TestNcCCLDist;

void test_nc_ccl_dist_new_model1_noNeff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model2_noNeff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model3_noNeff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model4_noNeff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model5_noNeff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model1_Neff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model2_Neff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model3_Neff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model4_Neff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model5_Neff_nomnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model1_noNeff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model2_noNeff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model3_noNeff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model4_noNeff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model5_noNeff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model1_Neff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model2_Neff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model3_Neff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model4_Neff_mnu  (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_new_model5_Neff_mnu  (TestNcCCLDist *test, gconstpointer pdata);

void test_nc_ccl_dist_free (TestNcCCLDist *test, gconstpointer pdata);

void test_nc_ccl_dist_cmp (TestNcCCLDist *test, gconstpointer pdata);

void test_nc_ccl_dist_traps (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_invalid_st (TestNcCCLDist *test, gconstpointer pdata);

typedef struct _ccl_params_data 
{
  gdouble Omega_c;
  gdouble Omega_b;
  gdouble h;
  gdouble A_s;
  gdouble n_s;
  gdouble sigma8;
  gdouble Neff[2];
  gdouble mnu[5][3];
  gdouble Omega_v[5];
  gdouble Omega_k[5];
  gdouble w_0[5];
  gdouble w_a[5];
} ccl_params_data;

static void
ccl_params_data_init (ccl_params_data *data)
{
  gdouble mnu[5][3]  = {{0.04, 0.0, 0.0}, {0.05, 0.01, 0.0}, {0.03, 0.02, 0.04}, {0.05, 0.0, 0.0}, {0.03, 0.02, 0.0}};
  gdouble Omega_v[5] = { 0.7,  0.7,  0.7,  0.65, 0.75};
  gdouble w_0[5]     = {-1.0, -0.9, -0.9, -0.9, -0.9};
  gdouble w_a[5]     = { 0.0,  0.0,  0.1,  0.1,  0.1};
	gint i;

	data->Omega_c  = 0.25;
  data->Omega_b  = 0.05;
  data->h        = 0.7;
  data->A_s      = 2.1e-9;
  data->n_s      = 0.96;
  data->sigma8   = 0.8;
	data->Neff[0]  = 0.0;
	data->Neff[1]  = 3.0;

	memcpy (data->Omega_v, Omega_v, sizeof (gdouble) * 5);
	memcpy (data->w_0,     w_0,     sizeof (gdouble) * 5);
	memcpy (data->w_a,     w_a,     sizeof (gdouble) * 5);
	memcpy (data->mnu,     mnu,     sizeof (gdouble) * 5 * 3);

	for (i = 0; i < 5; i++) 
  {
    data->Omega_k[i] = 1.0 - data->Omega_c - data->Omega_b - data->Omega_v[i];
  }
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();
  
  g_test_add ("/nc/ccl/model1/noNeff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model1_noNeff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model2/noNeff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model2_noNeff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model3/noNeff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model3_noNeff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model4/noNeff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model4_noNeff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model5/noNeff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model5_noNeff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model1/Neff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model1_Neff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model2/Neff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model2_Neff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model3/Neff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model3_Neff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model4/Neff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model4_Neff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model5/Neff/nomnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model5_Neff_nomnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);

  g_test_add ("/nc/ccl/model1/noNeff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model1_noNeff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model2/noNeff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model2_noNeff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model3/noNeff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model3_noNeff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model4/noNeff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model4_noNeff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model5/noNeff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model5_noNeff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model1/Neff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model1_Neff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model2/Neff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model2_Neff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model3/Neff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model3_Neff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model4/Neff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model4_Neff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);
  g_test_add ("/nc/ccl/model5/Neff/mnu/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_new_model5_Neff_mnu ,
              &test_nc_ccl_dist_cmp,
              &test_nc_ccl_dist_free);

/*
  g_test_add ("/nc/ccl/model2/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_comoving_new_model2,
              &test_nc_ccl_cmp_dist,
              &test_nc_ccl_free);
              
  g_test_add ("/nc/ccl/model3/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_comoving_new_model3,
              &test_nc_ccl_cmp_dist,
              &test_nc_ccl_free);
              
  g_test_add ("/nc/ccl/model4/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_comoving_new_model4,
              &test_nc_ccl_cmp_dist,
              &test_nc_ccl_free);
              
  g_test_add ("/nc/ccl/model5/dist/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_comoving_new_model5,
              &test_nc_ccl_cmp_dist,
              &test_nc_ccl_free);

  g_test_add ("/nc/ccl/model1/dist_modulus/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_modulus_new_model1,
              &test_nc_ccl_cmp_dist,
              &test_nc_ccl_free);

  g_test_add ("/nc/ccl/model5/dist_modulus/cmp", TestNcCCLDist, NULL,
              &test_nc_ccl_dist_modulus_new_model5,
              &test_nc_ccl_cmp_dist,
              &test_nc_ccl_free);
 
  g_test_add ("/nc/ccl/model1/bbks/traps", TestNcCCLDist, NULL,
              &test_nc_ccl_bbks_new_model1,
              &test_nc_ccl_bbks_traps,
              &test_nc_ccl_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/nc/ccl/model1/bbks/invalid/st/subprocess", TestNcCCLDist, NULL,
              &test_nc_ccl_bbks_new_model1,
              &test_nc_ccl_bbks_invalid_st,
              &test_nc_ccl_free);
#endif
*/
  g_test_run ();
}

static gdouble 
neff_from_n_ur_n_ncdm (gdouble N_ur, gdouble N_ncdm)
{
  gdouble Neff = N_ur + N_ncdm * pow (ccl_constants.TNCDM, 4.0) / pow (4.0/11.0, 4.0/3.0);
  return Neff;
}

static ccl_cosmology *
test_nc_create_ccl_cosmo (gint i_model, gint i_Neff, gboolean has_mnu)
{
  ccl_params_data *data    = g_new0 (ccl_params_data, 1);
  ccl_configuration config = default_config;
  gint status              = 0;
  ccl_parameters params;
  ccl_cosmology *cosmo;

  ccl_params_data_init (data);

	if (!has_mnu)
	{
		gdouble mnu = 0.0;
    params = ccl_parameters_create (data->Omega_c, data->Omega_b, data->Omega_k[i_model], 
                                    data->Neff[i_Neff], &mnu, ccl_mnu_sum, data->w_0[i_model], 
                                    data->w_a[i_model], data->h, data->A_s, data->n_s, 
                                    -1, -1, -1, -1, NULL, NULL, &status);
	}
	else
	{
		params = ccl_parameters_create (data->Omega_c, data->Omega_b, data->Omega_k[i_model], 
		                                data->Neff[i_Neff], data->mnu[i_model * 3], 
		                                ccl_mnu_list, data->w_0[i_model], data->w_a[i_model], 
		                                data->h, data->A_s, data->n_s, -1, -1, -1, -1, NULL, NULL, &status);
	}

	cosmo = ccl_cosmology_create (params, config);

  g_free (data);

  return cosmo;
}

#define MAX_Z (1200.0)

/* NOMNU */

void
test_nc_ccl_dist_new_model1_noNeff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (0, 0, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model2_noNeff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, 0, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model3_noNeff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, 0, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model4_noNeff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, 0, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model5_noNeff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (4, 0, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model1_Neff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (0, 1, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model2_Neff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, 1, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model3_Neff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, 1, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model4_Neff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, 1, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model5_Neff_nomnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (4, 1, FALSE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

/* MNU */

void
test_nc_ccl_dist_new_model1_noNeff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (0, 0, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model2_noNeff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, 0, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model3_noNeff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, 0, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model4_noNeff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, 0, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model5_noNeff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (4, 0, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model1_Neff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (0, 1, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model2_Neff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, 1, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model3_Neff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, 1, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model4_Neff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, 1, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_new_model5_Neff_mnu (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (4, 1, TRUE);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_dist_free (TestNcCCLDist *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
	NCM_TEST_FREE (nc_distance_free, test->dist);

  ccl_parameters_free (&test->ccl_cosmo->params);
  ccl_cosmology_free (test->ccl_cosmo);  
}

void
test_nc_ccl_dist_cmp (TestNcCCLDist *test, gconstpointer pdata)
{
  const gint ntests       = 10000;
  const gdouble z_max     = 1000.0;
  const gdouble tol       = 5.0e-5;
	const gdouble RH_Mpc    = nc_hicosmo_RH_Mpc (test->cosmo);
  gint status = 0;
  gint i;
	
  for (i = 0; i < ntests; i++)
  {
    const gdouble z     = z_max / (1.0 * ntests) * i;
    const gdouble a     = 1.0 / (1.0 + z);
    const gdouble E2    = nc_hicosmo_E2 (test->cosmo, z);
    const gdouble cclE2 = gsl_pow_2 (ccl_h_over_h0 (test->ccl_cosmo, a, &status));

    if (TRUE)
    {
      const gdouble E2Omega_x    = nc_hicosmo_de_E2Omega_de (NC_HICOSMO_DE (test->cosmo), z);
      const gdouble cclE2Omega_x = ccl_omega_x (test->ccl_cosmo, a, ccl_species_l_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_x, ==, cclE2Omega_x, tol, 0.0);
      
      /* printf ("(% 22.15g, % 22.15g) | OMEGA_X: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_x, cclE2Omega_x, E2Omega_x / cclE2Omega_x - 1.0); */
    }
		
    if (TRUE)
    {
      const gdouble E2Omega_m    = nc_hicosmo_E2Omega_m (test->cosmo, z);
      const gdouble cclE2Omega_m = ccl_omega_x (test->ccl_cosmo, a, ccl_species_m_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_m, ==, cclE2Omega_m, tol, 0.0);
      printf ("(% 22.15g, % 22.15g) | OMEGA_M: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_m, cclE2Omega_m, E2Omega_m / cclE2Omega_m - 1.0);
    }

    if (TRUE)
    {
      const gdouble E2Omega_g    = nc_hicosmo_E2Omega_g (test->cosmo, z);
      const gdouble cclE2Omega_g = ccl_omega_x (test->ccl_cosmo, a, ccl_species_g_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_g, ==, cclE2Omega_g, tol, 0.0);
      /*printf ("(% 22.15g, % 22.15g) | OMEGA_G: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_g, cclE2Omega_g, E2Omega_g / cclE2Omega_g - 1.0);*/    
    }

    if (TRUE)
    {
      const gdouble E2Omega_nu    = nc_hicosmo_E2Omega_nu (test->cosmo, z);
      const gdouble cclE2Omega_nu = ccl_omega_x (test->ccl_cosmo, a, ccl_species_ur_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_nu, ==, cclE2Omega_nu, tol, 0.0);
      /*printf ("(% 22.15g, % 22.15g) | OMEGA_NU: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_nu, cclE2Omega_nu, E2Omega_nu / cclE2Omega_nu - 1.0);*/    
    }
    
    if (TRUE)
    {
      const gdouble E2Omega_k    = nc_hicosmo_E2Omega_k (test->cosmo, z);
      const gdouble cclE2Omega_k = ccl_omega_x (test->ccl_cosmo, a, ccl_species_k_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_k, ==, cclE2Omega_k, tol, 0.0);
      /* printf ("(% 22.15g, % 22.15g) | OMEGA_K: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_k, cclE2Omega_k, E2Omega_k / cclE2Omega_k - 1.0); */
    }

    if (TRUE)
    {
      ncm_assert_cmpdouble_e (E2, ==, cclE2, tol, 0.0);
      /*printf ("(% 22.15g, % 22.15g) | E2 % 22.15g % 22.15g %17.10e\n", z, a, E2, cclE2, E2 / cclE2 - 1.0);*/ 
    }

    if (TRUE)
    {
			const gdouble chi    = RH_Mpc * nc_distance_comoving (test->dist, test->cosmo, z);
			const gdouble cclchi = ccl_comoving_radial_distance (test->ccl_cosmo, a, &status);

      ncm_assert_cmpdouble_e (chi, ==, cclchi, tol, 0.0);
      /*printf ("(% 22.15g, % 22.15g) | chi % 22.15g % 22.15g %17.10e\n", z, a, chi, cclchi, chi / cclchi - 1.0);*/ 
    }

    if (TRUE)
    {
			const gdouble Dt    = RH_Mpc * nc_distance_transverse (test->dist, test->cosmo, z);
			const gdouble cclDt = ccl_comoving_angular_distance (test->ccl_cosmo, a, &status);

      ncm_assert_cmpdouble_e (Dt, ==, cclDt, tol, 0.0);
      /* printf ("(% 22.15g, % 22.15g) | Dt % 22.15g % 22.15g %17.10e\n", z, a, Dt, cclDt, Dt / cclDt - 1.0); */
    }

    if (TRUE)
    {
			const gdouble Dl    = RH_Mpc * nc_distance_luminosity (test->dist, test->cosmo, z);
			const gdouble cclDl = ccl_luminosity_distance (test->ccl_cosmo, a, &status);

      ncm_assert_cmpdouble_e (Dl, ==, cclDl, tol, 0.0);
      /* printf ("(% 22.15g, % 22.15g) | Dl % 22.15g % 22.15g %17.10e\n", z, a, Dl, cclDl, Dl / cclDl - 1.0); */
    }

    if (a != 1.0)
    {
			const gdouble mu    = nc_distance_dmodulus (test->dist, test->cosmo, z) + 5.0 * log10 (RH_Mpc);
			const gdouble cclmu = ccl_distance_modulus (test->ccl_cosmo, a, &status);

      ncm_assert_cmpdouble_e (mu, ==, cclmu, tol, 0.0);
      /* printf ("(% 22.15g, % 22.15g) | mu % 22.15g % 22.15g %17.10e\n", z, a, mu, cclmu, mu / cclmu - 1.0); */
    }
	}
}

void
test_nc_ccl_bbks_traps (TestNcCCLDist *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/nc/ccl/model1/bbks/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_nc_ccl_bbks_invalid_st (TestNcCCLDist *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}
