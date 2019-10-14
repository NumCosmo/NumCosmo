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

void test_nc_ccl_dist_new (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_free (TestNcCCLDist *test, gconstpointer pdata);

void test_nc_ccl_dist_cmp (TestNcCCLDist *test, gconstpointer pdata);

void test_nc_ccl_dist_traps (TestNcCCLDist *test, gconstpointer pdata);
void test_nc_ccl_dist_invalid_st (TestNcCCLDist *test, gconstpointer pdata);

typedef struct _ccl_params_test_data
{
  const gchar *name; 
  gboolean use_Omega_k;
  gboolean use_Neff_func;
  gdouble Omega_c;
  gdouble Omega_b;
  gdouble h;
  gdouble A_s;
  gdouble n_s;
  gdouble Neff;
  gdouble N_ur;
  gdouble N_ncdm;
  gdouble mnu[3];
  gdouble Omega_k;
  gdouble Omega_v;
  gdouble w_0;
  gdouble w_a;
} ccl_params_test_data;

ccl_params_test_data test_data[29] = {
/* Name           use_Ok  use_f  Oc    Ob    h    As      ns    Neff Nur   Ncdm   m1    m2    m3      Ok    Ode   w0    wa */
  {"model_0",          FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.70, -1.0, 0.0},
  {"model_1",          FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.70, -0.9, 0.0},
  {"model_2",          FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.70, -0.9, 0.1},
  {"model_3",          FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.65, -0.9, 0.1},
  {"model_4",          FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.75, -0.9, 0.1},
  {"mnu/model_0",      FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.04, 0.0,  0.0 },  0.0,  0.70, -1.0, 0.0},
  {"mnu/model_1",      FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.05, 0.01, 0.0 },  0.0,  0.70, -0.9, 0.0},
  {"mnu/model_2",      FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.046, 0.0, 0.0, {0.03, 0.02, 0.04},  0.0,  0.70, -0.9, 0.1},
  {"mnu/model_3",      FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.05, 0.0,  0.0 },  0.0,  0.65, -0.9, 0.1},
  {"mnu/model_4",      FALSE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.03, 0.02, 0.0 },  0.0,  0.75, -0.9, 0.1},

  {"flat/nonu",         TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.0,  0.0 , 0.0 },  0.0,  0.00, -1.0, 0.0},
  {"pos_curv/nonu",     TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.0,  0.0 , 0.0 },  0.01, 0.00, -1.0, 0.0},
  {"neg_curv/nonu",     TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.0,   0.0, 0.0, {0.0,  0.0 , 0.0 }, -0.01, 0.00, -1.0, 0.0},
  {"flat/massnu1",      TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   2.0, 1.0, {0.0,  0.0 , 0.1 },  0.0,  0.00, -1.0, 0.0},
  {"flat/massnu2",      TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 3.0, {0.03, 0.03, 0.1 },  0.0,  0.00, -1.0, 0.0},
  {"flat/massnu3",      TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 3.0, {0.03, 0.05, 0.1 },  0.0,  0.00, -1.0, 0.0},
  {"flat_manynu1",      TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 6.0,   0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.00, -1.0, 0.0},
  {"neg_curv/massnu1",  TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   4.0, 2.0, {0.0,  0.03, 0.1 }, -0.01, 0.00, -1.0, 0.0},
  {"pos_curv/manynu1",  TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   3.0, 3.0, {0.03, 0.05, 0.1 },  0.01, 0.00, -1.0, 0.0},

  {"CCL1",              TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.046, 0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.00, -1.0, 0.0},
  {"CCL2",              TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.046, 0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.00, -0.9, 0.0},
  {"CCL3",              TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.046, 0.0, 0.0, {0.0,  0.0,  0.0 },  0.0,  0.00, -0.9, 0.1},
  {"CCL4",              TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.046, 0.0, 0.0, {0.0,  0.0,  0.0 },  0.05, 0.00, -0.9, 0.1},
  {"CCL5",              TRUE, FALSE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 3.046, 0.0, 0.0, {0.0,  0.0,  0.0 }, -0.05, 0.00, -0.9, 0.1},
  {"CCL7",              TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   2.0, 1.0, {0.04, 0.0,  0.0 },  0.0,  0.00, -1.0, 0.0},
  {"CCL8",              TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   1.0, 2.0, {0.05, 0.01, 0.0 },  0.0,  0.00, -0.9, 0.0},
  {"CCL9",              TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   0.0, 3.0, {0.03, 0.02, 0.04},  0.0,  0.00, -0.9, 0.1},
  {"CCL10",             TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   2.0, 1.0, {0.05, 0.0,  0.0 },  0.05, 0.00, -0.9, 0.1},
  {"CCL11",             TRUE,  TRUE, 0.25, 0.05, 0.7, 2.1e-9, 0.96, 0.0,   1.0, 2.0, {0.03, 0.02, 0.0 }, -0.05, 0.00, -0.9, 0.1},
};

gint
main (gint argc, gchar *argv[])
{
  gint i;
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < sizeof (test_data) / sizeof (test_data[0]); i++)
  {
    gchar *name = g_strdup_printf ("/nc/ccl/%s/dist/cmp", test_data[i].name);
    g_test_add (name, TestNcCCLDist, &test_data[i],
                &test_nc_ccl_dist_new,
                &test_nc_ccl_dist_cmp,
                &test_nc_ccl_dist_free);
    g_free (name);
  }

  g_test_add ("/nc/ccl/dist/traps", TestNcCCLDist, &test_data[0],
              &test_nc_ccl_dist_new,
              &test_nc_ccl_dist_traps,
              &test_nc_ccl_dist_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/nc/ccl/dist/invalid/st/subprocess", TestNcCCLDist, &test_data[0],
              &test_nc_ccl_dist_new,
              &test_nc_ccl_dist_invalid_st,
              &test_nc_ccl_dist_free);
#endif

  g_test_run ();
}

static gdouble 
neff_from_n_ur_n_ncdm (gdouble N_ur, gdouble N_ncdm)
{
  gdouble Neff = N_ur + N_ncdm * pow (ccl_constants.TNCDM, 4.0) / pow (4.0 / 11.0, 4.0 / 3.0);
  return Neff;
}

static ccl_cosmology *
test_nc_create_ccl_cosmo (const ccl_params_test_data *test_data)
{
  ccl_configuration config = default_config;
  gint status              = 0;
  ccl_parameters params;
  ccl_cosmology *cosmo;
  const gdouble Omega_k = test_data->use_Omega_k ? test_data->Omega_k : (1.0 - test_data->Omega_c - test_data->Omega_b - test_data->Omega_v);
  const gdouble Neff    = test_data->use_Neff_func ? neff_from_n_ur_n_ncdm (test_data->N_ur, test_data->N_ncdm) : test_data->Neff;
  
  params = ccl_parameters_create (test_data->Omega_c, test_data->Omega_b, Omega_k, 
                                  Neff, (gdouble *)test_data->mnu, ccl_mnu_list, test_data->w_0, 
                                  test_data->w_a, test_data->h, test_data->A_s, test_data->n_s, 
                                  -1, -1, -1, 0.0, 0.0, -1, NULL, NULL, &status);

	cosmo = ccl_cosmology_create (params, config);

  return cosmo;
}

#define MAX_Z (1200.0)

/* NOMNU */

void
test_nc_ccl_dist_new (TestNcCCLDist *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (pdata);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));
	test->dist      = nc_distance_new (MAX_Z);

  nc_distance_prepare (test->dist, test->cosmo);  
  
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
      
      /*printf ("(% 22.15g, % 22.15g) | OMEGA_X: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_x, cclE2Omega_x, E2Omega_x / cclE2Omega_x - 1.0);*/ 
    }

    if (TRUE)
    {
      const gdouble E2Omega_mnu    = nc_hicosmo_E2Omega_mnu (NC_HICOSMO (test->cosmo), z);
      const gdouble cclE2Omega_mnu = ccl_omega_x (test->ccl_cosmo, a, ccl_species_nu_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_mnu, ==, cclE2Omega_mnu, tol, 0.0);
      
      /*printf ("(% 22.15g, % 22.15g) | OMEGA_MNU: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_mnu, cclE2Omega_mnu, E2Omega_mnu / cclE2Omega_mnu - 1.0);*/
    }

    if (TRUE)
    {
      const gdouble E2Omega_m    = nc_hicosmo_E2Omega_m (test->cosmo, z) - nc_hicosmo_E2Omega_mnu (NC_HICOSMO (test->cosmo), z) + 3.0 * nc_hicosmo_E2Press_mnu (NC_HICOSMO (test->cosmo), z);
      const gdouble cclE2Omega_m = ccl_omega_x (test->ccl_cosmo, a, ccl_species_m_label, &status) * cclE2;

      ncm_assert_cmpdouble_e (E2Omega_m, ==, cclE2Omega_m, tol, 0.0);
      /*printf ("(% 22.15g, % 22.15g) | OMEGA_M: % 22.15g % 22.15g %17.10e;\n", z, a, E2Omega_m, cclE2Omega_m, E2Omega_m / cclE2Omega_m - 1.0);*/
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
test_nc_ccl_dist_traps (TestNcCCLDist *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/nc/ccl/dist/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_nc_ccl_dist_invalid_st (TestNcCCLDist *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}
