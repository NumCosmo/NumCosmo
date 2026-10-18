/***************************************************************************
 *            test_nc_hcosmo_bi_dist.c
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_hcosmo_bi_dist.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_integration.h>

typedef struct _TestNcHCosmoBIDist
{
  NcHCosmoBI *cosmo;
  NcHCosmoBIDist *dist;
} TestNcHCosmoBIDist;

typedef struct _TestNcHCosmoBIDistK
{
  NcHCosmoBI *cosmo;
  gdouble n2[3];
  guint i;
  guint j;
} TestNcHCosmoBIDistK;

void test_nc_hcosmo_bi_dist_new (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_free (TestNcHCosmoBIDist *test, gconstpointer pdata);

void test_nc_hcosmo_bi_dist_lcdm (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_quad (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_rotation (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_data (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_properties (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_zmax (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_zmax_unreached (TestNcHCosmoBIDist *test, gconstpointer pdata);
void test_nc_hcosmo_bi_dist_pitfalls (TestNcHCosmoBIDist *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/hcosmo_bi_dist/lcdm", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_lcdm,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/quad", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_quad,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/rotation", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_rotation,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/data", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_data,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/properties", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_properties,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/pitfalls", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_pitfalls,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/zmax/subprocess", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_zmax,
              &test_nc_hcosmo_bi_dist_free);

  g_test_add ("/nc/hcosmo_bi_dist/zmax_unreached/subprocess", TestNcHCosmoBIDist, NULL,
              &test_nc_hcosmo_bi_dist_new,
              &test_nc_hcosmo_bi_dist_zmax_unreached,
              &test_nc_hcosmo_bi_dist_free);

  g_test_run ();
}

static void
_test_nc_hcosmo_bi_dist_set (NcHCosmoBI *cosmo, const gdouble Omega_m0, const gdouble H1, const gdouble H2, const gdouble H3)
{
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_H1, H1);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_H2, H2);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_H3, H3);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HCOSMO_BI_OMEGA_M, Omega_m0);
}

/* Integrand of the flat LCDM comoving distance in lookback time, 1 / a(t_lb); the QAG references below are at reltol 1e-13 */
static gdouble
_test_nc_hcosmo_bi_dist_inv_a (gdouble t_lb, gpointer params)
{
  const gdouble *p = params; /* Omega_m0, Omega_L0, age */
  const gdouble r  = sqrt (p[1]);

  return 1.0 / (cbrt (p[0] / p[1]) * pow (sinh (1.5 * r * (p[2] - t_lb)), 2.0 / 3.0));
}

/* Integrand u_i u_j / E^3 of K_ij along a direction, from the closed form of the model */
static gdouble
_test_nc_hcosmo_bi_dist_K_integrand (gdouble t_lb, gpointer params)
{
  TestNcHCosmoBIDistK *K = params;
  gdouble ln_a[3], u[3], E2;

  nc_hcosmo_bi_ln_a (K->cosmo, t_lb, &ln_a[0], &ln_a[1], &ln_a[2]);
  u[0] = exp (-2.0 * ln_a[0]);
  u[1] = exp (-2.0 * ln_a[1]);
  u[2] = exp (-2.0 * ln_a[2]);
  E2   = K->n2[0] * u[0] + K->n2[1] * u[1] + K->n2[2] * u[2];

  return u[K->i] * u[K->j] / (E2 * sqrt (E2));
}

void
test_nc_hcosmo_bi_dist_new (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  test->cosmo = nc_hcosmo_bi_new ();
  test->dist  = nc_hcosmo_bi_dist_new (0.3);

  g_assert_true (NC_IS_HCOSMO_BI (test->cosmo));
  g_assert_true (NC_IS_HCOSMO_BI_DIST (test->dist));
}

void
test_nc_hcosmo_bi_dist_free (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hcosmo_bi_dist_free, test->dist);
  NCM_TEST_FREE (nc_hcosmo_bi_free, test->cosmo);
}

void
test_nc_hcosmo_bi_dist_lcdm (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  const gdouble Omega_m0        = 0.3;
  const gdouble Omega_L0        = 0.7;
  const gdouble r               = sqrt (Omega_L0);
  const gdouble age             = asinh (sqrt (Omega_L0 / Omega_m0)) / (1.5 * r);
  gdouble p[3]                  = {Omega_m0, Omega_L0, age};
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  guint j;

  F.function = &_test_nc_hcosmo_bi_dist_inv_a;
  F.params   = p;

  _test_nc_hcosmo_bi_dist_set (test->cosmo, Omega_m0, 70.0, 70.0, 70.0);
  nc_hcosmo_bi_dist_prepare (test->dist, test->cosmo);

  for (j = 1; j <= 10; j++)
  {
    const gdouble z   = 0.3 * j / 10.0;
    const gdouble a_e = 1.0 / (1.0 + z);
    const gdouble t_e = age - asinh (sqrt (Omega_L0 / Omega_m0) * pow (a_e, 1.5)) / (1.5 * r);
    gdouble chi, err;

    gsl_integration_qag (&F, 0.0, t_e, 0.0, 1.0e-13, 1000, GSL_INTEG_GAUSS61, ws, &chi, &err);

    ncm_assert_cmpdouble_e (nc_hcosmo_bi_dist_t_lb_emission (test->dist, z, 0.6, 0.0, 0.8), ==, t_e, 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (nc_hcosmo_bi_dist_luminosity (test->dist, z, 0.6, 0.0, 0.8), ==, (1.0 + z) * chi, 1.0e-13, 0.0);
  }

  gsl_integration_workspace_free (ws);
}

void
test_nc_hcosmo_bi_dist_quad (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  const gdouble dirs[3][3] = {
    {0.6, 0.0, 0.8},
    {0.48, 0.6, 0.64},
    {1.0, 0.0, 0.0},
  };
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (1000);
  TestNcHCosmoBIDistK K;
  gsl_function F;
  guint d, j;

  F.function = &_test_nc_hcosmo_bi_dist_K_integrand;
  F.params   = &K;
  K.cosmo    = test->cosmo;

  _test_nc_hcosmo_bi_dist_set (test->cosmo, 0.334, 0.85 * 73.04, 73.04, 1.15 * 73.04);
  nc_hcosmo_bi_dist_prepare (test->dist, test->cosmo);

  for (d = 0; d < 3; d++)
  {
    K.n2[0] = gsl_pow_2 (dirs[d][0]);
    K.n2[1] = gsl_pow_2 (dirs[d][1]);
    K.n2[2] = gsl_pow_2 (dirs[d][2]);

    for (j = 1; j <= 6; j++)
    {
      const gdouble z      = 0.3 * j / 6.0;
      const gdouble t_e    = nc_hcosmo_bi_dist_t_lb_emission (test->dist, z, dirs[d][0], dirs[d][1], dirs[d][2]);
      const guint ij[3][2] = {
        {
          0, 1
        }, {
          0, 2
        }, {
          1, 2
        }
      };
      gdouble Kij[3], ln_a[3], E2, err, DA2;
      guint l;

      nc_hcosmo_bi_ln_a (test->cosmo, t_e, &ln_a[0], &ln_a[1], &ln_a[2]);
      E2 = K.n2[0] * exp (-2.0 * ln_a[0]) + K.n2[1] * exp (-2.0 * ln_a[1]) + K.n2[2] * exp (-2.0 * ln_a[2]);
      ncm_assert_cmpdouble_e (sqrt (E2), ==, 1.0 + z, 1.0e-14, 0.0);

      for (l = 0; l < 3; l++)
      {
        K.i = ij[l][0];
        K.j = ij[l][1];
        gsl_integration_qag (&F, 0.0, t_e, 0.0, 1.0e-13, 1000, GSL_INTEG_GAUSS61, ws, &Kij[l], &err);
      }

      DA2 = exp (ln_a[0] + ln_a[1] + ln_a[2]) * (1.0 + z) *
            (K.n2[0] * Kij[0] * Kij[1] + K.n2[1] * Kij[0] * Kij[2] + K.n2[2] * Kij[1] * Kij[2]);

      ncm_assert_cmpdouble_e (nc_hcosmo_bi_dist_luminosity (test->dist, z, dirs[d][0], dirs[d][1], dirs[d][2]), ==,
                              gsl_pow_2 (1.0 + z) * sqrt (DA2), 1.0e-13, 0.0);
    }
  }

  gsl_integration_workspace_free (ws);
}

void
test_nc_hcosmo_bi_dist_rotation (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  gdouble n1, n2, n3;

  _test_nc_hcosmo_bi_dist_set (test->cosmo, 0.334, 0.85 * 73.04, 73.04, 1.15 * 73.04);

  /* Default orientation is the identity */
  nc_hcosmo_bi_dist_prepare (test->dist, test->cosmo);
  nc_hcosmo_bi_dist_sky_to_axes (test->dist, 0.48, 0.6, 0.64, &n1, &n2, &n3);
  ncm_assert_cmpdouble_e (n1, ==, 0.48, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (n2, ==, 0.6, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (n3, ==, 0.64, 1.0e-15, 0.0);

  /* u = (1, 0, 1/8): q = (cos(pi/4); 0, 0, sin(pi/4)), rotation by pi/2 around the z axis */
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HCOSMO_BI_U3, 0.125);
  nc_hcosmo_bi_dist_prepare_if_needed (test->dist, test->cosmo);
  nc_hcosmo_bi_dist_sky_to_axes (test->dist, 1.0, 0.0, 0.0, &n1, &n2, &n3);
  g_assert_cmpfloat (fabs (n1), <, 1.0e-15);
  ncm_assert_cmpdouble_e (n2, ==, 1.0, 1.0e-15, 0.0);
  g_assert_cmpfloat (fabs (n3), <, 1.0e-15);

  /* The sky distance modulus is the axes one at the rotated direction */
  nc_hcosmo_bi_dist_sky_to_axes (test->dist, 0.48, 0.6, 0.64, &n1, &n2, &n3);
  ncm_assert_cmpdouble_e (nc_hcosmo_bi_dist_dmodulus_sky (test->dist, 0.1, 0.48, 0.6, 0.64), ==,
                          5.0 * log10 (nc_hcosmo_bi_dist_luminosity (test->dist, 0.1, n1, n2, n3)) + 25.0, 1.0e-15, 0.0);
}

void
test_nc_hcosmo_bi_dist_data (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  NcDataHCosmoBIDistMu *dmu = nc_data_hcosmo_bi_dist_mu_new (test->dist);
  NcmMSet *mset             = ncm_mset_new (NCM_MODEL (test->cosmo), NULL, NULL);
  NcmVector *z              = ncm_vector_new (3);
  NcmVector *mu             = ncm_vector_new (3);
  NcmVector *sigma          = ncm_vector_new (3);
  NcmMatrix *p              = ncm_matrix_new (3, 3);
  const gdouble dirs[3][3]  = {
    {
      0.6, 0.0, 0.8
    }, {
      0.48, 0.6, 0.64
    }, {
      0.0, 1.0, 0.0
    }
  };
  gdouble m2lnL_ref = 0.0;
  gdouble term[3];
  gdouble m2lnL;
  guint i;

  _test_nc_hcosmo_bi_dist_set (test->cosmo, 0.3, 70.0, 72.0, 75.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HCOSMO_BI_U1, 0.3);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HCOSMO_BI_U2, 0.2);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HCOSMO_BI_U3, 0.7);
  nc_hcosmo_bi_dist_prepare (test->dist, test->cosmo);

  for (i = 0; i < 3; i++)
  {
    const gdouble zi   = 0.01 * (i + 1);
    const gdouble mu_t = nc_hcosmo_bi_dist_dmodulus_sky (test->dist, zi, dirs[i][0], dirs[i][1], dirs[i][2]) +
                         5.0 * log10 (nc_hcosmo_bi_RH_Mpc (test->cosmo));
    const gdouble mu_i = mu_t + 0.05 * (i + 1.0);
    const gdouble s_i  = 0.1 + 0.05 * i;

    ncm_vector_set (z, i, zi);
    ncm_vector_set (mu, i, mu_i);
    ncm_vector_set (sigma, i, s_i);
    ncm_matrix_set (p, i, 0, dirs[i][0]);
    ncm_matrix_set (p, i, 1, dirs[i][1]);
    ncm_matrix_set (p, i, 2, dirs[i][2]);

    term[i]    = gsl_pow_2 ((mu_i - mu_t) / s_i) + 2.0 * log (s_i) + ncm_c_ln2pi ();
    m2lnL_ref += term[i];
  }

  nc_data_hcosmo_bi_dist_mu_set_data (dmu, z, p, mu, sigma);
  g_assert_true (ncm_data_is_init (NCM_DATA (dmu)));
  g_assert_true (nc_data_hcosmo_bi_dist_mu_peek_z (dmu) != z);
  ncm_assert_cmpdouble (ncm_vector_get (nc_data_hcosmo_bi_dist_mu_peek_z (dmu), 2), ==, 0.03);
  ncm_assert_cmpdouble (ncm_matrix_get (nc_data_hcosmo_bi_dist_mu_peek_p (dmu), 1, 1), ==, 0.6);

  ncm_data_prepare (NCM_DATA (dmu), mset);
  ncm_data_m2lnL_val (NCM_DATA (dmu), mset, &m2lnL);
  ncm_assert_cmpdouble_e (m2lnL, ==, m2lnL_ref, 1.0e-13, 0.0);

  /* The catalog built from the data, read back with a cut keeping the first two objects */
  {
    NcDistMuObs *obs              = nc_data_hcosmo_bi_dist_mu_get_obs (dmu);
    NcDataHCosmoBIDistMu *dmu_obs = nc_data_hcosmo_bi_dist_mu_new (test->dist);
    NcmVector *z_obs;
    gdouble m2lnL_obs, m2lnL_two;
    GError *error = NULL;

    g_assert_cmpuint (ncm_catalog_len (NCM_CATALOG (obs)), ==, 3);
    ncm_assert_cmpdouble_e (ncm_catalog_get (NCM_CATALOG (obs), NC_DIST_MU_OBS_COL_MU, 1, NULL), ==, ncm_vector_get (mu, 1), 0.0, 0.0);

    nc_data_hcosmo_bi_dist_mu_set_obs (dmu_obs, obs, NULL, NULL, NULL, 0.0, 1.0, &error);
    g_assert_no_error (error);
    ncm_data_prepare (NCM_DATA (dmu_obs), mset);
    ncm_data_m2lnL_val (NCM_DATA (dmu_obs), mset, &m2lnL_obs);
    ncm_assert_cmpdouble_e (m2lnL_obs, ==, m2lnL, 1.0e-13, 0.0);

    nc_data_hcosmo_bi_dist_mu_set_obs (dmu_obs, obs, NULL, NULL, NULL, 0.0, 0.025, &error);
    g_assert_no_error (error);
    z_obs = nc_data_hcosmo_bi_dist_mu_peek_z (dmu_obs);
    g_assert_cmpuint (ncm_vector_len (z_obs), ==, 2);
    ncm_data_prepare (NCM_DATA (dmu_obs), mset);
    ncm_data_m2lnL_val (NCM_DATA (dmu_obs), mset, &m2lnL_two);
    ncm_assert_cmpdouble_e (m2lnL_two, ==, term[0] + term[1], 1.0e-13, 0.0);

    nc_data_hcosmo_bi_dist_mu_set_obs (dmu_obs, obs, "DM_av", NULL, NULL, 0.0, 1.0, &error);
    g_assert_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND);
    g_clear_error (&error);

    nc_data_hcosmo_bi_dist_mu_set_obs (dmu_obs, obs, NULL, "eDM_av", NULL, 0.0, 1.0, &error);
    g_assert_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND);
    g_clear_error (&error);

    /* Offsets and the velocity term through a NcDistMuNuisance in the model set */
    {
      NcDistMuNuisance *nuisance = nc_dist_mu_nuisance_new (2);
      NcmMSet *mset_n            = ncm_mset_new (NCM_MODEL (test->cosmo), NULL, nuisance, NULL);
      NcmVector *sample          = ncm_vector_new (3);
      const gdouble dmu_s[2]     = {0.1, -0.2};
      const gdouble sigma_v      = 300.0;
      gdouble m2lnL_n, m2lnL_ref_n = 0.0, m2lnL_again;

      ncm_vector_set (sample, 0, 0.0);
      ncm_vector_set (sample, 1, 1.0);
      ncm_vector_set (sample, 2, -1.0);
      nc_data_hcosmo_bi_dist_mu_set_samples (dmu, sample);
      ncm_assert_cmpdouble (ncm_vector_get (nc_data_hcosmo_bi_dist_mu_peek_samples (dmu), 1), ==, 1.0);

      ncm_model_orig_param_set (NCM_MODEL (nuisance), NC_DIST_MU_NUISANCE_SIGMA_V, sigma_v);
      ncm_model_orig_vparam_set (NCM_MODEL (nuisance), NC_DIST_MU_NUISANCE_DMU, 0, dmu_s[0]);
      ncm_model_orig_vparam_set (NCM_MODEL (nuisance), NC_DIST_MU_NUISANCE_DMU, 1, dmu_s[1]);
      g_assert_cmpuint (nc_dist_mu_nuisance_get_nsamples (nuisance), ==, 2);

      for (i = 0; i < 3; i++)
      {
        const gdouble zi    = ncm_vector_get (z, i);
        const gdouble mu_t  = ncm_vector_get (mu, i) - 0.05 * (i + 1.0) + ((i < 2) ? dmu_s[i] : 0.0);
        const gdouble s_v   = 5.0 / M_LN10 * (sigma_v / (1.0e-3 * ncm_c_c ())) * (1.0 + zi) / (zi * (1.0 + 0.5 * zi));
        const gdouble s_tot = hypot (ncm_vector_get (sigma, i), s_v);

        ncm_assert_cmpdouble_e (nc_dist_mu_nuisance_sigma_mu_v (nuisance, zi), ==, s_v, 1.0e-15, 0.0);
        m2lnL_ref_n += gsl_pow_2 ((ncm_vector_get (mu, i) - mu_t) / s_tot) + 2.0 * log (s_tot) + ncm_c_ln2pi ();
      }

      ncm_data_prepare (NCM_DATA (dmu), mset_n);
      ncm_data_m2lnL_val (NCM_DATA (dmu), mset_n, &m2lnL_n);
      ncm_assert_cmpdouble_e (m2lnL_n, ==, m2lnL_ref_n, 1.0e-13, 0.0);

      /* Unchanged sigma_v: the same value; everything back to zero: the original likelihood */
      ncm_data_m2lnL_val (NCM_DATA (dmu), mset_n, &m2lnL_again);
      ncm_assert_cmpdouble_e (m2lnL_again, ==, m2lnL_n, 0.0, 0.0);
      ncm_model_orig_param_set (NCM_MODEL (nuisance), NC_DIST_MU_NUISANCE_SIGMA_V, 0.0);
      ncm_model_orig_vparam_set (NCM_MODEL (nuisance), NC_DIST_MU_NUISANCE_DMU, 0, 0.0);
      ncm_model_orig_vparam_set (NCM_MODEL (nuisance), NC_DIST_MU_NUISANCE_DMU, 1, 0.0);
      ncm_data_m2lnL_val (NCM_DATA (dmu), mset_n, &m2lnL_again);
      ncm_assert_cmpdouble_e (m2lnL_again, ==, m2lnL, 1.0e-13, 0.0);

      /* The catalog written from the data carries the sample column */
      {
        NcDistMuObs *obs_s = nc_data_hcosmo_bi_dist_mu_get_obs (dmu);

        g_assert_cmpint (ncm_catalog_get_int (NCM_CATALOG (obs_s), NC_DIST_MU_OBS_COL_SAMPLE, 1, NULL), ==, 1);
        g_assert_cmpint (ncm_catalog_get_int (NCM_CATALOG (obs_s), NC_DIST_MU_OBS_COL_SAMPLE, 2, NULL), ==, -1);

        nc_data_hcosmo_bi_dist_mu_set_obs (dmu_obs, obs_s, NULL, NULL, NC_DIST_MU_OBS_COL_SAMPLE, 0.0, 1.0, &error);
        g_assert_no_error (error);
        ncm_assert_cmpdouble (ncm_vector_get (nc_data_hcosmo_bi_dist_mu_peek_samples (dmu_obs), 0), ==, 0.0);
        ncm_assert_cmpdouble (ncm_vector_get (nc_data_hcosmo_bi_dist_mu_peek_samples (dmu_obs), 2), ==, -1.0);

        nc_dist_mu_obs_free (obs_s);
      }

      ncm_vector_free (sample);
      ncm_mset_free (mset_n);
      nc_dist_mu_nuisance_free (nuisance);
    }

    nc_data_hcosmo_bi_dist_mu_free (dmu_obs);
    nc_dist_mu_obs_free (obs);
  }

  ncm_vector_free (z);
  ncm_vector_free (mu);
  ncm_vector_free (sigma);
  ncm_matrix_free (p);
  ncm_mset_free (mset);
  NCM_TEST_FREE (nc_data_hcosmo_bi_dist_mu_free, dmu);
}

void
test_nc_hcosmo_bi_dist_properties (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcHCosmoBIDist *dup;
  NcDataHCosmoBIDistMu *dmu;
  NcDataHCosmoBIDistMu *dmu_dup;
  gdouble dl_a, dl_b;

  ncm_assert_cmpdouble (nc_hcosmo_bi_dist_get_zmax (test->dist), ==, 0.3);
  ncm_assert_cmpdouble (nc_hcosmo_bi_dist_get_reltol (test->dist), ==, NC_HCOSMO_BI_DIST_DEFAULT_RELTOL);

  _test_nc_hcosmo_bi_dist_set (test->cosmo, 0.334, 0.85 * 73.04, 73.04, 1.15 * 73.04);
  nc_hcosmo_bi_dist_prepare_if_needed (test->dist, test->cosmo);
  dl_a = nc_hcosmo_bi_dist_luminosity (test->dist, 0.2, 0.6, 0.0, 0.8);

  /* A looser tolerance changes the node count and the result stays within it */
  nc_hcosmo_bi_dist_set_reltol (test->dist, 1.0e-6);
  nc_hcosmo_bi_dist_prepare_if_needed (test->dist, test->cosmo);
  g_assert_cmpuint (nc_hcosmo_bi_dist_get_n_nodes (test->dist), <=, 33);
  dl_b = nc_hcosmo_bi_dist_luminosity (test->dist, 0.2, 0.6, 0.0, 0.8);
  ncm_assert_cmpdouble_e (dl_b, ==, dl_a, 1.0e-6, 0.0);

  nc_hcosmo_bi_dist_set_zmax (test->dist, 0.5);
  nc_hcosmo_bi_dist_prepare_if_needed (test->dist, test->cosmo);
  g_assert_cmpfloat (nc_hcosmo_bi_dist_get_t_lb_max (test->dist), >, 0.0);

  /* Serialization keeps zmax and reltol; the data object keeps its calculator and tables */
  dup = NC_HCOSMO_BI_DIST (ncm_serialize_dup_obj (ser, G_OBJECT (test->dist)));
  ncm_assert_cmpdouble (nc_hcosmo_bi_dist_get_zmax (dup), ==, 0.5);
  ncm_assert_cmpdouble (nc_hcosmo_bi_dist_get_reltol (dup), ==, 1.0e-6);

  dmu = nc_data_hcosmo_bi_dist_mu_new (NULL);
  nc_data_hcosmo_bi_dist_mu_set_dist (dmu, dup);
  {
    NcmVector *z     = ncm_vector_new (2);
    NcmVector *mu    = ncm_vector_new (2);
    NcmVector *sigma = ncm_vector_new (2);
    NcmMatrix *p     = ncm_matrix_new (2, 3);

    ncm_vector_set_all (z, 0.1);
    ncm_vector_set_all (mu, 38.0);
    ncm_vector_set_all (sigma, 0.2);
    ncm_matrix_set_identity (p);
    nc_data_hcosmo_bi_dist_mu_set_data (dmu, z, p, mu, sigma);
    g_assert_cmpuint (ncm_matrix_nrows (nc_data_hcosmo_bi_dist_mu_peek_p (dmu)), ==, 2);

    ncm_vector_free (z);
    ncm_vector_free (mu);
    ncm_vector_free (sigma);
    ncm_matrix_free (p);
  }
  dmu_dup = NC_DATA_HCOSMO_BI_DIST_MU (ncm_serialize_dup_obj (ser, G_OBJECT (dmu)));
  g_assert_cmpuint (ncm_vector_len (nc_data_hcosmo_bi_dist_mu_peek_z (dmu_dup)), ==, 2);
  ncm_assert_cmpdouble (ncm_matrix_get (nc_data_hcosmo_bi_dist_mu_peek_p (dmu_dup), 1, 1), ==, 1.0);
  g_assert_true (ncm_data_is_init (NCM_DATA (dmu_dup)));

  nc_data_hcosmo_bi_dist_mu_free (dmu_dup);
  nc_data_hcosmo_bi_dist_mu_free (dmu);
  nc_hcosmo_bi_dist_free (dup);
  ncm_serialize_free (ser);
}

void
test_nc_hcosmo_bi_dist_pitfalls (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/nc/hcosmo_bi_dist/zmax/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*redshift*above zmax*");

  g_test_trap_subprocess ("/nc/hcosmo_bi_dist/zmax_unreached/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*an axis contracts*");
}

void
test_nc_hcosmo_bi_dist_zmax (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  _test_nc_hcosmo_bi_dist_set (test->cosmo, 0.3, 70.0, 70.0, 70.0);
  nc_hcosmo_bi_dist_prepare (test->dist, test->cosmo);
  nc_hcosmo_bi_dist_luminosity (test->dist, 0.31, 1.0, 0.0, 0.0);
}

void
test_nc_hcosmo_bi_dist_zmax_unreached (TestNcHCosmoBIDist *test, gconstpointer pdata)
{
  /* With this shear axis 1 contracts before 1 + z = 2 is reached along it */
  _test_nc_hcosmo_bi_dist_set (test->cosmo, 0.334, 0.85 * 73.04, 73.04, 1.15 * 73.04);
  nc_hcosmo_bi_dist_set_zmax (test->dist, 1.0);
  nc_hcosmo_bi_dist_prepare (test->dist, test->cosmo);
}

