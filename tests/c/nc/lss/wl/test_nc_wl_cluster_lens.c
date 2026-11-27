/***************************************************************************
 *            test_nc_wl_cluster_lens.c
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_wl_cluster_lens.c
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
/*
 * Minor contributions:
 * Caio Lima de Oliveira <caiolimadeoliveira@pm.me> (2026)
 * - /nc/wl_cluster_lens/hash: the geometry hash does not follow the halo center.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcWLClusterLens
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcHaloMassSummary *hms;
  NcHaloDensityProfile *dp;
  NcHaloPosition *hp;
  NcWLSurfaceMassDensity *smd;
  NcWLClusterLens *lens;
} TestNcWLClusterLens;

void test_nc_wl_cluster_lens_new (TestNcWLClusterLens *test, gconstpointer pdata);
void test_nc_wl_cluster_lens_free (TestNcWLClusterLens *test, gconstpointer pdata);

void test_nc_wl_cluster_lens_exact (TestNcWLClusterLens *test, gconstpointer pdata);
void test_nc_wl_cluster_lens_hash (TestNcWLClusterLens *test, gconstpointer pdata);
void test_nc_wl_cluster_lens_serialize (TestNcWLClusterLens *test, gconstpointer pdata);

/* The curvatures of the fixtures: flat, open and closed */
static const gdouble test_omega_k[3] = {0.0, 0.05, -0.05};

#define TEST_Z_CL (0.4)

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/wl_cluster_lens/flat/exact", TestNcWLClusterLens, &test_omega_k[0],
              &test_nc_wl_cluster_lens_new,
              &test_nc_wl_cluster_lens_exact,
              &test_nc_wl_cluster_lens_free);

  g_test_add ("/nc/wl_cluster_lens/Okp/exact", TestNcWLClusterLens, &test_omega_k[1],
              &test_nc_wl_cluster_lens_new,
              &test_nc_wl_cluster_lens_exact,
              &test_nc_wl_cluster_lens_free);

  g_test_add ("/nc/wl_cluster_lens/Okn/exact", TestNcWLClusterLens, &test_omega_k[2],
              &test_nc_wl_cluster_lens_new,
              &test_nc_wl_cluster_lens_exact,
              &test_nc_wl_cluster_lens_free);

  g_test_add ("/nc/wl_cluster_lens/hash", TestNcWLClusterLens, &test_omega_k[0],
              &test_nc_wl_cluster_lens_new,
              &test_nc_wl_cluster_lens_hash,
              &test_nc_wl_cluster_lens_free);

  g_test_add ("/nc/wl_cluster_lens/serialize", TestNcWLClusterLens, &test_omega_k[0],
              &test_nc_wl_cluster_lens_new,
              &test_nc_wl_cluster_lens_serialize,
              &test_nc_wl_cluster_lens_free);

  g_test_run ();

  return 0;
}

void
test_nc_wl_cluster_lens_new (TestNcWLClusterLens *test, gconstpointer pdata)
{
  const gdouble Omega_k = *((const gdouble *) pdata);

  test->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  test->dist  = nc_distance_new (3.0);
  test->hms   = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  test->dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (test->hms));
  test->hp    = nc_halo_position_new (test->dist);
  test->smd   = nc_wl_surface_mass_density_new (test->dist);
  test->lens  = nc_wl_cluster_lens_new (test->dist);

  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", Omega_k, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 14.7, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "cDelta", 4.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "z", TEST_Z_CL, NULL);

  nc_wl_surface_mass_density_prepare (test->smd, test->cosmo);
  nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);
}

void
test_nc_wl_cluster_lens_free (TestNcWLClusterLens *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_wl_cluster_lens_free, test->lens);
  nc_wl_surface_mass_density_clear (&test->smd);
  nc_halo_position_clear (&test->hp);
  nc_halo_density_profile_clear (&test->dp);
  nc_halo_mass_summary_clear (&test->hms);
  nc_distance_clear (&test->dist);
  nc_hicosmo_clear (&test->cosmo);
}

/*
 * The critical surface density, the surface densities and the reduced shear of the lens
 * equal, bit for bit, those of the cached functions of NcWLSurfaceMassDensity (the lens
 * context, the critical and sigma caches, the optimization context and their reduced
 * shears) over radii and source redshifts, and the reduced shear is 0 in front of the
 * cluster.
 */
void
test_nc_wl_cluster_lens_exact (TestNcWLClusterLens *test, gconstpointer pdata)
{
  NcWLSurfaceMassDensityLensCtx ctx;
  gdouble r_s, rho_s;
  guint i, j;

  nc_wl_surface_mass_density_lens_ctx_prep (&ctx, test->smd, test->cosmo, TEST_Z_CL);
  nc_halo_density_profile_r_s_rho_s (test->dp, test->cosmo, TEST_Z_CL, &r_s, &rho_s);

  g_assert_cmpfloat (nc_wl_cluster_lens_get_z_cl (test->lens), ==, TEST_Z_CL);

  for (i = 0; i < 20; i++)
  {
    const gdouble R = 0.1 * exp (i * log (100.0) / 19.0);
    NcWLSurfaceMassDensitySigmaCache sigma_cache;
    NcWLSurfaceMassDensityOptzs optzs;
    gdouble sigma, mean_sigma;

    nc_wl_surface_mass_density_reduced_shear_sigma_cache_prep_with_rs (test->dp, R, r_s, rho_s, &sigma_cache);
    nc_wl_surface_mass_density_reduced_shear_optzs_prep_with_lens_ctx (test->dp, &ctx, R, r_s, rho_s, &optzs);
    nc_wl_cluster_lens_sigma (test->lens, test->dp, R, &sigma, &mean_sigma);

    g_assert_cmpfloat (sigma, ==, sigma_cache.sigma);
    g_assert_cmpfloat (mean_sigma, ==, sigma_cache.mean_sigma);

    for (j = 0; j < 30; j++)
    {
      const gdouble z_s = 0.05 + 0.1 * j;
      NcWLSurfaceMassDensityCritCache crit_cache;
      const gdouble g_optzs = nc_wl_surface_mass_density_reduced_shear_optzs (test->smd, test->dp, test->cosmo, z_s, TEST_Z_CL, &optzs);
      const gdouble g_lens  = nc_wl_cluster_lens_reduced_shear_at (test->lens, test->cosmo, sigma, mean_sigma, z_s);

      nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx (test->smd, test->cosmo, &ctx, z_s, &crit_cache);

      g_assert_cmpfloat (g_lens, ==, g_optzs);
      g_assert_cmpfloat (g_lens, ==, nc_wl_surface_mass_density_reduced_shear_cache (&crit_cache, &sigma_cache));

      if (z_s < TEST_Z_CL)
      {
        g_assert_false (crit_cache.is_source);
        g_assert_cmpfloat (g_lens, ==, 0.0);
      }
      else
      {
        const gdouble sigma_crit = nc_wl_cluster_lens_sigma_crit (test->lens, test->cosmo, z_s);

        g_assert_cmpfloat (sigma_crit, ==, crit_cache.sigma_crit);
        g_assert_cmpfloat (nc_wl_cluster_lens_reduced_shear (sigma, mean_sigma, sigma_crit), ==, g_lens);
      }
    }
  }
}

/*
 * The geometry hash changes with the cosmology and the cluster redshift, also when the
 * redshift is set back, not with the center; the profile hash changes with the mass, the
 * center, the redshift and the cosmology. A prepare after a change computes the values of
 * the new models.
 */
void
test_nc_wl_cluster_lens_hash (TestNcWLClusterLens *test, gconstpointer pdata)
{
  const guint64 geometry_hash0 = nc_wl_cluster_lens_get_geometry_hash (test->lens);
  const guint64 profile_hash0  = nc_wl_cluster_lens_get_profile_hash (test->lens);
  const gdouble z_cl0          = nc_wl_cluster_lens_get_z_cl (test->lens);
  gdouble sigma0, mean_sigma0, sigma, mean_sigma;

  nc_wl_cluster_lens_sigma (test->lens, test->dp, 1.0, &sigma0, &mean_sigma0);

  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", 1.0, NULL);
  nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);
  g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), ==, geometry_hash0);
  g_assert_cmpuint (nc_wl_cluster_lens_get_profile_hash (test->lens), !=, profile_hash0);

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 14.9, NULL);
  nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);
  g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), ==, geometry_hash0);
  nc_wl_cluster_lens_sigma (test->lens, test->dp, 1.0, &sigma, &mean_sigma);
  g_assert_cmpfloat (sigma, >, sigma0);
  g_assert_cmpfloat (mean_sigma, >, mean_sigma0);

  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "z", 0.5, NULL);
  nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);
  g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), !=, geometry_hash0);
  g_assert_cmpfloat (nc_wl_cluster_lens_get_z_cl (test->lens), ==, 0.5);
  g_assert_cmpfloat (nc_wl_cluster_lens_reduced_shear_at (test->lens, test->cosmo, sigma, mean_sigma, 0.45), ==, 0.0);

  {
    const guint64 geometry_hash_z = nc_wl_cluster_lens_get_geometry_hash (test->lens);

    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "z", z_cl0, NULL);
    nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);
    g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), !=, geometry_hash_z);
    g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), !=, geometry_hash0);
  }

  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "z", 0.5, NULL);
  nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);

  {
    const guint64 geometry_hash1 = nc_wl_cluster_lens_get_geometry_hash (test->lens);

    ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegac", 0.3, NULL);
    nc_wl_cluster_lens_prepare (test->lens, test->cosmo, test->hp, test->dp);
    g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), !=, geometry_hash1);
  }

  /* Another halo position object with the same parameters changes the profile, not the geometry */
  {
    NcHaloPosition *hp2          = nc_halo_position_new (test->dist);
    const guint64 geometry_hash2 = nc_wl_cluster_lens_get_geometry_hash (test->lens);
    const guint64 profile_hash2  = nc_wl_cluster_lens_get_profile_hash (test->lens);

    ncm_model_params_copyto (NCM_MODEL (test->hp), NCM_MODEL (hp2));
    nc_wl_cluster_lens_prepare (test->lens, test->cosmo, hp2, test->dp);
    g_assert_cmpuint (nc_wl_cluster_lens_get_geometry_hash (test->lens), ==, geometry_hash2);
    g_assert_cmpuint (nc_wl_cluster_lens_get_profile_hash (test->lens), !=, profile_hash2);

    nc_halo_position_free (hp2);
  }
}

/* A lens read back from its serialization holds a distance and gives the same values */
void
test_nc_wl_cluster_lens_serialize (TestNcWLClusterLens *test, gconstpointer pdata)
{
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcWLClusterLens *dup = NC_WL_CLUSTER_LENS (ncm_serialize_dup_obj (ser, G_OBJECT (test->lens)));
  gdouble sigma, mean_sigma;

  g_assert_true (NC_IS_WL_CLUSTER_LENS (dup));
  g_assert_true (NC_IS_DISTANCE (nc_wl_cluster_lens_peek_distance (dup)));

  nc_wl_cluster_lens_prepare (dup, test->cosmo, test->hp, test->dp);
  nc_wl_cluster_lens_sigma (dup, test->dp, 1.0, &sigma, &mean_sigma);

  g_assert_cmpfloat (nc_wl_cluster_lens_reduced_shear_at (dup, test->cosmo, sigma, mean_sigma, 1.2), ==,
                     nc_wl_cluster_lens_reduced_shear_at (test->lens, test->cosmo, sigma, mean_sigma, 1.2));

  ncm_serialize_free (ser);
  NCM_TEST_FREE (nc_wl_cluster_lens_free, dup);
}

