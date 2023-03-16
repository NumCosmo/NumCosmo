/***************************************************************************
 *            test_ncm_stats_dist1d_epdf.c
 *
 *  Mon February 27 11:01:27 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br> & <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2023 <vitenti@uel.br>
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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>


typedef struct _TestNcGalaxyWLEllipticityBinned
{
  NcGalaxyWLEllipticityBinned *gebin;
  NcHICosmo *cosmo;
  NcHaloDensityProfile *dp;
  NcWLSurfaceMassDensity *smd;
  NcGalaxyRedshift *gz;
  gdouble z_cluster;
} TestNcGalaxyWLEllipticityBinned;

static void test_nc_galaxy_wl_ellipticity_binned_new (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_binned_rep (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_binned_n_bins(TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_binned_binning (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_binned_free (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  // g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/galaxy_wl_ellipticity_binned/binned_rep", TestNcGalaxyWLEllipticityBinned, NULL,
              &test_nc_galaxy_wl_ellipticity_binned_new,
              &test_nc_galaxy_wl_ellipticity_binned_rep,
              &test_nc_galaxy_wl_ellipticity_binned_free);

  g_test_add ("/nc/galaxy_wl_ellipticity_binned/n_bins", TestNcGalaxyWLEllipticityBinned, NULL,
              &test_nc_galaxy_wl_ellipticity_binned_new,
              &test_nc_galaxy_wl_ellipticity_binned_n_bins,
              &test_nc_galaxy_wl_ellipticity_binned_free);

  g_test_add ("/nc/galaxy_wl_ellipticity_binned/binning", TestNcGalaxyWLEllipticityBinned, NULL,
              &test_nc_galaxy_wl_ellipticity_binned_new,
              &test_nc_galaxy_wl_ellipticity_binned_binning,
              &test_nc_galaxy_wl_ellipticity_binned_free);

  g_test_run ();
}

static void 
test_nc_galaxy_wl_ellipticity_binned_new (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata)
{
  // NcGalaxyWLEllipticityBinned *gebin = nc_galaxy_wl_ellipticity_binned_new ();
  NcHICosmo *cosmo                   = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoLCDM");
  NcHaloDensityProfile *dp           = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL, 200.0));
  NcDistance *dist                   = nc_distance_new (3.0);
  NcWLSurfaceMassDensity *smd        = nc_wl_surface_mass_density_new (dist);
  NcmRNG *rng                        = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const gdouble n                    = 10000;
  NcmMatrix *data                    = ncm_matrix_new (n, 3);
  NcmVector *z_vec                   = ncm_vector_new (n);
  NcGalaxyRedshiftSpec *gzs          = nc_galaxy_redshift_spec_new ();
  const gdouble rl                   = g_test_rand_double_range (0.0002, 0.1300);
  const gdouble ru                   = g_test_rand_double_range (5.55, 5.65);
  const gint bin_n                   = g_test_rand_int_range (10, 100);
  NcmVector *bin_vector              = ncm_vector_new (bin_n+1);
  const gdouble bin_size             = (ru - rl)/bin_n;
  const gdouble mu_e                 = g_test_rand_double_range (0.015, 0.020);
  const gdouble sigma_e              = g_test_rand_double_range (0.050, 0.055);
  const gdouble e_noise              = g_test_rand_double_range (0.02, 0.08);
  const gdouble mu_z                 = g_test_rand_double_range (1.30, 1.35);
  const gdouble sigma_z              = g_test_rand_double_range (0.65, 0.75);
  const gdouble z_cluster            = g_test_rand_double_range (0.2, 0.4);
  gint i;
  gdouble j;

  for (j = 0; j < bin_n+1; j++)
  {
    ncm_vector_set (bin_vector, j, rl + j*bin_size);
  }

  for (i = 0; i < n; i++)
  {
    const gdouble r = g_test_rand_double_range (rl, ru);
    const gdouble e = ncm_rng_gaussian_gen (rng, mu_e, sigma_e);
    const gdouble z = ncm_rng_gaussian_gen (rng, mu_z, sigma_z);

    ncm_matrix_set (data, i, 0, r);
    ncm_matrix_set (data, i, 1, e);
    ncm_matrix_set (data, i, 2, e_noise);

    ncm_vector_set (z_vec, i, z);
  }

  nc_galaxy_redshift_spec_set_z (gzs, z_vec);
  nc_distance_prepare (dist, cosmo);

  test->gebin     = nc_galaxy_wl_ellipticity_binned_new ();
  test->cosmo     = cosmo;
  test->dp        = dp;
  test->smd       = smd;
  test->gz        = NC_GALAXY_REDSHIFT (gzs);
  test->z_cluster = z_cluster;

  nc_galaxy_wl_ellipticity_binned_set_binobs (test->gebin, data, bin_vector);

  g_assert_true (NC_IS_GALAXY_WL_ELLIPTICITY_BINNED (test->gebin));
}

static void
test_nc_galaxy_wl_ellipticity_binned_rep (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata)
{
  NcmObjArray *bin_obs = nc_galaxy_wl_ellipticity_binned_peek_binobs (test->gebin);
  NcmVector *bins = nc_galaxy_wl_ellipticity_binned_peek_bins (test->gebin);
  gint i;


  for (i = 0; i < ncm_vector_len (bins)-1; i++)
  {
    NcmMatrix *bin_data_i = NCM_MATRIX (ncm_obj_array_get (bin_obs, i));
    gint j;


    for (j = 0; j < ncm_matrix_nrows (bin_data_i); j++)
    {
      gdouble gal_j[3] = {ncm_matrix_get (bin_data_i, j, 0), ncm_matrix_get (bin_data_i, j, 1), ncm_matrix_get (bin_data_i, j, 2)};
      gint k;

      for (k = 0; k < ncm_vector_len (bins)-1; k++)
      {
        if (k != i)
        {
          NcmMatrix *bin_data_k = NCM_MATRIX (ncm_obj_array_get (bin_obs, k));
          gint l;

          for (l = 0; l < ncm_matrix_nrows (bin_data_k); l++)
          {
            gdouble gal_l[3] = {ncm_matrix_get (bin_data_k, l, 0), ncm_matrix_get (bin_data_k, l, 1), ncm_matrix_get (bin_data_k, l, 2)};

            g_assert_true (!((gal_j[0] == gal_l[0]) && (gal_j[1] == gal_l[1]) && (gal_j[2] == gal_l[2])));
          }
        }
      }
    }
  }
}

static void
test_nc_galaxy_wl_ellipticity_binned_n_bins (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata)
{
  NcmObjArray *bin_obs = nc_galaxy_wl_ellipticity_binned_peek_binobs (test->gebin);
  NcmVector *bins = nc_galaxy_wl_ellipticity_binned_peek_bins (test->gebin);

  g_assert_cmpfloat (ncm_obj_array_len (bin_obs) + 1, ==, ncm_vector_len (bins));
}

static void
test_nc_galaxy_wl_ellipticity_binned_binning (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata)
{
  NcmObjArray *bin_obs = nc_galaxy_wl_ellipticity_binned_peek_binobs (test->gebin);
  NcmVector *bins      = nc_galaxy_wl_ellipticity_binned_peek_bins (test->gebin);
  gint bin_i;

  for (bin_i = 0; bin_i < ncm_obj_array_len (bin_obs); bin_i++)
  {
    NcmMatrix *bin_data = NCM_MATRIX (ncm_obj_array_get (bin_obs, bin_i));
    gint j;
    printf ("%d\n", ncm_matrix_nrows (bin_data));

    for (j = 0; j < ncm_matrix_nrows (bin_data); j++)
    {
      gdouble gal_r = ncm_matrix_get (bin_data, j, 0);

      if (bin_i == ncm_obj_array_len (bin_obs)-1)
      {
        g_assert_true ((gal_r >= ncm_vector_get (bins, bin_i)) && (gal_r <= ncm_vector_get (bins, bin_i+1)));
      }
      else
      {
        g_assert_true ((gal_r >= ncm_vector_get (bins, bin_i)) && (gal_r < ncm_vector_get (bins, bin_i+1)));
      }
    }
  }
}

static void
test_nc_galaxy_wl_ellipticity_binned_free (TestNcGalaxyWLEllipticityBinned *test, gconstpointer pdata)
{
  g_assert_true (NC_IS_GALAXY_WL_ELLIPTICITY_BINNED (test->gebin));
  NCM_TEST_FREE (nc_galaxy_wl_ellipticity_binned_free, NC_GALAXY_WL_ELLIPTICITY_BINNED(test->gebin));
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
  NCM_TEST_FREE (nc_halo_density_profile_free, test->dp);
  NCM_TEST_FREE (nc_wl_surface_mass_density_free, test->smd);
  NCM_TEST_FREE (nc_galaxy_redshift_free, test->gz);


  // ncm_model_free (NCM_MODEL (test->cosmo));
  // ncm_model_free (NCM_MODEL (test->dp));

  // nc_cluster_redshift_free (NC_CLUSTER_REDSHIFT (test->gz));
  // nc_distance_free (test->dist);
  // ncm_vector_free (test->bin_vector);
  // ncm_vector_free (test->z_vec);
  // ncm_matrix_free (test->data);
  // ncm_rng_free (test->rng);
}
