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


typedef struct _TestNcGalaxyWLEllipticityKDE
{
  NcGalaxyWLEllipticityKDE *gekde;
  NcHICosmo *cosmo;
  NcHaloDensityProfile *dp;
  NcWLSurfaceMassDensity *smd;
  NcGalaxyRedshift *gz;
  gdouble z_cluster;
} TestNcGalaxyWLEllipticityKDE;

static void test_nc_galaxy_wl_ellipticity_kde_new (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_kde_reset (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_kde_e_vec (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_kde_m2lnP (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata);
static void test_nc_galaxy_wl_ellipticity_kde_free (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/galaxy_wl_dist_ellipticity_kde/reset", TestNcGalaxyWLEllipticityKDE, NULL,
              &test_nc_galaxy_wl_ellipticity_kde_new,
              &test_nc_galaxy_wl_ellipticity_kde_reset,
              &test_nc_galaxy_wl_ellipticity_kde_free);

  g_test_add ("/nc/galaxy_wl_dist_ellipticity_kde/e_vec", TestNcGalaxyWLEllipticityKDE, NULL,
              &test_nc_galaxy_wl_ellipticity_kde_new,
              &test_nc_galaxy_wl_ellipticity_kde_e_vec,
              &test_nc_galaxy_wl_ellipticity_kde_free);

  g_test_add ("/nc/galaxy_wl_dist_ellipticity_kde/m2lnp", TestNcGalaxyWLEllipticityKDE, NULL,
              &test_nc_galaxy_wl_ellipticity_kde_new,
              &test_nc_galaxy_wl_ellipticity_kde_m2lnP,
              &test_nc_galaxy_wl_ellipticity_kde_free);

  g_test_run ();
}

static void 
test_nc_galaxy_wl_ellipticity_kde_new (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata)
{
  NcGalaxyWLEllipticityKDE *gekde = nc_galaxy_wl_ellipticity_kde_new ();
  NcHICosmo *cosmo                = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcHaloDensityProfile *dp        = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL, 200.0));
  NcDistance *dist                = nc_distance_new (3.0);
  NcWLSurfaceMassDensity *smd     = nc_wl_surface_mass_density_new (dist);
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ndata               = 10000;
  NcmMatrix *data                 = ncm_matrix_new (ndata, 3);
  NcmVector *z_vec                = ncm_vector_new (ndata);
  NcGalaxyRedshiftSpec *gzs       = nc_galaxy_redshift_spec_new ();
  const gdouble rl                = g_test_rand_double_range (0.0002, 0.1300);
  const gdouble ru                = g_test_rand_double_range (5.55, 5.65);
  const gdouble mu_e              = g_test_rand_double_range (0.015, 0.020);
  const gdouble sigma_e           = g_test_rand_double_range (0.050, 0.055);
  const gdouble e_noise           = g_test_rand_double_range (0.02, 0.08);
  const gdouble mu_z              = g_test_rand_double_range (1.30, 1.35);
  const gdouble sigma_z           = g_test_rand_double_range (0.65, 0.75);
  const gdouble z_cluster         = g_test_rand_double_range (0.2, 0.4);
  guint i;

  for (i = 0; i < ndata; i++)
  {
    const gdouble r = ncm_rng_uniform_gen (rng, rl, ru);
    const gdouble e = ncm_rng_gaussian_gen (rng, mu_e, sigma_e);
    const gdouble z = ncm_rng_gaussian_gen (rng, mu_z, sigma_z);

    ncm_matrix_set (data, i, 0, r);
    ncm_matrix_set (data, i, 1, e);
    ncm_matrix_set (data, i, 2, e_noise);

    ncm_vector_set (z_vec, i, z);
  }

  nc_galaxy_wl_ellipticity_kde_set_obs (gekde, data);
  nc_galaxy_redshift_spec_set_z (gzs, z_vec);

  test->gekde     = nc_galaxy_wl_ellipticity_kde_new ();
  test->cosmo     = cosmo;
  test->dp        = dp;
  test->smd       = smd;
  test->gz        = NC_GALAXY_REDSHIFT (gzs);
  test->z_cluster = z_cluster;

  g_assert_true (NC_IS_GALAXY_WL_ELLIPTICITY_KDE (gekde));
}

static void
test_nc_galaxy_wl_ellipticity_kde_free (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_wl_ellipticity_kde_free, test->gekde);
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
  NCM_TEST_FREE (nc_halo_density_profile_free, test->dp);
  NCM_TEST_FREE (nc_wl_surface_mass_density_free, test->smd);
  NCM_TEST_FREE (nc_galaxy_redshift_free, test->gz);
}

static void
test_nc_galaxy_wl_ellipticity_kde_reset (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata)
{
  NcGalaxyWLEllipticityKDE *gekde = test->gekde;
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;
  const guint nruns = 100000;
  gint i;

  for (i = 0; i < nruns; i++)
  {
    const gdouble h;

    nc_galaxy_wl_dist_m2lnP_initial_prep (NC_GALAXY_WL_DIST (test->gekde), test->gz, test->cosmo, test->dp, test->smd, test->z_cluster);
    h = ncm_stats_dist1d_get_current_h (NCM_STATS_DIST1D (self->kde));

    g_assert_cmpfloat (h, <, 1);
  }
}

static void 
test_nc_galaxy_wl_ellipticity_kde_e_vec (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata)
{
  NcGalaxyWLEllipticityKDE *gekde = test->gekde;
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;
  const guint nruns = 100000;
  const guint ndata = 10000;
  gint i;

  for (i = 0; i < nruns; i++)
  {
    nc_galaxy_wl_dist_m2lnP_initial_prep (NC_GALAXY_WL_DIST (test->gekde), test->gz, test->cosmo, test->dp, test->smd, test->z_cluster);

    g_assert_cmpfloat (ncm_vector_len (self->e_vec), <=, ndata);
  }
}

static void 
test_nc_galaxy_wl_ellipticity_kde_m2lnP (TestNcGalaxyWLEllipticityKDE *test, gconstpointer pdata)
{
  NcGalaxyWLEllipticityKDE *gekde = test->gekde;
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;
  const guint nruns = 100000;
  gint i;

  for (i = 0; i < nruns; i++)
  {
    const gdouble p;
    const gdouble e_i;
    gint gal_i;

    nc_galaxy_wl_dist_m2lnP_initial_prep (NC_GALAXY_WL_DIST (test->gekde), test->gz, test->cosmo, test->dp, test->smd, test->z_cluster);
    gal_i = g_test_rand_int_range (0, 10000);
    e_i = ncm_vector_get (self->e_vec, gal_i);

    p = ncm_stats_dist1d_eval_p (NCM_STATS_DIST1D (self->kde), e_i);

    g_assert_cmpfloat (p, >, 0);
  }
}