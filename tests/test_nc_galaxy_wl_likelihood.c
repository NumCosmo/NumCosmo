/***************************************************************************
 *            test_nc_galaxy_wl_likelihood.c
 *
 *  thu December 07 13:54:21 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
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

typedef struct _TestNcGalaxyWLLikelihood
{
  NcGalaxyWLLikelihood *gwl;
} TestNcGalaxyWLLikelihood;


static void test_nc_galaxy_wl_likelihood_new (TestNcGalaxyWLLikelihood *test, gconstpointer pdata);
static void test_nc_galaxy_wl_likelihood_fit (TestNcGalaxyWLLikelihood *test, gconstpointer pdata);
static void test_nc_galaxy_wl_likelihood_free (TestNcGalaxyWLLikelihood *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_wl_likelihood/fit", TestNcGalaxyWLLikelihood, NULL,
              &test_nc_galaxy_wl_likelihood_new,
              &test_nc_galaxy_wl_likelihood_fit,
              &test_nc_galaxy_wl_likelihood_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_wl_likelihood_new (TestNcGalaxyWLLikelihood *test, gconstpointer pdata)
{
  NcGalaxySDPositionFlat *gsdpf = nc_galaxy_sd_position_flat_new ();
  NcGalaxySDZProxyGauss *gsdzpg = nc_galaxy_sd_z_proxy_gauss_new ();
  NcGalaxySDShapeGauss *gss     = nc_galaxy_sd_shape_gauss_new ();
  NcGalaxyWLLikelihood *gwl     = nc_galaxy_wl_likelihood_new (NC_GALAXY_SD_SHAPE (gss),
                                                               NC_GALAXY_SD_Z_PROXY (gsdzpg),
                                                               NC_GALAXY_SD_POSITION (gsdpf));

  {
    NcmVector *z_lim = ncm_vector_new (2);
    NcmVector *r_lim = ncm_vector_new (2);

    ncm_vector_set (z_lim, 0, 0.0);
    ncm_vector_set (z_lim, 1, 2.0);

    ncm_vector_set (r_lim, 0, 0.5);
    ncm_vector_set (r_lim, 1, 1.5);

    nc_galaxy_sd_position_flat_set_z_lim (gsdpf, z_lim);
    nc_galaxy_sd_position_flat_set_r_lim (gsdpf, r_lim);

    ncm_vector_free (z_lim);
    ncm_vector_free (r_lim);
  }

  {
    NcmVector *z_lim = ncm_vector_new (2);

    ncm_vector_set (z_lim, 0, 0.2);
    ncm_vector_set (z_lim, 1, 1.0);

    nc_galaxy_sd_z_proxy_gauss_set_z_lim (gsdzpg, z_lim);
    nc_galaxy_sd_z_proxy_gauss_set_sigma (gsdzpg, 0.03);

    ncm_vector_free (z_lim);
  }

  nc_galaxy_sd_shape_gauss_set_sigma (gss, 0.1);
  nc_galaxy_wl_likelihood_set_cut (gwl, 0.0, 2.0);

  test->gwl = gwl;

  nc_galaxy_sd_position_flat_free (gsdpf);
  nc_galaxy_sd_z_proxy_gauss_free (gsdzpg);
  nc_galaxy_sd_shape_gauss_free (gss);

  g_assert_true (NC_IS_GALAXY_WL_LIKELIHOOD (gwl));
}

static void
test_nc_galaxy_wl_likelihood_free (TestNcGalaxyWLLikelihood *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_wl_likelihood_free, test->gwl);
}

static void
test_nc_galaxy_wl_likelihood_fit (TestNcGalaxyWLLikelihood *test, gconstpointer pdata)
{
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  const guint ngals           = 10000;
  NcmVector *m2lnP_int_gal    = ncm_vector_new (ngals);
  NcmVector *m2lnP_kde_gal    = ncm_vector_new (ngals);
  NcmMatrix *gal_obs          = NULL;
  const gdouble z_cluster     = 0.4;

  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0);

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  nc_galaxy_wl_likelihood_gen_obs (test->gwl, cosmo, dp, smd, z_cluster, ngals, rng);

  {
    NcDataClusterWLL *data_wll = nc_data_cluster_wll_new ();
    NcmObjArray *obs           = ncm_obj_array_new ();

    ncm_obj_array_add (obs, G_OBJECT (test->gwl));
    nc_data_cluster_wll_set_kde (data_wll, FALSE);

    g_object_set (data_wll,
                  "galaxy-array", obs,
                  "z-cluster", z_cluster,
                  NULL);
    ncm_data_set_init (NCM_DATA (data_wll), TRUE);

    {
      NcmDataset *dataset       = ncm_dataset_new_list (data_wll, NULL);
      NcmLikelihood *likelihood = ncm_likelihood_new (dataset);
      NcmMSet *mset             = ncm_mset_new (cosmo, dp, smd, NULL);
      NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", likelihood, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);

      ncm_fit_log_start (fit);
      ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_SIMPLE, 1.0e-3, 0.0, NULL, NULL);
      ncm_fit_log_info (fit);

      ncm_fit_obs_fisher (fit);
      ncm_fit_log_covar (fit);

      ncm_dataset_free (dataset);
      ncm_likelihood_free (likelihood);
      ncm_mset_free (mset);
      ncm_fit_free (fit);
    }

    nc_data_cluster_wll_free (data_wll);
    ncm_obj_array_unref (obs);
  }

  ncm_vector_free (m2lnP_int_gal);
  ncm_vector_free (m2lnP_kde_gal);
  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}

