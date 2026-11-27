/***************************************************************************
 *            test_nc_galaxy_halo_wl.c
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_halo_wl.c
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
 * - /nc/galaxy_halo_wl/radius/moves and /nc/galaxy_halo_wl/shear_rotation, first written in
 *   Python for the shape factor.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_N_GAL (16)
#define TEST_RA (30.0)
#define TEST_DEC (-5.0)
#define TEST_Z_CL (0.3)

typedef struct _TestNcGalaxyHaloWL
{
  NcmMSet *mset;
  NcHICosmo *cosmo;
  NcHaloPosition *hp;
  NcHaloMassSummary *hms;
  NcHaloDensityProfile *dp;
  NcWLSurfaceMassDensity *smd;
  NcGalaxyPositionFactor *posf;
  NcGalaxyRedshiftFactor *zf;
  NcGalaxyShapeFactor *gsf;
  NcGalaxyHaloWL *hwl;
  GPtrArray *pos_data;
  GPtrArray *s_data;
  GPtrArray *h_data;
} TestNcGalaxyHaloWL;

void test_nc_galaxy_halo_wl_new (TestNcGalaxyHaloWL *test, gconstpointer pdata);
void test_nc_galaxy_halo_wl_free (TestNcGalaxyHaloWL *test, gconstpointer pdata);

void test_nc_galaxy_halo_wl_radius (TestNcGalaxyHaloWL *test, gconstpointer pdata);
void test_nc_galaxy_halo_wl_radius_moves (TestNcGalaxyHaloWL *test, gconstpointer pdata);
void test_nc_galaxy_halo_wl_shear_rotation (TestNcGalaxyHaloWL *test, gconstpointer pdata);
void test_nc_galaxy_halo_wl_profile (TestNcGalaxyHaloWL *test, gconstpointer pdata);
void test_nc_galaxy_halo_wl_hash (TestNcGalaxyHaloWL *test, gconstpointer pdata);
void test_nc_galaxy_halo_wl_bound (TestNcGalaxyHaloWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_halo_wl/radius", TestNcGalaxyHaloWL, NULL,
              &test_nc_galaxy_halo_wl_new,
              &test_nc_galaxy_halo_wl_radius,
              &test_nc_galaxy_halo_wl_free);

  g_test_add ("/nc/galaxy_halo_wl/radius/moves", TestNcGalaxyHaloWL, NULL,
              &test_nc_galaxy_halo_wl_new,
              &test_nc_galaxy_halo_wl_radius_moves,
              &test_nc_galaxy_halo_wl_free);

  g_test_add ("/nc/galaxy_halo_wl/shear_rotation", TestNcGalaxyHaloWL, NULL,
              &test_nc_galaxy_halo_wl_new,
              &test_nc_galaxy_halo_wl_shear_rotation,
              &test_nc_galaxy_halo_wl_free);

  g_test_add ("/nc/galaxy_halo_wl/profile", TestNcGalaxyHaloWL, NULL,
              &test_nc_galaxy_halo_wl_new,
              &test_nc_galaxy_halo_wl_profile,
              &test_nc_galaxy_halo_wl_free);

  g_test_add ("/nc/galaxy_halo_wl/hash", TestNcGalaxyHaloWL, NULL,
              &test_nc_galaxy_halo_wl_new,
              &test_nc_galaxy_halo_wl_hash,
              &test_nc_galaxy_halo_wl_free);

  g_test_add ("/nc/galaxy_halo_wl/bound", TestNcGalaxyHaloWL, NULL,
              &test_nc_galaxy_halo_wl_new,
              &test_nc_galaxy_halo_wl_bound,
              &test_nc_galaxy_halo_wl_free);

  g_test_run ();

  return 0;
}

/* Prepares the shape factor and the halo WL object and updates every galaxy */
static void
_test_update (TestNcGalaxyHaloWL *test)
{
  guint i;

  nc_galaxy_shape_factor_prepare (test->gsf, test->mset);
  nc_galaxy_halo_wl_prepare (test->hwl, test->mset);

  for (i = 0; i < TEST_N_GAL; i++)
  {
    NcGalaxyHaloWLData *h_data = g_ptr_array_index (test->h_data, i);

    nc_galaxy_halo_wl_update_data_radius (test->hwl, h_data, g_ptr_array_index (test->pos_data, i));
    nc_galaxy_halo_wl_update_data_profile (test->hwl, test->dp, h_data);
  }
}

void
test_nc_galaxy_halo_wl_new (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  NcDistance *dist = nc_distance_new (3.0);
  guint i;

  test->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  test->hms   = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  test->dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (test->hms));
  test->hp    = nc_halo_position_new (dist);
  test->smd   = nc_wl_surface_mass_density_new (dist);
  test->mset  = ncm_mset_empty_new ();

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 14.7, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", TEST_RA, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", TEST_DEC, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "z", TEST_Z_CL, NULL);

  ncm_mset_set (test->mset, NCM_MODEL (test->cosmo), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (test->dp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (test->hp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (test->smd), NULL);
  {
    NcGalaxyShapePop *pop             = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
    NcGalaxyRedshiftPopLSSTSRD *pop_z = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
    NcGalaxyRedshiftObsGauss *obs_z   = nc_galaxy_redshift_obs_gauss_new ();

    ncm_mset_set (test->mset, NCM_MODEL (pop), NULL);
    ncm_mset_set (test->mset, NCM_MODEL (pop_z), NULL);
    ncm_mset_set (test->mset, NCM_MODEL (obs_z), NULL);

    nc_galaxy_shape_pop_free (pop);
    nc_galaxy_redshift_pop_lsst_srd_free (pop_z);
    nc_galaxy_redshift_obs_gauss_free (obs_z);
  }

  test->posf = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (TEST_RA - 1.0, TEST_RA + 1.0, TEST_DEC - 1.0, TEST_DEC + 1.0));
  test->zf   = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 20.0));
  test->gsf  = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE));
  {
    NcWLClusterLens *lens = nc_wl_cluster_lens_new (dist);

    test->hwl = nc_galaxy_halo_wl_new (lens);
    nc_wl_cluster_lens_free (lens);
  }

  test->pos_data = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_position_factor_data_unref);
  test->s_data   = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_shape_factor_data_unref);
  test->h_data   = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_halo_wl_data_unref);

  for (i = 0; i < TEST_N_GAL; i++)
  {
    const gdouble alpha                = 0.3 + 2.0 * M_PI * i / TEST_N_GAL;
    const gdouble sep                  = 0.02 + 0.03 * i;
    NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->posf, test->mset);

    p_data->ra  = TEST_RA + sep * cos (alpha) / cos (TEST_DEC * M_PI / 180.0);
    p_data->dec = TEST_DEC + sep * sin (alpha);

    g_ptr_array_add (test->pos_data, p_data);
    g_ptr_array_add (test->s_data, nc_galaxy_shape_factor_data_new (test->gsf, test->mset));
    g_ptr_array_add (test->h_data, nc_galaxy_halo_wl_data_new ());
  }

  nc_distance_free (dist);

  _test_update (test);
}

void
test_nc_galaxy_halo_wl_free (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  g_ptr_array_unref (test->h_data);
  g_ptr_array_unref (test->s_data);
  g_ptr_array_unref (test->pos_data);

  NCM_TEST_FREE (nc_galaxy_halo_wl_free, test->hwl);
  nc_galaxy_shape_factor_clear (&test->gsf);
  nc_galaxy_redshift_factor_clear (&test->zf);
  nc_galaxy_position_factor_clear (&test->posf);
  ncm_mset_clear (&test->mset);
  nc_wl_surface_mass_density_clear (&test->smd);
  nc_halo_position_clear (&test->hp);
  nc_halo_density_profile_clear (&test->dp);
  nc_halo_mass_summary_clear (&test->hms);
  nc_hicosmo_clear (&test->cosmo);
}

/*
 * The projected radius and the rotation factor (cos 2 phi, sin 2 phi) agree to rounding with
 * those of the polar angles of nc_halo_position_polar_angles().
 */
void
test_nc_galaxy_halo_wl_radius (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  guint i;

  for (i = 0; i < TEST_N_GAL; i++)
  {
    NcGalaxyPositionFactorData *p_data = g_ptr_array_index (test->pos_data, i);
    NcGalaxyHaloWLData *h_data         = g_ptr_array_index (test->h_data, i);
    gdouble theta, phi, rot_re, rot_im;

    nc_halo_position_polar_angles (test->hp, p_data->ra, p_data->dec, &theta, &phi);
    ncm_assert_cmpdouble_e (nc_galaxy_halo_wl_data_get_radius (h_data), ==, nc_halo_position_projected_radius (test->hp, test->cosmo, theta), 1.0e-13, 0.0);
    nc_galaxy_halo_wl_data_get_rotation (h_data, &rot_re, &rot_im);

    ncm_assert_cmpdouble_e (rot_re, ==, cos (2.0 * phi), 0.0, 1.0e-12);
    ncm_assert_cmpdouble_e (rot_im, ==, sin (2.0 * phi), 0.0, 1.0e-12);
  }
}

/*
 * Offsets in degrees of the halo center and of the galaxy from (TEST_RA, TEST_DEC): the
 * halo moves with the galaxy fixed, then the galaxy moves, then the halo again.
 */
static const gdouble test_moves[4][4] = {
  {0.0, 0.0, 0.03, 0.02},
  {0.01, -0.02, 0.03, 0.02},
  {0.01, -0.02, -0.1, 0.15},
  {-0.05, 0.05, -0.1, 0.15},
};

/*
 * The projected radius of a galaxy is that of the polar angle of
 * nc_halo_position_polar_angles() after the halo center moves (the unit vector of the
 * galaxy is kept and rotated) and after the galaxy moves (the vector is computed again).
 */
void
test_nc_galaxy_halo_wl_radius_moves (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->posf, test->mset);
  NcGalaxyHaloWLData *h_data         = nc_galaxy_halo_wl_data_new ();
  guint k;

  for (k = 0; k < G_N_ELEMENTS (test_moves); k++)
  {
    gdouble theta, phi;

    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", TEST_RA + test_moves[k][0], NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", TEST_DEC + test_moves[k][1], NULL);
    p_data->ra  = TEST_RA + test_moves[k][2];
    p_data->dec = TEST_DEC + test_moves[k][3];

    nc_galaxy_halo_wl_prepare (test->hwl, test->mset);
    nc_galaxy_halo_wl_update_data_radius (test->hwl, h_data, p_data);

    nc_halo_position_polar_angles (test->hp, p_data->ra, p_data->dec, &theta, &phi);
    ncm_assert_cmpdouble_e (nc_galaxy_halo_wl_data_get_radius (h_data), ==, nc_halo_position_projected_radius (test->hp, test->cosmo, theta), 1.0e-13, 0.0);
  }

  nc_galaxy_halo_wl_data_unref (h_data);
  nc_galaxy_position_factor_data_unref (p_data);
}

/*
 * The marginal of a galaxy at a source redshift, evaluated at the reduced shear of
 * #NcWLClusterLens times the rotation factor of #NcGalaxyHaloWL (conjugated for CARTESIAN),
 * agrees with the marginal at the shear from the azimuth of nc_halo_position_polar_angles()
 * taken to the frame of the galaxy by nc_wl_ellipticity_celestial_to_frame_angle(), and the
 * reduced shear of the lens agrees with that of #NcWLSurfaceMassDensity. The other frame
 * gives a different value, so the frame conversion is seen, in both frames and after the
 * halo center and the galaxy move.
 */
void
test_nc_galaxy_halo_wl_shear_rotation (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  const NcWLEllipticityFrame frames[2] = {NC_WL_ELLIPTICITY_FRAME_CELESTIAL, NC_WL_ELLIPTICITY_FRAME_CARTESIAN};
  const gdouble e_1                    = 0.05;
  const gdouble e_2                    = -0.02;
  const gdouble z_s                    = 0.9;
  NcGalaxyShapePop *pop                = NC_GALAXY_SHAPE_POP (ncm_mset_peek (test->mset, nc_galaxy_shape_pop_id ()));
  NcWLClusterLens *lens                = nc_galaxy_halo_wl_peek_lens (test->hwl);
  NcGalaxyPositionFactorData *p_data   = nc_galaxy_position_factor_data_new (test->posf, test->mset);
  NcGalaxyHaloWLData *h_data           = nc_galaxy_halo_wl_data_new ();
  NcGalaxyShapeFactorData *s_data      = nc_galaxy_shape_factor_data_new (test->gsf, test->mset);
  GPtrArray *s_array                   = g_ptr_array_new ();
  guint f, k;

  g_ptr_array_add (s_array, s_data);

  for (f = 0; f < 2; f++)
  {
    nc_galaxy_shape_factor_data_set (test->gsf, s_data, e_1, e_2, 0.03, 0.0, 0.0, 0.0, frames[f]);

    for (k = 0; k < G_N_ELEMENTS (test_moves); k++)
    {
      gdouble sigma, mean_sigma, rot_re, rot_im, gt_lens, gt, theta, phi, ln_p[2];
      guint q;

      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", TEST_RA + test_moves[k][0], NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", TEST_DEC + test_moves[k][1], NULL);
      p_data->ra  = TEST_RA + test_moves[k][2];
      p_data->dec = TEST_DEC + test_moves[k][3];

      g_assert_true (nc_galaxy_shape_factor_prepare_data_array (test->gsf, test->mset, s_array));
      nc_galaxy_halo_wl_prepare (test->hwl, test->mset);
      nc_galaxy_halo_wl_update_data_radius (test->hwl, h_data, p_data);
      nc_galaxy_halo_wl_update_data_profile (test->hwl, test->dp, h_data);
      nc_galaxy_halo_wl_data_get_sigma (h_data, &sigma, &mean_sigma);
      nc_galaxy_halo_wl_data_get_rotation (h_data, &rot_re, &rot_im);

      if (frames[f] == NC_WL_ELLIPTICITY_FRAME_CARTESIAN)
        rot_im = -rot_im;

      gt_lens = nc_wl_cluster_lens_reduced_shear_at (lens, test->cosmo, sigma, mean_sigma, z_s);
      gt      = nc_wl_surface_mass_density_reduced_shear (test->smd, test->dp, test->cosmo, nc_galaxy_halo_wl_data_get_radius (h_data), z_s, TEST_Z_CL, TEST_Z_CL);
      ncm_assert_cmpdouble_e (gt_lens, ==, gt, 1.0e-12, 0.0);

      nc_halo_position_polar_angles (test->hp, p_data->ra, p_data->dec, &theta, &phi);

      for (q = 0; q < 2; q++)
      {
        const gdouble phi_q = nc_wl_ellipticity_celestial_to_frame_angle (frames[q], phi);

        ln_p[q] = nc_galaxy_shape_factor_eval_ln_marginal (test->gsf, pop, s_data, gt * cos (2.0 * phi_q), gt * sin (2.0 * phi_q), e_1, e_2);
      }

      ncm_assert_cmpdouble_e (nc_galaxy_shape_factor_eval_shear (test->gsf, pop, s_data, gt_lens * rot_re, gt_lens * rot_im, TRUE), ==, ln_p[f], 1.0e-9, 0.0);
      g_assert_cmpfloat (fabs (ln_p[1 - f] - ln_p[f]), >, 1.0e-3);
    }
  }

  g_ptr_array_unref (s_array);
  nc_galaxy_shape_factor_data_unref (s_data);
  nc_galaxy_halo_wl_data_unref (h_data);
  nc_galaxy_position_factor_data_unref (p_data);
}

/* The surface densities equal, bit for bit, those of the cached NcWLSurfaceMassDensity path */
void
test_nc_galaxy_halo_wl_profile (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  gdouble r_s, rho_s;
  guint i;

  nc_halo_density_profile_r_s_rho_s (test->dp, test->cosmo, TEST_Z_CL, &r_s, &rho_s);

  for (i = 0; i < TEST_N_GAL; i++)
  {
    NcGalaxyHaloWLData *h_data = g_ptr_array_index (test->h_data, i);
    NcWLSurfaceMassDensitySigmaCache sigma_cache;
    gdouble sigma, mean_sigma;

    nc_wl_surface_mass_density_reduced_shear_sigma_cache_prep_with_rs (test->dp, nc_galaxy_halo_wl_data_get_radius (h_data), r_s, rho_s, &sigma_cache);
    nc_galaxy_halo_wl_data_get_sigma (h_data, &sigma, &mean_sigma);

    g_assert_cmpfloat (sigma, ==, sigma_cache.sigma);
    g_assert_cmpfloat (mean_sigma, ==, sigma_cache.mean_sigma);
  }
}

/*
 * The radius hash changes with the halo center, not with the mass; the profile hash
 * changes with the mass. A new position recomputes the radius of that galaxy.
 */
void
test_nc_galaxy_halo_wl_hash (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  const guint64 radius_hash0  = nc_galaxy_halo_wl_get_radius_hash (test->hwl);
  const guint64 profile_hash0 = nc_galaxy_halo_wl_get_profile_hash (test->hwl);
  NcGalaxyHaloWLData *h_data  = g_ptr_array_index (test->h_data, 0);
  const gdouble R0            = nc_galaxy_halo_wl_data_get_radius (h_data);

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 14.9, NULL);
  nc_galaxy_halo_wl_prepare (test->hwl, test->mset);
  g_assert_cmpuint (nc_galaxy_halo_wl_get_radius_hash (test->hwl), ==, radius_hash0);
  g_assert_cmpuint (nc_galaxy_halo_wl_get_profile_hash (test->hwl), !=, profile_hash0);

  ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", TEST_RA + 0.01, NULL);
  _test_update (test);
  g_assert_cmpuint (nc_galaxy_halo_wl_get_radius_hash (test->hwl), !=, radius_hash0);
  g_assert_cmpfloat (nc_galaxy_halo_wl_data_get_radius (h_data), !=, R0);
  test_nc_galaxy_halo_wl_radius (test, pdata);

  {
    NcGalaxyPositionFactorData *p_data = g_ptr_array_index (test->pos_data, 0);

    p_data->dec += 0.05;
    _test_update (test);
    test_nc_galaxy_halo_wl_radius (test, pdata);
  }
}

/*
 * The bound holds when only the mass and the center are free; its corner evaluation
 * equals the convergence and the reduced shear of a profile with the mass at its upper
 * bound, and the smallest radius is R minus the largest displacement of the center. With a
 * profile of another class put in the mset, the corner is a copy of that profile.
 */
void
test_nc_galaxy_halo_wl_bound (TestNcGalaxyHaloWL *test, gconstpointer pdata)
{
  NcGalaxyHaloWLData *h_data = g_ptr_array_index (test->h_data, 3);
  const gdouble R            = nc_galaxy_halo_wl_data_get_radius (h_data);
  NcHaloMassSummary *hms_box = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  NcHaloDensityProfile *dp_box;
  gdouble kappa, g_t;

  ncm_model_param_set_ftype (NCM_MODEL (test->hms), NC_HALO_CM_PARAM_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  nc_galaxy_halo_wl_prepare_bound (test->hwl, test->mset);
  g_assert_true (nc_galaxy_halo_wl_bound_ready (test->hwl));
  g_assert_true (nc_galaxy_halo_wl_bound_holds (test->hwl));
  g_assert_cmpfloat (nc_galaxy_halo_wl_get_radius_min (test->hwl, h_data), ==, R);

  ncm_model_param_set_by_name (NCM_MODEL (hms_box), "log10MDelta", ncm_model_param_get_upper_bound (NCM_MODEL (test->hms), NC_HALO_CM_PARAM_LOG10M_DELTA), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (hms_box), "cDelta", ncm_model_param_get_by_name (NCM_MODEL (test->hms), "cDelta", NULL), NULL);
  dp_box = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms_box));

  nc_galaxy_halo_wl_eval_bound (test->hwl, test->smd, test->cosmo, R, 1.5, &kappa, &g_t);
  g_assert_cmpfloat (kappa, ==, nc_wl_surface_mass_density_convergence (test->smd, dp_box, test->cosmo, R, 1.5, TEST_Z_CL, TEST_Z_CL));
  g_assert_cmpfloat (g_t, ==, fabs (nc_wl_surface_mass_density_reduced_shear (test->smd, dp_box, test->cosmo, R, 1.5, TEST_Z_CL, TEST_Z_CL)));

  ncm_model_param_set_ftype (NCM_MODEL (test->hp), NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  nc_galaxy_halo_wl_prepare_bound (test->hwl, test->mset);
  g_assert_true (nc_galaxy_halo_wl_bound_holds (test->hwl));
  g_assert_cmpfloat (nc_galaxy_halo_wl_get_radius_min (test->hwl, h_data), <, R);

  ncm_model_param_set_ftype (NCM_MODEL (test->hms), NC_HALO_CM_PARAM_C_DELTA, NCM_PARAM_TYPE_FREE);
  nc_galaxy_halo_wl_prepare_bound (test->hwl, test->mset);
  g_assert_false (nc_galaxy_halo_wl_bound_holds (test->hwl));

  /* A profile of another class in the mset: the corner is a copy of it */
  {
    NcHaloMassSummary *hms_h     = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
    NcHaloMassSummary *hms_h_box = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
    NcHaloDensityProfile *dp_h   = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_hernquist_new (hms_h));
    NcHaloDensityProfile *dp_h_box;

    ncm_model_param_set_by_name (NCM_MODEL (hms_h), "log10MDelta", 14.2, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (hms_h), "cDelta", 6.0, NULL);
    ncm_mset_set (test->mset, NCM_MODEL (dp_h), NULL);
    nc_galaxy_halo_wl_prepare_bound (test->hwl, test->mset);
    g_assert_true (nc_galaxy_halo_wl_bound_ready (test->hwl));

    ncm_model_param_set_by_name (NCM_MODEL (hms_h_box), "log10MDelta", ncm_model_param_get_upper_bound (NCM_MODEL (hms_h), NC_HALO_CM_PARAM_LOG10M_DELTA), NULL);
    ncm_model_param_set_by_name (NCM_MODEL (hms_h_box), "cDelta", 6.0, NULL);
    dp_h_box = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_hernquist_new (hms_h_box));

    nc_galaxy_halo_wl_eval_bound (test->hwl, test->smd, test->cosmo, R, 1.5, &kappa, &g_t);
    g_assert_cmpfloat (kappa, ==, nc_wl_surface_mass_density_convergence (test->smd, dp_h_box, test->cosmo, R, 1.5, TEST_Z_CL, TEST_Z_CL));
    g_assert_cmpfloat (g_t, ==, fabs (nc_wl_surface_mass_density_reduced_shear (test->smd, dp_h_box, test->cosmo, R, 1.5, TEST_Z_CL, TEST_Z_CL)));

    nc_halo_density_profile_free (dp_h_box);
    nc_halo_density_profile_free (dp_h);
    nc_halo_mass_summary_free (hms_h_box);
    nc_halo_mass_summary_free (hms_h);
  }

  nc_halo_density_profile_free (dp_box);
  nc_halo_mass_summary_free (hms_box);
}

