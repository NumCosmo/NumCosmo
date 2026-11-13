/***************************************************************************
 *            test_nc_data_cluster_wl_factor_shape_interp.c
 *
 *  Thu Nov 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * Tests of nc_data_cluster_wl_factor_write_shape_interp(): a catalog prepared with the
 * interpolants of a NcGalaxyShapeFactorSpectral gives, after a serialization round trip, the
 * -2 ln L of the dataset that prepared it; another shape factor aborts.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define SHAPE_INTERP_NROWS (6)
#define SHAPE_INTERP_Z_CL (0.5)
#define SHAPE_INTERP_SEED (20261112)

typedef struct _TestNcDataClusterWLFactorShapeInterp
{
  NcmMSet *mset;
  NcDataClusterWLFactor *dcwlf;
  NcGalaxyWLObs *obs;
} TestNcDataClusterWLFactorShapeInterp;

void test_nc_data_cluster_wl_factor_shape_interp_new (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_shape_interp_free (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_shape_interp_write (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_shape_interp_shape_factor (void);
void test_nc_data_cluster_wl_factor_shape_interp_shape_factor_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_cluster_wl_factor/shape_interp/write", TestNcDataClusterWLFactorShapeInterp, GINT_TO_POINTER (TRUE),
              &test_nc_data_cluster_wl_factor_shape_interp_new,
              &test_nc_data_cluster_wl_factor_shape_interp_write,
              &test_nc_data_cluster_wl_factor_shape_interp_free);
  g_test_add_func ("/nc/data_cluster_wl_factor/shape_interp/shape_factor", &test_nc_data_cluster_wl_factor_shape_interp_shape_factor);
  g_test_add_func ("/nc/data_cluster_wl_factor/shape_interp/shape_factor/subprocess",
                   &test_nc_data_cluster_wl_factor_shape_interp_shape_factor_subprocess);

  g_test_run ();

  return 0;
}

/*
 * A dataset of SHAPE_INTERP_NROWS generated galaxies around a halo at SHAPE_INTERP_Z_CL,
 * with a NcGalaxyShapeFactorSpectral using its interpolants when pdata is TRUE and a
 * NcGalaxyShapeFactorVarAdd otherwise.
 */
void
test_nc_data_cluster_wl_factor_shape_interp_new (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata)
{
  NcmRNG *rng                             = ncm_rng_seeded_new (NULL, SHAPE_INTERP_SEED);
  NcHICosmo *cosmo                        = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                        = nc_distance_new (100.0);
  NcHaloMassSummary *hms                  = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp                = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                      = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd             = nc_wl_surface_mass_density_new (dist);
  NcGalaxyRedshiftPopLSSTSRD *pop_z       = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z         = nc_galaxy_redshift_obs_gauss_new ();
  NcGalaxyShapePop *pop_shape             = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  const NcGalaxyWLObsEllipConv ec         = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;
  const NcWLEllipticityFrame co           = NC_WL_ELLIPTICITY_FRAME_CELESTIAL;
  const gdouble min_radius                = 0.2;
  const gdouble max_radius                = 4.0;
  NcGalaxyPositionFactor *position_factor = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcGalaxyRedshiftFactor *redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0));
  NcGalaxyShapeFactor *shape_factor;
  GStrvBuilder *builder = g_strv_builder_new ();
  NcGalaxyPositionFactorData *p_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *s_data;
  GList *columns, *l;
  GStrv columns_strv;
  guint i;

  if (GPOINTER_TO_INT (pdata))
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (ec);

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
    shape_factor = NC_GALAXY_SHAPE_FACTOR (gsfs);
  }
  else
  {
    shape_factor = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (ec));
  }

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, 14.5);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_Z, SHAPE_INTERP_Z_CL);
  nc_halo_position_prepare (hp, cosmo);
  ncm_model_param_set (NCM_MODEL (pop_shape), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.25);

  test->mset = ncm_mset_empty_new ();
  ncm_mset_set (test->mset, NCM_MODEL (cosmo), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (dp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (hp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (smd), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (pop_z), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (obs_z), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (pop_shape), NULL);

  p_data = nc_galaxy_position_factor_data_new (position_factor, test->mset);
  z_data = nc_galaxy_redshift_factor_data_new (redshift_factor, test->mset);
  s_data = nc_galaxy_shape_factor_data_new (shape_factor, test->mset, p_data, z_data);

  columns = nc_galaxy_shape_factor_data_required_columns (s_data);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  test->obs    = nc_galaxy_wl_obs_new (ec, co, SHAPE_INTERP_NROWS, columns_strv);

  for (i = 0; i < SHAPE_INTERP_NROWS; i++)
  {
    gdouble radius;

    nc_galaxy_wl_obs_set (test->obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, ncm_rng_uniform_gen (rng, 0.02, 0.06), NULL);
    nc_galaxy_redshift_factor_data_read_row (z_data, test->obs, i);
    nc_galaxy_redshift_factor_gen (redshift_factor, test->mset, z_data, rng);

    do {
      nc_galaxy_position_factor_gen (position_factor, test->mset, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    nc_galaxy_shape_factor_data_set (shape_factor, s_data, 0.0, 0.0, ncm_rng_uniform_gen (rng, 0.02, 0.1), 0.0, 0.0, 0.0, co);
    nc_galaxy_shape_factor_gen (shape_factor, test->mset, s_data, rng);
    nc_galaxy_shape_factor_data_write_row (s_data, test->obs, i);
  }

  test->dcwlf = nc_data_cluster_wl_factor_new (position_factor, redshift_factor, shape_factor);
  nc_data_cluster_wl_factor_set_cut (test->dcwlf, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (test->dcwlf, test->obs);

  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_mass_summary_free (hms);
  nc_halo_density_profile_free (dp);
  nc_halo_position_free (hp);
  nc_wl_surface_mass_density_free (smd);
  nc_galaxy_redshift_pop_lsst_srd_free (pop_z);
  nc_galaxy_redshift_obs_gauss_free (obs_z);
  nc_galaxy_shape_pop_free (pop_shape);
  nc_galaxy_position_factor_free (position_factor);
  nc_galaxy_redshift_factor_free (redshift_factor);
  nc_galaxy_shape_factor_free (shape_factor);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_data_unref (s_data);
  g_strfreev (columns_strv);
  g_strv_builder_unref (builder);
  g_list_free_full (columns, g_free);
  ncm_rng_free (rng);
}

void
test_nc_data_cluster_wl_factor_shape_interp_free (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  nc_galaxy_wl_obs_clear (&test->obs);
  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, test->dcwlf);
}

/*
 * The interpolants written into the catalog of the dataset leave its columns unchanged, and the
 * catalog read back from its serialization gives the same -2 ln L.
 */
void
test_nc_data_cluster_wl_factor_shape_interp_write (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmData *data     = NCM_DATA (test->dcwlf);
  const guint64 key = ncm_catalog_get_key (NCM_CATALOG (test->obs));
  NcGalaxyWLObs *obs_dup;
  GVariant *obs_ser;
  gdouble m2lnL, m2lnL_dup;
  guint i;

  nc_data_cluster_wl_factor_write_shape_interp (test->dcwlf, test->mset, test->obs);
  g_assert_cmpuint (ncm_catalog_get_key (NCM_CATALOG (test->obs)), ==, key);
  g_assert_nonnull (nc_galaxy_wl_obs_peek_shape_interp_source (test->obs));

  for (i = 0; i < SHAPE_INTERP_NROWS; i++)
    g_assert_true (NCM_IS_SPECTRAL_POLAR_PANELS (nc_galaxy_wl_obs_peek_shape_interp (test->obs, i)));

  ncm_data_m2lnL_val (data, test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));

  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->obs));
  ncm_serialize_reset (ser, FALSE);
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));

  nc_data_cluster_wl_factor_set_obs (test->dcwlf, obs_dup);
  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL_dup);
  g_assert_cmpfloat (m2lnL_dup, ==, m2lnL);

  g_variant_unref (obs_ser);
  nc_galaxy_wl_obs_free (obs_dup);
  ncm_serialize_free (ser);
}

/* A shape factor other than NcGalaxyShapeFactorSpectral with use-interp aborts */
void
test_nc_data_cluster_wl_factor_shape_interp_shape_factor (void)
{
  g_test_trap_subprocess ("/nc/data_cluster_wl_factor/shape_interp/shape_factor/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the shape factor NcGalaxyShapeFactorVarAdd is not a NcGalaxyShapeFactorSpectral with use-interp*");
}

void
test_nc_data_cluster_wl_factor_shape_interp_shape_factor_subprocess (void)
{
  TestNcDataClusterWLFactorShapeInterp test;

  test_nc_data_cluster_wl_factor_shape_interp_new (&test, GINT_TO_POINTER (FALSE));
  nc_data_cluster_wl_factor_write_shape_interp (test.dcwlf, test.mset, test.obs);
}

