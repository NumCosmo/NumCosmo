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
 * -2 ln L of the dataset that prepared it, also with a free population parameter; another
 * shape factor aborts.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <glib/gstdio.h>

/* Every column of these catalogs is DOUBLE (NCM_CATALOG_COL_TYPE_DOUBLE is zero) */
static const NcmCatalogColType test_double_cols[64] = {NCM_CATALOG_COL_TYPE_DOUBLE};

#define SHAPE_INTERP_NROWS (6)
#define SHAPE_INTERP_Z_CL (0.5)
#define SHAPE_INTERP_SEED (20261112)
#define SHAPE_INTERP_SIGMA_LB (0.24)
#define SHAPE_INTERP_SIGMA_UB (0.30)
#define SHAPE_INTERP_POP_ABSTOL (1.0e-3)

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
void test_nc_data_cluster_wl_factor_shape_interp_free_pop (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_shape_interp_file (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_shape_interp_file_traps (void);
static void _test_nc_data_cluster_wl_factor_shape_interp_free_sigma (TestNcDataClusterWLFactorShapeInterp *test);
void test_nc_data_cluster_wl_factor_shape_interp_file_shape_factor_subprocess (void);
void test_nc_data_cluster_wl_factor_shape_interp_file_rows_subprocess (void);
void test_nc_data_cluster_wl_factor_shape_interp_file_create_subprocess (void);

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
  g_test_add ("/nc/data_cluster_wl_factor/shape_interp/free_pop", TestNcDataClusterWLFactorShapeInterp, GINT_TO_POINTER (TRUE),
              &test_nc_data_cluster_wl_factor_shape_interp_new,
              &test_nc_data_cluster_wl_factor_shape_interp_free_pop,
              &test_nc_data_cluster_wl_factor_shape_interp_free);
  g_test_add ("/nc/data_cluster_wl_factor/shape_interp/file", TestNcDataClusterWLFactorShapeInterp, GINT_TO_POINTER (TRUE),
              &test_nc_data_cluster_wl_factor_shape_interp_new,
              &test_nc_data_cluster_wl_factor_shape_interp_file,
              &test_nc_data_cluster_wl_factor_shape_interp_free);
  g_test_add_func ("/nc/data_cluster_wl_factor/shape_interp/file/traps", &test_nc_data_cluster_wl_factor_shape_interp_file_traps);
  g_test_add_func ("/nc/data_cluster_wl_factor/shape_interp/file/shape_factor/subprocess",
                   &test_nc_data_cluster_wl_factor_shape_interp_file_shape_factor_subprocess);
  g_test_add_func ("/nc/data_cluster_wl_factor/shape_interp/file/rows/subprocess",
                   &test_nc_data_cluster_wl_factor_shape_interp_file_rows_subprocess);
  g_test_add_func ("/nc/data_cluster_wl_factor/shape_interp/file/create/subprocess",
                   &test_nc_data_cluster_wl_factor_shape_interp_file_create_subprocess);

  g_test_run ();

  return 0;
}

/*
 * Makes sigma a variable of the interpolants on [SHAPE_INTERP_SIGMA_LB, SHAPE_INTERP_SIGMA_UB],
 * at interp-abstol 1e-3 and level cap 3. A table over this range builds in about a seventh of
 * the time of one over the default bounds; the accuracy over wide ranges is tested in
 * test_nc_galaxy_shape_factor_spectral.c.
 */
static void
_test_nc_data_cluster_wl_factor_shape_interp_free_sigma (TestNcDataClusterWLFactorShapeInterp *test)
{
  NcmModel *pop = ncm_mset_peek (test->mset, nc_galaxy_shape_pop_id ());
  NcGalaxyShapeFactorSpectral *gsfs;

  g_object_get (test->dcwlf, "shape-factor", &gsfs, NULL);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, SHAPE_INTERP_POP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs, 3);
  nc_galaxy_shape_factor_spectral_free (gsfs);

  ncm_model_param_set_lower_bound (pop, NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, SHAPE_INTERP_SIGMA_LB);
  ncm_model_param_set_upper_bound (pop, NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, SHAPE_INTERP_SIGMA_UB);
  ncm_mset_param_set_ftype (test->mset, nc_galaxy_shape_pop_id (), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
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
    shape_factor = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_spectral_new (ec));
  else
    shape_factor = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (ec));

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
  s_data = nc_galaxy_shape_factor_data_new (shape_factor, test->mset);

  test->dcwlf = nc_data_cluster_wl_factor_new (position_factor, redshift_factor, shape_factor);
  columns     = nc_data_cluster_wl_factor_required_columns (test->dcwlf, test->mset);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  test->obs    = nc_galaxy_wl_obs_new (ec, co, SHAPE_INTERP_NROWS, columns_strv, test_double_cols, g_strv_length (columns_strv));

  for (i = 0; i < SHAPE_INTERP_NROWS; i++)
  {
    gdouble radius;

    ncm_catalog_set (NCM_CATALOG (test->obs), NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, ncm_rng_uniform_gen (rng, 0.02, 0.06), NULL);
    nc_galaxy_redshift_factor_data_read_row (z_data, test->obs, i);
    nc_galaxy_redshift_factor_gen (redshift_factor, test->mset, z_data, rng);

    do {
      nc_galaxy_position_factor_gen (position_factor, test->mset, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    nc_galaxy_shape_factor_data_set (shape_factor, s_data, 0.0, 0.0, ncm_rng_uniform_gen (rng, 0.02, 0.1), 0.0, 0.0, 0.0, co);
    nc_galaxy_position_factor_data_write_row (p_data, test->obs, i);
    nc_galaxy_redshift_factor_data_write_row (z_data, test->obs, i);
    nc_galaxy_shape_factor_data_write_row (s_data, test->obs, i);
  }

  nc_data_cluster_wl_factor_set_cut (test->dcwlf, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (test->dcwlf, test->obs);
  nc_data_cluster_wl_factor_set_resample_flag (test->dcwlf, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_SHAPE);
  ncm_data_resample (NCM_DATA (test->dcwlf), test->mset, rng);

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
    g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (test->obs, i)));

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

/* A shape factor other than NcGalaxyShapeFactorSpectral aborts */
void
test_nc_data_cluster_wl_factor_shape_interp_shape_factor (void)
{
  g_test_trap_subprocess ("/nc/data_cluster_wl_factor/shape_interp/shape_factor/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the shape factor NcGalaxyShapeFactorVarAdd is not a NcGalaxyShapeFactorSpectral.*");
}

void
test_nc_data_cluster_wl_factor_shape_interp_shape_factor_subprocess (void)
{
  TestNcDataClusterWLFactorShapeInterp test;

  test_nc_data_cluster_wl_factor_shape_interp_new (&test, GINT_TO_POINTER (FALSE));
  nc_data_cluster_wl_factor_write_shape_interp (test.dcwlf, test.mset, test.obs);
}

/*
 * With sigma a variable (_test_nc_data_cluster_wl_factor_shape_interp_free_sigma()), the
 * interpolants written into the catalog have sigma as a third variable, and -2 ln L of the
 * dataset at sigma 0.25 and 0.29 agrees with the one of the same catalog without interpolants
 * (the integrals) within 2 times the number of galaxies times interp-abstol, the bound of the
 * change of each galaxy's ln P.
 */
void
test_nc_data_cluster_wl_factor_shape_interp_free_pop (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata)
{
  const gdouble sigmas[] = {0.25, 0.29};
  NcmSerialize *ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmData *data          = NCM_DATA (test->dcwlf);
  NcGalaxyWLObs *obs_plain;
  GVariant *obs_ser;
  guint i, k;

  _test_nc_data_cluster_wl_factor_shape_interp_free_sigma (test);

  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->obs));
  ncm_serialize_reset (ser, FALSE);
  obs_plain = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));

  nc_data_cluster_wl_factor_write_shape_interp (test->dcwlf, test->mset, test->obs);

  for (i = 0; i < SHAPE_INTERP_NROWS; i++)
    g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_peek_shape_interp (test->obs, i)), ==, 3);

  for (k = 0; k < G_N_ELEMENTS (sigmas); k++)
  {
    gdouble m2lnL, m2lnL_plain;

    ncm_mset_param_set (test->mset, nc_galaxy_shape_pop_id (), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, sigmas[k]);

    nc_data_cluster_wl_factor_set_obs (test->dcwlf, test->obs);
    ncm_data_m2lnL_val (data, test->mset, &m2lnL);
    nc_data_cluster_wl_factor_set_obs (test->dcwlf, obs_plain);
    ncm_data_m2lnL_val (data, test->mset, &m2lnL_plain);

    g_assert_true (gsl_finite (m2lnL));
    g_test_message ("sigma %g: -2lnL %.15g, without interpolants %.15g, difference %.3e", sigmas[k], m2lnL, m2lnL_plain, m2lnL - m2lnL_plain);
    g_assert_cmpfloat (fabs (m2lnL - m2lnL_plain), <=, 2.0 * SHAPE_INTERP_NROWS * SHAPE_INTERP_POP_ABSTOL);
  }

  g_variant_unref (obs_ser);
  nc_galaxy_wl_obs_free (obs_plain);
  ncm_serialize_free (ser);
}

/*
 * nc_data_cluster_wl_factor_write_shape_interp_file() against
 * nc_data_cluster_wl_factor_write_shape_interp(), with sigma a variable of the interpolants
 * (_test_nc_data_cluster_wl_factor_shape_interp_free_sigma()): the temporary records file is
 * removed and a records file left by an earlier call is overwritten; the catalog given to the
 * call refers to the interpolants file, with the record; the catalog written next to that
 * file and read back refers to it too and gives, row by row, interpolants whose tree equals
 * the one built in memory, and the -2 ln L of the catalog prepared in memory.
 */
void
test_nc_data_cluster_wl_factor_shape_interp_file (TestNcDataClusterWLFactorShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmData *data     = NCM_DATA (test->dcwlf);
  GError *error     = NULL;
  gchar *dir        = g_dir_make_tmp ("nc_shape_interp_XXXXXX", &error);
  gchar *tables, *catalog, *rec_name;
  NcGalaxyWLObs *obs_mem, *obs_file;
  NcGalaxyWLObsShapeInterp *si;
  gdouble m2lnL, m2lnL_file;
  guint i;

  g_assert_no_error (error);
  tables   = g_build_filename (dir, "catalog.shape_interp.gvar", NULL);
  catalog  = g_build_filename (dir, "catalog.gvar", NULL);
  rec_name = g_strconcat (tables, ".records", NULL);

  _test_nc_data_cluster_wl_factor_shape_interp_free_sigma (test);
  obs_mem = NC_GALAXY_WL_OBS (ncm_serialize_dup_obj (ser, G_OBJECT (test->obs)));
  ncm_serialize_reset (ser, TRUE);

  /* A file left by an interrupted call is not read */
  g_assert_true (g_file_set_contents (rec_name, "leftover", -1, &error));
  g_assert_no_error (error);

  nc_data_cluster_wl_factor_write_shape_interp_file (test->dcwlf, test->mset, test->obs, tables);
  g_assert_false (g_file_test (rec_name, G_FILE_TEST_EXISTS));
  g_assert_nonnull (nc_galaxy_wl_obs_peek_shape_interp_source (test->obs));
  g_object_get (test->obs, "shape-interp", &si, NULL);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, tables);
  nc_galaxy_wl_obs_shape_interp_free (si);

  ncm_serialize_to_binfile (ser, G_OBJECT (test->obs), catalog);
  ncm_serialize_reset (ser, TRUE);
  obs_file = NC_GALAXY_WL_OBS (ncm_serialize_from_binfile (ser, catalog));
  ncm_serialize_reset (ser, TRUE);
  g_assert_cmpuint (ncm_catalog_len (NCM_CATALOG (obs_file)), ==, SHAPE_INTERP_NROWS);
  g_assert_nonnull (nc_galaxy_wl_obs_peek_shape_interp_source (obs_file));
  g_object_get (obs_file, "shape-interp", &si, NULL);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, tables);
  nc_galaxy_wl_obs_shape_interp_free (si);

  nc_data_cluster_wl_factor_write_shape_interp (test->dcwlf, test->mset, obs_mem);

  for (i = 0; i < SHAPE_INTERP_NROWS; i++)
  {
    NcmSpectralTensor *st_file = nc_galaxy_wl_obs_peek_shape_interp (obs_file, i);
    NcmSpectralTensor *st_mem  = nc_galaxy_wl_obs_peek_shape_interp (obs_mem, i);
    GVariant *tree_file, *tree_mem;

    g_assert_true (NCM_IS_SPECTRAL_TENSOR (st_file));
    g_assert_true (NCM_IS_SPECTRAL_TENSOR (st_mem));
    g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (st_file), ==, 3);

    g_object_get (st_file, "tree", &tree_file, NULL);
    g_object_get (st_mem, "tree", &tree_mem, NULL);
    g_assert_true (g_variant_equal (tree_file, tree_mem));
    g_variant_unref (tree_file);
    g_variant_unref (tree_mem);
  }

  nc_data_cluster_wl_factor_set_obs (test->dcwlf, obs_mem);
  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));

  nc_data_cluster_wl_factor_set_obs (test->dcwlf, obs_file);
  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);
  ncm_data_m2lnL_val (data, test->mset, &m2lnL_file);
  g_assert_cmpfloat (m2lnL_file, ==, m2lnL);

  g_assert_cmpint (g_remove (tables), ==, 0);
  g_assert_cmpint (g_remove (catalog), ==, 0);
  g_assert_cmpint (g_rmdir (dir), ==, 0);

  nc_galaxy_wl_obs_free (obs_file);
  nc_galaxy_wl_obs_free (obs_mem);
  ncm_serialize_free (ser);
  g_free (rec_name);
  g_free (catalog);
  g_free (tables);
  g_free (dir);
}

/*
 * nc_data_cluster_wl_factor_write_shape_interp_file() aborts, before it builds any
 * interpolant, when the shape factor is not a NcGalaxyShapeFactorSpectral, when the catalog
 * has another number of rows than the dataset, and when the records file cannot be created.
 */
void
test_nc_data_cluster_wl_factor_shape_interp_file_traps (void)
{
  g_test_trap_subprocess ("/nc/data_cluster_wl_factor/shape_interp/file/shape_factor/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the shape factor NcGalaxyShapeFactorVarAdd is not a NcGalaxyShapeFactorSpectral.*");

  g_test_trap_subprocess ("/nc/data_cluster_wl_factor/shape_interp/file/rows/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the catalog has 3 rows, the dataset 6 galaxies.*");

  g_test_trap_subprocess ("/nc/data_cluster_wl_factor/shape_interp/file/create/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*cannot create *catalog.gvar.records*");
}

void
test_nc_data_cluster_wl_factor_shape_interp_file_shape_factor_subprocess (void)
{
  TestNcDataClusterWLFactorShapeInterp test;

  test_nc_data_cluster_wl_factor_shape_interp_new (&test, GINT_TO_POINTER (FALSE));
  nc_data_cluster_wl_factor_write_shape_interp_file (test.dcwlf, test.mset, test.obs, "catalog.gvar");
}

void
test_nc_data_cluster_wl_factor_shape_interp_file_rows_subprocess (void)
{
  TestNcDataClusterWLFactorShapeInterp test;
  NcGalaxyWLObs *obs;

  test_nc_data_cluster_wl_factor_shape_interp_new (&test, GINT_TO_POINTER (TRUE));
  obs = nc_galaxy_wl_obs_new (nc_galaxy_wl_obs_get_ellip_conv (test.obs), nc_galaxy_wl_obs_get_coord (test.obs), 3,
                              ncm_catalog_peek_columns (NCM_CATALOG (test.obs)), test_double_cols, g_strv_length (ncm_catalog_peek_columns (NCM_CATALOG (test.obs))));
  nc_data_cluster_wl_factor_write_shape_interp_file (test.dcwlf, test.mset, obs, "catalog.gvar");
}

void
test_nc_data_cluster_wl_factor_shape_interp_file_create_subprocess (void)
{
  TestNcDataClusterWLFactorShapeInterp test;
  gchar *filename = g_build_filename (g_get_tmp_dir (), "nc_shape_interp_no_such_directory", "catalog.gvar", NULL);

  test_nc_data_cluster_wl_factor_shape_interp_new (&test, GINT_TO_POINTER (TRUE));
  nc_data_cluster_wl_factor_write_shape_interp_file (test.dcwlf, test.mset, test.obs, filename);
}

