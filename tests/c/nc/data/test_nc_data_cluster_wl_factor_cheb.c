/***************************************************************************
 *            test_nc_data_cluster_wl_factor_cheb.c
 *
 *  Wed November 25 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_data_cluster_wl_factor_cheb.c
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
 * Tests of the CHEB_WEIGHTS integration method of NcDataClusterWLFactor that the truth test
 * does not cover: its properties, the statistic of the galaxies not converged, a change of
 * method on one object, and the grid shared between copies, on a catalog of tabulated p(z)
 * (NcGalaxyRedshiftFactorSpline, whose grids are shared).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <string.h>
#include <glib.h>
#include <glib-object.h>

/* Rows, lens redshift and seed of the catalog; each p(z) is a Gaussian of width
 * TEST_PZ_SIGMA tabulated every 0.01 on [TEST_PZ_Z0, TEST_PZ_Z1] */
#define TEST_NROWS (8)
#define TEST_Z_CL (0.5)
#define TEST_SEED (20261125)
#define TEST_PZ_SIGMA (0.08)
#define TEST_PZ_Z0 (0.05)
#define TEST_PZ_Z1 (2.0)

typedef struct _TestNcDataClusterWLFactorCheb
{
  NcmMSet *mset;
  NcDataClusterWLFactor *dcwlf;
  NcGalaxyWLObs *obs;
} TestNcDataClusterWLFactorCheb;

/* Half-width of the position box [deg] and the projected radius cut [Mpc] of the catalog */
typedef struct _TestNcDataClusterWLFactorChebCatalog
{
  gdouble box;
  gdouble min_radius;
  gdouble max_radius;
} TestNcDataClusterWLFactorChebCatalog;

/* Galaxies from 0.2 to 4 Mpc, where the shape factor is smooth in z */
static const TestNcDataClusterWLFactorChebCatalog test_catalog_default = {0.2, 0.2, 4.0};

/* Galaxies from 0.02 to 0.3 Mpc: at log10M 16 the convergence nears 1 and level 7 does not
 * resolve the shape factor in z */
static const TestNcDataClusterWLFactorChebCatalog test_catalog_center = {0.015, 0.02, 0.3};

void test_nc_data_cluster_wl_factor_cheb_new (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_free (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);

void test_nc_data_cluster_wl_factor_cheb_properties (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_converged (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_unconverged (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_switch (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_share (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_stored (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cheb_stored_errors (void);
void test_nc_data_cluster_wl_factor_cheb_stored_errors_subprocess (void);

void test_nc_data_cluster_wl_factor_cheb_level_min (void);
void test_nc_data_cluster_wl_factor_cheb_level_min_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_cluster_wl_factor/cheb/properties", TestNcDataClusterWLFactorCheb, &test_catalog_default,
              &test_nc_data_cluster_wl_factor_cheb_new,
              &test_nc_data_cluster_wl_factor_cheb_properties,
              &test_nc_data_cluster_wl_factor_cheb_free);
  g_test_add ("/nc/data_cluster_wl_factor/cheb/converged", TestNcDataClusterWLFactorCheb, &test_catalog_default,
              &test_nc_data_cluster_wl_factor_cheb_new,
              &test_nc_data_cluster_wl_factor_cheb_converged,
              &test_nc_data_cluster_wl_factor_cheb_free);
  g_test_add ("/nc/data_cluster_wl_factor/cheb/unconverged", TestNcDataClusterWLFactorCheb, &test_catalog_center,
              &test_nc_data_cluster_wl_factor_cheb_new,
              &test_nc_data_cluster_wl_factor_cheb_unconverged,
              &test_nc_data_cluster_wl_factor_cheb_free);
  g_test_add ("/nc/data_cluster_wl_factor/cheb/switch", TestNcDataClusterWLFactorCheb, &test_catalog_default,
              &test_nc_data_cluster_wl_factor_cheb_new,
              &test_nc_data_cluster_wl_factor_cheb_switch,
              &test_nc_data_cluster_wl_factor_cheb_free);
  g_test_add ("/nc/data_cluster_wl_factor/cheb/share", TestNcDataClusterWLFactorCheb, &test_catalog_default,
              &test_nc_data_cluster_wl_factor_cheb_new,
              &test_nc_data_cluster_wl_factor_cheb_share,
              &test_nc_data_cluster_wl_factor_cheb_free);

  g_test_add ("/nc/data_cluster_wl_factor/cheb/stored", TestNcDataClusterWLFactorCheb, &test_catalog_default,
              &test_nc_data_cluster_wl_factor_cheb_new,
              &test_nc_data_cluster_wl_factor_cheb_stored,
              &test_nc_data_cluster_wl_factor_cheb_free);
  g_test_add_func ("/nc/data_cluster_wl_factor/cheb/stored/errors", &test_nc_data_cluster_wl_factor_cheb_stored_errors);
  g_test_add_func ("/nc/data_cluster_wl_factor/cheb/stored/errors/subprocess", &test_nc_data_cluster_wl_factor_cheb_stored_errors_subprocess);
  g_test_add_func ("/nc/data_cluster_wl_factor/cheb/level_min", &test_nc_data_cluster_wl_factor_cheb_level_min);
  g_test_add_func ("/nc/data_cluster_wl_factor/cheb/level_min/subprocess", &test_nc_data_cluster_wl_factor_cheb_level_min_subprocess);

  g_test_run ();

  return 0;
}

/* The tabulated Gaussian of mean mu, the p(z) of one row */
static NcGalaxyPz *
_test_pz_new (const gdouble mu)
{
  const guint n = (guint) lround ((TEST_PZ_Z1 - TEST_PZ_Z0) / 0.01) + 1;
  NcmVector *zv = ncm_vector_new (n);
  NcmVector *pv = ncm_vector_new (n);
  NcGalaxyPz *pz;
  guint j;

  for (j = 0; j < n; j++)
  {
    const gdouble z = TEST_PZ_Z0 + 0.01 * j;

    ncm_vector_set (zv, j, z);
    ncm_vector_set (pv, j, exp (-0.5 * gsl_pow_2 ((z - mu) / TEST_PZ_SIGMA)));
  }

  pz = NC_GALAXY_PZ (nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv));

  ncm_vector_free (zv);
  ncm_vector_free (pv);

  return pz;
}

/* -2 ln P of every galaxy under @method */
static NcmVector *
_test_eval (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcDataClusterWLIntegMethod method)
{
  NcmVector *v = ncm_vector_new (TEST_NROWS);

  nc_data_cluster_wl_factor_set_integ_method (dcwlf, method);
  nc_data_cluster_wl_factor_eval_m2lnP_gal (dcwlf, mset, v);

  return v;
}

static void
_test_assert_equal (NcmVector *a, NcmVector *b)
{
  guint i;

  for (i = 0; i < TEST_NROWS; i++)
    g_assert_cmpfloat (ncm_vector_get (a, i), ==, ncm_vector_get (b, i));
}

/*
 * A catalog of TEST_NROWS galaxies around a lens at TEST_Z_CL, with p(z) means from 0.35 to
 * 1.2 (several straddle the lens), positions and radius cut from @pdata (a
 * TestNcDataClusterWLFactorChebCatalog), the VarAdd shape factor, and the dataset under
 * CHEB_WEIGHTS.
 */
void
test_nc_data_cluster_wl_factor_cheb_new (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  const TestNcDataClusterWLFactorChebCatalog *cat = pdata;
  NcmRNG *rng                                     = ncm_rng_seeded_new (NULL, TEST_SEED);
  NcHICosmo *cosmo                                = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                                = nc_distance_new (100.0);
  NcHaloMassSummary *hms                          = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp                        = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                              = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd                     = nc_wl_surface_mass_density_new (dist);
  NcGalaxyShapePop *pop_shape                     = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  const NcGalaxyWLObsEllipConv ec                 = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  const NcWLEllipticityFrame co                   = NC_WL_ELLIPTICITY_FRAME_CELESTIAL;
  const gdouble min_radius                        = cat->min_radius;
  const gdouble max_radius                        = cat->max_radius;
  NcGalaxyPositionFactor *position_factor         = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (-cat->box, cat->box, -cat->box, cat->box));
  NcGalaxyRedshiftFactor *redshift_factor         = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());
  NcGalaxyShapeFactor *shape_factor               = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (ec));
  GStrvBuilder *builder                           = g_strv_builder_new ();
  NcGalaxyPositionFactorData *p_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *s_data;
  GList *columns, *l;
  GStrv columns_strv;
  guint i;

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, 14.5);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_Z, TEST_Z_CL);
  nc_halo_position_prepare (hp, cosmo);
  ncm_model_param_set (NCM_MODEL (pop_shape), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.2);

  test->mset = ncm_mset_empty_new ();
  ncm_mset_set (test->mset, NCM_MODEL (cosmo), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (dp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (hp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (smd), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (pop_shape), NULL);

  p_data = nc_galaxy_position_factor_data_new (position_factor, test->mset);
  z_data = nc_galaxy_redshift_factor_data_new (redshift_factor, test->mset);
  s_data = nc_galaxy_shape_factor_data_new (shape_factor, test->mset, p_data, z_data);

  columns = nc_galaxy_shape_factor_data_required_columns (s_data);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  test->obs    = nc_galaxy_wl_obs_new (ec, co, TEST_NROWS, columns_strv);

  for (i = 0; i < TEST_NROWS; i++)
  {
    NcGalaxyPz *pz = _test_pz_new (0.35 + 0.85 * i / (TEST_NROWS - 1.0));
    gdouble radius;

    nc_galaxy_wl_obs_set_pz (test->obs, i, pz);
    nc_galaxy_redshift_factor_data_read_row (z_data, test->obs, i);
    nc_galaxy_redshift_factor_gen (redshift_factor, test->mset, z_data, rng);

    do {
      nc_galaxy_position_factor_gen (position_factor, test->mset, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    nc_galaxy_shape_factor_data_set (shape_factor, s_data, 0.0, 0.0, ncm_rng_uniform_gen (rng, 0.005, 0.02), 0.0, 0.0, 0.0, co);
    nc_galaxy_shape_factor_gen (shape_factor, test->mset, s_data, rng);
    nc_galaxy_shape_factor_data_write_row (s_data, test->obs, i);
    nc_galaxy_pz_free (pz);
  }

  test->dcwlf = nc_data_cluster_wl_factor_new (position_factor, redshift_factor, shape_factor);
  nc_data_cluster_wl_factor_set_cut (test->dcwlf, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (test->dcwlf, test->obs);
  nc_data_cluster_wl_factor_set_integ_method (test->dcwlf, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);

  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_mass_summary_free (hms);
  nc_halo_density_profile_free (dp);
  nc_halo_position_free (hp);
  nc_wl_surface_mass_density_free (smd);
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
test_nc_data_cluster_wl_factor_cheb_free (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  nc_galaxy_wl_obs_clear (&test->obs);
  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, test->dcwlf);
}

/* The defaults, the setters, the properties, and the description under CHEB_WEIGHTS */
void
test_nc_data_cluster_wl_factor_cheb_properties (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  gdouble tol;
  guint level;

  g_assert_cmpfloat (nc_data_cluster_wl_factor_get_cheb_tol (test->dcwlf), ==, 1.0e-8);
  g_assert_cmpuint (nc_data_cluster_wl_factor_get_cheb_level_min (test->dcwlf), ==, 2);
  g_assert_nonnull (strstr (ncm_data_get_desc (NCM_DATA (test->dcwlf)), "cheb_tol=1e-08, cheb_level_min=2"));

  nc_data_cluster_wl_factor_set_cheb_tol (test->dcwlf, 1.0e-10);
  nc_data_cluster_wl_factor_set_cheb_level_min (test->dcwlf, 5);
  g_assert_cmpfloat (nc_data_cluster_wl_factor_get_cheb_tol (test->dcwlf), ==, 1.0e-10);
  g_assert_cmpuint (nc_data_cluster_wl_factor_get_cheb_level_min (test->dcwlf), ==, 5);

  g_object_set (test->dcwlf, "cheb-tol", 1.0e-6, "cheb-level-min", 7, NULL);
  g_object_get (test->dcwlf, "cheb-tol", &tol, "cheb-level-min", &level, NULL);
  g_assert_cmpfloat (tol, ==, 1.0e-6);
  g_assert_cmpuint (level, ==, 7);
}

/* At the default tolerance every galaxy of the default catalog converges and the statistic is zero */
void
test_nc_data_cluster_wl_factor_cheb_converged (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  NcmVector *v = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);

  g_assert_cmpuint (nc_data_cluster_wl_factor_get_cheb_unconverged_count (test->dcwlf), ==, 0);
  g_assert_cmpfloat (nc_data_cluster_wl_factor_get_cheb_worst_relerr (test->dcwlf), ==, 0.0);
  ncm_vector_free (v);
}

/*
 * The central catalog evaluated at log10M 16, at the default tolerance: a galaxy reaches level 7
 * unconverged and is counted; the values stay finite. Its relative difference measures 1.7e-2,
 * far above the 1e-8 tolerance and the rounding of the sums; the bar is 1e-3.
 */
void
test_nc_data_cluster_wl_factor_cheb_unconverged (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcHaloMassSummary *hms   = nc_halo_density_profile_peek_mass_summary (dp);
  NcmVector *v;
  guint i;

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, 16.0);
  v = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  g_assert_cmpuint (nc_data_cluster_wl_factor_get_cheb_unconverged_count (test->dcwlf), >, 0);
  g_assert_cmpfloat (nc_data_cluster_wl_factor_get_cheb_worst_relerr (test->dcwlf), >, 1.0e-3);

  for (i = 0; i < TEST_NROWS; i++)
    g_assert_true (gsl_finite (ncm_vector_get (v, i)));

  ncm_vector_free (v);
}

/*
 * FIXED_NODES, CHEB_WEIGHTS and FIXED_NODES again on one object: a change of method rebuilds
 * the grid, so each method gives the same values whenever it runs, and the two methods differ
 * by their quadrature error only (FIXED_NODES 10 x 5 against the exact p-weighted levels).
 */
void
test_nc_data_cluster_wl_factor_cheb_switch (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  NcmVector *f1 = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES);
  NcmVector *c1 = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  NcmVector *f2 = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES);
  NcmVector *c2 = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  guint i;

  _test_assert_equal (f1, f2);
  _test_assert_equal (c1, c2);

  for (i = 0; i < TEST_NROWS; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (c1, i), ==, ncm_vector_get (f1, i), 1.0e-4, 0.0);

  ncm_vector_free (f1);
  ncm_vector_free (c1);
  ncm_vector_free (f2);
  ncm_vector_free (c2);
}

/*
 * A copy made through a serializer in which the dataset registered its catalog as shared adopts
 * the CHEB_WEIGHTS grid the dataset published and gives the same values. Switched to
 * FIXED_NODES, the copy does not take that grid: it builds its own and gives the values of the
 * dataset under FIXED_NODES.
 */
void
test_nc_data_cluster_wl_factor_cheb_share (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  NcmDataset *dset  = ncm_dataset_new ();
  NcDataClusterWLFactor *dup;
  NcmVector *c_orig, *c_dup, *f_dup, *f_orig;

  /* Registered before its first build, the dataset publishes the grid it builds */
  ncm_dataset_append_data (dset, NCM_DATA (test->dcwlf));
  ncm_dataset_register_shared (dset, ser);
  c_orig = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);

  dup = NC_DATA_CLUSTER_WL_FACTOR (ncm_serialize_dup_obj (ser, G_OBJECT (test->dcwlf)));
  g_assert_true (nc_data_cluster_wl_factor_peek_obs (dup) == nc_data_cluster_wl_factor_peek_obs (test->dcwlf));

  c_dup = _test_eval (dup, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  _test_assert_equal (c_dup, c_orig);

  f_dup  = _test_eval (dup, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES);
  f_orig = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES);
  _test_assert_equal (f_dup, f_orig);

  ncm_vector_free (c_orig);
  ncm_vector_free (c_dup);
  ncm_vector_free (f_dup);
  ncm_vector_free (f_orig);
  nc_data_cluster_wl_factor_free (dup);
  ncm_dataset_free (dset);
  ncm_serialize_free (ser);
}

void
test_nc_data_cluster_wl_factor_cheb_level_min (void)
{
  g_test_trap_subprocess ("/nc/data_cluster_wl_factor/cheb/level_min/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the level 8 is outside 2 to 7*");
}

void
test_nc_data_cluster_wl_factor_cheb_level_min_subprocess (void)
{
  NcGalaxyWLObsEllipConv ec = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  NcDataClusterWLFactor *dcwlf;

  dcwlf = nc_data_cluster_wl_factor_new (NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (-0.2, 0.2, -0.2, 0.2)),
                                         NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ()),
                                         NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (ec)));
  nc_data_cluster_wl_factor_set_cheb_level_min (dcwlf, 8);
}

/*
 * A dataset on a catalog that stores the grid takes it instead of building one: with the
 * grid intact its values equal those of the built grid, and with one normalization of the
 * stored grid changed its value for that galaxy changes.
 */
void
test_nc_data_cluster_wl_factor_cheb_stored (TestNcDataClusterWLFactorCheb *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmVector *built  = _test_eval (test->dcwlf, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  NcDataClusterWLFactor *dup;
  NcGalaxyWLObsChebGrid *grid;
  NcmVector *stored, *altered, *z_nodes = NULL;
  GArray *w = NULL;
  gdouble norm, bg_norm, fg_norm;
  guint n_total, i;

  nc_data_cluster_wl_factor_write_cheb_grid (test->dcwlf, test->mset, test->obs);
  grid = nc_galaxy_wl_obs_peek_cheb_grid (test->obs);
  g_assert_nonnull (grid);
  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_len (grid), ==, TEST_NROWS);
  g_assert_cmpfloat (nc_galaxy_wl_obs_cheb_grid_get_z_cl (grid), ==, TEST_Z_CL);
  g_assert_cmpstr (nc_galaxy_wl_obs_cheb_grid_peek_z_type (grid), ==, "NcGalaxyRedshiftFactorSpline");

  dup    = NC_DATA_CLUSTER_WL_FACTOR (ncm_serialize_dup_obj (ser, G_OBJECT (test->dcwlf)));
  stored = _test_eval (dup, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  _test_assert_equal (stored, built);
  nc_data_cluster_wl_factor_free (dup);

  /* A galaxy with nodes: its value is changed by a change of its stored normalization */
  for (i = 0; i < TEST_NROWS; i++)
  {
    n_total = nc_galaxy_wl_obs_cheb_grid_get_galaxy (grid, i, &w, &z_nodes, &norm, &bg_norm, &fg_norm);

    if (n_total > 0)
      break;

    ncm_vector_free (z_nodes);
  }

  g_assert_cmpuint (i, <, TEST_NROWS);
  nc_galaxy_wl_obs_cheb_grid_set_galaxy (grid, i, n_total, w, z_nodes, 2.0 * norm, bg_norm, fg_norm);
  ncm_serialize_reset (ser, TRUE);

  dup     = NC_DATA_CLUSTER_WL_FACTOR (ncm_serialize_dup_obj (ser, G_OBJECT (test->dcwlf)));
  altered = _test_eval (dup, test->mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
  g_assert_cmpfloat (ncm_vector_get (altered, i), !=, ncm_vector_get (built, i));

  g_array_unref (w);
  ncm_vector_free (z_nodes);
  ncm_vector_free (built);
  ncm_vector_free (stored);
  ncm_vector_free (altered);
  nc_data_cluster_wl_factor_free (dup);
  ncm_serialize_free (ser);
}

void
test_nc_data_cluster_wl_factor_cheb_stored_errors (void)
{
  const gchar *cases[][2] = {
    {"z_cl",   "*the Chebyshev grid of the catalog was built for the cluster redshift 0.5, the dataset has 0.6*"},
    {"method", "*nc_data_cluster_wl_factor_write_cheb_grid: the integration method is not CHEB_WEIGHTS*"},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    g_setenv ("TEST_NC_DATA_CLUSTER_WL_FACTOR_CHEB_ERROR", cases[i][0], TRUE);
    g_test_trap_subprocess ("/nc/data_cluster_wl_factor/cheb/stored/errors/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (cases[i][1]);
  }

  g_unsetenv ("TEST_NC_DATA_CLUSTER_WL_FACTOR_CHEB_ERROR");
}

void
test_nc_data_cluster_wl_factor_cheb_stored_errors_subprocess (void)
{
  const gchar *which = g_getenv ("TEST_NC_DATA_CLUSTER_WL_FACTOR_CHEB_ERROR");
  TestNcDataClusterWLFactorCheb test;

  test_nc_data_cluster_wl_factor_cheb_new (&test, &test_catalog_default);

  if (g_strcmp0 (which, "z_cl") == 0)
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    NcDataClusterWLFactor *dup;
    NcmVector *v;

    nc_data_cluster_wl_factor_set_integ_method (test.dcwlf, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
    nc_data_cluster_wl_factor_write_cheb_grid (test.dcwlf, test.mset, test.obs);
    dup = NC_DATA_CLUSTER_WL_FACTOR (ncm_serialize_dup_obj (ser, G_OBJECT (test.dcwlf)));
    ncm_model_param_set_by_name (ncm_mset_peek (test.mset, nc_halo_position_id ()), "z", 0.6, NULL);
    v = _test_eval (dup, test.mset, NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS);
    ncm_vector_free (v);
  }
  else
  {
    nc_data_cluster_wl_factor_set_integ_method (test.dcwlf, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES);
    nc_data_cluster_wl_factor_write_cheb_grid (test.dcwlf, test.mset, test.obs);
  }
}

