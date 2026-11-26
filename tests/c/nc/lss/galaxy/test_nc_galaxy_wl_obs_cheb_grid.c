/***************************************************************************
 *            test_nc_galaxy_wl_obs_cheb_grid.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_wl_obs_cheb_grid.c
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

#include <glib.h>
#include <glib-object.h>

/* A small grid: 3 galaxies, 4 weights and at most 3 redshift values per galaxy */
#define TEST_LEN (3)
#define TEST_N_WEIGHTS (4)
#define TEST_N_NODES (3)
#define TEST_Z_CL (0.25)
#define TEST_Z_TYPE "NcGalaxyRedshiftFactorSpline"
#define TEST_HASH "0123abcd"

typedef struct _TestNcGalaxyWLObsChebGrid
{
  NcGalaxyWLObsChebGrid *grid;
} TestNcGalaxyWLObsChebGrid;

void test_nc_galaxy_wl_obs_cheb_grid_new (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata);
void test_nc_galaxy_wl_obs_cheb_grid_free (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata);

void test_nc_galaxy_wl_obs_cheb_grid_record (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata);
void test_nc_galaxy_wl_obs_cheb_grid_galaxy (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata);
void test_nc_galaxy_wl_obs_cheb_grid_serialize (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata);
void test_nc_galaxy_wl_obs_cheb_grid_errors (void);
void test_nc_galaxy_wl_obs_cheb_grid_errors_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_wl_obs_cheb_grid/record", TestNcGalaxyWLObsChebGrid, NULL,
              &test_nc_galaxy_wl_obs_cheb_grid_new,
              &test_nc_galaxy_wl_obs_cheb_grid_record,
              &test_nc_galaxy_wl_obs_cheb_grid_free);

  g_test_add ("/nc/galaxy_wl_obs_cheb_grid/galaxy", TestNcGalaxyWLObsChebGrid, NULL,
              &test_nc_galaxy_wl_obs_cheb_grid_new,
              &test_nc_galaxy_wl_obs_cheb_grid_galaxy,
              &test_nc_galaxy_wl_obs_cheb_grid_free);

  g_test_add ("/nc/galaxy_wl_obs_cheb_grid/serialize", TestNcGalaxyWLObsChebGrid, NULL,
              &test_nc_galaxy_wl_obs_cheb_grid_new,
              &test_nc_galaxy_wl_obs_cheb_grid_serialize,
              &test_nc_galaxy_wl_obs_cheb_grid_free);

  g_test_add_func ("/nc/galaxy_wl_obs_cheb_grid/errors", &test_nc_galaxy_wl_obs_cheb_grid_errors);
  g_test_add_func ("/nc/galaxy_wl_obs_cheb_grid/errors/subprocess", &test_nc_galaxy_wl_obs_cheb_grid_errors_subprocess);

  g_test_run ();

  return 0;
}

/*
 * Fills galaxy 0 with two nodes and galaxy 2 with none (no background support); galaxy 1
 * keeps the zero values of a new grid.
 */
static void
_test_nc_galaxy_wl_obs_cheb_grid_fill (NcGalaxyWLObsChebGrid *grid)
{
  const gdouble w0[TEST_N_WEIGHTS] = {0.5, -0.25, 0.125, 2.0};
  GArray *weights                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *z0                    = ncm_vector_new (3);
  NcmVector *z2                    = ncm_vector_new (1);

  g_array_append_vals (weights, w0, TEST_N_WEIGHTS);
  ncm_vector_set (z0, 0, 0.1);
  ncm_vector_set (z0, 1, 0.4);
  ncm_vector_set (z0, 2, 1.2);
  ncm_vector_set (z2, 0, 0.05);

  nc_galaxy_wl_obs_cheb_grid_set_galaxy (grid, 0, 2, weights, z0, 1.0, 0.75, 0.25);
  nc_galaxy_wl_obs_cheb_grid_set_galaxy (grid, 2, 0, NULL, z2, 0.9, 0.0, 0.0);

  g_array_unref (weights);
  ncm_vector_free (z0);
  ncm_vector_free (z2);
}

static void
_test_nc_galaxy_wl_obs_cheb_grid_assert_filled (NcGalaxyWLObsChebGrid *grid)
{
  GArray *weights    = NULL;
  NcmVector *z_nodes = NULL;
  gdouble norm, bg_norm, fg_norm;

  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_get_galaxy (grid, 0, &weights, &z_nodes, &norm, &bg_norm, &fg_norm), ==, 2);
  g_assert_nonnull (weights);
  g_assert_cmpuint (weights->len, ==, TEST_N_WEIGHTS);
  g_assert_cmpfloat (g_array_index (weights, gdouble, 1), ==, -0.25);
  g_assert_cmpfloat (g_array_index (weights, gdouble, 3), ==, 2.0);
  g_assert_cmpuint (ncm_vector_len (z_nodes), ==, 3);
  g_assert_cmpfloat (ncm_vector_get (z_nodes, 0), ==, 0.1);
  g_assert_cmpfloat (ncm_vector_get (z_nodes, 2), ==, 1.2);
  g_assert_cmpfloat (norm, ==, 1.0);
  g_assert_cmpfloat (bg_norm, ==, 0.75);
  g_assert_cmpfloat (fg_norm, ==, 0.25);
  g_array_unref (weights);
  ncm_vector_free (z_nodes);

  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_get_galaxy (grid, 2, &weights, &z_nodes, &norm, &bg_norm, &fg_norm), ==, 0);
  g_assert_null (weights);
  g_assert_cmpuint (ncm_vector_len (z_nodes), ==, 1);
  g_assert_cmpfloat (ncm_vector_get (z_nodes, 0), ==, 0.05);
  g_assert_cmpfloat (norm, ==, 0.9);
  ncm_vector_free (z_nodes);
}

void
test_nc_galaxy_wl_obs_cheb_grid_new (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata)
{
  test->grid = nc_galaxy_wl_obs_cheb_grid_new (TEST_LEN, TEST_N_WEIGHTS, TEST_N_NODES, TEST_Z_CL, TEST_Z_TYPE, TEST_HASH);

  g_assert_true (NC_IS_GALAXY_WL_OBS_CHEB_GRID (test->grid));
}

void
test_nc_galaxy_wl_obs_cheb_grid_free (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_wl_obs_cheb_grid_free, test->grid);
}

void
test_nc_galaxy_wl_obs_cheb_grid_record (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata)
{
  NcGalaxyWLObsChebGrid *ref = nc_galaxy_wl_obs_cheb_grid_ref (test->grid);

  g_assert_true (ref == test->grid);
  nc_galaxy_wl_obs_cheb_grid_clear (&ref);
  g_assert_null (ref);

  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_len (test->grid), ==, TEST_LEN);
  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_get_n_weights (test->grid), ==, TEST_N_WEIGHTS);
  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_get_n_nodes (test->grid), ==, TEST_N_NODES);
  g_assert_cmpfloat (nc_galaxy_wl_obs_cheb_grid_get_z_cl (test->grid), ==, TEST_Z_CL);
  g_assert_cmpstr (nc_galaxy_wl_obs_cheb_grid_peek_z_type (test->grid), ==, TEST_Z_TYPE);
  g_assert_cmpstr (nc_galaxy_wl_obs_cheb_grid_peek_data_hash (test->grid), ==, TEST_HASH);
}

/*
 * A galaxy's grid is read back as it was set; a new grid has no node, only the zero lower
 * end, and every normalization zero.
 */
void
test_nc_galaxy_wl_obs_cheb_grid_galaxy (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata)
{
  GArray *weights    = NULL;
  NcmVector *z_nodes = NULL;
  gdouble norm, bg_norm, fg_norm;

  _test_nc_galaxy_wl_obs_cheb_grid_fill (test->grid);
  _test_nc_galaxy_wl_obs_cheb_grid_assert_filled (test->grid);

  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_get_galaxy (test->grid, 1, &weights, &z_nodes, &norm, &bg_norm, &fg_norm), ==, 0);
  g_assert_null (weights);
  g_assert_cmpuint (ncm_vector_len (z_nodes), ==, 1);
  g_assert_cmpfloat (ncm_vector_get (z_nodes, 0), ==, 0.0);
  g_assert_cmpfloat (norm, ==, 0.0);
  ncm_vector_free (z_nodes);
}

void
test_nc_galaxy_wl_obs_cheb_grid_serialize (TestNcGalaxyWLObsChebGrid *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObsChebGrid *dup;
  GVariant *var;

  _test_nc_galaxy_wl_obs_cheb_grid_fill (test->grid);

  var = ncm_serialize_to_variant (ser, G_OBJECT (test->grid));
  dup = NC_GALAXY_WL_OBS_CHEB_GRID (ncm_serialize_from_variant (ser, var));

  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_len (dup), ==, TEST_LEN);
  g_assert_cmpfloat (nc_galaxy_wl_obs_cheb_grid_get_z_cl (dup), ==, TEST_Z_CL);
  g_assert_cmpstr (nc_galaxy_wl_obs_cheb_grid_peek_z_type (dup), ==, TEST_Z_TYPE);
  g_assert_cmpstr (nc_galaxy_wl_obs_cheb_grid_peek_data_hash (dup), ==, TEST_HASH);
  _test_nc_galaxy_wl_obs_cheb_grid_assert_filled (dup);

  g_variant_unref (var);
  ncm_serialize_free (ser);
  NCM_TEST_FREE (nc_galaxy_wl_obs_cheb_grid_free, dup);
}

void
test_nc_galaxy_wl_obs_cheb_grid_errors (void)
{
  const gchar *cases[][2] = {
    {"no_nodes",    "*a grid needs at least one redshift value per galaxy*"},
    {"set_index",   "*nc_galaxy_wl_obs_cheb_grid_set_galaxy: galaxy 3 of 3*"},
    {"get_index",   "*nc_galaxy_wl_obs_cheb_grid_get_galaxy: galaxy 3 of 3*"},
    {"z_len",       "*2 redshift values for 2 nodes and at most 3 values*"},
    {"no_weights",  "*nc_galaxy_wl_obs_cheb_grid_set_galaxy: 0 weights, not 4*"},
    {"matrix",      "*the weights matrix is 3 x 5, not 3 x 4*"},
    {"n_total_len", "*2 node counts for 3 galaxies*"},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    g_setenv ("TEST_NC_GALAXY_WL_OBS_CHEB_GRID_ERROR", cases[i][0], TRUE);
    g_test_trap_subprocess ("/nc/galaxy_wl_obs_cheb_grid/errors/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (cases[i][1]);
  }

  g_unsetenv ("TEST_NC_GALAXY_WL_OBS_CHEB_GRID_ERROR");
}

void
test_nc_galaxy_wl_obs_cheb_grid_errors_subprocess (void)
{
  const gchar *which          = g_getenv ("TEST_NC_GALAXY_WL_OBS_CHEB_GRID_ERROR");
  NcGalaxyWLObsChebGrid *grid = NULL;
  NcmVector *z                = ncm_vector_new (2);
  GArray *w                   = NULL;
  gdouble norm, bg_norm, fg_norm;

  if (g_strcmp0 (which, "no_nodes") == 0)
  {
    nc_galaxy_wl_obs_cheb_grid_new (TEST_LEN, TEST_N_WEIGHTS, 0, TEST_Z_CL, TEST_Z_TYPE, TEST_HASH);

    return;
  }

  grid = nc_galaxy_wl_obs_cheb_grid_new (TEST_LEN, TEST_N_WEIGHTS, TEST_N_NODES, TEST_Z_CL, TEST_Z_TYPE, TEST_HASH);

  if (g_strcmp0 (which, "set_index") == 0)
  {
    nc_galaxy_wl_obs_cheb_grid_set_galaxy (grid, TEST_LEN, 1, NULL, z, 1.0, 1.0, 0.0);
  }
  else if (g_strcmp0 (which, "get_index") == 0)
  {
    nc_galaxy_wl_obs_cheb_grid_get_galaxy (grid, TEST_LEN, &w, &z, &norm, &bg_norm, &fg_norm);
  }
  else if (g_strcmp0 (which, "z_len") == 0)
  {
    nc_galaxy_wl_obs_cheb_grid_set_galaxy (grid, 0, 2, NULL, z, 1.0, 1.0, 0.0);
  }
  else if (g_strcmp0 (which, "no_weights") == 0)
  {
    nc_galaxy_wl_obs_cheb_grid_set_galaxy (grid, 0, 1, NULL, z, 1.0, 1.0, 0.0);
  }
  else if (g_strcmp0 (which, "matrix") == 0)
  {
    g_object_set (grid, "weights", ncm_matrix_new (TEST_LEN, TEST_N_WEIGHTS + 1), NULL);
  }
  else
  {
    const guint32 n_total[2] = {0, 0};

    g_object_set (grid, "n-total", g_variant_new_fixed_array (G_VARIANT_TYPE_UINT32, n_total, 2, sizeof (guint32)), NULL);
  }
}

