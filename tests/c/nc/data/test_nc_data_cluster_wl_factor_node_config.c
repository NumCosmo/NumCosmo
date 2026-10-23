/***************************************************************************
 *            test_nc_data_cluster_wl_factor_node_config.c
 *
 *  Fri Oct 23 2026
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
 * Tests of the NcDataClusterWLFactor:node-config property under FIXED_NODES:
 * the stored settings are the ones the per-galaxy pairs were calibrated
 * under, a new catalog drops the pairs, and the entries no calibration wrote
 * are zero.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define NODE_CONFIG_NROWS (8)
#define NODE_CONFIG_Z_CL (0.5)
#define NODE_CONFIG_SEED (20261023)
#define NODE_CONFIG_RELTOL_CALIB (1.0e-2)
#define NODE_CONFIG_RELTOL_LATER (1.0e-4)

typedef struct _TestNcDataClusterWLFactorNodeConfig
{
  NcmMSet *mset;
  NcDataClusterWLFactor *dcwlf;
  NcGalaxyWLObs *obs;
} TestNcDataClusterWLFactorNodeConfig;

void test_nc_data_cluster_wl_factor_node_config_new (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_node_config_free (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_node_config_calibrated_settings (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_node_config_set_obs (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_node_config_zero_fill (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/data_cluster_wl_factor/node_config/calibrated_settings", TestNcDataClusterWLFactorNodeConfig, NULL,
              &test_nc_data_cluster_wl_factor_node_config_new,
              &test_nc_data_cluster_wl_factor_node_config_calibrated_settings,
              &test_nc_data_cluster_wl_factor_node_config_free);

  g_test_add ("/nc/data_cluster_wl_factor/node_config/set_obs", TestNcDataClusterWLFactorNodeConfig, NULL,
              &test_nc_data_cluster_wl_factor_node_config_new,
              &test_nc_data_cluster_wl_factor_node_config_set_obs,
              &test_nc_data_cluster_wl_factor_node_config_free);

  g_test_add ("/nc/data_cluster_wl_factor/node_config/zero_fill", TestNcDataClusterWLFactorNodeConfig, NULL,
              &test_nc_data_cluster_wl_factor_node_config_new,
              &test_nc_data_cluster_wl_factor_node_config_zero_fill,
              &test_nc_data_cluster_wl_factor_node_config_free);

  g_test_run ();

  return 0;
}

/* Returns the node-config property (NULL when there is none) and, when
 * @n_nodes / @rule_n are not NULL, its per-galaxy arrays. */
static NcmVarDict *
_node_config_peek (NcDataClusterWLFactor *dcwlf, GArray **n_nodes, GArray **rule_n)
{
  NcmVarDict *vd = NULL;

  g_object_get (dcwlf, "node-config", &vd, NULL);

  if (vd != NULL)
  {
    if (n_nodes != NULL)
      g_assert_true (ncm_var_dict_get_int_array (vd, "n-nodes", n_nodes));

    if (rule_n != NULL)
      g_assert_true (ncm_var_dict_get_int_array (vd, "rule-n", rule_n));
  }

  return vd;
}

void
test_nc_data_cluster_wl_factor_node_config_new (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata)
{
  NcmRNG *rng                             = ncm_rng_seeded_new (NULL, NODE_CONFIG_SEED);
  NcHICosmo *cosmo                        = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                        = nc_distance_new (100.0);
  NcHaloMassSummary *hms                  = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp                = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                      = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd             = nc_wl_surface_mass_density_new (dist);
  NcGalaxyRedshiftPopLSSTSRD *pop_z       = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z         = nc_galaxy_redshift_obs_gauss_new ();
  NcGalaxyShapePop *pop_shape             = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  const NcGalaxyWLObsEllipConv ec         = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  const NcWLEllipticityFrame co           = NC_WL_ELLIPTICITY_FRAME_CELESTIAL;
  const gdouble min_radius                = 0.2;
  const gdouble max_radius                = 4.0;
  NcGalaxyPositionFactor *position_factor = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcGalaxyRedshiftFactor *redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0));
  NcGalaxyShapeFactor *shape_factor       = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (ec));
  GStrvBuilder *builder                   = g_strv_builder_new ();
  NcGalaxyPositionFactorData *p_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *s_data;
  GList *columns, *l;
  GStrv columns_strv;
  guint i;

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, 14.5);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_Z, NODE_CONFIG_Z_CL);
  nc_halo_position_prepare (hp, cosmo);
  ncm_model_param_set (NCM_MODEL (pop_shape), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.2);

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
  test->obs    = nc_galaxy_wl_obs_new (ec, co, NODE_CONFIG_NROWS, columns_strv);

  for (i = 0; i < NODE_CONFIG_NROWS; i++)
  {
    gdouble radius;

    /* Composed's gen() draws zp from the row's sigma0, which has to be on
     * the catalog first (see test_nc_data_cluster_wl_factor_truth.c). */
    nc_galaxy_wl_obs_set (test->obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, ncm_rng_uniform_gen (rng, 0.02, 0.06), NULL);
    nc_galaxy_redshift_factor_data_read_row (z_data, test->obs, i);
    nc_galaxy_redshift_factor_gen (redshift_factor, test->mset, z_data, rng);

    do {
      nc_galaxy_position_factor_gen (position_factor, test->mset, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    nc_galaxy_shape_factor_data_set (shape_factor, s_data, 0.0, 0.0, ncm_rng_uniform_gen (rng, 0.005, 0.02), 0.0, 0.0, 0.0, co);
    nc_galaxy_shape_factor_gen (shape_factor, test->mset, s_data, rng);
    nc_galaxy_shape_factor_data_write_row (s_data, test->obs, i);
  }

  test->dcwlf = nc_data_cluster_wl_factor_new (position_factor, redshift_factor, shape_factor);
  nc_data_cluster_wl_factor_set_cut (test->dcwlf, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (test->dcwlf, test->obs);
  nc_data_cluster_wl_factor_set_integ_method (test->dcwlf, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES);
  nc_data_cluster_wl_factor_set_auto_nodes (test->dcwlf, TRUE);
  nc_data_cluster_wl_factor_set_node_reltol (test->dcwlf, NODE_CONFIG_RELTOL_CALIB);

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
test_nc_data_cluster_wl_factor_node_config_free (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  nc_galaxy_wl_obs_clear (&test->obs);
  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, test->dcwlf);
}

/* A node-reltol change after the calibration does not change the settings
 * the node-config states; a copy carrying that config and the new
 * node-reltol discards it and calibrates again. */
void
test_nc_data_cluster_wl_factor_node_config_calibrated_settings (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcDataClusterWLFactor *dup;
  NcmVarDict *vd;
  gdouble reltol = 0.0;

  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);
  g_assert_cmpuint (nc_data_cluster_wl_factor_get_calib_count (test->dcwlf), >, 0);

  nc_data_cluster_wl_factor_set_node_reltol (test->dcwlf, NODE_CONFIG_RELTOL_LATER);

  vd = _node_config_peek (test->dcwlf, NULL, NULL);
  g_assert_nonnull (vd);
  g_assert_true (ncm_var_dict_get_double (vd, "node-reltol", &reltol));
  g_assert_cmpfloat (reltol, ==, NODE_CONFIG_RELTOL_CALIB);
  ncm_var_dict_unref (vd);

  dup = NC_DATA_CLUSTER_WL_FACTOR (ncm_serialize_dup_obj (ser, G_OBJECT (test->dcwlf)));
  g_assert_cmpfloat (nc_data_cluster_wl_factor_get_node_reltol (dup), ==, NODE_CONFIG_RELTOL_LATER);

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*stored node configuration does not match*");
  nc_data_cluster_wl_factor_data_prepare (dup, test->mset);
  g_test_assert_expected_messages ();

  g_assert_cmpuint (nc_data_cluster_wl_factor_get_calib_count (dup), >, 0);

  vd = _node_config_peek (dup, NULL, NULL);
  g_assert_nonnull (vd);
  g_assert_true (ncm_var_dict_get_double (vd, "node-reltol", &reltol));
  g_assert_cmpfloat (reltol, ==, NODE_CONFIG_RELTOL_LATER);
  ncm_var_dict_unref (vd);

  nc_data_cluster_wl_factor_free (dup);
  ncm_serialize_free (ser);
}

/* A new catalog drops the calibrated pairs: until the next prepare the
 * dataset is not prepared and reports no node-config. */
void
test_nc_data_cluster_wl_factor_node_config_set_obs (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata)
{
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcGalaxyWLObs *copy = NC_GALAXY_WL_OBS (ncm_serialize_dup_obj (ser, G_OBJECT (test->obs)));
  NcmVarDict *vd;

  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);
  g_assert_true (nc_data_cluster_wl_factor_is_data_prepared (test->dcwlf));

  nc_data_cluster_wl_factor_set_obs (test->dcwlf, copy);
  g_assert_false (nc_data_cluster_wl_factor_is_data_prepared (test->dcwlf));

  vd = _node_config_peek (test->dcwlf, NULL, NULL);
  g_assert_null (vd);

  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);
  g_assert_true (nc_data_cluster_wl_factor_is_data_prepared (test->dcwlf));

  nc_galaxy_wl_obs_free (copy);
  ncm_serialize_free (ser);
}

/* With auto-nodes off no galaxy is calibrated, so every node-config entry is
 * zero. The first prepare with auto-nodes on fills the arrays with nonzero
 * pairs, so entries the second prepare left unwritten would show them. */
void
test_nc_data_cluster_wl_factor_node_config_zero_fill (TestNcDataClusterWLFactorNodeConfig *test, gconstpointer pdata)
{
  GArray *n_nodes = NULL;
  GArray *rule_n  = NULL;
  gboolean any_nonzero;
  NcmVarDict *vd;
  guint i;

  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);

  vd = _node_config_peek (test->dcwlf, &n_nodes, NULL);
  g_assert_nonnull (vd);
  g_assert_cmpuint (n_nodes->len, ==, NODE_CONFIG_NROWS);

  any_nonzero = FALSE;

  for (i = 0; i < n_nodes->len; i++)
    any_nonzero = any_nonzero || (g_array_index (n_nodes, gint, i) != 0);

  g_assert_true (any_nonzero);
  g_array_unref (n_nodes);
  ncm_var_dict_unref (vd);

  nc_data_cluster_wl_factor_set_auto_nodes (test->dcwlf, FALSE);
  nc_data_cluster_wl_factor_data_prepare (test->dcwlf, test->mset);

  vd = _node_config_peek (test->dcwlf, &n_nodes, &rule_n);
  g_assert_nonnull (vd);
  g_assert_cmpuint (n_nodes->len, ==, NODE_CONFIG_NROWS);
  g_assert_cmpuint (rule_n->len, ==, NODE_CONFIG_NROWS);

  for (i = 0; i < NODE_CONFIG_NROWS; i++)
  {
    g_assert_cmpint (g_array_index (n_nodes, gint, i), ==, 0);
    g_assert_cmpint (g_array_index (rule_n, gint, i), ==, 0);
  }

  g_array_unref (n_nodes);
  g_array_unref (rule_n);
  ncm_var_dict_unref (vd);
}

