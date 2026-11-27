/***************************************************************************
 *            test_nc_data_cluster_wl_factor_cache.c
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_data_cluster_wl_factor_cache.c
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
 * Regression tests of the change values that gate the per-galaxy updates of
 * #NcDataClusterWLFactor. prepare() recomputes the projected radius and the rotation factor
 * of a galaxy when the radius value of #NcGalaxyHaloWL changes (the halo position and the
 * cosmology), its surface densities when the profile value of #NcWLClusterLens changes (the
 * cosmology, the density profile, its mass summary and the halo position), its critical
 * surface densities when the geometry value of the lens changes (the cosmology and the
 * cluster redshift), and its population fragment when the population value of the shape
 * factor changes. A submodel's pkey does not propagate to its parent; the #NcmModelCtrl of
 * the profile records its mass summary, so a mass change counts as a profile change.
 *
 * One dataset, reused over a sequence of model values that returns to earlier values, is
 * compared with a dataset built from new models at every step, which holds no cached
 * state. The two agree at every step when no change value misses a change of its models.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <string.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_N_GAL (8)
#define TEST_SEED (20261127)
#define TEST_ELLIP_CONV (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET)

/* Every column of these catalogs is DOUBLE (NCM_CATALOG_COL_TYPE_DOUBLE is zero) */
static const NcmCatalogColType test_double_cols[64] = {NCM_CATALOG_COL_TYPE_DOUBLE};

typedef struct _TestNcDataClusterWLFactorCacheModels
{
  NcmMSet *mset;
  NcHICosmo *cosmo;
  NcHaloMassSummary *hms;
  NcHaloPosition *hp;
  NcWLSurfaceMassDensity *smd;
  NcGalaxyShapePop *pop;
} TestNcDataClusterWLFactorCacheModels;

typedef struct _TestNcDataClusterWLFactorCache
{
  TestNcDataClusterWLFactorCacheModels models;
} TestNcDataClusterWLFactorCache;

/* One model state: log10MDelta, cDelta, the cluster redshift and the population sigma */
typedef struct _TestNcDataClusterWLFactorCacheConfig
{
  gdouble log10_mdelta;
  gdouble c_delta;
  gdouble z_cl;
  gdouble pop_sigma;
} TestNcDataClusterWLFactorCacheConfig;

/*
 * Each of the four values returns to an earlier value at least once: a change value that
 * misses a model gives the cached values of the earlier state on a return, which a
 * monotonic sweep does not show.
 */
static const TestNcDataClusterWLFactorCacheConfig test_sequence[] = {
  {14.0, 4.0, 0.2, 0.25},
  {14.5, 4.0, 0.2, 0.25},  /* mass */
  {14.0, 4.0, 0.2, 0.25},  /* mass returns */
  {14.0, 5.0, 0.2, 0.25},  /* concentration */
  {14.0, 4.0, 0.2, 0.25},  /* concentration returns */
  {14.0, 4.0, 0.35, 0.25}, /* cluster redshift */
  {14.0, 4.0, 0.2, 0.25},  /* cluster redshift returns */
  {14.0, 4.0, 0.2, 0.35},  /* population sigma */
  {14.0, 4.0, 0.2, 0.25},  /* population sigma returns */
};

static const NcDataClusterWLIntegMethod test_method_lnint = NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT;
static const NcDataClusterWLIntegMethod test_method_cheb  = NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS;

void test_nc_data_cluster_wl_factor_cache_new (TestNcDataClusterWLFactorCache *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cache_free (TestNcDataClusterWLFactorCache *test, gconstpointer pdata);

void test_nc_data_cluster_wl_factor_cache_revisits (TestNcDataClusterWLFactorCache *test, gconstpointer pdata);
void test_nc_data_cluster_wl_factor_cache_change_values (TestNcDataClusterWLFactorCache *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_cluster_wl_factor/cache/revisits/lnint", TestNcDataClusterWLFactorCache, &test_method_lnint,
              &test_nc_data_cluster_wl_factor_cache_new,
              &test_nc_data_cluster_wl_factor_cache_revisits,
              &test_nc_data_cluster_wl_factor_cache_free);
  g_test_add ("/nc/data_cluster_wl_factor/cache/revisits/cheb", TestNcDataClusterWLFactorCache, &test_method_cheb,
              &test_nc_data_cluster_wl_factor_cache_new,
              &test_nc_data_cluster_wl_factor_cache_revisits,
              &test_nc_data_cluster_wl_factor_cache_free);
  g_test_add ("/nc/data_cluster_wl_factor/cache/change_values", TestNcDataClusterWLFactorCache, NULL,
              &test_nc_data_cluster_wl_factor_cache_new,
              &test_nc_data_cluster_wl_factor_cache_change_values,
              &test_nc_data_cluster_wl_factor_cache_free);

  g_test_run ();

  return 0;
}

/* The models of the mset, built once: a fit changes the parameter values of the same models */
static void
_test_models_init (TestNcDataClusterWLFactorCacheModels *m)
{
  NcDistance *dist                  = nc_distance_new (100.0);
  NcHaloDensityProfile *dp          = NULL;
  NcGalaxyRedshiftPopLSSTSRD *pop_z = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z   = nc_galaxy_redshift_obs_gauss_new ();

  m->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  m->hms   = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  dp       = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (m->hms));
  m->hp    = nc_halo_position_new (dist);
  m->smd   = nc_wl_surface_mass_density_new (dist);
  m->pop   = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  m->mset  = ncm_mset_empty_new ();

  ncm_mset_set (m->mset, NCM_MODEL (m->cosmo), NULL);
  ncm_mset_set (m->mset, NCM_MODEL (dp), NULL);
  ncm_mset_set (m->mset, NCM_MODEL (m->hp), NULL);
  ncm_mset_set (m->mset, NCM_MODEL (m->smd), NULL);
  ncm_mset_set (m->mset, NCM_MODEL (m->pop), NULL);
  ncm_mset_set (m->mset, NCM_MODEL (pop_z), NULL);
  ncm_mset_set (m->mset, NCM_MODEL (obs_z), NULL);

  nc_halo_density_profile_free (dp);
  nc_galaxy_redshift_pop_lsst_srd_free (pop_z);
  nc_galaxy_redshift_obs_gauss_free (obs_z);
  nc_distance_free (dist);
}

static void
_test_models_clear (TestNcDataClusterWLFactorCacheModels *m)
{
  ncm_mset_clear (&m->mset);
  nc_hicosmo_clear (&m->cosmo);
  nc_halo_mass_summary_clear (&m->hms);
  nc_halo_position_clear (&m->hp);
  nc_wl_surface_mass_density_clear (&m->smd);
  nc_galaxy_shape_pop_clear (&m->pop);
}

/*
 * Sets @config, touching only the parameters that differ from @prev (NULL sets all).
 * ncm_model_param_set() changes the pkey on every call, also when it sets the value the
 * parameter already has, so setting all four at every step would change every value at
 * every step and a value that misses one of its models would pass.
 */
static void
_test_models_set (TestNcDataClusterWLFactorCacheModels *m, const TestNcDataClusterWLFactorCacheConfig *config, const TestNcDataClusterWLFactorCacheConfig *prev)
{
  if ((prev == NULL) || (config->log10_mdelta != prev->log10_mdelta))
    ncm_model_param_set_by_name (NCM_MODEL (m->hms), "log10MDelta", config->log10_mdelta, NULL);

  if ((prev == NULL) || (config->c_delta != prev->c_delta))
    ncm_model_param_set_by_name (NCM_MODEL (m->hms), "cDelta", config->c_delta, NULL);

  if ((prev == NULL) || (config->z_cl != prev->z_cl))
  {
    ncm_model_param_set_by_name (NCM_MODEL (m->hp), "z", config->z_cl, NULL);
    nc_halo_position_prepare (m->hp, m->cosmo);
  }

  if ((prev == NULL) || (config->pop_sigma != prev->pop_sigma))
    ncm_model_param_set_by_name (NCM_MODEL (m->pop), "sigma", config->pop_sigma, NULL);
}

/* The dataset of TEST_N_GAL galaxies drawn from TEST_SEED, the same for every call */
static NcDataClusterWLFactor *
_test_dataset_new (NcmMSet *mset, NcDataClusterWLIntegMethod method)
{
  NcGalaxyPositionFactor *posf = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcGalaxyRedshiftFactor *zf   = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0));
  NcGalaxyShapeFactor *gsf     = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (TEST_ELLIP_CONV));
  NcDataClusterWLFactor *dcwlf = nc_data_cluster_wl_factor_new (posf, zf, gsf);
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, TEST_SEED);
  GList *columns               = nc_data_cluster_wl_factor_required_columns (dcwlf, mset);
  GStrvBuilder *builder        = g_strv_builder_new ();
  NcGalaxyWLObs *obs;
  GStrv columns_strv;
  GList *l;
  guint i;

  for (l = columns; l != NULL; l = l->next)
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  obs          = nc_galaxy_wl_obs_new (TEST_ELLIP_CONV, NC_WL_ELLIPTICITY_FRAME_CELESTIAL, TEST_N_GAL, columns_strv, test_double_cols, g_strv_length (columns_strv));

  for (i = 0; i < TEST_N_GAL; i++)
  {
    NcmCatalog *cat = NCM_CATALOG (obs);

    for (l = columns; l != NULL; l = l->next)
      ncm_catalog_set (cat, l->data, i, 0.0, NULL);

    ncm_catalog_set (cat, "ra", i, ncm_rng_uniform_gen (rng, -0.15, 0.15), NULL);
    ncm_catalog_set (cat, "dec", i, ncm_rng_uniform_gen (rng, -0.15, 0.15), NULL);
    ncm_catalog_set (cat, "zp", i, ncm_rng_uniform_gen (rng, 0.3, 1.5), NULL);
    ncm_catalog_set (cat, "sigma0", i, ncm_rng_uniform_gen (rng, 0.02, 0.06), NULL);
    ncm_catalog_set (cat, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1, i, ncm_rng_gaussian_gen (rng, 0.0, 0.2), NULL);
    ncm_catalog_set (cat, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2, i, ncm_rng_gaussian_gen (rng, 0.0, 0.2), NULL);
    ncm_catalog_set (cat, NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE, i, ncm_rng_uniform_gen (rng, 0.03, 0.08), NULL);
  }

  nc_data_cluster_wl_factor_set_obs (dcwlf, obs);
  nc_data_cluster_wl_factor_set_integ_method (dcwlf, method);

  nc_galaxy_wl_obs_free (obs);
  g_strfreev (columns_strv);
  g_strv_builder_unref (builder);
  g_list_free_full (columns, g_free);
  ncm_rng_free (rng);
  nc_galaxy_position_factor_free (posf);
  nc_galaxy_redshift_factor_free (zf);
  nc_galaxy_shape_factor_free (gsf);

  return dcwlf;
}

/* The -2lnP of each galaxy */
static NcmVector *
_test_m2lnp (NcDataClusterWLFactor *dcwlf, NcmMSet *mset)
{
  NcmVector *m2lnP = ncm_vector_new (TEST_N_GAL);

  ncm_data_prepare (NCM_DATA (dcwlf), mset);
  nc_data_cluster_wl_factor_eval_m2lnP_gal (dcwlf, mset, m2lnP);

  return m2lnP;
}

void
test_nc_data_cluster_wl_factor_cache_new (TestNcDataClusterWLFactorCache *test, gconstpointer pdata)
{
  _test_models_init (&test->models);
}

void
test_nc_data_cluster_wl_factor_cache_free (TestNcDataClusterWLFactorCache *test, gconstpointer pdata)
{
  _test_models_clear (&test->models);
}

/*
 * A dataset reused over test_sequence, its models changed in place (only the parameters
 * that change at each step), gives the -2lnP of a dataset built with new models at every
 * step, under the integration method of @pdata.
 */
void
test_nc_data_cluster_wl_factor_cache_revisits (TestNcDataClusterWLFactorCache *test, gconstpointer pdata)
{
  const NcDataClusterWLIntegMethod method = *(const NcDataClusterWLIntegMethod *) pdata;
  NcDataClusterWLFactor *cached           = NULL;
  guint k, i;

  for (k = 0; k < G_N_ELEMENTS (test_sequence); k++)
  {
    TestNcDataClusterWLFactorCacheModels fresh_models;
    NcDataClusterWLFactor *fresh;
    NcmVector *v_cached, *v_fresh;

    _test_models_set (&test->models, &test_sequence[k], (k == 0) ? NULL : &test_sequence[k - 1]);

    if (cached == NULL)
      cached = _test_dataset_new (test->models.mset, method);

    _test_models_init (&fresh_models);
    _test_models_set (&fresh_models, &test_sequence[k], NULL);
    fresh = _test_dataset_new (fresh_models.mset, method);

    v_cached = _test_m2lnp (cached, test->models.mset);
    v_fresh  = _test_m2lnp (fresh, fresh_models.mset);

    for (i = 0; i < TEST_N_GAL; i++)
      ncm_assert_cmpdouble_e (ncm_vector_get (v_cached, i), ==, ncm_vector_get (v_fresh, i), 1.0e-12, 0.0);

    ncm_vector_free (v_cached);
    ncm_vector_free (v_fresh);
    nc_data_cluster_wl_factor_free (fresh);
    _test_models_clear (&fresh_models);
  }

  nc_data_cluster_wl_factor_free (cached);
}

/*
 * Each change value changes when its models differ from the previous prepare, including on
 * a return to an earlier state. Only this direction is checked: a value may also change
 * without need, which test_nc_data_cluster_wl_factor_cache_revisits() would not detect,
 * but it never misses a change.
 */
void
test_nc_data_cluster_wl_factor_cache_change_values (TestNcDataClusterWLFactorCache *test, gconstpointer pdata)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (TEST_ELLIP_CONV));
  NcWLClusterLens *lens    = nc_wl_cluster_lens_new (test->models.smd->dist);
  NcGalaxyHaloWL *hwl      = nc_galaxy_halo_wl_new (lens);
  guint64 prev[4]          = {0, 0, 0, 0};
  guint k;

  for (k = 0; k < G_N_ELEMENTS (test_sequence); k++)
  {
    const TestNcDataClusterWLFactorCacheConfig *c = &test_sequence[k];
    const TestNcDataClusterWLFactorCacheConfig *p = (k == 0) ? NULL : &test_sequence[k - 1];
    guint64 cur[4];

    _test_models_set (&test->models, c, p);
    nc_galaxy_shape_factor_prepare (gsf, test->models.mset);
    nc_galaxy_halo_wl_prepare (hwl, test->models.mset);

    cur[0] = nc_galaxy_halo_wl_get_radius_hash (hwl);
    cur[1] = nc_galaxy_halo_wl_get_profile_hash (hwl);
    cur[2] = nc_wl_cluster_lens_get_geometry_hash (lens);
    cur[3] = nc_galaxy_shape_factor_get_pop_hash (gsf);

    if (p != NULL)
    {
      if (c->z_cl != p->z_cl)
      {
        g_assert_cmpuint (cur[0], !=, prev[0]);
        g_assert_cmpuint (cur[2], !=, prev[2]);
      }

      if ((c->log10_mdelta != p->log10_mdelta) || (c->c_delta != p->c_delta) || (c->z_cl != p->z_cl))
        g_assert_cmpuint (cur[1], !=, prev[1]);

      if (c->pop_sigma != p->pop_sigma)
        g_assert_cmpuint (cur[3], !=, prev[3]);
    }

    memcpy (prev, cur, sizeof (cur));
  }

  nc_galaxy_halo_wl_free (hwl);
  nc_wl_cluster_lens_free (lens);
  nc_galaxy_shape_factor_free (gsf);
}

