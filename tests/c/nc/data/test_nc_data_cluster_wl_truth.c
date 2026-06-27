/***************************************************************************
 *            test_nc_data_cluster_wl_truth.c
 *
 *  Truth-table validation of the NcDataClusterWL redshift-integration methods.
 *
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
 * Deterministic "truth table" for the redshift-integration methods. For each
 * (shape, redshift, ellip-conv, coord) configuration the test regenerates - from
 * a fixed seed, so no sample data is stored in the repository - a curated set of
 * galaxies that spans the kinds of integrand that occur in practice: sources in
 * front of / straddling / behind the reduced-shear kink at the lens redshift
 * z_cl, plus tail galaxies. The reference per-galaxy -2lnP is recomputed each run
 * with FIXED_NODES at 320x7, which self-converges to ~1e-12; every cheaper method
 * is validated against it at a tolerance matched to that method's accuracy floor.
 *
 * The pz integrand is a cubic spline (only C2 at its knots) and the adaptive
 * methods find it hardest exactly when the photo-z bump straddles z_cl. Random
 * photo-z draws (z_avg > z_cl) rarely produce that, so the first few pz galaxies
 * are deliberately given spline centres around z_cl (STRADDLE_ZAVG). gauss is
 * analytic and spec is a delta, so neither has such a hard case.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TRUTH_NROWS (80)
#define TRUTH_Z_CL (0.5)
#define TRUTH_GOLDEN_NODES (320)
#define TRUTH_GOLDEN_RULE (7)
#define TRUTH_PROD_NODES (20)
#define TRUTH_PROD_RULE (5)
#define TRUTH_PREC (1.0e-6)

/* pz spline centres placed around z_cl so the photo-z bump straddles the kink. */
static const gdouble STRADDLE_ZAVG[] = {
  0.30, 0.35, 0.40, 0.44, 0.47, 0.50, 0.53, 0.56, 0.60, 0.65, 0.70, 0.80
};

typedef struct _TruthCfg
{
  const gchar *shape;
  const gchar *redshift;
  const gchar *conv;
  const gchar *coord;
  guint seed;
} TruthCfg;

typedef struct _TestNcDataClusterWLTruth
{
  NcmMSet *mset;
  NcDataClusterWL *dcwl;
  NcmVector *golden;
  const gchar *redshift;
} TestNcDataClusterWLTruth;

static void test_nc_data_cluster_wl_truth_new (TestNcDataClusterWLTruth *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_truth_free (TestNcDataClusterWLTruth *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_truth_methods (TestNcDataClusterWLTruth *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  static const gchar *shapes[]    = {"gauss_global", "gauss"};
  static const gchar *redshifts[] = {"spec", "gauss", "pz"};
  static const gchar *convs[]     = {"trace", "trace_det"};
  static const gchar *coords[]    = {"celestial", "cartesian"};
  static TruthCfg cfgs[2 * 3 * 2 * 2];
  guint n = 0;
  guint si, ri, ci, di;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  g_test_set_nonfatal_assertions ();

  for (si = 0; si < G_N_ELEMENTS (shapes); si++)
  {
    for (ri = 0; ri < G_N_ELEMENTS (redshifts); ri++)
    {
      for (ci = 0; ci < G_N_ELEMENTS (convs); ci++)
      {
        for (di = 0; di < G_N_ELEMENTS (coords); di++)
        {
          gchar *path = g_strdup_printf ("/nc/data_cluster_wl_truth/%s/%s/%s/%s",
                                         shapes[si], redshifts[ri], convs[ci], coords[di]);

          cfgs[n].shape    = shapes[si];
          cfgs[n].redshift = redshifts[ri];
          cfgs[n].conv     = convs[ci];
          cfgs[n].coord    = coords[di];
          cfgs[n].seed     = 1000 + n;

          g_test_add (path, TestNcDataClusterWLTruth, &cfgs[n],
                      &test_nc_data_cluster_wl_truth_new,
                      &test_nc_data_cluster_wl_truth_methods,
                      &test_nc_data_cluster_wl_truth_free);
          g_free (path);
          n++;
        }
      }
    }
  }

  g_test_run ();

  return 0;
}

static NcmVector *
_truth_eval (TestNcDataClusterWLTruth *test, NcDataClusterWLIntegMethod method, guint n_nodes, guint rule_n)
{
  NcmVector *v = ncm_vector_new (TRUTH_NROWS);

  g_object_set (test->dcwl, "n-nodes", n_nodes, "rule-n", rule_n, NULL);
  nc_data_cluster_wl_set_integ_method (test->dcwl, method);
  nc_data_cluster_wl_eval_m2lnP_gal (test->dcwl, test->mset, v);

  return v;
}

static void
_truth_gen_redshift (const gchar *redshift, NcGalaxySDObsRedshift *zd, NcmMSet *mset,
                     NcGalaxySDObsRedshiftData *z_data, NcmRNG *rng, gint straddle_i)
{
  if (g_strcmp0 (redshift, "gauss") == 0)
  {
    nc_galaxy_sd_obs_redshift_gauss_gen (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (zd), mset, z_data, 0.03, rng);
  }
  else if (g_strcmp0 (redshift, "spec") == 0)
  {
    nc_galaxy_sd_obs_redshift_spec_gen (NC_GALAXY_SD_OBS_REDSHIFT_SPEC (zd), mset, z_data, rng);
  }
  else /* pz */
  {
    const guint npoints = 20;
    const gdouble z_min = 0.01;
    const gdouble z_max = 5.0;
    const gdouble z_avg = (straddle_i >= 0) ? STRADDLE_ZAVG[straddle_i] : ncm_rng_uniform_gen (rng, z_min + 0.5, z_max - 0.5);
    const gdouble z_sd  = 0.03 * (1.0 + z_avg);
    const gdouble x_min = MAX (z_avg - 5.0 * z_sd, z_min);
    const gdouble x_max = MIN (z_avg + 5.0 * z_sd, z_max);
    NcmVector *xv       = ncm_vector_new (npoints);
    NcmVector *yv       = ncm_vector_new (npoints);
    NcmSpline *pz;
    guint j;

    for (j = 0; j < npoints; j++)
    {
      const gdouble x = x_min + (x_max - x_min) * j / (npoints - 1.0);
      const gdouble y = exp (-0.5 * gsl_pow_2 ((x - z_avg) / z_sd)) / (sqrt (2.0 * M_PI) * z_sd);

      ncm_vector_set (xv, j, x);
      ncm_vector_set (yv, j, y);
    }

    pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

    nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (zd), z_data, pz);
    nc_galaxy_sd_obs_redshift_pz_prepare (NC_GALAXY_SD_OBS_REDSHIFT_PZ (zd), z_data);
    nc_galaxy_sd_obs_redshift_pz_gen (NC_GALAXY_SD_OBS_REDSHIFT_PZ (zd), mset, z_data, rng);

    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_spline_free (pz);
  }
}

static void
_truth_gen_shape (const gchar *shape, NcGalaxySDShape *sd, NcmMSet *mset,
                  NcGalaxySDShapeData *s_data, NcGalaxyWLObsCoord coord, NcmRNG *rng)
{
  const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
  const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
  const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
  const gdouble std_noise = ncm_rng_uniform_gen (rng, 0.001, 0.003);

  if (g_strcmp0 (shape, "gauss_global") == 0)
  {
    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (sd), mset, s_data, std_noise, c1, c2, m, coord, rng);
  }
  else
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.05, 0.15);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (sd), mset, s_data, std_shape, std_noise, c1, c2, m, coord, rng);
  }
}

/* Regenerate the curated sample for one configuration (deterministic in @seed)
 * and recompute the FIXED-320x7 golden reference. Leaves dcwl at the production
 * fixed-node configuration. */
static void
test_nc_data_cluster_wl_truth_new (TestNcDataClusterWLTruth *test, gconstpointer pdata)
{
  const TruthCfg *cfg             = (const TruthCfg *) pdata;
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, cfg->seed);
  NcHICosmo *cosmo                = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                = nc_distance_new (100.0);
  NcHaloMassSummary *hms          = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp        = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp              = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd     = nc_wl_surface_mass_density_new (dist);
  NcGalaxySDTrueRedshift *z_true  = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  const NcGalaxyWLObsEllipConv ec = (g_strcmp0 (cfg->conv, "trace") == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  const NcGalaxyWLObsCoord co     = (g_strcmp0 (cfg->coord, "celestial") == 0) ? NC_GALAXY_WL_OBS_COORD_CELESTIAL : NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;
  const gdouble min_radius        = ncm_rng_uniform_gen (rng, 0.1, 0.5);
  const gdouble max_radius        = ncm_rng_uniform_gen (rng, 2.0, 5.0);
  const gdouble ra                = ncm_rng_uniform_gen (rng, -180, 180);
  const gdouble dec               = ncm_rng_uniform_gen (rng, -90, 90);
  const gint n_straddle           = (g_strcmp0 (cfg->redshift, "pz") == 0) ? (gint) G_N_ELEMENTS (STRADDLE_ZAVG) : 0;
  NcGalaxySDPosition *p_dist      = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (ra - 0.2, ra + 0.2, dec - 0.2, dec + 0.2));
  NcGalaxySDShape *s_dist;
  NcGalaxySDObsRedshift *z_dist;
  NcGalaxySDObsRedshiftData *z_data;
  NcGalaxySDPositionData *p_data;
  NcGalaxySDShapeData *s_data;
  NcGalaxyWLObs *obs;
  GStrvBuilder *builder = g_strv_builder_new ();
  GList *columns, *l;
  GStrv columns_strv;
  guint i;

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, ncm_rng_uniform_gen (rng, 13.5, 15.5));
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, ra);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, dec);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_Z, TRUTH_Z_CL);
  nc_halo_position_prepare (hp, cosmo);

  if (g_strcmp0 (cfg->shape, "gauss_global") == 0)
  {
    s_dist = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_hsm_gauss_global_new (ec));
    ncm_model_param_set (NCM_MODEL (s_dist), NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SIGMA, ncm_rng_uniform_gen (rng, 0.05, 0.15));
  }
  else
  {
    s_dist = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_hsm_gauss_new (ec));
  }

  if (g_strcmp0 (cfg->redshift, "spec") == 0)
    z_dist = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true, 0.0, 5.0));
  else if (g_strcmp0 (cfg->redshift, "gauss") == 0)
    z_dist = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_gauss_new (z_true, 0.1, 4.8));
  else
    z_dist = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_pz_new ());

  test->mset     = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);
  test->redshift = cfg->redshift;

  z_data = nc_galaxy_sd_obs_redshift_data_new (z_dist);
  p_data = nc_galaxy_sd_position_data_new (p_dist, z_data);
  s_data = nc_galaxy_sd_shape_data_new (s_dist, p_data);

  columns = nc_galaxy_sd_shape_data_required_columns (s_data);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  obs          = nc_galaxy_wl_obs_new (ec, co, TRUTH_NROWS, columns_strv);

  for (i = 0; i < TRUTH_NROWS; i++)
  {
    const gint straddle_i = ((gint) i < n_straddle) ? (gint) i : -1;
    gdouble radius;

    _truth_gen_redshift (cfg->redshift, z_dist, test->mset, z_data, rng, straddle_i);

    do {
      nc_galaxy_sd_position_gen (p_dist, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    _truth_gen_shape (cfg->shape, s_dist, test->mset, s_data, co, rng);
    nc_galaxy_sd_shape_data_write_row (s_data, obs, i);
  }

  test->dcwl = nc_data_cluster_wl_new ();
  nc_data_cluster_wl_set_cut (test->dcwl, min_radius, max_radius);
  nc_data_cluster_wl_set_obs (test->dcwl, obs);
  g_object_set (test->dcwl, "enable-parallel", FALSE, NULL);
  g_object_set (test->dcwl, "auto-nodes", FALSE, NULL);
  g_object_set (test->dcwl, "node_reltol", 1e-6, NULL);

  test->golden = _truth_eval (test, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES, TRUTH_GOLDEN_NODES, TRUTH_GOLDEN_RULE);

  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_mass_summary_free (hms);
  nc_halo_density_profile_free (dp);
  nc_halo_position_free (hp);
  nc_wl_surface_mass_density_free (smd);
  nc_galaxy_sd_true_redshift_free (z_true);
  nc_galaxy_sd_position_free (p_dist);
  nc_galaxy_sd_shape_free (s_dist);
  nc_galaxy_sd_obs_redshift_free (z_dist);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_wl_obs_free (obs);
  g_strfreev (columns_strv);
  g_strv_builder_unref (builder);
  g_list_free_full (columns, g_free);
  ncm_rng_free (rng);
}

static void
test_nc_data_cluster_wl_truth_free (TestNcDataClusterWLTruth *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  ncm_vector_clear (&test->golden);
  NCM_TEST_FREE (nc_data_cluster_wl_free, test->dcwl);
}

/* Per-galaxy comparison of @v against the golden -2lnP vector: fail (reporting
 * the closest-to-tolerance galaxy) if any galaxy exceeds reltol*max(|.|) + abstol.
 * Reports max|abs| separately, since the worst-excess galaxy (abstol-bound, near
 * zero -2lnP) is usually not the largest-error one. */
static void
_truth_cmp (const gchar *label, NcmVector *v, NcmVector *golden, gdouble reltol, gdouble abstol)
{
  const guint len      = ncm_vector_len (golden);
  gdouble worst_excess = -G_MAXDOUBLE;
  gdouble max_abs      = 0.0;
  guint worst_i        = 0;
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble ai     = ncm_vector_get (v, i);
    const gdouble bi     = ncm_vector_get (golden, i);
    const gdouble mean   = GSL_MAX (fabs (ai), fabs (bi));
    const gdouble adiff  = fabs (ai - bi);
    const gdouble excess = adiff - (reltol * mean + abstol);

    g_assert_true (gsl_finite (ai));

    max_abs = GSL_MAX (max_abs, adiff);

    if (excess > worst_excess)
    {
      worst_excess = excess;
      worst_i      = i;
    }
  }

  {
    const gdouble ai  = ncm_vector_get (v, worst_i);
    const gdouble bi  = ncm_vector_get (golden, worst_i);
    const gdouble rel = (fabs (bi) > 0.0) ? fabs (ai - bi) / fabs (bi) : 0.0;

    g_test_message ("%s: max|abs|=%.3e; closest-to-tol at %u: % .17g vs golden % .17g (abs %.3e, rel %.3e; reltol %.3e abstol %.3e)",
                    label, max_abs, worst_i, ai, bi, fabs (ai - bi), rel, reltol, abstol);

    if (worst_excess > 0.0)
      g_error ("%s disagrees with golden at galaxy %u: % .17g vs % .17g (abs %.3e, rel %.3e exceeds reltol %.3e abstol %.3e)",
               label, worst_i, ai, bi, fabs (ai - bi), rel, reltol, abstol);
  }
}

static void
test_nc_data_cluster_wl_truth_methods (TestNcDataClusterWLTruth *test, gconstpointer pdata)
{
  const gboolean is_spec = (g_strcmp0 (test->redshift, "spec") == 0);
  const gboolean is_pz   = (g_strcmp0 (test->redshift, "pz") == 0);

  /* Per-redshift tolerance vs the FIXED-320x7 golden over this frozen sample. The
   * production FIXED 20-node truncation is the largest error in each case and
   * sets the bar (adaptive at prec=1e-6 is well inside it): spec is a delta (all
   * methods exact); the analytic gauss kernel reaches ~1.6e-8; the pz cubic
   * spline (only C2 at its knots) reaches ~2e-6. The bars below keep a few x
   * margin - tight enough to catch dropping the z_cl split or the effective
   * support, or a node-count regression. */
  const gdouble abstol = is_spec ? 1.0e-9 : (is_pz ? 1.0e-5 : 1.0e-6);
  NcmVector *vF        = _truth_eval (test, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES, TRUTH_PROD_NODES, TRUTH_PROD_RULE);

  g_object_set (test->dcwl, "auto-nodes", TRUE, NULL);

  NcmVector *vA = _truth_eval (test, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES, TRUTH_PROD_NODES, TRUTH_PROD_RULE);
  NcmVector *vL = _truth_eval (test, NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT, TRUTH_PROD_NODES, TRUTH_PROD_RULE);
  NcmVector *vC = _truth_eval (test, NC_DATA_CLUSTER_WL_INTEG_METHOD_CUBATURE, TRUTH_PROD_NODES, TRUTH_PROD_RULE);

  _truth_cmp ("FIXED vs golden", vF, test->golden, 1.0e-6, abstol);
  _truth_cmp ("FIXED AUTO vs golden", vA, test->golden, 1.0e-6, abstol);
  _truth_cmp ("LNINT vs golden", vL, test->golden, 1.0e-6, abstol);
  _truth_cmp ("CUBATURE vs golden", vC, test->golden, 1.0e-6, abstol);

  ncm_vector_free (vF);
  ncm_vector_free (vL);
  ncm_vector_free (vC);
}

