/***************************************************************************
 *            test_nc_data_cluster_wl_factor_truth.c
 *
 *  Thu Jul 16 2026
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
 * Deterministic "truth table" for NcDataClusterWLFactor, self-contained: no
 * legacy oracle is involved. For each (shape scheme, population, redshift
 * scheme, ellip-conv, coord) configuration the test regenerates -- from a
 * fixed seed, so no sample data is stored in the repository -- a curated set
 * of galaxies spanning the hard integrand cases (sources in front of /
 * straddling / behind the reduced-shear kink at the lens redshift z_cl),
 * recomputes a self-converged FIXED_NODES-320x7 golden reference (already
 * shown to have full parity with LNINT/CUBATURE for the same integrand
 * family -- see the nc-data-cluster-wl-factor-integ-parity work), and
 * validates the cheaper production configuration (FIXED_NODES-20x5, LNINT,
 * CUBATURE) against it.
 *
 * The Factor pipeline has several independently-coded *exact*
 * intrinsic-ellipticity marginalization schemes (SeriesLensed, Quad,
 * Laplace, FixedQuad) and a non-Gaussian population (Beta); this file gives
 * each of them its own self-contained correctness net. The matrix below is
 * deliberately curated, not a full cross-product: each shape scheme is
 * exercised against the golden reference at least once, VarAdd+Gauss is
 * kept as the cheap "does the wiring still work" baseline, and
 * SeriesLensed/Beta and Quad/Beta each get one spline-straddle (hardest
 * integrand) case. A final, non-golden test cross-checks SeriesLensed
 * against Quad directly on the identical generated galaxy set: two
 * independently-coded exact schemes for the same physical marginal integral
 * must agree with each other.
 *
 * The Spline redshift scheme's pz is a cubic spline (only C2 at its knots),
 * the hardest case for the adaptive methods, exactly when the photo-z bump
 * straddles z_cl. The Composed scheme's joint density (population x
 * Gaussian photo-z kernel) is smooth, so its tolerance is tighter.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TRUTH_NROWS (14)
#define TRUTH_Z_CL (0.5)
#define TRUTH_GOLDEN_NODES (320)
#define TRUTH_GOLDEN_RULE (7)
#define TRUTH_PROD_NODES (20)
#define TRUTH_PROD_RULE (5)
#define TRUTH_PREC (1.0e-6)
#define TRUTH_SERIES_ORDER (4)

/* NcGalaxyShapeFactorQuad's *default* reltol (1e-7) is appropriate for
 * production single-galaxy work but pathologically expensive to call
 * O(golden-nodes x nrows) times in a truth table: measured at ~70ms per
 * eval_marginal() call (vs ~1us for SeriesLensed/VarAdd on the same input),
 * driven entirely by the 2D disc cubature's own adaptive refinement, not by
 * anything z-integration-related. Loosening to 1e-5 (measured ~1.6ms/call,
 * a ~45x speedup) changes eval_marginal()'s own output at the ~1e-5 relative
 * level -- negligible next to this file's ~1e-6/1e-5 z-integration
 * tolerances -- while keeping the Quad-specific configs' total runtime
 * reasonable. Quad configs also use fewer rows and a smaller golden node
 * count than the other shape schemes (see TRUTH_QUAD_NROWS below): the
 * z-integrand's smoothness is governed by the redshift scheme, not the
 * shape scheme, and that has already been validated at full resolution by
 * the (near-free) VarAdd/SeriesLensed configs above. */
#define TRUTH_QUAD_RELTOL (1.0e-5)
#define TRUTH_QUAD_NROWS (5)
#define TRUTH_QUAD_GOLDEN_NODES (60)
#define TRUTH_QUAD_GOLDEN_RULE (5)

/* pz spline centres placed around z_cl so the photo-z bump straddles the kink
 * (see test_nc_data_cluster_wl_truth.c's identical STRADDLE_ZAVG). */
static const gdouble STRADDLE_ZAVG[] = {
  0.32, 0.40, 0.46, 0.50, 0.54, 0.60, 0.68, 0.80
};

typedef struct _TruthCfg
{
  const gchar *shape;    /* "var_add" | "series_lensed" | "quad" */
  const gchar *pop;      /* "gauss" | "beta" */
  const gchar *redshift; /* "composed" | "spline" */
  const gchar *conv;
  const gchar *coord;
  guint seed;
  guint nrows;
  guint golden_nodes;
  guint golden_rule;
} TruthCfg;

typedef struct _TestNcDataClusterWLFactorTruth
{
  NcmMSet *mset;
  NcDataClusterWLFactor *dcwlf;
  NcmVector *golden;
  const gchar *redshift;
  const gchar *shape;
  guint nrows;
} TestNcDataClusterWLFactorTruth;

static void test_nc_data_cluster_wl_factor_truth_new (TestNcDataClusterWLFactorTruth *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_factor_truth_free (TestNcDataClusterWLFactorTruth *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_factor_truth_methods (TestNcDataClusterWLFactorTruth *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_factor_truth_cross_shape (void);

gint
main (gint argc, gchar *argv[])
{
  static const TruthCfg cfgs[] = {
    /* Cheap baseline: mirrors legacy's own gauss_global+gauss combo, checks
     * the orchestrator wiring (prepare/cache cascade, obs plumbing) is
     * unchanged for the closed-form approximation. */
    { "var_add",       "gauss", "composed", "trace_det", "celestial", 2001, TRUTH_NROWS, TRUTH_GOLDEN_NODES, TRUTH_GOLDEN_RULE },
    { "var_add",       "gauss", "spline",   "trace",     "cartesian", 2002, TRUTH_NROWS, TRUTH_GOLDEN_NODES, TRUTH_GOLDEN_RULE },

    /* SeriesLensed: exact truncated-series marginalization, both populations
     * it supports, one of each redshift scheme (Beta paired with the harder
     * Spline straddle case). Its own eval_marginal() is ~us-scale (a closed-
     * form Taylor-series evaluation), so it affords the same full resolution
     * as VarAdd. */
    { "series_lensed", "gauss", "composed", "trace_det", "celestial", 2003, TRUTH_NROWS, TRUTH_GOLDEN_NODES, TRUTH_GOLDEN_RULE },
    { "series_lensed", "beta",  "spline",   "trace_det", "cartesian", 2004, TRUTH_NROWS, TRUTH_GOLDEN_NODES, TRUTH_GOLDEN_RULE },

    /* Quad: exact 2D disc quadrature, same shape-scheme-coverage rationale,
     * but its own eval_marginal() is orders of magnitude more expensive (see
     * TRUTH_QUAD_RELTOL's comment) -- fewer rows, a smaller golden node
     * count, and a loosened internal reltol keep this file's total runtime
     * reasonable. */
    { "quad",          "gauss", "spline",   "trace",     "celestial", 2005, TRUTH_QUAD_NROWS, TRUTH_QUAD_GOLDEN_NODES, TRUTH_QUAD_GOLDEN_RULE },
    { "quad",          "beta",  "composed", "trace_det", "cartesian", 2006, TRUTH_QUAD_NROWS, TRUTH_QUAD_GOLDEN_NODES, TRUTH_QUAD_GOLDEN_RULE },
  };
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (cfgs); i++)
  {
    gchar *path = g_strdup_printf ("/nc/data_cluster_wl_factor_truth/%s/%s/%s/%s/%s",
                                   cfgs[i].shape, cfgs[i].pop, cfgs[i].redshift, cfgs[i].conv, cfgs[i].coord);

    g_test_add (path, TestNcDataClusterWLFactorTruth, &cfgs[i],
                &test_nc_data_cluster_wl_factor_truth_new,
                &test_nc_data_cluster_wl_factor_truth_methods,
                &test_nc_data_cluster_wl_factor_truth_free);
    g_free (path);
  }

  g_test_add_func ("/nc/data_cluster_wl_factor_truth/cross_shape/series_lensed_vs_quad",
                   &test_nc_data_cluster_wl_factor_truth_cross_shape);

  g_test_run ();

  return 0;
}

static NcmVector *
_truth_eval (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, guint nrows, NcDataClusterWLIntegMethod method, guint n_nodes, guint rule_n)
{
  NcmVector *v = ncm_vector_new (nrows);

  nc_data_cluster_wl_factor_set_n_nodes (dcwlf, n_nodes);
  nc_data_cluster_wl_factor_set_rule_n (dcwlf, rule_n);
  nc_data_cluster_wl_factor_set_integ_method (dcwlf, method);
  nc_data_cluster_wl_factor_eval_m2lnP_gal (dcwlf, mset, v);

  return v;
}

/* Builds a cubic-notaknot pz spline of a Gaussian bump: STRADDLE_ZAVG[straddle_i]
 * when @straddle_i >= 0 (deliberately centred near z_cl), otherwise a random
 * centre away from the domain edges. Mirrors
 * test_nc_data_cluster_wl_truth.c's pz-branch construction exactly, just
 * standalone (the Spline scheme takes the spline directly, no SD wrapper). */
static NcmSpline *
_truth_make_pz_spline (NcmRNG *rng, gint straddle_i)
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

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return pz;
}

static void
_truth_gen_redshift (const gchar *redshift, NcGalaxyRedshiftFactor *z_dist, NcmMSet *mset,
                     NcGalaxyRedshiftFactorData *z_data, NcmRNG *rng, gint straddle_i,
                     NcGalaxyWLObs *obs, guint i)
{
  if (g_strcmp0 (redshift, "spline") == 0)
  {
    NcmSpline *pz = _truth_make_pz_spline (rng, straddle_i);

    nc_galaxy_redshift_factor_spline_data_set (NC_GALAXY_REDSHIFT_FACTOR_SPLINE (z_dist), z_data, pz);
    nc_galaxy_redshift_factor_gen (z_dist, mset, z_data, rng);
    ncm_spline_free (pz);
  }
  else /* composed */
  {
    /* Composed's gen() only draws (z, zp) FROM a given photo-z scatter
     * sigma0 -- sigma0 itself is per-galaxy survey-depth data, set here
     * exactly like a real pipeline would (not something gen() invents). It
     * must be written onto @obs and loaded back via read_row() into
     * z_data's own ldata fragment before gen() can use it: there is no
     * public accessor into Composed's opaque per-galaxy fragment otherwise.
     * Without this, sigma0 stays at its zero default, gen() draws zp
     * exactly at z (sigmaz=0), and the resulting near-delta Gaussian kernel
     * makes the orchestrator's own z-integral fail to converge (GSL qag
     * error) -- a real trap for any caller driving Composed generatively. */
    const gdouble sigma0 = ncm_rng_uniform_gen (rng, 0.02, 0.06);

    nc_galaxy_wl_obs_set (obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, sigma0, NULL);
    nc_galaxy_redshift_factor_data_read_row (z_data, obs, i);
    nc_galaxy_redshift_factor_gen (z_dist, mset, z_data, rng);
  }
}

static void
_truth_gen_shape (NcGalaxyShapeFactor *s_dist, NcmMSet *mset, NcGalaxyShapeFactorData *s_data,
                  NcWLEllipticityFrame coord, NcmRNG *rng)
{
  const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
  const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
  const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
  const gdouble std_noise = ncm_rng_uniform_gen (rng, 0.005, 0.02);

  nc_galaxy_shape_factor_data_set (s_dist, s_data, 0.0, 0.0, std_noise, c1, c2, m, coord);
  nc_galaxy_shape_factor_gen (s_dist, mset, s_data, rng);
}

/* Regenerate the curated sample for one configuration (deterministic in @seed)
 * and recompute the FIXED-320x7 golden reference. Leaves dcwlf at the
 * production fixed-node configuration. */
static void
test_nc_data_cluster_wl_factor_truth_new (TestNcDataClusterWLFactorTruth *test, gconstpointer pdata)
{
  const TruthCfg *cfg                     = (const TruthCfg *) pdata;
  NcmRNG *rng                             = ncm_rng_seeded_new (NULL, cfg->seed);
  NcHICosmo *cosmo                        = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                        = nc_distance_new (100.0);
  NcHaloMassSummary *hms                  = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp                = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                      = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd             = nc_wl_surface_mass_density_new (dist);
  NcGalaxyRedshiftPopLSSTSRD *pop_z       = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z         = nc_galaxy_redshift_obs_gauss_new ();
  const NcGalaxyWLObsEllipConv ec         = (g_strcmp0 (cfg->conv, "trace") == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  const NcWLEllipticityFrame co           = (g_strcmp0 (cfg->coord, "celestial") == 0) ? NC_WL_ELLIPTICITY_FRAME_CELESTIAL : NC_WL_ELLIPTICITY_FRAME_CARTESIAN;
  const gdouble min_radius                = ncm_rng_uniform_gen (rng, 0.1, 0.5);
  const gdouble max_radius                = ncm_rng_uniform_gen (rng, 2.0, 5.0);
  const gdouble ra                        = ncm_rng_uniform_gen (rng, -180, 180);
  const gdouble dec                       = ncm_rng_uniform_gen (rng, -90, 90);
  const gint n_straddle                   = (g_strcmp0 (cfg->redshift, "spline") == 0) ? (gint) G_N_ELEMENTS (STRADDLE_ZAVG) : 0;
  NcGalaxyPositionFactor *position_factor = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (ra - 0.2, ra + 0.2, dec - 0.2, dec + 0.2));
  NcGalaxyRedshiftFactor *redshift_factor;
  NcGalaxyShapeFactor *shape_factor;
  NcGalaxyShapePop *pop_shape;
  NcGalaxyPositionFactorData *p_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *s_data;
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

  if (g_strcmp0 (cfg->pop, "gauss") == 0)
  {
    pop_shape = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
    ncm_model_param_set (NCM_MODEL (pop_shape), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, ncm_rng_uniform_gen (rng, 0.05, 0.2));
  }
  else /* beta: alpha=0.9, beta=4.1 (the class's own pre->=1-bound default,
        * still directly settable here) matches the population used
        * elsewhere to exercise the non-integer-exponent path (see
        * test_nc_galaxy_shape_pop_series.c). */
  {
    pop_shape = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_beta_new ());
    ncm_model_param_set_by_name (NCM_MODEL (pop_shape), "alpha", 0.9, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (pop_shape), "beta", 4.1, NULL);
  }

  if (g_strcmp0 (cfg->redshift, "composed") == 0)
    redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0));
  else
    redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_spline_new ());

  if (g_strcmp0 (cfg->shape, "var_add") == 0)
  {
    shape_factor = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (ec));
  }
  else if (g_strcmp0 (cfg->shape, "series_lensed") == 0)
  {
    shape_factor = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_series_lensed_new (ec, TRUTH_SERIES_ORDER));
  }
  else /* quad */
  {
    NcGalaxyShapeFactorQuad *gsfq = nc_galaxy_shape_factor_quad_new (ec);

    nc_galaxy_shape_factor_quad_set_reltol (gsfq, TRUTH_QUAD_RELTOL);
    shape_factor = NC_GALAXY_SHAPE_FACTOR (gsfq);
  }

  test->mset = ncm_mset_empty_new ();
  ncm_mset_set (test->mset, NCM_MODEL (cosmo), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (dp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (hp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (smd), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (pop_z), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (obs_z), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (pop_shape), NULL);

  test->redshift = cfg->redshift;
  test->shape    = cfg->shape;
  test->nrows    = cfg->nrows;

  p_data = nc_galaxy_position_factor_data_new (position_factor, test->mset);
  z_data = nc_galaxy_redshift_factor_data_new (redshift_factor, test->mset);
  s_data = nc_galaxy_shape_factor_data_new (shape_factor, test->mset, p_data, z_data);

  columns = nc_galaxy_shape_factor_data_required_columns (s_data);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  obs          = nc_galaxy_wl_obs_new (ec, co, cfg->nrows, columns_strv);

  for (i = 0; i < cfg->nrows; i++)
  {
    const gint straddle_i = ((gint) i < n_straddle) ? (gint) i : -1;
    gdouble radius;

    _truth_gen_redshift (cfg->redshift, redshift_factor, test->mset, z_data, rng, straddle_i, obs, i);

    do {
      nc_galaxy_position_factor_gen (position_factor, test->mset, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    _truth_gen_shape (shape_factor, test->mset, s_data, co, rng);
    nc_galaxy_shape_factor_data_write_row (s_data, obs, i);
  }

  test->dcwlf = nc_data_cluster_wl_factor_new (position_factor, redshift_factor, shape_factor);
  nc_data_cluster_wl_factor_set_cut (test->dcwlf, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (test->dcwlf, obs);
  nc_data_cluster_wl_factor_set_prec (test->dcwlf, TRUTH_PREC);

  test->golden = _truth_eval (test->dcwlf, test->mset, cfg->nrows, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES, cfg->golden_nodes, cfg->golden_rule);

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
  nc_galaxy_wl_obs_free (obs);
  g_strfreev (columns_strv);
  g_strv_builder_unref (builder);
  g_list_free_full (columns, g_free);
  ncm_rng_free (rng);
}

static void
test_nc_data_cluster_wl_factor_truth_free (TestNcDataClusterWLFactorTruth *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  ncm_vector_clear (&test->golden);
  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, test->dcwlf);
}

/* Per-galaxy comparison of @v against the golden -2lnP vector: fail (reporting
 * the closest-to-tolerance galaxy) if any galaxy exceeds reltol*max(|.|) + abstol.
 * Reports max|abs| separately, since the worst-excess galaxy (abstol-bound, near
 * zero -2lnP) is usually not the largest-error one. Verbatim from
 * test_nc_data_cluster_wl_truth.c's _truth_cmp. */
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
test_nc_data_cluster_wl_factor_truth_methods (TestNcDataClusterWLFactorTruth *test, gconstpointer pdata)
{
  /* Per-redshift-scheme tolerance vs the golden reference over this frozen
   * sample: Composed's joint density (population x Gaussian photo-z kernel)
   * is smooth in z, reaching ~1e-6; Spline's cubic pz (only C2 at its knots,
   * deliberately straddling z_cl for several galaxies) is the harder case,
   * reaching ~1e-5 -- both mirror test_nc_data_cluster_wl_truth.c's
   * gauss/pz bars. The bars below keep a few x margin -- tight enough to
   * catch dropping the z_cl split or the effective support, or a node-count
   * regression, for ANY of the shape schemes exercised across the matrix.
   * Quad's own marginal is only accurate to TRUTH_QUAD_RELTOL (see its
   * comment) and golden/production evaluate it at different z-nodes, so
   * their residual disagreement floor is set by that, not by the
   * z-integration itself -- the abstol is widened accordingly for it. */
  const gboolean is_spline = (g_strcmp0 (test->redshift, "spline") == 0);
  const gboolean is_quad   = (g_strcmp0 (test->shape, "quad") == 0);
  const gdouble abstol     = GSL_MAX (is_spline ? 1.0e-5 : 1.0e-6, is_quad ? 5.0e-4 : 0.0);
  NcmVector *vF            = _truth_eval (test->dcwlf, test->mset, test->nrows, NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES, TRUTH_PROD_NODES, TRUTH_PROD_RULE);
  NcmVector *vL            = _truth_eval (test->dcwlf, test->mset, test->nrows, NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT, TRUTH_PROD_NODES, TRUTH_PROD_RULE);
  NcmVector *vC            = _truth_eval (test->dcwlf, test->mset, test->nrows, NC_DATA_CLUSTER_WL_INTEG_METHOD_CUBATURE, TRUTH_PROD_NODES, TRUTH_PROD_RULE);

  _truth_cmp ("FIXED vs golden", vF, test->golden, 1.0e-6, abstol);
  _truth_cmp ("LNINT vs golden", vL, test->golden, 1.0e-6, abstol);
  _truth_cmp ("CUBATURE vs golden", vC, test->golden, 1.0e-6, abstol);

  ncm_vector_free (vF);
  ncm_vector_free (vL);
  ncm_vector_free (vC);
}

/* Direct shape-scheme-vs-shape-scheme cross-check: SeriesLensed (truncated
 * g-series) and Quad (2D disc cubature) are two independently-coded *exact*
 * marginalization schemes -- unlike the golden-vs-production tests above
 * (which only probe the z-integration's own convergence for a FIXED shape
 * scheme), this probes the shape marginalization itself, by feeding the
 * IDENTICAL generated galaxy set (same ra/dec/z/epsilon_obs/std_noise/c/m,
 * same Gauss population) through both and comparing -2lnP_gal directly. No
 * golden reference is needed: the two schemes' own physical answer must
 * agree with each other. Both are configured with LNINT at a tight
 * precision, so any disagreement beyond adaptive-quadrature-level tolerance
 * reflects a real discrepancy between the two marginalization
 * implementations. */
static void
test_nc_data_cluster_wl_factor_truth_cross_shape (void)
{
  NcmRNG *rng                             = ncm_rng_seeded_new (NULL, 3001);
  NcHICosmo *cosmo                        = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                        = nc_distance_new (100.0);
  NcHaloMassSummary *hms                  = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp                = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                      = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd             = nc_wl_surface_mass_density_new (dist);
  NcGalaxyRedshiftPopLSSTSRD *pop_z       = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z         = nc_galaxy_redshift_obs_gauss_new ();
  NcGalaxyShapePopGauss *pop_shape        = nc_galaxy_shape_pop_gauss_new ();
  const NcGalaxyWLObsEllipConv ec         = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  const NcWLEllipticityFrame co           = NC_WL_ELLIPTICITY_FRAME_CELESTIAL;
  const gdouble min_radius                = 0.2;
  const gdouble max_radius                = 3.0;
  const gdouble ra                        = 10.0;
  const gdouble dec                       = -20.0;
  NcGalaxyPositionFactor *position_factor = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (ra - 0.2, ra + 0.2, dec - 0.2, dec + 0.2));
  NcGalaxyRedshiftFactor *redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0));
  NcGalaxyShapeFactor *shape_factor_gen   = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_series_lensed_new (ec, TRUTH_SERIES_ORDER));
  NcGalaxyShapeFactor *shape_factor_lensed;
  NcGalaxyShapeFactor *shape_factor_quad;
  NcDataClusterWLFactor *dcwlf_lensed;
  NcDataClusterWLFactor *dcwlf_quad;
  NcmMSet *mset;
  NcGalaxyPositionFactorData *p_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *s_data;
  NcGalaxyWLObs *obs;
  GStrvBuilder *builder = g_strv_builder_new ();
  GList *columns, *l;
  GStrv columns_strv;
  NcmVector *v_lensed, *v_quad;
  guint i;

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, 14.3);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, ra);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, dec);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_Z, TRUTH_Z_CL);
  nc_halo_position_prepare (hp, cosmo);
  ncm_model_param_set (NCM_MODEL (pop_shape), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.28);

  mset = ncm_mset_empty_new ();
  ncm_mset_set (mset, NCM_MODEL (cosmo), NULL);
  ncm_mset_set (mset, NCM_MODEL (dp), NULL);
  ncm_mset_set (mset, NCM_MODEL (hp), NULL);
  ncm_mset_set (mset, NCM_MODEL (smd), NULL);
  ncm_mset_set (mset, NCM_MODEL (pop_z), NULL);
  ncm_mset_set (mset, NCM_MODEL (obs_z), NULL);
  ncm_mset_set (mset, NCM_MODEL (pop_shape), NULL);

  p_data = nc_galaxy_position_factor_data_new (position_factor, mset);
  z_data = nc_galaxy_redshift_factor_data_new (redshift_factor, mset);
  s_data = nc_galaxy_shape_factor_data_new (shape_factor_gen, mset, p_data, z_data);

  columns = nc_galaxy_shape_factor_data_required_columns (s_data);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  columns_strv = g_strv_builder_end (builder);
  obs          = nc_galaxy_wl_obs_new (ec, co, TRUTH_QUAD_NROWS, columns_strv);

  for (i = 0; i < TRUTH_QUAD_NROWS; i++)
  {
    gdouble radius;
    const gdouble sigma0 = ncm_rng_uniform_gen (rng, 0.02, 0.06);

    nc_galaxy_wl_obs_set (obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, sigma0, NULL);
    nc_galaxy_redshift_factor_data_read_row (z_data, obs, i);
    nc_galaxy_redshift_factor_gen (redshift_factor, mset, z_data, rng);

    do {
      nc_galaxy_position_factor_gen (position_factor, mset, p_data, rng);
      radius = nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, p_data->ra, p_data->dec);
    } while ((radius < min_radius) || (radius > max_radius));

    _truth_gen_shape (shape_factor_gen, mset, s_data, co, rng);
    nc_galaxy_shape_factor_data_write_row (s_data, obs, i);
  }

  /* Two fresh, independent Factor/orchestrator stacks sharing only @mset and
   * @obs: SeriesLensed and Quad. Quad's reltol is loosened for the same
   * runtime reason as the main matrix's Quad configs (see
   * TRUTH_QUAD_RELTOL's comment). */
  shape_factor_lensed = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_series_lensed_new (ec, TRUTH_SERIES_ORDER));
  shape_factor_quad   = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_quad_new (ec));
  nc_galaxy_shape_factor_quad_set_reltol (NC_GALAXY_SHAPE_FACTOR_QUAD (shape_factor_quad), TRUTH_QUAD_RELTOL);

  dcwlf_lensed = nc_data_cluster_wl_factor_new (
    NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (ra - 0.2, ra + 0.2, dec - 0.2, dec + 0.2)),
    NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0)),
    shape_factor_lensed);
  nc_data_cluster_wl_factor_set_cut (dcwlf_lensed, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (dcwlf_lensed, obs);
  nc_data_cluster_wl_factor_set_prec (dcwlf_lensed, TRUTH_PREC);
  v_lensed = _truth_eval (dcwlf_lensed, mset, TRUTH_QUAD_NROWS, NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT, TRUTH_PROD_NODES, TRUTH_PROD_RULE);

  dcwlf_quad = nc_data_cluster_wl_factor_new (
    NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (ra - 0.2, ra + 0.2, dec - 0.2, dec + 0.2)),
    NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0)),
    shape_factor_quad);
  nc_data_cluster_wl_factor_set_cut (dcwlf_quad, min_radius, max_radius);
  nc_data_cluster_wl_factor_set_obs (dcwlf_quad, obs);
  nc_data_cluster_wl_factor_set_prec (dcwlf_quad, TRUTH_PREC);
  v_quad = _truth_eval (dcwlf_quad, mset, TRUTH_QUAD_NROWS, NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT, TRUTH_PROD_NODES, TRUTH_PROD_RULE);

  _truth_cmp ("SeriesLensed vs Quad", v_lensed, v_quad, 1.0e-4, 5.0e-4);

  ncm_vector_free (v_lensed);
  ncm_vector_free (v_quad);
  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, dcwlf_lensed);
  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, dcwlf_quad);

  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_mass_summary_free (hms);
  nc_halo_density_profile_free (dp);
  nc_halo_position_free (hp);
  nc_wl_surface_mass_density_free (smd);
  nc_galaxy_redshift_pop_lsst_srd_free (pop_z);
  nc_galaxy_redshift_obs_gauss_free (obs_z);
  nc_galaxy_shape_pop_gauss_free (pop_shape);
  nc_galaxy_position_factor_free (position_factor);
  nc_galaxy_redshift_factor_free (redshift_factor);
  nc_galaxy_shape_factor_free (shape_factor_gen);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_data_unref (s_data);
  nc_galaxy_wl_obs_free (obs);
  g_strfreev (columns_strv);
  g_strv_builder_unref (builder);
  g_list_free_full (columns, g_free);
  ncm_mset_clear (&mset);
  ncm_rng_free (rng);
}

