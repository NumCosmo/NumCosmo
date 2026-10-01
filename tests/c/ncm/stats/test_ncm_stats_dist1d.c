/***************************************************************************
 *            test_ncm_stats_dist1d.c
 *
 *  Tue September 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_stats_dist1d.c
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_cdf.h>

/*
 * Toy distribution: the unnormalized Gaussian p(x) = A exp(-(x - mu)^2 / (2 sigma^2))
 * on [xi, xf], whose normalization, cumulative distribution and quantiles have closed
 * forms through the unit Gaussian P and Q.
 */

#define TEST_TYPE_GAUSS1D (test_gauss1d_get_type ())
G_DECLARE_FINAL_TYPE (TestGauss1d, test_gauss1d, TEST, GAUSS1D, NcmStatsDist1d)

struct _TestGauss1d
{
  NcmStatsDist1d parent_instance;
  gdouble A;
  gdouble mu;
  gdouble sigma;
};

G_DEFINE_TYPE (TestGauss1d, test_gauss1d, NCM_TYPE_STATS_DIST1D)

static void
test_gauss1d_init (TestGauss1d *tg)
{
  tg->A     = 1.0;
  tg->mu    = 0.0;
  tg->sigma = 1.0;
}

static gdouble
_test_gauss1d_p (NcmStatsDist1d *sd1, gdouble x)
{
  TestGauss1d *tg = TEST_GAUSS1D (sd1);

  return tg->A * exp (-0.5 * gsl_pow_2 ((x - tg->mu) / tg->sigma));
}

static gdouble
_test_gauss1d_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  TestGauss1d *tg = TEST_GAUSS1D (sd1);

  return gsl_pow_2 ((x - tg->mu) / tg->sigma) - 2.0 * log (tg->A);
}

static void
test_gauss1d_class_init (TestGauss1dClass *klass)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_CLASS (klass);

  sd1_class->p     = &_test_gauss1d_p;
  sd1_class->m2lnp = &_test_gauss1d_m2lnp;
}

static NcmStatsDist1d *
_test_gauss1d_new (gdouble A, gdouble mu, gdouble sigma, gdouble xi, gdouble xf, gdouble reltol)
{
  TestGauss1d *tg = g_object_new (TEST_TYPE_GAUSS1D,
                                  "xi", xi,
                                  "xf", xf,
                                  "reltol", reltol,
                                  NULL);

  tg->A     = A;
  tg->mu    = mu;
  tg->sigma = sigma;

  return NCM_STATS_DIST1D (tg);
}

/* Probability of [xi, x] under the truncated Gaussian */
static gdouble
_test_gauss1d_cdf (TestGauss1d *tg, gdouble xi, gdouble xf, gdouble x)
{
  const gdouble Pi = gsl_cdf_gaussian_P (xi - tg->mu, tg->sigma);
  const gdouble Pf = gsl_cdf_gaussian_P (xf - tg->mu, tg->sigma);

  return (gsl_cdf_gaussian_P (x - tg->mu, tg->sigma) - Pi) / (Pf - Pi);
}

/* Probability of [x, xf] under the truncated Gaussian, accurate in the upper tail */
static gdouble
_test_gauss1d_tail (TestGauss1d *tg, gdouble xi, gdouble xf, gdouble x)
{
  const gdouble Qi = gsl_cdf_gaussian_Q (xi - tg->mu, tg->sigma);
  const gdouble Qf = gsl_cdf_gaussian_Q (xf - tg->mu, tg->sigma);

  return (gsl_cdf_gaussian_Q (x - tg->mu, tg->sigma) - Qf) / (Qi - Qf);
}

/*
 * Toy distribution sampling a tabulated density: p(x) = max(P(x), 0) for a spline P, the
 * photometric-redshift densities of the HSC PDR1 catalogs in the tests below.
 */

#define TEST_TYPE_SPLINE1D (test_spline1d_get_type ())
G_DECLARE_FINAL_TYPE (TestSpline1d, test_spline1d, TEST, SPLINE1D, NcmStatsDist1d)

struct _TestSpline1d
{
  NcmStatsDist1d parent_instance;
  NcmSpline *P;
};

G_DEFINE_TYPE (TestSpline1d, test_spline1d, NCM_TYPE_STATS_DIST1D)

static void
test_spline1d_init (TestSpline1d *ts)
{
  ts->P = NULL;
}

static void
_test_spline1d_dispose (GObject *object)
{
  ncm_spline_clear (&TEST_SPLINE1D (object)->P);

  G_OBJECT_CLASS (test_spline1d_parent_class)->dispose (object);
}

static gdouble
_test_spline1d_p (NcmStatsDist1d *sd1, gdouble x)
{
  return GSL_MAX (ncm_spline_eval (TEST_SPLINE1D (sd1)->P, x), 0.0);
}

static gdouble
_test_spline1d_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  const gdouble p = _test_spline1d_p (sd1, x);

  return (p > 0.0) ? -2.0 * log (p) : GSL_POSINF;
}

static void
test_spline1d_class_init (TestSpline1dClass *klass)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_CLASS (klass);

  G_OBJECT_CLASS (klass)->dispose = &_test_spline1d_dispose;
  sd1_class->p                    = &_test_spline1d_p;
  sd1_class->m2lnp                = &_test_spline1d_m2lnp;
}

/*
 * Toy distribution with two bumps separated by an interval of zero density, on [0, 10]:
 * [1, 3] and [6.5, 7.5], either top hats of heights 1 and 3 or the smooth bumps
 * a (1 - t^2)^2, t = (x - c) / w, with the same supports and amplitudes 1 and 3.
 */

#define TEST_TYPE_BUMPS1D (test_bumps1d_get_type ())
G_DECLARE_FINAL_TYPE (TestBumps1d, test_bumps1d, TEST, BUMPS1D, NcmStatsDist1d)

struct _TestBumps1d
{
  NcmStatsDist1d parent_instance;
  gboolean smooth;
};

G_DEFINE_TYPE (TestBumps1d, test_bumps1d, NCM_TYPE_STATS_DIST1D)

static const gdouble test_bumps_lo[2] = {1.0, 6.5};
static const gdouble test_bumps_hi[2] = {3.0, 7.5};
static const gdouble test_bumps_a[2]  = {1.0, 3.0};

static void
test_bumps1d_init (TestBumps1d *tb)
{
  tb->smooth = FALSE;
}

static gdouble
_test_bumps1d_p (NcmStatsDist1d *sd1, gdouble x)
{
  TestBumps1d *tb = TEST_BUMPS1D (sd1);
  guint j;

  for (j = 0; j < 2; j++)
  {
    if ((x >= test_bumps_lo[j]) && (x <= test_bumps_hi[j]))
    {
      const gdouble w = 0.5 * (test_bumps_hi[j] - test_bumps_lo[j]);
      const gdouble t = (x - test_bumps_lo[j] - w) / w;

      return tb->smooth ? test_bumps_a[j] * gsl_pow_2 (1.0 - t * t) : test_bumps_a[j];
    }
  }

  return 0.0;
}

static gdouble
_test_bumps1d_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  const gdouble p = _test_bumps1d_p (sd1, x);

  return (p > 0.0) ? -2.0 * log (p) : GSL_POSINF;
}

static void
test_bumps1d_class_init (TestBumps1dClass *klass)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_CLASS (klass);

  sd1_class->p     = &_test_bumps1d_p;
  sd1_class->m2lnp = &_test_bumps1d_m2lnp;
}

/* Mass of bump j in [lo_j, x]; the smooth bump integrates to 16 w a / 15 */
static gdouble
_test_bumps1d_mass (TestBumps1d *tb, guint j, gdouble x)
{
  const gdouble w = 0.5 * (test_bumps_hi[j] - test_bumps_lo[j]);
  const gdouble t = GSL_MIN (GSL_MAX ((x - test_bumps_lo[j] - w) / w, -1.0), 1.0);

  if (tb->smooth)
    return test_bumps_a[j] * w * (t - 2.0 * gsl_pow_3 (t) / 3.0 + gsl_pow_5 (t) / 5.0 + 8.0 / 15.0);
  else
    return test_bumps_a[j] * w * (t + 1.0);
}

static gdouble
_test_bumps1d_cdf (TestBumps1d *tb, gdouble x)
{
  const gdouble M = _test_bumps1d_mass (tb, 0, 10.0) + _test_bumps1d_mass (tb, 1, 10.0);

  return (_test_bumps1d_mass (tb, 0, x) + _test_bumps1d_mass (tb, 1, x)) / M;
}

#define TEST_XI (-1.0)
#define TEST_XF 4.0
#define TEST_MU 0.7
#define TEST_SIGMA 0.9
#define TEST_NPOINTS 201

/*
 * At reltol 1e-8 the measured errors are 1.2e-9 (normalization, relative), 3.6e-9
 * (cumulative distribution) and 4.7e-9 (quantiles, in probability), the same for density
 * scales A from 1e-300 to 1e300 and for x scales from 1e-100 to 1e100, whenever the
 * normalization A sigma sqrt(2 pi) is a normal double.
 */
static void
test_ncm_stats_dist1d_gauss (void)
{
  const gdouble A_list[] = {1.0e-300, 1.0e-20, 1.0, 1.0e20, 1.0e300};
  const gdouble s_list[] = {1.0e-100, 1.0e-10, 1.0, 1.0e10, 1.0e100};
  guint a, b;

  for (a = 0; a < G_N_ELEMENTS (A_list); a++)
  {
    for (b = 0; b < G_N_ELEMENTS (s_list); b++)
    {
      const gdouble A     = A_list[a];
      const gdouble sc    = s_list[b];
      const gdouble xi    = TEST_XI * sc;
      const gdouble xf    = TEST_XF * sc;
      const gdouble mu    = TEST_MU * sc;
      const gdouble sigma = TEST_SIGMA * sc;
      NcmStatsDist1d *sd1 = _test_gauss1d_new (A, mu, sigma, xi, xf, 1.0e-8);
      TestGauss1d *tg     = TEST_GAUSS1D (sd1);
      const gdouble Z     = sqrt (2.0 * M_PI) * sigma *
                            (gsl_cdf_gaussian_P (xf - mu, sigma) - gsl_cdf_gaussian_P (xi - mu, sigma));
      guint i;

      if (fabs (log10 (A * sc)) > 300.0)
      {
        ncm_stats_dist1d_free (sd1);
        continue;
      }

      ncm_stats_dist1d_prepare (sd1);

      ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_norma (sd1) / Z, ==, A, 1.0e-8, 0.0);
      ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_mode (sd1), ==, mu, 0.0, 1.0e-6 * sc);

      for (i = 0; i < TEST_NPOINTS; i++)
      {
        const gdouble x = xi + (xf - xi) * i / (TEST_NPOINTS - 1.0);
        const gdouble u = 1.0e-6 + (1.0 - 2.0e-6) * i / (TEST_NPOINTS - 1.0);
        const gdouble v = pow (10.0, -12.0 + 11.0 * i / (TEST_NPOINTS - 1.0));

        ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_p (sd1, x) * Z, ==, exp (-0.5 * gsl_pow_2 ((x - mu) / sigma)), 0.0, 1.0e-8);
        ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_pdf (sd1, x), ==, _test_gauss1d_cdf (tg, xi, xf, x), 0.0, 5.0e-8);
        ncm_assert_cmpdouble_e (_test_gauss1d_cdf (tg, xi, xf, ncm_stats_dist1d_eval_inv_pdf (sd1, u)), ==, u, 0.0, 5.0e-8);
        ncm_assert_cmpdouble_e (_test_gauss1d_tail (tg, xi, xf, ncm_stats_dist1d_eval_inv_pdf_tail (sd1, v)), ==, v, 0.0, 5.0e-8);
      }

      g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf (sd1, 0.0), ==, xi);
      g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf (sd1, 1.0), ==, xf);
      g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf_tail (sd1, 0.0), ==, xf);
      g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf_tail (sd1, 1.0), ==, xi);

      ncm_stats_dist1d_free (sd1);
    }
  }
}

/* A support starting at zero with the default zero abstol */
static void
test_ncm_stats_dist1d_zero_xi (void)
{
  NcmStatsDist1d *sd1 = _test_gauss1d_new (1.0, TEST_MU, TEST_SIGMA, 0.0, TEST_XF, 1.0e-8);
  TestGauss1d *tg     = TEST_GAUSS1D (sd1);
  guint i;

  ncm_stats_dist1d_prepare (sd1);

  for (i = 1; i < TEST_NPOINTS - 1; i++)
  {
    const gdouble u = (gdouble) i / (TEST_NPOINTS - 1.0);

    ncm_assert_cmpdouble_e (_test_gauss1d_cdf (tg, 0.0, TEST_XF, ncm_stats_dist1d_eval_inv_pdf (sd1, u)), ==, u, 0.0, 5.0e-8);
  }

  ncm_stats_dist1d_free (sd1);
}

/* A maximum at a bound returns the bound */
static void
test_ncm_stats_dist1d_mode_edge (void)
{
  NcmStatsDist1d *lower = _test_gauss1d_new (1.0, -2.0, TEST_SIGMA, TEST_XI, TEST_XF, 1.0e-8);
  NcmStatsDist1d *upper = _test_gauss1d_new (1.0, 5.0, TEST_SIGMA, TEST_XI, TEST_XF, 1.0e-8);

  g_assert_cmpfloat (ncm_stats_dist1d_eval_mode (lower), ==, TEST_XI);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_mode (upper), ==, TEST_XF);

  ncm_stats_dist1d_prepare (lower);
  ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_pdf (lower, 0.0), ==, _test_gauss1d_cdf (TEST_GAUSS1D (lower), TEST_XI, TEST_XF, 0.0), 0.0, 5.0e-8);

  ncm_stats_dist1d_free (lower);
  ncm_stats_dist1d_free (upper);
}

/* Without compute-cdf the density is not normalized */
static void
test_ncm_stats_dist1d_no_cdf (void)
{
  NcmStatsDist1d *sd1 = _test_gauss1d_new (3.0, TEST_MU, TEST_SIGMA, TEST_XI, TEST_XF, 1.0e-8);

  ncm_stats_dist1d_set_compute_cdf (sd1, FALSE);
  g_assert_false (ncm_stats_dist1d_get_compute_cdf (sd1));
  ncm_stats_dist1d_prepare (sd1);

  g_assert_cmpfloat (ncm_stats_dist1d_eval_norma (sd1), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_p (sd1, TEST_MU), ==, 3.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (sd1, TEST_MU), ==, -2.0 * log (3.0));

  ncm_stats_dist1d_free (sd1);
}

/* xi = xf is a point mass */
static void
test_ncm_stats_dist1d_point_mass (void)
{
  NcmStatsDist1d *sd1 = _test_gauss1d_new (1.0, TEST_MU, TEST_SIGMA, 1.0, 1.0, 1.0e-8);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 1);

  ncm_stats_dist1d_prepare (sd1);

  g_assert_cmpfloat (ncm_stats_dist1d_eval_p (sd1, 1.0), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_p (sd1, 2.0), ==, 0.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (sd1, 1.0), ==, 0.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (sd1, 2.0), ==, GSL_POSINF);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_pdf (sd1, 0.5), ==, 0.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_pdf (sd1, 1.0), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_norma (sd1), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf (sd1, 0.3), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf_tail (sd1, 0.3), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_mode (sd1), ==, 1.0);
  g_assert_cmpfloat (ncm_stats_dist1d_gen (sd1, rng), ==, 1.0);

  ncm_rng_free (rng);
  ncm_stats_dist1d_free (sd1);
}

/* Half of the draws fall below the median */
static void
test_ncm_stats_dist1d_gen (void)
{
  const guint ndraws  = 100000;
  NcmStatsDist1d *sd1 = _test_gauss1d_new (1.0, TEST_MU, TEST_SIGMA, TEST_XI, TEST_XF, 1.0e-8);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 1);
  guint below         = 0;
  gdouble median;
  guint i;

  ncm_stats_dist1d_prepare (sd1);
  median = ncm_stats_dist1d_eval_inv_pdf (sd1, 0.5);

  for (i = 0; i < ndraws; i++)
  {
    const gdouble x = ncm_stats_dist1d_gen (sd1, rng);

    g_assert_cmpfloat (x, >=, TEST_XI);
    g_assert_cmpfloat (x, <=, TEST_XF);

    if (x < median)
      below++;
  }

  ncm_assert_cmpdouble_e (below, ==, 0.5 * ndraws, 0.0, 5.0 * sqrt (0.25 * ndraws));

  ncm_rng_free (rng);
  ncm_stats_dist1d_free (sd1);
}

/*
 * HSC PDR1 P(z) with five to nine peaks separated by zero density, at the reltol 1e-5 of
 * NcGalaxyRedshiftFactorSpline. The reference is the integral of the P(z) spline, which
 * keeps the small negative lobes between zero knots that p(x) = max(P(x), 0) removes; the
 * difference is part of the measured 1.0e-4 below.
 */
static void
test_ncm_stats_dist1d_hsc_pz (void)
{
  gchar *path        = ncm_cfg_get_data_filename ("hsc_pdr1_pz_multipeak.bin", TRUE);
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmObjDictStr *ods = ncm_serialize_dict_str_from_binfile (ser, path);
  GStrv keys         = ncm_obj_dict_str_keys (ods);
  gdouble max_err    = 0.0;
  guint k;

  g_assert_cmpuint (ncm_obj_dict_str_len (ods), ==, 10);

  for (k = 0; keys[k] != NULL; k++)
  {
    TestSpline1d *ts    = g_object_new (TEST_TYPE_SPLINE1D, "reltol", 1.0e-5, NULL);
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ts);
    gdouble zmin, zmax, norm;
    guint i;

    ts->P = NCM_SPLINE (ncm_obj_dict_str_get (ods, keys[k]));
    ncm_spline_prepare (ts->P);
    ncm_spline_get_bounds (ts->P, &zmin, &zmax);
    norm = ncm_spline_eval_integ (ts->P, zmin, zmax);

    ncm_stats_dist1d_set_xi (sd1, zmin);
    ncm_stats_dist1d_set_xf (sd1, zmax);
    ncm_stats_dist1d_prepare (sd1);

    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_norma (sd1), ==, norm, 1.0e-3, 0.0);

    for (i = 1; i < 1000; i++)
    {
      const gdouble u = i / 1000.0;
      const gdouble z = ncm_stats_dist1d_eval_inv_pdf (sd1, u);
      const gdouble e = fabs (ncm_spline_eval_integ (ts->P, zmin, z) / norm - u);

      g_assert_cmpfloat (z, >=, zmin);
      g_assert_cmpfloat (z, <=, zmax);
      max_err = GSL_MAX (max_err, e);
    }

    ncm_stats_dist1d_free (sd1);
  }

  g_test_message ("max |F(inv(u)) - u| = %e", max_err);
  g_assert_cmpfloat (max_err, <, 1.0e-3);

  g_free (keys);
  ncm_obj_dict_str_unref (ods);
  ncm_serialize_free (ser);
  g_free (path);
}

/*
 * Two bumps separated by a zero-density interval, at reltol 1e-8. Measured: 4.0e-9 (top
 * hats) and 5.9e-7 (smooth bumps) in probability, no quantile and no draw inside the gap,
 * the mode on the taller bump (anywhere on its top for the top hat).
 */
static void
test_ncm_stats_dist1d_zero_gap (void)
{
  const gdouble bound[2] = {5.0e-8, 5.0e-6};
  const guint ndraws     = 100000;
  guint smooth;

  for (smooth = 0; smooth < 2; smooth++)
  {
    TestBumps1d *tb     = g_object_new (TEST_TYPE_BUMPS1D, "xi", 0.0, "xf", 10.0, "reltol", 1.0e-8, NULL);
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (tb);
    NcmRNG *rng         = ncm_rng_seeded_new (NULL, 1);
    guint in_first      = 0;
    gdouble m1;
    guint i;

    tb->smooth = smooth;
    m1         = _test_bumps1d_cdf (tb, 3.0);
    ncm_stats_dist1d_prepare (sd1);

    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_norma (sd1), ==, _test_bumps1d_mass (tb, 0, 10.0) + _test_bumps1d_mass (tb, 1, 10.0), bound[smooth], 0.0);

    if (smooth)
    {
      ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_mode (sd1), ==, 7.0, 0.0, 1.0e-6);
    }
    else
    {
      g_assert_cmpfloat (ncm_stats_dist1d_eval_mode (sd1), >=, 6.5);
      g_assert_cmpfloat (ncm_stats_dist1d_eval_mode (sd1), <=, 7.5);
    }

    g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (sd1, 5.0), ==, GSL_POSINF);

    for (i = 1; i < 20000; i++)
    {
      const gdouble u = i / 20000.0;
      const gdouble x = ncm_stats_dist1d_eval_inv_pdf (sd1, u);

      ncm_assert_cmpdouble_e (_test_bumps1d_cdf (tb, x), ==, u, 0.0, bound[smooth]);
      g_assert_false ((x > 3.0 + 1.0e-6) && (x < 6.5 - 1.0e-6));
      g_assert_cmpfloat (x, >=, 1.0 - 1.0e-6);
      g_assert_cmpfloat (x, <=, 7.5 + 1.0e-6);
    }

    for (i = 0; i < ndraws; i++)
    {
      const gdouble x = ncm_stats_dist1d_gen (sd1, rng);

      g_assert_false ((x > 3.0 + 1.0e-6) && (x < 6.5 - 1.0e-6));

      if (x < 4.75)
        in_first++;
    }

    ncm_assert_cmpdouble_e (in_first, ==, m1 * ndraws, 0.0, 5.0 * sqrt (m1 * (1.0 - m1) * ndraws));

    /* At reltol 1e-5 a probability above the plateau by 1e-6 or more starts after the gap (measured 8.9e-8 top hats, 4.2e-7 smooth) */
    g_object_set (sd1, "reltol", 1.0e-5, NULL);
    ncm_stats_dist1d_prepare (sd1);

    {
      const gdouble U = ncm_stats_dist1d_eval_pdf (sd1, 5.0);

      for (i = 0; i <= 50; i++)
      {
        const gdouble x = ncm_stats_dist1d_eval_inv_pdf (sd1, U + pow (10.0, -6.0 + 5.0 * i / 50.0));

        g_assert_false ((x > 3.0 + 1.0e-6) && (x < 6.5 - 1.0e-6));
      }
    }

    ncm_rng_free (rng);
    ncm_stats_dist1d_free (sd1);
  }
}

static void
test_ncm_stats_dist1d_reversed_subprocess (void)
{
  NcmStatsDist1d *sd1 = _test_gauss1d_new (1.0, TEST_MU, TEST_SIGMA, TEST_XF, TEST_XI, 1.0e-8);

  ncm_stats_dist1d_prepare (sd1);
}

static void
test_ncm_stats_dist1d_zero_density_subprocess (void)
{
  NcmStatsDist1d *sd1 = _test_gauss1d_new (0.0, TEST_MU, TEST_SIGMA, TEST_XI, TEST_XF, 1.0e-8);

  ncm_stats_dist1d_prepare (sd1);
}

static void
test_ncm_stats_dist1d_no_h_subprocess (void)
{
  NcmStatsDist1d *sd1 = _test_gauss1d_new (1.0, TEST_MU, TEST_SIGMA, TEST_XI, TEST_XF, 1.0e-8);

  ncm_stats_dist1d_get_current_h (sd1);
}

static void
test_ncm_stats_dist1d_traps (void)
{
  g_test_trap_subprocess ("/ncm/stats/dist1d/reversed/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*`TestGauss1d' has xf =*below xi =*");

  g_test_trap_subprocess ("/ncm/stats/dist1d/zero_density/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*`TestGauss1d' has density*at its mode*");

  g_test_trap_subprocess ("/ncm/stats/dist1d/no_h/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*`TestGauss1d' does not implement get_current_h*");
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/stats/dist1d/gauss", &test_ncm_stats_dist1d_gauss);
  g_test_add_func ("/ncm/stats/dist1d/zero_xi", &test_ncm_stats_dist1d_zero_xi);
  g_test_add_func ("/ncm/stats/dist1d/mode_edge", &test_ncm_stats_dist1d_mode_edge);
  g_test_add_func ("/ncm/stats/dist1d/no_cdf", &test_ncm_stats_dist1d_no_cdf);
  g_test_add_func ("/ncm/stats/dist1d/point_mass", &test_ncm_stats_dist1d_point_mass);
  g_test_add_func ("/ncm/stats/dist1d/gen", &test_ncm_stats_dist1d_gen);
  g_test_add_func ("/ncm/stats/dist1d/hsc_pz", &test_ncm_stats_dist1d_hsc_pz);
  g_test_add_func ("/ncm/stats/dist1d/zero_gap", &test_ncm_stats_dist1d_zero_gap);
  g_test_add_func ("/ncm/stats/dist1d/traps", &test_ncm_stats_dist1d_traps);
  g_test_add_func ("/ncm/stats/dist1d/reversed/subprocess", &test_ncm_stats_dist1d_reversed_subprocess);
  g_test_add_func ("/ncm/stats/dist1d/zero_density/subprocess", &test_ncm_stats_dist1d_zero_density_subprocess);
  g_test_add_func ("/ncm/stats/dist1d/no_h/subprocess", &test_ncm_stats_dist1d_no_h_subprocess);

  g_test_run ();
}

