/***************************************************************************
 *            test_ncm_stats_dist1d_spline.c
 *
 *  Tue September 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_stats_dist1d_spline.c
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

#define TEST_NKNOTS 101

/* A cubic spline of f on n equally spaced knots in [a, b] */
static NcmSpline *
_test_spline_new (gdouble ( *f ) (gdouble), gdouble a, gdouble b, guint n)
{
  NcmVector *xv = ncm_vector_new (n);
  NcmVector *yv = ncm_vector_new (n);
  NcmSpline *s;
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble x = a + (b - a) * i / (n - 1.0);

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, f (x));
  }

  s = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return s;
}

static gdouble
_test_x2 (gdouble x)
{
  return x * x;
}

static gdouble
_test_x2_p1600 (gdouble x)
{
  return x * x + 1600.0;
}

static gdouble
_test_x2_m1600 (gdouble x)
{
  return x * x - 1600.0;
}

/* Concave at +-3: slope 0.6 outward, second derivative -3.4 */
static gdouble
_test_concave (gdouble x)
{
  return x * x - 0.05 * gsl_pow_4 (x);
}

/* On [0, 3] the right edge has outward slope -2 */
static gdouble
_test_inward (gdouble x)
{
  return gsl_pow_2 (x - 4.0);
}

static gdouble
_test_gauss_density (gdouble x)
{
  return exp (-0.5 * x * x);
}

/*
 * m2lnp = x^2 on [-5, 5], reproduced exactly by the cubic spline: the unit Gaussian
 * truncated to [-5, 5]. At reltol 1e-8 the normalization is 5.8e-9 from its value, and an
 * offset of +-1600 in m2lnp moves the quantiles by 4.1e-10 (rounding in m2lnp - m_min).
 */
static void
test_ncm_stats_dist1d_spline_gauss (void)
{
  gdouble (*f[3]) (gdouble) = {&_test_x2, &_test_x2_p1600, &_test_x2_m1600};

  const gdouble Z  = sqrt (2.0 * M_PI) * erf (5.0 / M_SQRT2);
  const gdouble P5 = gsl_cdf_ugaussian_P (-5.0);
  gdouble ref[3];
  guint k;

  for (k = 0; k < 3; k++)
  {
    NcmSpline *m2lnp          = _test_spline_new (f[k], -5.0, 5.0, TEST_NKNOTS);
    NcmStatsDist1dSpline *sds = ncm_stats_dist1d_spline_new (m2lnp);
    NcmStatsDist1d *sd1       = NCM_STATS_DIST1D (sds);
    guint i;

    g_object_set (sds, "reltol", 1.0e-8, NULL);
    ncm_stats_dist1d_set_xi (sd1, -100.0);
    ncm_stats_dist1d_prepare (sd1);

    /* The support is the knot range */
    g_assert_cmpfloat (ncm_stats_dist1d_get_xi (sd1), ==, -5.0);
    g_assert_cmpfloat (ncm_stats_dist1d_get_xf (sd1), ==, 5.0);

    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_norma (sd1), ==, Z, 5.0e-8, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_mode (sd1), ==, 0.0, 0.0, 1.0e-6);

    for (i = 1; i < 100; i++)
    {
      const gdouble u = i / 100.0;
      const gdouble x = ncm_stats_dist1d_eval_inv_pdf (sd1, u);

      ncm_assert_cmpdouble_e ((gsl_cdf_ugaussian_P (x) - P5) / (1.0 - 2.0 * P5), ==, u, 0.0, 5.0e-8);
    }

    ref[k] = ncm_stats_dist1d_eval_inv_pdf (sd1, 0.975);
    ncm_assert_cmpdouble_e (ref[k], ==, ref[0], 5.0e-9, 0.0);

    ncm_stats_dist1d_free (sd1);
    ncm_spline_free (m2lnp);
  }
}

/*
 * A mode exactly at zero: m2lnp = x^2 on three knot sets whose refinement, with no
 * absolute tolerance, crawled towards zero or stopped at a single unchanged bracket. The
 * mode is found within the absolute tolerance sqrt (reltol) dx, dx the grid spacing.
 */
static void
test_ncm_stats_dist1d_spline_mode_at_zero (void)
{
  const gdouble lo[3] = {-5.0, -7.0, -5.0};
  const guint n[3]    = {51, 101, 201};
  guint k;

  for (k = 0; k < 3; k++)
  {
    NcmSpline *m2lnp    = _test_spline_new (&_test_x2, lo[k], lo[k] + 10.0, n[k]);
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (m2lnp));

    g_object_set (sd1, "reltol", 1.0e-8, NULL);
    ncm_stats_dist1d_set_compute_cdf (sd1, FALSE);
    ncm_stats_dist1d_prepare (sd1);

    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_mode (sd1), ==, 0.0, 0.0, 1.0e-4 * 10.0 / 999.0);

    ncm_stats_dist1d_free (sd1);
    ncm_spline_free (m2lnp);
  }
}

/* Continuity of the value and of the outward slope at a bound */
static void
_test_assert_c1 (NcmStatsDist1d *sd1, gdouble xb, gdouble outward, gdouble L)
{
  const gdouble h     = 1.0e-6 * L;
  const gdouble m_b   = ncm_stats_dist1d_eval_m2lnp (sd1, xb);
  const gdouble s_in  = outward * (m_b - ncm_stats_dist1d_eval_m2lnp (sd1, xb - outward * h)) / h;
  const gdouble s_out = outward * (ncm_stats_dist1d_eval_m2lnp (sd1, xb + outward * h) - m_b) / h;

  ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_m2lnp (sd1, xb + outward * 1.0e-12 * L), ==, m_b, 0.0, 1.0e-9);
  ncm_assert_cmpdouble_e (s_out, ==, s_in, 1.0e-3, 1.0e-3);
}

/*
 * Tails: exact quadratic continuation of x^2; a concave edge continues with curvature
 * 1/L^2 and grows; an edge with inward slope dips by at most 1 and grows.
 */
static void
test_ncm_stats_dist1d_spline_tails (void)
{
  {
    NcmSpline *m2lnp    = _test_spline_new (&_test_x2, -5.0, 5.0, TEST_NKNOTS);
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (m2lnp));

    ncm_stats_dist1d_prepare (sd1);
    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_m2lnp (sd1, 7.0), ==, 49.0, 1.0e-10, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_m2lnp (sd1, -9.0), ==, 81.0, 1.0e-10, 0.0);

    ncm_stats_dist1d_free (sd1);
    ncm_spline_free (m2lnp);
  }

  {
    NcmSpline *m2lnp    = _test_spline_new (&_test_concave, -3.0, 3.0, TEST_NKNOTS);
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (m2lnp));
    const gdouble L     = 6.0;
    gdouble last;
    guint i;

    ncm_stats_dist1d_prepare (sd1);
    last = ncm_stats_dist1d_eval_m2lnp (sd1, 3.0);
    _test_assert_c1 (sd1, 3.0, +1.0, L);
    _test_assert_c1 (sd1, -3.0, -1.0, L);

    for (i = 1; i <= 300; i++)
    {
      const gdouble m = ncm_stats_dist1d_eval_m2lnp (sd1, 3.0 + 0.01 * L * i);

      g_assert_cmpfloat (m, >, last);
      last = m;
    }

    /* m2lnp(3 + d) = m2lnp(3) + m2lnp'(3) d + d^2 / (2 L^2), with the spline's slope */
    ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_m2lnp (sd1, 3.0 + L), ==, ncm_spline_eval (m2lnp, 3.0) + ncm_spline_eval_deriv (m2lnp, 3.0) * L + 0.5, 1.0e-12, 0.0);

    ncm_stats_dist1d_free (sd1);
    ncm_spline_free (m2lnp);
  }

  {
    NcmSpline *m2lnp    = _test_spline_new (&_test_inward, 0.0, 3.0, TEST_NKNOTS);
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (m2lnp));
    gdouble min_out     = GSL_POSINF;
    guint i;

    ncm_stats_dist1d_prepare (sd1);
    _test_assert_c1 (sd1, 3.0, +1.0, 3.0);

    for (i = 0; i <= 1000; i++)
      min_out = GSL_MIN (min_out, ncm_stats_dist1d_eval_m2lnp (sd1, 3.0 + 0.01 * i));

    g_assert_cmpfloat (min_out, >=, ncm_stats_dist1d_eval_m2lnp (sd1, 3.0) - 1.0 - 1.0e-10);
    g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (sd1, 13.0), >, ncm_stats_dist1d_eval_m2lnp (sd1, 3.0) + 10.0);

    ncm_stats_dist1d_free (sd1);
    ncm_spline_free (m2lnp);
  }
}

/* The density path: the truncated unit Gaussian from a spline of p, zero outside the knots */
static void
test_ncm_stats_dist1d_spline_density (void)
{
  NcmSpline *p        = _test_spline_new (&_test_gauss_density, -5.0, 5.0, 401);
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new_from_density (p));
  const gdouble P5    = gsl_cdf_ugaussian_P (-5.0);
  guint i;

  g_object_set (sd1, "reltol", 1.0e-8, NULL);
  ncm_stats_dist1d_prepare (sd1);

  ncm_assert_cmpdouble_e (ncm_stats_dist1d_eval_norma (sd1), ==, sqrt (2.0 * M_PI) * erf (5.0 / M_SQRT2), 1.0e-6, 0.0);

  for (i = 1; i < 100; i++)
  {
    const gdouble u = i / 100.0;
    const gdouble x = ncm_stats_dist1d_eval_inv_pdf (sd1, u);

    ncm_assert_cmpdouble_e ((gsl_cdf_ugaussian_P (x) - P5) / (1.0 - 2.0 * P5), ==, u, 0.0, 1.0e-6);
  }

  g_assert_cmpfloat (ncm_stats_dist1d_eval_p (sd1, 5.5), ==, 0.0);
  g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (sd1, -5.5), ==, GSL_POSINF);

  ncm_stats_dist1d_free (sd1);
  ncm_spline_free (p);
}

/*
 * The density path on the HSC PDR1 multi-peak P(z) at reltol 1e-5. The reference is the
 * integral of the P(z) spline, whose negative lobes between zero knots the density path
 * removes; measured 8.7e-5 in probability.
 */
static void
test_ncm_stats_dist1d_spline_hsc_pz (void)
{
  gchar *path        = ncm_cfg_get_data_filename ("hsc_pdr1_pz_multipeak.bin", TRUE);
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmObjDictStr *ods = ncm_serialize_dict_str_from_binfile (ser, path);
  GStrv keys         = ncm_obj_dict_str_keys (ods);
  gdouble max_err    = 0.0;
  guint k;

  for (k = 0; keys[k] != NULL; k++)
  {
    NcmSpline *pz       = NCM_SPLINE (ncm_obj_dict_str_get (ods, keys[k]));
    NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new_from_density (pz));
    gdouble zmin, zmax, norm;
    guint i;

    g_object_set (sd1, "reltol", 1.0e-5, NULL);
    ncm_stats_dist1d_prepare (sd1);
    ncm_spline_get_bounds (pz, &zmin, &zmax);
    norm = ncm_spline_eval_integ (pz, zmin, zmax);

    for (i = 1; i < 1000; i++)
    {
      const gdouble u = i / 1000.0;
      const gdouble z = ncm_stats_dist1d_eval_inv_pdf (sd1, u);

      max_err = GSL_MAX (max_err, fabs (ncm_spline_eval_integ (pz, zmin, z) / norm - u));
    }

    ncm_stats_dist1d_free (sd1);
    ncm_spline_free (pz);
  }

  g_test_message ("max |F(inv(u)) - u| = %e", max_err);
  g_assert_cmpfloat (max_err, <, 1.0e-3);

  g_free (keys);
  ncm_obj_dict_str_unref (ods);
  ncm_serialize_free (ser);
  g_free (path);
}

/* A serialized copy prepares and evaluates like the original */
static void
test_ncm_stats_dist1d_spline_serialize (void)
{
  NcmSpline *m2lnp          = _test_spline_new (&_test_x2, -5.0, 5.0, TEST_NKNOTS);
  NcmStatsDist1dSpline *sds = ncm_stats_dist1d_spline_new (m2lnp);
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmStatsDist1d *dup       = NCM_STATS_DIST1D (ncm_serialize_dup_obj (ser, G_OBJECT (sds)));

  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (sds));
  ncm_stats_dist1d_prepare (dup);

  g_assert_cmpfloat (ncm_stats_dist1d_eval_inv_pdf (dup, 0.3), ==, ncm_stats_dist1d_eval_inv_pdf (NCM_STATS_DIST1D (sds), 0.3));
  g_assert_cmpfloat (ncm_stats_dist1d_eval_m2lnp (dup, 6.0), ==, ncm_stats_dist1d_eval_m2lnp (NCM_STATS_DIST1D (sds), 6.0));

  ncm_stats_dist1d_free (dup);
  ncm_stats_dist1d_free (NCM_STATS_DIST1D (sds));
  ncm_serialize_free (ser);
  ncm_spline_free (m2lnp);
}

static void
test_ncm_stats_dist1d_spline_none_subprocess (void)
{
  NcmStatsDist1d *sd1 = g_object_new (NCM_TYPE_STATS_DIST1D_SPLINE, NULL);

  ncm_stats_dist1d_prepare (sd1);
}

static void
test_ncm_stats_dist1d_spline_both_subprocess (void)
{
  NcmSpline *s        = _test_spline_new (&_test_x2, -5.0, 5.0, TEST_NKNOTS);
  NcmStatsDist1d *sd1 = g_object_new (NCM_TYPE_STATS_DIST1D_SPLINE, "m2lnp", s, "density", s, NULL);

  ncm_stats_dist1d_prepare (sd1);
}

static void
test_ncm_stats_dist1d_spline_traps (void)
{
  g_test_trap_subprocess ("/ncm/stats/dist1d_spline/none/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exactly one of the m2lnp and density splines must be set*");

  g_test_trap_subprocess ("/ncm/stats/dist1d_spline/both/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exactly one of the m2lnp and density splines must be set*");
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/stats/dist1d_spline/gauss", &test_ncm_stats_dist1d_spline_gauss);
  g_test_add_func ("/ncm/stats/dist1d_spline/mode_at_zero", &test_ncm_stats_dist1d_spline_mode_at_zero);
  g_test_add_func ("/ncm/stats/dist1d_spline/tails", &test_ncm_stats_dist1d_spline_tails);
  g_test_add_func ("/ncm/stats/dist1d_spline/density", &test_ncm_stats_dist1d_spline_density);
  g_test_add_func ("/ncm/stats/dist1d_spline/hsc_pz", &test_ncm_stats_dist1d_spline_hsc_pz);
  g_test_add_func ("/ncm/stats/dist1d_spline/serialize", &test_ncm_stats_dist1d_spline_serialize);
  g_test_add_func ("/ncm/stats/dist1d_spline/traps", &test_ncm_stats_dist1d_spline_traps);
  g_test_add_func ("/ncm/stats/dist1d_spline/none/subprocess", &test_ncm_stats_dist1d_spline_none_subprocess);
  g_test_add_func ("/ncm/stats/dist1d_spline/both/subprocess", &test_ncm_stats_dist1d_spline_both_subprocess);

  g_test_run ();
}

