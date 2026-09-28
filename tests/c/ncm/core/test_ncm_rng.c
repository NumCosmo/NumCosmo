/***************************************************************************
 *            test_ncm_rng.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_rng.c
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
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_rng.h>

#define TEST_RNG_SEED 123456
#define TEST_RNG_NDRAWS 100000

/* Mean of the unit Gaussian restricted to x > 1: 2 phi(1) / erfc(1 / sqrt(2)) */
#define TEST_RNG_TAIL_MEAN (M_SQRT2 / M_SQRTPI * exp (-0.5) / gsl_sf_erfc (M_SQRT1_2))

static void
test_ncm_rng_seed_state (void)
{
  NcmRNG *rng1 = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);
  NcmRNG *rng2 = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);
  gdouble draws[5];
  gchar *state;
  guint i;

  g_assert_cmpstr (ncm_rng_get_algo (rng1), ==, "mt19937");
  g_assert_cmpuint (ncm_rng_get_seed (rng1), ==, TEST_RNG_SEED);
  g_assert_false (ncm_rng_check_seed (rng1, TEST_RNG_SEED));

  /* Same algorithm and seed, same sequence */
  for (i = 0; i < 5; i++)
    g_assert_cmpuint (ncm_rng_gen_ulong (rng1), ==, ncm_rng_gen_ulong (rng2));

  /* Restoring a state repeats the sequence from it */
  state = ncm_rng_get_state (rng1);

  for (i = 0; i < 5; i++)
    draws[i] = ncm_rng_uniform01_gen (rng1);

  ncm_rng_set_state (rng2, state);

  for (i = 0; i < 5; i++)
    g_assert_cmpfloat (ncm_rng_uniform01_gen (rng2), ==, draws[i]);

  {
    gchar *state_prop;
    gchar *state1 = ncm_rng_get_state (rng1);

    g_object_get (rng2, "state", &state_prop, NULL);
    g_assert_cmpstr (state_prop, ==, state1);
    g_free (state_prop);
    g_free (state1);
  }

  ncm_rng_set_random_seed (rng2, FALSE);
  g_assert_false (ncm_rng_check_seed (rng2, ncm_rng_get_seed (rng2)));

  ncm_rng_lock (rng1);
  ncm_rng_unlock (rng1);

  g_free (state);
  ncm_rng_free (rng1);
  ncm_rng_free (rng2);
}

static void
test_ncm_rng_pool (void)
{
  NcmRNG *rng_a  = ncm_rng_pool_get ("test_ncm_rng_pool_a");
  NcmRNG *rng_a2 = ncm_rng_pool_get ("test_ncm_rng_pool_a");
  NcmRNG *rng_b  = ncm_rng_pool_get ("test_ncm_rng_pool_b");

  g_assert_true (rng_a == rng_a2);
  g_assert_true (rng_a != rng_b);

  ncm_rng_free (rng_a);
  ncm_rng_free (rng_a2);
  ncm_rng_free (rng_b);
}

typedef gdouble (*TestNcmRNGDraw) (NcmRNG *rng);

static gdouble
_test_gaussian (NcmRNG *rng)
{
  return ncm_rng_gaussian_gen (rng, 1.5, 2.0);
}

static gdouble
_test_ugaussian (NcmRNG *rng)
{
  return ncm_rng_ugaussian_gen (rng);
}

static gdouble
_test_uniform (NcmRNG *rng)
{
  return ncm_rng_uniform_gen (rng, -1.0, 3.0);
}

static gdouble
_test_exponential (NcmRNG *rng)
{
  return ncm_rng_exponential_gen (rng, 2.5);
}

static gdouble
_test_laplace (NcmRNG *rng)
{
  return ncm_rng_laplace_gen (rng, 1.5);
}

static gdouble
_test_exppow (NcmRNG *rng)
{
  return ncm_rng_exppow_gen (rng, 1.0, 2.0);
}

static gdouble
_test_beta (NcmRNG *rng)
{
  return ncm_rng_beta_gen (rng, 2.0, 3.0);
}

static gdouble
_test_gamma (NcmRNG *rng)
{
  return ncm_rng_gamma_gen (rng, 3.0, 0.5);
}

static gdouble
_test_chisq (NcmRNG *rng)
{
  return ncm_rng_chisq_gen (rng, 4.0);
}

static gdouble
_test_poisson (NcmRNG *rng)
{
  return ncm_rng_poisson_gen (rng, 3.5);
}

static gdouble
_test_rayleigh (NcmRNG *rng)
{
  return ncm_rng_rayleigh_gen (rng, 2.0);
}

static gdouble
_test_gaussian_tail (NcmRNG *rng)
{
  return ncm_rng_gaussian_tail_gen (rng, 1.0, 1.0);
}

static void
test_ncm_rng_distributions (void)
{
  NcmRNG *rng = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);

  /* {draw, mean, variance, lower bound} */
  const struct
  {
    TestNcmRNGDraw draw;
    gdouble mean;
    gdouble var;
    gdouble lb;
  } cases[] = {
    {&_test_gaussian, 1.5, 4.0, -GSL_POSINF},
    {&_test_ugaussian, 0.0, 1.0, -GSL_POSINF},
    {&_test_uniform, 1.0, 16.0 / 12.0, -1.0},
    {&_test_exponential, 2.5, 6.25, 0.0},
    {&_test_laplace, 0.0, 2.0 * 1.5 * 1.5, -GSL_POSINF},
    {&_test_exppow, 0.0, 0.5, -GSL_POSINF},
    {&_test_beta, 0.4, 6.0 / (25.0 * 6.0), 0.0},
    {&_test_gamma, 1.5, 0.75, 0.0},
    {&_test_chisq, 4.0, 8.0, 0.0},
    {&_test_poisson, 3.5, 3.5, 0.0},
    {&_test_rayleigh, 2.0 * M_SQRTPI / M_SQRT2, (4.0 - M_PI) / 2.0 * 4.0, 0.0},
    {&_test_gaussian_tail, TEST_RNG_TAIL_MEAN, 1.0 + TEST_RNG_TAIL_MEAN - TEST_RNG_TAIL_MEAN * TEST_RNG_TAIL_MEAN, 1.0},
  };

  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    gdouble sum = 0.0;
    gdouble sd;

    for (j = 0; j < TEST_RNG_NDRAWS; j++)
    {
      const gdouble x = cases[i].draw (rng);

      g_assert_cmpfloat (x, >=, cases[i].lb);
      sum += x;
    }

    sd = sqrt (cases[i].var);
    ncm_assert_cmpdouble_e (sum / TEST_RNG_NDRAWS, ==, cases[i].mean, 0.0, 5.0 * sd / sqrt (TEST_RNG_NDRAWS));
  }

  ncm_rng_free (rng);
}

static void
test_ncm_rng_uniform_ranges (void)
{
  NcmRNG *rng = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);
  guint j;

  for (j = 0; j < TEST_RNG_NDRAWS; j++)
  {
    const gdouble u  = ncm_rng_uniform01_gen (rng);
    const gdouble up = ncm_rng_uniform01_pos_gen (rng);

    g_assert_cmpfloat (u, >=, 0.0);
    g_assert_cmpfloat (u, <, 1.0);
    g_assert_cmpfloat (up, >, 0.0);
    g_assert_cmpfloat (up, <, 1.0);
    g_assert_cmpuint (ncm_rng_uniform_int_gen (rng, 7), <, 7);
  }

  ncm_rng_free (rng);
}

static void
test_ncm_rng_discrete (void)
{
  NcmRNG *rng                   = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);
  const gdouble weights[]       = {1.0, 2.0, 7.0};
  NcmRNGDiscrete *rng_discrete  = ncm_rng_discrete_new (weights, 3);
  NcmRNGDiscrete *rng_discrete2 = ncm_rng_discrete_copy (rng_discrete);
  guint counts[3]               = {0, 0, 0};
  guint i, j;

  for (j = 0; j < TEST_RNG_NDRAWS; j++)
  {
    const gsize k = ncm_rng_discrete_gen (rng, (j % 2) ? rng_discrete : rng_discrete2);

    g_assert_cmpuint (k, <, 3);
    counts[k]++;
  }

  for (i = 0; i < 3; i++)
  {
    const gdouble p = weights[i] / 10.0;

    ncm_assert_cmpdouble_e (counts[i] * 1.0 / TEST_RNG_NDRAWS, ==, p, 0.0, 5.0 * sqrt (p * (1.0 - p) / TEST_RNG_NDRAWS));
  }

  ncm_rng_discrete_free (rng_discrete);
  ncm_rng_discrete_free (rng_discrete2);
  ncm_rng_free (rng);
}

static void
test_ncm_rng_sample_choose (void)
{
  NcmRNG *rng  = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);
  gint src[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  gint dest[20];
  guint i, j;

  /* Without replacement, in the order of src */
  for (j = 0; j < 100; j++)
  {
    ncm_rng_choose (rng, dest, 6, src, 10, sizeof (gint));

    for (i = 1; i < 6; i++)
      g_assert_cmpint (dest[i], >, dest[i - 1]);
  }

  /* With replacement: more draws than elements */
  ncm_rng_sample (rng, dest, 20, src, 10, sizeof (gint));

  for (i = 0; i < 20; i++)
  {
    g_assert_cmpint (dest[i], >=, 0);
    g_assert_cmpint (dest[i], <, 10);
  }

  ncm_rng_free (rng);
}

static void
test_ncm_rng_multivariate (void)
{
  NcmRNG *rng       = ncm_rng_seeded_new ("mt19937", TEST_RNG_SEED);
  const gdouble p[] = {0.2, 0.3, 0.5};
  const gdouble rho = 0.6;
  guint n[3];
  gdouble sxy = 0.0;
  guint j;

  ncm_rng_multinomial (rng, 3, 1000, p, n);
  g_assert_cmpuint (n[0] + n[1] + n[2], ==, 1000);

  for (j = 0; j < TEST_RNG_NDRAWS; j++)
  {
    gdouble x, y;

    ncm_rng_bivariate_gaussian_gen (rng, 2.0, 0.5, rho, &x, &y);
    sxy += x * y;
  }

  /* For unit variances, Var(x y) = 1 + rho^2 */
  ncm_assert_cmpdouble_e (sxy / TEST_RNG_NDRAWS / (2.0 * 0.5), ==, rho, 0.0, 5.0 * sqrt (1.0 + rho * rho) / sqrt (TEST_RNG_NDRAWS));

  ncm_rng_free (rng);
}

static void
test_ncm_rng_invalid_algo_subprocess (void)
{
  NcmRNG *rng = ncm_rng_new (NULL);

  ncm_rng_set_algo (rng, "not-a-gsl-algorithm");
}

static void
test_ncm_rng_traps (void)
{
  g_test_trap_subprocess ("/ncm/rng/invalid_algo/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

/* Changing the algorithm keeps the seed; NULL selects the GSL default */
static void
test_ncm_rng_set_algo (void)
{
  NcmRNG *rng   = ncm_rng_new ("mt19937");
  NcmRNG *fresh = ncm_rng_new ("ranlxd2");
  guint i;

  ncm_rng_set_seed (rng, 1234);
  ncm_rng_set_algo (rng, "ranlxd2");
  ncm_rng_set_seed (fresh, 1234);

  g_assert_cmpstr (ncm_rng_get_algo (rng), ==, "ranlxd2");
  g_assert_cmpuint (ncm_rng_get_seed (rng), ==, 1234);

  for (i = 0; i < 10; i++)
    g_assert_cmpfloat (ncm_rng_uniform01_gen (rng), ==, ncm_rng_uniform01_gen (fresh));

  ncm_rng_set_algo (rng, NULL);
  g_assert_cmpstr (ncm_rng_get_algo (rng), ==, gsl_rng_default->name);
  g_assert_cmpuint (ncm_rng_get_seed (rng), ==, 1234);

  ncm_rng_free (rng);
  ncm_rng_free (fresh);
}

/* The used-seed table keeps the whole seed */
static void
test_ncm_rng_check_seed_width (void)
{
  NcmRNG *rng       = ncm_rng_new (NULL);
  const gulong seed = 987654321UL;

  ncm_rng_set_seed (rng, seed);
  g_assert_false (ncm_rng_check_seed (rng, seed));

  if (sizeof (gulong) > 4)
    g_assert_true (ncm_rng_check_seed (rng, seed + (((gulong) 1) << 32)));

  ncm_rng_free (rng);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/rng/seed_state", &test_ncm_rng_seed_state);
  g_test_add_func ("/ncm/rng/set_algo", &test_ncm_rng_set_algo);
  g_test_add_func ("/ncm/rng/check_seed_width", &test_ncm_rng_check_seed_width);
  g_test_add_func ("/ncm/rng/pool", &test_ncm_rng_pool);
  g_test_add_func ("/ncm/rng/distributions", &test_ncm_rng_distributions);
  g_test_add_func ("/ncm/rng/uniform_ranges", &test_ncm_rng_uniform_ranges);
  g_test_add_func ("/ncm/rng/discrete", &test_ncm_rng_discrete);
  g_test_add_func ("/ncm/rng/sample_choose", &test_ncm_rng_sample_choose);
  g_test_add_func ("/ncm/rng/multivariate", &test_ncm_rng_multivariate);
  g_test_add_func ("/ncm/rng/traps", &test_ncm_rng_traps);
  g_test_add_func ("/ncm/rng/invalid_algo/subprocess", &test_ncm_rng_invalid_algo_subprocess);

  g_test_run ();
}

