/***************************************************************************
 *            test_ncm_spectral.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spectral.c
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gegenbauer.h>

#define _assert_small(err, tol) g_assert_cmpfloat ((err), <, (tol))

/* T_k at t in [-1, 1] */
static gdouble
_T (guint k, gdouble t)
{
  return cos (k * acos (t));
}

/* T_k' = k U_{k-1}, with U from its recurrence */
static gdouble
_dT (guint k, gdouble t)
{
  gdouble Um1 = 0.0, U0 = 1.0;
  guint i;

  if (k == 0)
    return 0.0;

  for (i = 1; i < k; i++)
  {
    const gdouble Up1 = 2.0 * t * U0 - Um1;

    Um1 = U0;
    U0  = Up1;
  }

  return k * U0;
}

/* T_k'' = 2 k C^(2)_{k-2} */
static gdouble
_d2T (guint k, gdouble t)
{
  if (k < 2)
    return 0.0;

  return 2.0 * k * gsl_sf_gegenpoly_n (k - 2, 2.0, t);
}

/* T_k(s) for any s, by the three-term recurrence */
static gdouble
_T_any (guint k, gdouble s)
{
  gdouble Tm1 = 1.0, T0 = s;
  guint i;

  if (k == 0)
    return 1.0;

  for (i = 1; i < k; i++)
  {
    const gdouble Tp1 = 2.0 * s * T0 - Tm1;

    Tm1 = T0;
    T0  = Tp1;
  }

  return T0;
}

static GArray *
_array_new (guint n)
{
  GArray *a = g_array_sized_new (FALSE, TRUE, sizeof (gdouble), n);

  g_array_set_size (a, n);

  return a;
}

/* Decaying coefficients with no special structure */
static GArray *
_pattern_new (guint n)
{
  GArray *a = _array_new (n);
  guint k;

  for (k = 0; k < n; k++)
    g_array_index (a, gdouble, k) = cos (2.1 * k + 0.3) / (1.0 + k);

  return a;
}

typedef gdouble (*_Basis) (guint k, gdouble t);

static gdouble
_sum (GArray *c, _Basis basis, gdouble t)
{
  gdouble s = 0.0;
  guint k;

  for (k = 0; k < c->len; k++)
    s += g_array_index (c, gdouble, k) * basis (k, t);

  return s;
}

static gdouble
_gegen_sum (GArray *c, gdouble lambda, gdouble t)
{
  gdouble s = 0.0;
  guint k;

  for (k = 0; k < c->len; k++)
    s += g_array_index (c, gdouble, k) * gsl_sf_gegenpoly_n (k, lambda, t);

  return s;
}

static const gdouble test_t[] = {-1.0, -1.0 + 1.0e-10, -0.95, -0.5, 0.0, 0.3, 0.9, 0.97, 1.0 - 1.0e-10, 1.0};

typedef struct _Counted
{
  guint calls;
} Counted;

static gdouble
_gauss (gpointer user_data, gdouble x)
{
  Counted *c = user_data;

  if (c != NULL)
    c->calls++;

  return exp (-(x - 2.0) * (x - 2.0));
}

static gdouble
_exp (gpointer user_data, gdouble x)
{
  return exp (x);
}

static gdouble
_tiny (gpointer user_data, gdouble x)
{
  return 1.0e-8 * exp (-x * x);
}

static gdouble
_sin_hard (gpointer user_data, gdouble x)
{
  return sin (40.0 * x);
}

static gdouble
_cheb_mode (gpointer user_data, gdouble x)
{
  const guint m = GPOINTER_TO_UINT (user_data);

  return _T (m, CLAMP (ncm_spectral_x_to_s (0.5, 3.0, x), -1.0, 1.0));
}

static void
test_ncm_spectral_max_level (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectral *dup;
  guint max_level;

  g_assert_cmpuint (ncm_spectral_get_max_level (spectral), ==, 16);

  ncm_spectral_set_max_level (spectral, 7);
  g_object_get (spectral, "max-level", &max_level, NULL);
  g_assert_cmpuint (max_level, ==, 7);

  dup = NCM_SPECTRAL (ncm_serialize_dup_obj (ser, G_OBJECT (spectral)));
  g_assert_cmpuint (ncm_spectral_get_max_level (dup), ==, 7);

  ncm_spectral_clear (&dup);
  ncm_serialize_free (ser);
  ncm_spectral_free (spectral);

  spectral = ncm_spectral_new_with_max_level (9);
  g_assert_cmpuint (ncm_spectral_get_max_level (spectral), ==, 9);
  g_assert_true (ncm_spectral_ref (spectral) == spectral);
  ncm_spectral_free (spectral);
  ncm_spectral_clear (&spectral);
  g_assert_null (spectral);
}

static void
test_ncm_spectral_x_t (void)
{
  const gdouble a = -0.7, b = 4.2;
  gdouble x;

  _assert_small (fabs (ncm_spectral_x_to_s (a, b, a) + 1.0), 1.0e-15);
  _assert_small (fabs (ncm_spectral_x_to_s (a, b, b) - 1.0), 1.0e-15);
  _assert_small (fabs (ncm_spectral_s_to_x (a, b, -1.0) - a), 1.0e-15);
  _assert_small (fabs (ncm_spectral_s_to_x (a, b, 1.0) - b), 1.0e-15);

  for (x = a; x <= b; x += 0.37)
    _assert_small (fabs (ncm_spectral_s_to_x (a, b, ncm_spectral_x_to_s (a, b, x)) - x), 1.0e-15);
}

/* On N Lobatto nodes, T_m for m < N is reproduced exactly: its coefficients are delta_{mk} */
static void
test_ncm_spectral_fixed_modes (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *coeffs        = NULL;
  const guint N         = 17;
  guint m, k;

  for (m = 0; m < N; m++)
  {
    gdouble err = 0.0;

    ncm_spectral_compute_chebyshev_coeffs (spectral, _cheb_mode, 0.5, 3.0, N, &coeffs, GUINT_TO_POINTER (m));
    g_assert_cmpuint (coeffs->len, ==, N);

    for (k = 0; k < N; k++)
      err = MAX (err, fabs (g_array_index (coeffs, gdouble, k) - (k == m)));

    _assert_small (err, 1.0e-14);
  }

  g_array_unref (coeffs);
  ncm_spectral_free (spectral);
}

/* e^t = I_0(1) + 2 sum_k I_k(1) T_k(t); the order changes between calls */
static void
test_ncm_spectral_fixed_exp (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *coeffs        = NULL;
  const guint Ns[]      = {24, 9, 24};
  guint i, k;

  for (i = 0; i < G_N_ELEMENTS (Ns); i++)
  {
    const guint N = Ns[i];
    gdouble err   = 0.0;

    ncm_spectral_compute_chebyshev_coeffs (spectral, _exp, -1.0, 1.0, N, &coeffs, NULL);

    for (k = 0; k < N - 1; k++)
      err = MAX (err, fabs (g_array_index (coeffs, gdouble, k) - (k == 0 ? 1.0 : 2.0) * gsl_sf_bessel_In (k, 1.0)));

    _assert_small (err, (N == 24) ? 3.0e-15 : 1.0e-7);
  }

  g_array_unref (coeffs);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_adaptive (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  Counted counted       = {0};
  GArray *coeffs        = NULL;
  GArray *fixed         = NULL;
  gdouble err           = 0.0;
  guint k, N, i;

  k = ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, _gauss, 0.0, 4.0, 3, 1.0e-10, &coeffs, &counted);
  N = (1 << k) + 1;

  g_assert_cmpuint (coeffs->len, ==, N);
  g_assert_cmpuint (counted.calls, ==, N);

  for (i = 0; i <= 100; i++)
  {
    const gdouble x = 4.0 * i / 100.0;

    err = MAX (err, fabs (ncm_spectral_chebyshev_eval_x (coeffs, 0.0, 4.0, x) - _gauss (NULL, x)));
  }

  _assert_small (err, 3.0e-15);

  ncm_spectral_compute_chebyshev_coeffs (spectral, _gauss, 0.0, 4.0, N, &fixed, NULL);
  err = 0.0;

  for (i = 0; i < N; i++)
    err = MAX (err, fabs (g_array_index (coeffs, gdouble, i) - g_array_index (fixed, gdouble, i)));

  _assert_small (err, 1.0e-16);

  g_array_unref (coeffs);
  g_array_unref (fixed);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_adaptive_full (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *c_rel         = NULL;
  GArray *c_full        = NULL;
  guint k_rel, k_full, i;

  k_rel  = ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, _tiny, -3.0, 3.0, 2, 1.0e-12, &c_rel, NULL);
  k_full = ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, _tiny, -3.0, 3.0, 2, 1.0e-12, 0.0, &c_full, NULL);

  g_assert_cmpuint (k_rel, ==, k_full);

  for (i = 0; i < c_rel->len; i++)
    g_assert_cmpfloat (g_array_index (c_rel, gdouble, i), ==, g_array_index (c_full, gdouble, i));

  k_full = ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, _tiny, -3.0, 3.0, 2, 1.0e-12, 1.0e-6, &c_full, NULL);
  g_assert_cmpuint (k_full, ==, 3);
  g_assert_cmpuint (k_full, <, k_rel);

  g_array_unref (c_rel);
  g_array_unref (c_full);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_adaptive_try (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *c_try         = NULL;
  GArray *c_full        = NULL;
  gboolean converged    = TRUE;
  guint k, i;

  k = ncm_spectral_compute_chebyshev_coeffs_adaptive_try (spectral, _sin_hard, -1.0, 1.0, 2, 4, 1.0e-10, 0.0, &c_try, NULL, &converged);
  g_assert_false (converged);
  g_assert_cmpuint (k, ==, 4);
  g_assert_cmpuint (c_try->len, ==, 17);

  k = ncm_spectral_compute_chebyshev_coeffs_adaptive_try (spectral, _exp, -1.0, 1.0, 2, 30, 1.0e-12, 0.0, &c_try, NULL, &converged);
  g_assert_true (converged);
  g_assert_cmpuint (k, ==, ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, _exp, -1.0, 1.0, 2, 1.0e-12, 0.0, &c_full, NULL));

  for (i = 0; i < c_try->len; i++)
    g_assert_cmpfloat (g_array_index (c_try, gdouble, i), ==, g_array_index (c_full, gdouble, i));

  g_array_unref (c_try);
  g_array_unref (c_full);
  ncm_spectral_free (spectral);
}

/* With level_min equal to the cap there is no doubling: the level_min coefficients come back,
 * also into a reused array, and the tolerance is reported as not met */
static void
test_ncm_spectral_adaptive_single_level (void)
{
  NcmSpectral *spectral = ncm_spectral_new_with_max_level (4);
  GArray *coeffs        = NULL;
  GArray *fixed         = NULL;
  gboolean converged    = TRUE;
  guint pass, k, i;

  ncm_spectral_compute_chebyshev_coeffs (spectral, _exp, -1.0, 1.0, 17, &fixed, NULL);

  for (pass = 0; pass < 2; pass++)
  {
    if (pass == 1)
    {
      g_array_set_size (coeffs, 3);
      g_array_index (coeffs, gdouble, 0) = 7.0;
    }

    k = ncm_spectral_compute_chebyshev_coeffs_adaptive_try (spectral, _exp, -1.0, 1.0, 4, 4, 1.0e-12, 0.0, &coeffs, NULL, &converged);

    g_assert_cmpuint (k, ==, 4);
    g_assert_false (converged);
    g_assert_cmpuint (coeffs->len, ==, 17);

    for (i = 0; i < 17; i++)
      g_assert_cmpfloat (g_array_index (coeffs, gdouble, i), ==, g_array_index (fixed, gdouble, i));
  }

  g_array_unref (coeffs);
  g_array_unref (fixed);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_invalid_arguments (void)
{
  g_test_trap_subprocess ("/ncm/spectral/invalid/max_level/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*max_level <= *");

  g_test_trap_subprocess ("/ncm/spectral/invalid/order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*N >= 2*");
}

static void
test_ncm_spectral_invalid_max_level_subprocess (void)
{
  NcmSpectral *spectral = ncm_spectral_new_with_max_level (4);

  ncm_spectral_set_max_level (spectral, 31);
}

static void
test_ncm_spectral_invalid_order_subprocess (void)
{
  NcmSpectral *spectral = ncm_spectral_new_with_max_level (4);
  GArray *coeffs        = NULL;

  ncm_spectral_compute_chebyshev_coeffs (spectral, _exp, -1.0, 1.0, 1, &coeffs, NULL);
}

static void
test_ncm_spectral_adaptive_fatal (void)
{
  g_test_trap_subprocess ("/ncm/spectral/adaptive/fatal/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*without converging*");
}

static void
test_ncm_spectral_adaptive_fatal_subprocess (void)
{
  NcmSpectral *spectral = ncm_spectral_new_with_max_level (5);
  GArray *coeffs        = NULL;

  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, _sin_hard, -1.0, 1.0, 2, 1.0e-12, &coeffs, NULL);
}

/* Batch: exp(x), cos(5x), 1/(1+x^2) on [-1, 2] */
static void
_batch3 (gpointer user_data, gdouble x, NcmVector *y)
{
  Counted *c = user_data;

  c->calls++;
  ncm_vector_set (y, 0, exp (x));
  ncm_vector_set (y, 1, cos (5.0 * x));

  if (ncm_vector_len (y) > 2)
    ncm_vector_set (y, 2, 1.0 / (1.0 + x * x));
}

static gdouble
_batch3_comp (gpointer user_data, gdouble x)
{
  switch (GPOINTER_TO_UINT (user_data))
  {
    case 0:
      return exp (x);

    case 1:
      return cos (5.0 * x);

    default:
      return 1.0 / (1.0 + x * x);
  }
}

static void
_batch_hard (gpointer user_data, gdouble x, NcmVector *y)
{
  ncm_vector_set (y, 0, exp (x));
  ncm_vector_set (y, 1, sin (60.0 * x));
}

static void
test_ncm_spectral_batch (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  NcmMatrix *coeffs     = NULL;
  GArray *fixed         = NULL;
  const guint n_comps[] = {3, 2};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (n_comps); i++)
  {
    const guint n_comp = n_comps[i];
    Counted counted    = {0};
    guint k, N, c, j;

    k = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (spectral, _batch3, n_comp, -1.0, 2.0, 3, 1.0e-12, 0.0, &coeffs, &counted);
    N = (1 << k) + 1;

    g_assert_cmpuint (counted.calls, ==, N);
    g_assert_cmpuint (ncm_matrix_nrows (coeffs), ==, n_comp);
    g_assert_cmpuint (ncm_matrix_ncols (coeffs), ==, N);

    for (c = 0; c < n_comp; c++)
    {
      gdouble err = 0.0;

      ncm_spectral_compute_chebyshev_coeffs (spectral, _batch3_comp, -1.0, 2.0, N, &fixed, GUINT_TO_POINTER (c));

      for (j = 0; j < N; j++)
        err = MAX (err, fabs (ncm_matrix_get (coeffs, c, j) - g_array_index (fixed, gdouble, j)));

      _assert_small (err, 1.0e-15);
    }
  }

  ncm_matrix_free (coeffs);
  g_array_unref (fixed);
  ncm_spectral_free (spectral);
}

/* 1 + 2x, a Gaussian and a component six orders below them */
static void
_batch_scales (gpointer user_data, gdouble x, NcmVector *y)
{
  ncm_vector_set (y, 0, 1.0 + 2.0 * x);
  ncm_vector_set (y, 1, exp (-(x - 2.0) * (x - 2.0)));
  ncm_vector_set (y, 2, 1.0e-6 * cos (30.0 * x));
}

/* Convergence is per component: the small one is resolved relative to its own size */
static void
test_ncm_spectral_batch_small_component (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  NcmMatrix *coeffs     = NULL;
  NcmVector *y          = ncm_vector_new (3);
  const gdouble a       = 0.3, b = 4.7;
  gdouble err[3] = {0.0, 0.0, 0.0};
  gdouble high   = 0.0;
  guint k, N, c, i;

  k = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (spectral, _batch_scales, 3, a, b, 3, 1.0e-12, 0.0, &coeffs, NULL);
  N = (1 << k) + 1;

  for (i = 0; i <= 200; i++)
  {
    const gdouble x = a + (b - a) * i / 200.0;
    GArray *row     = _array_new (N);

    _batch_scales (NULL, x, y);

    for (c = 0; c < 3; c++)
    {
      memcpy (row->data, ncm_matrix_ptr (coeffs, c, 0), N * sizeof (gdouble));
      err[c] = MAX (err[c], fabs (ncm_spectral_chebyshev_eval_x (row, a, b, x) - ncm_vector_get (y, c)));
    }

    g_array_unref (row);
  }

  for (i = 2; i < N; i++)
    high = MAX (high, fabs (ncm_matrix_get (coeffs, 0, i)));

  _assert_small (err[0], 2.0e-14);
  _assert_small (err[1], 5.0e-15);
  _assert_small (err[2] / 1.0e-6, 2.0e-13);
  _assert_small (high, 5.0e-15);

  ncm_vector_free (y);
  ncm_matrix_free (coeffs);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_batch_cap (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  NcmMatrix *coeffs     = ncm_matrix_new (2, 3);
  NcmMatrix *orig       = coeffs;
  guint k;

  ncm_matrix_set_all (coeffs, 7.0);

  k = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (spectral, _batch_hard, 2, -1.0, 1.0, 2, 4, 1.0e-10, 0.0, FALSE, &coeffs, NULL);

  g_assert_cmpuint (k, ==, 0);
  g_assert_true (coeffs == orig);
  g_assert_cmpfloat (ncm_matrix_get (coeffs, 1, 2), ==, 7.0);

  ncm_matrix_free (coeffs);
  ncm_spectral_free (spectral);
}

/* Below four intervals the coefficients have no two bands to extrapolate from, so the
 * last doubling the cap allows is always tried: every node of the next level is
 * evaluated, even though it cannot converge at these tolerances. */
static void
test_ncm_spectral_batch_no_prediction_below_four (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  NcmMatrix *coeffs     = NULL;
  guint level_min;

  for (level_min = 0; level_min <= 1; level_min++)
  {
    Counted calls = {0};
    guint k;

    k = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (spectral, _batch3, 3, -1.0, 2.0, level_min, level_min + 1, 1.0e-10, 0.0, FALSE, &coeffs, &calls);

    g_assert_cmpuint (k, ==, 0);
    g_assert_cmpuint (calls.calls, ==, (1u << (level_min + 1)) + 1u);
  }

  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_batch_fatal (void)
{
  g_test_trap_subprocess ("/ncm/spectral/batch/fatal/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*without converging*");
}

static void
test_ncm_spectral_batch_fatal_subprocess (void)
{
  NcmSpectral *spectral = ncm_spectral_new_with_max_level (5);
  NcmMatrix *coeffs     = NULL;

  ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (spectral, _batch_hard, 2, -1.0, 1.0, 2, 1.0e-12, 0.0, &coeffs, NULL);
}

/* Every expansion path gives the same coefficients after the buffers, plans and node
 * tables are released and recreated. */
static void
test_ncm_spectral_free_buffers (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *fixed_a       = NULL;
  GArray *fixed_b       = NULL;
  GArray *adapt_a       = NULL;
  GArray *adapt_b       = NULL;
  NcmMatrix *batch_a    = NULL;
  NcmMatrix *batch_b    = NULL;
  Counted counted_a     = {0};
  Counted counted_b     = {0};
  guint k_a, k_b, i, c;

  ncm_spectral_compute_chebyshev_coeffs (spectral, _exp, 0.5, 2.0, 24, &fixed_a, NULL);
  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, _exp, -1.0, 1.0, 2, 1.0e-12, &adapt_a, NULL);
  k_a = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (spectral, _batch3, 3, -1.0, 2.0, 3, 1.0e-12, 0.0, &batch_a, &counted_a);

  ncm_spectral_free_buffers (spectral);

  ncm_spectral_compute_chebyshev_coeffs (spectral, _exp, 0.5, 2.0, 24, &fixed_b, NULL);
  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, _exp, -1.0, 1.0, 2, 1.0e-12, &adapt_b, NULL);
  k_b = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (spectral, _batch3, 3, -1.0, 2.0, 3, 1.0e-12, 0.0, &batch_b, &counted_b);

  g_assert_cmpuint (fixed_a->len, ==, fixed_b->len);
  g_assert_cmpuint (adapt_a->len, ==, adapt_b->len);
  g_assert_cmpuint (k_a, ==, k_b);

  for (i = 0; i < fixed_a->len; i++)
    g_assert_cmpfloat (g_array_index (fixed_a, gdouble, i), ==, g_array_index (fixed_b, gdouble, i));

  for (i = 0; i < adapt_a->len; i++)
    g_assert_cmpfloat (g_array_index (adapt_a, gdouble, i), ==, g_array_index (adapt_b, gdouble, i));

  for (c = 0; c < 3; c++)
    for (i = 0; i < ncm_matrix_ncols (batch_a); i++)
      g_assert_cmpfloat (ncm_matrix_get (batch_a, c, i), ==, ncm_matrix_get (batch_b, c, i));

  /* Releasing twice, and releasing an instance that never expanded, is allowed */
  ncm_spectral_free_buffers (spectral);
  ncm_spectral_free_buffers (spectral);

  g_array_unref (fixed_a);
  g_array_unref (fixed_b);
  g_array_unref (adapt_a);
  g_array_unref (adapt_b);
  ncm_matrix_free (batch_a);
  ncm_matrix_free (batch_b);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_eval_deriv (void)
{
  GArray *c = _pattern_new (30);
  const gdouble a   = -2.0, b = 1.5;
  gdouble err_f     = 0.0, err_df = 0.0, err_x = 0.0, scale_x = 0.0;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_t); i++)
  {
    const gdouble t     = test_t[i];
    const gdouble x     = ncm_spectral_s_to_x (a, b, t);
    const gdouble f_x   = ncm_spectral_chebyshev_eval (c, ncm_spectral_x_to_s (a, b, x));
    const gdouble df_dx = ncm_spectral_chebyshev_deriv (c, ncm_spectral_x_to_s (a, b, x)) * 2.0 / (b - a);

    err_f   = MAX (err_f, fabs (ncm_spectral_chebyshev_eval (c, t) - _sum (c, _T, t)));
    err_df  = MAX (err_df, fabs (ncm_spectral_chebyshev_deriv (c, t) - _sum (c, _dT, t)));
    err_x   = MAX (err_x, fabs (ncm_spectral_chebyshev_eval_x (c, a, b, x) - f_x));
    err_x   = MAX (err_x, fabs (ncm_spectral_chebyshev_deriv_x (c, a, b, x) - df_dx));
    scale_x = MAX (scale_x, MAX (fabs (f_x), fabs (df_dx)));
  }

  _assert_small (err_f, 1.0e-14);
  _assert_small (err_df, 5.0e-13);

  /* The same computation on both sides, so the two agree to rounding; exactly on x86-64,
   * within 1.02e-15 on macOS arm64, where the compiler may fuse x_to_s differently. */
  _assert_small (err_x, 8.0 * GSL_DBL_EPSILON * scale_x);

  g_array_set_size (c, 1);
  g_assert_cmpfloat (ncm_spectral_chebyshev_eval (c, 0.3), ==, g_array_index (c, gdouble, 0));
  g_assert_cmpfloat (ncm_spectral_chebyshev_deriv (c, 0.3), ==, 0.0);
  g_array_set_size (c, 0);
  g_assert_cmpfloat (ncm_spectral_chebyshev_eval (c, 0.3), ==, 0.0);

  g_array_unref (c);
}

static void
test_ncm_spectral_integrate (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *coeffs        = NULL;
  GArray *mode          = _array_new (12);
  guint k;

  ncm_spectral_compute_chebyshev_coeffs (spectral, _exp, 0.5, 2.0, 24, &coeffs, NULL);
  _assert_small (fabs (ncm_spectral_chebyshev_integrate (coeffs, 0.5, 2.0) / (exp (2.0) - exp (0.5)) - 1.0), 1.0e-15);

  for (k = 0; k < mode->len; k++)
  {
    g_array_index (mode, gdouble, k) = 1.0;
    _assert_small (fabs (ncm_spectral_chebyshev_integrate (mode, -1.0, 1.0) - ((k % 2) ? 0.0 : 2.0 / (1.0 - k * k))), 1.0e-16);
    g_array_index (mode, gdouble, k) = 0.0;
  }

  g_array_unref (coeffs);
  g_array_unref (mode);
  ncm_spectral_free (spectral);
}

static void
test_ncm_spectral_gegenbauer_eval (void)
{
  GArray *c = _pattern_new (12);
  gdouble err1 = 0.0, err2 = 0.0, err_x = 0.0;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_t); i++)
  {
    const gdouble t = test_t[i];
    const gdouble x = ncm_spectral_s_to_x (1.0, 3.0, t);

    err1  = MAX (err1, fabs (ncm_spectral_gegenbauer_alpha1_eval (c, t) - _gegen_sum (c, 1.0, t)));
    err2  = MAX (err2, fabs (ncm_spectral_gegenbauer_alpha2_eval (c, t) - _gegen_sum (c, 2.0, t)));
    err_x = MAX (err_x, fabs (ncm_spectral_gegenbauer_alpha1_eval_x (c, 1.0, 3.0, x) - ncm_spectral_gegenbauer_alpha1_eval (c, ncm_spectral_x_to_s (1.0, 3.0, x))));
    err_x = MAX (err_x, fabs (ncm_spectral_gegenbauer_alpha2_eval_x (c, 1.0, 3.0, x) - ncm_spectral_gegenbauer_alpha2_eval (c, ncm_spectral_x_to_s (1.0, 3.0, x))));
  }

  _assert_small (err1, 1.0e-14);
  _assert_small (err2, 3.0e-14);
  _assert_small (err_x, 1.0e-15);

  g_array_unref (c);
}

/* Conversions: the C^(1) and C^(2) series of c, and of its derivatives, evaluated against the T_k sums */
static void
test_ncm_spectral_gegenbauer_conversions (void)
{
  GArray *c = _pattern_new (20);
  GArray *g1  = NULL, *g2 = NULL, *gd = NULL, *gd2 = NULL, *gx = NULL;
  gdouble err = 0.0, err_d = 0.0, err_d2 = 0.0, err_x = 0.0;
  guint i;

  ncm_spectral_chebT_to_gegenbauer_alpha1 (c, &g1);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (c, &g2);
  ncm_spectral_chebT_deriv_to_gegenbauer_alpha2 (c, &gd);
  ncm_spectral_chebT_deriv2_to_gegenbauer_alpha2 (c, &gd2);
  ncm_spectral_gegenbauer_alpha2_mul_affine (g2, 0.7, -0.4, &gx);

  g_assert_cmpuint (g1->len, ==, 20);
  g_assert_cmpuint (g2->len, ==, 20);
  g_assert_cmpuint (gd->len, ==, 19);
  g_assert_cmpuint (gd2->len, ==, 18);
  g_assert_cmpuint (gx->len, ==, 21);

  for (i = 0; i < G_N_ELEMENTS (test_t); i++)
  {
    const gdouble t = test_t[i];
    const gdouble f = _sum (c, _T, t);

    err    = MAX (err, fabs (ncm_spectral_gegenbauer_alpha1_eval (g1, t) - f));
    err    = MAX (err, fabs (ncm_spectral_gegenbauer_alpha2_eval (g2, t) - f));
    err_d  = MAX (err_d, fabs (ncm_spectral_gegenbauer_alpha2_eval (gd, t) - _sum (c, _dT, t)));
    err_x  = MAX (err_x, fabs (ncm_spectral_gegenbauer_alpha2_eval (gx, t) - (0.7 * t - 0.4) * f));
    err_d2 = MAX (err_d2, fabs (ncm_spectral_gegenbauer_alpha2_eval (gd2, t) - _sum (c, _d2T, t)));
  }

  _assert_small (err, 1.0e-14);
  _assert_small (err_d, 5.0e-13);
  _assert_small (err_d2, 3.0e-12);
  _assert_small (err_x, 5.0e-15);

  g_array_set_size (c, 1);
  ncm_spectral_chebT_deriv_to_gegenbauer_alpha2 (c, &gd);
  ncm_spectral_chebT_deriv2_to_gegenbauer_alpha2 (c, &gd2);
  g_assert_cmpuint (gd->len, ==, 1);
  g_assert_cmpuint (gd2->len, ==, 1);
  g_assert_cmpfloat (g_array_index (gd, gdouble, 0), ==, 0.0);
  g_assert_cmpfloat (g_array_index (gd2, gdouble, 0), ==, 0.0);

  g_array_unref (c);
  g_array_unref (g1);
  g_array_unref (g2);
  g_array_unref (gd);
  g_array_unref (gd2);
  g_array_unref (gx);
}

typedef NcmMatrix *(*_OpMatrix) (guint N);
typedef void (*_OpRow) (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);

static gdouble
_op_f (guint op, GArray *c, gdouble t)
{
  switch (op)
  {
    case 0:
      return _sum (c, _T, t);

    case 1:
      return t * _sum (c, _T, t);

    case 2:
      return t * t * _sum (c, _T, t);

    case 3:
      return _sum (c, _dT, t);

    case 4:
      return t * _sum (c, _dT, t);

    case 5:
      return _sum (c, _d2T, t);

    case 6:
      return t * _sum (c, _d2T, t);

    default:
      return t * t * _sum (c, _d2T, t);
  }
}

static void
_d_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  ncm_spectral_compute_d_row (row_data, offset, coeff);
}

/* Each operator matrix applied to c, evaluated as a C^(2) series, against the operator
 * applied to the T_k sum. c ends 8 below N, so the truncated columns hold nothing. */
static void
test_ncm_spectral_operators (void)
{
  const _OpMatrix mats[] = {
    ncm_spectral_get_proj_matrix, ncm_spectral_get_s_matrix, ncm_spectral_get_s2_matrix,
    ncm_spectral_get_d_matrix, ncm_spectral_get_s_d_matrix, ncm_spectral_get_d2_matrix,
    ncm_spectral_get_s_d2_matrix, ncm_spectral_get_s2_d2_matrix
  };
  const _OpRow rows[] = {
    ncm_spectral_compute_proj_row, ncm_spectral_compute_s_row, ncm_spectral_compute_s2_row,
    _d_row, ncm_spectral_compute_s_d_row, ncm_spectral_compute_d2_row,
    ncm_spectral_compute_s_d2_row, ncm_spectral_compute_s2_d2_row
  };
  const gdouble tols[] = {1.0e-14, 1.0e-14, 1.0e-14, 2.0e-13, 2.0e-13, 3.0e-12, 3.0e-12, 3.0e-12};
  const guint N        = 24;
  GArray *c            = _pattern_new (N - 8);
  GArray *g            = _array_new (N);
  guint op, i, j;

  for (op = 0; op < G_N_ELEMENTS (mats); op++)
  {
    NcmMatrix *M = mats[op](N);
    gdouble err  = 0.0;

    for (i = 0; i < N; i++)
    {
      gdouble gi = 0.0;

      for (j = 0; j < c->len; j++)
        gi += ncm_matrix_get (M, i, j) * g_array_index (c, gdouble, j);

      g_array_index (g, gdouble, i) = gi;
    }

    for (i = 0; i < G_N_ELEMENTS (test_t); i++)
    {
      const gdouble t = test_t[i];

      err = MAX (err, fabs (ncm_spectral_gegenbauer_alpha2_eval (g, t) - _op_f (op, c, t)));
    }

    _assert_small (err, tols[op]);

    /* A row is linear in its coefficient and adds to what is there */
    for (i = 0; i < 4; i++)
    {
      gdouble r1[12] = {0.0}, r2[12] = {0.0};

      rows[op](r1, i, 2, 1.0);
      rows[op](r2, i, 2, 1.5);
      rows[op](r2, i, 2, 1.0);

      for (j = 0; j < 12; j++)
        _assert_small (fabs (r2[j] - 2.5 * r1[j]), 1.0e-15);
    }

    ncm_matrix_free (M);
  }

  g_array_unref (c);
  g_array_unref (g);
}

static void
test_ncm_spectral_rebase (void)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *c             = _array_new (3);
  GArray *r             = NULL;
  GArray *r5            = NULL;
  GArray *c5;
  gdouble norm, err;
  guint i;

  /* T_2(s) on [0, 2] is (T_2(t) + 1)/4 - T_1(t) - 1/2 on [0, 1] */
  g_array_index (c, gdouble, 2) = 1.0;
  norm                          = ncm_spectral_chebyshev_rebase (spectral, c, 0, 0.0, 2.0, 0.0, 1.0, &r);
  g_assert_cmpuint (r->len, ==, 3);
  _assert_small (fabs (g_array_index (r, gdouble, 0) + 0.25), 1.0e-15);
  _assert_small (fabs (g_array_index (r, gdouble, 1) + 1.0), 1.0e-15);
  _assert_small (fabs (g_array_index (r, gdouble, 2) - 0.25), 1.0e-15);
  _assert_small (fabs (norm - 1.5), 1.0e-15);
  g_array_unref (c);

  /* The same interval gives the same series; another gives the same polynomial */
  c = _pattern_new (14);
  ncm_spectral_chebyshev_rebase (spectral, c, 0, -1.0, 3.0, -1.0, 3.0, &r);
  err = 0.0;

  for (i = 0; i < c->len; i++)
    err = MAX (err, fabs (g_array_index (r, gdouble, i) - g_array_index (c, gdouble, i)));

  _assert_small (err, 1.0e-15);

  norm = ncm_spectral_chebyshev_rebase (spectral, c, 0, -1.0, 3.0, 0.5, 5.0, &r);
  err  = 0.0;

  {
    gdouble fmax = 0.0, sum_abs = 0.0;

    for (i = 0; i <= 50; i++)
    {
      const gdouble x = 0.5 + 4.5 * i / 50.0;
      const gdouble f = _sum (c, _T_any, ncm_spectral_x_to_s (-1.0, 3.0, x));

      err  = MAX (err, fabs (_sum (r, _T_any, ncm_spectral_x_to_s (0.5, 5.0, x)) - f));
      fmax = MAX (fmax, fabs (f));
    }

    for (i = 0; i < r->len; i++)
      sum_abs += fabs (g_array_index (r, gdouble, i));

    _assert_small (err / fmax, 2.0e-15);
    g_assert_cmpfloat (fmax, <=, norm);
    g_assert_cmpfloat (norm, ==, sum_abs);
  }

  /* Extending past the source interval grows the norm */
  {
    const gdouble norm_in = ncm_spectral_chebyshev_rebase (spectral, c, 0, -1.0, 3.0, -1.0, 3.0, &r);

    g_assert_cmpfloat (ncm_spectral_chebyshev_rebase (spectral, c, 0, -1.0, 3.0, -1.0, 30.0, &r), >, 100.0 * norm_in);
  }

  /* The rows variant equals the scalar rebase on every row */
  {
    const guint n_rows = 4;
    const guint ncol   = c->len;
    NcmMatrix *cm      = ncm_matrix_new (n_rows, ncol);
    NcmMatrix *rm      = ncm_matrix_new (n_rows, ncol);
    gdouble norm_rows, norm_max = 0.0;
    guint row;

    for (row = 0; row < n_rows; row++)
      for (i = 0; i < ncol; i++)
        ncm_matrix_set (cm, row, i, g_array_index (c, gdouble, i) * (row + 1.0) * ((i % 3 == row % 3) ? -1.0 : 1.0));

    norm_rows = ncm_spectral_chebyshev_rebase_rows (spectral, cm, -1.0, 3.0, 0.5, 5.0, rm);

    for (row = 0; row < n_rows; row++)
    {
      GArray *crow = _array_new (ncol);

      for (i = 0; i < ncol; i++)
        g_array_index (crow, gdouble, i) = ncm_matrix_get (cm, row, i);

      norm_max = MAX (norm_max, ncm_spectral_chebyshev_rebase (spectral, crow, 0, -1.0, 3.0, 0.5, 5.0, &r));

      for (i = 0; i < ncol; i++)
        g_assert_cmpfloat (ncm_matrix_get (rm, row, i), ==, g_array_index (r, gdouble, i));

      g_array_unref (crow);
    }

    g_assert_cmpfloat (norm_rows, ==, norm_max);

    ncm_matrix_free (cm);
    ncm_matrix_free (rm);
  }

  /* len keeps the leading coefficients */
  c5 = _array_new (5);
  memcpy (c5->data, c->data, 5 * sizeof (gdouble));
  ncm_spectral_chebyshev_rebase (spectral, c, 5, -1.0, 3.0, 0.5, 5.0, &r);
  ncm_spectral_chebyshev_rebase (spectral, c5, 0, -1.0, 3.0, 0.5, 5.0, &r5);
  g_assert_cmpuint (r->len, ==, 5);

  for (i = 0; i < 5; i++)
    g_assert_cmpfloat (g_array_index (r, gdouble, i), ==, g_array_index (r5, gdouble, i));

  g_array_unref (c);
  g_array_unref (c5);

  /* Overflow is reported as infinity, an empty series as zero */
  c = _array_new (400);

  for (i = 0; i < c->len; i++)
    g_array_index (c, gdouble, i) = 1.0;

  g_assert_true (isinf (ncm_spectral_chebyshev_rebase (spectral, c, 0, -1.0, 1.0, -1.0, 1.0e3, &r)));

  g_array_set_size (c, 0);
  g_assert_cmpfloat (ncm_spectral_chebyshev_rebase (spectral, c, 0, -1.0, 1.0, 0.0, 1.0, &r), ==, 0.0);
  g_assert_cmpuint (r->len, ==, 0);

  g_array_unref (c);
  g_array_unref (r);
  g_array_unref (r5);
  ncm_spectral_free (spectral);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/spectral/max_level", &test_ncm_spectral_max_level);
  g_test_add_func ("/ncm/spectral/x_t", &test_ncm_spectral_x_t);
  g_test_add_func ("/ncm/spectral/fixed/modes", &test_ncm_spectral_fixed_modes);
  g_test_add_func ("/ncm/spectral/fixed/exp", &test_ncm_spectral_fixed_exp);
  g_test_add_func ("/ncm/spectral/adaptive", &test_ncm_spectral_adaptive);
  g_test_add_func ("/ncm/spectral/adaptive/full", &test_ncm_spectral_adaptive_full);
  g_test_add_func ("/ncm/spectral/adaptive/try", &test_ncm_spectral_adaptive_try);
  g_test_add_func ("/ncm/spectral/adaptive/single_level", &test_ncm_spectral_adaptive_single_level);
  g_test_add_func ("/ncm/spectral/adaptive/fatal", &test_ncm_spectral_adaptive_fatal);
  g_test_add_func ("/ncm/spectral/adaptive/fatal/subprocess", &test_ncm_spectral_adaptive_fatal_subprocess);
  g_test_add_func ("/ncm/spectral/batch", &test_ncm_spectral_batch);
  g_test_add_func ("/ncm/spectral/batch/small_component", &test_ncm_spectral_batch_small_component);
  g_test_add_func ("/ncm/spectral/batch/cap", &test_ncm_spectral_batch_cap);
  g_test_add_func ("/ncm/spectral/batch/no_prediction_below_four", &test_ncm_spectral_batch_no_prediction_below_four);
  g_test_add_func ("/ncm/spectral/batch/fatal", &test_ncm_spectral_batch_fatal);
  g_test_add_func ("/ncm/spectral/batch/fatal/subprocess", &test_ncm_spectral_batch_fatal_subprocess);
  g_test_add_func ("/ncm/spectral/free_buffers", &test_ncm_spectral_free_buffers);
  g_test_add_func ("/ncm/spectral/eval_deriv", &test_ncm_spectral_eval_deriv);
  g_test_add_func ("/ncm/spectral/integrate", &test_ncm_spectral_integrate);
  g_test_add_func ("/ncm/spectral/gegenbauer/eval", &test_ncm_spectral_gegenbauer_eval);
  g_test_add_func ("/ncm/spectral/gegenbauer/conversions", &test_ncm_spectral_gegenbauer_conversions);
  g_test_add_func ("/ncm/spectral/operators", &test_ncm_spectral_operators);
  g_test_add_func ("/ncm/spectral/rebase", &test_ncm_spectral_rebase);
  g_test_add_func ("/ncm/spectral/invalid/arguments", &test_ncm_spectral_invalid_arguments);
  g_test_add_func ("/ncm/spectral/invalid/max_level/subprocess", &test_ncm_spectral_invalid_max_level_subprocess);
  g_test_add_func ("/ncm/spectral/invalid/order/subprocess", &test_ncm_spectral_invalid_order_subprocess);

  g_test_run ();
}

