/***************************************************************************
 *            test_nc_xcor_cov_gauss.c
 *
 *  Tue November 03 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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
 * #NcXcorCovGauss on three tracers with analytic spectra C^ab_l = A_ab (10 / l) and constant
 * noise. Every tracer pair appears once in the data vector, one of them as (b, a), with
 * top-hat windows and one spectrum with overlapping windows. The covariance is compared
 * with the formula evaluated term by term in this file.
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

#define TEST_N_TRACERS 3
#define TEST_ELL_MIN 2
#define TEST_ELL_MAX 31
#define TEST_N_ELL (TEST_ELL_MAX - TEST_ELL_MIN + 1)
#define TEST_F_SKY 0.3
#define TEST_N_SPECTRA 6

typedef struct _TestNcXcorCovGauss
{
  NcXcorCovGauss *cov;
  NcmMatrix *cl;
  NcmMatrix *noise;
  NcmMatrix *cov_mat;
  gdouble amp[TEST_N_TRACERS][TEST_N_TRACERS];
  gdouble noise_level[TEST_N_TRACERS];
  guint pair_a[TEST_N_SPECTRA];
  guint pair_b[TEST_N_SPECTRA];
  NcmMatrix *window[TEST_N_SPECTRA];
} TestNcXcorCovGauss;

void test_nc_xcor_cov_gauss_new (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_free (TestNcXcorCovGauss *test, gconstpointer pdata);

void test_nc_xcor_cov_gauss_pair_index (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_formula (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_negative_cross (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_noise_only (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_f_sky (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_positive_definite (TestNcXcorCovGauss *test, gconstpointer pdata);
void test_nc_xcor_cov_gauss_serialize (TestNcXcorCovGauss *test, gconstpointer pdata);

void test_nc_xcor_cov_gauss_invalid (void);
void test_nc_xcor_cov_gauss_invalid_tracer_subprocess (void);
void test_nc_xcor_cov_gauss_invalid_window_subprocess (void);
void test_nc_xcor_cov_gauss_invalid_cl_subprocess (void);
void test_nc_xcor_cov_gauss_invalid_f_sky_subprocess (void);
void test_nc_xcor_cov_gauss_invalid_ell_range_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/xcor/cov_gauss/pair_index", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_pair_index,
              &test_nc_xcor_cov_gauss_free);

  g_test_add ("/nc/xcor/cov_gauss/formula", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_formula,
              &test_nc_xcor_cov_gauss_free);

  g_test_add ("/nc/xcor/cov_gauss/negative_cross", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_negative_cross,
              &test_nc_xcor_cov_gauss_free);

  g_test_add ("/nc/xcor/cov_gauss/noise_only", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_noise_only,
              &test_nc_xcor_cov_gauss_free);

  g_test_add ("/nc/xcor/cov_gauss/f_sky", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_f_sky,
              &test_nc_xcor_cov_gauss_free);

  g_test_add ("/nc/xcor/cov_gauss/positive_definite", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_positive_definite,
              &test_nc_xcor_cov_gauss_free);

  g_test_add ("/nc/xcor/cov_gauss/serialize", TestNcXcorCovGauss, NULL,
              &test_nc_xcor_cov_gauss_new,
              &test_nc_xcor_cov_gauss_serialize,
              &test_nc_xcor_cov_gauss_free);

  g_test_add_func ("/nc/xcor/cov_gauss/invalid", &test_nc_xcor_cov_gauss_invalid);
  g_test_add_func ("/nc/xcor/cov_gauss/invalid/tracer/subprocess", &test_nc_xcor_cov_gauss_invalid_tracer_subprocess);
  g_test_add_func ("/nc/xcor/cov_gauss/invalid/window/subprocess", &test_nc_xcor_cov_gauss_invalid_window_subprocess);
  g_test_add_func ("/nc/xcor/cov_gauss/invalid/cl/subprocess", &test_nc_xcor_cov_gauss_invalid_cl_subprocess);
  g_test_add_func ("/nc/xcor/cov_gauss/invalid/f_sky/subprocess", &test_nc_xcor_cov_gauss_invalid_f_sky_subprocess);
  g_test_add_func ("/nc/xcor/cov_gauss/invalid/ell_range/subprocess", &test_nc_xcor_cov_gauss_invalid_ell_range_subprocess);

  g_test_run ();
}

/* Window of @n_bp top-hat band powers of @width multipoles each, from the first multipole on. */
static NcmMatrix *
_test_nc_xcor_cov_gauss_window_tophat (guint n_bp, guint width)
{
  NcmMatrix *window = ncm_matrix_new (n_bp, TEST_N_ELL);
  guint i, l;

  ncm_matrix_set_zero (window);

  for (i = 0; i < n_bp; i++)
    for (l = i * width; l < (i + 1) * width; l++)
      ncm_matrix_set (window, i, l, 1.0 / width);

  return window;
}

/* Two triangular band powers over multipoles 0-19 and 10-29, overlapping on 10-19. */
static NcmMatrix *
_test_nc_xcor_cov_gauss_window_overlap (void)
{
  NcmMatrix *window = ncm_matrix_new (2, TEST_N_ELL);
  guint i, l;

  ncm_matrix_set_zero (window);

  for (i = 0; i < 2; i++)
    for (l = 0; l < 20; l++)
      ncm_matrix_set (window, i, 10 * i + l, (l < 10) ? (l + 1.0) / 110.0 : (20.0 - l) / 110.0);

  return window;
}

/* Signal plus noise of the tracers (@x, @y) at multipole @ell, from the fixture's amplitudes. */
static gdouble
_test_nc_xcor_cov_gauss_cl_tot (TestNcXcorCovGauss *test, gboolean with_signal, guint x, guint y, guint ell)
{
  const gdouble signal = with_signal ? test->amp[x][y] * 10.0 / ell : 0.0;

  return signal + ((x == y) ? test->noise_level[x] : 0.0);
}

/* Covariance of band power @i of spectrum @s and band power @j of spectrum @t, term by term. */
static gdouble
_test_nc_xcor_cov_gauss_reference (TestNcXcorCovGauss *test, gboolean with_signal, gdouble f_sky, guint s, guint i, guint t, guint j)
{
  const guint a = test->pair_a[s];
  const guint b = test->pair_b[s];
  const guint c = test->pair_a[t];
  const guint d = test->pair_b[t];
  gdouble res   = 0.0;
  guint l;

  for (l = 0; l < TEST_N_ELL; l++)
  {
    const guint ell  = TEST_ELL_MIN + l;
    const gdouble ac = _test_nc_xcor_cov_gauss_cl_tot (test, with_signal, a, c, ell);
    const gdouble bd = _test_nc_xcor_cov_gauss_cl_tot (test, with_signal, b, d, ell);
    const gdouble ad = _test_nc_xcor_cov_gauss_cl_tot (test, with_signal, a, d, ell);
    const gdouble bc = _test_nc_xcor_cov_gauss_cl_tot (test, with_signal, b, c, ell);

    res += ncm_matrix_get (test->window[s], i, l) * ncm_matrix_get (test->window[t], j, l) * (ac * bd + ad * bc) / ((2.0 * ell + 1.0) * f_sky);
  }

  return res;
}

/* Compares @cov_mat with the reference at every entry, to 1e-13 relative to the largest entry. */
static void
_test_nc_xcor_cov_gauss_check (TestNcXcorCovGauss *test, NcmMatrix *cov_mat, gboolean with_signal, gdouble f_sky)
{
  const guint n_data = nc_xcor_cov_gauss_get_n_data (test->cov);
  NcmMatrix *ref     = ncm_matrix_new (n_data, n_data);
  guint row_s        = 0;
  guint s, t, i, j;

  for (s = 0; s < TEST_N_SPECTRA; s++)
  {
    guint col_t = 0;

    for (t = 0; t < TEST_N_SPECTRA; t++)
    {
      for (i = 0; i < ncm_matrix_nrows (test->window[s]); i++)
        for (j = 0; j < ncm_matrix_nrows (test->window[t]); j++)
          ncm_matrix_set (ref, row_s + i, col_t + j, _test_nc_xcor_cov_gauss_reference (test, with_signal, f_sky, s, i, t, j));

      col_t += ncm_matrix_nrows (test->window[t]);
    }

    row_s += ncm_matrix_nrows (test->window[s]);
  }

  {
    gdouble scale = 0.0;
    guint k;

    for (i = 0; i < n_data; i++)
      for (k = 0; k < n_data; k++)
        scale = GSL_MAX (scale, fabs (ncm_matrix_get (ref, i, k)));

    for (i = 0; i < n_data; i++)
      for (k = 0; k < n_data; k++)
        ncm_assert_cmpdouble_e (ncm_matrix_get (cov_mat, i, k), ==, ncm_matrix_get (ref, i, k), 1.0e-13, 1.0e-15 * scale);
  }

  ncm_matrix_free (ref);
}

/* Fills the signal matrix of the fixture, with or without the signal. */
static void
_test_nc_xcor_cov_gauss_fill (TestNcXcorCovGauss *test, gboolean with_signal)
{
  guint a, b, l;

  for (a = 0; a < TEST_N_TRACERS; a++)
  {
    for (b = a; b < TEST_N_TRACERS; b++)
    {
      const guint row = nc_xcor_cov_gauss_get_pair_index (test->cov, a, b);

      for (l = 0; l < TEST_N_ELL; l++)
        ncm_matrix_set (test->cl, row, l, with_signal ? test->amp[a][b] * 10.0 / (TEST_ELL_MIN + l) : 0.0);
    }
  }
}

void
test_nc_xcor_cov_gauss_new (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  const gdouble amp[TEST_N_TRACERS][TEST_N_TRACERS] = {
    { 2.0,  0.6, -0.4},
    { 0.6,  3.0,  0.5},
    {-0.4,  0.5,  1.5},
  };
  const gdouble noise_level[TEST_N_TRACERS] = {0.1, 0.2, 0.05};
  const guint pair_a[TEST_N_SPECTRA]        = {0, 1, 0, 1, 2, 2};
  const guint pair_b[TEST_N_SPECTRA]        = {0, 0, 2, 1, 1, 2};
  guint s, a, l;

  memcpy (test->amp, amp, sizeof (amp));
  memcpy (test->noise_level, noise_level, sizeof (noise_level));
  memcpy (test->pair_a, pair_a, sizeof (pair_a));
  memcpy (test->pair_b, pair_b, sizeof (pair_b));

  test->window[0] = _test_nc_xcor_cov_gauss_window_tophat (3, 10);
  test->window[1] = _test_nc_xcor_cov_gauss_window_tophat (3, 10);
  test->window[2] = _test_nc_xcor_cov_gauss_window_overlap ();
  test->window[3] = _test_nc_xcor_cov_gauss_window_tophat (3, 10);
  test->window[4] = _test_nc_xcor_cov_gauss_window_tophat (2, 15);
  test->window[5] = _test_nc_xcor_cov_gauss_window_tophat (3, 10);

  test->cov = nc_xcor_cov_gauss_new (TEST_N_TRACERS, TEST_ELL_MIN, TEST_ELL_MAX, TEST_F_SKY);

  for (s = 0; s < TEST_N_SPECTRA; s++)
    nc_xcor_cov_gauss_add_spectrum (test->cov, test->pair_a[s], test->pair_b[s], test->window[s]);

  test->cl      = ncm_matrix_new (nc_xcor_cov_gauss_get_n_pairs (test->cov), TEST_N_ELL);
  test->noise   = ncm_matrix_new (TEST_N_TRACERS, TEST_N_ELL);
  test->cov_mat = ncm_matrix_new (nc_xcor_cov_gauss_get_n_data (test->cov), nc_xcor_cov_gauss_get_n_data (test->cov));

  _test_nc_xcor_cov_gauss_fill (test, TRUE);

  for (a = 0; a < TEST_N_TRACERS; a++)
    for (l = 0; l < TEST_N_ELL; l++)
      ncm_matrix_set (test->noise, a, l, test->noise_level[a]);
}

void
test_nc_xcor_cov_gauss_free (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  guint s;

  for (s = 0; s < TEST_N_SPECTRA; s++)
    ncm_matrix_free (test->window[s]);

  ncm_matrix_free (test->cl);
  ncm_matrix_free (test->noise);
  ncm_matrix_free (test->cov_mat);
  NCM_TEST_FREE (nc_xcor_cov_gauss_free, test->cov);
}

void
test_nc_xcor_cov_gauss_pair_index (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  guint a, b, row = 0;

  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_pairs (test->cov), ==, 6);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_ell (test->cov), ==, TEST_N_ELL);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_spectra (test->cov), ==, TEST_N_SPECTRA);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_data (test->cov), ==, 16);

  for (a = 0; a < TEST_N_TRACERS; a++)
  {
    for (b = a; b < TEST_N_TRACERS; b++)
    {
      g_assert_cmpuint (nc_xcor_cov_gauss_get_pair_index (test->cov, a, b), ==, row);
      g_assert_cmpuint (nc_xcor_cov_gauss_get_pair_index (test->cov, b, a), ==, row);
      row++;
    }
  }
}

void
test_nc_xcor_cov_gauss_formula (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, test->cov_mat);
  _test_nc_xcor_cov_gauss_check (test, test->cov_mat, TRUE, TEST_F_SKY);
}

/*
 * The first band power of (0, 0) against the first of (0, 2): the term is
 * 2 Ct^00 C^02 with C^02 < 0, so the covariance is negative. A form built from |C^ad C^bc|
 * gives a positive value here.
 */
void
test_nc_xcor_cov_gauss_negative_cross (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  const guint offset_02 = 6;

  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, test->cov_mat);

  g_assert_cmpfloat (ncm_matrix_get (test->cov_mat, 0, offset_02), <, 0.0);
  ncm_assert_cmpdouble_e (ncm_matrix_get (test->cov_mat, 0, offset_02), ==,
                          _test_nc_xcor_cov_gauss_reference (test, TRUE, TEST_F_SKY, 0, 0, 2, 0), 1.0e-13, 0.0);
}

/*
 * With no signal only the noise of the auto spectra remains: (0, 0) against (1, 1) is
 * exactly zero, and (1, 0) against itself is N^0 N^1 summed with the windows.
 */
void
test_nc_xcor_cov_gauss_noise_only (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  const guint offset_10 = 3;
  const guint offset_11 = 8;

  _test_nc_xcor_cov_gauss_fill (test, FALSE);
  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, test->cov_mat);

  _test_nc_xcor_cov_gauss_check (test, test->cov_mat, FALSE, TEST_F_SKY);
  g_assert_cmpfloat (ncm_matrix_get (test->cov_mat, 0, offset_11), ==, 0.0);
  g_assert_cmpfloat (ncm_matrix_get (test->cov_mat, offset_10, offset_10), >, 0.0);
}

/* Halving f_sky doubles every entry; both factors are powers of two, so the check is exact. */
void
test_nc_xcor_cov_gauss_f_sky (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  const guint n_data = nc_xcor_cov_gauss_get_n_data (test->cov);
  NcmMatrix *half    = ncm_matrix_new (n_data, n_data);
  guint i, j;

  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, test->cov_mat);

  nc_xcor_cov_gauss_set_f_sky (test->cov, 0.5 * TEST_F_SKY);
  g_assert_cmpfloat (nc_xcor_cov_gauss_get_f_sky (test->cov), ==, 0.5 * TEST_F_SKY);
  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, half);

  for (i = 0; i < n_data; i++)
    for (j = 0; j < n_data; j++)
      g_assert_cmpfloat (ncm_matrix_get (half, i, j), ==, 2.0 * ncm_matrix_get (test->cov_mat, i, j));

  ncm_matrix_free (half);
}

/* The signal plus noise matrix is positive definite at every multipole, so is the covariance. */
void
test_nc_xcor_cov_gauss_positive_definite (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  const guint n_data = nc_xcor_cov_gauss_get_n_data (test->cov);
  guint i, j;

  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, test->cov_mat);

  for (i = 0; i < n_data; i++)
    for (j = 0; j < n_data; j++)
      g_assert_cmpfloat (ncm_matrix_get (test->cov_mat, i, j), ==, ncm_matrix_get (test->cov_mat, j, i));

  g_assert_cmpint (ncm_matrix_cholesky_decomp (test->cov_mat, 'U'), ==, 0);
}

void
test_nc_xcor_cov_gauss_serialize (TestNcXcorCovGauss *test, gconstpointer pdata)
{
  NcXcorCovGauss *dup = NC_XCOR_COV_GAUSS (ncm_serialize_global_dup_obj (G_OBJECT (test->cov)));
  const guint n_data  = nc_xcor_cov_gauss_get_n_data (test->cov);
  NcmMatrix *cov_dup  = ncm_matrix_new (n_data, n_data);
  guint i, j;

  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_tracers (dup), ==, TEST_N_TRACERS);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_ell_min (dup), ==, TEST_ELL_MIN);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_ell_max (dup), ==, TEST_ELL_MAX);
  g_assert_cmpfloat (nc_xcor_cov_gauss_get_f_sky (dup), ==, TEST_F_SKY);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_spectra (dup), ==, TEST_N_SPECTRA);
  g_assert_cmpuint (nc_xcor_cov_gauss_get_n_data (dup), ==, n_data);

  nc_xcor_cov_gauss_compute (test->cov, test->cl, test->noise, test->cov_mat);
  nc_xcor_cov_gauss_compute (dup, test->cl, test->noise, cov_dup);

  for (i = 0; i < n_data; i++)
    for (j = 0; j < n_data; j++)
      g_assert_cmpfloat (ncm_matrix_get (cov_dup, i, j), ==, ncm_matrix_get (test->cov_mat, i, j));

  ncm_matrix_free (cov_dup);
  NCM_TEST_FREE (nc_xcor_cov_gauss_free, dup);
}

void
test_nc_xcor_cov_gauss_invalid (void)
{
  g_test_trap_subprocess ("/nc/xcor/cov_gauss/invalid/tracer/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nc_xcor_cov_gauss_add_spectrum: tracer pair (0, 3) out of range, 3 tracers.*");

  g_test_trap_subprocess ("/nc/xcor/cov_gauss/invalid/window/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nc_xcor_cov_gauss_add_spectrum: the window has 29 columns, the multipole range [2, 31] has 30.*");

  g_test_trap_subprocess ("/nc/xcor/cov_gauss/invalid/cl/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nc_xcor_cov_gauss_compute: the signal matrix is 5 x 30, expected 6 x 30.*");

  g_test_trap_subprocess ("/nc/xcor/cov_gauss/invalid/f_sky/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nc_xcor_cov_gauss_set_f_sky: f-sky must be in (0, 1], got 1.5.*");

  g_test_trap_subprocess ("/nc/xcor/cov_gauss/invalid/ell_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*_nc_xcor_cov_gauss_constructed: ell-max (9) is smaller than ell-min (10).*");
}

void
test_nc_xcor_cov_gauss_invalid_tracer_subprocess (void)
{
  NcXcorCovGauss *cov = nc_xcor_cov_gauss_new (TEST_N_TRACERS, TEST_ELL_MIN, TEST_ELL_MAX, TEST_F_SKY);
  NcmMatrix *window   = _test_nc_xcor_cov_gauss_window_tophat (3, 10);

  nc_xcor_cov_gauss_add_spectrum (cov, 0, TEST_N_TRACERS, window);
}

void
test_nc_xcor_cov_gauss_invalid_window_subprocess (void)
{
  NcXcorCovGauss *cov = nc_xcor_cov_gauss_new (TEST_N_TRACERS, TEST_ELL_MIN, TEST_ELL_MAX, TEST_F_SKY);
  NcmMatrix *window   = ncm_matrix_new (3, TEST_N_ELL - 1);

  nc_xcor_cov_gauss_add_spectrum (cov, 0, 1, window);
}

void
test_nc_xcor_cov_gauss_invalid_cl_subprocess (void)
{
  NcXcorCovGauss *cov = nc_xcor_cov_gauss_new (TEST_N_TRACERS, TEST_ELL_MIN, TEST_ELL_MAX, TEST_F_SKY);
  NcmMatrix *window   = _test_nc_xcor_cov_gauss_window_tophat (3, 10);
  NcmMatrix *cl       = ncm_matrix_new (5, TEST_N_ELL);
  NcmMatrix *noise    = ncm_matrix_new (TEST_N_TRACERS, TEST_N_ELL);
  NcmMatrix *cov_mat  = ncm_matrix_new (3, 3);

  nc_xcor_cov_gauss_add_spectrum (cov, 0, 1, window);
  nc_xcor_cov_gauss_compute (cov, cl, noise, cov_mat);
}

void
test_nc_xcor_cov_gauss_invalid_f_sky_subprocess (void)
{
  NcXcorCovGauss *cov = nc_xcor_cov_gauss_new (TEST_N_TRACERS, TEST_ELL_MIN, TEST_ELL_MAX, TEST_F_SKY);

  nc_xcor_cov_gauss_set_f_sky (cov, 1.5);
}

void
test_nc_xcor_cov_gauss_invalid_ell_range_subprocess (void)
{
  NcXcorCovGauss *cov = nc_xcor_cov_gauss_new (TEST_N_TRACERS, 10, 9, TEST_F_SKY);

  nc_xcor_cov_gauss_free (cov);
}

