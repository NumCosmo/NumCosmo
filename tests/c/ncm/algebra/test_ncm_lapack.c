/***************************************************************************
 *            test_ncm_lapack.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_lapack.c
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
#include <string.h>

/* A symmetric positive definite 3x3 matrix, row-major */
static const gdouble A_spd[9] = {4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0};

/* A general 3x3 matrix, row-major */
static const gdouble A_gen[9] = {2.0, 1.0, 0.0, -1.0, 3.0, 1.0, 0.5, 0.0, 1.0};

static void
_test_lapack_matvec (const gdouble *A, const gdouble *x, gdouble *y, gboolean trans)
{
  guint i, j;

  for (i = 0; i < 3; i++)
  {
    y[i] = 0.0;

    for (j = 0; j < 3; j++)
      y[i] += (trans ? A[j * 3 + i] : A[i * 3 + j]) * x[j];
  }
}

static void
test_ncm_lapack_dsyevd (void)
{
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gdouble a[9], w[3];
  guint k;

  memcpy (a, A_spd, sizeof (a));
  g_assert_cmpint (ncm_lapack_dsyevd ('V', 'U', 3, a, 3, w, ws), ==, 0);

  g_assert_cmpfloat (w[0], <=, w[1]);
  g_assert_cmpfloat (w[1], <=, w[2]);

  /* The eigenvectors are the rows of a */
  for (k = 0; k < 3; k++)
  {
    gdouble Av[3];
    guint i;

    _test_lapack_matvec (A_spd, &a[k * 3], Av, FALSE);

    for (i = 0; i < 3; i++)
      ncm_assert_cmpdouble_e (Av[i], ==, w[k] * a[k * 3 + i], 1.0e-13, 1.0e-14);
  }

  ncm_lapack_ws_free (ws);
}

static void
test_ncm_lapack_dgeev (void)
{
  gdouble a[9], wr[3], wi[3], vl[9], vr[9], work[64];
  guint k;

  memcpy (a, A_gen, sizeof (a));
  g_assert_cmpint (ncm_lapack_dgeev ('V', 'V', 3, a, 3, wr, wi, vl, 3, vr, 3, work, 64), ==, 0);

  /* For the real eigenvalues, the rows of vr are right eigenvectors of the row-major matrix */
  for (k = 0; k < 3; k++)
  {
    gdouble Av[3];
    guint i;

    if (wi[k] != 0.0)
      continue;

    _test_lapack_matvec (A_gen, &vr[k * 3], Av, FALSE);

    for (i = 0; i < 3; i++)
      ncm_assert_cmpdouble_e (Av[i], ==, wr[k] * vr[k * 3 + i], 1.0e-12, 1.0e-13);
  }
}

static void
test_ncm_lapack_dgesv (void)
{
  const gdouble b0[3] = {1.0, -2.0, 0.5};
  gdouble a[9], b[3], Ax[3];
  gint ipiv[3];
  guint i;

  /* The array is read in column-major order: given row-major A it solves A^T x = b */
  memcpy (a, A_gen, sizeof (a));
  memcpy (b, b0, sizeof (b));
  g_assert_cmpint (ncm_lapack_dgesv (3, 1, a, 3, ipiv, b, 3), ==, 0);

  _test_lapack_matvec (A_gen, b, Ax, TRUE);

  for (i = 0; i < 3; i++)
    ncm_assert_cmpdouble_e (Ax[i], ==, b0[i], 1.0e-14, 1.0e-15);
}

static void
test_ncm_lapack_dsytrf_dsytrs (void)
{
  NcmLapackWS *ws     = ncm_lapack_ws_new ();
  const gdouble b0[3] = {1.0, -2.0, 0.5};
  gdouble a[9], b[3], Ax[3];
  gint ipiv[3];
  guint i;

  memcpy (a, A_spd, sizeof (a));
  memcpy (b, b0, sizeof (b));
  g_assert_cmpint (ncm_lapack_dsytrf ('U', 3, a, 3, ipiv, ws), ==, 0);
  g_assert_cmpint (ncm_lapack_dsytrs ('U', 3, 1, a, 3, ipiv, b, 3), ==, 0);

  _test_lapack_matvec (A_spd, b, Ax, FALSE);

  for (i = 0; i < 3; i++)
    ncm_assert_cmpdouble_e (Ax[i], ==, b0[i], 1.0e-14, 1.0e-15);

  ncm_lapack_ws_free (ws);
}

static void
test_ncm_lapack_dptsv (void)
{
  /* Tridiagonal: diagonal 4, off-diagonal 1 */
  gdouble d[4]        = {4.0, 4.0, 4.0, 4.0};
  gdouble e[3]        = {1.0, 1.0, 1.0};
  const gdouble b0[4] = {1.0, 2.0, 3.0, 4.0};
  gdouble b[4], x[4];
  guint i;

  memcpy (b, b0, sizeof (b));
  g_assert_cmpint (ncm_lapack_dptsv (d, e, b, x, 4), ==, 0);

  for (i = 0; i < 4; i++)
  {
    const gdouble Ax = 4.0 * x[i] + ((i > 0) ? x[i - 1] : 0.0) + ((i < 3) ? x[i + 1] : 0.0);

    ncm_assert_cmpdouble_e (Ax, ==, b0[i], 1.0e-14, 1.0e-15);
  }
}

static void
test_ncm_lapack_dgeqrf (void)
{
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gdouble a[9], tau[3];
  guint i, j, k;

  /* QR of the row-major matrix: R is the upper triangle and A^T A = R^T R */
  memcpy (a, A_gen, sizeof (a));
  g_assert_cmpint (ncm_lapack_dgeqrf (3, 3, a, 3, tau, ws), ==, 0);

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      gdouble AtA = 0.0, RtR = 0.0;

      for (k = 0; k < 3; k++)
      {
        AtA += A_gen[k * 3 + i] * A_gen[k * 3 + j];

        if ((k <= i) && (k <= j))
          RtR += a[k * 3 + i] * a[k * 3 + j];
      }

      ncm_assert_cmpdouble_e (RtR, ==, AtA, 1.0e-13, 1.0e-14);
    }
  }

  ncm_lapack_ws_free (ws);
}

static void
test_ncm_lapack_dgels (void)
{
  /* Least squares fit of y = c0 + c1 t to four points, row-major 4x2 design matrix */
  const gdouble t[4] = {0.0, 1.0, 2.0, 3.0};
  const gdouble y[4] = {1.0, 2.9, 5.1, 7.0};
  gdouble a[8], b[4], work[64];
  gdouble S0 = 0.0, S1 = 0.0, S2 = 0.0, Y0 = 0.0, Y1 = 0.0;
  guint i;

  for (i = 0; i < 4; i++)
  {
    a[i * 2 + 0] = 1.0;
    a[i * 2 + 1] = t[i];
    b[i]         = y[i];
    S0          += 1.0;
    S1          += t[i];
    S2          += t[i] * t[i];
    Y0          += y[i];
    Y1          += t[i] * y[i];
  }

  g_assert_cmpint (ncm_lapack_dgels ('N', 4, 2, 1, a, 2, b, 4, work, 64), ==, 0);

  /* Normal equations */
  {
    const gdouble det = S0 * S2 - S1 * S1;
    const gdouble c0  = (S2 * Y0 - S1 * Y1) / det;
    const gdouble c1  = (S0 * Y1 - S1 * Y0) / det;

    ncm_assert_cmpdouble_e (b[0], ==, c0, 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (b[1], ==, c1, 1.0e-13, 0.0);
  }
}

static void
test_ncm_lapack_ws (void)
{
  NcmLapackWS *ws  = ncm_lapack_ws_new ();
  NcmLapackWS *ws2 = ncm_lapack_ws_dup (ws);

  ncm_lapack_ws_free (ws);
  ncm_lapack_ws_clear (&ws2);
  g_assert_null (ws2);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/lapack/ws", &test_ncm_lapack_ws);
#ifdef HAVE_LAPACK
  g_test_add_func ("/ncm/lapack/dsyevd", &test_ncm_lapack_dsyevd);
  g_test_add_func ("/ncm/lapack/dgeev", &test_ncm_lapack_dgeev);
  g_test_add_func ("/ncm/lapack/dgesv", &test_ncm_lapack_dgesv);
  g_test_add_func ("/ncm/lapack/dsytrf_dsytrs", &test_ncm_lapack_dsytrf_dsytrs);
  g_test_add_func ("/ncm/lapack/dptsv", &test_ncm_lapack_dptsv);
  g_test_add_func ("/ncm/lapack/dgeqrf", &test_ncm_lapack_dgeqrf);
  g_test_add_func ("/ncm/lapack/dgels", &test_ncm_lapack_dgels);
#endif /* HAVE_LAPACK */

  g_test_run ();
}

