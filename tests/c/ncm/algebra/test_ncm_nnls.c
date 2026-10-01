/***************************************************************************
 *            test_ncm_nnls.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_nnls.c
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
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define _assert_small(err, tol) g_assert_cmpfloat ((err), <, (tol))

typedef struct _Problem
{
  NcmMatrix *A;
  NcmVector *x_true;
  NcmVector *f;
} Problem;

/* A Gaussian, x_true >= 0 with every other entry zero, f = A x_true + noise */
static Problem
_problem_new (NcmRNG *rng, guint m, guint n, gdouble noise)
{
  Problem p = {ncm_matrix_new (m, n), ncm_vector_new (n), ncm_vector_new (m)};
  guint i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      ncm_matrix_set (p.A, i, j, ncm_rng_gaussian_gen (rng, 0.0, 1.0));

  for (j = 0; j < n; j++)
    ncm_vector_set (p.x_true, j, (j % 2) ? 0.0 : 0.5 + ncm_rng_uniform01_gen (rng));

  ncm_matrix_update_vector (p.A, 'N', 1.0, p.x_true, 0.0, p.f);

  for (i = 0; i < m; i++)
    ncm_vector_addto (p.f, i, noise * ncm_rng_gaussian_gen (rng, 0.0, 1.0));

  return p;
}

static void
_problem_clear (Problem *p)
{
  ncm_matrix_clear (&p->A);
  ncm_vector_clear (&p->x_true);
  ncm_vector_clear (&p->f);
}

/* |f - A x| and the largest KKT violation, relative to |A^T| |f| */
static gdouble
_kkt_violation (Problem *p, NcmVector *x, gdouble *rnorm)
{
  const guint m = ncm_matrix_nrows (p->A), n = ncm_matrix_ncols (p->A);
  NcmVector *r = ncm_vector_dup (p->f);
  NcmVector *g = ncm_vector_new (n);
  gdouble viol  = 0.0, scale;
  guint j;

  ncm_matrix_update_vector (p->A, 'N', -1.0, x, 1.0, r);
  ncm_matrix_update_vector (p->A, 'T', -1.0, r, 0.0, g);
  *rnorm = ncm_vector_dnrm2 (r);
  scale  = sqrt ((gdouble) m) * ncm_vector_dnrm2 (p->f);

  for (j = 0; j < n; j++)
  {
    const gdouble xj = ncm_vector_get (x, j), gj = ncm_vector_get (g, j);

    viol = MAX (viol, MAX (0.0, -xj));
    viol = MAX (viol, (xj > 0.0) ? fabs (gj) / scale : MAX (0.0, -gj / scale));
  }

  ncm_vector_free (r);
  ncm_vector_free (g);

  return viol;
}

static void
test_ncm_nnls_properties (void)
{
  NcmNNLS *nnls = ncm_nnls_new (7, 3);

  g_assert_cmpuint (ncm_nnls_get_nrows (nnls), ==, 7);
  g_assert_cmpuint (ncm_nnls_get_ncols (nnls), ==, 3);
  g_assert_cmpint (ncm_nnls_get_umethod (nnls), ==, NCM_NNLS_UMETHOD_NORMAL);
  g_assert_cmpfloat (ncm_nnls_get_reltol (nnls), ==, GSL_DBL_EPSILON);

  ncm_nnls_set_umethod (nnls, NCM_NNLS_UMETHOD_QR);
  ncm_nnls_set_reltol (nnls, 1.0e-10);
  g_assert_cmpint (ncm_nnls_get_umethod (nnls), ==, NCM_NNLS_UMETHOD_QR);
  g_assert_cmpfloat (ncm_nnls_get_reltol (nnls), ==, 1.0e-10);

  g_assert_true (ncm_nnls_ref (nnls) == nnls);
  ncm_nnls_free (nnls);
  ncm_nnls_clear (&nnls);
  g_assert_null (nnls);
}

/* Every unconstrained method gives a KKT point, and the stored residuals are f - A x */
static void
test_ncm_nnls_solve_kkt (void)
{
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, 1234);
  const guint sizes[][2] = {
    {
      20, 6
    }, {
      40, 15
    }, {
      8, 8
    }
  };
  NcmNNLSUMethod um;
  guint s;

  for (s = 0; s < G_N_ELEMENTS (sizes); s++)
  {
    Problem p = _problem_new (rng, sizes[s][0], sizes[s][1], 0.3);

    for (um = NCM_NNLS_UMETHOD_NORMAL; um < NCM_NNLS_UMETHOD_LEN; um++)
    {
      NcmNNLS *nnls   = ncm_nnls_new (sizes[s][0], sizes[s][1]);
      NcmVector *x    = ncm_vector_new (sizes[s][1]);
      NcmVector *diff = ncm_vector_dup (p.f);
      gdouble rnorm, rnorm_chk, viol;

      ncm_nnls_set_umethod (nnls, um);
      rnorm = ncm_nnls_solve (nnls, p.A, x, p.f);
      viol  = _kkt_violation (&p, x, &rnorm_chk);

      _assert_small (viol, 5.0e-15);
      _assert_small (fabs (rnorm - rnorm_chk) / rnorm_chk, 1.0e-14);

      ncm_matrix_update_vector (p.A, 'N', -1.0, x, 1.0, diff);
      ncm_vector_sub (diff, ncm_nnls_get_residuals (nnls));
      _assert_small (ncm_vector_dnrm2 (diff) / rnorm_chk, 1.0e-14);

      ncm_vector_free (diff);
      ncm_vector_free (x);
      ncm_nnls_free (nnls);
    }

    _problem_clear (&p);
  }

  ncm_rng_free (rng);
}

/* Without noise the solution is x_true */
static void
test_ncm_nnls_solve_exact (void)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, 4321);
  Problem p   = _problem_new (rng, 30, 10, 0.0);
  NcmNNLSUMethod um;

  for (um = NCM_NNLS_UMETHOD_NORMAL; um < NCM_NNLS_UMETHOD_LEN; um++)
  {
    NcmNNLS *nnls = ncm_nnls_new (30, 10);
    NcmVector *x  = ncm_vector_new (10);
    gdouble err   = 0.0;
    guint j;

    ncm_nnls_set_umethod (nnls, um);
    ncm_nnls_solve (nnls, p.A, x, p.f);

    for (j = 0; j < 10; j++)
      err = MAX (err, fabs (ncm_vector_get (x, j) - ncm_vector_get (p.x_true, j)));

    _assert_small (err, 5.0e-15);

    ncm_vector_free (x);
    ncm_nnls_free (nnls);
  }

  _problem_clear (&p);
  ncm_rng_free (rng);
}

/* Lawson-Hanson agrees with the active-set solver, leaves A alone and stores the residuals */
static void
test_ncm_nnls_solve_LH (void)
{
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, 55);
  const guint sizes[][2] = {
    {
      20, 6
    }, {
      40, 15
    }, {
      8, 8
    }
  };
  guint s;

  for (s = 0; s < G_N_ELEMENTS (sizes); s++)
  {
    const guint m   = sizes[s][0], n = sizes[s][1];
    Problem p       = _problem_new (rng, m, n, 0.3);
    NcmMatrix *A0   = ncm_matrix_dup (p.A);
    NcmNNLS *nnls   = ncm_nnls_new (m, n);
    NcmNNLS *ref    = ncm_nnls_new (m, n);
    NcmVector *x    = ncm_vector_new (n);
    NcmVector *x_r  = ncm_vector_new (n);
    NcmVector *diff = ncm_vector_dup (p.f);
    gdouble rnorm, rnorm_chk, err = 0.0;
    guint i, j;

    rnorm = ncm_nnls_solve_LH (nnls, p.A, x, p.f);
    ncm_nnls_solve (ref, p.A, x_r, p.f);

    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        g_assert_cmpfloat (ncm_matrix_get (p.A, i, j), ==, ncm_matrix_get (A0, i, j));

    for (j = 0; j < n; j++)
      err = MAX (err, fabs (ncm_vector_get (x, j) - ncm_vector_get (x_r, j)));

    _assert_small (err, 5.0e-14);
    _assert_small (_kkt_violation (&p, x, &rnorm_chk), 5.0e-15);
    _assert_small (fabs (rnorm - rnorm_chk) / rnorm_chk, 1.0e-14);

    ncm_matrix_update_vector (p.A, 'N', -1.0, x, 1.0, diff);
    ncm_vector_sub (diff, ncm_nnls_get_residuals (nnls));
    _assert_small (ncm_vector_dnrm2 (diff) / rnorm_chk, 1.0e-14);

    ncm_vector_free (diff);
    ncm_vector_free (x);
    ncm_vector_free (x_r);
    ncm_nnls_free (nnls);
    ncm_nnls_free (ref);
    ncm_matrix_free (A0);
    _problem_clear (&p);
  }

  ncm_rng_free (rng);
}

static void
test_ncm_nnls_solve_lowrankqp (void)
{
  NcmRNG *rng    = ncm_rng_seeded_new (NULL, 99);
  Problem p      = _problem_new (rng, 20, 6, 0.3);
  NcmNNLS *nnls  = ncm_nnls_new (20, 6);
  NcmNNLS *ref   = ncm_nnls_new (20, 6);
  NcmVector *x   = ncm_vector_new (6);
  NcmVector *x_r = ncm_vector_new (6);
  gdouble rnorm, rnorm_chk, err = 0.0;
  guint j;

  rnorm = ncm_nnls_solve_lowrankqp (nnls, p.A, x, p.f);
  ncm_nnls_solve (ref, p.A, x_r, p.f);
  _kkt_violation (&p, x, &rnorm_chk);

  for (j = 0; j < 6; j++)
    err = MAX (err, fabs (ncm_vector_get (x, j) - ncm_vector_get (x_r, j)));

  _assert_small (err, 3.0e-9);
  _assert_small (fabs (rnorm - rnorm_chk) / rnorm_chk, 1.0e-14);

  ncm_vector_free (x);
  ncm_vector_free (x_r);
  ncm_nnls_free (nnls);
  ncm_nnls_free (ref);
  _problem_clear (&p);
  ncm_rng_free (rng);
}

/* KKT of min |A x - f| with x >= 0 and sum x fixed: g = A^T (A x - f) equals mu on the
 * positive entries and is at least mu on the others. Returns the violation relative to
 * |A^T f| and sets mu and sum x. */
static gdouble
_simplex_kkt_violation (Problem *p, NcmVector *x, gdouble *mu, gdouble *sum)
{
  const guint n = ncm_matrix_ncols (p->A);
  NcmVector *r  = ncm_vector_dup (p->f);
  NcmVector *g  = ncm_vector_new (n);
  NcmVector *b  = ncm_vector_new (n);
  gdouble g_min = GSL_POSINF, g_max_pos = GSL_NEGINF, scale;
  guint j;

  ncm_matrix_update_vector (p->A, 'N', -1.0, x, 1.0, r);
  ncm_matrix_update_vector (p->A, 'T', -1.0, r, 0.0, g);
  ncm_matrix_update_vector (p->A, 'T', 1.0, p->f, 0.0, b);
  scale = ncm_vector_dnrm2 (b);
  *sum  = 0.0;

  for (j = 0; j < n; j++)
  {
    const gdouble xj = ncm_vector_get (x, j), gj = ncm_vector_get (g, j);

    g_assert_cmpfloat (xj, >=, 0.0);
    *sum += xj;
    g_min = MIN (g_min, gj);

    if (xj > 0.0)
      g_max_pos = MAX (g_max_pos, gj);
  }

  *mu = g_max_pos;

  ncm_vector_free (r);
  ncm_vector_free (g);
  ncm_vector_free (b);

  return (g_max_pos - g_min) / scale;
}

/* gsmo: sum x = 1; splx: sum x <= 1, with mu <= 0, and the NNLS solution when that is inside */
static void
test_ncm_nnls_solve_simplex (void)
{
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, 7);
  const guint sizes[][2] = {
    {
      20, 6
    }, {
      40, 15
    }
  };
  guint s;

  for (s = 0; s < G_N_ELEMENTS (sizes); s++)
  {
    const guint m = sizes[s][0], n = sizes[s][1];
    Problem p     = _problem_new (rng, m, n, 0.3);
    NcmNNLS *nnls = ncm_nnls_new (m, n);
    NcmVector *x  = ncm_vector_new (n);
    gdouble mu, sum;

    ncm_nnls_solve_gsmo (nnls, p.A, x, p.f);
    _assert_small (_simplex_kkt_violation (&p, x, &mu, &sum), 2.0e-15);
    _assert_small (fabs (sum - 1.0), 1.0e-14);

    ncm_nnls_solve_splx (nnls, p.A, x, p.f);
    _assert_small (_simplex_kkt_violation (&p, x, &mu, &sum), 2.0e-15);
    _assert_small (fabs (sum - 1.0), 1.0e-14);
    g_assert_cmpfloat (mu, <=, 0.0);

    ncm_vector_free (x);
    ncm_nnls_free (nnls);
    _problem_clear (&p);
  }

  /* x_true sums to 0.25, so the constraint of splx is inactive */
  {
    Problem p      = _problem_new (rng, 30, 8, 0.0);
    NcmNNLS *nnls  = ncm_nnls_new (30, 8);
    NcmNNLS *ref   = ncm_nnls_new (30, 8);
    NcmVector *x   = ncm_vector_new (8);
    NcmVector *x_r = ncm_vector_new (8);
    gdouble err    = 0.0;
    guint j;

    ncm_vector_scale (p.x_true, 0.25 / ncm_vector_sum_cpts (p.x_true));
    ncm_matrix_update_vector (p.A, 'N', 1.0, p.x_true, 0.0, p.f);

    ncm_nnls_solve_splx (nnls, p.A, x, p.f);
    ncm_nnls_solve (ref, p.A, x_r, p.f);

    for (j = 0; j < 8; j++)
      err = MAX (err, fabs (ncm_vector_get (x, j) - ncm_vector_get (x_r, j)));

    _assert_small (err, 1.0e-15);

    ncm_vector_free (x);
    ncm_vector_free (x_r);
    ncm_nnls_free (nnls);
    ncm_nnls_free (ref);
    _problem_clear (&p);
  }

  ncm_rng_free (rng);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/nnls/properties", &test_ncm_nnls_properties);
  g_test_add_func ("/ncm/nnls/solve/kkt", &test_ncm_nnls_solve_kkt);
  g_test_add_func ("/ncm/nnls/solve/exact", &test_ncm_nnls_solve_exact);
  g_test_add_func ("/ncm/nnls/solve/LH", &test_ncm_nnls_solve_LH);
  g_test_add_func ("/ncm/nnls/solve/lowrankqp", &test_ncm_nnls_solve_lowrankqp);
  g_test_add_func ("/ncm/nnls/solve/simplex", &test_ncm_nnls_solve_simplex);

  g_test_run ();
}

