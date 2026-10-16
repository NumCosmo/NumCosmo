/***************************************************************************
 *            test_ncm_spectral_panels_product.c
 *
 *  Sat Oct 10 2026
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/* The closure error has to stand well above the rounding of the exact integrals for the
 * error estimate to be what is checked, so the panels are built at 1e-8, not near 1e-16 */
#define TEST_RELTOL (1.0e-8)

/*
 * A side: on [lo, hi] with one breakpoint bp, row r is factor[r][i] exp (rate[r] x) on the
 * continuous interval i, so it jumps by the ratio of the factors at bp.
 */
typedef struct _Side
{
  gdouble lo;
  gdouble hi;
  gdouble bp;
  gdouble rate[2];
  gdouble factor[2][2];
} Side;

typedef struct _ThreadProduct
{
  NcmSpectralPanelsProduct *prod;
  NcmSpectralPanels *sp_P;
  NcmSpectralPanels *sp_Q;
  NcmVector *s;
  NcmVector *result;
} ThreadProduct;

static void test_ncm_spectral_panels_product_p (void);
static void test_ncm_spectral_panels_product_cross (void);
static void test_ncm_spectral_panels_product_auto (void);
static void test_ncm_spectral_panels_product_scaled (void);
static void test_ncm_spectral_panels_product_mixed (void);
static void test_ncm_spectral_panels_product_large_tables (void);
static void test_ncm_spectral_panels_product_disjoint (void);
static void test_ncm_spectral_panels_product_narrow (void);
static void test_ncm_spectral_panels_product_threads (void);
static void test_ncm_spectral_panels_product_fdim (void);
static void test_ncm_spectral_panels_product_fdim_subprocess (void);
static void test_ncm_spectral_panels_product_lengths (void);
static void test_ncm_spectral_panels_product_lengths_subprocess (void);
static void test_ncm_spectral_panels_product_rows (void);
static void test_ncm_spectral_panels_product_rows_fit (void);
static void test_ncm_spectral_panels_product_rows_fit_subprocess (void);

static const Side side_P = {
  0.5, 3.0, 1.3, {
    1.0, 0.5
  }, {
    {
      1.0, 3.0
    }, {
      1.0, 2.0
    }
  }
};

static const Side side_Q = {
  0.2, 2.6, 2.1, {
    -2.0, 1.0
  }, {
    {
      1.0, 0.5
    }, {
      1.0, 0.25
    }
  }
};

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/spectral_panels_product/p", &test_ncm_spectral_panels_product_p);
  g_test_add_func ("/ncm/spectral_panels_product/cross", &test_ncm_spectral_panels_product_cross);
  g_test_add_func ("/ncm/spectral_panels_product/auto", &test_ncm_spectral_panels_product_auto);
  g_test_add_func ("/ncm/spectral_panels_product/scaled", &test_ncm_spectral_panels_product_scaled);
  g_test_add_func ("/ncm/spectral_panels_product/mixed", &test_ncm_spectral_panels_product_mixed);
  g_test_add_func ("/ncm/spectral_panels_product/large_tables", &test_ncm_spectral_panels_product_large_tables);
  g_test_add_func ("/ncm/spectral_panels_product/disjoint", &test_ncm_spectral_panels_product_disjoint);
  g_test_add_func ("/ncm/spectral_panels_product/narrow", &test_ncm_spectral_panels_product_narrow);
  g_test_add_func ("/ncm/spectral_panels_product/threads", &test_ncm_spectral_panels_product_threads);
  g_test_add_func ("/ncm/spectral_panels_product/fdim", &test_ncm_spectral_panels_product_fdim);
  g_test_add_func ("/ncm/spectral_panels_product/fdim/subprocess", &test_ncm_spectral_panels_product_fdim_subprocess);
  g_test_add_func ("/ncm/spectral_panels_product/lengths", &test_ncm_spectral_panels_product_lengths);
  g_test_add_func ("/ncm/spectral_panels_product/lengths/subprocess", &test_ncm_spectral_panels_product_lengths_subprocess);
  g_test_add_func ("/ncm/spectral_panels_product/rows", &test_ncm_spectral_panels_product_rows);
  g_test_add_func ("/ncm/spectral_panels_product/rows_fit", &test_ncm_spectral_panels_product_rows_fit);
  g_test_add_func ("/ncm/spectral_panels_product/rows_fit/subprocess", &test_ncm_spectral_panels_product_rows_fit_subprocess);

  g_test_run ();

  return 0;
}

static void
_side_F (gpointer user_data, guint i, gdouble x, NcmVector *fval)
{
  const Side *side = user_data;
  guint r;

  for (r = 0; r < 2; r++)
    ncm_vector_set (fval, r, side->factor[r][i] * exp (side->rate[r] * x));
}

static NcmSpectralPanels *
_side_panels_new (const Side *side, guint level_min, guint level_cap)
{
  NcmSpectralPanelsBuilder *builder = ncm_spectral_panels_builder_new ();
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (2);
  NcmVector *breakpoints            = ncm_vector_new (1);

  ncm_vector_set (breakpoints, 0, side->bp);
  ncm_spectral_panels_builder_set_level_cap (builder, level_cap);
  ncm_spectral_panels_builder_set_level_min (builder, level_min);
  ncm_spectral_panels_builder_build (builder, _side_F, (gpointer) side, side->lo, side->hi, breakpoints, NULL, NULL,
                                     TEST_RELTOL, 0.0, sp);

  ncm_vector_free (breakpoints);
  ncm_spectral_panels_builder_free (builder);

  return sp;
}

/* int_{x0}^{x1} x^p exp (a x) dx from the antiderivative exp (a x) sum_j (-1)^j p! / (p - j)! x^(p - j) / a^(j + 1) */
static gdouble
_int_xp_exp (guint p, gdouble a, gdouble x0, gdouble x1)
{
  gdouble F0 = 0.0, F1 = 0.0, coef = 1.0 / a;
  guint j;

  for (j = 0; j <= p; j++)
  {
    F0  += coef * gsl_pow_int (x0, p - j);
    F1  += coef * gsl_pow_int (x1, p - j);
    coef = -coef * (p - j) / a;
  }

  return exp (a * x1) * F1 - exp (a * x0) * F0;
}

static gint
_cmp_gdouble (gconstpointer pa, gconstpointer pb)
{
  const gdouble a = *(const gdouble *) pa;
  const gdouble b = *(const gdouble *) pb;

  return (a > b) - (a < b);
}

/*
 * The exact int k^p P_r1(k / s1) Q_r2(k / s2) dk over the overlap, piece by piece between
 * the breakpoints in k; with the largest |P| and |Q| there and int k^p dk, for the bound
 * of the error estimate.
 */
static gdouble
_exact (const Side *P, gdouble s1, guint r1, const Side *Q, gdouble s2, guint r2, guint p, gdouble *max_P, gdouble *max_Q, gdouble *int_kp)
{
  const gdouble k_lo = MAX (P->lo * s1, Q->lo * s2);
  const gdouble k_hi = MIN (P->hi * s1, Q->hi * s2);
  const gdouble a    = P->rate[r1] / s1 + Q->rate[r2] / s2;
  GArray *cuts       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble bp[] = {P->bp * s1, Q->bp * s2};
  gdouble total      = 0.0;
  guint j;

  *max_P  = 0.0;
  *max_Q  = 0.0;
  *int_kp = 0.0;

  if (k_lo >= k_hi)
  {
    g_array_unref (cuts);

    return 0.0;
  }

  g_array_append_val (cuts, k_lo);

  for (j = 0; j < 2; j++)
    if ((bp[j] > k_lo) && (bp[j] < k_hi))
      g_array_append_val (cuts, bp[j]);

  g_array_append_val (cuts, k_hi);
  g_array_sort (cuts, _cmp_gdouble);

  for (j = 0; j + 1 < cuts->len; j++)
  {
    const gdouble k0  = g_array_index (cuts, gdouble, j);
    const gdouble k1  = g_array_index (cuts, gdouble, j + 1);
    const gdouble km  = 0.5 * (k0 + k1);
    const gdouble f_P = P->factor[r1][(km / s1 < P->bp) ? 0 : 1];
    const gdouble f_Q = Q->factor[r2][(km / s2 < Q->bp) ? 0 : 1];

    total += f_P * f_Q * _int_xp_exp (p, a, k0, k1);

    /* Each piece is monotone, so its largest value is at an end */
    *max_P = MAX (*max_P, MAX (fabs (f_P * exp (P->rate[r1] * k0 / s1)), fabs (f_P * exp (P->rate[r1] * k1 / s1))));
    *max_Q = MAX (*max_Q, MAX (fabs (f_Q * exp (Q->rate[r2] * k0 / s2)), fabs (f_Q * exp (Q->rate[r2] * k1 / s2))));
  }

  *int_kp = (gsl_pow_int (k_hi, p + 1) - gsl_pow_int (k_lo, p + 1)) / (p + 1);

  g_array_unref (cuts);

  return total;
}

static NcmVector *
_scales_new (gdouble s0, gdouble s1)
{
  NcmVector *s = ncm_vector_new (2);

  ncm_vector_set (s, 0, s0);
  ncm_vector_set (s, 1, s1);

  return s;
}

/*
 * Integrates with every p in [0, 3] and checks each row against the exact integral: the
 * error within the estimate, and the estimate within 2 reltol max|P| max|Q| int k^p dk
 */
static void
_check_pair (NcmSpectralPanels *sp_P, const Side *P, NcmVector *s1, NcmSpectralPanels *sp_Q, const Side *Q, NcmVector *s2)
{
  NcmVector *result = ncm_vector_new (2);
  NcmVector *err    = ncm_vector_new (2);
  guint p, r;

  for (p = 0; p <= 3; p++)
  {
    NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (p);

    ncm_spectral_panels_product_integrate (prod, sp_P, s1, sp_Q, s2, result, err);

    for (r = 0; r < 2; r++)
    {
      const gdouble sr1 = (s1 != NULL) ? ncm_vector_get (s1, r) : 1.0;
      const gdouble sr2 = (s2 != NULL) ? ncm_vector_get (s2, r) : 1.0;
      gdouble max_P, max_Q, int_kp;
      const gdouble exact = _exact (P, sr1, r, Q, sr2, r, p, &max_P, &max_Q, &int_kp);

      g_assert_cmpfloat (fabs (ncm_vector_get (result, r) - exact), <=, ncm_vector_get (err, r));
      g_assert_cmpfloat (ncm_vector_get (err, r), <=, 2.0 * TEST_RELTOL * max_P * max_Q * int_kp);
    }

    ncm_spectral_panels_product_free (prod);
  }

  ncm_vector_free (err);
  ncm_vector_free (result);
}

static gpointer
_thread_product (gpointer data)
{
  ThreadProduct *tp = data;
  guint k;

  for (k = 0; k < 50; k++)
    ncm_spectral_panels_product_integrate (tp->prod, tp->sp_P, tp->s, tp->sp_Q, NULL, tp->result, NULL);

  return NULL;
}

static void
test_ncm_spectral_panels_product_p (void)
{
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (3);
  guint p;

  g_assert_cmpuint (ncm_spectral_panels_product_get_p (prod), ==, 3);

  g_object_get (prod, "p", &p, NULL);
  g_assert_cmpuint (p, ==, 3);

  ncm_spectral_panels_product_free (prod);
}

/* Two panel sets with different edges and breakpoints, in the same variable */
static void
test_ncm_spectral_panels_product_cross (void)
{
  NcmSpectralPanels *sp_P = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_Q = _side_panels_new (&side_Q, 3, 5);

  _check_pair (sp_P, &side_P, NULL, sp_Q, &side_Q, NULL);

  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

/* A panel set with itself is the square; the same as the cross product with a copy */
static void
test_ncm_spectral_panels_product_auto (void)
{
  NcmSpectralPanels *sp_P    = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_copy = _side_panels_new (&side_P, 3, 5);
  NcmVector *res_auto        = ncm_vector_new (2);
  NcmVector *res_cross       = ncm_vector_new (2);
  NcmVector *err_auto        = ncm_vector_new (2);
  NcmVector *err_cross       = ncm_vector_new (2);
  guint r;

  _check_pair (sp_P, &side_P, NULL, sp_P, &side_P, NULL);

  {
    NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (2);

    ncm_spectral_panels_product_integrate (prod, sp_P, NULL, sp_P, NULL, res_auto, err_auto);
    ncm_spectral_panels_product_integrate (prod, sp_P, NULL, sp_copy, NULL, res_cross, err_cross);
    ncm_spectral_panels_product_free (prod);
  }

  for (r = 0; r < 2; r++)
  {
    g_assert_cmpfloat (ncm_vector_get (res_auto, r), ==, ncm_vector_get (res_cross, r));
    g_assert_cmpfloat (fabs (ncm_vector_get (err_auto, r) / ncm_vector_get (err_cross, r) - 1.0), <=, 1.0e-14);
  }

  ncm_vector_free (err_cross);
  ncm_vector_free (err_auto);
  ncm_vector_free (res_cross);
  ncm_vector_free (res_auto);
  ncm_spectral_panels_free (sp_copy);
  ncm_spectral_panels_free (sp_P);
}

/* Both sides scaled by the same s_r: one refinement in x = k / s_r */
static void
test_ncm_spectral_panels_product_scaled (void)
{
  NcmSpectralPanels *sp_P = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_Q = _side_panels_new (&side_Q, 3, 5);
  NcmVector *s            = _scales_new (2.5, 3.5);

  _check_pair (sp_P, &side_P, s, sp_Q, &side_Q, s);
  _check_pair (sp_P, &side_P, s, sp_P, &side_P, s);

  ncm_vector_free (s);
  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

/* One side scaled, either one, or both differently: one refinement per row. The scales keep
 * the exponent of every product away from zero, where the exact integral has another form */
static void
test_ncm_spectral_panels_product_mixed (void)
{
  NcmSpectralPanels *sp_P = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_Q = _side_panels_new (&side_Q, 3, 5);
  NcmVector *s            = _scales_new (0.6, 0.75);
  NcmVector *s_other      = _scales_new (0.8, 0.6);

  _check_pair (sp_P, &side_P, s, sp_Q, &side_Q, NULL);
  _check_pair (sp_P, &side_P, NULL, sp_Q, &side_Q, s);
  _check_pair (sp_P, &side_P, s, sp_Q, &side_Q, s_other);

  ncm_vector_free (s_other);
  ncm_vector_free (s);
  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

/* Panels with more coefficients than the object's tables use tables built for the call */
static void
test_ncm_spectral_panels_product_large_tables (void)
{
  NcmSpectralPanels *sp_P = _side_panels_new (&side_P, 7, 8);
  NcmSpectralPanels *sp_Q = _side_panels_new (&side_Q, 3, 5);
  gdouble a, b;
  guint N;

  ncm_spectral_panels_peek_panel (sp_P, 0, &a, &b, &N);
  g_assert_cmpuint (N, >, 65);

  _check_pair (sp_P, &side_P, NULL, sp_Q, &side_Q, NULL);

  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

/* Ranges that do not overlap give zero */
static void
test_ncm_spectral_panels_product_disjoint (void)
{
  NcmSpectralPanels *sp_P        = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_Q        = _side_panels_new (&side_Q, 3, 5);
  NcmVector *s                   = _scales_new (20.0, 20.0);
  NcmVector *result              = ncm_vector_new (2);
  NcmVector *err                 = ncm_vector_new (2);
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (2);
  guint r;

  ncm_vector_set_all (result, 1.0);
  ncm_vector_set_all (err, 1.0);
  ncm_spectral_panels_product_integrate (prod, sp_P, NULL, sp_Q, s, result, err);

  for (r = 0; r < 2; r++)
  {
    g_assert_cmpfloat (ncm_vector_get (result, r), ==, 0.0);
    g_assert_cmpfloat (ncm_vector_get (err, r), ==, 0.0);
  }

  ncm_spectral_panels_product_free (prod);
  ncm_vector_free (err);
  ncm_vector_free (result);
  ncm_vector_free (s);
  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

/*
 * Integrates two constant sides, one in k with a panel of relative width 1e-5 ending at
 * b, the other in u = k / nu, so that b is mapped to u and back. The round trip can
 * differ from b by one unit in the last place, and the cell must still be found inside
 * the panel ending at b. Both sides are 1 on [0.1, 0.5] in k, so the integral is 0.4 for
 * every end.
 */
static void
test_ncm_spectral_panels_product_narrow (void)
{
  const gdouble nu               = 100.5;
  const gdouble ends[]           = {0.123, 0.3, 0.31, 0.29, 0.191, 0.333, 0.271};
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (0);
  NcmMatrix *one                 = ncm_matrix_new (1, 1);
  NcmVector *s                   = ncm_vector_new (1);
  NcmVector *result              = ncm_vector_new (1);
  NcmSpectralPanels *sp_u        = ncm_spectral_panels_new (1);
  guint i;

  ncm_matrix_set (one, 0, 0, 1.0);
  ncm_vector_set (s, 0, nu);
  ncm_spectral_panels_append (sp_u, 0.1 / nu, 0.5 / nu, one);

  for (i = 0; i < G_N_ELEMENTS (ends); i++)
  {
    const gdouble b         = ends[i];
    const gdouble w         = 1.0e-5 * b;
    NcmSpectralPanels *sp_k = ncm_spectral_panels_new (1);

    ncm_spectral_panels_append (sp_k, 0.1, b - w, one);
    ncm_spectral_panels_append (sp_k, b - w, b, one);
    ncm_spectral_panels_append (sp_k, b, 0.5, one);

    ncm_spectral_panels_product_integrate (prod, sp_u, s, sp_k, NULL, result, NULL);
    ncm_assert_cmpdouble_e (ncm_vector_get (result, 0), ==, 0.4, 1.0e-13, 0.0);

    ncm_spectral_panels_product_integrate (prod, sp_k, NULL, sp_u, s, result, NULL);
    ncm_assert_cmpdouble_e (ncm_vector_get (result, 0), ==, 0.4, 1.0e-13, 0.0);

    ncm_spectral_panels_free (sp_k);
  }

  ncm_spectral_panels_free (sp_u);
  ncm_vector_free (result);
  ncm_vector_free (s);
  ncm_matrix_free (one);
  ncm_spectral_panels_product_free (prod);
}

/* Calls in two threads on one object give the serial result */
static void
test_ncm_spectral_panels_product_threads (void)
{
  NcmSpectralPanels *sp_P        = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_Q        = _side_panels_new (&side_Q, 3, 5);
  NcmVector *s                   = _scales_new (0.5, 0.75);
  NcmVector *serial              = ncm_vector_new (2);
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (2);
  ThreadProduct tp[2]            = {
    {
      prod, sp_P, sp_Q, s, ncm_vector_new (2)
    }, {
      prod, sp_P, sp_Q, s, ncm_vector_new (2)
    }
  };
  GThread *threads[2];
  guint t, r;

  ncm_spectral_panels_product_integrate (prod, sp_P, s, sp_Q, NULL, serial, NULL);

  for (t = 0; t < 2; t++)
    threads[t] = g_thread_new ("product", _thread_product, &tp[t]);

  for (t = 0; t < 2; t++)
  {
    g_thread_join (threads[t]);

    for (r = 0; r < 2; r++)
      g_assert_cmpfloat (ncm_vector_get (tp[t].result, r), ==, ncm_vector_get (serial, r));

    ncm_vector_free (tp[t].result);
  }

  ncm_spectral_panels_product_free (prod);
  ncm_vector_free (serial);
  ncm_vector_free (s);
  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

static void
test_ncm_spectral_panels_product_fdim (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_product/fdim/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*have 2 and 1 rows*");
}

static void
test_ncm_spectral_panels_product_fdim_subprocess (void)
{
  NcmSpectralPanels *sp_P        = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_1        = ncm_spectral_panels_new (1);
  NcmVector *result              = ncm_vector_new (2);
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (2);

  ncm_spectral_panels_product_integrate (prod, sp_P, NULL, sp_1, NULL, result, NULL);
}

static void
test_ncm_spectral_panels_product_lengths (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_product/lengths/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*must have fdim = 2 entries*");
}

static void
test_ncm_spectral_panels_product_lengths_subprocess (void)
{
  NcmSpectralPanels *sp_P        = _side_panels_new (&side_P, 3, 5);
  NcmVector *result              = ncm_vector_new (3);
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (2);

  ncm_spectral_panels_product_integrate (prod, sp_P, NULL, sp_P, NULL, result, NULL);
}

/*
 * One row of P against one row of Q at every pair of offsets, and P against its own other
 * row with other scales, each checked against the exact integral as in _check_pair()
 */
static void
test_ncm_spectral_panels_product_rows (void)
{
  NcmSpectralPanels *sp_P = _side_panels_new (&side_P, 3, 5);
  NcmSpectralPanels *sp_Q = _side_panels_new (&side_Q, 3, 5);
  NcmVector *s1           = _scales_new (1.5, 2.5);
  NcmVector *s2           = _scales_new (0.75, 1.25);
  NcmVector *result       = ncm_vector_new (1);
  NcmVector *err          = ncm_vector_new (1);
  guint p, o1, o2, self;

  for (p = 0; p <= 3; p++)
  {
    NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (p);

    for (self = 0; self < 2; self++)
    {
      NcmSpectralPanels *sp_2 = self ? sp_P : sp_Q;
      const Side *side_2      = self ? &side_P : &side_Q;

      for (o1 = 0; o1 < 2; o1++)
      {
        for (o2 = 0; o2 < 2; o2++)
        {
          gdouble max_P, max_Q, int_kp;
          const gdouble exact = _exact (&side_P, ncm_vector_get (s1, o1), o1, side_2, ncm_vector_get (s2, o2), o2, p, &max_P, &max_Q, &int_kp);

          ncm_spectral_panels_product_integrate_rows (prod, sp_P, o1, s1, sp_2, o2, s2, result, err);

          g_assert_cmpfloat (fabs (ncm_vector_get (result, 0) - exact), <=, ncm_vector_get (err, 0));
          g_assert_cmpfloat (ncm_vector_get (err, 0), <=, 2.0 * TEST_RELTOL * max_P * max_Q * int_kp);
        }
      }
    }

    ncm_spectral_panels_product_free (prod);
  }

  ncm_vector_free (err);
  ncm_vector_free (result);
  ncm_vector_free (s2);
  ncm_vector_free (s1);
  ncm_spectral_panels_free (sp_Q);
  ncm_spectral_panels_free (sp_P);
}

static void
test_ncm_spectral_panels_product_rows_fit (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_product/rows_fit/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*do not fit the panel sets of 2 and 2 rows*");
}

static void
test_ncm_spectral_panels_product_rows_fit_subprocess (void)
{
  NcmSpectralPanels *sp_P        = _side_panels_new (&side_P, 3, 5);
  NcmVector *result              = ncm_vector_new (2);
  NcmSpectralPanelsProduct *prod = ncm_spectral_panels_product_new (2);

  ncm_spectral_panels_product_integrate_rows (prod, sp_P, 0, NULL, sp_P, 1, NULL, result, NULL);
}

