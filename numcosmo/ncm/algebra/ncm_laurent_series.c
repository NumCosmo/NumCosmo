/***************************************************************************
 *            ncm_laurent_series.c
 *
 *  Tue Jul 8 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_laurent_series.c
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

/**
 * NcmLaurentSeries:
 *
 * Laurent polynomial with complex coefficients,
 * $$a(w) = \sum_{h = h_\mathrm{min}}^{h_\mathrm{max}} c_h w^h.$$
 *
 * Supports sums, scaling, products, the conjugation $a(w) \to \overline{a(1/\bar w)}$, which is
 * $\overline{a(w)}$ on the unit circle, evaluation, and the Jacobi-Anger reduction of
 * ncm_laurent_series_jacobi_anger_reduce(). #NcmLaurentSeriesTPS is a truncated power
 * series in a second variable $g$ whose coefficients are Laurent polynomials. They are used
 * by the weak-lensing shape series, see <a
 * href="../../theory/nc/lss/galaxy/wl_shape_marginalization_series.html">A Small-Shear Series
 * Marginalization for the Shape Likelihood</a>.
 *
 * Functions taking or returning a complex value by value are skipped by introspection;
 * their `_ptr` variants pass it by pointer and write a result into an #NcmComplex the
 * caller provides. The `_into` variants write into an existing series, resized by
 * ncm_laurent_series_reset() without shrinking its buffer, instead of allocating one;
 * their output must not alias an input. The boxed copy of both types is a new reference
 * sharing the data; ncm_laurent_series_copy() makes an independent copy.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_laurent_series.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_BOXED_TYPE (NcmLaurentSeries, ncm_laurent_series, ncm_laurent_series_ref, ncm_laurent_series_free)

/**
 * ncm_laurent_series_new:
 * @hmin: lowest power
 * @hmax: highest power, at least @hmin
 *
 * Creates a series with every coefficient of $w^{h_\mathrm{min}}, \dots, w^{h_\mathrm{max}}$ zero.
 *
 * Returns: (transfer full): a new #NcmLaurentSeries.
 */
NcmLaurentSeries *
ncm_laurent_series_new (gint hmin, gint hmax)
{
  NcmLaurentSeries *a = g_new (NcmLaurentSeries, 1);
  gint n              = hmax - hmin + 1;

  g_assert_cmpint (hmin, <=, hmax);

  a->hmin  = hmin;
  a->hmax  = hmax;
  a->c_cap = n;
  a->c     = g_new0 (NcmComplex, n);
  g_atomic_ref_count_init (&a->ref_count);

  return a;
}

/**
 * ncm_laurent_series_reset: (skip)
 * @a: a #NcmLaurentSeries
 * @hmin: new lowest power
 * @hmax: new highest power, at least @hmin
 *
 * Resizes @a to the powers @hmin to @hmax and sets every coefficient to zero. The buffer is
 * reallocated only when it is too small, so a reused series stops allocating.
 */
void
ncm_laurent_series_reset (NcmLaurentSeries *a, gint hmin, gint hmax)
{
  gint n = hmax - hmin + 1;

  g_assert_cmpint (hmin, <=, hmax);

  if (n > a->c_cap)
  {
    g_free (a->c);
    a->c     = g_new (NcmComplex, n);
    a->c_cap = n;
  }

  memset (a->c, 0, n * sizeof (NcmComplex));
  a->hmin = hmin;
  a->hmax = hmax;
}

/**
 * ncm_laurent_series_copy:
 * @a: a #NcmLaurentSeries
 *
 * Returns: (transfer full): an independent copy of @a.
 */
NcmLaurentSeries *
ncm_laurent_series_copy (const NcmLaurentSeries *a)
{
  NcmLaurentSeries *out = ncm_laurent_series_new (a->hmin, a->hmax);
  gint n                = a->hmax - a->hmin + 1;
  gint i;

  for (i = 0; i < n; i++)
    out->c[i] = a->c[i];

  return out;
}

/**
 * ncm_laurent_series_ref:
 * @a: a #NcmLaurentSeries
 *
 * Increases the reference count of @a by one; this is also the boxed copy function.
 *
 * Returns: (transfer full): @a.
 */
NcmLaurentSeries *
ncm_laurent_series_ref (NcmLaurentSeries *a)
{
  g_atomic_ref_count_inc (&a->ref_count);

  return a;
}

/**
 * ncm_laurent_series_free:
 * @a: a #NcmLaurentSeries
 *
 * Decreases the reference count of @a by one.
 */
void
ncm_laurent_series_free (NcmLaurentSeries *a)
{
  if (g_atomic_ref_count_dec (&a->ref_count))
  {
    g_free (a->c);
    g_free (a);
  }
}

/**
 * ncm_laurent_series_clear:
 * @a: a #NcmLaurentSeries
 *
 * If *@a is not %NULL, decreases its reference count by one and sets *@a to %NULL.
 */
void
ncm_laurent_series_clear (NcmLaurentSeries **a)
{
  if (*a != NULL)
  {
    ncm_laurent_series_free (*a);
    *a = NULL;
  }
}

/**
 * ncm_laurent_series_get_hmin:
 * @a: a #NcmLaurentSeries
 *
 * Returns: the lowest power $h_\mathrm{min}$.
 */
gint
ncm_laurent_series_get_hmin (const NcmLaurentSeries *a)
{
  return a->hmin;
}

/**
 * ncm_laurent_series_get_hmax:
 * @a: a #NcmLaurentSeries
 *
 * Returns: the highest power $h_\mathrm{max}$.
 */
gint
ncm_laurent_series_get_hmax (const NcmLaurentSeries *a)
{
  return a->hmax;
}

/**
 * ncm_laurent_series_new_single: (skip)
 * @h: the power
 * @val: the coefficient
 *
 * Creates the series $\mathrm{val}\,w^h$.
 *
 * Returns: (transfer full): a new #NcmLaurentSeries.
 */
NcmLaurentSeries *
ncm_laurent_series_new_single (gint h, NcmComplex val)
{
  NcmLaurentSeries *a = ncm_laurent_series_new (h, h);

  a->c[0] = val;

  return a;
}

/**
 * ncm_laurent_series_set_single_into: (skip)
 * @out: a #NcmLaurentSeries
 * @h: the power
 * @val: the coefficient
 *
 * Sets @out to $\mathrm{val}\,w^h$, as ncm_laurent_series_new_single() without allocating.
 */
void
ncm_laurent_series_set_single_into (NcmLaurentSeries *out, gint h, NcmComplex val)
{
  ncm_laurent_series_reset (out, h, h);
  out->c[0] = val;
}

/**
 * ncm_laurent_series_get: (skip)
 * @a: a #NcmLaurentSeries
 * @h: the power
 *
 * Returns: the coefficient $c_h$, zero outside $[h_\mathrm{min}, h_\mathrm{max}]$.
 */
NcmComplex
ncm_laurent_series_get (const NcmLaurentSeries *a, gint h)
{
  if ((h < a->hmin) || (h > a->hmax))
    return 0.0;

  return a->c[h - a->hmin];
}

/**
 * ncm_laurent_series_set: (skip)
 * @a: a #NcmLaurentSeries
 * @h: the power
 * @val: the coefficient
 *
 * Sets $c_h$ to @val. Aborts if $h$ is outside $[h_\mathrm{min}, h_\mathrm{max}]$.
 */
void
ncm_laurent_series_set (NcmLaurentSeries *a, gint h, NcmComplex val)
{
  g_assert_cmpint (h, >=, a->hmin);
  g_assert_cmpint (h, <=, a->hmax);

  a->c[h - a->hmin] = val;
}

/**
 * ncm_laurent_series_get_ptr:
 * @a: a #NcmLaurentSeries
 * @h: the power
 * @out: the coefficient $c_h$
 *
 * Same as ncm_laurent_series_get().
 */
void
ncm_laurent_series_get_ptr (const NcmLaurentSeries *a, gint h, NcmComplex *out)
{
  ncm_complex_set_c (out, ncm_laurent_series_get (a, h));
}

/**
 * ncm_laurent_series_set_ptr:
 * @a: a #NcmLaurentSeries
 * @h: the power
 * @val: the coefficient
 *
 * Same as ncm_laurent_series_set().
 */
void
ncm_laurent_series_set_ptr (NcmLaurentSeries *a, gint h, const NcmComplex *val)
{
  ncm_laurent_series_set (a, h, ncm_complex_c (val));
}

/**
 * ncm_laurent_series_add:
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 * @sb: factor $s_b$
 *
 * Returns: (transfer full): a new series equal to $a + s_b b$, over the union of the two ranges.
 */
NcmLaurentSeries *
ncm_laurent_series_add (const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb)
{
  gint hmin             = MIN (a->hmin, b->hmin);
  gint hmax             = MAX (a->hmax, b->hmax);
  NcmLaurentSeries *out = ncm_laurent_series_new (hmin, hmax);
  gint h;

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[h - hmin] += ncm_laurent_series_get (a, h);

  for (h = b->hmin; h <= b->hmax; h++)
    out->c[h - hmin] += sb * ncm_laurent_series_get (b, h);

  return out;
}

/**
 * ncm_laurent_series_add_into: (skip)
 * @out: a #NcmLaurentSeries
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 * @sb: factor $s_b$
 *
 * Sets @out to $a + s_b b$, as ncm_laurent_series_add() without allocating.
 */
void
ncm_laurent_series_add_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb)
{
  gint hmin = MIN (a->hmin, b->hmin);
  gint hmax = MAX (a->hmax, b->hmax);
  gint h;

  ncm_laurent_series_reset (out, hmin, hmax);

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[h - hmin] += ncm_laurent_series_get (a, h);

  for (h = b->hmin; h <= b->hmax; h++)
    out->c[h - hmin] += sb * ncm_laurent_series_get (b, h);
}

/**
 * ncm_laurent_series_scale: (skip)
 * @a: a #NcmLaurentSeries
 * @s: the factor
 *
 * Returns: (transfer full): a new series equal to $s\,a$.
 */
NcmLaurentSeries *
ncm_laurent_series_scale (const NcmLaurentSeries *a, NcmComplex s)
{
  NcmLaurentSeries *out = ncm_laurent_series_new (a->hmin, a->hmax);
  gint i, n             = a->hmax - a->hmin + 1;

  for (i = 0; i < n; i++)
    out->c[i] = s * a->c[i];

  return out;
}

/**
 * ncm_laurent_series_scale_into: (skip)
 * @out: a #NcmLaurentSeries
 * @a: a #NcmLaurentSeries
 * @s: the factor
 *
 * Sets @out to $s\,a$, as ncm_laurent_series_scale() without allocating.
 */
void
ncm_laurent_series_scale_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, NcmComplex s)
{
  gint i, n = a->hmax - a->hmin + 1;

  ncm_laurent_series_reset (out, a->hmin, a->hmax);

  for (i = 0; i < n; i++)
    out->c[i] = s * a->c[i];
}

/**
 * ncm_laurent_series_scale_ptr:
 * @a: a #NcmLaurentSeries
 * @s: the factor
 *
 * Same as ncm_laurent_series_scale().
 *
 * Returns: (transfer full): a new series equal to $s\,a$.
 */
NcmLaurentSeries *
ncm_laurent_series_scale_ptr (const NcmLaurentSeries *a, const NcmComplex *s)
{
  return ncm_laurent_series_scale (a, ncm_complex_c (s));
}

/**
 * ncm_laurent_series_conv:
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 *
 * Returns: (transfer full): a new series equal to the product $a\,b$, with powers from
 * $h_{a,\mathrm{min}} + h_{b,\mathrm{min}}$ to $h_{a,\mathrm{max}} + h_{b,\mathrm{max}}$.
 */
NcmLaurentSeries *
ncm_laurent_series_conv (const NcmLaurentSeries *a, const NcmLaurentSeries *b)
{
  gint hmin             = a->hmin + b->hmin;
  gint hmax             = a->hmax + b->hmax;
  NcmLaurentSeries *out = ncm_laurent_series_new (hmin, hmax);
  gint h1, h2;

  for (h1 = a->hmin; h1 <= a->hmax; h1++)
  {
    NcmComplex v1 = a->c[h1 - a->hmin];

    if (v1 == 0.0)
      continue;

    for (h2 = b->hmin; h2 <= b->hmax; h2++)
      out->c[(h1 + h2) - hmin] += v1 * b->c[h2 - b->hmin];
  }

  return out;
}

/**
 * ncm_laurent_series_conv_into: (skip)
 * @out: a #NcmLaurentSeries
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 *
 * Sets @out to $a\,b$, as ncm_laurent_series_conv() without allocating.
 */
void
ncm_laurent_series_conv_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b)
{
  gint hmin = a->hmin + b->hmin;
  gint hmax = a->hmax + b->hmax;
  gint h1, h2;

  ncm_laurent_series_reset (out, hmin, hmax);

  for (h1 = a->hmin; h1 <= a->hmax; h1++)
  {
    NcmComplex v1 = a->c[h1 - a->hmin];

    if (v1 == 0.0)
      continue;

    for (h2 = b->hmin; h2 <= b->hmax; h2++)
      out->c[(h1 + h2) - hmin] += v1 * b->c[h2 - b->hmin];
  }
}

/**
 * ncm_laurent_series_conj:
 * @a: a #NcmLaurentSeries
 *
 * Returns: (transfer full): a new series equal to $\overline{a(1/\bar w)} = \sum_h \bar c_h w^{-h}$,
 * which is $\overline{a(w)}$ for $|w| = 1$.
 */
NcmLaurentSeries *
ncm_laurent_series_conj (const NcmLaurentSeries *a)
{
  NcmLaurentSeries *out = ncm_laurent_series_new (-a->hmax, -a->hmin);
  gint h;

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[(-h) - out->hmin] = conj (ncm_laurent_series_get (a, h));

  return out;
}

/**
 * ncm_laurent_series_conj_into: (skip)
 * @out: a #NcmLaurentSeries
 * @a: a #NcmLaurentSeries
 *
 * Sets @out to $\overline{a(1/\bar w)}$, as ncm_laurent_series_conj() without allocating.
 */
void
ncm_laurent_series_conj_into (NcmLaurentSeries *out, const NcmLaurentSeries *a)
{
  gint h;

  ncm_laurent_series_reset (out, -a->hmax, -a->hmin);

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[(-h) - out->hmin] = conj (ncm_laurent_series_get (a, h));
}

/**
 * ncm_laurent_series_eval: (skip)
 * @a: a #NcmLaurentSeries
 * @w: the point $w$
 *
 * Returns: $a(w)$.
 */
NcmComplex
ncm_laurent_series_eval (const NcmLaurentSeries *a, NcmComplex w)
{
  NcmComplex result = 0.0;
  gint h;

  /* Horner's scheme on the powers h - hmin >= 0, then the factor w^hmin */
  for (h = a->hmax; h >= a->hmin; h--)
    result = result * w + ncm_laurent_series_get (a, h);

  return result * cpow (w, a->hmin);
}

/**
 * ncm_laurent_series_eval_ptr:
 * @a: a #NcmLaurentSeries
 * @w: the point $w$
 * @out: $a(w)$
 *
 * Same as ncm_laurent_series_eval().
 */
void
ncm_laurent_series_eval_ptr (const NcmLaurentSeries *a, const NcmComplex *w, NcmComplex *out)
{
  ncm_complex_set_c (out, ncm_laurent_series_eval (a, ncm_complex_c (w)));
}

/**
 * NcmLaurentSeriesTPS:
 *
 * Truncated power series of order $N$ in $g$ with Laurent polynomial coefficients,
 * $$\sum_{n=0}^{N} L_n(w)\,g^n \mod g^{N+1}.$$
 *
 * The order is fixed by ncm_laurent_series_tps_new(). Products are truncated at order $N$,
 * and every operation requires its operands and output to have the same order. The object is
 * meant to be kept and refilled; ncm_laurent_series_tps_conv() uses scratch owned by its
 * output, so a reused output does not allocate.
 */

/* conv_acc and conv_term: scratch of ncm_laurent_series_tps_conv() when this series is its @out */
struct _NcmLaurentSeriesTPS
{
  GPtrArray *coeffs;
  NcmLaurentSeries *conv_acc[2];
  NcmLaurentSeries *conv_term;
  gatomicrefcount ref_count;
};

G_DEFINE_BOXED_TYPE (NcmLaurentSeriesTPS, ncm_laurent_series_tps, ncm_laurent_series_tps_ref, ncm_laurent_series_tps_unref)

/**
 * ncm_laurent_series_tps_new:
 * @order: the order $N$
 *
 * Creates a series of order @order with every coefficient zero.
 *
 * Returns: (transfer full): a new #NcmLaurentSeriesTPS.
 */
NcmLaurentSeriesTPS *
ncm_laurent_series_tps_new (guint order)
{
  NcmLaurentSeriesTPS *tps = g_new (NcmLaurentSeriesTPS, 1);
  guint i;

  tps->coeffs = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_laurent_series_free);

  for (i = 0; i <= order; i++)
    g_ptr_array_add (tps->coeffs, ncm_laurent_series_new (0, 0));

  tps->conv_acc[0] = ncm_laurent_series_new (0, 0);
  tps->conv_acc[1] = ncm_laurent_series_new (0, 0);
  tps->conv_term   = ncm_laurent_series_new (0, 0);

  g_atomic_ref_count_init (&tps->ref_count);

  return tps;
}

/**
 * ncm_laurent_series_tps_ref:
 * @tps: a #NcmLaurentSeriesTPS
 *
 * Increases the reference count of @tps by one; this is also the boxed copy function, so a
 * copy shares the coefficients.
 *
 * Returns: (transfer full): @tps.
 */
NcmLaurentSeriesTPS *
ncm_laurent_series_tps_ref (NcmLaurentSeriesTPS *tps)
{
  g_atomic_ref_count_inc (&tps->ref_count);

  return tps;
}

/**
 * ncm_laurent_series_tps_unref:
 * @tps: a #NcmLaurentSeriesTPS
 *
 * Decreases the reference count of @tps by one.
 */
void
ncm_laurent_series_tps_unref (NcmLaurentSeriesTPS *tps)
{
  if (g_atomic_ref_count_dec (&tps->ref_count))
  {
    g_ptr_array_unref (tps->coeffs);
    ncm_laurent_series_free (tps->conv_acc[0]);
    ncm_laurent_series_free (tps->conv_acc[1]);
    ncm_laurent_series_free (tps->conv_term);
    g_free (tps);
  }
}

/**
 * ncm_laurent_series_tps_clear:
 * @tps: a #NcmLaurentSeriesTPS
 *
 * If *@tps is not %NULL, decreases its reference count by one and sets *@tps to %NULL.
 */
void
ncm_laurent_series_tps_clear (NcmLaurentSeriesTPS **tps)
{
  if (*tps != NULL)
  {
    ncm_laurent_series_tps_unref (*tps);
    *tps = NULL;
  }
}

/**
 * ncm_laurent_series_tps_order:
 * @tps: a #NcmLaurentSeriesTPS
 *
 * Returns: the order $N$.
 */
guint
ncm_laurent_series_tps_order (const NcmLaurentSeriesTPS *tps)
{
  return tps->coeffs->len - 1;
}

/**
 * ncm_laurent_series_tps_get:
 * @tps: a #NcmLaurentSeriesTPS
 * @n: index, at most the order
 *
 * Returns: (transfer none): the coefficient $L_n$.
 */
NcmLaurentSeries *
ncm_laurent_series_tps_get (const NcmLaurentSeriesTPS *tps, guint n)
{
  g_assert_cmpuint (n, <, tps->coeffs->len);

  return g_ptr_array_index (tps->coeffs, n);
}

/**
 * ncm_laurent_series_tps_conv:
 * @out: a #NcmLaurentSeriesTPS of the same order, not aliasing the inputs
 * @a: a #NcmLaurentSeriesTPS
 * @b: a #NcmLaurentSeriesTPS
 *
 * Sets $\mathrm{out}_m = \sum_{k=0}^{m} a_k\,b_{m-k}$ for $m \le N$, the product truncated at order $N$.
 */
void
ncm_laurent_series_tps_conv (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b)
{
  const guint order = ncm_laurent_series_tps_order (a);
  guint m;

  g_assert_cmpuint (ncm_laurent_series_tps_order (b), ==, order);
  g_assert_cmpuint (ncm_laurent_series_tps_order (out), ==, order);

  for (m = 0; m <= order; m++)
  {
    NcmLaurentSeries *acc = out->conv_acc[0];
    guint k;

    ncm_laurent_series_reset (acc, 0, 0);

    for (k = 0; k <= m; k++)
    {
      NcmLaurentSeries *acc2 = out->conv_acc[(k + 1) % 2];

      ncm_laurent_series_conv_into (out->conv_term, ncm_laurent_series_tps_get (a, k), ncm_laurent_series_tps_get (b, m - k));
      ncm_laurent_series_add_into (acc2, acc, out->conv_term, 1.0);
      acc = acc2;
    }

    ncm_laurent_series_scale_into (g_ptr_array_index (out->coeffs, m), acc, 1.0);
  }
}

/**
 * ncm_laurent_series_tps_conj:
 * @out: a #NcmLaurentSeriesTPS of the same order, not aliasing the inputs
 * @a: a #NcmLaurentSeriesTPS
 *
 * Sets each coefficient of @out to the conjugate, ncm_laurent_series_conj(), of that of @a.
 */
void
ncm_laurent_series_tps_conj (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a)
{
  const guint order = ncm_laurent_series_tps_order (a);
  guint m;

  g_assert_cmpuint (ncm_laurent_series_tps_order (out), ==, order);

  for (m = 0; m <= order; m++)
    ncm_laurent_series_conj_into (g_ptr_array_index (out->coeffs, m), ncm_laurent_series_tps_get (a, m));
}

/**
 * ncm_laurent_series_tps_add:
 * @out: a #NcmLaurentSeriesTPS of the same order, not aliasing the inputs
 * @a: a #NcmLaurentSeriesTPS
 * @b: a #NcmLaurentSeriesTPS
 * @sb: factor $s_b$
 *
 * Sets $\mathrm{out}_m = a_m + s_b b_m$.
 */
void
ncm_laurent_series_tps_add (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b, gdouble sb)
{
  const guint order = ncm_laurent_series_tps_order (a);
  guint m;

  g_assert_cmpuint (ncm_laurent_series_tps_order (b), ==, order);
  g_assert_cmpuint (ncm_laurent_series_tps_order (out), ==, order);

  for (m = 0; m <= order; m++)
    ncm_laurent_series_add_into (g_ptr_array_index (out->coeffs, m), ncm_laurent_series_tps_get (a, m), ncm_laurent_series_tps_get (b, m), sb);
}

/**
 * ncm_laurent_series_tps_scale: (skip)
 * @out: a #NcmLaurentSeriesTPS of the same order, not aliasing the inputs
 * @a: a #NcmLaurentSeriesTPS
 * @s: the factor
 *
 * Sets $\mathrm{out}_m = s\,a_m$.
 */
void
ncm_laurent_series_tps_scale (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, NcmComplex s)
{
  const guint order = ncm_laurent_series_tps_order (a);
  guint m;

  g_assert_cmpuint (ncm_laurent_series_tps_order (out), ==, order);

  for (m = 0; m <= order; m++)
    ncm_laurent_series_scale_into (g_ptr_array_index (out->coeffs, m), ncm_laurent_series_tps_get (a, m), s);
}

/**
 * ncm_laurent_series_tps_eval: (skip)
 * @tps: a #NcmLaurentSeriesTPS
 * @w: the point $w$
 * @g: the point $g$
 *
 * Returns: $\sum_{n=0}^{N} L_n(w)\,g^n$.
 */
NcmComplex
ncm_laurent_series_tps_eval (const NcmLaurentSeriesTPS *tps, NcmComplex w, NcmComplex g)
{
  const guint order = ncm_laurent_series_tps_order (tps);
  NcmComplex result = 0.0;
  gint n;

  for (n = (gint) order; n >= 0; n--)
    result = result * g + ncm_laurent_series_eval (ncm_laurent_series_tps_get (tps, (guint) n), w);

  return result;
}

/**
 * ncm_laurent_series_tps_eval_ptr:
 * @tps: a #NcmLaurentSeriesTPS
 * @w: the point $w$
 * @g: the point $g$
 * @out: the value
 *
 * Same as ncm_laurent_series_tps_eval().
 */
void
ncm_laurent_series_tps_eval_ptr (const NcmLaurentSeriesTPS *tps, const NcmComplex *w, const NcmComplex *g, NcmComplex *out)
{
  ncm_complex_set_c (out, ncm_laurent_series_tps_eval (tps, ncm_complex_c (w), ncm_complex_c (g)));
}

/**
 * ncm_laurent_series_tps_pow:
 * @out: a #NcmLaurentSeriesTPS of the same order, not aliasing the inputs
 * @a: a #NcmLaurentSeriesTPS
 * @p: the exponent
 *
 * Sets @out to $a(g)^p \mod g^{N+1}$ for a real @p. With $a = a_0 (1 + u)$ and $u(0) = 0$,
 * the coefficients of $(1 + u)^p = \sum_n c_n g^n$ follow from $(1 + u) F' = p\,u' F$,
 * $$c_0 = 1, \qquad n\,c_n = \sum_{k=1}^{n} \left[k p - (n - k)\right] u_k\,c_{n-k},$$
 * and $\mathrm{out} = a_0^p \sum_n c_n g^n$, with the principal branch of $a_0^p$.
 *
 * $L_0$ must be the constant $a_0 \ne 0$, with the single power $w^0$; aborts otherwise.
 */
void
ncm_laurent_series_tps_pow (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, gdouble p)
{
  const guint order             = ncm_laurent_series_tps_order (a);
  NcmLaurentSeriesTPS *u        = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *c        = ncm_laurent_series_tps_new (order);
  NcmLaurentSeries *fold_acc[2] = {ncm_laurent_series_new (0, 0), ncm_laurent_series_new (0, 0)};
  NcmLaurentSeries *fold_term   = ncm_laurent_series_new (0, 0);
  NcmComplex a0;
  guint n;

  g_assert_cmpuint (ncm_laurent_series_tps_order (out), ==, order);

  {
    const NcmLaurentSeries *L0 = ncm_laurent_series_tps_get (a, 0);

    g_assert_cmpint (ncm_laurent_series_get_hmin (L0), ==, 0);
    g_assert_cmpint (ncm_laurent_series_get_hmax (L0), ==, 0);

    a0 = ncm_laurent_series_get (L0, 0);
    g_assert (a0 != 0.0);
  }

  for (n = 1; n <= order; n++)
    ncm_laurent_series_scale_into (ncm_laurent_series_tps_get (u, n), ncm_laurent_series_tps_get (a, n), 1.0 / a0);

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (c, 0), 0, 1.0);

  for (n = 1; n <= order; n++)
  {
    NcmLaurentSeries *acc = fold_acc[0];
    guint k;

    ncm_laurent_series_reset (acc, 0, 0);

    for (k = 1; k <= n; k++)
    {
      NcmLaurentSeries *acc2 = fold_acc[k % 2];
      const gdouble weight   = (gdouble) k * p - (gdouble) (n - k);

      ncm_laurent_series_conv_into (fold_term, ncm_laurent_series_tps_get (u, k), ncm_laurent_series_tps_get (c, n - k));
      ncm_laurent_series_add_into (acc2, acc, fold_term, weight);
      acc = acc2;
    }

    ncm_laurent_series_scale_into (ncm_laurent_series_tps_get (c, n), acc, 1.0 / n);
  }

  ncm_laurent_series_tps_scale (out, c, cpow (a0, p));

  ncm_laurent_series_tps_unref (u);
  ncm_laurent_series_tps_unref (c);
  ncm_laurent_series_free (fold_acc[0]);
  ncm_laurent_series_free (fold_acc[1]);
  ncm_laurent_series_free (fold_term);
}

/**
 * ncm_laurent_series_jacobi_anger_reduce:
 * @cm: a #NcmLaurentSeries $c(\theta) = \sum_h c_h e^{ih\theta}$
 * @phi: the phase $\phi$
 * @Ik: (array length=n_Ik) (element-type gdouble): $e^{-z} I_k(z)$ for $k = 0, \dots, n_{I_k} - 1$
 * @n_Ik: length of @Ik, at least one
 *
 * Computes
 * $$\int_0^{2\pi} c(\theta)\,e^{z[\cos(\theta - \phi) - 1]}\,\mathrm{d}\theta
 * = 2\pi \sum_h c_h\,e^{-z} I_{|h|}(z)\,e^{ih\phi}$$
 * from the Jacobi-Anger expansion, for a real $c(\theta)$, that is $c_{-h} = \bar c_h$: only the
 * powers $h \ge 0$ are read, as $c_0 I_0 + 2 \sum_{h \ge 1} I_h\,\mathrm{Re}(c_h e^{ih\phi})$.
 * Powers $h \ge n_{I_k}$ are not included.
 *
 * Returns: the integral.
 */
gdouble
ncm_laurent_series_jacobi_anger_reduce (const NcmLaurentSeries *cm, gdouble phi, const gdouble *Ik, gint n_Ik)
{
  gdouble term;
  gint k;

  g_assert_cmpint (n_Ik, >, 0);

  term = creal (ncm_laurent_series_get (cm, 0)) * Ik[0];

  for (k = 1; k < n_Ik; k++)
  {
    NcmComplex v = ncm_laurent_series_get (cm, k);

    if (v != 0.0)
      term += 2.0 * Ik[k] * creal (v * cexp (I * k * phi));
  }

  return 2.0 * M_PI * term;
}

