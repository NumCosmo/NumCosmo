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
 * Complex Laurent-polynomial arithmetic (add/scale/convolve/conjugate) plus
 * the Jacobi-Anger reduction of a truncated Fourier series against a
 * von-Mises-shaped kernel to scaled modified Bessel functions, and
 * #NcmLaurentSeriesTPS, a truncated power series (in a second, formal
 * variable) whose coefficients are #NcmLaurentSeries. Generic
 * complex-analysis machinery with no physics content of its own, used and
 * independently tested by `nc_wl_ellipticity_series.c` (see
 * docs/theory/wl_shape_marginalization_series.qmd for the physics context
 * and derivation).
 *
 * Two parallel calling conventions for every #NcmLaurentSeries operation,
 * matching nc_wl_ellipticity.h's own introspectable/native dual interface:
 *
 * - an introspectable, #NcmComplex-based API (the plain function names),
 *   usable and testable from Python;
 * - a native, `complex double`-based secondary interface (`_c`-suffixed or,
 *   for functions with no complex-valued parameters at all, just direct
 *   struct-field access), declared behind `#ifndef NUMCOSMO_GIR_SCAN` since
 *   C99 complex types are not introspectable -- what the hot loop in
 *   nc_wl_ellipticity_series.c actually uses.
 *
 * Every #NcmLaurentSeries is independently heap-allocated: callers own
 * whatever they create and must free it themselves (or track a short-lived
 * batch, e.g. a #GPtrArray with ncm_laurent_series_free() as its free
 * function, freed in one sweep at the end of a computation).
 *
 * For hot loops that would otherwise allocate and free many of these per
 * call, ncm_laurent_series_reset() plus the `_into`-suffixed counterparts
 * of add/conv/scale_c/conj/new_single write into a caller-supplied,
 * already-allocated #NcmLaurentSeries instead of allocating a new one
 * (grow-only, so a reused instance never reallocates once it reaches its
 * steady-state size) -- #NcmLaurentSeriesTPS's own coefficients and
 * ncm_laurent_series_tps_conv()'s private scratch are exactly such
 * long-lived, reused instances.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_laurent_series.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#include <math.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_BOXED_TYPE (NcmLaurentSeries, ncm_laurent_series, ncm_laurent_series_copy, ncm_laurent_series_free)

/**
 * ncm_laurent_series_new:
 * @hmin: lowest harmonic
 * @hmax: highest harmonic
 *
 * Creates a new zero-initialized #NcmLaurentSeries with harmonics
 * $h\in[hmin,hmax]$.
 *
 * Returns: (transfer full): a new #NcmLaurentSeries
 */
NcmLaurentSeries *
ncm_laurent_series_new (gint hmin, gint hmax)
{
  NcmLaurentSeries *a = g_new (NcmLaurentSeries, 1);
  gint n              = hmax - hmin + 1;

  a->hmin  = hmin;
  a->hmax  = hmax;
  a->c_cap = n;
  a->c     = g_new0 (complex double, n);

  return a;
}

/**
 * ncm_laurent_series_reset: (skip)
 * @a: a #NcmLaurentSeries
 * @hmin: new lowest harmonic
 * @hmax: new highest harmonic
 *
 * Resizes @a in place to harmonics $h\in[hmin,hmax]$, zeroing every
 * coefficient. Grow-only: reuses @a's existing buffer (no allocation) when
 * it is already big enough, reallocates bigger otherwise -- never shrinks
 * the underlying allocation. Lets a hot loop reuse one #NcmLaurentSeries
 * across many calls instead of allocating fresh every time (see
 * nc_galaxy_shape_factor_series_lensed.c).
 */
void
ncm_laurent_series_reset (NcmLaurentSeries *a, gint hmin, gint hmax)
{
  gint n = hmax - hmin + 1;

  if (n > a->c_cap)
  {
    g_free (a->c);
    a->c     = g_new (complex double, n);
    a->c_cap = n;
  }

  memset (a->c, 0, n * sizeof (complex double));
  a->hmin = hmin;
  a->hmax = hmax;
}

/**
 * ncm_laurent_series_copy:
 * @a: a #NcmLaurentSeries
 *
 * Returns: (transfer full): a new, independent copy of @a
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
 * ncm_laurent_series_free:
 * @a: a #NcmLaurentSeries
 *
 * Frees @a.
 */
void
ncm_laurent_series_free (NcmLaurentSeries *a)
{
  g_free (a->c);
  g_free (a);
}

/**
 * ncm_laurent_series_clear:
 * @a: a #NcmLaurentSeries
 *
 * Frees *@a and sets it to %NULL.
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
 * ncm_laurent_series_hmin:
 * @a: a #NcmLaurentSeries
 *
 * Returns: @a's lowest harmonic
 */
gint
ncm_laurent_series_hmin (const NcmLaurentSeries *a)
{
  return a->hmin;
}

/**
 * ncm_laurent_series_hmax:
 * @a: a #NcmLaurentSeries
 *
 * Returns: @a's highest harmonic
 */
gint
ncm_laurent_series_hmax (const NcmLaurentSeries *a)
{
  return a->hmax;
}

/**
 * ncm_laurent_series_new_single: (skip)
 * @h: the single nonzero harmonic
 * @val: its coefficient
 *
 * Creates a new #NcmLaurentSeries with a single nonzero term $val\cdot w^h$.
 *
 * Returns: (transfer full): the new series
 */
NcmLaurentSeries *
ncm_laurent_series_new_single (gint h, complex double val)
{
  NcmLaurentSeries *a = ncm_laurent_series_new (h, h);

  a->c[0] = val;

  return a;
}

/**
 * ncm_laurent_series_set_single_into: (skip)
 * @out: a #NcmLaurentSeries, resized in place
 * @h: the single nonzero harmonic
 * @val: its coefficient
 *
 * Reuse counterpart of ncm_laurent_series_new_single(): resets @out to
 * $[h,h]$ and sets its one coefficient, instead of allocating.
 */
void
ncm_laurent_series_set_single_into (NcmLaurentSeries *out, gint h, complex double val)
{
  ncm_laurent_series_reset (out, h, h);
  out->c[0] = val;
}

/**
 * ncm_laurent_series_get_c:
 * @a: a #NcmLaurentSeries
 * @h: the harmonic to query
 *
 * Returns: the coefficient of $w^h$ in @a, or 0 if $h$ is outside @a's range
 */
complex double
ncm_laurent_series_get_c (const NcmLaurentSeries *a, gint h)
{
  if ((h < a->hmin) || (h > a->hmax))
    return 0.0;

  return a->c[h - a->hmin];
}

/**
 * ncm_laurent_series_set_c: (skip)
 * @a: a #NcmLaurentSeries
 * @h: the harmonic to set, must already be within @a's $[hmin,hmax]$ range
 * @val: the new coefficient of $w^h$
 */
void
ncm_laurent_series_set_c (NcmLaurentSeries *a, gint h, complex double val)
{
  g_assert_cmpint (h, >=, a->hmin);
  g_assert_cmpint (h, <=, a->hmax);

  a->c[h - a->hmin] = val;
}

/**
 * ncm_laurent_series_get:
 * @a: a #NcmLaurentSeries
 * @h: the harmonic to query
 * @out: output #NcmComplex, the coefficient of $w^h$, or 0 if $h$ is
 * outside @a's range
 */
void
ncm_laurent_series_get (const NcmLaurentSeries *a, gint h, NcmComplex *out)
{
  ncm_complex_set_c (out, ncm_laurent_series_get_c (a, h));
}

/**
 * ncm_laurent_series_set:
 * @a: a #NcmLaurentSeries
 * @h: the harmonic to set, must already be within @a's $[hmin,hmax]$ range
 * @val: the new coefficient of $w^h$
 */
void
ncm_laurent_series_set (NcmLaurentSeries *a, gint h, const NcmComplex *val)
{
  ncm_laurent_series_set_c (a, h, ncm_complex_c (val));
}

/**
 * ncm_laurent_series_add:
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 * @sb: scale factor applied to @b
 *
 * Returns: (transfer full): a new series equal to $a+sb\cdot b$
 */
NcmLaurentSeries *
ncm_laurent_series_add (const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb)
{
  gint hmin             = MIN (a->hmin, b->hmin);
  gint hmax             = MAX (a->hmax, b->hmax);
  NcmLaurentSeries *out = ncm_laurent_series_new (hmin, hmax);
  gint h;

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[h - hmin] += ncm_laurent_series_get_c (a, h);

  for (h = b->hmin; h <= b->hmax; h++)
    out->c[h - hmin] += sb * ncm_laurent_series_get_c (b, h);

  return out;
}

/**
 * ncm_laurent_series_add_into: (skip)
 * @out: a #NcmLaurentSeries, resized in place, must not alias @a or @b
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 * @sb: scale factor applied to @b
 *
 * Reuse counterpart of ncm_laurent_series_add(): writes $a+sb\cdot b$ into
 * @out instead of allocating a new series.
 */
void
ncm_laurent_series_add_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb)
{
  gint hmin = MIN (a->hmin, b->hmin);
  gint hmax = MAX (a->hmax, b->hmax);
  gint h;

  ncm_laurent_series_reset (out, hmin, hmax);

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[h - hmin] += ncm_laurent_series_get_c (a, h);

  for (h = b->hmin; h <= b->hmax; h++)
    out->c[h - hmin] += sb * ncm_laurent_series_get_c (b, h);
}

/**
 * ncm_laurent_series_scale_c: (skip)
 * @a: a #NcmLaurentSeries
 * @s: scale factor
 *
 * Returns: (transfer full): a new series equal to $s\cdot a$
 */
NcmLaurentSeries *
ncm_laurent_series_scale_c (const NcmLaurentSeries *a, complex double s)
{
  NcmLaurentSeries *out = ncm_laurent_series_new (a->hmin, a->hmax);
  gint i, n             = a->hmax - a->hmin + 1;

  for (i = 0; i < n; i++)
    out->c[i] = s * a->c[i];

  return out;
}

/**
 * ncm_laurent_series_scale_c_into: (skip)
 * @out: a #NcmLaurentSeries, resized in place, must not alias @a
 * @a: a #NcmLaurentSeries
 * @s: scale factor
 *
 * Reuse counterpart of ncm_laurent_series_scale_c(): writes $s\cdot a$ into
 * @out instead of allocating a new series.
 */
void
ncm_laurent_series_scale_c_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, complex double s)
{
  gint i, n = a->hmax - a->hmin + 1;

  ncm_laurent_series_reset (out, a->hmin, a->hmax);

  for (i = 0; i < n; i++)
    out->c[i] = s * a->c[i];
}

/**
 * ncm_laurent_series_scale:
 * @a: a #NcmLaurentSeries
 * @s: scale factor
 *
 * Returns: (transfer full): a new series equal to $s\cdot a$
 */
NcmLaurentSeries *
ncm_laurent_series_scale (const NcmLaurentSeries *a, const NcmComplex *s)
{
  return ncm_laurent_series_scale_c (a, ncm_complex_c (s));
}

/**
 * ncm_laurent_series_conv:
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 *
 * Returns: (transfer full): the new series equal to the Laurent-polynomial
 * product $a\cdot b$
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
    complex double v1 = a->c[h1 - a->hmin];

    if (v1 == 0.0)
      continue;

    for (h2 = b->hmin; h2 <= b->hmax; h2++)
      out->c[(h1 + h2) - hmin] += v1 * b->c[h2 - b->hmin];
  }

  return out;
}

/**
 * ncm_laurent_series_conv_into: (skip)
 * @out: a #NcmLaurentSeries, resized in place, must not alias @a or @b
 * @a: a #NcmLaurentSeries
 * @b: a #NcmLaurentSeries
 *
 * Reuse counterpart of ncm_laurent_series_conv(): writes the
 * Laurent-polynomial product $a\cdot b$ into @out instead of allocating a
 * new series.
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
    complex double v1 = a->c[h1 - a->hmin];

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
 * Returns: (transfer full): the new series equal to $\overline{a(1/w)}$,
 * i.e. the harmonic-by-harmonic complex conjugate (harmonic $h$'s
 * coefficient becomes harmonic $-h$'s, conjugated)
 */
NcmLaurentSeries *
ncm_laurent_series_conj (const NcmLaurentSeries *a)
{
  NcmLaurentSeries *out = ncm_laurent_series_new (-a->hmax, -a->hmin);
  gint h;

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[(-h) - out->hmin] = conj (ncm_laurent_series_get_c (a, h));

  return out;
}

/**
 * ncm_laurent_series_conj_into: (skip)
 * @out: a #NcmLaurentSeries, resized in place, must not alias @a
 * @a: a #NcmLaurentSeries
 *
 * Reuse counterpart of ncm_laurent_series_conj(): writes
 * $\overline{a(1/w)}$ into @out instead of allocating a new series.
 */
void
ncm_laurent_series_conj_into (NcmLaurentSeries *out, const NcmLaurentSeries *a)
{
  gint h;

  ncm_laurent_series_reset (out, -a->hmax, -a->hmin);

  for (h = a->hmin; h <= a->hmax; h++)
    out->c[(-h) - out->hmin] = conj (ncm_laurent_series_get_c (a, h));
}

/**
 * ncm_laurent_series_eval_c: (skip)
 * @a: a #NcmLaurentSeries
 * @w: the point to evaluate at
 *
 * Returns: $a(w)=\sum_{h=hmin}^{hmax} c_h w^h$
 */
complex double
ncm_laurent_series_eval_c (const NcmLaurentSeries *a, complex double w)
{
  complex double result = 0.0;
  gint h;

  /* Horner over the *shifted*, all-non-negative exponents h-hmin (a plain
   * Laurent range can start negative, where Horner's scheme doesn't apply
   * directly); multiplying the shifted result by w^hmin at the end recovers
   * sum_h c_h w^h exactly. */
  for (h = a->hmax; h >= a->hmin; h--)
    result = result * w + ncm_laurent_series_get_c (a, h);

  return result * cpow (w, a->hmin);
}

/**
 * ncm_laurent_series_eval:
 * @a: a #NcmLaurentSeries
 * @w: the point to evaluate at
 * @out: (out): $a(w)$
 */
void
ncm_laurent_series_eval (const NcmLaurentSeries *a, const NcmComplex *w, NcmComplex *out)
{
  ncm_complex_set_c (out, ncm_laurent_series_eval_c (a, ncm_complex_c (w)));
}

/* @coeffs: owned directly, length order+1
 * @conv_acc: private ping-pong scratch for ncm_laurent_series_tps_conv()
 * when this TPS is the @out of a conv() call -- see that function's own
 * comment for why exactly two suffice regardless of order
 * @conv_term: private scratch for ncm_laurent_series_tps_conv() */
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
 * @order: the truncation order $N$
 *
 * Creates a new #NcmLaurentSeriesTPS of order @order ($N+1$ zero
 * coefficients, each an independently owned #NcmLaurentSeries) plus its own
 * private conv() scratch. Order is immutable for the object's whole life --
 * meant to be constructed once and refilled many times (see
 * nc_wl_ellipticity_series.c), not a short-lived per-call temporary.
 *
 * Returns: (transfer full): a new #NcmLaurentSeriesTPS
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
 * Increases the reference count of @tps by one. This is the boxed type's
 * "copy" function: it shares the same underlying storage rather than
 * deep-copying it (matching #NcGalaxySDShapeData and siblings), so later
 * mutations through any reference (e.g. the owning evaluator's own compute
 * step) are visible through every other outstanding reference.
 *
 * Returns: (transfer full): @tps
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
 * Decreases the reference count of @tps by one, freeing it once the count
 * reaches zero.
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
 * Unrefs *@tps and sets it to %NULL.
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
 * Returns: @tps's order $N$ (number of coefficients minus one)
 */
guint
ncm_laurent_series_tps_order (const NcmLaurentSeriesTPS *tps)
{
  return tps->coeffs->len - 1;
}

/**
 * ncm_laurent_series_tps_get:
 * @tps: a #NcmLaurentSeriesTPS
 * @n: the coefficient index, $0\le n\le$ @tps's order
 *
 * Returns: (transfer none): @tps's coefficient of $g^n$
 */
NcmLaurentSeries *
ncm_laurent_series_tps_get (const NcmLaurentSeriesTPS *tps, guint n)
{
  g_assert_cmpuint (n, <, tps->coeffs->len);

  return g_ptr_array_index (tps->coeffs, n);
}

/**
 * ncm_laurent_series_tps_conv: (skip)
 * @out: a #NcmLaurentSeriesTPS, same order as @a and @b, must not alias either
 * @a: a #NcmLaurentSeriesTPS
 * @b: a #NcmLaurentSeriesTPS, same order as @a
 *
 * Truncated Cauchy product $out_m=\sum_{k=0}^m a_k b_{m-k}$ for
 * $m=0..N$ -- truncated at the shared order $N$, not extended to
 * $\deg(a)+\deg(b)$ (this is a truncated power series, not a polynomial
 * ring element). The inner fold only ever needs two non-aliasing
 * accumulator buffers (ncm_laurent_series_add_into() forbids @out aliasing
 * an input, so the running sum alternates between them: at step $k$, the
 * accumulator from step $k-2$ is already fully consumed and free to reuse
 * for step $k$'s result) plus one transient product buffer -- three fixed
 * buffers regardless of order, drawn from @out's own private
 * @conv_acc/@conv_term rather than any external pool.
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

    ncm_laurent_series_scale_c_into (g_ptr_array_index (out->coeffs, m), acc, 1.0);
  }
}

/**
 * ncm_laurent_series_tps_conj: (skip)
 * @out: a #NcmLaurentSeriesTPS, same order as @a, must not alias @a
 * @a: a #NcmLaurentSeriesTPS
 *
 * Harmonic-by-harmonic conjugate of every coefficient of @a (see
 * ncm_laurent_series_conj_into()), term by term in $g$. Writes straight
 * into @out's existing slots, no scratch #NcmLaurentSeries needed.
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
 * ncm_laurent_series_tps_add: (skip)
 * @out: a #NcmLaurentSeriesTPS, same order as @a and @b, must not alias either
 * @a: a #NcmLaurentSeriesTPS
 * @b: a #NcmLaurentSeriesTPS, same order as @a
 * @sb: scale factor applied to @b
 *
 * Term-by-term $out_m=a_m+sb\cdot b_m$. Writes straight into @out's
 * existing slots, no scratch #NcmLaurentSeries needed.
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
 * @out: a #NcmLaurentSeriesTPS, same order as @a, must not alias @a
 * @a: a #NcmLaurentSeriesTPS
 * @s: scale factor
 *
 * Term-by-term $out_m=s\cdot a_m$. Writes straight into @out's existing
 * slots, no scratch #NcmLaurentSeries needed.
 */
void
ncm_laurent_series_tps_scale (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, complex double s)
{
  const guint order = ncm_laurent_series_tps_order (a);
  guint m;

  g_assert_cmpuint (ncm_laurent_series_tps_order (out), ==, order);

  for (m = 0; m <= order; m++)
    ncm_laurent_series_scale_c_into (g_ptr_array_index (out->coeffs, m), ncm_laurent_series_tps_get (a, m), s);
}

/**
 * ncm_laurent_series_tps_eval_c: (skip)
 * @tps: a #NcmLaurentSeriesTPS
 * @w: the coefficients' own formal variable, evaluated at this point
 * @g: the truncation variable, evaluated at this point
 *
 * Returns: $\mathrm{tps}(w,g)=\sum_{n=0}^N L_n(w)\,g^n$
 */
complex double
ncm_laurent_series_tps_eval_c (const NcmLaurentSeriesTPS *tps, complex double w, complex double g)
{
  const guint order     = ncm_laurent_series_tps_order (tps);
  complex double result = 0.0;
  gint n;

  for (n = (gint) order; n >= 0; n--)
    result = result * g + ncm_laurent_series_eval_c (ncm_laurent_series_tps_get (tps, (guint) n), w);

  return result;
}

/**
 * ncm_laurent_series_tps_eval:
 * @tps: a #NcmLaurentSeriesTPS
 * @w: the coefficients' own formal variable, evaluated at this point
 * @g: the truncation variable, evaluated at this point
 * @out: (out): $\mathrm{tps}(w,g)$
 */
void
ncm_laurent_series_tps_eval (const NcmLaurentSeriesTPS *tps, const NcmComplex *w, const NcmComplex *g, NcmComplex *out)
{
  ncm_complex_set_c (out, ncm_laurent_series_tps_eval_c (tps, ncm_complex_c (w), ncm_complex_c (g)));
}

/**
 * ncm_laurent_series_jacobi_anger_reduce:
 * @cm: a #NcmLaurentSeries
 * @phi: the kernel's own phase offset
 * @Ik: (array length=n_Ik) (element-type gdouble): scaled Bessel values
 * $\exp(-z)I_h(z)$ for $h=0..n\_Ik-1$
 * @n_Ik: length of @Ik
 *
 * Exact reduction of
 * $\int_0^{2\pi} cm(\theta)\exp(z(\cos(\theta-\phi)-1))\,d\theta$
 * via the Jacobi-Anger identity, given @cm's own harmonic content and
 * precomputed scaled Bessel values -- no numerical quadrature over $\theta$
 * at all. Verified against direct numerical theta-integration
 * (test_ncm_laurent_series.c); note this includes the overall $2\pi$ factor
 * from the $\theta$-integral itself, so the return value is meaningful on
 * its own without relying on a caller-side normalization to supply it.
 *
 * Returns: the reduced value
 */
gdouble
ncm_laurent_series_jacobi_anger_reduce (const NcmLaurentSeries *cm, gdouble phi, const gdouble *Ik, gint n_Ik)
{
  gdouble term = creal (ncm_laurent_series_get_c (cm, 0)) * Ik[0];
  gint k;

  for (k = 1; k < n_Ik; k++)
  {
    complex double v = ncm_laurent_series_get_c (cm, k);

    if (v != 0.0)
      term += 2.0 * Ik[k] * creal (v * cexp (I * k * phi));
  }

  return 2.0 * M_PI * term;
}

