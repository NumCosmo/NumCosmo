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
 * von-Mises-shaped kernel to scaled modified Bessel functions. Generic
 * complex-analysis machinery with no physics content of its own, used and
 * independently tested by `nc_galaxy_shape_factor_series_lensed.c` (see
 * docs/theory/wl_shape_marginalization_series.qmd for the physics context
 * and derivation).
 *
 * Two parallel calling conventions for every operation, matching
 * nc_wl_ellipticity.h's own introspectable/native dual interface:
 *
 * - an introspectable, #NcmComplex-based API (the plain function names),
 *   usable and testable from Python;
 * - a native, `complex double`-based secondary interface (`_c`-suffixed or,
 *   for functions with no complex-valued parameters at all, just direct
 *   struct-field access), declared behind `#ifndef NUMCOSMO_GIR_SCAN` since
 *   C99 complex types are not introspectable -- what the hot loop in
 *   nc_galaxy_shape_factor_series_lensed.c actually uses.
 *
 * Every #NcmLaurentSeries is independently heap-allocated (no shared arena):
 * callers own whatever they create and must free it themselves (or track a
 * short-lived batch, e.g. a #GPtrArray with ncm_laurent_series_free() as its
 * free function, freed in one sweep at the end of a computation).
 *
 * For hot loops that would otherwise allocate and free many of these per
 * call, ncm_laurent_series_reset() plus the `_into`-suffixed counterparts
 * of add/conv/scale_c/conj/new_single write into a caller-supplied,
 * already-allocated #NcmLaurentSeries instead of allocating a new one
 * (grow-only, so a reused instance never reallocates once it reaches its
 * steady-state size) -- see nc_galaxy_shape_factor_series_lensed.c for the
 * intended usage (a per-instance pool of reusable series).
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

