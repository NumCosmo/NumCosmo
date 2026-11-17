/***************************************************************************
 *            ncm_spectral_clenshaw_private.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * The Clenshaw recurrence of a Chebyshev series, $\sum_{k = 0}^{n - 1} c_k T_k (s)$, shared
 * by the evaluators of the series in the library. Each step is written as
 * $b_k = 2 s b_{k + 1} + (c_k - b_{k + 2})$: $c_k - b_{k + 2}$ does not depend on the
 * previous step, so the chain of dependent operations is one multiply-add per step.
 * The two-dimensional series of #NcmSpectralPolarPanels use the same step
 * (ncm_spectral_polar_panels.c). ncm_spectral_chebyshev_eval() is a different method
 * (Reinsch's modification near the ends) and does not use these.
 */

#ifndef _NCM_SPECTRAL_CLENSHAW_PRIVATE_H_
#define _NCM_SPECTRAL_CLENSHAW_PRIVATE_H_

#include <glib.h>

G_BEGIN_DECLS

/*
 * ncm_spectral_clenshaw:
 * @c: the coefficients, $c_k$ at @c[k @stride]
 * @n: the number of coefficients, at least one
 * @stride: the distance between consecutive coefficients
 * @s: the point, in $[-1, 1]$
 *
 * Returns: $\sum_{k = 0}^{n - 1} c_k T_k (s)$.
 */
static inline gdouble
ncm_spectral_clenshaw (const gdouble *c, const guint n, const gsize stride, const gdouble s)
{
  const gdouble s2 = 2.0 * s;
  gdouble b1       = 0.0;
  gdouble b2       = 0.0;
  guint k;

  for (k = n - 1; k >= 1; k--)
  {
    const gdouble b0 = s2 * b1 + (c[(gsize) k * stride] - b2);

    b2 = b1;
    b1 = b0;
  }

  return s * b1 + (c[0] - b2);
}

/* G_ALWAYS_INLINE is GLib 2.74; NumCosmo needs 2.54 */
#if defined (__GNUC__) || defined (__clang__)
#define NCM_SPECTRAL_CLENSHAW_ALWAYS_INLINE __attribute__ ((__always_inline__))
#else
#define NCM_SPECTRAL_CLENSHAW_ALWAYS_INLINE
#endif

/* Most series of ncm_spectral_clenshaw_series() */
#define NCM_SPECTRAL_CLENSHAW_SERIES_MAX (4)

/*
 * ncm_spectral_clenshaw_series:
 * @c: the coefficients of @m series, interleaved: $c_k$ of series a at @c[k @m + a]
 * @n: the number of coefficients of each series, at least one
 * @m: the number of series, at most %NCM_SPECTRAL_CLENSHAW_SERIES_MAX
 * @s: the point, in $[-1, 1]$
 * @out: the @m values
 *
 * Sets @out[a] to ncm_spectral_clenshaw() of series a at @s, the @m recurrences in one
 * pass over the coefficients: they are independent, so each step has @m chains. @m is
 * meant to be a constant at the call site (the function is always inlined), so the loop
 * over the series unrolls.
 */
NCM_SPECTRAL_CLENSHAW_ALWAYS_INLINE static inline void
ncm_spectral_clenshaw_series (const gdouble *c, const guint n, const guint m, const gdouble s, gdouble *out)
{
  const gdouble s2 = 2.0 * s;
  gdouble b1[NCM_SPECTRAL_CLENSHAW_SERIES_MAX];
  gdouble b2[NCM_SPECTRAL_CLENSHAW_SERIES_MAX];
  guint a, k;

  g_assert_cmpuint (m, <=, NCM_SPECTRAL_CLENSHAW_SERIES_MAX);

  for (a = 0; a < m; a++)
  {
    b1[a] = 0.0;
    b2[a] = 0.0;
  }

  for (k = n - 1; k >= 1; k--)
  {
    const gdouble *ck = &c[(gsize) k * m];

    for (a = 0; a < m; a++)
    {
      const gdouble b0 = s2 * b1[a] + (ck[a] - b2[a]);

      b2[a] = b1[a];
      b1[a] = b0;
    }
  }

  for (a = 0; a < m; a++)
    out[a] = s * b1[a] + (c[a] - b2[a]);
}

/*
 * ncm_spectral_clenshaw_points:
 * @c: the @n coefficients, contiguous
 * @n: the number of coefficients, at least one
 * @m: the number of points
 * @s: the @m points, in $[-1, 1]$
 * @b1: scratch of @m entries
 * @b2: scratch of @m entries
 * @out: the @m values, may be @s
 *
 * Sets @out[j] to ncm_spectral_clenshaw() at @s[j]. One pass over the coefficients
 * carries the @m recurrences, which are independent, so the loop over the points
 * vectorizes and each coefficient is loaded once.
 */
static inline void
ncm_spectral_clenshaw_points (const gdouble *c, const guint n, const guint m, const gdouble *s,
                              gdouble *b1, gdouble *b2, gdouble *out)
{
  guint j, k;

  #pragma omp simd

  for (j = 0; j < m; j++)
  {
    b1[j] = 0.0;
    b2[j] = 0.0;
  }

  for (k = n - 1; k >= 1; k--)
  {
    const gdouble ck = c[k];

    #pragma omp simd

    for (j = 0; j < m; j++)
    {
      const gdouble b0 = 2.0 * s[j] * b1[j] + (ck - b2[j]);

      b2[j] = b1[j];
      b1[j] = b0;
    }
  }

  #pragma omp simd

  for (j = 0; j < m; j++)
    out[j] = s[j] * b1[j] + (c[0] - b2[j]);
}

G_END_DECLS

#endif /* _NCM_SPECTRAL_CLENSHAW_PRIVATE_H_ */

