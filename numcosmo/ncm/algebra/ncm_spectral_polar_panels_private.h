/***************************************************************************
 *            ncm_spectral_polar_panels_private.h
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _NCM_SPECTRAL_POLAR_PANELS_PRIVATE_H_
#define _NCM_SPECTRAL_POLAR_PANELS_PRIVATE_H_

#include <glib.h>
#include "ncm/algebra/ncm_spectral_polar_panels.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

/*
 * The layout of NcmSpectralPolarPanels and the evaluation of a panel, shared by
 * ncm_spectral_polar_panels.c and ncm_spectral_polar_panels_array.c, which evaluates the
 * panels it combines with the same code. Not part of the public API.
 */

/* Largest level of a panel in either direction */
#define NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX (20)

/* G_ALWAYS_INLINE is GLib 2.74; NumCosmo needs 2.54 */
#if defined (__GNUC__) || defined (__clang__)
#define NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE __attribute__ ((__always_inline__))
#else
#define NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE
#endif

/* A panel, its numbers of radial and angular coefficients and the position of its
 * coefficients in the coefficient array */
typedef struct _NcmSpectralPolarPanelsPanel
{
  gdouble ra;
  gdouble rb;
  gdouble ta;
  gdouble tb;
  guint n_r;
  guint n_theta;
  gsize offset;
} NcmSpectralPolarPanelsPanel;

struct _NcmSpectralPolarPanels
{
  /*< private >*/
  GObject parent_instance;
  GArray *r_edges;       /* gdouble, the radial edges of the initial panels */
  GArray *theta_edges;   /* gdouble, their angular edges */
  GArray *first;         /* guint, n_r n_theta + 1 entries: initial panel c holds the panels
                          * first[c], ..., first[c + 1] - 1 */
  GArray *panels;        /* NcmSpectralPolarPanelsPanel */
  GArray *coeffs;        /* gdouble, the coefficients of every panel one after the other */
  GVariant *init_panels; /* the construct-only property, consumed by constructed() */
};

/*
 * The index i of the interval [e_i, e_{i+1}) of the len edges e that holds x, by
 * bisection. A point below the first edge gives the first interval, and a point on or above
 * the last edge the last one.
 */
static inline guint
_ncm_spectral_polar_panels_edge_index (const gdouble *e, guint len, gdouble x)
{
  guint lo = 0;
  guint hi = len - 2;

  while (lo < hi)
  {
    const guint mid = (lo + hi + 1) / 2;

    if (x >= e[mid])
      lo = mid;
    else
      hi = mid - 1;
  }

  return lo;
}

/* ncm_spectral_polar_panels_find_panel() on the arrays of a set or of a pack */
NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline guint
_ncm_spectral_polar_panels_find_panel_in (const gdouble *r_edges, guint n_r_edges, const gdouble *theta_edges, guint n_theta_edges,
                                          const guint *first, const NcmSpectralPolarPanelsPanel *panels, gdouble r, gdouble theta)
{
  const guint i    = _ncm_spectral_polar_panels_edge_index (r_edges, n_r_edges, r);
  const guint j    = _ncm_spectral_polar_panels_edge_index (theta_edges, n_theta_edges, theta);
  const guint cell = i * (n_theta_edges - 1) + j;
  const guint p0   = first[cell];
  const guint p1   = first[cell + 1];

  g_assert_cmpuint (p0, <, p1);

  if (p0 + 1 == p1)
    return p0;

  {
    const gdouble cell_ra = r_edges[i];
    const gdouble cell_rb = r_edges[i + 1];
    const gdouble cell_ta = theta_edges[j];
    const gdouble cell_tb = theta_edges[j + 1];
    guint p;

    /* A panel on a side of the initial panel also takes the points beyond that side */
    for (p = p0; p + 1 < p1; p++)
    {
      const NcmSpectralPolarPanelsPanel *panel = &panels[p];

      if (((r >= panel->ra) || (panel->ra == cell_ra)) && ((r < panel->rb) || (panel->rb == cell_rb)) &&
          ((theta >= panel->ta) || (panel->ta == cell_ta)) && ((theta < panel->tb) || (panel->tb == cell_tb)))
        return p;
    }

    return p;
  }
}

/* The series of panel with the coefficients c at (r, theta), see ncm_spectral_polar_panels_eval_on() */
static inline gdouble
_ncm_spectral_polar_panels_eval_on_c (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, gdouble r, gdouble theta)
{
  const guint n     = panel->n_r - 1;
  const guint m     = panel->n_theta - 1;
  const gsize ld    = panel->n_r;
  const gdouble sr  = (2.0 * r - panel->ra - panel->rb) / (panel->rb - panel->ra);
  const gdouble st  = (2.0 * theta - panel->ta - panel->tb) / (panel->tb - panel->ta);
  const gdouble sr2 = 2.0 * sr;
  const gdouble st2 = 2.0 * st;
  gdouble B1        = 0.0;
  gdouble B2        = 0.0;
  gint l            = (gint) m;
  guint k;

  /*
   * Each step is written as 2 s b_1 + (c_k - b_2): c_k - b_2 does not depend
   * on the previous step, so the chain of dependent operations is one
   * multiply-add per step instead of a multiply-add followed by an add. The
   * radial recurrences of different angular indices are independent, so
   * they run four at a time, from the highest index down, and their results
   * enter the angular recurrence in the same order.
   */
  for ( ; l >= 3; l -= 4)
  {
    const gdouble *c0 = c + (gsize) l * ld;
    const gdouble *c1 = c0 - ld;
    const gdouble *c2 = c1 - ld;
    const gdouble *c3 = c2 - ld;
    gdouble a1 = 0.0, a2 = 0.0, b1 = 0.0, b2 = 0.0, d1 = 0.0, d2 = 0.0, e1 = 0.0, e2 = 0.0;
    gdouble R[4];
    guint q;

    for (k = n; k > 0; k--)
    {
      const gdouble a0 = sr2 * a1 + (c0[k] - a2);
      const gdouble b0 = sr2 * b1 + (c1[k] - b2);
      const gdouble d0 = sr2 * d1 + (c2[k] - d2);
      const gdouble e0 = sr2 * e1 + (c3[k] - e2);

      a2 = a1;
      a1 = a0;
      b2 = b1;
      b1 = b0;
      d2 = d1;
      d1 = d0;
      e2 = e1;
      e1 = e0;
    }

    R[0] = sr * a1 + (c0[0] - a2);
    R[1] = sr * b1 + (c1[0] - b2);
    R[2] = sr * d1 + (c2[0] - d2);
    R[3] = sr * e1 + (c3[0] - e2);

    for (q = 0; q < 4; q++)
    {
      if (l - (gint) q > 0)
      {
        const gdouble B0 = st2 * B1 + (R[q] - B2);

        B2 = B1;
        B1 = B0;
      }
      else
      {
        return st * B1 + (R[q] - B2);
      }
    }
  }

  for ( ; l >= 0; l--)
  {
    const gdouble *row = c + (gsize) l * ld;
    gdouble b1         = 0.0;
    gdouble b2         = 0.0;
    gdouble R;

    for (k = n; k > 0; k--)
    {
      const gdouble b0 = sr2 * b1 + (row[k] - b2);

      b2 = b1;
      b1 = b0;
    }

    R = sr * b1 + (row[0] - b2);

    if (l > 0)
    {
      const gdouble B0 = st2 * B1 + (R - B2);

      B2 = B1;
      B1 = B0;
    }
    else
    {
      return st * B1 + (R - B2);
    }
  }

  g_assert_not_reached (); /* LCOV_EXCL_LINE */

  return 0.0; /* LCOV_EXCL_LINE */
}

/*
 * Lanes of ncm_spectral_polar_panels_eval_on_array(): the series of one panel at W points
 * at once, W a constant at each call site (the function is always inlined), so that the
 * loops over the points vectorize. Each coefficient is loaded once and used by every
 * point, and two angular rows run together, so that each point has two independent
 * radial recurrences. Every point goes through the operations of
 * ncm_spectral_polar_panels_eval_on() in the same order, so the two differ only where the
 * compiler contracts them into fma differently (ncm_spectral_polar_panels_eval_on_array()).
 * r and theta are read before out is written, so out may be r or theta.
 */
#define NCM_SPECTRAL_POLAR_PANELS_LANES_MAX 8

NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline void
_ncm_spectral_polar_panels_eval_lanes (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, const guint W,
                                       const gdouble *r, const gdouble *theta, gdouble *out)
{
  const guint n    = panel->n_r - 1;
  const gsize ld   = panel->n_r;
  const gdouble ra = panel->ra;
  const gdouble rb = panel->rb;
  const gdouble ta = panel->ta;
  const gdouble tb = panel->tb;
  gdouble sr[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], st[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  gdouble sr2[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], st2[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  gdouble B1[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], B2[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  gdouble a1[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], a2[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  gdouble b1[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], b2[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  gint l = (gint) panel->n_theta - 1;
  guint k, w;

  #pragma omp simd

  for (w = 0; w < W; w++)
  {
    /* As in ncm_spectral_polar_panels_eval_on(), operation by operation */
    sr[w]  = (2.0 * r[w] - ra - rb) / (rb - ra);
    st[w]  = (2.0 * theta[w] - ta - tb) / (tb - ta);
    sr2[w] = 2.0 * sr[w];
    st2[w] = 2.0 * st[w];
    B1[w]  = 0.0;
    B2[w]  = 0.0;
  }

  /* Two rows, l and l - 1, while both enter the angular recurrence before its last step */
  for ( ; l >= 2; l -= 2)
  {
    const gdouble *c0 = c + (gsize) l * ld;
    const gdouble *c1 = c0 - ld;

    #pragma omp simd

    for (w = 0; w < W; w++)
    {
      a1[w] = 0.0;
      a2[w] = 0.0;
      b1[w] = 0.0;
      b2[w] = 0.0;
    }

    for (k = n; k > 0; k--)
    {
      const gdouble x0 = c0[k];
      const gdouble x1 = c1[k];

      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble a0 = sr2[w] * a1[w] + (x0 - a2[w]);
        const gdouble b0 = sr2[w] * b1[w] + (x1 - b2[w]);

        a2[w] = a1[w];
        a1[w] = a0;
        b2[w] = b1[w];
        b1[w] = b0;
      }
    }

    #pragma omp simd

    for (w = 0; w < W; w++)
    {
      const gdouble R0  = sr[w] * a1[w] + (c0[0] - a2[w]);
      const gdouble R1  = sr[w] * b1[w] + (c1[0] - b2[w]);
      const gdouble B0  = st2[w] * B1[w] + (R0 - B2[w]);
      const gdouble B0m = st2[w] * B0 + (R1 - B1[w]);

      B2[w] = B0;
      B1[w] = B0m;
    }
  }

  /* The last one or two rows, ending with row 0 */
  for ( ; l >= 0; l--)
  {
    const gdouble *row = c + (gsize) l * ld;

    #pragma omp simd

    for (w = 0; w < W; w++)
    {
      b1[w] = 0.0;
      b2[w] = 0.0;
    }

    for (k = n; k > 0; k--)
    {
      const gdouble x = row[k];

      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble b0 = sr2[w] * b1[w] + (x - b2[w]);

        b2[w] = b1[w];
        b1[w] = b0;
      }
    }

    if (l > 0)
    {
      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble R  = sr[w] * b1[w] + (row[0] - b2[w]);
        const gdouble B0 = st2[w] * B1[w] + (R - B2[w]);

        B2[w] = B1[w];
        B1[w] = B0;
      }
    }
    else
    {
      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble R = sr[w] * b1[w] + (row[0] - b2[w]);

        out[w] = st[w] * B1[w] + (R - B2[w]);
      }
    }
  }
}

/* W lanes for the last nl < W points: the missing lanes repeat the last point */
NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline void
_ncm_spectral_polar_panels_eval_lanes_tail (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, const guint W, const guint nl,
                                            const gdouble *r, const gdouble *theta, gdouble *out)
{
  gdouble rr[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], tt[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], oo[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  guint w;

  for (w = 0; w < W; w++)
  {
    rr[w] = r[MIN (w, nl - 1)];
    tt[w] = theta[MIN (w, nl - 1)];
  }

  _ncm_spectral_polar_panels_eval_lanes (panel, c, W, rr, tt, oo);

  for (w = 0; w < nl; w++)
    out[w] = oo[w];
}

/*
 * The series of panel, with the coefficients c, at the n points (r_j, theta_j). Eight
 * points per pass with AVX-512 and four otherwise, then four for a rest of two or more and
 * the scalar recurrence for a last point, where the lanes cost more: at one point four
 * lanes take 1.3 times the scalar time, at two they already take less.
 */
NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline void
_ncm_spectral_polar_panels_eval_points (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, guint n,
                                        const gdouble *r, const gdouble *theta, gdouble *out)
{
  guint j = 0;

#ifdef __AVX512F__

  for ( ; j + 8 <= n; j += 8)
    _ncm_spectral_polar_panels_eval_lanes (panel, c, 8, &r[j], &theta[j], &out[j]);

#endif /* __AVX512F__ */

  for ( ; j + 4 <= n; j += 4)
    _ncm_spectral_polar_panels_eval_lanes (panel, c, 4, &r[j], &theta[j], &out[j]);

  if (n - j >= 2)
    _ncm_spectral_polar_panels_eval_lanes_tail (panel, c, 4, n - j, &r[j], &theta[j], &out[j]);
  else if (n - j == 1)
    out[j] = _ncm_spectral_polar_panels_eval_on_c (panel, c, r[j], theta[j]);
}

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_PANELS_PRIVATE_H_ */

