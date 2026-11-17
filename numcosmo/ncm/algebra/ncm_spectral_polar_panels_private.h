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

/*
 * The series of one panel at W points, R angular rows at a time: the kernel of both
 * ncm_spectral_polar_panels_eval_on() (W = 1) and ncm_spectral_polar_panels_eval_on_array()
 * (W points). W and R are constants at each call site (the function is always inlined),
 * so the loops over them unroll and the loop over the points vectorizes. Each step of a
 * recurrence is 2 s b_1 + (c_k - b_2), see ncm_spectral_clenshaw_private.h. The radial
 * recurrences of different angular rows are independent: R of them run together, which
 * gives each point R independent chains, and each coefficient is loaded once for the W
 * points. Their results enter the angular recurrence one at a time, from the highest row
 * down, so the operations on a point do not depend on W or R; the two entry points then
 * differ only where the compiler contracts them into fma differently. r and theta are
 * read before out is written, so out may be r or theta.
 */
#define NCM_SPECTRAL_POLAR_PANELS_LANES_MAX 8
#define NCM_SPECTRAL_POLAR_PANELS_ROWS_MAX 4

/*
 * Rows l down to l - Q + 1 of the kernel below, Q a constant: their radial recurrences at
 * the W points, together, then their results into the angular recurrence (B1, B2), one
 * row at a time; at row 0 the angular recurrence ends and out is written.
 */
NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline void
_ncm_spectral_polar_panels_eval_rows (const gdouble *c, const gsize ld, const guint n, const guint W, const guint Q, const gint l,
                                      const gdouble *sr, const gdouble *sr2, const gdouble *st, const gdouble *st2,
                                      gdouble *B1, gdouble *B2, gdouble *out)
{
  gdouble b1[NCM_SPECTRAL_POLAR_PANELS_ROWS_MAX][NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  gdouble b2[NCM_SPECTRAL_POLAR_PANELS_ROWS_MAX][NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  guint k, q, w;

  for (q = 0; q < Q; q++)
  {
    #pragma omp simd

    for (w = 0; w < W; w++)
    {
      b1[q][w] = 0.0;
      b2[q][w] = 0.0;
    }
  }

  for (k = n; k > 0; k--)
  {
    for (q = 0; q < Q; q++)
    {
      const gdouble x = c[(gsize) (l - (gint) q) * ld + k];

      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble b0 = sr2[w] * b1[q][w] + (x - b2[q][w]);

        b2[q][w] = b1[q][w];
        b1[q][w] = b0;
      }
    }
  }

  for (q = 0; q < Q; q++)
  {
    const gdouble x0 = c[(gsize) (l - (gint) q) * ld];

    if (l - (gint) q > 0)
    {
      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble Rl = sr[w] * b1[q][w] + (x0 - b2[q][w]);
        const gdouble B0 = st2[w] * B1[w] + (Rl - B2[w]);

        B2[w] = B1[w];
        B1[w] = B0;
      }
    }
    else
    {
      #pragma omp simd

      for (w = 0; w < W; w++)
      {
        const gdouble Rl = sr[w] * b1[q][w] + (x0 - b2[q][w]);

        out[w] = st[w] * B1[w] + (Rl - B2[w]);
      }
    }
  }
}

NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline void
_ncm_spectral_polar_panels_eval_kernel (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, const guint W, const guint R,
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
  gint l = (gint) panel->n_theta - 1;
  guint w;

  #pragma omp simd

  for (w = 0; w < W; w++)
  {
    sr[w]  = (2.0 * r[w] - ra - rb) / (rb - ra);
    st[w]  = (2.0 * theta[w] - ta - tb) / (tb - ta);
    sr2[w] = 2.0 * sr[w];
    st2[w] = 2.0 * st[w];
    B1[w]  = 0.0;
    B2[w]  = 0.0;
  }

  /* Blocks of R rows, the last one possibly ending with row 0, then the rows left one at
   * a time */
  for ( ; l + 1 >= (gint) R; l -= (gint) R)
    _ncm_spectral_polar_panels_eval_rows (c, ld, n, W, R, l, sr, sr2, st, st2, B1, B2, out);

  for ( ; l >= 0; l--)
    _ncm_spectral_polar_panels_eval_rows (c, ld, n, W, 1, l, sr, sr2, st, st2, B1, B2, out);
}

/* W points for the last nl < W: the missing lanes repeat the last point */
NCM_SPECTRAL_POLAR_PANELS_ALWAYS_INLINE static inline void
_ncm_spectral_polar_panels_eval_kernel_tail (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, const guint W, const guint nl,
                                             const gdouble *r, const gdouble *theta, gdouble *out)
{
  gdouble rr[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], tt[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX], oo[NCM_SPECTRAL_POLAR_PANELS_LANES_MAX];
  guint w;

  for (w = 0; w < W; w++)
  {
    rr[w] = r[MIN (w, nl - 1)];
    tt[w] = theta[MIN (w, nl - 1)];
  }

  _ncm_spectral_polar_panels_eval_kernel (panel, c, W, 2, rr, tt, oo);

  for (w = 0; w < nl; w++)
    out[w] = oo[w];
}

/* The series of panel with the coefficients c at (r, theta), see ncm_spectral_polar_panels_eval_on() */
static inline gdouble
_ncm_spectral_polar_panels_eval_on_c (const NcmSpectralPolarPanelsPanel *panel, const gdouble *c, gdouble r, gdouble theta)
{
  gdouble out;

  /* Four radial recurrences at a time, the independent chains of one point */
  _ncm_spectral_polar_panels_eval_kernel (panel, c, 1, 4, &r, &theta, &out);

  return out;
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
    _ncm_spectral_polar_panels_eval_kernel (panel, c, 8, 2, &r[j], &theta[j], &out[j]);

#endif /* __AVX512F__ */

  for ( ; j + 4 <= n; j += 4)
    _ncm_spectral_polar_panels_eval_kernel (panel, c, 4, 2, &r[j], &theta[j], &out[j]);

  if (n - j >= 2)
    _ncm_spectral_polar_panels_eval_kernel_tail (panel, c, 4, n - j, &r[j], &theta[j], &out[j]);
  else if (n - j == 1)
    out[j] = _ncm_spectral_polar_panels_eval_on_c (panel, c, r[j], theta[j]);
}

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_PANELS_PRIVATE_H_ */

