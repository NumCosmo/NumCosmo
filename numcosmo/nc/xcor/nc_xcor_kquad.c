/***************************************************************************
 *            nc_xcor_kquad.c
 *
 *  Sat August 29 2026
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
 * The outer k quadrature: everything the kernel-space method
 * %NC_XCOR_METHOD_KERNEL_EXACT does with a pair of k-space closures, from the
 * block integrators and their GL(5) sweep through the knot merge.
 *
 * Every entry point here takes closures, or the kernels to build them from,
 * and none of them knows about the Limber tier or the multipole policy that
 * chose it -- that lives in nc_xcor.c.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "nc/xcor/nc_xcor.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "nc/xcor/nc_xcor_priv.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * Five-node Gauss-Legendre rule on [-1, 1], exact through degree 9. The outer
 * integrand k^2 W_i W_j is degree 8 on a merged knot panel of two cubic
 * closures, so five nodes are the smallest exact choice: four cap at degree 7.
 * A closure that is not piecewise cubic is handled by
 * _nc_xcor_kernel_integrate_block_spectral() instead.
 */
#define NC_XCOR_GL5_N 5

static const gdouble _nc_xcor_gl5_x[NC_XCOR_GL5_N] = {
  -0.9061798459386640, -0.5384693101056831, 0.0,
  0.5384693101056831, 0.9061798459386640
};

static const gdouble _nc_xcor_gl5_w[NC_XCOR_GL5_N] = {
  0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
  0.4786286704993665, 0.2369268850561891
};

/*
 * Accumulators for the closure-error estimate.
 *
 * The estimate propagates d(W1 W2) = |W1| dW2 + |W2| dW1, with dW_i the
 * closure fit's own error. A cell whose achieved residual was recorded by
 * nc_xcor_kernel_integrand_peek_closure_error() uses that residual and
 * accumulates into @res. A cell with no record -- tracking off, or a
 * refinement never accepted -- accumulates into @unk_i instead and is closed
 * afterwards with the requested tolerance times the peak.
 *
 * The estimate is a property of the closures, not of the quadrature above
 * them. What a cell is does change with the representation -- a merged knot
 * interval, or a cell of the common refinement of two panel sets -- but the
 * propagation above does not. Both exact paths therefore fill these
 * accumulators with the same sweep, and @vp_err means one thing whichever
 * representations the pair carries.
 *
 * The sweep samples |W| at the GL(5) nodes of every cell. These are error
 * terms, not the answer: |W| is not a polynomial, so no rule is exact on them.
 * Asking for @vp_err therefore costs one extra sweep of closure evaluations,
 * a fraction of a pair's total cost, which the closure build dominates.
 */
typedef struct _NcXcorClosureErr
{
  gdouble *res;          /* int k^2 (|W1| dW2 + |W2| dW1), cells with a record   */
  gdouble *unk1;         /* int k^2 |W2|    over cells where W1 has no record    */
  gdouble *unk2;         /* int k^2 |W1|    over cells where W2 has no record    */
  gdouble *prod1;        /* int k^2 |W1 W2| over cells where W1 has no record    */
  gdouble *prod2;        /* int k^2 |W1 W2| over cells where W2 has no record    */
  gdouble *peak1;        /* max |W1|                                             */
  gdouble *peak2;        /* max |W2|                                             */
  NcmMatrix *residuals1; /* achieved residuals, or NULL                   */
  NcmMatrix *residuals2;
  GArray *rows1; /* cell -> row of @residuals1 (guint)                   */
  GArray *rows2;
  gdouble *dW1; /* per-cell scratch: residual, or 0 where unknown       */
  gdouble *dW2;
  gdouble *m1; /* per-cell scratch: 1.0 where unknown, else 0.0        */
  gdouble *m2;

  /* One double per multipole per accumulator, and a block is capped at
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK by the closure builder -- 512 bytes each. The
   * pointers above index into this, aliased in pairs on an auto spectrum,
   * where one closure means one set of residuals to accumulate. */
  gdouble store[11][NC_XCOR_KERNEL_MAX_ELL_BLOCK];
} NcXcorClosureErr;

/*
 * Fills the per-cell dW/mask scratch for one side from its recorded
 * residuals. A NaN entry means the interval was never accepted, and is treated
 * exactly as no record at all.
 */
static void
_nc_xcor_closure_cell_residual (NcmMatrix *residuals, GArray *rows, const guint ie, const guint nell, gdouble *dW, gdouble *m)
{
  guint il;

  if (residuals == NULL)
  {
    for (il = 0; il < nell; il++)
    {
      dW[il] = 0.0;
      m[il]  = 1.0;
    }

    return;
  }

  {
    const guint row = g_array_index (rows, guint, ie);

    for (il = 0; il < nell; il++)
    {
      const gdouble d = ncm_matrix_get (residuals, row, il);

      dW[il] = gsl_finite (d) ? d : 0.0;
      m[il]  = gsl_finite (d) ? 0.0 : 1.0;
    }
  }
}

/*
 * Maps each cell of @edges onto the row of the closure's residual record that
 * covers it -- the knot interval of a spline closure, the panel of a spectral
 * one. Both the cells and the closure's own breakpoints are sorted and every
 * cell lies inside one of those intervals, so a single marching index does it.
 */
static GArray *
_nc_xcor_closure_cell_rows (NcXcorKernelIntegrand *xclki, GArray *edges, const gdouble c)
{
  const guint n_panels = nc_xcor_kernel_integrand_get_n_panels (xclki);
  GArray *rows         = g_array_sized_new (FALSE, FALSE, sizeof (guint), edges->len);
  guint ie, j = 0;

  if (n_panels > 0)
  {
    for (ie = 0; ie + 1 < edges->len; ie++)
    {
      const gdouble cell_lo = c * g_array_index (edges, gdouble, ie);

      while (j + 1 < n_panels)
      {
        NcmMatrix *ignored = NULL;
        gdouble a, b;

        nc_xcor_kernel_integrand_peek_panel (xclki, j, &ignored, &a, &b);

        if (b > cell_lo)
          break;

        j++;
      }

      g_array_append_val (rows, j);
    }
  }
  else
  {
    NcmVector *knots   = nc_xcor_kernel_integrand_peek_knots (xclki);
    const guint nknots = ncm_vector_len (knots);

    for (ie = 0; ie + 1 < edges->len; ie++)
    {
      const gdouble cell_lo = c * g_array_index (edges, gdouble, ie);

      while ((j + 2 < nknots) && (ncm_vector_get (knots, j + 1) <= cell_lo))
        j++;

      g_array_append_val (rows, j);
    }
  }

  return rows;
}

static void
_nc_xcor_gl5_sweep_auto (NcXcorKernelIntegrand *xclki, GArray *edges, guint nell, gdouble *W, gdouble *sum)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);
    const gdouble panel_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid      = 0.5 * (panel_lo + panel_hi);
    const gdouble half     = 0.5 * (panel_hi - panel_lo);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki, k, W);

      for (il = 0; il < nell; il++)
        sum[il] += w * W[il] * W[il];
    }
  }
}

static void
_nc_xcor_gl5_sweep_cross (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, GArray *edges, guint nell, gdouble *W1, gdouble *W2, gdouble *sum)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);
    const gdouble panel_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid      = 0.5 * (panel_lo + panel_hi);
    const gdouble half     = 0.5 * (panel_hi - panel_lo);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki1, k, W1);
      nc_xcor_kernel_integrand_eval (xclki2, k, W2);

      for (il = 0; il < nell; il++)
        sum[il] += w * W1[il] * W2[il];
    }
  }
}

/*
 * Wires the accumulators onto their storage and reads each closure's residual
 * record. The auto case aliases every pair onto one buffer: one closure means
 * one set of residuals to accumulate.
 */
static void
_nc_xcor_closure_err_init (NcXcorClosureErr *err, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gboolean isauto)
{
  guint i;

  for (i = 0; i < G_N_ELEMENTS (err->store); i++)
    memset (err->store[i], 0, sizeof (err->store[i]));

  err->res   = err->store[0];
  err->unk1  = err->store[1];
  err->unk2  = isauto ? err->store[1] : err->store[2];
  err->prod1 = err->store[3];
  err->prod2 = isauto ? err->store[3] : err->store[4];
  err->peak1 = err->store[5];
  err->peak2 = isauto ? err->store[5] : err->store[6];
  err->dW1   = err->store[7];
  err->dW2   = isauto ? err->store[7] : err->store[8];
  err->m1    = err->store[9];
  err->m2    = isauto ? err->store[9] : err->store[10];

  err->residuals1 = nc_xcor_kernel_integrand_peek_closure_error (xclki1);
  err->residuals2 = isauto ? err->residuals1 : nc_xcor_kernel_integrand_peek_closure_error (xclki2);
  err->rows1      = NULL;
  err->rows2      = NULL;
}

/* The cell-to-row maps of one edge set, in each closure's own variable x_j = c_j x */
static void
_nc_xcor_closure_err_set_rows (NcXcorClosureErr *err, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2,
                               gboolean isauto, GArray *edges, const gdouble c1, const gdouble c2)
{
  g_clear_pointer (&err->rows1, g_array_unref);

  if (!isauto)
    g_clear_pointer (&err->rows2, g_array_unref);

  err->rows1 = (err->residuals1 != NULL) ? _nc_xcor_closure_cell_rows (xclki1, edges, c1) : NULL;
  err->rows2 = isauto ? err->rows1 :
               ((err->residuals2 != NULL) ? _nc_xcor_closure_cell_rows (xclki2, edges, c2) : NULL);
}

/*
 * The multipoles [il0, il0 + n) at the node x of a cell: k = sigma[il] x, one
 * component lookup each, or one evaluation of the whole vector when sigma is
 * NULL, which stands for every scale equal to 1.
 */
static void
_nc_xcor_closure_eval_group (NcXcorKernelIntegrand *xclki, const gdouble x, const guint il0, const guint n,
                             const gdouble *sigma, gdouble *W)
{
  guint il;

  if (sigma == NULL)
  {
    nc_xcor_kernel_integrand_eval (xclki, x, W);

    return;
  }

  for (il = il0; il < il0 + n; il++)
    nc_xcor_kernel_integrand_eval_comps (xclki, sigma[il] * x, il, 1, W);
}

/*
 * Samples |W| at the GL(5) nodes of every cell and accumulates the pair's
 * error terms. It does not compute the integral: that is left to whichever
 * exact path matches the representations.
 *
 * The auto and cross cases are separate sweeps rather than one with a test
 * inside, because the aliasing that makes them one set of accumulators also
 * means the auto case must add each term once and let the assembly double it.
 */
static void
_nc_xcor_closure_err_sweep_auto (NcXcorKernelIntegrand *xclki, GArray *edges, const guint il0, const guint n,
                                 const gdouble *sigma, gdouble *W, NcXcorClosureErr *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble cell_lo = g_array_index (edges, gdouble, ie);
    const gdouble cell_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid     = 0.5 * (cell_lo + cell_hi);
    const gdouble half    = 0.5 * (cell_hi - cell_lo);

    _nc_xcor_closure_cell_residual (err->residuals1, err->rows1, ie, il0 + n, err->dW1, err->m1);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble x = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * x * x;

      _nc_xcor_closure_eval_group (xclki, x, il0, n, sigma, W);

      for (il = il0; il < il0 + n; il++)
      {
        const gdouble s3   = (sigma != NULL) ? gsl_pow_3 (sigma[il]) : 1.0;
        const gdouble absW = fabs (W[il]);

        /* d(W^2) = 2 |W| dW, and the aliased unk2/prod2/peak2 supply the
         * second half of the unknown-cell term in the assembly. */
        err->res[il]   += 2.0 * s3 * w * absW * err->dW1[il];
        err->unk1[il]  += s3 * w * absW * err->m1[il];
        err->prod1[il] += fabs (s3 * w * W[il] * W[il]) * err->m1[il];
        err->peak1[il]  = GSL_MAX (err->peak1[il], absW);
      }
    }
  }
}

static void
_nc_xcor_closure_err_sweep_cross (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, GArray *edges,
                                  const guint il0, const guint n, const gdouble *sigma,
                                  gdouble *W1, gdouble *W2, NcXcorClosureErr *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble cell_lo = g_array_index (edges, gdouble, ie);
    const gdouble cell_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid     = 0.5 * (cell_lo + cell_hi);
    const gdouble half    = 0.5 * (cell_hi - cell_lo);

    _nc_xcor_closure_cell_residual (err->residuals1, err->rows1, ie, il0 + n, err->dW1, err->m1);
    _nc_xcor_closure_cell_residual (err->residuals2, err->rows2, ie, il0 + n, err->dW2, err->m2);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble x = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * x * x;

      _nc_xcor_closure_eval_group (xclki1, x, il0, n, sigma, W1);
      _nc_xcor_closure_eval_group (xclki2, x, il0, n, sigma, W2);

      for (il = il0; il < il0 + n; il++)
      {
        const gdouble s3    = (sigma != NULL) ? gsl_pow_3 (sigma[il]) : 1.0;
        const gdouble term  = s3 * w * W1[il] * W2[il];
        const gdouble absW1 = fabs (W1[il]);
        const gdouble absW2 = fabs (W2[il]);

        err->res[il]   += s3 * w * (absW1 * err->dW2[il] + absW2 * err->dW1[il]);
        err->unk1[il]  += s3 * w * absW2 * err->m1[il];
        err->unk2[il]  += s3 * w * absW1 * err->m2[il];
        err->prod1[il] += fabs (term) * err->m1[il];
        err->prod2[il] += fabs (term) * err->m2[il];
        err->peak1[il]  = GSL_MAX (err->peak1[il], absW1);
        err->peak2[il]  = GSL_MAX (err->peak2[il], absW2);
      }
    }
  }
}

/*
 * Closes the estimate and writes it into @vp_err.
 *
 * The quadrature is exact on both paths, so the only error is the closures'
 * own, propagated through d(W1 W2) = |W1| dW2 + |W2| dW1. Where a closure
 * recorded what its fit achieved, the sweep has already integrated that; what
 * is left is to close the cells that carry no record with the tolerance the fit
 * was asked for. That fallback keeps the two halves of the criterion apart the
 * way the criterion does -- the relative one riding on the product, the
 * peak-scaled floor against the other closure's amplitude -- so with
 * #NcXcorKernel:track-closure-error off it is the whole estimate, and is then
 * exactly the tolerance-only bound.
 */
static void
_nc_xcor_closure_err_assemble (NcXcorClosureErr *err, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gboolean isauto, guint nell, gdouble const_factor, NcmVector *vp_err)
{
  const gdouble reltol1 = nc_xcor_kernel_integrand_get_reltol (xclki1);
  const gdouble reltol2 = nc_xcor_kernel_integrand_get_reltol (xclki2);
  const gdouble sabs1   = nc_xcor_kernel_integrand_get_peak_epsilon (xclki1);
  const gdouble sabs2   = nc_xcor_kernel_integrand_get_peak_epsilon (xclki2);
  guint il;

  for (il = 0; il < nell; il++)
  {
    const gdouble unk_term = reltol1 * err->prod1[il] + sabs1 * err->peak1[il] * err->unk1[il] +
                             reltol2 * err->prod2[il] + sabs2 * err->peak2[il] * err->unk2[il];

    ncm_vector_set (vp_err, il, const_factor * (err->res[il] + unk_term));
  }

  g_clear_pointer (&err->rows1, g_array_unref);

  if (!isauto)
    g_clear_pointer (&err->rows2, g_array_unref);
}

/*
 * Common refinement of two knot sets, clipped to [@k_min, @k_max]. Both are
 * sorted, so this is a linear merge; duplicates are dropped so no zero-width
 * panel survives. The result is pre-sized to the exact upper bound of a merge,
 * so the append loop never reallocates: the union of two sorted sets cannot
 * exceed their combined length.
 */
static GArray *
_nc_xcor_merge_knots (NcmVector *knots1, NcmVector *knots2, gdouble k_min, gdouble k_max)
{
  const guint len1 = ncm_vector_len (knots1);
  const guint len2 = ncm_vector_len (knots2);
  GArray *edges    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len1 + len2);
  guint i1         = 0;
  guint i2         = 0;

  while ((i1 < len1) || (i2 < len2))
  {
    const gdouble x1 = (i1 < len1) ? ncm_vector_get (knots1, i1) : GSL_POSINF;
    const gdouble x2 = (i2 < len2) ? ncm_vector_get (knots2, i2) : GSL_POSINF;
    const gdouble x  = GSL_MIN (x1, x2);

    if (x1 <= x2)
      i1++;

    if (x2 <= x1)
      i2++;

    if ((x < k_min) || (x > k_max))
      continue;

    if ((edges->len > 0) && (x <= g_array_index (edges, gdouble, edges->len - 1)))
      continue;

    g_array_append_val (edges, x);
  }

  return edges;
}

/*
 * Exact outer quadrature for a pair on the union of both knot sets. Each
 * closure is cubic on a merged panel, so k^2 W_1 W_2 has degree at most 8 and
 * GL(5) is exact.
 *
 * The integration range is the intersection of the two fitted domains,
 * because NcmSpline does not range-check and an out-of-domain evaluation
 * returns an extrapolation rather than a small number.
 *
 * Refining every panel fourfold moves the result by 1e-15 to 1e-12, so an
 * embedded rule (Kronrod, or GL(5) against GL(9)) would report machine zero
 * on every call. Do not add one. The remaining error is the closure fit's,
 * amplified by cancellation in C_ell: two disjoint Gaussian bins cancel by a
 * factor 1.4e4 at ell = 9, so a closure good to 1e-8 leaves 1e-4 on C_ell.
 * @vp_err reports that product; see nc_xcor_compute_full().
 */
static gint
_nc_xcor_cmp_edge (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/*
 * Appends the breakpoints of one closure that fall strictly inside
 * (k_min, k_max): the panel edges of a spectral closure, the knots of a spline
 * one. Between two consecutive breakpoints a closure is a single polynomial,
 * which is the only property the common refinement below needs -- so the two
 * representations enter it on the same footing.
 */
static void
_nc_xcor_append_breakpoints (NcXcorKernelIntegrand *xclki, gdouble x_min, gdouble x_max, const gdouble c, GArray *edges)
{
  const guint n_panels = nc_xcor_kernel_integrand_get_n_panels (xclki);
  guint i;

  if (n_panels > 0)
  {
    for (i = 0; i < n_panels; i++)
    {
      NcmMatrix *ignored = NULL;
      gdouble a, b;

      nc_xcor_kernel_integrand_peek_panel (xclki, i, &ignored, &a, &b);
      b /= c;

      if ((b > x_min) && (b < x_max))
        g_array_append_val (edges, b);
    }
  }
  else
  {
    NcmVector *knots   = nc_xcor_kernel_integrand_peek_knots (xclki);
    const guint nknots = (knots != NULL) ? ncm_vector_len (knots) : 0;

    for (i = 0; i < nknots; i++)
    {
      const gdouble knot = ncm_vector_get (knots, i) / c;

      if ((knot > x_min) && (knot < x_max))
        g_array_append_val (edges, knot);
    }
  }
}

/*
 * Merges two closures' breakpoints over [k_min, k_max]. The result is the
 * common refinement, on each cell of which both closures are a single
 * polynomial -- the same argument the merged knot sets make for two splines,
 * stated so that a spline and a panel set can meet on it.
 */
static GArray *
_nc_xcor_merge_panel_edges (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2,
                            gboolean isauto, gdouble x_min, gdouble x_max, const gdouble c1, const gdouble c2)
{
  GArray *edges = g_array_new (FALSE, FALSE, sizeof (gdouble));

  g_array_append_val (edges, x_min);

  _nc_xcor_append_breakpoints (xclki1, x_min, x_max, c1, edges);

  if (!isauto)
    _nc_xcor_append_breakpoints (xclki2, x_min, x_max, c2, edges);

  g_array_sort (edges, _nc_xcor_cmp_edge);
  g_array_append_val (edges, x_max);

  /* Drop duplicates: two closures often break at the same place. */
  {
    guint w = 1;
    guint i;

    for (i = 1; i < edges->len; i++)
      if (g_array_index (edges, gdouble, i) > g_array_index (edges, gdouble, w - 1))
        g_array_index (edges, gdouble, w++) = g_array_index (edges, gdouble, i);

    g_array_set_size (edges, w);
  }

  return edges;
}

/*
 * Chebyshev coefficients of one closure on the cell [a, b] of the common
 * refinement, as an @nell by N matrix.
 *
 * A spectral closure restricts its panel's polynomial onto the cell, which is
 * a change of basis. A spline closure is a single cubic there -- the merged
 * breakpoints carry its knots, so the cell lies inside one knot interval -- and
 * interpolation at N + 1 Chebyshev-Lobatto nodes is exact for a polynomial of
 * degree N, so four nodes reproduce that cubic exactly. This is what lets a
 * mixed pair be integrated exactly rather than refused: the two
 * representations meet in coefficient space, on cells where each is one
 * polynomial.
 *
 * The four-node transform is written out instead of taken from a DCT because
 * at four nodes it is four sums.
 */
static void
_nc_xcor_cell_coeffs (NcXcorKernelIntegrand *xclki, gdouble a, gdouble b, guint nell, const gchar *side, NcmMatrix **coeffs_ptr)
{
  NcmMatrix *coeffs = *coeffs_ptr;

  if (nc_xcor_kernel_integrand_get_n_panels (xclki) > 0)
  {
    /* Every cell of the common refinement lies inside one panel -- its edges
     * are the panels' own, so the containment test in restrict() compares
     * identical doubles. A failure here would mean the refinement and the
     * panels disagree, and dropping the cell would return a quietly wrong
     * C_ell. */
    if (!nc_xcor_kernel_integrand_restrict (xclki, a, b, coeffs_ptr))
      g_error ("_nc_xcor_kernel_integrate_block_spectral: cell [%.17g, %.17g] "
               "is not inside a single panel of the %s closure.", a, b, side);

    return;
  }

  {
    /* The Chebyshev-Lobatto nodes of order three, cos(j pi / 3). */
    static const gdouble node_t[4] = { 1.0, 0.5, -0.5, -1.0 };
    gdouble f[4][NC_XCOR_KERNEL_MAX_ELL_BLOCK];
    const gdouble mid  = 0.5 * (a + b);
    const gdouble half = 0.5 * (b - a);
    guint j, il;

    for (j = 0; j < 4; j++)
      nc_xcor_kernel_integrand_eval (xclki, mid + half * node_t[j], f[j]);

    if ((coeffs != NULL) && ((ncm_matrix_nrows (coeffs) != nell) || (ncm_matrix_ncols (coeffs) != 4)))
      ncm_matrix_clear (&coeffs);

    if (coeffs == NULL)
      coeffs = ncm_matrix_new (nell, 4);

    *coeffs_ptr = coeffs;

    for (il = 0; il < nell; il++)
    {
      const gdouble f0 = f[0][il];
      const gdouble f1 = f[1][il];
      const gdouble f2 = f[2][il];
      const gdouble f3 = f[3][il];

      ncm_matrix_set (coeffs, il, 0, (0.5 * f0 + f1 + f2 + 0.5 * f3) / 3.0);
      ncm_matrix_set (coeffs, il, 1, (f0 + f1 - f2 - f3) / 3.0);
      ncm_matrix_set (coeffs, il, 2, (f0 - f1 - f2 + f3) / 3.0);
      ncm_matrix_set (coeffs, il, 3, (0.5 * f0 - f1 + f2 - 0.5 * f3) / 3.0);
    }

    return;
  }
}

/*
 * Scratch of the spectral integration of one block. The tables A_k[i][j] =
 * int_{-1}^{1} T_i T_j T_k dt, k = 0, 1, 2, depend only on (i, j); up to
 * NC_XCOR_CHEB_PRODUCT_TABLE_N coefficients they are built once per process,
 * above it once per call for the largest count seen. G holds one cell's
 * bilinear form.
 */
typedef struct _NcXcorSpectralScratch
{
  GArray *A;            /* tables built here, for counts above the process-wide ones */
  GArray *G;            /* n1 n2 doubles */
  const gdouble *A_ptr; /* the tables in use: A_0, A_1, A_2, row-major in (i, j) */
  guint nmax;           /* their row stride */
  guint nmax_local;     /* the count A was built for */
} NcXcorSpectralScratch;

/* inv(m) = 1 / (1 - m^2) for even m, 0 for odd m: half the integral of T_m over [-1, 1] */
static inline gdouble
_nc_xcor_cheb_half_int (const guint m)
{
  return (m % 2 == 0) ? 1.0 / (1.0 - (gdouble) m * (gdouble) m) : 0.0;
}

/*
 * Fills the tables A_k for coefficient counts n1 and n2, from
 * T_i T_j = (T_{i+j} + T_{|i-j|}) / 2 applied twice:
 * A_k[i][j] = (inv(p + k) + inv(|p - k|) + inv(q + k) + inv(|q - k|)) / 2, with
 * p = i + j and q = |i - j|.
 */

/* Fills the three tables for coefficient counts below nmax, row stride nmax */
static void
_nc_xcor_cheb_product_tables_fill (gdouble *A, const guint nmax)
{
  guint i, j, k;

  for (k = 0; k < 3; k++)
  {
    for (i = 0; i < nmax; i++)
    {
      for (j = 0; j < nmax; j++)
      {
        const guint p = i + j;
        const guint q = (i > j) ? i - j : j - i;

        A[(k * nmax + i) * nmax + j] = 0.5 * (_nc_xcor_cheb_half_int (p + k) +
                                              _nc_xcor_cheb_half_int ((p > k) ? p - k : k - p) +
                                              _nc_xcor_cheb_half_int (q + k) +
                                              _nc_xcor_cheb_half_int ((q > k) ? q - k : k - q));
      }
    }
  }
}

/* The tables are constants: up to this count they are built once per process. */
#define NC_XCOR_CHEB_PRODUCT_TABLE_N (65)

static const gdouble *
_nc_xcor_cheb_product_tables_static (void)
{
  static gsize init = 0;
  static gdouble *A = NULL;

  if (g_once_init_enter (&init))
  {
    A = g_new (gdouble, 3 * NC_XCOR_CHEB_PRODUCT_TABLE_N * NC_XCOR_CHEB_PRODUCT_TABLE_N);
    _nc_xcor_cheb_product_tables_fill (A, NC_XCOR_CHEB_PRODUCT_TABLE_N);
    g_once_init_leave (&init, 1);
  }

  return A;
}

/*
 * Points sc->A_ptr at tables valid for counts n1 and n2, with row stride
 * sc->nmax: the process-wide tables when both counts fit, otherwise tables
 * built here for the larger count.
 */
static void
_nc_xcor_spectral_scratch_prepare (NcXcorSpectralScratch *sc, const guint n1, const guint n2)
{
  const guint nmax = MAX (n1, n2);

  g_array_set_size (sc->G, n1 * n2);

  if (nmax <= NC_XCOR_CHEB_PRODUCT_TABLE_N)
  {
    sc->A_ptr = _nc_xcor_cheb_product_tables_static ();
    sc->nmax  = NC_XCOR_CHEB_PRODUCT_TABLE_N;

    return;
  }

  if (nmax > sc->nmax_local)
  {
    g_array_set_size (sc->A, 3 * nmax * nmax);
    _nc_xcor_cheb_product_tables_fill (&g_array_index (sc->A, gdouble, 0), nmax);
    sc->nmax_local = nmax;
  }

  sc->A_ptr = &g_array_index (sc->A, gdouble, 0);
  sc->nmax  = sc->nmax_local;
}

/* A closure's range in its own variable: the span of its panels, or its k range */
static void
_nc_xcor_closure_var_range (NcXcorKernelIntegrand *xclki, gdouble *lo, gdouble *hi)
{
  const guint n_panels = nc_xcor_kernel_integrand_get_n_panels (xclki);

  if (n_panels > 0)
  {
    NcmMatrix *ignored = NULL;
    gdouble a, b;

    nc_xcor_kernel_integrand_peek_panel (xclki, 0, &ignored, lo, &b);
    nc_xcor_kernel_integrand_peek_panel (xclki, n_panels - 1, &ignored, &a, hi);
  }
  else
  {
    nc_xcor_kernel_integrand_get_range (xclki, lo, hi);
  }
}

/*
 * Integrates the multipoles [il0, il0 + n) of a pair on the common refinement
 * of the two closures' breakpoints in a variable x, with k = sigma[il] x for
 * multipole il (sigma NULL: every scale 1) and closure j in x_j = c_j x. On
 * each cell both closures are polynomials in the cell's variable t, so the
 * product is integrated as a bilinear form of their coefficients, with k^2 dk
 * = sigma^3 x^2 dx. Adds sigma^3 times the integral in x to sum[il], and the
 * closure error terms to err when it is not NULL.
 */
static void
_nc_xcor_spectral_group (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gboolean isauto,
                         const guint il0, const guint n, const gdouble *sigma, const gdouble c1, const gdouble c2,
                         gdouble *sum, gdouble *W1, gdouble *W2, NcXcorClosureErr *err, NcXcorSpectralScratch *sc)
{
  const guint len1 = nc_xcor_kernel_integrand_get_len (xclki1);
  const guint len2 = nc_xcor_kernel_integrand_get_len (xclki2);
  NcmMatrix *cell1 = NULL;
  NcmMatrix *cell2 = NULL;
  gdouble lo1, hi1, lo2, hi2, x_min, x_max;
  GArray *edges;
  guint ie, il;

  _nc_xcor_closure_var_range (xclki1, &lo1, &hi1);
  _nc_xcor_closure_var_range (xclki2, &lo2, &hi2);

  x_min = GSL_MAX (lo1 / c1, lo2 / c2);
  x_max = GSL_MIN (hi1 / c1, hi2 / c2);

  if (x_min >= x_max)
    return;

  edges = _nc_xcor_merge_panel_edges (xclki1, xclki2, isauto, x_min, x_max, c1, c2);

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble a    = g_array_index (edges, gdouble, ie);
    const gdouble b    = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid  = 0.5 * (a + b);
    const gdouble half = 0.5 * (b - a);
    NcmMatrix *cm1, *cm2;

    _nc_xcor_cell_coeffs (xclki1, a * c1, b * c1, len1, "first", &cell1);
    cm1 = cell1;

    if (isauto)
    {
      cm2 = cell1;
    }
    else
    {
      _nc_xcor_cell_coeffs (xclki2, a * c2, b * c2, len2, "second", &cell2);
      cm2 = cell2;
    }

    {
      /* x^2 = w0 T_0 + w1 T_1 + w2 T_2 in the cell's variable x = mid + half t, so
       * int x^2 p q dx = half a^T G b with G = w0 A_0 + w1 A_1 + w2 A_2, the
       * same G for every multipole of the cell. */
      const gdouble w0   = mid * mid + 0.5 * half * half;
      const gdouble w1   = 2.0 * mid * half;
      const gdouble w2   = 0.5 * half * half;
      const guint n1     = ncm_matrix_ncols (cm1);
      const guint n2     = ncm_matrix_ncols (cm2);
      const guint tda1   = ncm_matrix_tda (cm1);
      const guint tda2   = ncm_matrix_tda (cm2);
      const gdouble *dc1 = ncm_matrix_data (cm1);
      const gdouble *dc2 = ncm_matrix_data (cm2);
      const gdouble *A0, *A1, *A2;
      gdouble *G;
      guint i, j;

      _nc_xcor_spectral_scratch_prepare (sc, n1, n2);

      A0 = sc->A_ptr;
      A1 = A0 + sc->nmax * sc->nmax;
      A2 = A1 + sc->nmax * sc->nmax;
      G  = &g_array_index (sc->G, gdouble, 0);

      for (i = 0; i < n1; i++)
      {
        const guint row = i * sc->nmax;

        for (j = 0; j < n2; j++)
          G[i * n2 + j] = w0 * A0[row + j] + w1 * A1[row + j] + w2 * A2[row + j];
      }

      for (il = il0; il < il0 + n; il++)
      {
        const gdouble s3 = (sigma != NULL) ? gsl_pow_3 (sigma[il]) : 1.0;
        const gdouble *a = &dc1[il * tda1];
        const gdouble *b = &dc2[il * tda2];
        gdouble acc      = 0.0;

        for (i = 0; i < n1; i++)
        {
          const gdouble *Gi = &G[i * n2];
          gdouble r         = 0.0;

          for (j = 0; j < n2; j++)
            r += Gi[j] * b[j];

          acc += a[i] * r;
        }

        sum[il] += s3 * half * acc;
      }
    }
  }

  ncm_matrix_clear (&cell1);
  ncm_matrix_clear (&cell2);

  if (err != NULL)
  {
    _nc_xcor_closure_err_set_rows (err, xclki1, xclki2, isauto, edges, c1, c2);

    if (isauto)
      _nc_xcor_closure_err_sweep_auto (xclki1, edges, il0, n, sigma, W1, err);
    else
      _nc_xcor_closure_err_sweep_cross (xclki1, xclki2, edges, il0, n, sigma, W1, W2, err);
  }

  g_array_unref (edges);
}

/*
 * Integrates a pair with a spectral closure on at least one side. Two closures
 * in k, or a spectral closure with a spline, share one edge set for the block.
 * Two closures in u = k / nu share one edge set too, with sigma = nu per
 * multipole. A closure in u paired with one in k has no common edge set, and
 * each multipole is integrated on its own refinement, the k closure's
 * breakpoints divided by that multipole's nu.
 */
static void
_nc_xcor_kernel_integrate_block_spectral (NcXcor *xc, NcXcorKernelIntegrand *xclki1,
                                          NcXcorKernelIntegrand *xclki2, guint lmin,
                                          guint lmax, gboolean isauto, NcmVector *vp,
                                          NcmVector *vp_err)
{
  const guint nell           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcXcorClosureErr err_acc;
  gboolean scaled1 = FALSE;
  gboolean scaled2 = FALSE;

  /* One accumulator per multipole; a block is capped at
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK by the closure builder. */
  gdouble sum[NC_XCOR_KERNEL_MAX_ELL_BLOCK] = { 0.0 };
  gdouble s1[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble s2[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble W1_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble W2_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  NcXcorSpectralScratch sc;
  guint il;

  ncm_vector_set_zero (vp);

  if (vp_err != NULL)
  {
    ncm_vector_set_zero (vp_err);
    _nc_xcor_closure_err_init (&err_acc, xclki1, xclki2, isauto);
  }

  for (il = 0; il < nell; il++)
  {
    s1[il]   = nc_xcor_kernel_integrand_get_scale (xclki1, il);
    s2[il]   = nc_xcor_kernel_integrand_get_scale (xclki2, il);
    scaled1 |= (s1[il] != 1.0);
    scaled2 |= (s2[il] != 1.0);
  }

  /* Local rather than kept on @xc: the solver's OpenMP team shares @xc across threads. */
  sc.A          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sc.G          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sc.nmax       = 0;
  sc.nmax_local = 0;
  sc.A_ptr      = NULL;

  if (!scaled1 && !scaled2)
  {
    _nc_xcor_spectral_group (xclki1, xclki2, isauto, 0, nell, NULL, 1.0, 1.0,
                             sum, W1_store, W2_store, (vp_err != NULL) ? &err_acc : NULL, &sc);
  }
  else if (scaled1 && scaled2)
  {
    for (il = 0; il < nell; il++)
      if (s1[il] != s2[il])
        g_error ("_nc_xcor_kernel_integrate_block_spectral: the two closures scale multipole %u "
                 "differently (%.17g and %.17g).", lmin + il, s1[il], s2[il]);

    _nc_xcor_spectral_group (xclki1, xclki2, isauto, 0, nell, s1, 1.0, 1.0,
                             sum, W1_store, W2_store, (vp_err != NULL) ? &err_acc : NULL, &sc);
  }
  else
  {
    const gdouble *sigma = scaled1 ? s1 : s2;

    for (il = 0; il < nell; il++)
      _nc_xcor_spectral_group (xclki1, xclki2, isauto, il, 1, sigma,
                               scaled1 ? 1.0 : sigma[il], scaled2 ? 1.0 : sigma[il],
                               sum, W1_store, W2_store, (vp_err != NULL) ? &err_acc : NULL, &sc);
  }

  for (il = 0; il < nell; il++)
    ncm_vector_set (vp, il, const_factor * sum[il]);

  /* The same estimate the merged-knot path reports, on the same cells: the
   * bilinear form is exact, so what is left is the closures' own fit error. */
  if (vp_err != NULL)
    _nc_xcor_closure_err_assemble (&err_acc, xclki1, xclki2, isauto, nell, const_factor, vp_err);

  g_array_unref (sc.A);
  g_array_unref (sc.G);
}

void
_nc_xcor_kernel_integrate_block_exact (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell             = lmax - lmin + 1;
  const gdouble const_factor   = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcXcorKernelIntegrand *side2 = isauto ? xclki1 : xclki2;
  NcmVector *knots1            = nc_xcor_kernel_integrand_peek_knots (xclki1);
  NcmVector *knots2            = nc_xcor_kernel_integrand_peek_knots (side2);
  gdouble k_min1, k_max1, k_min2, k_max2, k_min, k_max;
  NcXcorClosureErr err_acc;

  /* One double per multipole, and a block is capped at
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK by the closure builder -- 512 bytes each. On
   * the stack they need no allocation, and no exit path from this function has
   * anything of theirs to free. */
  gdouble sum[NC_XCOR_KERNEL_MAX_ELL_BLOCK] = { 0.0 };
  gdouble W1_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble W2_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble *W1, *W2;
  GArray *edges;
  guint il;

  if (ncm_vector_len (vp) != nell)
    g_error ("_nc_xcor_kernel_integrate_block_exact: vector size does not match multipole limits");

  if ((vp_err != NULL) && (ncm_vector_len (vp_err) != nell))
    g_error ("_nc_xcor_kernel_integrate_block_exact: error vector size does not match multipole limits");

  /* The scratch here and in the spectral path is sized at the cap
   * nc_xcor_kernel_get_eval_vectorized_full() enforces; an integrand built
   * some other way must still respect it. */
  if ((nell > NC_XCOR_KERNEL_MAX_ELL_BLOCK) ||
      (nc_xcor_kernel_integrand_get_len (xclki1) > NC_XCOR_KERNEL_MAX_ELL_BLOCK) ||
      (nc_xcor_kernel_integrand_get_len (side2) > NC_XCOR_KERNEL_MAX_ELL_BLOCK))
    g_error ("_nc_xcor_kernel_integrate_block_exact: block of %u multipoles exceeds "
             "NC_XCOR_KERNEL_MAX_ELL_BLOCK (%d).", nell, NC_XCOR_KERNEL_MAX_ELL_BLOCK);

  /* Chosen here rather than by the callers: NcXcorSolver and
   * _nc_xcor_kernel_space_compute() both enter through this function, and a
   * choice made in one of them is a choice the other silently misses.
   *
   * A pair with a panel-backed closure on at least one side goes to the common
   * refinement of the two breakpoint sets, where the product is a polynomial
   * and the integral is a bilinear form in the coefficients. Two splines take
   * the merged-knot GL(5) sweep below: on a merged interval the product is
   * two cubics times k^2, degree 8, and five-point Gauss-Legendre is exact
   * through degree 9. Either way the closures handed in are integrated
   * exactly. */
  if ((nc_xcor_kernel_integrand_get_n_panels (xclki1) > 0) ||
      (nc_xcor_kernel_integrand_get_n_panels (side2) > 0))
  {
    _nc_xcor_kernel_integrate_block_spectral (xc, xclki1, side2, lmin, lmax, isauto, vp, vp_err);

    return;
  }

  if ((knots1 == NULL) || (knots2 == NULL))
    g_error ("_nc_xcor_kernel_integrate_block_exact: %s method requires spline-backed "
             "integrands, which report their knots.", "NC_XCOR_METHOD_KERNEL_EXACT");

  nc_xcor_kernel_integrand_get_range (xclki1, &k_min1, &k_max1);
  nc_xcor_kernel_integrand_get_range (xclki2, &k_min2, &k_max2);

  k_min = GSL_MAX (k_min1, k_min2);
  k_max = GSL_MIN (k_max1, k_max2);

  ncm_vector_set_zero (vp);

  if (vp_err != NULL)
    ncm_vector_set_zero (vp_err);

  if (k_min >= k_max)
    return;

  edges = _nc_xcor_merge_knots (knots1, knots2, k_min, k_max);

  if (edges->len < 2)
  {
    g_array_unref (edges);

    return;
  }

  W1 = W1_store;
  W2 = isauto ? W1_store : W2_store;

  /* The auto/cross distinction is fixed for the whole sweep, so it is resolved
   * once here rather than tested at every quadrature node. */
  if (isauto)
    _nc_xcor_gl5_sweep_auto (xclki1, edges, nell, W1, sum);
  else
    _nc_xcor_gl5_sweep_cross (xclki1, side2, edges, nell, W1, W2, sum);

  for (il = 0; il < nell; il++)
    ncm_vector_set (vp, il, const_factor * sum[il]);

  if (vp_err != NULL)
  {
    _nc_xcor_closure_err_init (&err_acc, xclki1, side2, isauto);
    _nc_xcor_closure_err_set_rows (&err_acc, xclki1, side2, isauto, edges, 1.0, 1.0);

    if (isauto)
      _nc_xcor_closure_err_sweep_auto (xclki1, edges, 0, nell, NULL, W1, &err_acc);
    else
      _nc_xcor_closure_err_sweep_cross (xclki1, side2, edges, 0, nell, NULL, W1, W2, &err_acc);

    _nc_xcor_closure_err_assemble (&err_acc, xclki1, side2, isauto, nell, const_factor, vp_err);
  }

  g_array_unref (edges);
}

/*
 * Aborts when NcXcor:reltol asks for more precision than the kernel's k-space
 * closure is built to, the looser of its Levin integrator's reltol and
 * cheb-reltol. The exact quadrature carries no tolerance of its own, so this is
 * a guard on the caller's stated intent: a C_ell cannot carry more precision
 * than its integrand.
 */
void
_nc_xcor_check_kernel_tolerance (NcXcor *xc, NcXcorKernel *xclk)
{
  NcmSBesselIntegrator *sbi = nc_xcor_kernel_peek_integrator (xclk);
  gdouble closure_reltol;

  if ((sbi == NULL) || !NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbi))
    return;

  {
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

    closure_reltol = GSL_MAX (ncm_sbessel_integrator_levin_get_reltol (sbilv),
                              ncm_sbessel_integrator_levin_get_cheb_reltol (sbilv));
  }

  if (xc->reltol < closure_reltol)
    g_error ("_nc_xcor_check_kernel_tolerance: NcXcor:reltol is %.17g but kernel %s builds "
             "its k-space closure only to %.17g (the looser of the integrator's reltol and "
             "cheb-reltol). The outer integral cannot converge to more precision than the "
             "integrand carries. Loosen NcXcor:reltol to at least %.17g, or construct the "
             "integrator with tighter tolerances.",
             xc->reltol, G_OBJECT_TYPE_NAME (xclk), closure_reltol, closure_reltol);
}

/*
 * The kernel-space methods, one entry each. Adding a quadrature is a line here
 * plus its block function; nothing else selects on the method.
 */
static const NcXcorKQuad _nc_xcor_kquad_table[] = {
  { _nc_xcor_kernel_integrate_block_exact, TRUE, "NC_XCOR_METHOD_KERNEL_EXACT" },
};

const NcXcorKQuad *
_nc_xcor_kquad_for_method (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_EXACT:
      return &_nc_xcor_kquad_table[0];

    default:
      return NULL;
  }
}

void
_nc_xcor_kernel_space_run (NcXcor *xc, const NcXcorKQuad *kquad, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  NcmSBesselIntegrator *sbi1 = nc_xcor_kernel_peek_integrator (xclk1);
  NcmSBesselIntegrator *sbi2 = nc_xcor_kernel_peek_integrator (xclk2);
  const guint size           = lmax - lmin + 1;
  const guint block          = xc->ell_batch_size;
  guint i;

  _nc_xcor_check_kernel_tolerance (xc, xclk1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, xclk2);

  /* Either kernel's integrator serves a kernel that carries none of its own. */
  if (sbi1 == NULL)
    sbi1 = sbi2;

  if (sbi2 == NULL)
    sbi2 = sbi1;

  /* Batched by NcXcor:ell-batch-size: one k-space closure per kernel per
   * batch, built here and handed to the same per-block integrator
   * #NcXcorSolver drives with closures of its own. The batching is not merely
   * an optimization -- a single closure spanning more than
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK multipoles is a hard error in
   * nc_xcor_kernel_get_eval_vectorized_full(), so an unbatched sweep aborted
   * on any range wider than that. */
  for (i = 0; i < size; i += block)
  {
    const guint nells      = MIN (block, size - i);
    const guint block_lmin = lmin + i;
    const guint block_lmax = block_lmin + nells - 1;
    NcmVector *vp_i        = ncm_vector_get_subvector (vp, i, nells);
    NcmVector *vp_err_i    = ((vp_err != NULL) && kquad->has_err) ? ncm_vector_get_subvector (vp_err, i, nells) : NULL;
    NcXcorKernelIntegrand *xclki1;
    NcXcorKernelIntegrand *xclki2;

    xclki1 = nc_xcor_kernel_get_eval_vectorized_full (xclk1, cosmo, block_lmin, block_lmax, sbi1, xc->closure_type);
    xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized_full (xclk2, cosmo, block_lmin, block_lmax, sbi2, xc->closure_type);

    kquad->block (xc, xclki1, isauto ? xclki1 : xclki2, block_lmin, block_lmax, isauto, vp_i, vp_err_i);

    nc_xcor_kernel_integrand_unref (xclki1);

    if (xclki2 != NULL)
      nc_xcor_kernel_integrand_unref (xclki2);

    ncm_vector_free (vp_i);
    ncm_vector_clear (&vp_err_i);
  }
}

/**
 * nc_xcor_integrate_block:
 * @xc: a #NcXcor
 * @xclki1: a #NcXcorKernelIntegrand covering [@lmin, @lmax]
 * @xclki2: (nullable): the same for the second kernel, %NULL for an auto
 * @lmin: minimum multipole, matching @xclki1's (and @xclki2's) own range
 * @lmax: maximum multipole, matching @xclki1's (and @xclki2's) own range
 * @isauto: %TRUE for an auto spectrum, in which case @xclki2 is ignored
 * @meth: the quadrature to run, which need not be @xc's own #NcXcor:meth
 * @vp: a #NcmVector of length (@lmax - @lmin + 1), filled with the result
 * @vp_err: (nullable): a #NcmVector of the same length for the error
 * estimate, or %NULL
 *
 * Runs one kernel-space quadrature over closures the caller already holds,
 * instead of building them from kernels the way nc_xcor_compute() does.
 *
 * Everything else about a $C_\ell$ -- fitting $W_\ell(k)$, choosing the tier,
 * batching the multipoles -- is an order of magnitude more expensive than the
 * outer integral, so a timing taken through nc_xcor_compute() is a timing of
 * the closure. This is the entry point that separates them: build the pair
 * once with nc_xcor_kernel_get_eval_vectorized_full(), then integrate it here,
 * as #NcXcorSolver does with closures it shares across requests.
 *
 * @meth is given here rather than read from @xc for the same reason. It must
 * be the kernel-space method, %NC_XCOR_METHOD_KERNEL_EXACT. The redshift-space
 * Limber methods have no block quadrature, being a different approximation
 * rather than a different quadrature.
 *
 * @xc still supplies #NcXcor:reltol, #NcXcor:closure-type and the $2/(\pi
 * R_H^3)$ factor, so nc_xcor_prepare() must have been called for the cosmology
 * the closures were built at.
 *
 * @vp_err is filled only by the methods nc_xcor_method_has_error_estimate()
 * reports; the rest fill it with NaN, as nc_xcor_compute_full() does, since a
 * zero there would read as "no error". Ask that question before reading it.
 */
void
nc_xcor_integrate_block (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcXcorMethod meth, NcmVector *vp, NcmVector *vp_err)
{
  const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (meth);

  g_return_if_fail (NC_IS_XCOR (xc));
  g_return_if_fail (xclki1 != NULL);
  g_return_if_fail (isauto || (xclki2 != NULL));
  g_return_if_fail (lmax >= lmin);
  g_return_if_fail (vp != NULL);

  if (kquad == NULL)
    g_error ("nc_xcor_integrate_block: %s has no block quadrature; it is not a "
             "kernel-space method.",
             nc_xcor_method_get_name (meth));

  if ((vp_err != NULL) && !kquad->has_err)
    ncm_vector_set_all (vp_err, GSL_NAN);

  kquad->block (xc, xclki1, isauto ? xclki1 : xclki2, lmin, lmax, isauto, vp, kquad->has_err ? vp_err : NULL);
}

/**
 * nc_xcor_method_get_name:
 * @meth: a #NcXcorMethod
 *
 * Returns: (transfer none): @meth's enum name, for labelling and error messages.
 */
const gchar *
nc_xcor_method_get_name (NcXcorMethod meth)
{
  const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (meth);

  if (kquad != NULL)
    return kquad->name;

  switch (meth)
  {
    case NC_XCOR_METHOD_LIMBER_Z_GSL:
      return "NC_XCOR_METHOD_LIMBER_Z_GSL";

    case NC_XCOR_METHOD_LIMBER_Z_CUBATURE:
      return "NC_XCOR_METHOD_LIMBER_Z_CUBATURE";

    default:
      return "NC_XCOR_METHOD_INVALID";
  }
}

/**
 * nc_xcor_method_has_error_estimate:
 * @meth: a #NcXcorMethod
 *
 * Whether @meth fills the @vp_err of nc_xcor_compute_full() and
 * nc_xcor_integrate_block(). Every other method leaves it alone rather than
 * writing a zero that would read as "no error", so a caller that wants an
 * estimate has to ask this first.
 *
 * Returns: %TRUE when @meth reports an error estimate.
 */
gboolean
nc_xcor_method_has_error_estimate (NcXcorMethod meth)
{
  const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (meth);

  return (kquad != NULL) && kquad->has_err;
}

/**
 * nc_xcor_method_is_kernel_space:
 * @meth: a #NcXcorMethod
 *
 * Whether @meth integrates over $k$ with a pair of fitted closures, as opposed
 * to the redshift-space Limber tier.
 *
 * Returns: %TRUE for the kernel-space methods.
 */
gboolean
nc_xcor_method_is_kernel_space (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_EXACT:
      return TRUE;

    default:
      return FALSE;
  }
}

