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
 * The outer k quadrature: everything the kernel-space methods
 * (%NC_XCOR_METHOD_KERNEL_GSL, %NC_XCOR_METHOD_KERNEL_CUBATURE and
 * %NC_XCOR_METHOD_KERNEL_EXACT) do with a pair of k-space closures, from the
 * block integrators and their GL(5) sweep through the knot merge, the qagp
 * breakpoints and the cubature range splitting.
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
#include "ncm/integration/ncm_integral_nd.h"
#include "nc/xcor/nc_xcor.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "nc/xcor/nc_xcor_priv.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

static void nc_xcor_kernel_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_kernel_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_xcor_kernel_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_kernel_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_KERNEL_AUTO, NcXcorKernelAuto, nc_xcor_kernel_auto, nc_xcor_kernel_auto_dim, nc_xcor_kernel_auto_integ, NcXcorArg);
NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_KERNEL_CROSS, NcXcorKernelCross, nc_xcor_kernel_cross, nc_xcor_kernel_cross_dim, nc_xcor_kernel_cross_integ, NcXcorArg);

static void
nc_xcor_kernel_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorKernelAuto *xcor_kernel_auto = NC_XCOR_KERNEL_AUTO (intnd);
  NcXcorArg *xcor_arg                = &xcor_kernel_auto->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_kernel_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorKernelCross *xcor_kernel_cross = NC_XCOR_KERNEL_CROSS (intnd);
  NcXcorArg *xcor_arg                  = &xcor_kernel_cross->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_kernel_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorKernelAuto *xcor_kernel_auto = NC_XCOR_KERNEL_AUTO (intnd);
  NcXcorArg *xcor_arg                = &xcor_kernel_auto->data;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble lnk = ncm_vector_fast_get (x, i);
    const gdouble k   = exp (lnk);
    const gdouble k3  = gsl_pow_3 (k);
    guint j;

    /* Evaluate all multipoles at once using pre-computed integrand. The
     * outputs are the block's components starting at comp_offset: a run of
     * multipoles sharing one k range, not necessarily the whole block. */
    nc_xcor_kernel_integrand_eval_comps (xcor_arg->xclki1, k, xcor_arg->comp_offset, fdim, xcor_arg->W1);

    for (j = 0; j < fdim; j++)
    {
      const gdouble W1  = xcor_arg->W1[xcor_arg->comp_offset + j];
      const gdouble res = k3 * W1 * W1;

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

static void
nc_xcor_kernel_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorKernelCross *xcor_kernel_cross = NC_XCOR_KERNEL_CROSS (intnd);
  NcXcorArg *xcor_arg                  = &xcor_kernel_cross->data;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble lnk = ncm_vector_fast_get (x, i);
    const gdouble k   = exp (lnk);
    const gdouble k3  = gsl_pow_3 (k);
    guint j;

    /* Evaluate all multipoles at once using pre-computed integrands. The
     * outputs are the block's components starting at comp_offset: a run of
     * multipoles sharing one k range, not necessarily the whole block. */
    nc_xcor_kernel_integrand_eval_comps (xcor_arg->xclki1, k, xcor_arg->comp_offset, fdim, xcor_arg->W1);
    nc_xcor_kernel_integrand_eval_comps (xcor_arg->xclki2, k, xcor_arg->comp_offset, fdim, xcor_arg->W2);

    for (j = 0; j < fdim; j++)
    {
      const gdouble res = k3 * xcor_arg->W1[xcor_arg->comp_offset + j] * xcor_arg->W2[xcor_arg->comp_offset + j];

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

static gdouble
_xcor_kernel_gsl_cross_int (gdouble lnk, gpointer ptr)
{
  NcXcorKernelIntegrand **xclki = (NcXcorKernelIntegrand **) ptr;
  const gdouble k               = exp (lnk);
  gdouble W1[1], W2[1];

  nc_xcor_kernel_integrand_eval (xclki[0], k, W1);
  nc_xcor_kernel_integrand_eval (xclki[1], k, W2);

  return gsl_pow_3 (k) * W1[0] * W2[0];
}

static gdouble
_xcor_kernel_gsl_auto_int (gdouble lnk, gpointer ptr)
{
  NcXcorKernelIntegrand **xclki = (NcXcorKernelIntegrand **) ptr;
  const gdouble k               = exp (lnk);
  gdouble W[1];

  nc_xcor_kernel_integrand_eval (xclki[0], k, W);

  return gsl_pow_3 (k) * W[0] * W[0];
}

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
 * nc_xcor_kernel_integrand_peek_residuals() uses that residual and
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
_nc_xcor_closure_cell_rows (NcXcorKernelIntegrand *xclki, GArray *edges)
{
  const guint n_panels = nc_xcor_kernel_integrand_get_n_panels (xclki);
  GArray *rows         = g_array_sized_new (FALSE, FALSE, sizeof (guint), edges->len);
  guint ie, j = 0;

  if (n_panels > 0)
  {
    for (ie = 0; ie + 1 < edges->len; ie++)
    {
      const gdouble cell_lo = g_array_index (edges, gdouble, ie);

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
      const gdouble cell_lo = g_array_index (edges, gdouble, ie);

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
_nc_xcor_closure_err_init (NcXcorClosureErr *err, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gboolean isauto, GArray *edges)
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

  err->residuals1 = nc_xcor_kernel_integrand_peek_residuals (xclki1);
  err->residuals2 = isauto ? err->residuals1 : nc_xcor_kernel_integrand_peek_residuals (xclki2);
  err->rows1      = (err->residuals1 != NULL) ? _nc_xcor_closure_cell_rows (xclki1, edges) : NULL;
  err->rows2      = isauto ? err->rows1 :
                    ((err->residuals2 != NULL) ? _nc_xcor_closure_cell_rows (xclki2, edges) : NULL);
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
_nc_xcor_closure_err_sweep_auto (NcXcorKernelIntegrand *xclki, GArray *edges, guint nell, gdouble *W, NcXcorClosureErr *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble cell_lo = g_array_index (edges, gdouble, ie);
    const gdouble cell_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid     = 0.5 * (cell_lo + cell_hi);
    const gdouble half    = 0.5 * (cell_hi - cell_lo);

    _nc_xcor_closure_cell_residual (err->residuals1, err->rows1, ie, nell, err->dW1, err->m1);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki, k, W);

      for (il = 0; il < nell; il++)
      {
        const gdouble absW = fabs (W[il]);

        /* d(W^2) = 2 |W| dW, and the aliased unk2/prod2/peak2 supply the
         * second half of the unknown-cell term in the assembly. */
        err->res[il]   += 2.0 * w * absW * err->dW1[il];
        err->unk1[il]  += w * absW * err->m1[il];
        err->prod1[il] += fabs (w * W[il] * W[il]) * err->m1[il];
        err->peak1[il]  = GSL_MAX (err->peak1[il], absW);
      }
    }
  }
}

static void
_nc_xcor_closure_err_sweep_cross (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, GArray *edges, guint nell, gdouble *W1, gdouble *W2, NcXcorClosureErr *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble cell_lo = g_array_index (edges, gdouble, ie);
    const gdouble cell_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid     = 0.5 * (cell_lo + cell_hi);
    const gdouble half    = 0.5 * (cell_hi - cell_lo);

    _nc_xcor_closure_cell_residual (err->residuals1, err->rows1, ie, nell, err->dW1, err->m1);
    _nc_xcor_closure_cell_residual (err->residuals2, err->rows2, ie, nell, err->dW2, err->m2);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki1, k, W1);
      nc_xcor_kernel_integrand_eval (xclki2, k, W2);

      for (il = 0; il < nell; il++)
      {
        const gdouble term  = w * W1[il] * W2[il];
        const gdouble absW1 = fabs (W1[il]);
        const gdouble absW2 = fabs (W2[il]);

        err->res[il]   += w * (absW1 * err->dW2[il] + absW2 * err->dW1[il]);
        err->unk1[il]  += w * absW2 * err->m1[il];
        err->unk2[il]  += w * absW1 * err->m2[il];
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
 * #NcXcorKernel:track-fit-residual off it is the whole estimate, and is then
 * exactly the tolerance-only bound.
 */
static void
_nc_xcor_closure_err_assemble (NcXcorClosureErr *err, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gboolean isauto, guint nell, gdouble const_factor, NcmVector *vp_err)
{
  const gdouble reltol1 = nc_xcor_kernel_integrand_get_reltol (xclki1);
  const gdouble reltol2 = nc_xcor_kernel_integrand_get_reltol (xclki2);
  const gdouble sabs1   = nc_xcor_kernel_integrand_get_scaled_abstol (xclki1);
  const gdouble sabs2   = nc_xcor_kernel_integrand_get_scaled_abstol (xclki2);
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
_nc_xcor_append_breakpoints (NcXcorKernelIntegrand *xclki, gdouble k_min, gdouble k_max, GArray *edges)
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

      if ((b > k_min) && (b < k_max))
        g_array_append_val (edges, b);
    }
  }
  else
  {
    NcmVector *knots   = nc_xcor_kernel_integrand_peek_knots (xclki);
    const guint nknots = (knots != NULL) ? ncm_vector_len (knots) : 0;

    for (i = 0; i < nknots; i++)
    {
      const gdouble knot = ncm_vector_get (knots, i);

      if ((knot > k_min) && (knot < k_max))
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
                            gboolean isauto, gdouble k_min, gdouble k_max)
{
  GArray *edges = g_array_new (FALSE, FALSE, sizeof (gdouble));

  g_array_append_val (edges, k_min);

  _nc_xcor_append_breakpoints (xclki1, k_min, k_max, edges);

  if (!isauto)
    _nc_xcor_append_breakpoints (xclki2, k_min, k_max, edges);

  g_array_sort (edges, _nc_xcor_cmp_edge);
  g_array_append_val (edges, k_max);

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
 * INT_{-1}^{1} T_i T_j dt, from T_i T_j = (T_{i+j} + T_{|i-j|}) / 2 and
 * INT T_n dt = 2 / (1 - n^2) for even n, zero for odd.
 */
static inline gdouble
_nc_xcor_cheb_TT_integral (const guint i, const guint j)
{
  const guint sum  = i + j;
  const guint diff = (i > j) ? i - j : j - i;

  if (sum % 2 != 0)
    return 0.0;

  return 1.0 / (1.0 - (gdouble) (sum * sum)) + 1.0 / (1.0 - (gdouble) (diff * diff));
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
static NcmMatrix *
_nc_xcor_cell_coeffs (NcXcorKernelIntegrand *xclki, gdouble a, gdouble b, guint nell, const gchar *side)
{
  NcmMatrix *coeffs = NULL;

  if (nc_xcor_kernel_integrand_get_n_panels (xclki) > 0)
  {
    /* Every cell of the common refinement lies inside one panel -- its edges
     * are the panels' own, so the containment test in restrict() compares
     * identical doubles. A failure here would mean the refinement and the
     * panels disagree, and dropping the cell would return a quietly wrong
     * C_ell. */
    if (!nc_xcor_kernel_integrand_restrict (xclki, a, b, &coeffs))
      g_error ("_nc_xcor_kernel_integrate_block_spectral: cell [%.17g, %.17g] "
               "is not inside a single panel of the %s closure.", a, b, side);

    return coeffs;
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

    coeffs = ncm_matrix_new (nell, 4);

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

    return coeffs;
  }
}

/*
 * Exact outer integral on the common refinement of two closures' breakpoints.
 *
 * On each cell both closures are a single polynomial over the same interval, so
 * k^2 W_i W_j is a polynomial there and its integral is a fixed bilinear form
 * in the two coefficient sets -- no nodes, no adaptivity, no tolerance. The
 * k^2 weight is itself degree two in the cell's own variable and is folded into
 * one side.
 *
 * This is what a Chebyshev closure gains over feeding it to an adaptive rule:
 * the quadrature does not rediscover per pair what the closure already
 * describes. Either side may equally be a spline, at four coefficients per
 * cell, which is what a pair split across the Limber threshold produces.
 */
static void
_nc_xcor_kernel_integrate_block_spectral (NcXcor *xc, NcXcorKernelIntegrand *xclki1,
                                          NcXcorKernelIntegrand *xclki2, guint lmin,
                                          guint lmax, gboolean isauto, NcmVector *vp,
                                          NcmVector *vp_err)
{
  const guint nell           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  gdouble k_min1, k_max1, k_min2, k_max2, k_min, k_max;
  NcXcorClosureErr err_acc;
  GArray *edges;

  /* One accumulator per multipole; a block is capped at
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK by the closure builder. */
  gdouble sum[NC_XCOR_KERNEL_MAX_ELL_BLOCK] = { 0.0 };
  gdouble W1_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble W2_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  GArray *folded_a;
  guint ie, il;

  nc_xcor_kernel_integrand_get_range (xclki1, &k_min1, &k_max1);
  nc_xcor_kernel_integrand_get_range (xclki2, &k_min2, &k_max2);

  k_min = GSL_MAX (k_min1, k_min2);
  k_max = GSL_MIN (k_max1, k_max2);

  ncm_vector_set_zero (vp);

  if (vp_err != NULL)
    ncm_vector_set_zero (vp_err);

  if (k_min >= k_max)
    return;

  edges = _nc_xcor_merge_panel_edges (xclki1, xclki2, isauto, k_min, k_max);

  /* One buffer for every cell, grown to fit. Its length is a Chebyshev
   * coefficient count, so it is the one piece of scratch here with no bound to
   * size a fixed array by; it is local rather than kept on @xc because the
   * solver's OpenMP team shares @xc across threads. */
  folded_a = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 0);

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble a    = g_array_index (edges, gdouble, ie);
    const gdouble b    = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid  = 0.5 * (a + b);
    const gdouble half = 0.5 * (b - a);
    NcmMatrix *c1      = _nc_xcor_cell_coeffs (xclki1, a, b, nell, "first");
    NcmMatrix *c2      = isauto ? ncm_matrix_ref (c1) :
                         _nc_xcor_cell_coeffs (xclki2, a, b, nell, "second");

    {
      /* k^2 in the cell's variable: k = mid + half t. */
      const gdouble w0 = mid * mid + 0.5 * half * half;
      const gdouble w1 = 2.0 * mid * half;
      const gdouble w2 = 0.5 * half * half;
      const guint n1   = ncm_matrix_ncols (c1);
      const guint n2   = ncm_matrix_ncols (c2);
      gdouble *folded;

      g_array_set_size (folded_a, n2 + 2);
      folded = &g_array_index (folded_a, gdouble, 0);

      for (il = 0; il < nell; il++)
      {
        guint i, j;

        /* Fold k^2 into the second factor, once per multipole. */
        for (j = 0; j < n2 + 2; j++)
          folded[j] = 0.0;

        for (j = 0; j < n2; j++)
        {
          const gdouble bj = ncm_matrix_get (c2, il, j);

          folded[j]     += w0 * bj;
          folded[j + 1] += 0.5 * w1 * bj;
          folded[j + 2] += 0.5 * w2 * bj;

          if (j >= 1)
            folded[j - 1] += 0.5 * w1 * bj;
          else
            folded[1] += 0.5 * w1 * bj;

          if (j >= 2)
            folded[j - 2] += 0.5 * w2 * bj;
          else
            folded[2 - j] += 0.5 * w2 * bj;
        }

        for (i = 0; i < n1; i++)
        {
          const gdouble ai = ncm_matrix_get (c1, il, i);
          gdouble acc      = 0.0;

          if (ai == 0.0)
            continue;

          for (j = 0; j < n2 + 2; j++)
            acc += folded[j] * _nc_xcor_cheb_TT_integral (i, j);

          sum[il] += half * ai * acc;
        }
      }
    }

    ncm_matrix_clear (&c1);
    ncm_matrix_clear (&c2);
  }

  for (il = 0; il < nell; il++)
    ncm_vector_set (vp, il, const_factor * sum[il]);

  /* The same estimate the merged-knot path reports, on the same cells: the
   * bilinear form is exact, so what is left is the closures' own fit error. */
  if (vp_err != NULL)
  {
    _nc_xcor_closure_err_init (&err_acc, xclki1, xclki2, isauto, edges);

    if (isauto)
      _nc_xcor_closure_err_sweep_auto (xclki1, edges, nell, W1_store, &err_acc);
    else
      _nc_xcor_closure_err_sweep_cross (xclki1, xclki2, edges, nell, W1_store, W2_store, &err_acc);

    _nc_xcor_closure_err_assemble (&err_acc, xclki1, xclki2, isauto, nell, const_factor, vp_err);
  }

  g_array_unref (folded_a);
  g_array_unref (edges);
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
    _nc_xcor_closure_err_init (&err_acc, xclki1, side2, isauto, edges);

    if (isauto)
      _nc_xcor_closure_err_sweep_auto (xclki1, edges, nell, W1, &err_acc);
    else
      _nc_xcor_closure_err_sweep_cross (xclki1, side2, edges, nell, W1, W2, &err_acc);

    _nc_xcor_closure_err_assemble (&err_acc, xclki1, side2, isauto, nell, const_factor, vp_err);
  }

  g_array_unref (edges);
}

/*
 * Subinterval breakpoints in ln k for one pair, from the integrands' own
 * knots, with the integration limits as the first and last entry. %NULL when
 * either integrand is not spline-backed and so reports no knots.
 */
static GArray *
_nc_xcor_kernel_gsl_breakpoints (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gdouble k_min, gdouble k_max)
{
  NcmVector *knots1 = nc_xcor_kernel_integrand_peek_knots (xclki1);
  NcmVector *knots2 = (xclki2 == NULL) ? knots1 : nc_xcor_kernel_integrand_peek_knots (xclki2);
  GArray *edges;
  GArray *pts;
  guint i;

  if ((knots1 == NULL) || (knots2 == NULL))
    return NULL;

  edges = _nc_xcor_merge_knots (knots1, knots2, k_min, k_max);
  pts   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), edges->len + 2);

  {
    const gdouble ln_k_min = log (k_min);

    g_array_append_val (pts, ln_k_min);
  }

  for (i = 0; i < edges->len; i++)
  {
    const gdouble ln_k = log (g_array_index (edges, gdouble, i));

    if (ln_k > g_array_index (pts, gdouble, pts->len - 1))
      g_array_append_val (pts, ln_k);
  }

  {
    const gdouble ln_k_max = log (k_max);

    if (ln_k_max > g_array_index (pts, gdouble, pts->len - 1))
      g_array_append_val (pts, ln_k_max);
    else
      g_array_index (pts, gdouble, pts->len - 1) = ln_k_max;
  }

  g_array_unref (edges);

  /* Degenerate merge (a domain narrower than one panel): nothing to break on. */
  if ((pts->len < 2) || (g_array_index (pts, gdouble, pts->len - 1) <= g_array_index (pts, gdouble, pts->len - 2)))
  {
    g_array_unref (pts);

    return NULL;
  }

  return pts;
}

void
_nc_xcor_kernel_gsl (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint nell              = ncm_vector_len (vp);
  const gdouble const_factor    = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  NcXcorKernelIntegrand *xclki_array[2];
  gsl_function F;
  guint i;
  gint ret;

  if (nell != lmax - lmin + 1)
    g_error ("_nc_xcor_kernel_gsl: vector size does not match multipole limits");

  if (lmax < lmin)
    g_error ("_nc_xcor_kernel_gsl: lmax < lmin");

  if (isauto)
    F.function = &_xcor_kernel_gsl_auto_int;
  else
    F.function = &_xcor_kernel_gsl_cross_int;

  F.params = xclki_array;

  for (i = 0; i < nell; i++)
  {
    const guint ell = lmin + i;
    gdouble k_min, k_max, result, err;
    GArray *breakpoints;

    /* Build the integrand(s) first, then read the outer bound off their own
     * fitted domain (get_range) -- NOT the independent Limber-band formula
     * from nc_xcor_kernel_get_k_range(), which has no guarantee of matching
     * it (see plan doc dev-notes/xcor_ultralevin_batching_plan.md sec. 3). */
    xclki_array[0] = nc_xcor_kernel_get_eval (xclk1, cosmo, ell, xc->closure_type);
    nc_xcor_kernel_integrand_get_range (xclki_array[0], &k_min, &k_max);

    if (isauto)
    {
      xclki_array[1] = NULL;
    }
    else
    {
      gdouble k2_min, k2_max;

      xclki_array[1] = nc_xcor_kernel_get_eval (xclk2, cosmo, ell, xc->closure_type);
      nc_xcor_kernel_integrand_get_range (xclki_array[1], &k2_min, &k2_max);

      k_min = GSL_MAX (k_min, k2_min);
      k_max = GSL_MIN (k_max, k2_max);
    }

    /* Integrated on the integrands' own knots when they have them, and over
     * the bare interval when they do not. A spline-backed integrand is
     * piecewise cubic: its third derivative jumps at every knot, which a
     * Gauss-Kronrod rule spanning several knots cannot see. Its error estimate
     * then saturates -- measured on the CMB ISW kernel, a panel holding four
     * or five knots reports about 80 times the error the same region reports
     * split on them -- and bisection stops improving it, which is exactly what
     * QUADPACK reports as roundoff (GSL_EROUND). Breaking on the knots leaves
     * every subinterval a single cubic piece: the same integral converges to
     * machine precision instead of stalling near 1e-6, for about twice the
     * evaluations. No safety margin is applied over NcXcor:reltol: a margin
     * of (reltol * 1e-2) was measured to move the stall rather than remove it.
     */
    breakpoints = _nc_xcor_kernel_gsl_breakpoints (xclki_array[0], isauto ? NULL : xclki_array[1], k_min, k_max);

    if (breakpoints != NULL)
      ret = gsl_integration_qagp (&F, (gdouble *) breakpoints->data, breakpoints->len, 0.0, xc->reltol, NCM_INTEGRAL_PARTITION, *w, &result, &err);
    else
      ret = gsl_integration_qag (&F, log (k_min), log (k_max), 0.0, xc->reltol, NCM_INTEGRAL_PARTITION, 6, *w, &result, &err);

    _nc_xcor_check_qag_status ("_nc_xcor_kernel_gsl", ret, xc->reltol, result, err);

    if (breakpoints != NULL)
      g_array_unref (breakpoints);

    ncm_vector_set (vp, i, const_factor * result);

    nc_xcor_kernel_integrand_unref (xclki_array[0]);

    if (!isauto)
      nc_xcor_kernel_integrand_unref (xclki_array[1]);
  }

  ncm_memory_pool_return (w);
}

/*
 * One multipole of a block closure, in ln k. @il selects which, so the whole
 * block shares one closure and one evaluation buffer while qagp is driven once
 * per multipole -- the difference between %NC_XCOR_METHOD_KERNEL_GSL_BLOCK and
 * %NC_XCOR_METHOD_KERNEL_GSL is entirely in that sharing, not in the rule.
 */
typedef struct _NcXcorGSLBlockArg
{
  NcXcorKernelIntegrand *xclki1;
  NcXcorKernelIntegrand *xclki2; /* %NULL for an auto-correlation */
  gdouble *W1;
  gdouble *W2;
  guint il;
} NcXcorGSLBlockArg;

static gdouble
_xcor_kernel_gsl_block_auto_int (gdouble lnk, gpointer ptr)
{
  NcXcorGSLBlockArg *arg = (NcXcorGSLBlockArg *) ptr;
  const gdouble k        = exp (lnk);

  nc_xcor_kernel_integrand_eval (arg->xclki1, k, arg->W1);

  return gsl_pow_3 (k) * arg->W1[arg->il] * arg->W1[arg->il];
}

static gdouble
_xcor_kernel_gsl_block_cross_int (gdouble lnk, gpointer ptr)
{
  NcXcorGSLBlockArg *arg = (NcXcorGSLBlockArg *) ptr;
  const gdouble k        = exp (lnk);

  nc_xcor_kernel_integrand_eval (arg->xclki1, k, arg->W1);
  nc_xcor_kernel_integrand_eval (arg->xclki2, k, arg->W2);

  return gsl_pow_3 (k) * arg->W1[arg->il] * arg->W2[arg->il];
}

/*
 * _nc_xcor_kernel_integrate_block_gsl:
 * @xc: a #NcXcor
 * @xclki1: a pre-built #NcXcorKernelIntegrand covering [@lmin, @lmax]
 * @xclki2: (nullable): the same for the second kernel, %NULL for an auto
 * @lmin: minimum multipole, matching @xclki1's own range
 * @lmax: maximum multipole, matching @xclki1's own range
 * @isauto: %TRUE for an auto-correlation (only @xclki1 is used)
 * @vp: (out): output vector of length (@lmax - @lmin + 1)
 * @vp_err: unused -- qagp's own estimate is a statement about the quadrature,
 * which is not what nc_xcor_compute_full() reports; see there
 *
 * %NC_XCOR_METHOD_KERNEL_GSL_BLOCK's block quadrature: QUADPACK's qagp broken
 * on the merged knots, exactly the rule %NC_XCOR_METHOD_KERNEL_GSL runs, over
 * the block closure %NC_XCOR_METHOD_KERNEL_EXACT and
 * %NC_XCOR_METHOD_KERNEL_CUBATURE integrate. Its reason to exist is that
 * comparison: on one shared pair of closures it isolates the outer quadrature,
 * which no pair of methods could do while one of them fitted its own closure
 * per multipole.
 *
 * Each multipole is integrated over its component range because the Limber
 * integrand vanishes outside that range.
 *
 * This is a diagnostic method, not a production one. On spline closures it
 * costs 5.1 ms against %NC_XCOR_METHOD_KERNEL_EXACT's 0.13 ms at the same
 * accuracy: QUADPACK is scalar, so it calls a block closure once per
 * multipole and discards the rest of the block, and it subdivides where
 * GL(5) is already exact.
 *
 * The caller retains ownership of @xclki1/@xclki2.
 */
void
_nc_xcor_kernel_integrate_block_gsl (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell              = lmax - lmin + 1;
  const gdouble const_factor    = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  /* One double per multipole, and a block is capped at
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK by the closure builder. */
  gdouble W1_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble W2_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  NcXcorGSLBlockArg arg;
  gsl_function F;
  guint il;

  NCM_UNUSED (vp_err);

  if (ncm_vector_len (vp) != nell)
    g_error ("_nc_xcor_kernel_integrate_block_gsl: vector size does not match multipole limits");

  if (nell > NC_XCOR_KERNEL_MAX_ELL_BLOCK)
    g_error ("_nc_xcor_kernel_integrate_block_gsl: block of %u multipoles exceeds "
             "NC_XCOR_KERNEL_MAX_ELL_BLOCK (%d).", nell, NC_XCOR_KERNEL_MAX_ELL_BLOCK);

  arg.xclki1 = xclki1;
  arg.xclki2 = isauto ? NULL : xclki2;
  arg.W1     = W1_store;
  arg.W2     = isauto ? NULL : W2_store;
  arg.il     = 0;

  F.function = isauto ? &_xcor_kernel_gsl_block_auto_int : &_xcor_kernel_gsl_block_cross_int;
  F.params   = &arg;

  for (il = 0; il < nell; il++)
  {
    gdouble k_min, k_max, result, err;
    GArray *breakpoints;
    gint ret;

    nc_xcor_kernel_integrand_get_range_comp (xclki1, il, &k_min, &k_max);

    if (!isauto)
    {
      gdouble k2_min, k2_max;

      nc_xcor_kernel_integrand_get_range_comp (xclki2, il, &k2_min, &k2_max);

      k_min = GSL_MAX (k_min, k2_min);
      k_max = GSL_MIN (k_max, k2_max);
    }

    if (k_max <= k_min)
    {
      ncm_vector_set (vp, il, 0.0);
      continue;
    }

    arg.il      = il;
    breakpoints = _nc_xcor_kernel_gsl_breakpoints (xclki1, isauto ? NULL : xclki2, k_min, k_max);

    if (breakpoints != NULL)
      ret = gsl_integration_qagp (&F, (gdouble *) breakpoints->data, breakpoints->len, 0.0, xc->reltol, NCM_INTEGRAL_PARTITION, *w, &result, &err);
    else
      ret = gsl_integration_qag (&F, log (k_min), log (k_max), 0.0, xc->reltol, NCM_INTEGRAL_PARTITION, 6, *w, &result, &err);

    _nc_xcor_check_qag_status ("_nc_xcor_kernel_integrate_block_gsl", ret, xc->reltol, result, err);

    if (breakpoints != NULL)
      g_array_unref (breakpoints);

    ncm_vector_set (vp, il, const_factor * result);
  }

  ncm_memory_pool_return (w);
}

/*
 * The outer k-integral cannot resolve the closure more finely than the closure
 * itself is built. pcubature answers an impossible tolerance by exhausting its
 * Clenshaw-Curtis levels and reporting failure, far from the cause, so catch
 * the mismatch here where both numbers are in view.
 *
 * That reasoning is about %NC_XCOR_METHOD_KERNEL_CUBATURE. This is also called
 * from the exact path, where NcXcor:reltol governs nothing -- GL(5) on the knot
 * union carries no tolerance and cannot fail to converge -- so there the check
 * is a consistency guard on the caller's stated intent rather than a failure it
 * is preventing.
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
 * The k range one multipole of a pair is supported on: both integrands' own
 * ranges for that component, intersected.
 */
static void
_nc_xcor_kernel_cubature_comp_range (NcXcorArg *xcor_arg, gboolean isauto, guint i, gdouble *k_min, gdouble *k_max)
{
  nc_xcor_kernel_integrand_get_range_comp (xcor_arg->xclki1, i, k_min, k_max);

  if (!isauto)
  {
    gdouble k2_min, k2_max;

    nc_xcor_kernel_integrand_get_range_comp (xcor_arg->xclki2, i, &k2_min, &k2_max);

    *k_min = GSL_MAX (*k_min, k2_min);
    *k_max = GSL_MIN (*k_max, k2_max);
  }
}

/*
 * One ell block, integrated in runs of consecutive multipoles that share a k
 * range.
 *
 * A block's multipoles share one fitted domain but not necessarily one
 * support: under the Limber approximation each is confined to its own band in
 * k and vanishes outside it, so integrating the whole block over the shared
 * domain hands the quadrature a step per multipole. An adaptive rule does not
 * merely work harder there, it can miss the feature outright -- a region whose
 * nodes all land where the component vanishes reports no error and is never
 * subdivided -- and a kernel whose window does not tend to zero at its own
 * band edge (CMB ISW, whose radial integral truncates at recombination) then
 * loses a finite part of the integral under a converged status. Splitting on
 * the range puts each step on an integration limit, where it is not a
 * discontinuity of the integrand at all.
 *
 * Where every multipole shares one range -- the non-Limber case, whose support
 * is the fitted domain itself -- this is one call over the whole block, as
 * before.
 */
static void
_nc_xcor_kernel_cubature_runs (NcmIntegralND *xcor_int_nd, NcXcorArg *xcor_arg, gboolean isauto, guint nells, NcmVector *vp_i, NcmVector *err_i)
{
  NcmVector *lnk_min = ncm_vector_new (1);
  NcmVector *lnk_max = ncm_vector_new (1);
  guint j            = 0;

  while (j < nells)
  {
    gdouble k_min, k_max;
    guint run = 1;

    _nc_xcor_kernel_cubature_comp_range (xcor_arg, isauto, j, &k_min, &k_max);

    while (j + run < nells)
    {
      gdouble k_min_run, k_max_run;

      _nc_xcor_kernel_cubature_comp_range (xcor_arg, isauto, j + run, &k_min_run, &k_max_run);

      if ((k_min_run != k_min) || (k_max_run != k_max))
        break;

      run++;
    }

    if (k_min < k_max)
    {
      NcmVector *vp_j  = ncm_vector_get_subvector (vp_i, j, run);
      NcmVector *err_j = ncm_vector_get_subvector (err_i, j, run);

      xcor_arg->comp_offset = j;
      xcor_arg->nells       = run;

      ncm_vector_set (lnk_min, 0, log (k_min));
      ncm_vector_set (lnk_max, 0, log (k_max));

      ncm_integral_nd_eval (xcor_int_nd, lnk_min, lnk_max, vp_j, err_j);

      ncm_vector_free (vp_j);
      ncm_vector_free (err_j);
    }
    else
    {
      /* Empty support: the pair's bands do not meet inside the domain. */
      guint m;

      for (m = 0; m < run; m++)
      {
        ncm_vector_set (vp_i, j + m, 0.0);
        ncm_vector_set (err_i, j + m, 0.0);
      }
    }

    j += run;
  }

  ncm_vector_free (lnk_min);
  ncm_vector_free (lnk_max);
}

/*
 * _nc_xcor_kernel_integrate_block_cubature:
 * @xc: a #NcXcor
 * @xclki1: a pre-built #NcXcorKernelIntegrand for the first kernel,
 * covering multipoles [@lmin, @lmax]
 * @xclki2: (nullable): a pre-built #NcXcorKernelIntegrand for the second
 * kernel, covering the same range, or %NULL for an auto-correlation
 * @lmin: minimum multipole, matching @xclki1's (and @xclki2's) own range
 * @lmax: maximum multipole, matching @xclki1's (and @xclki2's) own range
 * @isauto: %TRUE for an auto-correlation (only @xclki1 is used)
 * @vp: (out): output vector of length (@lmax - @lmin + 1)
 * @vp_err: unused -- pcubature's estimate is a quadrature error, which is not
 * what nc_xcor_compute_full() reports; see there
 *
 * %NC_XCOR_METHOD_KERNEL_CUBATURE's block quadrature: adaptive pcubature over
 * integrand(s) the caller already built, rather than ones built internally.
 * This is #NcXcorSolver's building block for sharing one kernel's k-space
 * closure across every request that needs it in a given block, instead of
 * rebuilding it once per pair (see plan doc
 * dev-notes/xcor_ultralevin_batching_plan.md sec. 5-6). Reached through the
 * #NcXcorKQuad table, publicly through nc_xcor_integrate_block().
 *
 * The caller retains ownership of @xclki1/@xclki2 -- they are not unreffed
 * here.
 */
void
_nc_xcor_kernel_integrate_block_cubature (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const guint size           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcmVector *err             = ncm_vector_new (size);
  gdouble W1_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  gdouble W2_store[NC_XCOR_KERNEL_MAX_ELL_BLOCK];
  NcmIntegralND *xcor_int_nd;
  NcXcorArg *xcor_arg;

  NCM_UNUSED (vp_err);

  if (ncm_vector_len (vp) != size)
    g_error ("_nc_xcor_kernel_integrate_block_cubature: vector size does not match multipole limits");

  if (size > NC_XCOR_KERNEL_MAX_ELL_BLOCK)
    g_error ("_nc_xcor_kernel_integrate_block_cubature: block of %u multipoles exceeds "
             "NC_XCOR_KERNEL_MAX_ELL_BLOCK (%d).", size, NC_XCOR_KERNEL_MAX_ELL_BLOCK);

  if (isauto)
  {
    NcXcorKernelAuto *xcor_kernel_auto = g_object_new (nc_xcor_kernel_auto_get_type (), NULL);

    xcor_int_nd = NCM_INTEGRAL_ND (xcor_kernel_auto);
    xcor_arg    = &xcor_kernel_auto->data;
  }
  else
  {
    NcXcorKernelCross *xcor_kernel_cross = g_object_new (nc_xcor_kernel_cross_get_type (), NULL);

    xcor_int_nd = NCM_INTEGRAL_ND (xcor_kernel_cross);
    xcor_arg    = &xcor_kernel_cross->data;
  }

  xcor_arg->xc     = xc;
  xcor_arg->RH     = xc->RH;
  xcor_arg->xclki1 = xclki1;

  /* The integrand object is created and unreffed inside this function, so the
   * scratch it points at only has to outlive the runs below. */
  xcor_arg->W1 = W1_store;

  if (!isauto)
  {
    xcor_arg->xclki2 = xclki2;
    xcor_arg->W2     = W2_store;
  }

  ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
  ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
  ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

  _nc_xcor_kernel_cubature_runs (xcor_int_nd, xcor_arg, isauto, size, vp, err);

  ncm_vector_scale (vp, const_factor);

  g_object_unref (xcor_int_nd);
  ncm_vector_free (err);
}

/*
 * The kernel-space methods, one entry each. Adding a fifth quadrature is a
 * line here plus its block function; nothing else selects on the method.
 *
 * %NC_XCOR_METHOD_KERNEL_GSL is absent on purpose: it fits a closure per
 * multipole rather than per block, so it has no block quadrature to name and
 * keeps its own runner. See #NcXcorMethod.
 */
static const NcXcorKQuad _nc_xcor_kquad_table[] = {
  { _nc_xcor_kernel_integrate_block_cubature, FALSE, FALSE, "NC_XCOR_METHOD_KERNEL_CUBATURE"  },
  { _nc_xcor_kernel_integrate_block_exact,    TRUE,  TRUE,  "NC_XCOR_METHOD_KERNEL_EXACT"     },
  { _nc_xcor_kernel_integrate_block_gsl,      TRUE,  FALSE, "NC_XCOR_METHOD_KERNEL_GSL_BLOCK" },
};

const NcXcorKQuad *
_nc_xcor_kquad_for_method (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_CUBATURE:
      return &_nc_xcor_kquad_table[0];

    case NC_XCOR_METHOD_KERNEL_EXACT:
      return &_nc_xcor_kquad_table[1];

    case NC_XCOR_METHOD_KERNEL_GSL_BLOCK:
      return &_nc_xcor_kquad_table[2];

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

    if (kquad->needs_integrator)
    {
      xclki1 = nc_xcor_kernel_get_eval_vectorized_full (xclk1, cosmo, block_lmin, block_lmax, sbi1, xc->closure_type);
      xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized_full (xclk2, cosmo, block_lmin, block_lmax, sbi2, xc->closure_type);
    }
    else
    {
      xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, block_lmin, block_lmax, xc->closure_type);
      xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized (xclk2, cosmo, block_lmin, block_lmax, xc->closure_type);
    }

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
 * once with nc_xcor_kernel_get_eval_vectorized_full(), then drive every method
 * over it. That also makes a method comparison exact rather than approximate,
 * since the two runs then share an integrand instead of each fitting its own.
 *
 * @meth is given here rather than read from @xc for the same reason. It must
 * be a kernel-space method with a block quadrature --
 * %NC_XCOR_METHOD_KERNEL_CUBATURE, %NC_XCOR_METHOD_KERNEL_EXACT or
 * %NC_XCOR_METHOD_KERNEL_GSL_BLOCK. %NC_XCOR_METHOD_KERNEL_GSL has none: it
 * fits its closure one multipole at a time, so there is nothing for a caller
 * to hand it. The redshift-space Limber methods have none either, being a
 * different approximation rather than a different quadrature.
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
             "kernel-space method, or it builds its closure per multipole.",
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

    case NC_XCOR_METHOD_KERNEL_GSL:
      return "NC_XCOR_METHOD_KERNEL_GSL";

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
 * to the redshift-space Limber tier. %NC_XCOR_METHOD_KERNEL_GSL is kernel-space
 * but has no block quadrature, so this is the wider of the two questions --
 * nc_xcor_integrate_block() answers the narrower one by refusing.
 *
 * Returns: %TRUE for the kernel-space methods.
 */
gboolean
nc_xcor_method_is_kernel_space (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_GSL:
    case NC_XCOR_METHOD_KERNEL_CUBATURE:
    case NC_XCOR_METHOD_KERNEL_EXACT:
    case NC_XCOR_METHOD_KERNEL_GSL_BLOCK:
      return TRUE;

    default:
      return FALSE;
  }
}

