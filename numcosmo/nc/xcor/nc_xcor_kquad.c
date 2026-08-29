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
 * Five-node Gauss-Legendre rule on [-1, 1]. Exact through degree 9, which is
 * one more than the degree-8 polynomial the outer integrand k^2 W_i W_j is on
 * each knot panel of a cubic spline: 2 from k^2 plus 3 from each spline.
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
 * The two GL(5) sweeps over a pair's merged panels. Auto and cross differ only
 * in whether the second integrand is evaluated at all, which is fixed for the
 * whole sweep -- so they are separate functions and the inner loops carry no
 * branch. Each is a flat loop over panel x node x multipole.
 */
/*
 * Everything the error estimate of _nc_xcor_kernel_integrate_block_exact ()
 * needs, accumulated in the same pass as the integral itself.
 *
 * The estimate propagates d(W1 W2) = |W1| dW2 + |W2| dW1 with dW_i the closure
 * fit's own error, so what it needs from the sweep is that product integrated
 * against k^2. Where the closure recorded the residual it *achieved* on a
 * panel (nc_xcor_kernel_integrand_peek_residuals()), dW_i is that residual and
 * the term lands in @res. Where it did not -- tracking off, or a panel
 * refinement never accepted -- the panel lands in @unk_i instead, to be closed
 * afterwards with the tolerance the fit was *asked* for, times the peak. The
 * masks are per panel, so the inner loop multiplies rather than branches.
 *
 * The peaks enter only as multipliers of the accumulated integrals, so a
 * running maximum is enough and no second sweep is required.
 */
typedef struct _NcXcorGL5Err
{
  gdouble *res;          /* int k^2 (|W1| dW2 + |W2| dW1), panels with a record  */
  gdouble *unk1;         /* int k^2 |W2|    over panels where W1 has no record   */
  gdouble *unk2;         /* int k^2 |W1|    over panels where W2 has no record   */
  gdouble *prod1;        /* int k^2 |W1 W2| over panels where W1 has no record   */
  gdouble *prod2;        /* int k^2 |W1 W2| over panels where W2 has no record   */
  gdouble *peak1;        /* max |W1|                                             */
  gdouble *peak2;        /* max |W2|                                             */
  NcmMatrix *residuals1; /* achieved residuals, or NULL                   */
  NcmMatrix *residuals2;
  GArray *rows1; /* panel -> row of @residuals1 (guint)                  */
  GArray *rows2;
  gdouble *dW1; /* per-panel scratch: residual, or 0 where unknown      */
  gdouble *dW2;
  gdouble *m1; /* per-panel scratch: 1.0 where unknown, else 0.0       */
  gdouble *m2;
} NcXcorGL5Err;

/*
 * Fills the per-panel dW/mask scratch for one side from its recorded
 * residuals. A NaN entry means the interval was never accepted, and is treated
 * exactly as no record at all.
 */
static void
_nc_xcor_gl5_panel_residual (NcmMatrix *residuals, GArray *rows, const guint ie, const guint nell, gdouble *dW, gdouble *m)
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
 * Maps each panel of @edges onto the interval of @knots that contains it. Both
 * are sorted and every edge is a knot of one side or the other, so a single
 * marching index does it.
 */
static GArray *
_nc_xcor_gl5_panel_rows (NcmVector *knots, GArray *edges)
{
  const guint nknots = ncm_vector_len (knots);
  GArray *rows       = g_array_sized_new (FALSE, FALSE, sizeof (guint), edges->len);
  guint ie, j = 0;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);

    while ((j + 2 < nknots) && (ncm_vector_get (knots, j + 1) <= panel_lo))
      j++;

    g_array_append_val (rows, j);
  }

  return rows;
}

static void
_nc_xcor_gl5_sweep_auto (NcXcorKernelIntegrand *xclki, GArray *edges, guint nell, gdouble *W, gdouble *sum, NcXcorGL5Err *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);
    const gdouble panel_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid      = 0.5 * (panel_lo + panel_hi);
    const gdouble half     = 0.5 * (panel_hi - panel_lo);

    if (err != NULL)
      _nc_xcor_gl5_panel_residual (err->residuals1, err->rows1, ie, nell, err->dW1, err->m1);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki, k, W);

      for (il = 0; il < nell; il++)
      {
        const gdouble term = w * W[il] * W[il];

        sum[il] += term;

        if (err != NULL)
        {
          const gdouble absW = fabs (W[il]);

          /* d(W^2) = 2 |W| dW, and the aliased unk2/peak2 supply the second
           * half of the unknown-panel term the same way. */
          err->res[il]   += 2.0 * w * absW * err->dW1[il];
          err->unk1[il]  += w * absW * err->m1[il];
          err->prod1[il] += fabs (term) * err->m1[il];
          err->peak1[il]  = GSL_MAX (err->peak1[il], absW);
        }
      }
    }
  }
}

static void
_nc_xcor_gl5_sweep_cross (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, GArray *edges, guint nell, gdouble *W1, gdouble *W2, gdouble *sum, NcXcorGL5Err *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);
    const gdouble panel_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid      = 0.5 * (panel_lo + panel_hi);
    const gdouble half     = 0.5 * (panel_hi - panel_lo);

    if (err != NULL)
    {
      _nc_xcor_gl5_panel_residual (err->residuals1, err->rows1, ie, nell, err->dW1, err->m1);
      _nc_xcor_gl5_panel_residual (err->residuals2, err->rows2, ie, nell, err->dW2, err->m2);
    }

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki1, k, W1);
      nc_xcor_kernel_integrand_eval (xclki2, k, W2);

      for (il = 0; il < nell; il++)
      {
        const gdouble term = w * W1[il] * W2[il];

        sum[il] += term;

        if (err != NULL)
        {
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
 * Exact outer quadrature for one pair, on the union of the two kernels' own
 * knot sets.
 *
 * Each kernel is sampled independently -- the same per-kernel closures
 * KERNEL_CUBATURE builds and NcXcorSolver caches -- so the two splines live on
 * different abscissas. On the common refinement of those abscissas each spline
 * is still a single cubic piece per panel, so k^2 W_i W_j is a degree-8
 * polynomial there and GL(5) integrates it exactly. Merging two knot sets is
 * all the coupling the exactness argument needs; a shared abscissa built by
 * sampling the kernels jointly is not required, and costs about twice as much
 * to produce.
 *
 * The range is the intersection of the two integrands' fitted domains, exactly
 * as _nc_xcor_kernel_integrate_block_cubature() uses, because NcmSpline does
 * not range-check and an out-of-domain evaluation returns extrapolation rather
 * than a small number.
 *
 * ## What this means for error control
 *
 * Exact is meant literally: refining every panel fourfold moves the result by
 * 1e-15 to 1e-12, which is rounding. So there is nothing for an embedded
 * quadrature rule to measure here -- a Kronrod extension, or GL(5) against
 * GL(9), would report machine zero on every call and never fire. Do not add
 * one; it would be false confidence rather than error control.
 *
 * The error is entirely in the closure: a spline is a fit, and $C_\ell$ can be
 * far smaller than the integral of the absolute integrand, which amplifies
 * that fit's error. Two disjoint Gaussian bins already reach a cancellation of
 * 1.4e4 by $\ell = 9$, so a closure good to 1e-8 leaves 1e-4 on $C_\ell$.
 * @vp_err reports that product; see nc_xcor_compute_full().
 *
 * The same C^2 kinks are why the adaptive kernel-space methods struggle on
 * this integrand and this one does not. A cubic spline's third derivative
 * jumps at every knot, so an adaptive scheme subdivides at each of them
 * forever, chasing a relative criterion that the fit's own error puts out of
 * reach. Here the knots *are* the panel edges, so the kinks fall between
 * panels and never enter a rule's smoothness assumption.
 */
static gint
_nc_xcor_cmp_edge (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/*
 * Merges two panel edge sets over [k_min, k_max]. The result is the common
 * refinement, on each cell of which both closures are a single polynomial --
 * the same argument the merged knot sets make for splines, one level up.
 */
static GArray *
_nc_xcor_merge_panel_edges (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2,
                            gboolean isauto, gdouble k_min, gdouble k_max)
{
  GArray *edges  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint n1 = nc_xcor_kernel_integrand_get_n_panels (xclki1);
  const guint n2 = isauto ? 0 : nc_xcor_kernel_integrand_get_n_panels (xclki2);
  guint i;

  g_array_append_val (edges, k_min);

  for (i = 0; i < n1 + n2; i++)
  {
    NcmMatrix *ignored = NULL;
    gdouble a, b, edge;

    if (i < n1)
      nc_xcor_kernel_integrand_peek_panel (xclki1, i, &ignored, &a, &b);
    else
      nc_xcor_kernel_integrand_peek_panel (xclki2, i - n1, &ignored, &a, &b);

    edge = b;

    if ((edge > k_min) && (edge < k_max))
      g_array_append_val (edges, edge);
  }

  g_array_sort (edges, _nc_xcor_cmp_edge);
  g_array_append_val (edges, k_max);

  /* Drop duplicates: two closures often break at the same place. */
  {
    guint w = 1;

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
 * Exact outer integral for a pair of spectral closures.
 *
 * On each cell of the common refinement of the two panel edge sets both
 * closures are a single polynomial over the same interval, so k^2 W_i W_j is a
 * polynomial there and its integral is a fixed bilinear form in the two
 * coefficient sets -- no nodes, no adaptivity, no tolerance. The k^2 weight is
 * itself degree two in the cell's own variable and is folded into one side.
 *
 * This is what a Chebyshev closure buys over feeding it to an adaptive rule:
 * the quadrature stops rediscovering per pair what the closure already knows.
 */
static void
_nc_xcor_kernel_integrate_block_spectral (NcXcor *xc, NcXcorKernelIntegrand *xclki1,
                                          NcXcorKernelIntegrand *xclki2, guint lmin,
                                          guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint nell           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  gdouble k_min1, k_max1, k_min2, k_max2, k_min, k_max;
  GArray *edges;
  gdouble *sum;
  guint ie, il;

  nc_xcor_kernel_integrand_get_range (xclki1, &k_min1, &k_max1);
  nc_xcor_kernel_integrand_get_range (isauto ? xclki1 : xclki2, &k_min2, &k_max2);

  k_min = GSL_MAX (k_min1, k_min2);
  k_max = GSL_MIN (k_max1, k_max2);

  ncm_vector_set_zero (vp);

  if (k_min >= k_max)
    return;

  edges = _nc_xcor_merge_panel_edges (xclki1, xclki2, isauto, k_min, k_max);
  sum   = g_new0 (gdouble, nell);

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble a    = g_array_index (edges, gdouble, ie);
    const gdouble b    = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid  = 0.5 * (a + b);
    const gdouble half = 0.5 * (b - a);
    NcmMatrix *c1      = NULL;
    NcmMatrix *c2      = NULL;

    /* Every cell of the common refinement lies inside one panel of each
     * closure -- its edges are the panels' own, so the containment test in
     * restrict() compares identical doubles. A failure here would mean the
     * refinement and the panels disagree, and dropping the cell would return a
     * quietly wrong C_ell. */
    if (!nc_xcor_kernel_integrand_restrict (xclki1, a, b, &c1))
      g_error ("_nc_xcor_kernel_integrate_block_spectral: cell [%.17g, %.17g] "
               "is not inside a single panel of the first closure.", a, b);

    if (isauto)
      c2 = ncm_matrix_ref (c1);
    else if (!nc_xcor_kernel_integrand_restrict (xclki2, a, b, &c2))
      g_error ("_nc_xcor_kernel_integrate_block_spectral: cell [%.17g, %.17g] "
               "is not inside a single panel of the second closure.", a, b);

    {
      /* k^2 in the cell's variable: k = mid + half t. */
      const gdouble w0 = mid * mid + 0.5 * half * half;
      const gdouble w1 = 2.0 * mid * half;
      const gdouble w2 = 0.5 * half * half;
      const guint n1   = ncm_matrix_ncols (c1);
      const guint n2   = ncm_matrix_ncols (c2);
      gdouble *folded  = g_new0 (gdouble, n2 + 2);

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

      g_free (folded);
    }

    ncm_matrix_clear (&c1);
    ncm_matrix_clear (&c2);
  }

  for (il = 0; il < nell; il++)
    ncm_vector_set (vp, il, const_factor * sum[il]);

  g_free (sum);
  g_array_unref (edges);
}

void
_nc_xcor_kernel_integrate_block_exact (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcmVector *knots1          = nc_xcor_kernel_integrand_peek_knots (xclki1);
  NcmVector *knots2          = isauto ? knots1 : nc_xcor_kernel_integrand_peek_knots (xclki2);
  const gdouble reltol1      = nc_xcor_kernel_integrand_get_reltol (xclki1);
  const gdouble reltol2      = nc_xcor_kernel_integrand_get_reltol (xclki2);
  const gdouble sabs1        = nc_xcor_kernel_integrand_get_scaled_abstol (xclki1);
  const gdouble sabs2        = nc_xcor_kernel_integrand_get_scaled_abstol (xclki2);
  gdouble k_min1, k_max1, k_min2, k_max2, k_min, k_max;
  NcXcorGL5Err err_acc, *err = NULL;
  gdouble *sum, *W1, *W2;
  GArray *edges;
  guint il;

  /* Chosen here rather than by the callers: NcXcorSolver and
   * _nc_xcor_kernel_space_compute() both enter through this function, and a
   * choice made in one of them is a choice the other silently misses. A
   * spectral pair goes to the common refinement of its panels, where the
   * product is a polynomial and the integral is exact in closed form; splines
   * take the merged-knot GL(5) sweep below. Either way the method integrates
   * the closures it is handed exactly. */
  if ((nc_xcor_kernel_integrand_get_n_panels (xclki1) > 0) &&
      (nc_xcor_kernel_integrand_get_n_panels (isauto ? xclki1 : xclki2) > 0))
  {
    _nc_xcor_kernel_integrate_block_spectral (xc, xclki1, isauto ? xclki1 : xclki2,
                                              lmin, lmax, isauto, vp);

    if (vp_err != NULL)
      ncm_vector_set_all (vp_err, GSL_NAN);

    return;
  }

  if (ncm_vector_len (vp) != nell)
    g_error ("_nc_xcor_kernel_integrate_block_exact: vector size does not match multipole limits");

  if ((vp_err != NULL) && (ncm_vector_len (vp_err) != nell))
    g_error ("_nc_xcor_kernel_integrate_block_exact: error vector size does not match multipole limits");

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

  sum = g_new0 (gdouble, nell);
  W1  = g_new0 (gdouble, nc_xcor_kernel_integrand_get_len (xclki1));
  W2  = isauto ? W1 : g_new0 (gdouble, nc_xcor_kernel_integrand_get_len (xclki2));

  if (vp_err != NULL)
  {
    err_acc.res   = g_new0 (gdouble, nell);
    err_acc.unk1  = g_new0 (gdouble, nell);
    err_acc.unk2  = isauto ? err_acc.unk1 : g_new0 (gdouble, nell);
    err_acc.prod1 = g_new0 (gdouble, nell);
    err_acc.prod2 = isauto ? err_acc.prod1 : g_new0 (gdouble, nell);
    err_acc.peak1 = g_new0 (gdouble, nell);
    err_acc.peak2 = isauto ? err_acc.peak1 : g_new0 (gdouble, nell);

    err_acc.residuals1 = nc_xcor_kernel_integrand_peek_residuals (xclki1);
    err_acc.residuals2 = isauto ? err_acc.residuals1 : nc_xcor_kernel_integrand_peek_residuals (xclki2);
    err_acc.rows1      = (err_acc.residuals1 != NULL) ? _nc_xcor_gl5_panel_rows (knots1, edges) : NULL;
    err_acc.rows2      = isauto ? err_acc.rows1 :
                         ((err_acc.residuals2 != NULL) ? _nc_xcor_gl5_panel_rows (knots2, edges) : NULL);

    err_acc.dW1 = g_new0 (gdouble, nell);
    err_acc.dW2 = isauto ? err_acc.dW1 : g_new0 (gdouble, nell);
    err_acc.m1  = g_new0 (gdouble, nell);
    err_acc.m2  = isauto ? err_acc.m1 : g_new0 (gdouble, nell);

    err = &err_acc;
  }

  /* The auto/cross distinction is fixed for the whole sweep, so it is resolved
   * once here rather than tested at every quadrature node. */
  if (isauto)
    _nc_xcor_gl5_sweep_auto (xclki1, edges, nell, W1, sum, err);
  else
    _nc_xcor_gl5_sweep_cross (xclki1, xclki2, edges, nell, W1, W2, sum, err);

  for (il = 0; il < nell; il++)
    ncm_vector_set (vp, il, const_factor * sum[il]);

  /* The quadrature is exact, so the only error is the closures' own, propagated
   * through d(W1 W2) = |W1| dW2 + |W2| dW1. Where a closure recorded what its
   * fit achieved, the sweep has already integrated that; what is left is to
   * close the panels that carry no record with the tolerance the fit was asked
   * for. That fallback keeps the two halves of the criterion apart the way the
   * criterion does -- the relative one riding on the product, the peak-scaled
   * floor against the other closure's amplitude -- so with
   * #NcXcorKernel:track-fit-residual off it is the whole estimate, and is then
   * exactly the tolerance-only bound. */
  if (vp_err != NULL)
  {
    for (il = 0; il < nell; il++)
    {
      const gdouble unk_term = reltol1 * err->prod1[il] + sabs1 * err->peak1[il] * err->unk1[il] +
                               reltol2 * err->prod2[il] + sabs2 * err->peak2[il] * err->unk2[il];

      ncm_vector_set (vp_err, il, const_factor * (err->res[il] + unk_term));
    }

    g_free (err->res);
    g_free (err->unk1);
    g_free (err->prod1);
    g_free (err->peak1);
    g_free (err->dW1);
    g_free (err->m1);
    g_clear_pointer (&err->rows1, g_array_unref);

    if (!isauto)
    {
      g_free (err->unk2);
      g_free (err->prod2);
      g_free (err->peak2);
      g_free (err->dW2);
      g_free (err->m2);
      g_clear_pointer (&err->rows2, g_array_unref);
    }
  }

  g_free (sum);
  g_free (W1);

  if (!isauto)
    g_free (W2);

  g_array_unref (edges);
}

void
_nc_xcor_kernel_exact (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
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

  /* Batched by NcXcor:ell-batch-size exactly like _nc_xcor_kernel_cubature():
   * one k-space closure per kernel per batch. The batching is not merely an
   * optimization here -- a single closure spanning more than
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK multipoles is a hard error in
   * nc_xcor_kernel_get_eval_vectorized_full(), so an unbatched sweep aborted on
   * any range wider than that. */
  for (i = 0; i < size; i += block)
  {
    const guint block_lmin = lmin + i;
    const guint block_lmax = MIN (block_lmin + block - 1, lmax);
    NcmVector *vp_i        = ncm_vector_get_subvector (vp, i, block_lmax - block_lmin + 1);
    NcmVector *vp_err_i    = (vp_err != NULL) ? ncm_vector_get_subvector (vp_err, i, block_lmax - block_lmin + 1) : NULL;
    NcXcorKernelIntegrand *xclki1;
    NcXcorKernelIntegrand *xclki2;

    xclki1 = nc_xcor_kernel_get_eval_vectorized_full (xclk1, cosmo, block_lmin, block_lmax, sbi1, xc->closure_type);
    xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized_full (xclk2, cosmo, block_lmin, block_lmax, sbi2, xc->closure_type);

    _nc_xcor_kernel_integrate_block_exact (xc, xclki1, isauto ? xclki1 : xclki2, block_lmin, block_lmax, isauto, vp_i, vp_err_i);

    nc_xcor_kernel_integrand_unref (xclki1);

    if (xclki2 != NULL)
      nc_xcor_kernel_integrand_unref (xclki2);

    ncm_vector_free (vp_i);
    ncm_vector_clear (&vp_err_i);
  }
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
     * it (see plan doc dev-notes/xcor_ultralevin_batching_plan.md §3). */
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
     * evaluations. No safety margin is applied over NcXcor:reltol; the margin
     * that used to be (reltol * 1e-2) only moved the stall.
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

void
_nc_xcor_kernel_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint size  = lmax - lmin + 1;
  const guint block = xc->ell_batch_size;
  guint i;

  _nc_xcor_check_kernel_tolerance (xc, xclk1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, xclk2);

  /* Batched by NcXcor:ell-batch-size exactly like _nc_xcor_kernel_exact():
   * one k-space closure per kernel per batch, built here and handed to the
   * same per-block integrator #NcXcorSolver drives with closures of its own. */
  for (i = 0; i < size; i += block)
  {
    const guint nells      = MIN (block, size - i);
    const guint block_lmin = lmin + i;
    const guint block_lmax = block_lmin + nells - 1;
    NcmVector *vp_i        = ncm_vector_get_subvector (vp, i, nells);
    NcXcorKernelIntegrand *xclki1;
    NcXcorKernelIntegrand *xclki2;

    xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, block_lmin, block_lmax, xc->closure_type);
    xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized (xclk2, cosmo, block_lmin, block_lmax, xc->closure_type);

    _nc_xcor_kernel_integrate_block_cubature (xc, xclki1, xclki2, block_lmin, block_lmax, isauto, vp_i);

    nc_xcor_kernel_integrand_unref (xclki1);

    if (xclki2 != NULL)
      nc_xcor_kernel_integrand_unref (xclki2);

    ncm_vector_free (vp_i);
  }
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
 *
 * Computes the outer k-integral for one ell-block from integrand(s) the
 * caller already built, instead of building them internally the way
 * _nc_xcor_kernel_cubature() does. This is #NcXcorSolver's building block
 * for sharing one kernel's k-space closure across every request that needs
 * it in a given block, instead of rebuilding it once per pair (see plan
 * doc dev-notes/xcor_ultralevin_batching_plan.md §5-§6). Not public: see
 * nc_xcor_priv.h.
 *
 * The caller retains ownership of @xclki1/@xclki2 -- they are not unreffed
 * here.
 */
void
_nc_xcor_kernel_integrate_block_cubature (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint size           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcmVector *err             = ncm_vector_new (size);
  NcmIntegralND *xcor_int_nd;
  NcXcorArg *xcor_arg;

  if (ncm_vector_len (vp) != size)
    g_error ("_nc_xcor_kernel_integrate_block_cubature: vector size does not match multipole limits");

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
  xcor_arg->W1     = g_new (gdouble, size);

  if (!isauto)
  {
    xcor_arg->xclki2 = xclki2;
    xcor_arg->W2     = g_new (gdouble, size);
  }

  ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
  ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
  ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

  _nc_xcor_kernel_cubature_runs (xcor_int_nd, xcor_arg, isauto, size, vp, err);

  ncm_vector_scale (vp, const_factor);

  g_free (xcor_arg->W1);

  if (!isauto)
    g_free (xcor_arg->W2);

  g_object_unref (xcor_int_nd);
  ncm_vector_free (err);
}

