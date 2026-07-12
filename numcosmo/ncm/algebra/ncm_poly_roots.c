/***************************************************************************
 *            ncm_poly_roots.c
 *
 *  Mon Jul 6 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_poly_roots.c
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
 * NcmPolyRoots:
 *
 * Real roots of low-degree (2-4) polynomials, via bracket-and-bisect on the
 * derivative chain rather than a general eigenvalue-based solver.
 *
 * A real root of a degree-$n$ polynomial can only hide inside one of the
 * monotonic intervals bounded by its own derivative's real roots (Rolle's
 * theorem), so recursing one derivative at a time down to a quadratic
 * (closed-form, numerically safe -- no cancellation-prone quartic/cubic
 * formula is ever used) turns "find every real root" into a handful of
 * polynomial *evaluations* (Horner, $O(\text{degree})$ flops each) plus a
 * safeguarded Newton/bisection per bracket -- never a matrix decomposition.
 *
 * This was built to replace `gsl_poly_complex_solve()` in a hot loop where
 * its general eigenvalue-based approach (full complex root set, balancing,
 * Hessenberg reduction, QR iteration) measured slower overall than the
 * finite-difference method it was meant to speed up; this bracket-based
 * approach, needing only real roots, measured 1.6-1.8x faster than that
 * same baseline. See `nc_galaxy_shape_intrinsic_mode.c` for the physics
 * context this was extracted from.
 *
 * Each degree has its own hardcoded (not looped-over-a-runtime-degree)
 * evaluator: these are meant for hot inner loops (called up to ~100 times
 * per bracket), so every polynomial evaluation is unrolled straight-line
 * arithmetic rather than depending on the degree being known at compile
 * time for the optimizer to specialize.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_poly_roots.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

static void
_ncm_poly_eval_deriv3 (const gdouble *a, gdouble t, gdouble *p, gdouble *dp)
{
  *p  = ((a[3] * t + a[2]) * t + a[1]) * t + a[0];
  *dp = (3.0 * a[3] * t + 2.0 * a[2]) * t + a[1];
}

static void
_ncm_poly_eval_deriv4 (const gdouble *a, gdouble t, gdouble *p, gdouble *dp)
{
  *p  = (((a[4] * t + a[3]) * t + a[2]) * t + a[1]) * t + a[0];
  *dp = ((4.0 * a[4] * t + 3.0 * a[3]) * t + 2.0 * a[2]) * t + a[1];
}

static gdouble
_ncm_poly_cauchy_bound3 (const gdouble *a)
{
  return 1.0 + fmax (fabs (a[0] / a[3]), fmax (fabs (a[1] / a[3]), fabs (a[2] / a[3])));
}

static gdouble
_ncm_poly_cauchy_bound4 (const gdouble *a)
{
  return 1.0 + fmax (fmax (fabs (a[0] / a[4]), fabs (a[1] / a[4])), fmax (fabs (a[2] / a[4]), fabs (a[3] / a[4])));
}

static gdouble
_ncm_poly_bracketed_root3 (const gdouble *a, gdouble lo, gdouble hi)
{
  gdouble flo, dtmp, t;
  gint it;

  _ncm_poly_eval_deriv3 (a, lo, &flo, &dtmp);
  t = 0.5 * (lo + hi);

  for (it = 0; it < 100; it++)
  {
    gdouble f, df, tn;

    _ncm_poly_eval_deriv3 (a, t, &f, &df);

    if (f == 0.0)
      return t;

    if (flo * f < 0.0)
    {
      hi = t;
    }
    else
    {
      lo  = t;
      flo = f;
    }

    tn = (df != 0.0) ? t - f / df : G_MAXDOUBLE;
    t  = ((tn > lo) && (tn < hi)) ? tn : 0.5 * (lo + hi);

    if (fabs (hi - lo) < 1.0e-13 * (1.0 + fabs (t)))
      break;
  }

  return t;
}

static gdouble
_ncm_poly_bracketed_root4 (const gdouble *a, gdouble lo, gdouble hi)
{
  gdouble flo, dtmp, t;
  gint it;

  _ncm_poly_eval_deriv4 (a, lo, &flo, &dtmp);
  t = 0.5 * (lo + hi);

  for (it = 0; it < 100; it++)
  {
    gdouble f, df, tn;

    _ncm_poly_eval_deriv4 (a, t, &f, &df);

    if (f == 0.0)
      return t;

    if (flo * f < 0.0)
    {
      hi = t;
    }
    else
    {
      lo  = t;
      flo = f;
    }

    tn = (df != 0.0) ? t - f / df : G_MAXDOUBLE;
    t  = ((tn > lo) && (tn < hi)) ? tn : 0.5 * (lo + hi);

    if (fabs (hi - lo) < 1.0e-13 * (1.0 + fabs (t)))
      break;
  }

  return t;
}

static gint
_ncm_poly_bracket_and_collect3 (const gdouble *a, const gdouble *crit, gint n_crit, gdouble *roots)
{
  gdouble bp[4];
  gint n_bp = 0, n = 0, i;
  const gdouble bound = _ncm_poly_cauchy_bound3 (a);

  bp[n_bp++] = -bound * 1.01;

  for (i = 0; i < n_crit; i++)
    bp[n_bp++] = crit[i];

  bp[n_bp++] = bound * 1.01;

  for (i = 0; i < n_bp - 1; i++)
  {
    gdouble flo, fhi, dtmp;

    _ncm_poly_eval_deriv3 (a, bp[i], &flo, &dtmp);
    _ncm_poly_eval_deriv3 (a, bp[i + 1], &fhi, &dtmp);

    if (flo == 0.0)
      roots[n++] = bp[i];
    else if (flo * fhi < 0.0)
      roots[n++] = _ncm_poly_bracketed_root3 (a, bp[i], bp[i + 1]);
  }

  return n;
}

static gint
_ncm_poly_bracket_and_collect4 (const gdouble *a, const gdouble *crit, gint n_crit, gdouble *roots)
{
  gdouble bp[5];
  gint n_bp = 0, n = 0, i;
  const gdouble bound = _ncm_poly_cauchy_bound4 (a);

  bp[n_bp++] = -bound * 1.01;

  for (i = 0; i < n_crit; i++)
    bp[n_bp++] = crit[i];

  bp[n_bp++] = bound * 1.01;

  for (i = 0; i < n_bp - 1; i++)
  {
    gdouble flo, fhi, dtmp;

    _ncm_poly_eval_deriv4 (a, bp[i], &flo, &dtmp);
    _ncm_poly_eval_deriv4 (a, bp[i + 1], &fhi, &dtmp);

    if (flo == 0.0)
      roots[n++] = bp[i];
    else if (flo * fhi < 0.0)
      roots[n++] = _ncm_poly_bracketed_root4 (a, bp[i], bp[i + 1]);
  }

  return n;
}

/**
 * ncm_poly_roots_real_quadratic:
 * @a: (array fixed-size=3): ascending coefficients $a_0+a_1x+a_2x^2$, $a_2\neq0$
 * @roots: (out caller-allocates) (array fixed-size=2): the real roots found
 *
 * Real roots of a genuine quadratic ($a_2\neq0$ required -- callers with a
 * possibly-degenerate leading coefficient should check for that themselves,
 * or use ncm_poly_roots_real_quartic_or_lower() further down the chain).
 * Uses the numerically stable form of the quadratic formula (computing the
 * larger root directly, the other via the product-of-roots identity) to
 * avoid cancellation.
 *
 * Returns: the number of real roots found (0, 1, or 2)
 */
gint
ncm_poly_roots_real_quadratic (const gdouble a[3], gdouble roots[2])
{
  const gdouble disc = a[1] * a[1] - 4.0 * a[2] * a[0];
  gint n             = 0;

  if (disc < 0.0)
    return 0;

  {
    const gdouble sq = sqrt (disc);
    const gdouble q  = (a[1] >= 0.0) ? -0.5 * (a[1] + sq) : -0.5 * (a[1] - sq);

    if (q != 0.0)
    {
      roots[n++] = q / a[2];
      roots[n++] = a[0] / q;
    }
    else
    {
      roots[n++] = 0.0;
    }
  }

  return n;
}

/**
 * ncm_poly_roots_real_cubic:
 * @a: (array fixed-size=4): ascending coefficients $a_0+a_1x+\dots+a_3x^3$, $a_3\neq0$
 * @roots: (out caller-allocates) (array fixed-size=3): the real roots found
 *
 * Real roots of a genuine cubic ($a_3\neq0$ required), via bracket-and-
 * bisect on its derivative's (quadratic) real roots -- see #NcmPolyRoots.
 *
 * Returns: the number of real roots found (1, 2, or 3)
 */
gint
ncm_poly_roots_real_cubic (const gdouble a[4], gdouble roots[3])
{
  const gdouble deriv[3] = {a[1], 2.0 * a[2], 3.0 * a[3]};
  gdouble crit[2];
  gint n_crit = ncm_poly_roots_real_quadratic (deriv, crit);

  if ((n_crit == 2) && (crit[0] > crit[1]))
  {
    const gdouble tmp = crit[0];

    crit[0] = crit[1];
    crit[1] = tmp;
  }

  return _ncm_poly_bracket_and_collect3 (a, crit, n_crit, roots);
}

/**
 * ncm_poly_roots_real_quartic:
 * @a: (array fixed-size=5): ascending coefficients $a_0+a_1x+\dots+a_4x^4$, $a_4\neq0$
 * @roots: (out caller-allocates) (array fixed-size=4): the real roots found
 *
 * Real roots of a genuine quartic ($a_4\neq0$ required), via bracket-and-
 * bisect on its derivative's (cubic) real roots -- see #NcmPolyRoots.
 *
 * Returns: the number of real roots found (0, 1, 2, 3, or 4)
 */
gint
ncm_poly_roots_real_quartic (const gdouble a[5], gdouble roots[4])
{
  const gdouble deriv[4] = {a[1], 2.0 * a[2], 3.0 * a[3], 4.0 * a[4]};
  gdouble crit[3];
  gint n_crit = ncm_poly_roots_real_cubic (deriv, crit);
  gint i, j;

  /* insertion sort, n_crit <= 3 */
  for (i = 1; i < n_crit; i++)
  {
    const gdouble key = crit[i];

    for (j = i - 1; (j >= 0) && (crit[j] > key); j--)
      crit[j + 1] = crit[j];

    crit[j + 1] = key;
  }

  return _ncm_poly_bracket_and_collect4 (a, crit, n_crit, roots);
}

/**
 * ncm_poly_roots_real_quartic_or_lower:
 * @a: (array fixed-size=5): ascending coefficients $a_0+a_1x+\dots+a_4x^4$
 * @roots: (out caller-allocates) (array fixed-size=4): the real roots found
 *
 * Same as ncm_poly_roots_real_quartic(), but safe to call with a leading
 * coefficient that turns out to be (numerically) zero: trims from $a_4$
 * downward, relative to the scale of all five coefficients together, and
 * dispatches to whichever degree is actually genuine. A polynomial with all
 * five coefficients indistinguishable from zero returns 0 roots rather than
 * dividing by zero.
 *
 * Returns: the number of real roots found (0 to 4)
 */
gint
ncm_poly_roots_real_quartic_or_lower (const gdouble a[5], gdouble roots[4])
{
  const gdouble scale = fabs (a[0]) + fabs (a[1]) + fabs (a[2]) + fabs (a[3]) + fabs (a[4]) + 1.0e-300;
  gint deg;

  for (deg = 4; deg >= 1; deg--)
  {
    if (fabs (a[deg]) > 1.0e-14 * scale)
      break;
  }

  switch (deg)
  {
    case 4:
      return ncm_poly_roots_real_quartic (a, roots);

    case 3:
      return ncm_poly_roots_real_cubic (a, roots);

    case 2:
      return ncm_poly_roots_real_quadratic (a, roots);

    case 1:
      roots[0] = -a[0] / a[1];

      return 1;

    default:
      return 0; /* LCOV_EXCL_LINE : all five coefficients numerically zero */
  }
}

