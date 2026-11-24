/***************************************************************************
 *            ncm_bspline_private.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_bspline_private.h
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

#ifndef _NCM_BSPLINE_PRIVATE_H_
#define _NCM_BSPLINE_PRIVATE_H_

#include <glib.h>
#include "ncm/spline/ncm_spline_bspline.h"

G_BEGIN_DECLS

/*
 * Evaluates at x the B-spline of order k with the ncontrol + k knots t and the ncontrol
 * control points c, by de Boor's recursion (PPPACK bsplvb, as in GSL) on stack scratch.
 * It reads t and c and writes nothing outside its frame. Outside [t[k - 1], t[ncontrol]]
 * the clamped edge span extrapolates its polynomial, matching gsl_bspline_calc(). The
 * order k is at most NCM_SPLINE_BSPLINE_MAX_ORDER; it is not checked here.
 *
 * It is the one evaluation of NcmSplineBSpline and NcmBSpline, compiled into each.
 */
static inline gdouble
_ncm_bspline_eval (const gsize k, const gsize ncontrol, const gdouble *t, const gdouble *c, const gdouble x)
{
  gdouble deltal[NCM_SPLINE_BSPLINE_MAX_ORDER];
  gdouble deltar[NCM_SPLINE_BSPLINE_MAX_ORDER];
  gdouble B[NCM_SPLINE_BSPLINE_MAX_ORDER];
  gsize l, i, j;

  /* Largest span index l in [k - 1, ncontrol - 1] with t[l] <= x */
  if (x < t[k - 1])
  {
    l = k - 1;
  }
  else if (x >= t[ncontrol])
  {
    l = ncontrol - 1;
  }
  else
  {
    gsize lo = k - 1;
    gsize hi = ncontrol;

    while (hi - lo > 1)
    {
      const gsize mid = (lo + hi) / 2;

      if (x < t[mid])
        hi = mid;
      else
        lo = mid;
    }

    l = lo;
  }

  B[0] = 1.0;

  for (j = 0; j + 1 < k; j++)
  {
    gdouble saved = 0.0;

    deltar[j] = t[l + j + 1] - x;
    deltal[j] = x - t[l - j];

    for (i = 0; i <= j; i++)
    {
      const gdouble term = B[i] / (deltar[i] + deltal[j - i]);

      B[i]  = saved + deltar[i] * term;
      saved = deltal[j - i] * term;
    }

    B[j + 1] = saved;
  }

  {
    gdouble res = 0.0;

    for (i = 0; i < k; i++)
      res += B[i] * c[l - (k - 1) + i];

    return res;
  }
}

G_END_DECLS

#endif /* _NCM_BSPLINE_PRIVATE_H_ */

