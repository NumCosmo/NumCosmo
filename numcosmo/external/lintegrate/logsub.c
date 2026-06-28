/* Copyright (C) 2017 Matthew Pitkin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <float.h>

/* Numerically stable log(exp(x) - exp(y)) for any x, y
 * Returns -INFINITY if x == y
 */
static inline double
logsubexp (double x, double y)
{
  double a, b;
  double tmp;

  if (x > y)
  {
    a = x;
    b = y;
  }
  else
  {
    a = y;
    b = x;
  }

  tmp = a - b;

  if (tmp > M_LN2)
    /* Large difference: exp(-tmp) small, use log1p */
    return a + log1p (-exp (-tmp));
  else
    /* Small difference: avoid cancellation */
    return a + log (-expm1 (-tmp));
}

/* Macro for convenience, equivalent to LOGDIFF(x, y) */
#define LOGDIFF(x, y) logsubexp (x, y)

