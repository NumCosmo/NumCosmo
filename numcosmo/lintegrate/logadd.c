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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>

#include <math.h>
#include <float.h>

/* Numerically stable log(exp(x) + exp(y)) */
static inline double
logaddexp (double x, double y)
{
  double tmp = x - y;

  if (fabs (tmp) < 1e-12)
    return x + M_LN2;
  else if (tmp > 0)
    return x + log1p (exp (-tmp));
  else
    return y + log1p (exp (tmp));
}

