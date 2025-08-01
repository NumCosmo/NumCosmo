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

/* function to perform log(exp(lna) + exp(lnb)) maintaining numerical precision */
static double logaddexp(const double x, const double y);

/* function to perform log(exp(lna) + exp(lnb)) maintaining numerical precision */
double logaddexp(const double x, const double y){
  double tmp = x - y;
  if ( x == y || fabs(tmp) < 1e3*GSL_DBL_EPSILON ){ return x + M_LN2; } /* require the x == y to deal with cases when x and y are both -inf */
  else{
    if ( tmp > 0. ){
      return x + gsl_sf_log_1plusx(exp(-tmp));
    }
    else if ( tmp <= 0. ){
      return y + gsl_sf_log_1plusx(exp(tmp));
    }
    else{
      return tmp;
    }
  }
}
