/* err.c
 *
 * Copyright (C) 2017 Matthew Pitkin
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

/* slightly modified GSL rescale error function from http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/err.c */
static double rescale_error (double err, const double result_abs, const double result_asc);

static double rescale_error (double err, const double result_abs, const double result_asc) {
  if (gsl_isinf(result_asc) != -1 && gsl_isinf(err) != -1) {
    double scale = 1.5*(log(200.) + err - result_asc);

    if (scale < 0.) {
      err = result_asc + scale;
    }
    else {
      err = result_asc;
    }
  }

  if (result_abs > log(GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))) {
    double min_err = log(50 * GSL_DBL_EPSILON) + result_abs;

    if (min_err > err) {
      err = min_err;
    }
  }

  return err;
}

