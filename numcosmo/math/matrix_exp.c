/***************************************************************************
 *            matrix_exp.c
 *
 *  Thu Jun 18 12:50:32 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:matrix_exp
 * @title: Matrix Exponential
 * @short_description: Simple functions to calculate matrix exponential (only 2x2 to date)
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <string.h>
#include <glib.h>

/**
 * ncm_matrix_exp_2x2: (skip)
 * @B: FIXME
 * @exp_B: FIXME
 *
 * FIXME 
*/
void 
ncm_matrix_exp_2x2 (gsl_matrix *B, gsl_matrix *exp_B)
{
  gdouble delta2, a_md, exp_a_d_2;

  g_assert (B->size1 == 2 && B->size2 == 2);
#define a gsl_matrix_get (B, 0, 0)
#define b gsl_matrix_get (B, 0, 1)
#define c gsl_matrix_get (B, 1, 0)
#define d gsl_matrix_get (B, 1, 1)

  exp_a_d_2 = exp ((a + d) / 2.0);
  a_md = a - d;
  delta2 = gsl_pow_2 (a_md) + 4.0 * b * c;

  printf ("# why %.15g,%.15g = %.15g\n",a*a, b*c, a*a + b*c);

  
  if (delta2 == 0)
  {
    gsl_matrix_set (exp_B, 0, 0, exp_a_d_2 * (1.0 + a_md / 2.0));
    gsl_matrix_set (exp_B, 0, 1, exp_a_d_2 * (b));
    gsl_matrix_set (exp_B, 1, 0, exp_a_d_2 * (c));
    gsl_matrix_set (exp_B, 1, 1, exp_a_d_2 * (1.0 - a_md / 2.0));
  }
  else if (delta2 < 0)
  {
    gdouble mi_delta = sqrt(fabs(delta2));
    gdouble cos_delta_2 = cos (mi_delta / 2.0);
    gdouble sin_delta_2 = sin (mi_delta / 2.0);

    printf ("<<%.15g|%.15g>> == %g %g %g\n", mi_delta / 2.0, cos (mi_delta / 2.0), sin_delta_2, mi_delta / (2.0 * b), mi_delta/2.0);

    gsl_matrix_set (exp_B, 0, 0, exp_a_d_2 * (mi_delta * cos_delta_2 + a_md * sin_delta_2) / mi_delta);
    gsl_matrix_set (exp_B, 0, 1, exp_a_d_2 * (2.0 * b * sin_delta_2) / mi_delta);
    gsl_matrix_set (exp_B, 1, 0, exp_a_d_2 * (2.0 * c * sin_delta_2) / mi_delta);
    gsl_matrix_set (exp_B, 1, 1, exp_a_d_2 * (mi_delta * cos_delta_2 - a_md * sin_delta_2) / mi_delta);
  }
  else
  {
    gdouble delta = sqrt(delta2);
    gdouble cosh_delta_2 = cosh (delta / 2.0);
    gdouble sinh_delta_2 = sinh (delta / 2.0);

printf ("# [%.15g]EXXXX: %.15g %.15g %.15g | %.15g\n", delta/2.0, exp_a_d_2, 2.0 * b * sinh_delta_2, delta, exp_a_d_2 * (2.0 * b * sinh_delta_2) / delta);
    
    gsl_matrix_set (exp_B, 0, 0, exp_a_d_2 * (delta * cosh_delta_2 + a_md * sinh_delta_2) / delta);
    gsl_matrix_set (exp_B, 0, 1, exp_a_d_2 * (2.0 * b * sinh_delta_2) / delta);
    gsl_matrix_set (exp_B, 1, 0, exp_a_d_2 * (2.0 * c * sinh_delta_2) / delta);
    gsl_matrix_set (exp_B, 1, 1, exp_a_d_2 * (delta * cosh_delta_2 - a_md * sinh_delta_2) / delta);
  }
}

#undef a
#undef b
#undef c
#undef d
