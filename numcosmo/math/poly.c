/***************************************************************************
 *            poly.c
 *
 *  Wed Nov 21 14:39:45 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:poly
 * @title: Polynomials
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/poly.h"

#include <gsl/gsl_sf_log.h>

/**
 * ncm_poly_new: (skip)
 * @degree: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/ 
gsl_vector *
ncm_poly_new (gint degree)
{
  gsl_vector *poly = gsl_vector_calloc (degree + 1);
  return poly;
}

/**
 * ncm_poly_eval: (skip)
 * @poly: FIXME
 * @x: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_poly_eval (gsl_vector *poly, gdouble x)
{
  gint i;
  gdouble res = gsl_vector_get (poly, poly->size - 1);
  for (i = poly->size - 1; i > 0; i--)
    res = gsl_vector_get (poly, i - 1) + x * res / i;
  return res;
}

/**
 * ncm_poly_eval_diff: (skip)
 * @poly: FIXME
 * @n: FIXME 
 * @x: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_poly_eval_diff (gsl_vector *poly, gint n, gdouble x)
{
  gint i;
  gdouble res;
  if (n >= poly->size)
    return 0.0;
  res = gsl_vector_get (poly, poly->size - 1);
  for (i = poly->size - 1; i > n; i--)
    res = gsl_vector_get (poly, i - 1) + x * res / (i - n);
  return res;
}

/**
 * ncm_poly_eval_int: (skip)
 * @poly: FIXME
 * @x: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_poly_eval_int (gsl_vector *poly, gdouble x)
{
  gint i;
  gdouble res;
  res = gsl_vector_get (poly, poly->size - 1);
  for (i = poly->size - 1; i > 0; i--)
    res = gsl_vector_get (poly, i - 1) + x * res / (i + 1);
  return res * x;
}

/**
 * ncm_poly_eval_int_P_over_x: (skip)
 * @poly: FIXME
 * @x: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_poly_eval_int_P_over_x (gsl_vector *poly, gdouble x)
{
  gint i;
  gdouble res, res1;
  if (poly->size == 1)
    return gsl_sf_log (x);
  res = gsl_vector_get (poly, poly->size - 1) / (poly->size - 1.0);
  res1 = res;
  for (i = poly->size - 1; i > 1; i--)
  {
    res = gsl_vector_get (poly, i - 1) / (i - 1.0) + x * res / i;
    res1 = gsl_vector_get (poly, i - 1) / (i - 1.0) + res1 / i; 
  }
  
  return x * res + gsl_vector_get (poly, 0) * gsl_sf_log (x) - res1;
}
