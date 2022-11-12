/***************************************************************************
 *            ncm_sf_sbessel.c
 *
 *  Wed Mar 10 17:15:25 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_sf_sbessel
 * @title: NcmSFSBessel
 * @short_description: Double precision spherical bessel implementation.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sf_sbessel.h"
#include "math/ncm_mpsf_sbessel.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <mpfr.h>
#endif /* NUMCOSMO_GIR_SCAN */

/**
 * ncm_sf_sbessel:
 * @l: Spherical Bessel order $\ell$
 * @x: Spherical Bessel argument $x$
 *
 * Computes Spherical Bessel function $j_\ell(x)$.
 *
 * Returns: the value of $j_\ell(x)$.
*/
gdouble
ncm_sf_sbessel (gulong l, gdouble x)
{
  MPFR_DECL_INIT (res, 53); /* Should it be 53? FIXME */
  gdouble res_d;
  ncm_mpsf_sbessel_d (l, x, res, GMP_RNDN);
  res_d = mpfr_get_d (res, GMP_RNDN);
  return res_d;
}

static void
_taylor_jl (const glong l, const gdouble x, const gdouble x2, const gdouble x3, const gdouble jl, const gdouble jlp1, gdouble *deriv)
{
  const gdouble llm1 = l * (l - 1.0);
  const gdouble llm1lm2 = llm1 * (l - 2.0);

  deriv[0] = jl;
  deriv[1] = (l * jl - x * jlp1) / x;
  deriv[2] = (((llm1 - x2) * jl + 2.0 * x * jlp1) / x2) / (1.0 * 2.0);
  deriv[3] = (((llm1lm2 - (l - 2.0) * x2) * jl - x * (l * (l + 1.0) + 6.0 - x2) * jlp1) / x3) / (1.0 * 2.0 * 3.0);
}

/**
 * ncm_sf_sbessel_taylor:
 * @l: Spherical Bessel order $\ell$
 * @x: Spherical Bessel argument $x$
 * @djl: (out) (array fixed-size=4): Output power series coefficients
 *
 * Computes Spherical Bessel function power series
 * coefficients up to order three, i.e.,
 * $$\left(j_\ell(x),\; j'_\ell(x), \frac{j''_\ell(x)}{2!}, \frac{j'''_\ell(x)}{3!}\right).$$
*/
void
ncm_sf_sbessel_taylor (gulong l, gdouble x, gdouble *djl)
{
  const gdouble jl = ncm_sf_sbessel (l, x);
  const gdouble jlp1 = ncm_sf_sbessel (l + 1, x);
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;

  _taylor_jl (l, x, x2, x3, jl, jlp1, djl);
  return;
}

static gdouble
_ncm_sf_sbessel_spline_calc (gdouble x, gpointer data)
{
	gulong *l = (gulong *)data;
	return ncm_sf_sbessel (*l, x);
}

/**
 * ncm_sf_sbessel_spline:
 * @l: Spherical Bessel order $\ell$.
 * @xi: Spherical Bessel interval lower-bound $x_i$.
 * @xf: Spherical Bessel interval lower-bound $x_f$.
 * @reltol: Interpolation error tolerance.
 *
 * Computes a spline approximation of the Spherical Bessel
 * $j_\ell$ in the interval $[x_i, x_f]$.
 *
 * Returns: (transfer full): A #NcmSpline with the Spherical Bessel approximation.
 */
NcmSpline *
ncm_sf_sbessel_spline (gulong l, gdouble xi, gdouble xf, gdouble reltol)
{
  NcmSpline *s = ncm_spline_cubic_notaknot_new ();
  gsl_function F;

  F.function = &_ncm_sf_sbessel_spline_calc;
  F.params = &l;

  ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, xi, xf, 0, reltol);
  return s;
}
