/***************************************************************************
 *            ncm_spline_cubic.c
 *
 *  Wed Nov 21 19:09:20 2007
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
 * SECTION:ncm_spline_cubic
 * @title: NcmSplineCubic
 * @short_description: Abstract class for implementing cubic splines.
 *
 * This class implements the functions which use a polynomial interpolation
 * method of third degree.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_ABSTRACT_TYPE (NcmSplineCubic, ncm_spline_cubic, NCM_TYPE_SPLINE);

static void
ncm_spline_cubic_init (NcmSplineCubic *sc)
{
	sc->init    = FALSE;

	sc->b       = NULL;
	sc->c       = NULL;
	sc->d       = NULL;

	sc->g       = NULL;
	sc->diag    = NULL;
	sc->offdiag = NULL;
}

static void _ncm_spline_cubic_free (NcmSplineCubic *sc);

static void
ncm_spline_cubic_finalize (GObject *object)
{
	NcmSplineCubic *sc = NCM_SPLINE_CUBIC (object);
	_ncm_spline_cubic_free (sc);

  /* Chain up : end */
	G_OBJECT_CLASS (ncm_spline_cubic_parent_class)->finalize (object);
}

static void _ncm_spline_cubic_reset (NcmSpline *s);
static gdouble _ncm_spline_cubic_eval (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_cubic_deriv (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_cubic_deriv2 (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_cubic_deriv_nmax (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_cubic_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);

static void
ncm_spline_cubic_class_init (NcmSplineCubicClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcmSplineClass* s_class = NCM_SPLINE_CLASS (klass);

	object_class->finalize = ncm_spline_cubic_finalize;

  s_class->reset        = &_ncm_spline_cubic_reset;
	s_class->eval         = &_ncm_spline_cubic_eval;
	s_class->deriv        = &_ncm_spline_cubic_deriv;
	s_class->deriv2       = &_ncm_spline_cubic_deriv2;
  s_class->deriv_nmax   = &_ncm_spline_cubic_deriv_nmax;
	s_class->integ        = &_ncm_spline_cubic_integ;
}

static void
_ncm_spline_cubic_alloc (NcmSplineCubic *sc, gsize n)
{
	g_assert (!sc->init);

	sc->b = ncm_vector_new (n);
	sc->c = ncm_vector_new (n);
	sc->d = ncm_vector_new (n);

	sc->g = ncm_vector_new (n);
	sc->diag = ncm_vector_new (n);
	sc->offdiag = ncm_vector_new (n);

	sc->init = TRUE;
	sc->len = n;

	return;
}

static void
_ncm_spline_cubic_free (NcmSplineCubic *sc)
{
	if (sc->init)
	{
		ncm_vector_free (sc->b);
		ncm_vector_free (sc->c);
		ncm_vector_free (sc->d);
		sc->b = NULL;
		sc->c = NULL;
		sc->d = NULL;

		ncm_vector_free (sc->g);
		ncm_vector_free (sc->diag);
		ncm_vector_free (sc->offdiag);
		sc->g = NULL;
		sc->diag = NULL;
		sc->offdiag = NULL;

		sc->init = FALSE;
	}
}

static void
_ncm_spline_cubic_reset (NcmSpline *s)
{
	NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);

	if (sc->init)
	{
		if (sc->len != s->len)
		{
			_ncm_spline_cubic_free (sc);
			_ncm_spline_cubic_alloc (sc, s->len);
		}
	}
	else
		_ncm_spline_cubic_alloc (sc, s->len);
}

static gdouble
_ncm_spline_cubic_eval (const NcmSpline *s, const gdouble x)
{
	const NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	const size_t i = ncm_spline_get_index (s, x);
	{
		const gdouble delx = x - ncm_vector_get (s->xv, i);
    const gdouble a_i  = ncm_vector_get (s->yv, i);
		const gdouble b_i  = ncm_vector_fast_get (sc->b, i);
		const gdouble c_i  = ncm_vector_fast_get (sc->c, i);
		const gdouble d_i  = ncm_vector_fast_get (sc->d, i);		
#ifdef HAVE_FMA
    return fma (fma (fma (d_i, delx, c_i), delx, b_i), delx, a_i);
#else
    return a_i + delx * (b_i + delx * (c_i + delx * d_i));
#endif /* HAVE_FMA */
	}
}

static gdouble
_ncm_spline_cubic_deriv (const NcmSpline *s, const gdouble x)
{
	const NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	const size_t i = ncm_spline_get_index (s, x);

	{
		const gdouble delx = x - ncm_vector_get (s->xv, i);
		const gdouble b_i  = ncm_vector_fast_get (sc->b, i);
		const gdouble c2_i = 2.0 * ncm_vector_fast_get (sc->c, i);
		const gdouble d3_i = 3.0 * ncm_vector_fast_get (sc->d, i);

#ifdef HAVE_FMA
    return fma (fma (delx, d3_i, c2_i), delx, b_i);
#else
		return b_i + delx * (c2_i + delx * d3_i);
#endif /* HAVE_FMA */
	}
}

static gdouble
_ncm_spline_cubic_deriv2 (const NcmSpline *s, const gdouble x)
{
	const NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	const size_t i = ncm_spline_get_index (s, x);

	{
		const gdouble delx = x - ncm_vector_get (s->xv, i);
		const gdouble c2_i = 2.0 * ncm_vector_fast_get (sc->c, i);
		const gdouble d6_i = 6.0 * ncm_vector_fast_get (sc->d, i);

#ifdef HAVE_FMA
    return fma (delx, d6_i, c2_i);
#else
		return c2_i + delx * d6_i;
#endif /* HAVE_FMA */
	}
}

static gdouble
_ncm_spline_cubic_deriv_nmax (const NcmSpline *s, const gdouble x)
{
	const NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	const size_t i = ncm_spline_get_index (s, x);

	{
		const gdouble d_i = ncm_vector_fast_get (sc->d, i);
		return 6.0 * d_i;
	}
}

static gdouble
_ncm_spline_cubic_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
	const NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	size_t i;
	const size_t index_a = ncm_spline_get_index (s, x0);
	const size_t index_b = ncm_spline_get_index (s, x1);
	gdouble result = 0.0;

	for (i = index_a; i <= index_b; i++)
	{
		const gdouble x_hi = ncm_vector_get (s->xv, i + 1);
		const gdouble x_lo = ncm_vector_get (s->xv, i);
		const gdouble y_lo = ncm_vector_get (s->yv, i);

		{
			const gdouble b_i = ncm_vector_fast_get (sc->b, i);
			const gdouble c_i = ncm_vector_fast_get (sc->c, i);
			const gdouble d_i = ncm_vector_fast_get (sc->d, i);

			if (i == index_a || i == index_b)
			{
				gdouble a = (i == index_a) ? x0 : x_lo;
				gdouble b = (i == index_b) ? x1 : x_hi;
				result += _ncm_spline_util_integ_eval (y_lo, b_i, c_i, d_i, x_lo, a, b);
			}
			else
			{
				const gdouble dx = x_hi - x_lo;
				result += dx * (y_lo + dx * (0.5 * b_i + dx * (c_i / 3.0 + 0.25 * d_i * dx)));
			}
		}
	}

	return result;
}
