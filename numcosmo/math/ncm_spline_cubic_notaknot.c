/***************************************************************************
 *            ncm_spline_cubic_notaknot.c
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
 * SECTION:ncm_spline_cubic_notaknot
 * @title: Notaknot Cubic Spline
 * @short_description: Notaknot boundary conditions.
 *
 * This object implements the necessary functions to compute a cubic spline with
 * boundary conditions obtained with the notaknot method.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "nc_macros.h"

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

G_DEFINE_TYPE (NcmSplineCubicNotaknot, ncm_spline_cubic_notaknot, NCM_TYPE_SPLINE_CUBIC);

/**
 * ncm_spline_cubic_notaknot_new:
 *
 * This function returns a new cubic #NcmSpline.
 *
 * Returns: a new #NcmSpline.
 */
NcmSpline *
ncm_spline_cubic_notaknot_new ()
{
	return g_object_new (NCM_TYPE_SPLINE_CUBIC_NOTAKNOT, NULL);
}

/**
 * ncm_spline_cubic_notaknot_new_full:
 * @xv: #NcmVector of knots.
 * @yv: #NcmVector of the values of the function, to be interpolated, computed at @xv.
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it.
 *
 * This function returns a new #NcmSpline setting all its members.
 *
 * Returns: a new #NcmSpline.
 */
NcmSpline *
ncm_spline_cubic_notaknot_new_full (NcmVector *xv, NcmVector *yv, gboolean init)
{
	NcmSpline *s = ncm_spline_cubic_notaknot_new ();
	ncm_spline_set (s, xv, yv, init);
	return s;
}

static NcmSpline *
_ncm_spline_cubic_notaknot_copy_empty (const NcmSpline *s)
{
	return ncm_spline_cubic_notaknot_new ();
}

gsize
_ncm_spline_notaknot_min_size (const NcmSpline *s)
{
	return 6;
}

static void
_ncm_spline_notaknot_prepare_base (NcmSpline *s)
{
	NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	const size_t size = s->len;
	const size_t n = size - 1;
	const size_t sys_size = size - 2;
	const size_t nm1 = n - 1;
	const size_t nm2 = nm1 - 1;
	const size_t nm3 = nm2 - 1;
	const gdouble h_0 = ncm_vector_get (s->xv, 1) - ncm_vector_get (s->xv, 0);
	const gdouble h_1 = ncm_vector_get (s->xv, 2) - ncm_vector_get (s->xv, 1);
	const gdouble h_nm1 = ncm_vector_get (s->xv, n)   - ncm_vector_get (s->xv, nm1);
	const gdouble h_nm2 = ncm_vector_get (s->xv, nm1) - ncm_vector_get (s->xv, nm2);
#ifdef HAVE_LAPACKA
	NcmVector *g = sc->g;
	sc->g = sc->c;
#endif

	size_t start_i = 0, pad_i = 0;
	size_t i;

	g_assert (sys_size > 1);

	for (i = 1; i <= nm3; i++)
	{
		const gdouble h_i       = ncm_vector_get (s->xv, i + 1) - ncm_vector_get (s->xv, i);
		const gdouble h_ip1     = ncm_vector_get (s->xv, i + 2) - ncm_vector_get (s->xv, i + 1);
		const gdouble ydiff_i   = ncm_vector_get (s->yv, i + 1) - ncm_vector_get (s->yv, i);
		const gdouble ydiff_ip1 = ncm_vector_get (s->yv, i + 2) - ncm_vector_get (s->yv, i + 1);
		const gdouble g_i = 1.0 / h_i;
		const gdouble g_ip1 = 1.0 / h_ip1;

		ncm_vector_fast_set (sc->offdiag, i, h_ip1);
		ncm_vector_fast_set (sc->diag, i,    2.0 * (h_ip1 + h_i));
		ncm_vector_fast_set (sc->g, 1 + i,   3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i));
	}

	{
		const gdouble ydiff_0 = ncm_vector_get (s->yv, 1) - ncm_vector_get (s->yv, 0);
		const gdouble ydiff_1 = ncm_vector_get (s->yv, 2) - ncm_vector_get (s->yv, 1);
		if (fabs((h_0 - h_1) / h_0) < GSL_SQRT_DBL_EPSILON)
		{
			start_i = 1;
			ncm_vector_fast_set (sc->c, 1, 3.0 * (ydiff_1 -  ydiff_0 * h_1 / h_0) / ((h_1 + h_0) * (2.0 * h_1 + h_0)));
			ncm_vector_fast_subfrom (sc->g, 1 + 1, h_1 * ncm_vector_fast_get (sc->c, 1));
		}
		else
		{
			ncm_vector_fast_set (sc->offdiag, 0, h_1);
			ncm_vector_fast_set (sc->diag, 0,    h_1 * (2.0 * h_1 + h_0) / (h_1 - h_0));
			ncm_vector_fast_set (sc->g, 1 + 0,   3.0 * h_1 * (ydiff_1 -  ydiff_0 * h_1 / h_0) / (h_1 * h_1 - h_0 * h_0));
		}
	}

	{
		const gdouble ydiff_nm1 = ncm_vector_get (s->yv, n)   - ncm_vector_get (s->yv, nm1);
		const gdouble ydiff_nm2 = ncm_vector_get (s->yv, nm1) - ncm_vector_get (s->yv, nm2);

		if (fabs((h_nm2 - h_nm1) / h_nm2) < GSL_SQRT_DBL_EPSILON)
		{
			pad_i = 1;
			ncm_vector_fast_set (sc->c, nm1, 3.0 * (ydiff_nm1 * h_nm2 / h_nm1 - ydiff_nm2) / ((h_nm2 + h_nm1) * (2.0 * h_nm2 + h_nm1)));
			ncm_vector_fast_subfrom (sc->g, 1 + nm1 - 2, h_nm2 * ncm_vector_fast_get (sc->c, nm1));
		}
		else
		{
			ncm_vector_fast_set (sc->offdiag, nm2, h_nm1);
			ncm_vector_fast_set (sc->diag,    nm2, h_nm2 * (2.0 * h_nm2 + h_nm1) / (h_nm2 - h_nm1));
			ncm_vector_fast_set (sc->g,   1 + nm2, 3.0 * h_nm2 * (ydiff_nm1 * h_nm2 / h_nm1 - ydiff_nm2) / (h_nm2 * h_nm2 - h_nm1 * h_nm1));
		}
	}

	{
		gsize loc_sys_size = sys_size - start_i - pad_i;
#ifdef HAVE_LAPACKA
		{
			printf ("# nhocC[%zu] % 20.15g\n", nm2, ncm_vector_fast_get (sc->diag, nm2));
			gint info = ncm_lapack_dptsv (&NCM_VECTOR_DATA(sc->diag)[start_i],
			                              &NCM_VECTOR_DATA(sc->offdiag)[start_i],
			                              &NCM_VECTOR_DATA(sc->g)[start_i + 1],
			                              loc_sys_size);
			sc->g = g;
			if (info != 0)
			{
				for (i = 0; i < s->len; i++)
				{
					printf ("[%zu/%zu] [%zu, %zu] % 20.15g %20.15g [% 20.15g % 20.15g]\n", i, loc_sys_size, start_i, pad_i, ncm_vector_get (s->xv, i), ncm_vector_get (s->yv, i),
					        ncm_vector_fast_get (sc->diag, i),
					        ncm_vector_fast_get (sc->offdiag, i));
				}

			}
			NCM_LAPACK_CHECK_INFO ("dptsv", info);
		}
#else
		{
			gsl_vector_view g_vec = gsl_vector_subvector (ncm_vector_gsl (sc->g), 1 + start_i, loc_sys_size);
			gsl_vector_view diag_vec = gsl_vector_subvector (ncm_vector_gsl (sc->diag), start_i, loc_sys_size);
			gsl_vector_view offdiag_vec = gsl_vector_subvector (ncm_vector_gsl (sc->offdiag), start_i, loc_sys_size - 1);
			gsl_vector_view solution_vec = gsl_vector_subvector (ncm_vector_gsl (sc->c), start_i + 1, loc_sys_size);
			gint status = gsl_linalg_solve_symm_tridiag (&diag_vec.vector,
			                                             &offdiag_vec.vector,
			                                             &g_vec.vector,
			                                             &solution_vec.vector);
			if (status != GSL_SUCCESS)
			{
				gint i;
				for (i = 0; i < s->len; i++)
					printf ("x= % 20.8g y = % 20.8g\n", ncm_vector_get (s->xv, i), ncm_vector_get (s->yv, i));
			}

			NC_TEST_GSL_RESULT ("_ncm_spline_notaknot_prepare[gsl_linalg_solve_symm_tridiag]", status);
		}
#endif
		{
			const gdouble c1   = ncm_vector_fast_get (sc->c, 1);
			const gdouble c2   = ncm_vector_fast_get (sc->c, 2);
			const gdouble cnm1 = ncm_vector_fast_get (sc->c, nm1);
			const gdouble cnm2 = ncm_vector_fast_get (sc->c, nm2);

			ncm_vector_fast_set (sc->c, 0, c1 + h_0 * (c1 - c2) / h_1);
			ncm_vector_fast_set (sc->c, n, cnm1 + h_nm1 * (cnm1 - cnm2) / h_nm2);
		}

		return;
	}
}

static void
_ncm_spline_notaknot_prepare (NcmSpline *s)
{
	NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
	const size_t size = s->len;
	const size_t n = size - 1;
	size_t i;

	_ncm_spline_notaknot_prepare_base (s);

	for (i = 0; i < n; i++)
	{
		const gdouble dx = ncm_vector_get (s->xv, i + 1) - ncm_vector_get (s->xv, i);
		const gdouble dy = ncm_vector_get (s->yv, i + 1) - ncm_vector_get (s->yv, i);
		const gdouble c_ip1 = ncm_vector_fast_get (sc->c, i + 1);
		const gdouble c_i = ncm_vector_fast_get (sc->c, i);

		ncm_vector_fast_set (sc->b, i, (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0);
		ncm_vector_fast_set (sc->d, i, (c_ip1 - c_i) / (3.0 * dx));
	}

	return;
}

static void
ncm_spline_cubic_notaknot_init (NcmSplineCubicNotaknot *object)
{
}

static void
ncm_spline_cubic_notaknot_finalize (GObject *object)
{

  /* Chain up : end */
	G_OBJECT_CLASS (ncm_spline_cubic_notaknot_parent_class)->finalize (object);
}

static void
ncm_spline_cubic_notaknot_class_init (NcmSplineCubicNotaknotClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcmSplineClass *s_class = NCM_SPLINE_CLASS (klass);

	s_class->name         = "NcmSplineCubicNotaknot";
	s_class->prepare      = &_ncm_spline_notaknot_prepare;
	s_class->prepare_base = &_ncm_spline_notaknot_prepare_base;
	s_class->min_size     = &_ncm_spline_notaknot_min_size;
	s_class->copy_empty   = &_ncm_spline_cubic_notaknot_copy_empty;

	object_class->finalize = ncm_spline_cubic_notaknot_finalize;
}
