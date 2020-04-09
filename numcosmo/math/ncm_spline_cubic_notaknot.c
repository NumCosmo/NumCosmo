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
 * @title: NcmSplineCubicNotaknot
 * @short_description: Cubic spline implementation with 'not a knot' boundary conditions.
 * @stability: Stable
 * @include: numcosmo/math/ncm_spline_cubic_notaknot.h
 *
 * This object implements the necessary functions to compute a cubic spline with
 * boundary conditions obtained with the 'not a knot' method.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcmSplineCubicNotaknot, ncm_spline_cubic_notaknot, NCM_TYPE_SPLINE_CUBIC);

static void
ncm_spline_cubic_notaknot_init (NcmSplineCubicNotaknot *s)
{
  NCM_UNUSED (s);
}

static void
_ncm_spline_cubic_notaknot_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_cubic_notaknot_parent_class)->finalize (object);
}

static const gchar *_ncm_spline_cubic_notaknot_name (NcmSpline *s);
static void _ncm_spline_cubic_notaknot_prepare_base (NcmSpline *s);
static void _ncm_spline_cubic_notaknot_prepare (NcmSpline *s);
static gsize _ncm_spline_cubic_notaknot_min_size (const NcmSpline *s);
static NcmSpline *_ncm_spline_cubic_notaknot_copy_empty (const NcmSpline *s);

static void
ncm_spline_cubic_notaknot_class_init (NcmSplineCubicNotaknotClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmSplineClass *s_class    = NCM_SPLINE_CLASS (klass);
  
  object_class->finalize = &_ncm_spline_cubic_notaknot_finalize;
  
  s_class->name         = &_ncm_spline_cubic_notaknot_name;
  s_class->prepare      = &_ncm_spline_cubic_notaknot_prepare;
  s_class->prepare_base = &_ncm_spline_cubic_notaknot_prepare_base;
  s_class->min_size     = &_ncm_spline_cubic_notaknot_min_size;
  s_class->copy_empty   = &_ncm_spline_cubic_notaknot_copy_empty;
}

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
  NCM_UNUSED (s);
  
  return ncm_spline_cubic_notaknot_new ();
}

static gsize
_ncm_spline_cubic_notaknot_min_size (const NcmSpline *s)
{
  NCM_UNUSED (s);
  
  return 6;
}

static const gchar *
_ncm_spline_cubic_notaknot_name (NcmSpline *s)
{
  NCM_UNUSED (s);
  
  return "NcmSplineCubicNotaknot";
}

#define NUMCOSMO_CHECK_SPLINE_NODES 1

static void
_ncm_spline_cubic_notaknot_prepare_base (NcmSpline *s)
{
  NcmSplineCubic *sc           = NCM_SPLINE_CUBIC (s);
  const size_t size            = s->len;
  const size_t n               = size - 1;
  const size_t sys_size        = size - 2;
  const size_t nm1             = n - 1;
  const size_t nm2             = nm1 - 1;
  const size_t nm3             = nm2 - 1;
  const gdouble h_0            = ncm_vector_get (s->xv, 1) - ncm_vector_get (s->xv, 0);
  const gdouble h_1            = ncm_vector_get (s->xv, 2) - ncm_vector_get (s->xv, 1);
  const gdouble h_nm1          = ncm_vector_get (s->xv, n)   - ncm_vector_get (s->xv, nm1);
  const gdouble h_nm2          = ncm_vector_get (s->xv, nm1) - ncm_vector_get (s->xv, nm2);
  const gdouble h_0_p_2h_1     = h_0   + 2.0 * h_1;
  const gdouble h_nm1_p_2h_nm2 = h_nm1 + 2.0 * h_nm2;
  const gdouble ydiff_0        = ncm_vector_get (s->yv, 1) - ncm_vector_get (s->yv, 0);
  const gdouble ydiff_1        = ncm_vector_get (s->yv, 2) - ncm_vector_get (s->yv, 1);
  register gdouble h_i         = h_0;
  register gdouble h_ip1       = h_1;
  register gdouble ydiff_i     = ydiff_0;
  register gdouble ydiff_ip1   = ydiff_1;
  
  NcmVector *g = sc->g;
  
  sc->g = sc->c;
  
  size_t i;
  
  g_assert (sys_size > 1);
#ifdef NUMCOSMO_CHECK_SPLINE_NODES
  
  if ((h_0 < 0.0) || (h_1 < 0.0) || (h_nm1 < 0.0))
    g_error ("_ncm_spline_cubic_notaknot_prepare_base: in node [0, 1, 2, %zu, %zu] (% 20.15g, % 20.15g, % 20.15g, % 20.15g, % 20.15g), (% 20.15g, % 20.15g, % 20.15g, % 20.15g, % 20.15g).", nm1, n,
             ncm_vector_get (s->xv, 0), ncm_vector_get (s->xv, 1), ncm_vector_get (s->xv, 2), ncm_vector_get (s->xv, nm1), ncm_vector_get (s->xv, n),
             ncm_vector_get (s->yv, 0), ncm_vector_get (s->yv, 1), ncm_vector_get (s->yv, 2), ncm_vector_get (s->yv, nm1), ncm_vector_get (s->yv, n));

#endif
  
  for (i = 1; i <= nm3; i++)
  {
    h_i       = h_ip1;
    h_ip1     = ncm_vector_get (s->xv, i + 2) - ncm_vector_get (s->xv, i + 1);
    ydiff_i   = ydiff_ip1;
    ydiff_ip1 = ncm_vector_get (s->yv, i + 2) - ncm_vector_get (s->yv, i + 1);
    
#ifdef NUMCOSMO_CHECK_SPLINE_NODES
    
    if (h_ip1 <= 0.0)
      g_error ("_ncm_spline_cubic_notaknot_prepare_base: in node %zd of %zd x_ip2 = % 20.15g and x_ip1 = % 20.15g, y_ip2 = % 20.15g and y_ip1 = % 20.15g.",
               i, n,
               ncm_vector_get (s->xv, i + 2), ncm_vector_get (s->xv, i + 1),
               ncm_vector_get (s->yv, i + 2), ncm_vector_get (s->yv, i + 1));

#endif
    
    {
      const gdouble g_i   = 1.0 / h_i;
      const gdouble g_ip1 = 1.0 / h_ip1;
      
      ncm_vector_fast_set (sc->offdiag, i, h_ip1);
      ncm_vector_fast_set (sc->diag,    i, 2.0 * (h_ip1 + h_i));
      ncm_vector_fast_set (sc->g,   1 + i, 3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i));
    }
  }
  
  
  {
    const gdouble g_0 = 1.0 / h_0;
    const gdouble g_1 = 1.0 / h_1;
    
    ncm_vector_fast_set (sc->c,         0, 3.0 * (ydiff_1 * g_1 -  ydiff_0 * g_0) / h_0_p_2h_1);
    ncm_vector_fast_set (sc->c,         1, ncm_vector_fast_get (sc->c, 0) * h_1 / (h_1 + h_0));
    
    ncm_vector_fast_addto (sc->diag,    1, (h_0 - h_1) * h_1 / h_0_p_2h_1);
    ncm_vector_fast_subfrom (sc->g,     2, h_1 * ncm_vector_fast_get (sc->c, 1));
  }
  
  {
    const gdouble ydiff_nm2 = ncm_vector_get (s->yv, nm1) - ncm_vector_get (s->yv, nm2);
    const gdouble ydiff_nm1 = ncm_vector_get (s->yv, n)   - ncm_vector_get (s->yv, nm1);
    const gdouble g_nm2     = 1.0 / h_nm2;
    const gdouble g_nm1     = 1.0 / h_nm1;
    
    ncm_vector_fast_set (sc->c,   n, 3.0 * (ydiff_nm1 * g_nm1 - ydiff_nm2 * g_nm2) / h_nm1_p_2h_nm2);
    ncm_vector_fast_set (sc->c, nm1, ncm_vector_fast_get (sc->c, n) * h_nm2 / (h_nm2 + h_nm1));
    
    ncm_vector_fast_addto (sc->diag,    nm1 - 2, (h_nm1 - h_nm2) * h_nm2 / h_nm1_p_2h_nm2);
    ncm_vector_fast_subfrom (sc->g, 1 + nm1 - 2, h_nm2 * ncm_vector_fast_get (sc->c, nm1));
  }
  
  {
    gsize loc_sys_size = sys_size - 2;
    gint info          = ncm_lapack_dptsv (ncm_vector_ptr (sc->diag, 1),
                                           ncm_vector_ptr (sc->offdiag, 1),
                                           ncm_vector_ptr (sc->g, 2),
                                           ncm_vector_ptr (sc->g, 2),
                                           loc_sys_size);
    
    sc->g = g;
    NCM_LAPACK_CHECK_INFO ("dptsv", info);
  }
  
  {
    const gdouble c_2   = ncm_vector_fast_get (sc->c, 2);
    const gdouble c_nm2 = ncm_vector_fast_get (sc->c, nm2);
    
    ncm_vector_fast_subfrom (sc->c, 0, c_2 * (2.0 * h_0 + h_1) / h_0_p_2h_1);
    ncm_vector_fast_subfrom (sc->c, 1, c_2 * (h_1 - h_0) / h_0_p_2h_1);
    
    ncm_vector_fast_subfrom (sc->c, nm1, c_nm2 * (h_nm2 - h_nm1) / h_nm1_p_2h_nm2);
    ncm_vector_fast_subfrom (sc->c,   n, c_nm2 * (h_nm2 + 2.0 * h_nm1) / h_nm1_p_2h_nm2);
  }
  
  return;
}

static void
_ncm_spline_cubic_notaknot_prepare (NcmSpline *s)
{
  NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
  const size_t size  = s->len;
  const size_t n     = size - 1;
  size_t i;
  
  _ncm_spline_cubic_notaknot_prepare_base (s);
  
  for (i = 0; i < n; i++)
  {
    const gdouble dx    = ncm_vector_get (s->xv, i + 1) - ncm_vector_get (s->xv, i);
    const gdouble dy    = ncm_vector_get (s->yv, i + 1) - ncm_vector_get (s->yv, i);
    const gdouble c_ip1 = ncm_vector_fast_get (sc->c, i + 1);
    const gdouble c_i   = ncm_vector_fast_get (sc->c, i);
    
    ncm_vector_fast_set (sc->b, i, (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0);
    ncm_vector_fast_set (sc->d, i, (c_ip1 - c_i) / (3.0 * dx));
  }
  
  return;
}

