/***************************************************************************
 *            nc_window.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
 * 
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
 * SECTION:nc_window
 * @title: Window Function Abstract Class
 * @short_description: Defines the prototype of the #NcWindow object. 
 * 
 * This module comprises the set of functions to compute the window function in both real 
 * and Fourier spaces as well as its derivative with respect to the scale $R$ in Fourier space. 
 * 
 * In order to study the statistical properties of the density fluctuation field at a certain
 * scale $R$, we use the window function. As an example, to compute the variance of the density
 * contrast at scale $R$, we convolve the window function in the Fourier space with the power spectrum.
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_window.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcWindow, nc_window, G_TYPE_OBJECT);

/**
 * nc_window_new_from_name:
 * @window_name: "#NcWindowTophat" or "#NcWindowGaussian".
 * 
 * This function returns a new #NcWindow whose type is defined by @window_name string.
 * 
 * Returns: A new #NcWindow.
 */
NcWindow *
nc_window_new_from_name (gchar *window_name)
{
  GObject *obj = ncm_cfg_create_from_string (window_name);
  GType window_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (window_type, NC_TYPE_WINDOW))
	g_error ("nc_window_new_from_name: NcWindow %s do not descend from %s\n", window_name, g_type_name (NC_TYPE_WINDOW));
  return NC_WINDOW (obj);
}

/**
 * nc_window_volume:
 * @wf: a #NcWindow.
 *  
 * This function returns the volume of the region (with radius 1) defined by the window function.
 * 
 * Top-hat volume: #NC_WINDOW_VOLUME_TOPHAT.
 * 
 * Gaussian volume: #NC_WINDOW_VOLUME_GAUSSIAN.
 *
 * Returns: The volume (with radius 1) defined by @wf.
 */
gdouble
nc_window_volume (NcWindow *wf)
{
  return NC_WINDOW_GET_CLASS (wf)->volume;
}

/**
 * nc_window_eval_fourier:
 * @wf: a #NcWindow.
 * @k: mode. 
 * @R: scale.
 *  
 * This function computes the window function in the Fourier space. 
 *
 * Returns: The value of the window function in the Fourier space at scale @R.
 */
gdouble 
nc_window_eval_fourier (const NcWindow *wf, const gdouble k, const gdouble R) 
{ 
  return NC_WINDOW_GET_CLASS (wf)->eval_fourier (wf, k, R);
}

/**
 * nc_window_deriv_fourier:
 * @wf: a #NcWindow. 
 * @k: mode.  
 * @R: scale. 
 *
 * This function returns the derivative with respect to R of the window function 
 * in the Fourier space.
 * 
 * Returns: The value of the first derivative of the window function in the Fourier space at scale @R.
 */
gdouble 
nc_window_deriv_fourier (const NcWindow *wf, const gdouble k, const gdouble R) 
{ 
  return NC_WINDOW_GET_CLASS (wf)->deriv_fourier (wf, k, R);
}

/**
 * nc_window_eval_realspace:
 * @wf: a #NcWindow.  
 * @r: distance module to the center point of the filtered region.
 * @R: scale.
 *
 * This function computes the window function in real space. 
 *
 * Returns: The value of the window function in the real space at scale @R.
 */
gdouble 
nc_window_eval_realspace (const NcWindow *wf, const gdouble r, const gdouble R) 
{ 
  return NC_WINDOW_GET_CLASS (wf)->eval_real (wf, r, R);
}

/**
 * nc_window_free: 
 * @wf: a #NcWindow.
 *
 * Atomically decrements the reference count of @wf by one. If the reference count drops to 0, 
 * all memory allocated by @wf is released.
 *
 */
void 
nc_window_free (NcWindow *wf)
{
  g_object_unref (wf);
}

/**
 * nc_window_clear: 
 * @wf: a #NcWindow.
 *
 * Atomically decrements the reference count of @wf by one. If the reference count drops to 0, 
 * all memory allocated by @wf is released. Set the pointer to NULL;
 *
 */
void 
nc_window_clear (NcWindow **wf)
{
  g_clear_object (wf);
}

static void
nc_window_init (NcWindow *wf)
{
}

static void
_nc_window_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_window_parent_class)->finalize (object);
}

static void
nc_window_class_init (NcWindowClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = _nc_window_finalize;
}

