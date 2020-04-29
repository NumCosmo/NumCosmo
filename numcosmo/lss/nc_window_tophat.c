/***************************************************************************
 *            nc_window_tophat.c
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
 * SECTION:nc_window_tophat
 * @title: NcWindowTophat
 * @short_description: A top-hat window function.
 * @stability: Stable
 * @include: numcosmo/lss/nc_window_tophat.h
 *
 * This object implements the #NcWindow class for a top-hat window function.
 *
 * The top-hat window function in the real space is defined as, 
 * \begin{equation*}
 *   W_{TH}(r, R) = \frac{3}{4\pi R^3}, \,\,\,\, \mathrm{for} \,\,\, r \leq R \,\, ,
 * \end{equation*}
 * and 0 otherwise. The mass enclosed within the volume selected by this window function is
 * $$M_{TH}(R)= \frac{4\pi}{3}\bar{\rho} R^3 \, ,$$
 * where $\bar{\rho}(z)$ is the mean density of the universe at redshift $z$.
 *
 * When the function nc_window_eval_fourier() is applied,
 * it returns the top-hat window function in the Fourier space, 
 * \begin{equation*}
 *  W_{th}(k, R) = \frac{3}{(kR)^3} \left[ \sin (kR) - (kR)\cos (kR)\right] = \frac{3}{(kR)} j_1(kR),
 * \end{equation*}
 * where $j_\nu(kR)$ is the spherical Bessel function.
 * The function nc_window_deriv_fourier() returns the derivative with respect to R in Fourier space,
 * \begin{equation*}
 * \frac{\mathrm{d} W_{TH}(k, R)}{\mathrm{d} R} = -\frac{9}{k^3 R^4} \left[ \sin (kR) - (kR)\cos (kR) \right] + \frac{3}{k R^2} \sin (kR) \, .
 * \end{equation*}
 * The derivative with respect to R in real space is performed by nc_window_eval_realspace().
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_window_tophat.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcWindowTophat, nc_window_tophat, NC_TYPE_WINDOW);

/**
 * nc_window_tophat_new:
 *
 * This function returns a #NcWindow with a #NcWindowTophat implementation.
 *
 * Returns: A new #NcWindow.
 */
NcWindow *
nc_window_tophat_new ()
{
  return g_object_new (NC_TYPE_WINDOW_TOPHAT, NULL);
}

static gdouble
_nc_window_tophat_eval_fourier (const NcWindow *wp, const gdouble k, const gdouble R)
{
  gdouble kR = k * R;
  
  NCM_UNUSED (wp);
  
  if (kR == 0.0)
    return 1.0;
  
  return 3.0 * gsl_sf_bessel_j1 (kR) / kR;
}

static gdouble
_nc_window_tophat_deriv_fourier (const NcWindow *wp, const gdouble k, const gdouble R)
{
  gdouble dWT = -3.0 * gsl_sf_bessel_j2 (k * R) / R;
  
  NCM_UNUSED (wp);
  
  return dWT;
}

static gdouble
_nc_window_tophat_eval_realspace (const NcWindow *wp, const gdouble r, const gdouble R)
{
  gdouble WT_realspace;
  gdouble R3 = R * R * R;
  
  NCM_UNUSED (wp);
  
  if (r <= R)
    WT_realspace = 3.0 / (4.0 * M_PI * R3);
  else
    WT_realspace = 0.0;
  
  return WT_realspace;
}

static void
nc_window_tophat_init (NcWindowTophat *nc_window_tophat)
{
  NCM_UNUSED (nc_window_tophat);
}

static void
nc_window_tophat_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_window_tophat_parent_class)->finalize (object);
}

static void
nc_window_tophat_class_init (NcWindowTophatClass *klass)
{
  GObjectClass *object_class  = G_OBJECT_CLASS (klass);
  NcWindowClass *parent_class = NC_WINDOW_CLASS (klass);
  
  parent_class->volume        = NC_WINDOW_VOLUME_TOPHAT;
  parent_class->eval_fourier  = &_nc_window_tophat_eval_fourier;
  parent_class->deriv_fourier = &_nc_window_tophat_deriv_fourier;
  parent_class->eval_real     = &_nc_window_tophat_eval_realspace;
  
  object_class->finalize = nc_window_tophat_finalize;
}

