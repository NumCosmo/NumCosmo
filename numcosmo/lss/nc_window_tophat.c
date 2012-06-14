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
 * @title: Top-hat Window Function
 * @short_description: Provides a #NcWindow of top-hat type filter.
 * 
 * This object implements the #NcWindow abstract class for a top-hat window function.
 * See also <link linkend="sec_wf_th">Top-hat</link> for more details.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>


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
_nc_window_tophat_eval_fourier (const NcWindow * wp, const gdouble k, const gdouble R)
{
  gdouble kR = k * R;
  if (kR == 0.0)
    return 1.0;
  return 3.0 * gsl_sf_bessel_j1 (kR) / kR;
} 

static gdouble
_nc_window_tophat_deriv_fourier (const NcWindow * wp, const gdouble k, const gdouble R)
{
  gdouble dWT = -3.0 * gsl_sf_bessel_j2 (k * R) / R;

  return dWT;
}

static gdouble
_nc_window_tophat_eval_realspace (const NcWindow * wp, const gdouble r, const gdouble R)
{
  gdouble WT_realspace;  
  gdouble R3 = R * R * R;

  if (r <= R)
    WT_realspace = 3.0 / (4.0 * M_PI * R3);
  else
    WT_realspace = 0.0;

  return WT_realspace;
} 

static void
nc_window_tophat_init (NcWindowTophat *nc_window_tophat)
{
  /* TODO: Add initialization code here */ 
}

static void
nc_window_tophat_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_window_tophat_parent_class)->finalize (object);
}

static void
nc_window_tophat_class_init (NcWindowTophatClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcWindowClass *parent_class = NC_WINDOW_CLASS (klass);

  parent_class->volume = NC_WINDOW_VOLUME_TOPHAT;
  parent_class->eval_fourier = &_nc_window_tophat_eval_fourier;
  parent_class->deriv_fourier = &_nc_window_tophat_deriv_fourier;
  parent_class->eval_real = &_nc_window_tophat_eval_realspace;

  object_class->finalize = nc_window_tophat_finalize;
}

