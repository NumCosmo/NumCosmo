/***************************************************************************
 *            nc_window_gaussian.c
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
 * SECTION:nc_window_gaussian
 * @title: Gaussian window function
 * @short_description: Provides a #NcWindow of Gaussian type filter.  
 * 
 * This object implements the #NcWindow class for a Gaussian window function.
 * 
 * This function returns the gaussian window function in the real space,
 * \begin{equation}
 * W_G(r, R) = (2 \pi R^2)^{-3/2}\exp \left( \frac{-r^2}{2 R^2} \right).
 * \end{equation}
 * The mass enclosed within the volume selected by this window function is
 * $M_G(R) = (2\pi)^{3/2}\overline{\rho}(z) R^3$, where $\overline{\rho}(z)$
 * is the mean density of the universe at redshift $z$.
 * 
 * This function returns the gaussian window function in the Fourier space,
 * \begin{equation}
 * W_G(k, R) = \exp \left( \frac{-k^2 R^2}{2} \right).
 * \end{equation}
 * 
 * This function returns the derivative with respect to R of the gaussian window function
 * in the real space,
 * \begin{equation}
 * \frac{dW_G(k, R)}{dR} = -k^2 R \exp \left( \frac{-k^2 R^2}{2} \right).
 * \end{equation}
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_window_gaussian.h"
#include <math.h>
#include <gsl/gsl_math.h>

G_DEFINE_TYPE (NcWindowGaussian, nc_window_gaussian, NC_TYPE_WINDOW);

/**
 * nc_window_gaussian_new:
 *  
 * This function returns a #NcWindow with a #NcWindowGaussian implementation.  
 *
 * Returns: A new #NcWindow.
 */
NcWindow *
nc_window_gaussian_new ()
{
  return g_object_new (NC_TYPE_WINDOW_GAUSSIAN, NULL);
}

static gdouble
_nc_window_gaussian_eval_fourier (const NcWindow * wp, const gdouble k, const gdouble R)
{
  gdouble kR = k * R;
  gdouble kR2 = kR * kR;
  gdouble WG = exp (-kR2 / 2.0);

  return WG;
}

static gdouble
_nc_window_gaussian_deriv_fourier (const NcWindow * wp, const gdouble k, const gdouble R)
{
  gdouble kR = k * R;
  gdouble kR2 = kR * kR;
  gdouble k2R = kR2 / R;
  gdouble dWG = -k2R * exp (-kR2 / 2.0);

  return dWG;
}

static gdouble
_nc_window_gaussian_eval_realspace (const NcWindow * wp, const gdouble r, const gdouble R)
{
  gdouble r_R2 = r * r / (R * R);
  gdouble WG_realspace = 1.0 / gsl_pow_3 (sqrt (2 * M_PI * R * R)) * exp (-r_R2 / 2.0);

  return WG_realspace;
}

static void
nc_window_gaussian_init (NcWindowGaussian *nc_window_gaussian)
{
  /* TODO: Add initialization code here */
}

static void
nc_window_gaussian_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_window_gaussian_parent_class)->finalize (object);
}

static void
nc_window_gaussian_class_init (NcWindowGaussianClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcWindowClass* parent_class = NC_WINDOW_CLASS (klass);

  parent_class->volume = NC_WINDOW_VOLUME_GAUSSIAN;
  parent_class->eval_fourier = &_nc_window_gaussian_eval_fourier;
  parent_class->deriv_fourier = &_nc_window_gaussian_deriv_fourier;
  parent_class->eval_real = &_nc_window_gaussian_eval_realspace;

  object_class->finalize = nc_window_gaussian_finalize;
}

