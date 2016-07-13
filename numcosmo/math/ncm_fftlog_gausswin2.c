/***************************************************************************
 *            ncm_fftlog_gausswin2.c
 *
 *  Mon July 21 19:59:38 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/* excerpt from: */
/***************************************************************************
 *            nc_window_gaussian.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_fftlog_gausswin2.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fftlog_gausswin2
 * @title: NcmFftlogGausswin2
 * @short_description: Logarithm fast fourier transform for a kernel given by the square of a Gaussian window function.
 *
 *
 * This object computes the function (see #NcmFftlog)
 * $$Y_n = \int_0^\infty t^{\frac{2\pi i n}{L}} K(t) dt,$$
 * where the kernel is the square of the Gaussian window function $K(t) = W(t)^2$,
 * \begin{equation}
 * W(t) = \exp \left( \frac{-t^2}{2} \right).
 * \end{equation}
 *  
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog_gausswin2.h"
#include "math/ncm_cfg.h"

#include <math.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>

G_DEFINE_TYPE (NcmFftlogGausswin2, ncm_fftlog_gausswin2, NCM_TYPE_FFTLOG);

static void
ncm_fftlog_gausswin2_init (NcmFftlogGausswin2 *gwin2)
{
}

static void
_ncm_fftlog_gausswin2_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_gausswin2_parent_class)->finalize (object);
}

static void _ncm_fftlog_gausswin2_get_Ym (NcmFftlog *fftlog, gpointer Ym_0);

static void
ncm_fftlog_gausswin2_class_init (NcmFftlogGausswin2Class *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->finalize = &_ncm_fftlog_gausswin2_finalize;

  fftlog_class->name        = "gaussian_window_2";
  fftlog_class->get_Ym      = &_ncm_fftlog_gausswin2_get_Ym;
}

static void 
_ncm_fftlog_gausswin2_get_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
  const gdouble twopi_Lt  = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gint Nf           = ncm_fftlog_get_full_size (fftlog);
#ifdef NUMCOSMO_HAVE_FFTW3
  fftw_complex *Ym_base = (fftw_complex *) Ym_0;
  gint i;

  for (i = 0; i < Nf; i++)
  {
    const gint phys_i            = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double a       = twopi_Lt * phys_i * I;
    const complex double A       = a + 0.0/*fftlog->nu*/;
    const complex double onepA_2 = 0.5 * (1.0 + A);
    complex double U;
    gsl_sf_result lngamma_rho, lngamma_theta;

    gsl_sf_lngamma_complex_e (creal (onepA_2), cimag (onepA_2), &lngamma_rho, &lngamma_theta);
    U = 0.5 * cexp (lngamma_rho.val + I * lngamma_theta.val);

    Ym_base[i] = U;
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_gausswin2_new:
 * @lnr0: output center $\ln(r_0)$
 * @lnk0: input center $\ln(k_0)$
 * @Lk: input/output interval size
 * @N: number of knots
 * 
 * Creates a new fftlog Gaussian window squared object.
 * 
 * Returns: (transfer full): a new #NcmFftlogGausswin2
 */
NcmFftlogGausswin2 *
ncm_fftlog_gausswin2_new (gdouble lnr0, gdouble lnk0, gdouble Lk, guint N)
{
  NcmFftlogGausswin2 *fftlog = g_object_new (NCM_TYPE_FFTLOG_GAUSSWIN2, 
                                             "lnr0", lnr0,
                                             "lnk0", lnk0,
                                             "Lk", Lk,
                                             "N", N,
                                             NULL);
  return fftlog;
}
