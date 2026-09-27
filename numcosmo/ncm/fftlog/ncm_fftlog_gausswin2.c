/***************************************************************************
 *            ncm_fftlog_gausswin2.c
 *
 *  Mon July 21 19:59:38 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
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
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmFftlogGausswin2:
 *
 * Logarithm fast fourier transform for a kernel given by the square of a Gaussian
 * window function.
 *
 *
 * This object computes the coefficients (see #NcmFftlog)
 * $$Y_n = \int_0^\infty t^{A_n} K(t)\,\mathrm{d}t, \qquad A_n = b + \frac{2\pi i n}{L_T},$$
 * with $b$ the bias (#NcmFftlog:bias) and $L_T$ the full period, where the kernel is the
 * square of the Gaussian window function $K(t) = W(t)^2$,
 * \begin{equation}
 * W(t) = \exp \left( \frac{-t^2}{2} \right).
 * \end{equation}
 * The integral converges for $b > -1$ and is
 * $Y_n = \Gamma\left(\frac{1 + A_n}{2}\right) / 2$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fftlog/ncm_fftlog_gausswin2.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#include <fftw3.h>

#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmFftlogGausswin2
{
  NcmFftlog parent_instance;
};

G_DEFINE_TYPE (NcmFftlogGausswin2, ncm_fftlog_gausswin2, NCM_TYPE_FFTLOG)

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

static void _ncm_fftlog_gausswin2_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0);
static void _ncm_fftlog_gausswin2_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max);

static void
ncm_fftlog_gausswin2_class_init (NcmFftlogGausswin2Class *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->finalize = &_ncm_fftlog_gausswin2_finalize;

  fftlog_class->name           = "gaussian_window_2";
  fftlog_class->compute_Ym     = &_ncm_fftlog_gausswin2_compute_Ym;
  fftlog_class->get_bias_range = &_ncm_fftlog_gausswin2_get_bias_range;
}

static void
_ncm_fftlog_gausswin2_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
  const gdouble twopi_Lt = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gdouble bias     = ncm_fftlog_get_bias (fftlog);
  const gint Nf          = ncm_fftlog_get_full_size (fftlog);

  fftw_complex *Ym_base = (fftw_complex *) Ym_0;
  gint i;

  for (i = 0; i < Nf; i++)
  {
    const gint phys_i            = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double a       = twopi_Lt * phys_i * I;
    const complex double A       = a + bias;
    const complex double onepA_2 = 0.5 * (1.0 + A);
    complex double U;
    gsl_sf_result lngamma_rho, lngamma_theta;

    gsl_sf_lngamma_complex_e (creal (onepA_2), cimag (onepA_2), &lngamma_rho, &lngamma_theta);
    U = 0.5 * cexp (lngamma_rho.val + I * lngamma_theta.val);

    Ym_base[i] = U;
  }
}

static void
_ncm_fftlog_gausswin2_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max)
{
  *bias_min = -1.0;
  *bias_max = GSL_POSINF;
}

/**
 * ncm_fftlog_gausswin2_new:
 * @lnr0: output center $\ln(r_0)$
 * @lnk0: input center $\ln(k_0)$
 * @Lk: length $L$ of the fundamental interval in $\ln k$, see #NcmFftlog:Lk
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

