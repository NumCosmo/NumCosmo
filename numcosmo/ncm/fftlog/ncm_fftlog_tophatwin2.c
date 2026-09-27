/***************************************************************************
 *            ncm_fftlog_tophatwin2.c
 *
 *  Mon July 21 19:59:38 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/* excerpt from: */

/***************************************************************************
 *            nc_window_tophat.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/

/*
 * ncm_fftlog_tophatwin2.c
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
 * NcmFftlogTophatwin2:
 *
 * Logarithm fast fourier transform for a kernel given by the square of the spherical
 * Bessel function of order one.
 *
 *
 * This object computes the coefficients (see #NcmFftlog)
 * $$Y_n = \int_0^\infty t^{A_n} K(t)\,\mathrm{d}t, \qquad A_n = b + \frac{2\pi i n}{L_T},$$
 * with $b$ the bias (#NcmFftlog:bias) and $L_T$ the full period, where the kernel is the
 * square of the top hat window function in the Fourier space $K(t) = W(t)^2$,
 * \begin{eqnarray}
 * W(t) &=& \frac{3}{t^3}(\sin t - t \cos t) \\
 * &=& \frac{3}{t} j_1(t),
 * \end{eqnarray}
 * and $j_\nu(t)$ is the spherical Bessel function of the first kind. The integral
 * converges for $-1 < b < 3$, since $K(t) \to 1$ at $t \to 0$ and falls as $t^{-4}$ at
 * $t \to \infty$, and is
 * $$Y_n = \frac{9\sqrt{\pi}\,\Gamma\left(\frac{1 + A_n}{2}\right)}{(A_n - 3)(A_n - 5)\,\Gamma\left(2 - \frac{A_n}{2}\right)},$$
 * which is regular over that whole range.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fftlog/ncm_fftlog_tophatwin2.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmFftlogTophatwin2
{
  NcmFftlog parent_instance;
};

G_DEFINE_TYPE (NcmFftlogTophatwin2, ncm_fftlog_tophatwin2, NCM_TYPE_FFTLOG)

static void
ncm_fftlog_tophatwin2_init (NcmFftlogTophatwin2 *j1pow2)
{
}

static void
_ncm_fftlog_tophatwin2_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_tophatwin2_parent_class)->finalize (object);
}

static void _ncm_fftlog_tophatwin2_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0);
static void _ncm_fftlog_tophatwin2_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max);

static void
ncm_fftlog_tophatwin2_class_init (NcmFftlogTophatwin2Class *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->finalize = &_ncm_fftlog_tophatwin2_finalize;

  fftlog_class->name           = "tophat_window_2";
  fftlog_class->compute_Ym     = &_ncm_fftlog_tophatwin2_compute_Ym;
  fftlog_class->get_bias_range = &_ncm_fftlog_tophatwin2_get_bias_range;
}

/* The textbook form, -36 (A - 1) Gamma (A - 3) sin (pi A / 2) / (2^A (A - 5)), has
 * removable singularities at A = 0, 1, 2; the duplication and reflection formulas turn it
 * into the ratio of Gamma functions of the class documentation, whose arguments keep
 * positive real parts for -1 < b < 3. */
static void
_ncm_fftlog_tophatwin2_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
  const gdouble twopi_Lt = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gdouble bias     = ncm_fftlog_get_bias (fftlog);
  const gint Nf          = ncm_fftlog_get_full_size (fftlog);
  const gdouble c        = 9.0 * sqrt (M_PI);
  fftw_complex *Ym_base  = (fftw_complex *) Ym_0;
  gint i;

  for (i = 0; i < Nf; i++)
  {
    const gint phys_i       = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double A  = bias + twopi_Lt * phys_i * I;
    const complex double up = 0.5 * (1.0 + A);
    const complex double dw = 2.0 - 0.5 * A;
    gsl_sf_result lng_up_rho, lng_up_theta, lng_dw_rho, lng_dw_theta;

    gsl_sf_lngamma_complex_e (creal (up), cimag (up), &lng_up_rho, &lng_up_theta);
    gsl_sf_lngamma_complex_e (creal (dw), cimag (dw), &lng_dw_rho, &lng_dw_theta);

    Ym_base[i] = c * cexp ((lng_up_rho.val - lng_dw_rho.val) + I * (lng_up_theta.val - lng_dw_theta.val)) / ((A - 3.0) * (A - 5.0));
  }
}

static void
_ncm_fftlog_tophatwin2_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max)
{
  *bias_min = -1.0;
  *bias_max = 3.0;
}

/**
 * ncm_fftlog_tophatwin2_new:
 * @lnr0: output center $\ln(r_0)$
 * @lnk0: input center $\ln(k_0)$
 * @Lk: input/output interval size
 * @N: number of knots
 *
 * Creates a new fftlog top hat window squared object.
 *
 * Returns: (transfer full): a new #NcmFftlogTophatwin2
 */
NcmFftlogTophatwin2 *
ncm_fftlog_tophatwin2_new (gdouble lnr0, gdouble lnk0, gdouble Lk, guint N)
{
  NcmFftlogTophatwin2 *fftlog = g_object_new (NCM_TYPE_FFTLOG_TOPHATWIN2,
                                              "lnr0", lnr0,
                                              "lnk0", lnk0,
                                              "Lk", Lk,
                                              "N", N,
                                              NULL);

  return fftlog;
}

