/***************************************************************************
 *            ncm_fftlog_tophatwin2.c
 *
 *  Mon July 21 19:59:38 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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
 * SECTION:ncm_fftlog_tophatwin2
 * @title: NcmFftlogTophatwin2
 * @short_description: Logarithm fast fourier transform for a kernel given by the square of the spherical Bessel function of order one.
 *
 *
 * This object computes the function (see #NcmFftlog)
 * $$Y_n = \int_0^\infty t^{\frac{2\pi i n}{L}} K(t) dt,$$
 * where the kernel is the square of the top hat window function in the Fourier space $K(t) = W(t)^2$,
 * \begin{eqnarray}
 * W(t) &=& \frac{3}{t^3}(\sin t - t \cos t) \\
 * &=& \frac{3}{t} j_1(t),
 * \end{eqnarray}
 * and $j_\nu(t)$ is the spherical Bessel function of the first kind.
 *  
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog_tophatwin2.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#include <math.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>

G_DEFINE_TYPE (NcmFftlogTophatwin2, ncm_fftlog_tophatwin2, NCM_TYPE_FFTLOG);

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

static void _ncm_fftlog_tophatwin2_get_Ym (NcmFftlog *fftlog);
static void _ncm_fftlog_tophatwin2_generate_Gr (NcmFftlog *fftlog);

static void
ncm_fftlog_tophatwin2_class_init (NcmFftlogTophatwin2Class *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->finalize = &_ncm_fftlog_tophatwin2_finalize;

  fftlog_class->name        = "tophat_window_2";
  fftlog_class->ncomp       = 2;
  fftlog_class->get_Ym      = &_ncm_fftlog_tophatwin2_get_Ym;
  fftlog_class->generate_Gr = &_ncm_fftlog_tophatwin2_generate_Gr;
}

static void 
_ncm_fftlog_tophatwin2_get_Ym (NcmFftlog *fftlog)
{
  const gdouble twopi_Lt  = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gint Nf           = ncm_fftlog_get_full_size (fftlog);
  gint i;

#ifdef NUMCOSMO_HAVE_FFTW3
  i = 0;
  {
    const gint phys_i          = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double a     = twopi_Lt * phys_i * I;
    const complex double A     = a + 0.0/*fftlog->nu*/;
    complex double U           = - 36.0 * (A - 1.0) / (cexp (M_LN2 * A) * (A - 5.0));

    U *= - M_PI / 12.0;

    fftlog->Ym[0][i] = U * cexp (- a * (fftlog->lnk0 + fftlog->lnr0));
    fftlog->Ym[1][i] = -(1.0 + a) * fftlog->Ym[0][i];    
  }
  for (i = 1; i < Nf; i++)
  {
    const gint phys_i          = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double a     = twopi_Lt * phys_i * I;
    const complex double A     = a + 0.0/*fftlog->nu*/;
    complex double U           = - 36.0 * (A - 1.0) / (cexp (M_LN2 * A) * (A - 5.0));
    const complex double Api_2 = A * M_PI * 0.5;
    gsl_sf_result lngamma_rho, lngamma_theta, lnsin_rho, lnsin_theta;

    gsl_sf_lngamma_complex_e (creal (A) - 3.0, cimag (A), &lngamma_rho, &lngamma_theta);
    gsl_sf_complex_logsin_e (creal (Api_2), cimag (Api_2), &lnsin_rho, &lnsin_theta);
    U *= cexp (lngamma_rho.val + lnsin_rho.val + I * (lngamma_theta.val + lnsin_theta.val));

    fftlog->Ym[0][i] = U * cexp (- a * (fftlog->lnk0 + fftlog->lnr0));
    fftlog->Ym[1][i] = -(1.0 + a) * fftlog->Ym[0][i];
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

static void 
_ncm_fftlog_tophatwin2_generate_Gr (NcmFftlog *fftlog)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  const gint N          = ncm_fftlog_get_size (fftlog);
  const gdouble norma   = ncm_fftlog_get_norma (fftlog);
  gint i;
  
  for (i = 0; i < N; i++)
  {
    const gint phys_i        = i - fftlog->N_2;
    const gdouble lnr        = fftlog->lnr0 + phys_i * fftlog->Lk_N;

    const gdouble rGr        = fabs (creal (fftlog->Gr[0][i + fftlog->pad]) / norma);
    const gdouble lnGr       = log (rGr) - lnr;
    const gdouble dlnGr_dlnr = creal (fftlog->Gr[1][i + fftlog->pad]) / (rGr * norma);

    ncm_vector_set (fftlog->lnr_vec, i, lnr);
    ncm_vector_set (fftlog->Gr_vec[0], i, lnGr);
    ncm_vector_set (fftlog->Gr_vec[1], i, dlnGr_dlnr);
  }

#endif /* NUMCOSMO_HAVE_FFTW3 */
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
