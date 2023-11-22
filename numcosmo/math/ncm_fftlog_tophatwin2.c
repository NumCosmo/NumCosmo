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
 * SECTION:ncm_fftlog_tophatwin2
 * @title: NcmFftlogTophatwin2
 * @short_description: Logarithm fast fourier transform for a kernel given by the square of the spherical Bessel function of order one.
 * @stability: Stable
 * @include: numcosmo/math/ncm_fftlog_tophatwin2.h
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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */

#ifdef HAVE_ACB_H
#include <acb.h>
#endif /* HAVE_ACB_H */
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

static void _ncm_fftlog_tophatwin2_get_Ym (NcmFftlog *fftlog, gpointer Ym_0);

static void
ncm_fftlog_tophatwin2_class_init (NcmFftlogTophatwin2Class *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->finalize = &_ncm_fftlog_tophatwin2_finalize;

  fftlog_class->name   = "tophat_window_2";
  fftlog_class->get_Ym = &_ncm_fftlog_tophatwin2_get_Ym;
}

static void
_ncm_fftlog_tophatwin2_get_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
#ifdef NUMCOSMO_HAVE_FFTW3
#if defined (HAVE_ACB_H) && defined (NCM_FFTLOG_USE_ACB)
  const guint prec      = 120;
  const gint Nf         = ncm_fftlog_get_full_size (fftlog);
  fftw_complex *Ym_base = (fftw_complex *) Ym_0;
  acb_t twopi_Lt, Lt, a_i, Un_i, Ud_i, pi, two_a_i, a_i_pi_2, a_i_m3;
  gint i;

  acb_init (pi);
  acb_init (Lt);
  acb_init (twopi_Lt);
  acb_init (a_i);
  acb_init (Un_i);
  acb_init (Ud_i);
  acb_init (two_a_i);
  acb_init (a_i_pi_2);
  acb_init (a_i_m3);

  acb_set_d (Lt, ncm_fftlog_get_full_length (fftlog));

  acb_const_pi (pi, prec);

  acb_set (twopi_Lt, pi);
  acb_mul_ui (twopi_Lt, twopi_Lt, 2, prec);
  acb_div (twopi_Lt, twopi_Lt, Lt, prec);

  i = 0;
  {
    const gint phys_i = ncm_fftlog_get_mode_index (fftlog, i);

    acb_mul_si (a_i, twopi_Lt, phys_i, prec); /* a_i = twopi_Lt * phys_i */
    acb_mul_onei (a_i, a_i);                  /* a_i = a_i * I */

    acb_sub_ui (Un_i, a_i, 1, prec);  /* Un_i = a_i - 1   */
    acb_mul_si (Un_i, Un_i, 3, prec); /* Un_i = 3 * Un_i  */
    acb_mul (Un_i, Un_i, pi, prec);   /* Un_i = pi * Un_i */

    acb_sub_ui (Ud_i, a_i, 5, prec);       /* Ud_i    = a_i - 5 */
    acb_set_ui (two_a_i, 2);               /* two_a_i = 2 */
    acb_pow (two_a_i, two_a_i, a_i, prec); /* two_a_i = pow (two_a_i, a_i) */

    acb_mul (Ud_i, Ud_i, two_a_i, prec); /* Ud_i = Ud_i * two_a_i */
    acb_div (Un_i, Un_i, Ud_i, prec);    /* Un_i = Un_i / Ud_i */

    Ym_base[i] = ncm_acb_get_complex (Un_i);
    /*printf ("%d % 20.15g % 20.15g\n", i, creal (Ym_base[i]), cimag (Ym_base[i]));*/
  }

  for (i = 1; i < Nf; i++)
  {
    const gint phys_i = ncm_fftlog_get_mode_index (fftlog, i);

    acb_mul_si (a_i, twopi_Lt, phys_i, prec); /* a_i = twopi_Lt * phys_i */
    acb_mul_onei (a_i, a_i);                  /* a_i = a_i * I */

    acb_sub_ui (Un_i, a_i, 1, prec);    /* Un_i = a_i - 1   */
    acb_mul_si (Un_i, Un_i, -36, prec); /* Un_i = -36 * Un_i  */

    acb_sub_ui (Ud_i, a_i, 5, prec);       /* Ud_i    = a_i - 5 */
    acb_set_ui (two_a_i, 2);               /* two_a_i = 2 */
    acb_pow (two_a_i, two_a_i, a_i, prec); /* two_a_i = pow (two_a_i, a_i) */

    acb_mul (Ud_i, Ud_i, two_a_i, prec); /* Ud_i = Ud_i * two_a_i */
    acb_div (Un_i, Un_i, Ud_i, prec);    /* Un_i = Un_i / Ud_i */

    acb_mul (a_i_pi_2, a_i, pi, prec);        /* a_i_pi_2 = a_i * pi */
    acb_div_ui (a_i_pi_2, a_i_pi_2, 2, prec); /* a_i_pi_2 = a_i_pi_2 / 2 */
    acb_sub_ui (a_i_m3, a_i, 3, prec);        /* a_i_m3   = a_i - 3 */

    acb_sin (a_i_pi_2, a_i_pi_2, prec); /* a_i_pi_2 = sin (a_i_pi_2) */
    acb_gamma (a_i_m3, a_i_m3, prec);   /* a_i_m3 = gamma (a_i_m3) */

    acb_mul (Un_i, Un_i, a_i_pi_2, prec); /* Un_i = Un_i * a_i_pi_2 */
    acb_mul (Un_i, Un_i, a_i_m3, prec);   /* Un_i = Un_i * a_i_m3 */

    Ym_base[i] = ncm_acb_get_complex (Un_i);
    /*printf ("%d % 20.15g % 20.15g\n", i, creal (Ym_base[i]), cimag (Ym_base[i]));*/
  }

  acb_clear (pi);
  acb_clear (Lt);
  acb_clear (twopi_Lt);
  acb_clear (a_i);
  acb_clear (Un_i);
  acb_clear (Ud_i);
  acb_clear (two_a_i);
  acb_clear (a_i_pi_2);
  acb_clear (a_i_m3);

#else /* HAVE_ACB_H */
  const gdouble twopi_Lt = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gint Nf          = ncm_fftlog_get_full_size (fftlog);
  fftw_complex *Ym_base  = (fftw_complex *) Ym_0;
  gint i;

  i = 0;
  {
    const gint phys_i      = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double a = twopi_Lt * phys_i * I;
    complex double U       = -36.0 * (a - 1.0) / (cexp (M_LN2 * a) * (a - 5.0));

    U *= -M_PI / 12.0;

    Ym_base[i] = U;
    /*printf ("%d % 20.15g % 20.15g\n", i, creal (Ym_base[i]), cimag (Ym_base[i]));*/
  }

  for (i = 1; i < Nf; i++)
  {
    const gint phys_i          = ncm_fftlog_get_mode_index (fftlog, i);
    const complex double a     = twopi_Lt * phys_i * I;
    complex double U           = -36.0 * (a - 1.0) / (cexp (M_LN2 * a) * (a - 5.0));
    const complex double Api_2 = a * M_PI * 0.5;
    gsl_sf_result lngamma_rho, lngamma_theta, lnsin_rho, lnsin_theta;

    gsl_sf_lngamma_complex_e (creal (a) - 3.0, cimag (a), &lngamma_rho, &lngamma_theta);
    gsl_sf_complex_logsin_e (creal (Api_2), cimag (Api_2), &lnsin_rho, &lnsin_theta);
    U *= cexp (lngamma_rho.val + lnsin_rho.val + I * (lngamma_theta.val + lnsin_theta.val));

    Ym_base[i] = U;
    /*printf ("%d % 20.15g % 20.15g\n", i, creal (Ym_base[i]), cimag (Ym_base[i]));*/
  }

#endif /* HAVE_ACB_H */
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

