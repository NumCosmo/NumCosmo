/***************************************************************************
 *            ncm_fftlog_sbessel_j.c
 *
 *  Wed July 19 10:00:26 2017
 *  Copyright  2017  Fernando de Simoni
 *  <fernando.saliby@gmail.com>
 ****************************************************************************/

/***************************************************************************
 *            ncm_fftlog_sbessel_j.c
 *
 *  Sat September 02 18:11:00 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_fftlog_sbessel_j.c
 *
 * Copyright (C) 2017 - Fernando de Simoni
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmFftlogSBesselJ:
 *
 * Logarithm fast fourier transform for a kernel given by the spatial correlation
 * function multipoles.
 *
 * This object computes the coefficients (see #NcmFftlog)
 * $$Y_n = \int_0^\infty t^{A_n} j_\ell(t)\,\mathrm{d}t, \qquad A_n = b + \frac{2\pi i n}{L_T},$$
 * with $b$ the bias (#NcmFftlog:bias), $L_T$ the full period and $j_\ell$ the spherical
 * Bessel function of the first kind, of integer order $\ell$ (#NcmFftlogSBesselJ:ell).
 * The integral converges for $-\ell - 1 < b < 1$ and is
 * $$Y_n = \sqrt{\pi}\, 2^{A_n - 1} \frac{\Gamma\left(\frac{1 + \ell + A_n}{2}\right)}{\Gamma\left(\frac{2 + \ell - A_n}{2}\right)}.$$
 *
 * The spatial correlation function multipoles, $\xi_{\ell}^{(n)}(r)$, can be defined as
 * (see [Matsubara (2004)](https://arxiv.org/abs/astro-ph/0408349) [[arXiv](https://arxiv.org/abs/astro-ph/0408349)])
 *
 * \begin{equation}\label{eq:xi_multipoles}
 * \xi_{\ell}^{(n)}(r) = \frac{(-1)^{n+\ell}}{r^{2n-\ell}} \int_{0}^{\infty} \frac{\mathrm{d} k}{2\pi^2} \frac{k^2}{k^{2n-\ell}} j_{\ell}(kr) P(k) \,\, .
 * \end{equation}
 * Where, $P(k)$ is the power spectrum (see #NcmPowspec). The integral in it is the
 * transform of $F(k) = k^{2 - 2n + \ell} P(k) / (2\pi^2)$ with this kernel. The bias does
 * not change the result: it only chooses which function, $F(k) k^{-b}$, the discrete
 * transform represents; ncm_fftlog_get_best_bias() chooses one from the log-slopes of $F$
 * at the ends of the interval.
 *
 * The #NcmPowspecCorr3d object already evaluates Eq. \eqref{eq:xi_multipoles}
 * for the case of the monopole, $n=\ell=0$, with support for redshift evolution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fftlog/ncm_fftlog_sbessel_j.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#ifdef HAVE_ACB_H
#ifdef HAVE_FLINT_ACB_H
#include <flint/acb.h>
#else /* HAVE_FLINT_ACB_H */
#include <acb.h>
#endif /* HAVE_FLINT_ACB_H */
#endif /* HAVE_ACB_H */
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmFftlogSBesselJPrivate
{
  guint ell;
} NcmFftlogSBesselJPrivate;

enum
{
  PROP_0,
  PROP_ELL,
  PROP_SIZE,
};

struct _NcmFftlogSBesselJ
{
  NcmFftlog parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFftlogSBesselJ, ncm_fftlog_sbessel_j, NCM_TYPE_FFTLOG)

static void
ncm_fftlog_sbessel_j_init (NcmFftlogSBesselJ *fftlog_jl)
{
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);

  self->ell = 0;
}

static void
_ncm_fftlog_sbessel_j_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFftlogSBesselJ *fftlog_jl = NCM_FFTLOG_SBESSEL_J (object);

  g_return_if_fail (NCM_IS_FFTLOG_SBESSEL_J (object));

  switch (prop_id)
  {
    case PROP_ELL:
      ncm_fftlog_sbessel_j_set_ell (fftlog_jl, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fftlog_sbessel_j_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFftlogSBesselJ *fftlog_jl = NCM_FFTLOG_SBESSEL_J (object);

  g_return_if_fail (NCM_IS_FFTLOG_SBESSEL_J (object));

  switch (prop_id)
  {
    case PROP_ELL:
      g_value_set_uint (value, ncm_fftlog_sbessel_j_get_ell (fftlog_jl));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fftlog_sbessel_j_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_sbessel_j_parent_class)->finalize (object);
}

static void _ncm_fftlog_sbessel_j_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0);
static void _ncm_fftlog_sbessel_j_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max);

static void
ncm_fftlog_sbessel_j_class_init (NcmFftlogSBesselJClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->set_property = &_ncm_fftlog_sbessel_j_set_property;
  object_class->get_property = &_ncm_fftlog_sbessel_j_get_property;
  object_class->finalize     = &_ncm_fftlog_sbessel_j_finalize;

  /**
   * NcmFftlogSBesselJ:ell:
   *
   * The spherical Bessel integer order.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL,
                                   g_param_spec_uint ("ell",
                                                      NULL,
                                                      "Spherical Bessel integer order",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  fftlog_class->name           = "sbessel_j";
  fftlog_class->compute_Ym     = &_ncm_fftlog_sbessel_j_compute_Ym;
  fftlog_class->get_bias_range = &_ncm_fftlog_sbessel_j_get_bias_range;
}

/* At b = 1/2 the two Gamma arguments are complex conjugates, so the ratio is the phase
 * exp (2 i arg Gamma) and one lngamma per mode is enough. */
static void
_ncm_fftlog_sbessel_j_compute_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
  NcmFftlogSBesselJ *fftlog_jl          = NCM_FFTLOG_SBESSEL_J (fftlog);
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);

  const gdouble pi_sqrt  = sqrt (M_PI);
  const gdouble twopi_Lt = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gdouble bias     = ncm_fftlog_get_bias (fftlog);
  const gint Nf          = ncm_fftlog_get_full_size (fftlog);

  fftw_complex *Ym_base = (fftw_complex *) Ym_0;
  gint i;

  if (bias == 0.5)
  {
    for (i = 0; i < Nf; i++)
    {
      const gint phys_i             = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double A        = twopi_Lt * phys_i * I + 0.5;
      const complex double xup      = 0.5 * (1.0 + 1.0 * self->ell + A);
      const complex double two_x_m1 = cpow (2.0, A - 1.0);
      gsl_sf_result lngamma_rho_up, lngamma_theta_up;

      gsl_sf_lngamma_complex_e (creal (xup), cimag (xup), &lngamma_rho_up, &lngamma_theta_up);

      Ym_base[i] = pi_sqrt * two_x_m1 * cexp (2.0 * I * lngamma_theta_up.val);
    }
  }
  else
  {
    for (i = 0; i < Nf; i++)
    {
      const gint phys_i             = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double A        = twopi_Lt * phys_i * I + bias;
      const complex double xup      = 0.5 * (1.0 + 1.0 * self->ell + A);
      const complex double xdw      = 0.5 * (2.0 + 1.0 * self->ell - A);
      const complex double two_x_m1 = cpow (2.0, A - 1.0);
      gsl_sf_result lngamma_rho_up, lngamma_theta_up;
      gsl_sf_result lngamma_rho_dw, lngamma_theta_dw;

      gsl_sf_lngamma_complex_e (creal (xup), cimag (xup), &lngamma_rho_up, &lngamma_theta_up);
      gsl_sf_lngamma_complex_e (creal (xdw), cimag (xdw), &lngamma_rho_dw, &lngamma_theta_dw);

      Ym_base[i] = pi_sqrt * two_x_m1 * cexp ((lngamma_rho_up.val - lngamma_rho_dw.val) + I * (lngamma_theta_up.val - lngamma_theta_dw.val));
    }
  }
}

/* t^b j_l(t) goes as t^(b + l) at t -> 0 and oscillates as t^(b - 1) at t -> infinity. */
static void
_ncm_fftlog_sbessel_j_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max)
{
  NcmFftlogSBesselJ *fftlog_jl          = NCM_FFTLOG_SBESSEL_J (fftlog);
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);

  *bias_min = -1.0 - self->ell;
  *bias_max = 1.0;
}

/**
 * ncm_fftlog_sbessel_j_new:
 * @ell: Spherical Bessel Integer order
 * @lnr0: output center $\ln(r_0)$
 * @lnk0: input center $\ln(k_0)$
 * @Lk: length $L$ of the fundamental interval in $\ln k$, see #NcmFftlog:Lk
 * @N: number of knots
 *
 * Creates a new fftlog Spherical Bessel J object.
 *
 * Returns: (transfer full): a new #NcmFftlogSBesselJ
 */
NcmFftlogSBesselJ *
ncm_fftlog_sbessel_j_new (guint ell, gdouble lnr0, gdouble lnk0, gdouble Lk, guint N)
{
  NcmFftlogSBesselJ *fftlog_jl = g_object_new (NCM_TYPE_FFTLOG_SBESSEL_J,
                                               "ell",  ell,
                                               "lnr0", lnr0,
                                               "lnk0", lnk0,
                                               "Lk",   Lk,
                                               "N",    N,
                                               NULL);

  return fftlog_jl;
}

/**
 * ncm_fftlog_sbessel_j_set_ell:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 * @ell: Spherical Bessel integer order $\ell$
 *
 * Sets @ell as the Spherical Bessel integer order $\ell$.
 *
 */
void
ncm_fftlog_sbessel_j_set_ell (NcmFftlogSBesselJ *fftlog_jl, const guint ell)
{
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);

  if (self->ell != ell)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jl);

    self->ell = ell;
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_sbessel_j_get_ell:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 *
 * Returns: the current Spherical Bessel integer order $\ell$.
 */
guint
ncm_fftlog_sbessel_j_get_ell (NcmFftlogSBesselJ *fftlog_jl)
{
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);

  return self->ell;
}

/**
 * ncm_fftlog_sbessel_j_set_best_lnr0:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 *
 * Sets the value of $\ln(r_0)$ which gives the best results for
 * the transformation based on the current value of $\ln(k_0)$,
 * this is based in the rule of thumb $\mathrm{max}_{x^*}(j_l)$
 * where $ x^* \approx l + 1$.
 *
 */
void
ncm_fftlog_sbessel_j_set_best_lnr0 (NcmFftlogSBesselJ *fftlog_jl)
{
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);
  NcmFftlog *fftlog                     = NCM_FFTLOG (fftlog_jl);

  gint signp = 0;

  const gdouble lnk0      = ncm_fftlog_get_lnk0 (fftlog);
  const gdouble Lk        = ncm_fftlog_get_length (fftlog);
  const gdouble ell       = self->ell;
  const gdouble lnc0      = (ell == 0) ? 0.0 : ((ell - 1.0) * Lk + 2.0 * (ell + 1.0) * M_LN2 - ncm_c_lnpi () + 2.0 * lgamma_r (1.5 + ell, &signp)) / (2.0 * (1.0 + ell));
  const gdouble best_lnr0 = -lnk0 + lnc0;

  ncm_fftlog_set_lnr0 (fftlog, best_lnr0);
}

/**
 * ncm_fftlog_sbessel_j_set_best_lnk0:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 *
 * Sets the value of $\ln(k_0)$ which gives the best results for
 * the transformation based on the current value of $\ln(r_0)$,
 * this is based in the rule of thumb $\mathrm{max}_{x^*}(j_l)$
 * where $ x^* \approx l + 1$.
 *
 */
void
ncm_fftlog_sbessel_j_set_best_lnk0 (NcmFftlogSBesselJ *fftlog_jl)
{
  NcmFftlogSBesselJPrivate * const self = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);
  NcmFftlog *fftlog                     = NCM_FFTLOG (fftlog_jl);

  gint signp = 0;

  const gdouble lnr0      = ncm_fftlog_get_lnr0 (fftlog);
  const gdouble Lk        = ncm_fftlog_get_length (fftlog);
  const gdouble ell       = self->ell;
  const gdouble lnc0      = (ell == 0) ? 0.0 : ((ell - 1.0) * Lk + 2.0 * (ell + 1.0) * M_LN2 - ncm_c_lnpi () + 2.0 * lgamma_r (1.5 + ell, &signp)) / (2.0 * (1.0 + ell));
  const gdouble best_lnk0 = -lnr0 + lnc0;

  ncm_fftlog_set_lnk0 (fftlog, best_lnk0);
}

