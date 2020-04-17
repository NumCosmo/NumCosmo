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
 *  <sandro@isoftware.com.br>
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
 * SECTION:ncm_fftlog_sbessel_j
 * @title: NcmFftlogSBesselJ
 * @short_description: Logarithm fast fourier transform for a kernel given by the spatial correlation function multipoles.
 *
 *
 * This object computes the function (see #NcmFftlog)
 * $$Y_n = \int_0^\infty t^{\frac{2\pi i n}{L}} K(t) dt,$$
 * where the kernel are the spherical bessel function of the first kind $K(t) = t^q j_{\ell}(t)$.
 * 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog_sbessel_j.h"
#include "math/ncm_cfg.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmFftlogSBesselJPrivate
{
  guint ell;
  gdouble q;
};

enum
{
  PROP_0,
  PROP_ELL,
  PROP_Q,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFftlogSBesselJ, ncm_fftlog_sbessel_j, NCM_TYPE_FFTLOG);

static void
ncm_fftlog_sbessel_j_init (NcmFftlogSBesselJ *fftlog_jl)
{
  fftlog_jl->priv      = ncm_fftlog_sbessel_j_get_instance_private (fftlog_jl);
  fftlog_jl->priv->ell = 0;
  fftlog_jl->priv->q   = 0.0;
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
    case PROP_Q:
      ncm_fftlog_sbessel_j_set_q (fftlog_jl, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    case PROP_Q:
      g_value_set_double (value, ncm_fftlog_sbessel_j_get_q (fftlog_jl));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fftlog_sbessel_j_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_sbessel_j_parent_class)->finalize (object);
}

static void _ncm_fftlog_sbessel_j_get_Ym (NcmFftlog *fftlog, gpointer Ym_0);

static void
ncm_fftlog_sbessel_j_class_init (NcmFftlogSBesselJClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);
  
  object_class->set_property = &_ncm_fftlog_sbessel_j_set_property;
  object_class->get_property = &_ncm_fftlog_sbessel_j_get_property;
  object_class->finalize     = &_ncm_fftlog_sbessel_j_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ELL,
                                   g_param_spec_uint ("ell",
                                                      NULL,
                                                      "Spherical Bessel integer order",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Q,
                                   g_param_spec_double ("q",
                                                        NULL,
                                                        "Spherical Bessel power",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  fftlog_class->name     = "sbessel_j";
  fftlog_class->get_Ym   = &_ncm_fftlog_sbessel_j_get_Ym;
}

static void 
_ncm_fftlog_sbessel_j_get_Ym (NcmFftlog *fftlog, gpointer Ym_0)
{
  NcmFftlogSBesselJ *fftlog_jl = NCM_FFTLOG_SBESSEL_J (fftlog);
  const gdouble pi_sqrt        = sqrt (M_PI);
  const gdouble twopi_Lt       = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
  const gint Nf                = ncm_fftlog_get_full_size (fftlog);

#ifdef NUMCOSMO_HAVE_FFTW3
  fftw_complex *Ym_base = (fftw_complex *) Ym_0;
  gint i;

  if (fftlog_jl->priv->q == 0.5)
  {
    for (i = 0; i < Nf; i++)
    {
      const gint phys_i             = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double a        = twopi_Lt * phys_i * I;
      const complex double A        = a + 0.5;
      const complex double xup      = 0.5 * (1.0 + 1.0 * fftlog_jl->priv->ell + A);
      const complex double two_x_m1 = cpow (2.0, A - 1.0);
      complex double U;

      gsl_sf_result lngamma_rho_up, lngamma_theta_up;

      gsl_sf_lngamma_complex_e (creal (xup), cimag (xup), &lngamma_rho_up, &lngamma_theta_up);

      U = cexp (2.0 * I * lngamma_theta_up.val);

      Ym_base[i] = pi_sqrt * two_x_m1 * U;
    }
  }
  else
  {
    const gdouble q = fftlog_jl->priv->q;
    for (i = 0; i < Nf; i++)
    {
      const gint phys_i             = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double a        = twopi_Lt * phys_i * I;
      const complex double A        = a + q;
      const complex double xup      = 0.5 * (1.0 + 1.0 * fftlog_jl->priv->ell + A);
      const complex double xdw      = 0.5 * (2.0 + 1.0 * fftlog_jl->priv->ell - A);
      const complex double two_x_m1 = cpow (2.0, A - 1.0);
      complex double U;

      gsl_sf_result lngamma_rho_up, lngamma_theta_up;
      gsl_sf_result lngamma_rho_dw, lngamma_theta_dw;

      gsl_sf_lngamma_complex_e (creal (xup), cimag (xup), &lngamma_rho_up, &lngamma_theta_up);
      gsl_sf_lngamma_complex_e (creal (xdw), cimag (xdw), &lngamma_rho_dw, &lngamma_theta_dw);

      U = cexp ((lngamma_rho_up.val - lngamma_rho_dw.val) + I * (lngamma_theta_up.val - lngamma_theta_dw.val));

      /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g\n", lngamma_rho_up.val, lngamma_rho_dw.val, lngamma_theta_up.val, lngamma_theta_dw.val);*/

      Ym_base[i] = pi_sqrt * two_x_m1 * U;
    }
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_sbessel_j_new:
 * @ell: Spherical Bessel Integer order
 * @lnr0: output center $\ln(r_0)$
 * @lnk0: input center $\ln(k_0)$
 * @Lk: input/output interval size
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
 * ncm_fftlog_set_ell:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 * @ell: Spherical Bessel integer order $\ell$
 * 
 * Sets @ell as the Spherical Bessel integer order $\ell$.
 * 
 */
void 
ncm_fftlog_sbessel_j_set_ell (NcmFftlogSBesselJ *fftlog_jl, const guint ell)
{
  if (fftlog_jl->priv->ell != ell)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jl);

    fftlog_jl->priv->ell = ell;
    fftlog->prepared     = FALSE;
    fftlog->evaluated    = FALSE;
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
  return fftlog_jl->priv->ell;
}

/**
 * ncm_fftlog_set_q:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 * @ell: Spherical Bessel integer order $\ell$
 * 
 * Sets @q as the Spherical Bessel power $q$.
 * 
 */
void 
ncm_fftlog_sbessel_j_set_q (NcmFftlogSBesselJ *fftlog_jl, const gdouble q)
{
  if (fftlog_jl->priv->q != q)
  {
    NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jl);

    fftlog_jl->priv->q = q;
    fftlog->prepared   = FALSE;
    fftlog->evaluated  = FALSE;
  }
}

/**
 * ncm_fftlog_sbessel_j_get_q:
 * @fftlog_jl: a #NcmFftlogSBesselJ
 * 
 * Returns: the current Spherical Bessel power $q$.
 */
gdouble 
ncm_fftlog_sbessel_j_get_q (NcmFftlogSBesselJ *fftlog_jl)
{
  return fftlog_jl->priv->q;
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
  NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jl);

  gint signp = 0;
  
  const gdouble lnk0      = ncm_fftlog_get_lnk0 (fftlog);
  const gdouble Lk        = ncm_fftlog_get_length (fftlog);
  const gdouble ell       = fftlog_jl->priv->ell;
  const gdouble lnc0      = (ell == 0) ? 0.0 : ((ell - 1.0) * Lk + 2.0 * (ell + 1.0) * M_LN2 - ncm_c_lnpi () + 2.0 * lgamma_r (1.5 + ell, &signp)) / (2.0 * (1.0 + ell));
  const gdouble best_lnr0 = - lnk0 + lnc0;

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
  NcmFftlog *fftlog = NCM_FFTLOG (fftlog_jl);

  gint signp = 0;
  
  const gdouble lnr0      = ncm_fftlog_get_lnr0 (fftlog);
  const gdouble Lk        = ncm_fftlog_get_length (fftlog);
  const gdouble ell       = fftlog_jl->priv->ell;
  const gdouble lnc0      = (ell == 0) ? 0.0 : ((ell - 1.0) * Lk + 2.0 * (ell + 1.0) * M_LN2 - ncm_c_lnpi () + 2.0 * lgamma_r (1.5 + ell, &signp)) / (2.0 * (1.0 + ell));
  const gdouble best_lnk0 = - lnr0 + lnc0;

  ncm_fftlog_set_lnk0 (fftlog, best_lnk0);
}
