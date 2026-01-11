/***************************************************************************
 *            ncm_sbessel_integrator_fftl.c
 *
 *  Thu January 09 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_fftl.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmSBesselIntegratorFFTL:
 *
 * FFT-Legendre based spherical Bessel function integrator.
 *
 * This class implements integration of functions multiplied by spherical
 * Bessel functions using FFT and Legendre polynomials.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_integrator_fftl.h"
#include "math/ncm_sf_sbessel.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmSBesselIntegratorFFTL
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  gdouble oversample;
  guint Ny;
  gdouble *f_samp;
  fftw_complex *f_fft;
  fftw_plan plan_forward;
};

enum
{
  PROP_0,
  PROP_OVERSAMPLE,
};

G_DEFINE_TYPE (NcmSBesselIntegratorFFTL, ncm_sbessel_integrator_fftl, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_fftl_init (NcmSBesselIntegratorFFTL *sbilf)
{
  sbilf->oversample   = 2.0;
  sbilf->Ny           = 0;
  sbilf->f_samp       = NULL;
  sbilf->f_fft        = NULL;
  sbilf->plan_forward = NULL;
}

static void
_ncm_sbessel_integrator_fftl_dispose (GObject *object)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (object);

  if (sbilf->plan_forward != NULL)
  {
    fftw_destroy_plan (sbilf->plan_forward);
    sbilf->plan_forward = NULL;
  }

  if (sbilf->f_samp != NULL)
  {
    fftw_free (sbilf->f_samp);
    sbilf->f_samp = NULL;
  }

  if (sbilf->f_fft != NULL)
  {
    fftw_free (sbilf->f_fft);
    sbilf->f_fft = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_fftl_parent_class)->dispose (object);
}

static void
_ncm_sbessel_integrator_fftl_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_fftl_parent_class)->finalize (object);
}

static void _ncm_sbessel_integrator_fftl_prepare (NcmSBesselIntegrator *sbi);
static gdouble _ncm_sbessel_integrator_fftl_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gint ell, gpointer user_data);

static void
_ncm_sbessel_integrator_fftl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_FFTL (object));

  switch (prop_id)
  {
    case PROP_OVERSAMPLE:
      sbilf->oversample = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_integrator_fftl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_FFTL (object));

  switch (prop_id)
  {
    case PROP_OVERSAMPLE:
      g_value_set_double (value, sbilf->oversample);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_sbessel_integrator_fftl_class_init (NcmSBesselIntegratorFFTLClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmSBesselIntegratorClass *parent_class = NCM_SBESSEL_INTEGRATOR_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_integrator_fftl_set_property;
  object_class->get_property = &_ncm_sbessel_integrator_fftl_get_property;
  object_class->dispose      = &_ncm_sbessel_integrator_fftl_dispose;
  object_class->finalize     = &_ncm_sbessel_integrator_fftl_finalize;

  /**
   * NcmSBesselIntegratorFFTL:oversample:
   *
   * Oversampling factor for FFT grid (Ny ~ oversample * (b-a)/pi).
   */
  g_object_class_install_property (object_class,
                                   PROP_OVERSAMPLE,
                                   g_param_spec_double ("oversample",
                                                        NULL,
                                                        "Oversampling factor",
                                                        1.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->prepare       = &_ncm_sbessel_integrator_fftl_prepare;
  parent_class->integrate_ell = &_ncm_sbessel_integrator_fftl_integrate_ell;
}

static void
_ncm_sbessel_integrator_fftl_prepare (NcmSBesselIntegrator *sbi)
{
  /* NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (sbi); */

  /* Preparation logic */
}

static gdouble
_ncm_sbessel_integrator_fftl_integrate_ell (NcmSBesselIntegrator *sbi,
                                            NcmSBesselIntegratorF F,
                                            gdouble a, gdouble b,
                                            gint ell,
                                            gpointer user_data)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (sbi);
  const gdouble L_phys            = b - a;
  const gdouble L_phys_2pi        = L_phys / (2.0 * M_PI);
  gint alpha_epsilon              = 5; /* controls accuracy epsilon = 2^(-alpha_epsilon)/kappa */
  const gint alpha_J              = 5; /* controls accuracy the mu term */
  const gint J                    = 1 << alpha_J;
  gdouble L;
  gint alpha_L;
  guint i, n, Ny;

  frexp (0.5 * L_phys_2pi, &alpha_L); /* This is kappa/2 */
  alpha_epsilon = GSL_MAX (alpha_epsilon, alpha_L + 1);
  L             = 2.0 * M_PI * J;
  Ny            = 1 << (alpha_epsilon + alpha_J - alpha_L);

  /* Reallocate arrays if needed */
  if (Ny != sbilf->Ny)
  {
    if (sbilf->plan_forward != NULL)
      fftw_destroy_plan (sbilf->plan_forward);

    if (sbilf->f_samp != NULL)
      fftw_free (sbilf->f_samp);

    if (sbilf->f_fft != NULL)
      fftw_free (sbilf->f_fft);

    sbilf->Ny     = Ny;
    sbilf->f_samp = fftw_malloc (sizeof (gdouble) * Ny);
    sbilf->f_fft  = fftw_malloc (sizeof (fftw_complex) * (Ny / 2 + 1));

    sbilf->plan_forward = fftw_plan_dft_r2c_1d (Ny, sbilf->f_samp, sbilf->f_fft, FFTW_ESTIMATE);
  }

  /* Step 3: Sample f(x) on grid */
  const gdouble dx = L / Ny;

  memset (sbilf->f_samp, 0, sizeof (gdouble) * Ny);

  for (i = 0; i < Ny; i++)
  {
    const gdouble x = a + (i + 0.5) * dx;

    if (x > b)
      sbilf->f_samp[i] = 0.0;
    else
      sbilf->f_samp[i] = F (user_data, x);
  }

  /* Step 5: Compute FFT */
  fftw_execute (sbilf->plan_forward);

  /* Steps 6-9: Build mu-grid, compute Legendre polynomials, and sum */
  gdouble sum            = 0.0;
  guint m                = ell / 2;
  gboolean ell_is_even   = (ell % 2 == 0);
  gboolean m_is_even     = (m % 2 == 0);
  complex double phase   = cexp (I * (a + 0.5 * dx) * (-2.0 * M_PI / L));
  complex double phase_n = phase;

  if (ell_is_even)
  {
    sum = gsl_sf_legendre_Pl (ell, 0.0) * creal (sbilf->f_fft[0]);

    for (n = 1; n <= J; n++)
    {
      const gdouble mu       = -2.0 * M_PI * n / L; /* positive  frequency */
      const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
      const gdouble fft_real = creal (phase_n * sbilf->f_fft[n]);
      const gdouble w        = (n == J) ? 1.0 : ((n % 2 == 1) ? 4.0 : 2.0); /* Simpson's rule weight */

      sum     += w * P_ell * fft_real;
      phase_n *= phase;
    }
  }
  else
  {
    for (n = 1; n <= J; n++)
    {
      double mu              = -2.0 * M_PI * n / L; /* positive  frequency */
      gdouble P_ell          = gsl_sf_legendre_Pl (ell, mu);
      const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[n]);
      const gdouble w        = (n == J) ? 1.0 : ((n % 2 == 1) ? 4.0 : 2.0); /* Simpson's rule weight */

      sum     += w * P_ell * fft_imag;
      phase_n *= phase;
    }
  }

  sum = sum / 3.0;
  sum = (m_is_even ? sum : -sum) * (2.0 * M_PI / Ny);

  return sum;
}

/**
 * ncm_sbessel_integrator_fftl_new:
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 *
 * Creates a new #NcmSBesselIntegratorFFTL.
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorFFTL
 */
NcmSBesselIntegratorFFTL *
ncm_sbessel_integrator_fftl_new (guint lmin, guint lmax)
{
  NcmSBesselIntegratorFFTL *sbilf = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_FFTL,
                                                  "lmin", lmin,
                                                  "lmax", lmax,
                                                  NULL);

  return sbilf;
}

/**
 * ncm_sbessel_integrator_fftl_ref:
 * @sbilf: a #NcmSBesselIntegratorFFTL
 *
 * Increases the reference count of @sbilf by one.
 *
 * Returns: (transfer full): @sbilf
 */
NcmSBesselIntegratorFFTL *
ncm_sbessel_integrator_fftl_ref (NcmSBesselIntegratorFFTL *sbilf)
{
  return g_object_ref (sbilf);
}

/**
 * ncm_sbessel_integrator_fftl_free:
 * @sbilf: a #NcmSBesselIntegratorFFTL
 *
 * Decreases the reference count of @sbilf by one.
 *
 */
void
ncm_sbessel_integrator_fftl_free (NcmSBesselIntegratorFFTL *sbilf)
{
  g_object_unref (sbilf);
}

/**
 * ncm_sbessel_integrator_fftl_clear:
 * @sbilf: a #NcmSBesselIntegratorFFTL
 *
 * If @sbilf is different from NULL, decreases the reference count of
 * @sbilf by one and sets @sbilf to NULL.
 *
 */
void
ncm_sbessel_integrator_fftl_clear (NcmSBesselIntegratorFFTL **sbilf)
{
  g_clear_object (sbilf);
}

