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
  guint lmax_alloc;
  gdouble *f_samp;
  fftw_complex *f_fft;
  fftw_plan plan_forward;
  gdouble *jl_arr;
  gdouble *Pl_arr;
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
  sbilf->lmax_alloc   = 0;
  sbilf->f_samp       = NULL;
  sbilf->f_fft        = NULL;
  sbilf->plan_forward = NULL;
  sbilf->jl_arr       = NULL;
  sbilf->Pl_arr       = NULL;
}

static void
_ncm_sbessel_integrator_fftl_dispose (GObject *object)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (object);

  g_clear_pointer (&sbilf->plan_forward, fftw_destroy_plan);
  g_clear_pointer (&sbilf->f_samp, fftw_free);
  g_clear_pointer (&sbilf->f_fft, fftw_free);
  g_clear_pointer (&sbilf->jl_arr, g_free);
  g_clear_pointer (&sbilf->Pl_arr, g_free);

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
static void _ncm_sbessel_integrator_fftl_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, NcmVector *result, gpointer user_data);

typedef struct _NcmSBesselIntegratorFFTLParams
{
  gdouble L;
  gdouble dx;
  guint Ny;
  guint J;
} NcmSBesselIntegratorFFTLParams;

static gdouble
_ncm_sbessel_integrator_fftl_simpson_direct (NcmSBesselIntegratorF F, const gdouble a, const gdouble b, const gint ell, gpointer user_data)
{
  const guint N    = 2 * GSL_MAX (256, 4 * ell);
  const gdouble dx = (b - a) / N;
  gdouble result   = 0.0;
  guint i;

  /* First term */
  result = F (user_data, a) * gsl_sf_bessel_jl (ell, a);

  /* Interior terms with alternating weights 4 and 2 */
  for (i = 1; i < N; i++)
  {
    const gdouble x      = a + i * dx;
    const gdouble weight = (i % 2 == 1) ? 4.0 : 2.0;

    result += weight * F (user_data, x) * gsl_sf_bessel_jl (ell, x);
  }

  /* Last term */
  result += F (user_data, b) * gsl_sf_bessel_jl (ell, b);

  return result * dx / 3.0;
}

static void
_ncm_sbessel_integrator_fftl_compute_params (const gdouble L_phys, NcmSBesselIntegratorFFTLParams *params)
{
  const gdouble L_phys_2pi = L_phys / (2.0 * M_PI);
  gint alpha_epsilon       = 10; /* controls accuracy epsilon = 2^(-alpha_epsilon)/kappa */
  gint alpha_J             = 10; /* controls accuracy the mu term */
  const gint min_alpha     = 13;
  gint alpha_L;

  frexp (0.5 * L_phys_2pi, &alpha_L); /* This is kappa/2 */
  alpha_epsilon = GSL_MAX (alpha_epsilon, alpha_L + 1);

  /* FFT of at least Ny = 2^min_alpha */
  if (alpha_epsilon + alpha_J - alpha_L < min_alpha)
  {
    gint delta = min_alpha - (alpha_epsilon + alpha_J - alpha_L);

    alpha_J += delta / 2;
    alpha_L += delta / 2 + (delta % 2);
  }

  params->J  = 1 << alpha_J;
  params->L  = 2.0 * M_PI * params->J;
  params->Ny = 1 << (alpha_epsilon + alpha_J - alpha_L);
  params->dx = params->L / params->Ny;
/*#define DEBUG */
#ifdef DEBUG
  g_print ("[FFTL] Computed params: L = %.3e, J = %u, Ny = %u, dx = %.3e\n",
           params->L, params->J, params->Ny, params->dx);
#endif /* DEBUG */
}

static void
_ncm_sbessel_integrator_fftl_sample (NcmSBesselIntegratorFFTL *sbilf, NcmSBesselIntegratorF F, const gdouble a, const gdouble b, const NcmSBesselIntegratorFFTLParams *params, gpointer user_data)
{
  guint i;

  memset (sbilf->f_samp, 0, sizeof (gdouble) * params->Ny);

  for (i = 0; i < params->Ny; i++)
  {
    const gdouble x = a + (i + 0.5) * params->dx;

    if (x > b)
      sbilf->f_samp[i] = 0.0;
    else
      sbilf->f_samp[i] = F (user_data, x);
  }
}

static gdouble
_ncm_sbessel_integrator_fftl_sum_even_ell (NcmSBesselIntegratorFFTL *sbilf, const gint ell, const gdouble a, const NcmSBesselIntegratorFFTLParams *params)
{
  gdouble sum                = 0.0;
  const complex double phase = cexp (I * (a + 0.5 * params->dx) * (-2.0 * M_PI / params->L));
  complex double phase_n     = phase;
  const gdouble mu_over_n    = -2.0 * M_PI / params->L;
  guint n;

  /* n = 0 term */
  sum = gsl_sf_legendre_Pl (ell, 0.0) * creal (sbilf->f_fft[0]);

  /* Loop over pairs of indices with Simpson weights (4,2) */
  /* Stop before the last two terms which need special weights (4,1) */
  for (n = 1; n + 1 < params->J; n += 2)
  {
    /* n = odd (weight 4.0) */
    {
      const gdouble mu       = mu_over_n * n;
      const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
      const gdouble fft_real = creal (phase_n * sbilf->f_fft[n]);

      sum     += 4.0 * P_ell * fft_real;
      phase_n *= phase;
    }

    /* n+1 = even (weight 2.0) */
    {
      const gdouble mu       = mu_over_n * (n + 1);
      const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
      const gdouble fft_real = creal (phase_n * sbilf->f_fft[n + 1]);

      sum     += 2.0 * P_ell * fft_real;
      phase_n *= phase;
    }
  }

  /* Second to last term n = J-1 (weight 4.0) */
  {
    const gdouble mu       = mu_over_n * (params->J - 1);
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_real = creal (phase_n * sbilf->f_fft[params->J - 1]);

    sum     += 4.0 * P_ell * fft_real;
    phase_n *= phase;
  }

  /* Last term n = J (weight 1.0) */
  {
    const gdouble mu       = mu_over_n * params->J;
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_real = creal (phase_n * sbilf->f_fft[params->J]);

    sum += 1.0 * P_ell * fft_real;
  }

  return sum;
}

static gdouble
_ncm_sbessel_integrator_fftl_sum_odd_ell (NcmSBesselIntegratorFFTL *sbilf, const gint ell, const gdouble a, const NcmSBesselIntegratorFFTLParams *params)
{
  gdouble sum                = 0.0;
  const complex double phase = cexp (I * (a + 0.5 * params->dx) * (-2.0 * M_PI / params->L));
  complex double phase_n     = phase;
  const gdouble mu_over_n    = -2.0 * M_PI / params->L;
  guint n;

  /* Loop over pairs of indices with Simpson weights (4,2) */
  /* Stop before the last two terms which need special weights (4,1) */
  for (n = 1; n + 1 < params->J; n += 2)
  {
    /* n = odd (weight 4.0) */
    {
      const gdouble mu       = mu_over_n * n;
      const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
      const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[n]);

      sum     += 4.0 * P_ell * fft_imag;
      phase_n *= phase;
    }

    /* n+1 = even (weight 2.0) */
    {
      const gdouble mu       = mu_over_n * (n + 1);
      const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
      const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[n + 1]);

      sum     += 2.0 * P_ell * fft_imag;
      phase_n *= phase;
    }
  }

  /* Second to last term n = J-1 (weight 4.0) */
  {
    const gdouble mu       = mu_over_n * (params->J - 1);
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[params->J - 1]);

    sum     += 4.0 * P_ell * fft_imag;
    phase_n *= phase;
  }

  /* Last term n = J (weight 1.0) */
  {
    const gdouble mu       = mu_over_n * params->J;
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[params->J]);

    sum += 1.0 * P_ell * fft_imag;
  }

  return sum;
}

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
  parent_class->integrate     = &_ncm_sbessel_integrator_fftl_integrate;
}

static void
_ncm_sbessel_integrator_fftl_prepare (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (sbi);
  const guint lmax                = ncm_sbessel_integrator_get_lmax (sbi);

  /* Reallocate arrays if lmax changed */
  if (lmax != sbilf->lmax_alloc)
  {
    g_clear_pointer (&sbilf->jl_arr, g_free);
    g_clear_pointer (&sbilf->Pl_arr, g_free);

    sbilf->lmax_alloc = lmax;
    sbilf->jl_arr     = g_new (gdouble, lmax + 1);
    sbilf->Pl_arr     = g_new (gdouble, lmax + 1);
  }
}

static gdouble
_ncm_sbessel_integrator_fftl_integrate_ell (NcmSBesselIntegrator *sbi,
                                            NcmSBesselIntegratorF F,
                                            gdouble a, gdouble b,
                                            gint ell,
                                            gpointer user_data)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (sbi);
  NcmSBesselIntegratorFFTLParams params;
  const gdouble L_phys = b - a;
  gdouble sum;
  guint m;
  gboolean ell_is_even;
  gboolean m_is_even;

  /* For small intervals (b < ell), use direct Simpson's rule */
  if (b < ell)
  {
    printf ("# DIRECT SIMPSON METHOD FOR ELL %d\n", ell);

    return _ncm_sbessel_integrator_fftl_simpson_direct (F, a, b, ell, user_data);
  }

  /* Compute integration parameters */
  _ncm_sbessel_integrator_fftl_compute_params (L_phys, &params);

  /* Reallocate arrays if needed */
  if (params.Ny != sbilf->Ny)
  {
    if (sbilf->plan_forward != NULL)
      fftw_destroy_plan (sbilf->plan_forward);

    if (sbilf->f_samp != NULL)
      fftw_free (sbilf->f_samp);

    if (sbilf->f_fft != NULL)
      fftw_free (sbilf->f_fft);

    sbilf->Ny     = params.Ny;
    sbilf->f_samp = fftw_malloc (sizeof (gdouble) * params.Ny);
    sbilf->f_fft  = fftw_malloc (sizeof (fftw_complex) * (params.Ny / 2 + 1));

    sbilf->plan_forward = fftw_plan_dft_r2c_1d (params.Ny, sbilf->f_samp, sbilf->f_fft, FFTW_ESTIMATE);
  }

  /* Sample function on grid */
  _ncm_sbessel_integrator_fftl_sample (sbilf, F, a, b, &params, user_data);

  /* Compute FFT */
  fftw_execute (sbilf->plan_forward);

  /* Compute Legendre sum (different for even/odd ell) */
  ell_is_even = (ell % 2 == 0);

  if (ell_is_even)
    sum = _ncm_sbessel_integrator_fftl_sum_even_ell (sbilf, ell, a, &params);
  else
    sum = _ncm_sbessel_integrator_fftl_sum_odd_ell (sbilf, ell, a, &params);

  /* Apply Simpson's rule factor and final scaling */
  m         = ell / 2;
  m_is_even = (m % 2 == 0);
  sum       = sum / 3.0;
  sum       = (m_is_even ? sum : -sum) * (2.0 * M_PI / params.Ny);

  return sum;
}

static void
_ncm_sbessel_integrator_fftl_integrate_direct (NcmSBesselIntegratorFFTL *sbilf, const guint lmin, const guint ell_direct_min, const guint ell_direct_max,
                                               NcmSBesselIntegratorF F, const gdouble a, const gdouble b,
                                               NcmVector *result, gpointer user_data)
{
  const guint N    = 2 * GSL_MAX (256, 4 * ell_direct_max);
  const gdouble dx = (b - a) / N;
  guint i, ell;

  /* Initialize direct results to zero */
  for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
    ncm_vector_set (result, ell - lmin, 0.0);

  /* First term */
  {
    const gdouble fa = F (user_data, a);

    gsl_sf_bessel_jl_steed_array (ell_direct_max, a, sbilf->jl_arr);

    for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
      ncm_vector_set (result, ell - lmin, fa * sbilf->jl_arr[ell]);
  }

  /* Interior terms with alternating weights 4 and 2 */
  for (i = 1; i < N; i++)
  {
    const gdouble x      = a + i * dx;
    const gdouble weight = (i % 2 == 1) ? 4.0 : 2.0;
    const gdouble fx     = F (user_data, x);

    gsl_sf_bessel_jl_steed_array (ell_direct_max, x, sbilf->jl_arr);

    for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
      ncm_vector_addto (result, ell - lmin, weight * fx * sbilf->jl_arr[ell]);
  }

  /* Last term */
  {
    const gdouble fb = F (user_data, b);

    gsl_sf_bessel_jl_steed_array (ell_direct_max, b, sbilf->jl_arr);

    for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
      ncm_vector_addto (result, ell - lmin, fb * sbilf->jl_arr[ell]);
  }

  /* Apply Simpson's rule factor */
  for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
    ncm_vector_set (result, ell - lmin, ncm_vector_get (result, ell - lmin) * dx / 3.0);
}

static void
_ncm_sbessel_integrator_fftl_integrate_fft (NcmSBesselIntegratorFFTL *sbilf, const guint lmin, const guint ell_fft_min, const guint ell_fft_max,
                                            NcmSBesselIntegratorF F, const gdouble a, const gdouble b, const gdouble L_phys,
                                            NcmVector *result, gpointer user_data)
{
  NcmSBesselIntegratorFFTLParams params;
  guint ell, n;

  /* Compute integration parameters */
  _ncm_sbessel_integrator_fftl_compute_params (L_phys, &params);

  /* Reallocate arrays if needed */
  if (params.Ny != sbilf->Ny)
  {
    if (sbilf->plan_forward != NULL)
      fftw_destroy_plan (sbilf->plan_forward);

    if (sbilf->f_samp != NULL)
      fftw_free (sbilf->f_samp);

    if (sbilf->f_fft != NULL)
      fftw_free (sbilf->f_fft);

    sbilf->Ny     = params.Ny;
    sbilf->f_samp = fftw_malloc (sizeof (gdouble) * params.Ny);
    sbilf->f_fft  = fftw_malloc (sizeof (fftw_complex) * (params.Ny / 2 + 1));

    sbilf->plan_forward = fftw_plan_dft_r2c_1d (params.Ny, sbilf->f_samp, sbilf->f_fft, FFTW_ESTIMATE);
  }

  /* Sample function on grid */
  _ncm_sbessel_integrator_fftl_sample (sbilf, F, a, b, &params, user_data);

  /* Compute FFT */
  fftw_execute (sbilf->plan_forward);

  /* Initialize FFT results to zero */
  for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    ncm_vector_set (result, ell - lmin, 0.0);

  /* Compute phase factor */
  const complex double phase = cexp (I * (a + 0.5 * params.dx) * (-2.0 * M_PI / params.L));
  const gdouble mu_over_n    = -2.0 * M_PI / params.L;

  /* n = 0 term */
  {
    gsl_sf_legendre_Pl_array (ell_fft_max, 0.0, sbilf->Pl_arr);

    const gdouble fft_real = creal (sbilf->f_fft[0]);

    for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    {
      if (ell % 2 == 0) /* even ell uses real part */
        ncm_vector_set (result, ell - lmin, sbilf->Pl_arr[ell] * fft_real);
    }
  }

  /* Loop over pairs with Simpson weights */
  complex double phase_n = phase;

  for (n = 1; n + 1 < params.J; n += 2)
  {
    /* n = odd (weight 4.0) */
    {
      const gdouble mu = mu_over_n * n;

      gsl_sf_legendre_Pl_array (ell_fft_max, mu, sbilf->Pl_arr);

      for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
      {
        const gdouble P_ell = sbilf->Pl_arr[ell];

        if (ell % 2 == 0)
          ncm_vector_addto (result, ell - lmin, 4.0 * P_ell * creal (phase_n * sbilf->f_fft[n]));
        else
          ncm_vector_addto (result, ell - lmin, 4.0 * P_ell * cimag (phase_n * sbilf->f_fft[n]));
      }

      phase_n *= phase;
    }

    /* n+1 = even (weight 2.0) */
    {
      const gdouble mu = mu_over_n * (n + 1);

      gsl_sf_legendre_Pl_array (ell_fft_max, mu, sbilf->Pl_arr);

      for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
      {
        const gdouble P_ell = sbilf->Pl_arr[ell];

        if (ell % 2 == 0)
          ncm_vector_addto (result, ell - lmin, 2.0 * P_ell * creal (phase_n * sbilf->f_fft[n + 1]));
        else
          ncm_vector_addto (result, ell - lmin, 2.0 * P_ell * cimag (phase_n * sbilf->f_fft[n + 1]));
      }

      phase_n *= phase;
    }
  }

  /* Second to last term n = J-1 (weight 4.0) */
  {
    const gdouble mu = mu_over_n * (params.J - 1);

    gsl_sf_legendre_Pl_array (ell_fft_max, mu, sbilf->Pl_arr);

    for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    {
      const gdouble P_ell = sbilf->Pl_arr[ell];

      if (ell % 2 == 0)
        ncm_vector_addto (result, ell - lmin, 4.0 * P_ell * creal (phase_n * sbilf->f_fft[params.J - 1]));
      else
        ncm_vector_addto (result, ell - lmin, 4.0 * P_ell * cimag (phase_n * sbilf->f_fft[params.J - 1]));
    }

    phase_n *= phase;
  }

  /* Last term n = J (weight 1.0) */
  {
    const gdouble mu = mu_over_n * params.J;

    gsl_sf_legendre_Pl_array (ell_fft_max, mu, sbilf->Pl_arr);

    for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    {
      const gdouble P_ell = sbilf->Pl_arr[ell];

      if (ell % 2 == 0)
        ncm_vector_addto (result, ell - lmin, 1.0 * P_ell * creal (phase_n * sbilf->f_fft[params.J]));
      else
        ncm_vector_addto (result, ell - lmin, 1.0 * P_ell * cimag (phase_n * sbilf->f_fft[params.J]));
    }
  }

  /* Apply Simpson's rule factor and final scaling */
  for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
  {
    const guint m         = ell / 2;
    const gboolean m_even = (m % 2 == 0);
    const gdouble val     = ncm_vector_get (result, ell - lmin);

    ncm_vector_set (result, ell - lmin, (m_even ? val : -val) * (2.0 * M_PI / params.Ny) / 3.0);
  }
}

static void
_ncm_sbessel_integrator_fftl_integrate (NcmSBesselIntegrator *sbi,
                                        NcmSBesselIntegratorF F,
                                        gdouble a, gdouble b,
                                        NcmVector *result,
                                        gpointer user_data)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (sbi);
  const guint lmin                = ncm_sbessel_integrator_get_lmin (sbi);
  const guint lmax                = ncm_sbessel_integrator_get_lmax (sbi);
  const gdouble L_phys            = b - a;
  const guint n_ell               = lmax - lmin + 1;

  g_assert_cmpuint (ncm_vector_len (result), ==, n_ell);

  /* Find the threshold ell where b < ell (use direct Simpson's for ell > b) */
  const guint ell_b = (guint) ceil (b);

  /* Compute ells that need direct Simpson's rule (where b < ell) */
  if (ell_b <= lmax)
  {
    const guint ell_direct_min = GSL_MAX (ell_b, lmin);
    const guint ell_direct_max = lmax;

    if (ell_direct_min <= ell_direct_max)
      _ncm_sbessel_integrator_fftl_integrate_direct (sbilf, lmin, ell_direct_min, ell_direct_max, F, a, b, result, user_data);
  }

  /* Compute ells that need FFT-Legendre (where b >= ell) */
  if (ell_b > lmin)
  {
    const guint ell_fft_min = lmin;
    const guint ell_fft_max = GSL_MIN (ell_b - 1, lmax);

    if (ell_fft_min <= ell_fft_max)
      _ncm_sbessel_integrator_fftl_integrate_fft (sbilf, lmin, ell_fft_min, ell_fft_max, F, a, b, L_phys, result, user_data);
  }
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

