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
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_c.h"
#include "ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_roots.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NCM_SBESSEL_INTEGRATOR_FFTL_MAX_ALPHA_J 20

typedef struct _NcmSBesselIntegratorFFTLCache
{
  guint lmax;      /* lmax for which this cache was computed */
  guint J;         /* J = 1 << alpha_J */
  gdouble **Pl_mu; /* Array of size (J+1), each entry is array of size (lmax+1) with P_ell(mu_n) */
} NcmSBesselIntegratorFFTLCache;

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
  guint N_direct;
  guint alpha_J;
  NcmSBesselIntegratorFFTLCache *cache[NCM_SBESSEL_INTEGRATOR_FFTL_MAX_ALPHA_J]; /* Cache indexed by alpha_J */
  NcmSFSBesselArray *sba;                                                        /* Spherical Bessel array with cutoff */
  NcmSpline *wrap_f;
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
  guint i;

  sbilf->oversample   = 2.0;
  sbilf->Ny           = 0;
  sbilf->lmax_alloc   = 0;
  sbilf->f_samp       = NULL;
  sbilf->f_fft        = NULL;
  sbilf->plan_forward = NULL;
  sbilf->jl_arr       = NULL;
  sbilf->Pl_arr       = NULL;
  sbilf->N_direct     = 0;
  sbilf->alpha_J      = 0;
  sbilf->sba          = ncm_sf_sbessel_array_new ();

  for (i = 0; i < NCM_SBESSEL_INTEGRATOR_FFTL_MAX_ALPHA_J; i++)
    sbilf->cache[i] = NULL;

  sbilf->wrap_f = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  {
    NcmVector *x_vals = ncm_vector_new (10);
    NcmVector *y_vals = ncm_vector_new (10);

    ncm_spline_set (sbilf->wrap_f, x_vals, y_vals, FALSE);
  }
}

static void
_ncm_sbessel_integrator_fftl_dispose (GObject *object)
{
  NcmSBesselIntegratorFFTL *sbilf = NCM_SBESSEL_INTEGRATOR_FFTL (object);
  guint i;

  g_clear_pointer (&sbilf->plan_forward, fftw_destroy_plan);
  g_clear_pointer (&sbilf->f_samp, fftw_free);
  g_clear_pointer (&sbilf->f_fft, fftw_free);
  g_clear_pointer (&sbilf->jl_arr, g_free);
  g_clear_pointer (&sbilf->Pl_arr, g_free);

  ncm_sf_sbessel_array_clear (&sbilf->sba);

  for (i = 0; i < NCM_SBESSEL_INTEGRATOR_FFTL_MAX_ALPHA_J; i++)
  {
    if (sbilf->cache[i] != NULL)
    {
      guint n;

      for (n = 0; n <= sbilf->cache[i]->J; n++)
        g_clear_pointer (&sbilf->cache[i]->Pl_mu[n], g_free);

      g_clear_pointer (&sbilf->cache[i]->Pl_mu, g_free);
      g_clear_pointer (&sbilf->cache[i], g_free);
    }
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

typedef struct _NcmSBesselIntegratorFFTLIntegrationData
{
  /* Function information */
  NcmSBesselIntegratorF F;
  gpointer user_data;
  gdouble a;
  gdouble b;
  gdouble L_phys;      /* b - a */
  gdouble x_peak;      /* Location of peak */
  gdouble peak_val;    /* |f(x_peak)| */
  gdouble support_eps; /* Epsilon threshold for support (default 1e-3) */
  gdouble x_left;      /* Left edge where |f| drops to epsilon*peak_val */
  gdouble x_right;     /* Right edge where |f| drops to epsilon*peak_val */
  gboolean peak_found;
  gboolean support_found;
  /* FFT parameters */
  gdouble L;
  gdouble dx;
  guint Ny;
  guint J;
  guint alpha_J;
} NcmSBesselIntegratorFFTLIntegrationData;

static void _ncm_sbessel_integrator_fftl_prepare (NcmSBesselIntegrator *sbi);
static guint _ncm_sbessel_integrator_fftl_get_ell_threshold (NcmSBesselIntegratorFFTL *sbilf, NcmSBesselIntegratorFFTLIntegrationData *data);
static void _ncm_sbessel_integrator_fftl_find_peak (NcmSBesselIntegratorFFTLIntegrationData *data);
static gdouble _ncm_sbessel_integrator_fftl_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gint ell, gpointer user_data);
static void _ncm_sbessel_integrator_fftl_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, NcmVector *result, gpointer user_data);

/* Wrapper for GSL minimization - we minimize -|f(x)| to find max|f(x)| */
static gdouble
_ncm_sbessel_integrator_fftl_func_neg_abs (gdouble x, gpointer params)
{
  NcmSBesselIntegratorFFTLIntegrationData *data = (NcmSBesselIntegratorFFTLIntegrationData *) params;
  const gdouble fx                              = data->F (data->user_data, x);

  return -fabs (fx);
}

/* Helper to compute |f(x)| - threshold for root finding */
static gdouble
_ncm_sbessel_integrator_fftl_func_minus_threshold (gdouble x, gpointer params)
{
  NcmSBesselIntegratorFFTLIntegrationData *data = (NcmSBesselIntegratorFFTLIntegrationData *) params;
  const gdouble fx                              = data->F (data->user_data, x);
  const gdouble threshold                       = data->support_eps * data->peak_val;

  return fabs (fx) - threshold;
}

/* Find the peak (maximum of |f|) in [a, b] using GSL's Brent method */
static void
_ncm_sbessel_integrator_fftl_find_peak (NcmSBesselIntegratorFFTLIntegrationData *data)
{
  gsl_function F;
  gsl_min_fminimizer *minimizer;
  const gint max_iter = 100;
  gint iter           = 0;
  gint status;
  gdouble x_lower, x_upper, x_minimum;

  if (data->peak_found)
    return;

  F.function = &_ncm_sbessel_integrator_fftl_func_neg_abs;
  F.params   = data;

  minimizer = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);

  /* Initial guess: midpoint */
  x_minimum = 0.5 * (data->a + data->b);
  x_lower   = data->a;
  x_upper   = data->b;

  /* Set up the minimizer */
  gsl_min_fminimizer_set (minimizer, &F, x_minimum, x_lower, x_upper);

  /* Iterate to find minimum */
  do {
    iter++;
    status = gsl_min_fminimizer_iterate (minimizer);

    x_minimum = gsl_min_fminimizer_x_minimum (minimizer);
    x_lower   = gsl_min_fminimizer_x_lower (minimizer);
    x_upper   = gsl_min_fminimizer_x_upper (minimizer);

    status = gsl_min_test_interval (x_lower, x_upper, 0.0, 1.0e-3);
  } while (status == GSL_CONTINUE && iter < max_iter);

  data->x_peak     = x_minimum;
  data->peak_val   = fabs (data->F (data->user_data, x_minimum));
  data->peak_found = TRUE;

  gsl_min_fminimizer_free (minimizer);
}

/* Find support edges where |f| drops to epsilon * peak_val */
static void
_ncm_sbessel_integrator_fftl_find_support (NcmSBesselIntegratorFFTLIntegrationData *data)
{
  gsl_function F;
  gsl_root_fsolver *solver;
  const gint max_iter = 100;
  gint status;

  if (data->support_found || !data->peak_found)
    return;

  /* If function value at boundaries is already below threshold, use boundaries */
  const gdouble threshold = data->support_eps * data->peak_val;
  const gdouble f_a       = fabs (data->F (data->user_data, data->a));
  const gdouble f_b       = fabs (data->F (data->user_data, data->b));

  F.function = &_ncm_sbessel_integrator_fftl_func_minus_threshold;
  F.params   = data;

  /* Find left edge: search between a and x_peak */
  if (f_a < threshold)
  {
    data->x_left = data->a;
  }
  else
  {
    gdouble x_lo = data->a;
    gdouble x_hi = data->x_peak;
    gint iter    = 0;

    solver = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
    gsl_root_fsolver_set (solver, &F, x_lo, x_hi);

    do {
      iter++;
      status       = gsl_root_fsolver_iterate (solver);
      data->x_left = gsl_root_fsolver_root (solver);
      x_lo         = gsl_root_fsolver_x_lower (solver);
      x_hi         = gsl_root_fsolver_x_upper (solver);
      status       = gsl_root_test_interval (x_lo, x_hi, 0.0, 1.0e-3);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (solver);
  }

  /* Find right edge: search between x_peak and b */
  if (f_b < threshold)
  {
    data->x_right = data->b;
  }
  else
  {
    gdouble x_lo = data->x_peak;
    gdouble x_hi = data->b;
    gint iter    = 0;

    solver = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
    gsl_root_fsolver_set (solver, &F, x_lo, x_hi);

    do {
      iter++;
      status        = gsl_root_fsolver_iterate (solver);
      data->x_right = gsl_root_fsolver_root (solver);
      x_lo          = gsl_root_fsolver_x_lower (solver);
      x_hi          = gsl_root_fsolver_x_upper (solver);
      status        = gsl_root_test_interval (x_lo, x_hi, 0.0, 1.0e-3);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (solver);
  }

  data->support_found = TRUE;
}

static gdouble
_ncm_sbessel_integrator_fftl_simpson_direct (NcmSBesselIntegratorFFTL *sbilf, NcmSBesselIntegratorF F, const gdouble a, const gdouble b, const guint ell, gpointer user_data)
{
  const guint N    = sbilf->N_direct;
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

static NcmSBesselIntegratorFFTLCache *
_ncm_sbessel_integrator_fftl_get_cache (NcmSBesselIntegratorFFTL *sbilf, guint alpha_J, guint lmax)
{
  g_assert_cmpuint (alpha_J, <, NCM_SBESSEL_INTEGRATOR_FFTL_MAX_ALPHA_J);

  /* Check if cache exists and is valid */
  if ((sbilf->cache[alpha_J] != NULL) && (sbilf->cache[alpha_J]->lmax >= lmax))
    return sbilf->cache[alpha_J];

  /* Free old cache if it exists */
  if (sbilf->cache[alpha_J] != NULL)
  {
    guint n;

    for (n = 0; n <= sbilf->cache[alpha_J]->J; n++)
      g_clear_pointer (&sbilf->cache[alpha_J]->Pl_mu[n], g_free);

    g_clear_pointer (&sbilf->cache[alpha_J]->Pl_mu, g_free);
    g_clear_pointer (&sbilf->cache[alpha_J], g_free);
  }

  /* Create new cache */
  {
    NcmSBesselIntegratorFFTLCache *cache = g_new (NcmSBesselIntegratorFFTLCache, 1);
    const guint J                        = 1 << alpha_J;
    guint n;

    cache->lmax  = lmax;
    cache->J     = J;
    cache->Pl_mu = g_new (gdouble *, J + 1);

    /* Precompute Legendre polynomials for all mu knots */
    for (n = 0; n <= J; n++)
    {
      const gdouble mu = -((gdouble) n) / J;

      cache->Pl_mu[n] = g_new (gdouble, lmax + 1);
      gsl_sf_legendre_Pl_array (lmax, mu, cache->Pl_mu[n]);
    }

    sbilf->cache[alpha_J] = cache;

    return cache;
  }
}

static void
_ncm_sbessel_integrator_fftl_compute_params (NcmSBesselIntegratorFFTLIntegrationData *data)
{
  const gint alpha_L_phys = ceil (log2 (data->L_phys / M_PI));
  gint alpha_epsilon      = 10; /* controls accuracy epsilon = 2^(-alpha_epsilon)/kappa */
  gint alpha_J            = 10; /* controls accuracy the mu term */
  gint alpha0;

  if (alpha_epsilon < alpha_L_phys + 1)
    alpha_epsilon = alpha_L_phys + 1;

  if (alpha_J < alpha_L_phys + 1)
    alpha_J = alpha_L_phys + 1;

  alpha0 = alpha_epsilon - (alpha_L_phys + 1);

  data->alpha_J = alpha_J;
  data->J       = 1 << alpha_J;
  data->L       = 2.0 * M_PI * data->J;
  data->Ny      = 1 << (alpha_J + 2 + alpha0);
  data->dx      = data->L / data->Ny;
/* #define DEBUG */
#ifdef DEBUG
  g_print ("[FFTL] Computed alphas: alpha_L_phys = %d, alpha_J = %d, alpha0 = %d, alpha_Ny = %d\n", alpha_L_phys, alpha_J, alpha0, alpha_J + 2 + alpha0);
  g_print ("[FFTL] Computed params: L = %.3e, L_phys = %.3e, J = %u, Ny = %u, dx = %.3e, L_phys / dx = %.3e\n",
           data->L, data->L_phys, data->J, data->Ny, data->dx, data->L_phys / data->dx);
#endif /* DEBUG */
}

static void
_ncm_sbessel_integrator_fftl_sample (NcmSBesselIntegratorFFTL *sbilf, const NcmSBesselIntegratorFFTLIntegrationData *data)
{
  guint i;
  NcmVector *x_vals = ncm_spline_peek_xv (sbilf->wrap_f);
  NcmVector *y_vals = ncm_spline_peek_yv (sbilf->wrap_f);

  for (i = 0; i < 5; i++)
  {
    const gdouble x = data->b - (4.0 - i) * data->dx;
    const gdouble y = data->F (data->user_data, x);

    /* printf ("[FFTL] Pre-sample  %u: x=% 22.15e, y=% 22.15e\n", i, x, y); */

    ncm_vector_fast_set (x_vals, i, x);
    ncm_vector_fast_set (y_vals, i, y);
  }

  for (i = 0; i < 5; i++)
  {
    const gdouble x   = data->a + i * data->dx;
    const gdouble x_L = x + data->L;
    const gdouble y   = data->F (data->user_data, x);

    /* printf ("[FFTL] Post-sample %u: x=% 22.15e, y=% 22.15e\n", i + 5, x_L, y); */

    ncm_vector_fast_set (x_vals, i + 5, x_L);
    ncm_vector_fast_set (y_vals, i + 5, y);
  }

  ncm_spline_prepare (sbilf->wrap_f);

  memset (sbilf->f_samp, 0, sizeof (gdouble) * data->Ny);

  for (i = 0; i < data->Ny; i++)
  {
    const gdouble x = data->a + (i + 0.5) * data->dx;

    if (x > data->b)
      sbilf->f_samp[i] = ncm_spline_eval (sbilf->wrap_f, x);
    else
      sbilf->f_samp[i] = data->F (data->user_data, x);
  }
}

static gdouble
_ncm_sbessel_integrator_fftl_sum_even_ell (NcmSBesselIntegratorFFTL *sbilf, const guint ell, const NcmSBesselIntegratorFFTLIntegrationData *data)
{
  gdouble sum                = 0.0;
  const complex double phase = cexp (I * (data->a + 0.5 * data->dx) * (-2.0 * M_PI / data->L));
  complex double phase_n     = phase;
  const gdouble mu_over_n    = -2.0 * M_PI / data->L;
  guint n;

  /* n = 0 term */
  sum = gsl_sf_legendre_Pl (ell, 0.0) * creal (sbilf->f_fft[0]);

  /* Loop over pairs of indices with Simpson weights (4,2) */
  /* Stop before the last two terms which need special weights (4,1) */
  for (n = 1; n + 1 < data->J; n += 2)
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
    const gdouble mu       = mu_over_n * (data->J - 1);
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_real = creal (phase_n * sbilf->f_fft[data->J - 1]);

    sum     += 4.0 * P_ell * fft_real;
    phase_n *= phase;
  }

  /* Last term n = J (weight 1.0) */
  {
    const gdouble mu       = mu_over_n * data->J;
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_real = creal (phase_n * sbilf->f_fft[data->J]);

    sum += 1.0 * P_ell * fft_real;
  }

  return sum;
}

static gdouble
_ncm_sbessel_integrator_fftl_sum_odd_ell (NcmSBesselIntegratorFFTL *sbilf, const guint ell, const NcmSBesselIntegratorFFTLIntegrationData *data)
{
  gdouble sum                = 0.0;
  const complex double phase = cexp (I * (data->a + 0.5 * data->dx) * (-2.0 * M_PI / data->L));
  complex double phase_n     = phase;
  const gdouble mu_over_n    = -2.0 * M_PI / data->L;
  guint n;

  /* Loop over pairs of indices with Simpson weights (4,2) */
  /* Stop before the last two terms which need special weights (4,1) */
  for (n = 1; n + 1 < data->J; n += 2)
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
    const gdouble mu       = mu_over_n * (data->J - 1);
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[data->J - 1]);

    sum     += 4.0 * P_ell * fft_imag;
    phase_n *= phase;
  }

  /* Last term n = J (weight 1.0) */
  {
    const gdouble mu       = mu_over_n * data->J;
    const gdouble P_ell    = gsl_sf_legendre_Pl (ell, mu);
    const gdouble fft_imag = cimag (phase_n * sbilf->f_fft[data->J]);

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

  /* Load FFTW wisdom */
  ncm_cfg_load_fftw_wisdom ("ncm_sbessel_integrator_fftl");

  /* Reallocate arrays if lmax changed */
  if (lmax != sbilf->lmax_alloc)
  {
    g_clear_pointer (&sbilf->jl_arr, g_free);
    g_clear_pointer (&sbilf->Pl_arr, g_free);

    sbilf->lmax_alloc = lmax;
    sbilf->jl_arr     = g_new (gdouble, lmax + 1);
    sbilf->Pl_arr     = g_new (gdouble, lmax + 1);
  }

  sbilf->N_direct = GSL_MAX (256, 4 * lmax);
}

static guint
_ncm_sbessel_integrator_fftl_get_ell_threshold (NcmSBesselIntegratorFFTL *sbilf, NcmSBesselIntegratorFFTLIntegrationData *data)
{
  /* Return the ell value where we switch from FFT to direct Simpson's */
  /* Use FFT for ell < ell_threshold, direct Simpson's for ell >= ell_threshold */

  /* For small L_phys, just use direct Simpson's */
  if (data->L_phys < 1.0e3)
    return 0;

  /* Find the peak if not already found */
  _ncm_sbessel_integrator_fftl_find_peak (data);
  _ncm_sbessel_integrator_fftl_find_support (data);

  /* Optionally find support width for more accurate threshold */
  /* _ncm_sbessel_integrator_fftl_find_support (data); */
  /* gdouble support_width = data->x_right - data->x_left; */

  /* Use the peak location to determine threshold */
  return (guint) floor (data->x_peak);
}

static gdouble
_ncm_sbessel_integrator_fftl_integrate_ell (NcmSBesselIntegrator *sbi,
                                            NcmSBesselIntegratorF F,
                                            gdouble a, gdouble b,
                                            gint ell,
                                            gpointer user_data)
{
  NcmSBesselIntegratorFFTL *sbilf              = NCM_SBESSEL_INTEGRATOR_FFTL (sbi);
  NcmSBesselIntegratorFFTLIntegrationData data = {
    .F             = F,
    .user_data     = user_data,
    .a             = a,
    .b             = b,
    .L_phys        = b - a,
    .x_peak        = 0.0,
    .peak_val      = 0.0,
    .support_eps   = 1.0e-1,
    .x_left        = a,
    .x_right       = b,
    .peak_found    = FALSE,
    .support_found = FALSE,
    .L             = 0.0,
    .dx            = 0.0,
    .Ny            = 0,
    .J             = 0,
    .alpha_J       = 0
  };
  const guint ell_threshold = _ncm_sbessel_integrator_fftl_get_ell_threshold (sbilf, &data);
  gboolean ell_is_even;
  gboolean m_is_even;
  gdouble sum;
  guint m;

  /* For ell >= threshold, use direct Simpson's rule */
  if ((guint) ell >= ell_threshold)
    return _ncm_sbessel_integrator_fftl_simpson_direct (sbilf, F, a, b, ell, user_data);

  /* Compute integration parameters */
  _ncm_sbessel_integrator_fftl_compute_params (&data);

  /* Reallocate arrays if needed */
  if (data.Ny != sbilf->Ny)
  {
    if (sbilf->plan_forward != NULL)
      fftw_destroy_plan (sbilf->plan_forward);

    if (sbilf->f_samp != NULL)
      fftw_free (sbilf->f_samp);

    if (sbilf->f_fft != NULL)
      fftw_free (sbilf->f_fft);

    sbilf->Ny     = data.Ny;
    sbilf->f_samp = fftw_malloc (sizeof (gdouble) * data.Ny);
    sbilf->f_fft  = fftw_malloc (sizeof (fftw_complex) * (data.Ny / 2 + 1));

    ncm_cfg_load_fftw_wisdom ("ncm_sbessel_integrator_fftl");
    ncm_cfg_lock_plan_fftw ();
    sbilf->plan_forward = fftw_plan_dft_r2c_1d (data.Ny, sbilf->f_samp, sbilf->f_fft, ncm_cfg_get_fftw_default_flag ());
    ncm_cfg_unlock_plan_fftw ();
    printf ("Created FFT plan for Ny = %d\n", data.Ny);
    ncm_cfg_save_fftw_wisdom ("ncm_sbessel_integrator_fftl");
  }

  /* Sample function on grid */
  _ncm_sbessel_integrator_fftl_sample (sbilf, &data);

  /* Compute FFT */
  fftw_execute (sbilf->plan_forward);

  /* Compute Legendre sum (different for even/odd ell) */
  ell_is_even = (ell % 2 == 0);

  if (ell_is_even)
    sum = _ncm_sbessel_integrator_fftl_sum_even_ell (sbilf, ell, &data);
  else
    sum = _ncm_sbessel_integrator_fftl_sum_odd_ell (sbilf, ell, &data);

  /* Apply Simpson's rule factor and final scaling */
  m         = ell / 2;
  m_is_even = (m % 2 == 0);
  sum       = sum / 3.0;
  sum       = (m_is_even ? sum : -sum) * (2.0 * M_PI / data.Ny);

  return sum;
}

static void
_ncm_sbessel_integrator_fftl_integrate_direct (NcmSBesselIntegratorFFTL *sbilf, const guint lmin, const guint ell_direct_min, guint ell_direct_max,
                                               NcmSBesselIntegratorF F, const gdouble a, const gdouble b,
                                               NcmVector *result, gpointer user_data)
{
  const guint N                 = sbilf->N_direct;
  const gdouble dx              = (b - a) / N;
  gdouble * restrict result_ptr = ncm_vector_data (result);
  guint i, ell;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);
  /* Initialize direct results to zero */
  memset (result_ptr, 0, sizeof (gdouble) * (ell_direct_max - ell_direct_min + 1));
  ell_direct_max = GSL_MIN (ell_direct_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilf->sba, b));

  /* First term */
  {
    const gdouble fa = F (user_data, a);

    ncm_sf_sbessel_array_eval (sbilf->sba, ell_direct_max, a, sbilf->jl_arr);

    for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
    {
      result_ptr[ell - lmin] = fa * sbilf->jl_arr[ell];
    }
  }

  /* Interior terms with alternating weights 4 and 2 */
  for (i = 1; i < N; i++)
  {
    const gdouble x      = a + i * dx;
    const gdouble weight = (i % 2 == 1) ? 4.0 : 2.0;
    const gdouble fx     = F (user_data, x);

    ncm_sf_sbessel_array_eval (sbilf->sba, ell_direct_max, x, sbilf->jl_arr);

    for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
    {
      result_ptr[ell - lmin] += weight * fx * sbilf->jl_arr[ell];
    }
  }

  /* Last term */
  {
    const gdouble fb = F (user_data, b);

    ncm_sf_sbessel_array_eval (sbilf->sba, ell_direct_max, b, sbilf->jl_arr);

    for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
    {
      result_ptr[ell - lmin] += fb * sbilf->jl_arr[ell];
    }
  }

  /* Apply Simpson's rule factor */
  for (ell = ell_direct_min; ell <= ell_direct_max; ell++)
    result_ptr[ell - lmin] *= dx / 3.0;
}

static void
_ncm_sbessel_integrator_fftl_integrate_fft (NcmSBesselIntegratorFFTL *sbilf, const guint lmin, const guint ell_fft_min, const guint ell_fft_max,
                                            NcmSBesselIntegratorF F, const gdouble a, const gdouble b, const gdouble L_phys,
                                            NcmVector *result, gpointer user_data)
{
  NcmSBesselIntegratorFFTLIntegrationData data = {
    .F             = F,
    .user_data     = user_data,
    .a             = a,
    .b             = b,
    .L_phys        = L_phys,
    .x_peak        = 0.0,
    .peak_val      = 0.0,
    .support_eps   = 1.0e-3,
    .x_left        = a,
    .x_right       = b,
    .peak_found    = FALSE,
    .support_found = FALSE,
    .L             = 0.0,
    .dx            = 0.0,
    .Ny            = 0,
    .J             = 0,
    .alpha_J       = 0
  };
  NcmSBesselIntegratorFFTLCache *cache;
  gdouble * restrict result_ptr = ncm_vector_data (result);
  guint ell, n;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);

  /* Compute integration parameters */
  _ncm_sbessel_integrator_fftl_compute_params (&data);

  /* Get or create cache for this alpha_J */
  cache = _ncm_sbessel_integrator_fftl_get_cache (sbilf, data.alpha_J, ell_fft_max);

  /* Reallocate arrays if needed */
  if (data.Ny != sbilf->Ny)
  {
    if (sbilf->plan_forward != NULL)
      fftw_destroy_plan (sbilf->plan_forward);

    if (sbilf->f_samp != NULL)
      fftw_free (sbilf->f_samp);

    if (sbilf->f_fft != NULL)
      fftw_free (sbilf->f_fft);

    sbilf->Ny     = data.Ny;
    sbilf->f_samp = fftw_malloc (sizeof (gdouble) * data.Ny);
    sbilf->f_fft  = fftw_malloc (sizeof (fftw_complex) * (data.Ny / 2 + 1));

    ncm_cfg_load_fftw_wisdom ("ncm_sbessel_integrator_fftl");
    ncm_cfg_lock_plan_fftw ();
    sbilf->plan_forward = fftw_plan_dft_r2c_1d (data.Ny, sbilf->f_samp, sbilf->f_fft, ncm_cfg_get_fftw_default_flag ());
    ncm_cfg_unlock_plan_fftw ();
    ncm_cfg_save_fftw_wisdom ("ncm_sbessel_integrator_fftl");
  }

  /* Sample function on grid */
  _ncm_sbessel_integrator_fftl_sample (sbilf, &data);

  /* Compute FFT */
  fftw_execute (sbilf->plan_forward);

  /* Initialize FFT results to zero */
  for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    result_ptr[ell - lmin] = 0.0;

  /* Compute phase factor */
  const complex double phase = cexp (I * (a + 0.5 * data.dx) * (-2.0 * M_PI / data.L));

  /* n = 0 term */
  {
    const gdouble fft_real = creal (sbilf->f_fft[0]);

    for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    {
      if (ell % 2 == 0) /* even ell uses real part */
        result_ptr[ell - lmin] = cache->Pl_mu[0][ell] * fft_real;
    }
  }

  /* Loop over pairs with Simpson weights */
  complex double phase_n = phase;

  for (n = 1; n + 1 < data.J; n += 2)
  {
    /* n = odd (weight 4.0) */
    {
      for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
      {
        const gdouble P_ell = cache->Pl_mu[n][ell];

        if (ell % 2 == 0)
          result_ptr[ell - lmin] += 4.0 * P_ell * creal (phase_n * sbilf->f_fft[n]);
        else
          result_ptr[ell - lmin] += 4.0 * P_ell * cimag (phase_n * sbilf->f_fft[n]);
      }

      phase_n *= phase;
    }

    /* n+1 = even (weight 2.0) */
    {
      for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
      {
        const gdouble P_ell = cache->Pl_mu[n + 1][ell];

        if (ell % 2 == 0)
          result_ptr[ell - lmin] += 2.0 * P_ell * creal (phase_n * sbilf->f_fft[n + 1]);
        else
          result_ptr[ell - lmin] += 2.0 * P_ell * cimag (phase_n * sbilf->f_fft[n + 1]);
      }

      phase_n *= phase;
    }
  }

  /* Second to last term n = J-1 (weight 4.0) */
  {
    for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    {
      const gdouble P_ell = cache->Pl_mu[data.J - 1][ell];

      if (ell % 2 == 0)
        result_ptr[ell - lmin] += 4.0 * P_ell * creal (phase_n * sbilf->f_fft[data.J - 1]);
      else
        result_ptr[ell - lmin] += 4.0 * P_ell * cimag (phase_n * sbilf->f_fft[data.J - 1]);
    }

    phase_n *= phase;
  }

  /* Last term n = J (weight 1.0) */
  {
    for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
    {
      const gdouble P_ell = cache->Pl_mu[data.J][ell];

      if (ell % 2 == 0)
        result_ptr[ell - lmin] += 1.0 * P_ell * creal (phase_n * sbilf->f_fft[data.J]);
      else
        result_ptr[ell - lmin] += 1.0 * P_ell * cimag (phase_n * sbilf->f_fft[data.J]);
    }
  }

  /* Apply Simpson's rule factor and final scaling */
  for (ell = ell_fft_min; ell <= ell_fft_max; ell++)
  {
    const guint m         = ell / 2;
    const gboolean m_even = (m % 2 == 0);
    const gdouble val     = result_ptr[ell - lmin];

    result_ptr[ell - lmin] = (m_even ? val : -val) * (2.0 * M_PI / data.Ny) / 3.0;
  }
}

static void
_ncm_sbessel_integrator_fftl_integrate (NcmSBesselIntegrator *sbi,
                                        NcmSBesselIntegratorF F,
                                        gdouble a, gdouble b,
                                        NcmVector *result,
                                        gpointer user_data)
{
  NcmSBesselIntegratorFFTL *sbilf              = NCM_SBESSEL_INTEGRATOR_FFTL (sbi);
  const guint lmin                             = ncm_sbessel_integrator_get_lmin (sbi);
  const guint lmax                             = ncm_sbessel_integrator_get_lmax (sbi);
  const guint n_ell                            = lmax - lmin + 1;
  NcmSBesselIntegratorFFTLIntegrationData data = {
    .F             = F,
    .user_data     = user_data,
    .a             = a,
    .b             = b,
    .L_phys        = b - a,
    .x_peak        = 0.0,
    .peak_val      = 0.0,
    .support_eps   = 1.0e-3,
    .x_left        = a,
    .x_right       = b,
    .peak_found    = FALSE,
    .support_found = FALSE,
    .L             = 0.0,
    .dx            = 0.0,
    .Ny            = 0,
    .J             = 0,
    .alpha_J       = 0
  };
  const guint ell_threshold = _ncm_sbessel_integrator_fftl_get_ell_threshold (sbilf, &data); /* Find the threshold ell where we switch methods */

  g_assert_cmpuint (ncm_vector_len (result), ==, n_ell);

  /* printf ("[FFTL] Integrating from ell = %u to ell = %u with threshold at ell = %u (peak at x=%.3e)\n", lmin, lmax, ell_threshold, data.x_peak); */

  if (lmax >= ell_threshold)
  {
    const guint ell_direct_max = lmax;
    const guint ell_direct_min = GSL_MAX (lmin, ell_threshold);

    /* printf ("[FFTL]  Direct integration for ell = %u to ell = %u\n", ell_direct_min, ell_direct_max); */

    if (ell_direct_min <= ell_direct_max)
      _ncm_sbessel_integrator_fftl_integrate_direct (sbilf, lmin, ell_direct_min, ell_direct_max, F, a, b, result, user_data);
  }

  if (lmin < ell_threshold)
  {
    const guint ell_fft_min = lmin;
    const guint ell_fft_max = GSL_MIN (ell_threshold - 1, lmax);

    /* printf ("[FFTL]  FFT integration for ell = %u to ell = %u\n", ell_fft_min, ell_fft_max); */

    if (ell_fft_min <= ell_fft_max)
      _ncm_sbessel_integrator_fftl_integrate_fft (sbilf, lmin, ell_fft_min, ell_fft_max, F, a, b, data.L_phys, result, user_data);
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

